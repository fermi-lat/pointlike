"""A module for classes that perform spectral fitting.

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/specfitter.py,v 1.7 2009/07/04 15:33:22 burnett Exp $

    author: Matthew Kerr
"""
import numpy as N

class SpectralModelFitter(object):
   """A host class for the spectral fitting methods.  All called statically."""

   @staticmethod
   def least_squares(pslw,model,quiet=False):
      """Perform a least squares fit of the spectral parameters.  Appropriate for very bright sources
         and for finding seed positions for the more accurate Poisson fitter."""
      
      #Products for fit
      photons = N.fromiter( (sl.sl.photons() for sl in pslw ), int )
      alphas  = N.array([ [sl.sl.alpha(),sl.sl.sigma_alpha()] for sl in pslw ]).transpose()
      errors  = (alphas[0,:]**2*photons+alphas[1,:]**2*photons.astype(float)**2)**0.5
      signals = N.column_stack([alphas[0,:]*photons,errors]).transpose()

      def chi(parameters,*args):
         """Routine for use in least squares routine."""
         model,photons,signals,sl_wrappers = args
         model.set_parameters(parameters)
         expected = N.fromiter( (sl.expected(model) for sl in sl_wrappers ), float)
         return ((expected - signals[0,:])/signals[1,:])[photons > 0]

      from scipy.optimize import leastsq
      try:
         fit = leastsq(chi,model.get_parameters(),args=(model,photons,signals,pslw.sl_wrappers),full_output=1)
      except:
         fit = [0]*5

      if fit[4] == 1:

         model.good_fit=True #Save parameters to model
         model.cov_matrix=fit[1] #Save covariance matrix
         model.set_parameters(fit[0])
         vals = chi(fit[0],model,photons,signals,pslw.sl_wrappers)
         model.dof=len(vals)
         model.chi_sq=(vals**2).sum() #Save chi-sq information

         if not pslw.quiet and not quiet:
            print '\nFit converged!  Chi-squared at best fit: %.4f'%model.chi_sq
            print str(model)+'\n'
      return model

   @staticmethod
   def poisson(pslw,model,prefit=True):
      """Fit a spectral model using Poisson statistics."""
      
      if prefit: SpectralModelFitter.least_squares(pslw,model,quiet=True)

      def logLikelihood(parameters,*args):
         """Routine for use in minimum Poisson likelihood."""
         model = args[0]
         model.set_parameters(parameters)
         ret = sum( (sl.logLikelihood(sl.expected(model)) for sl in pslw.sl_wrappers) )
         #print 'parameters, value:' , parameters, ret
         #print 'expected: ', [sl.expected(model) for sl in pslw.sl_wrappers]
         return ret
      
      from scipy.optimize import fmin
      #fit = fmin(logLikelihood,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=2000)
      fit = fmin(logLikelihood,model.get_parameters(),args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=2000)
      warnflag = (fit[4]==1 or fit[4]==2)
      
      if not warnflag: #Good fit (claimed, anyway!)      
         try:
            from numpy.linalg import inv
            model.good_fit   = True
            model.logl       = -fit[1] #Store the log likelihood at best fit
            model.set_parameters(fit[0])
            model.set_cov_matrix(inv(SpectralModelFitter.hessian(model,logLikelihood)))
            if not pslw.quiet:
               print '\nFit converged!  Function value at minimum: %.4f'%fit[1]
               print str(model)+'\n'
         except:
            print 'Hessian inversion failed!'
      return model

   @staticmethod
   def hessian(m,mf,*args):
      """Calculate the Hessian; mf is the minimizing function, m is the model,args additional arguments for mf."""
      #p = m.p.copy()
      p  = m.get_parameters().copy()
      np = len(p)
      deltas   = N.abs(0.01 * p) #initial guess
      hessian  = N.zeros([np,np])
      bad_mask = N.asarray([False] * np)
      return_code = N.zeros(np)

      l0 = mf(p,m,*args)

      #find good values with which to estimate the covariance matrix -- look at diagonal deviations
      #iterate until change in function consistent with ~1 sigma conditional error
      for i in xrange(np):
         #print 'Working on parameter %d'%(i)
         #if p[i] < 0: deltas[i] *= -1 #necessary?
         h,l = p.copy(),p.copy()
         for j in xrange(10):
            h[:] = p[:]; l[:] = p[:];
            h[i] += deltas[i]
            l[i] -= deltas[i]

            delta_f_1 = mf(h,m,*args) - l0
            delta_f_2 = mf(l,m,*args) - l0
            delta_f = delta_f_1 + delta_f_2 #twice difference, really
            deltas[i] /= max(delta_f**0.5,0.33) # can change by half decade

            #print delta_f,delta_f_1,delta_f_2
            

            if delta_f < 5 and delta_f > 0.5: break

         if delta_f < 5e-3:
            # no constraint on parameter -- ignore it in further fittingor :
            bad_mask[i] = True
            return_code[i] = 1
            #print 'BAAAAAAD 1'
         if (delta_f_1/delta_f_2 > 10 or delta_f_1/delta_f_2 < 1./10):
            # significant asymmetry in likelihood          
            #print 'BAAAAAAD 2'
            bad_mask[i] = True
            return_code[i] = 2
         if (delta_f_2 < 5e-3 and delta_f_1 > 0.5):
            # not actually at maximum of likelihood -- upper limit condition
            #print 'BAAAAAAD 3'
            bad_mask[i] = True
            return_code[i] = 3

      #print deltas
         
      for i in xrange(np):
         if bad_mask[i]:
            hessian[i,:] = 0 #no correlation?
            hessian[:,i] = 0
            continue
         for j in xrange(i,np): #Second partials by finite difference
            
            xhyh,xhyl,xlyh,xlyl=p.copy(),p.copy(),p.copy(),p.copy()
            #xdelt = delt if p[i] >= 0 else -delt
            #ydelt = delt if p[j] >= 0 else -delt
            xdelt = deltas[i]
            ydelt = deltas[j]

            xhyh[i] += xdelt
            xhyh[j] += ydelt

            xhyl[i] += xdelt
            xhyl[j] -= ydelt

            xlyh[i] -= xdelt
            xlyh[j] += ydelt

            xlyl[i] -= xdelt
            xlyl[j] -= ydelt

            vals = N.asarray([mf(xhyh,m,*args),-mf(xhyl,m,*args),-mf(xlyh,m,*args),mf(xlyl,m,*args)])
            #if i == j:
            #print (N.abs(vals).max() - N.abs(vals).min())

            hessian[i][j]=hessian[j][i]=(mf(xhyh,m,*args)-mf(xhyl,m,*args)-mf(xlyh,m,*args)+mf(xlyl,m,*args))/\
                                          (4*xdelt*ydelt)

      mf(p,m,*args) #call likelihood with original values; this resets model and any other values that might be used later
      
      return hessian,return_code

   @staticmethod
   def gradient(m,mf,*args):
      """Calculate the gradient; mf is the minimizing function, m is the model, args additional arguments for mf."""
      p = m.get_parameters().copy()
      delt = 0.01
      gradient = N.zeros([len(p)])
      for i in xrange(len(p)):
         xh,xl = p.copy(),p.copy()
         xdelt = delt if p[i] > 0 else -delt

         xh[i]*=(1+xdelt)
         xl[i]*=(1-xdelt)

         gradient[i] = (mf(xh,m,*args)-mf(xl,m,*args))/(p[i]*2*delt)

      m.set_parameters(p)
      return gradient

def hessian(*args):
   return SpectralModelFitter.hessian(*args)