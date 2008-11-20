"""A module for classes that perform spectral fitting.

    $Header$
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
         model.p = parameters
         expected = N.fromiter( (sl.expected(model) for sl in sl_wrappers ), float)
         return ((expected - signals[0,:])/signals[1,:])[photons > 0]

      from scipy.optimize import leastsq
      try:
         fit = leastsq(chi,model.p,args=(model,photons,signals,pslw.sl_wrappers),full_output=1)
      except:
         fit = [0]*5

      if fit[4] == 1:

         model.good_fit=True #Save parameters to model
         model.cov_matrix=fit[1] #Save covariance matrix
         model.p=fit[0]
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
         model.p = parameters
         return sum( (sl.logLikelihood(sl.expected(model)) for sl in pslw.sl_wrappers) ) #for speed
      
      from scipy.optimize import fmin
      fit = fmin(logLikelihood,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=2000)
      warnflag = (fit[4]==1 or fit[4]==2)
      
      if not warnflag: #Good fit (claimed, anyway!)      
         try:
            from numpy.linalg import inv
            model.cov_matrix = inv(SpectralModelFitter.hessian(model,logLikelihood))
            model.good_fit   = True
            model.p          = fit[0]
            model.logl       = -fit[1] #Store the log likelihood at best fit
            if not pslw.quiet:
               print '\nFit converged!  Function value at minimum: %.4f'%fit[1]
               print str(model)+'\n'
         except:
            print 'Hessian inversion failed!'
      return model

   @staticmethod
   def hessian(m,mf,*args):
      """Calculate the Hessian; f is the minimizing function, m is the model,args additional arguments for mf."""
      delt=0.01
      p = m.p.copy()
      hessian=N.zeros([len(p),len(p)])
      for i in xrange(len(p)):
         for j in xrange(i,len(p)): #Second partials by finite difference
            
            xhyh,xhyl,xlyh,xlyl=p.copy(),p.copy(),p.copy(),p.copy()

            xhyh[i]*=(1+delt)
            xhyh[j]*=(1+delt)

            xhyl[i]*=(1+delt)
            xhyl[j]*=(1-delt)

            xlyh[i]*=(1-delt)
            xlyh[j]*=(1+delt)

            xlyl[i]*=(1-delt)
            xlyl[j]*=(1-delt)

            hessian[i][j]=hessian[j][i]=(mf(xhyh,m,*args)-mf(xhyl,m,*args)-mf(xlyh,m,*args)+mf(xlyl,m,*args))/\
                                          (p[i]*p[j]*4*delt**2)

      m.p = p #Restore parameters
      return hessian