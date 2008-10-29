"""A module containing classes to do spectral fitting with pointlike.

"""

import numpy as N
import math as M
from Models import *
from types import *
from scipy.optimize import leastsq,fmin,brute,fmin_powell,fsolve,brentq,fminbound,newton,brent
from scipy.integrate import quad
from scipy.stats import poisson
from numpy.linalg import inv

from math import exp

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
class EnergyBands(object):
   """Keep track of energy bands in use by source and any additional bands."""

   def __init__(self,s_bands,r_bands=None,l_bands=None):
      """s_bands -- left bin edges associated with Pointlike fit
         r_bands -- right bin edges associated with bins *above* Pointlike energies
                    The leftmost of these is interpreted as "infinity".
         l_bands -- left bin edges associated with bins *below* Pointlike energies"""
      if r_bands is None:
         r_bands=[s_bands[-1]**2/s_bands[-2]] #Set final bin width equal to previous one
      if l_bands is None:
         l_bands=[s_bands[0]**2/s_bands[1]] #Set initial width equal to first
      self.s_bands=s_bands
      self.r_bands=r_bands
      self.l_bands=l_bands
      self.joint=l_bands+s_bands+r_bands

   def __getitem__(self,slice_or_int):
      return N.array(self.joint[slice_or_int])

   def __call__(self,infinity=True):
      """Return the left band edges and (optionally) rightmost edge of last band."""
      if infinity: return N.append(self.s_bands,self.r_bands[0])
      else: return N.array(self.s_bands)

   def psf_c(self,event_class=-1):
      if event_class < 0: return N.append(self.psf_corrections,self.psf_corrections)
      else: return self.psf_corrections

   def centers(self):
      """Return the geometric bin centers."""
      b=self(infinity=True)
      return (b[:-1]*b[1:])**0.5

   def diffs(self,from_center=False):
      """Return the bin widths."""
      b=self(infinity=True)
      if from_center:
         c=self.centers()
         return N.array([c-b[:-1],b[1:]-c])
      return b[1:]-b[:-1]

   def all(self):
      return N.array(self.joint)

   def trim(self):
      return slice(len(self.l_bands),len(self.l_bands)+len(self.s_bands))

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
class Mask(object):
   """Implement a mask based on user-specified energies to exclude."""
   
   def __init__(self,front_ranges=[[30,300000]],back_ranges=[[30,300000]]):
      """Pass a list of ranges that will be allowed."""
      self.front_ranges=N.asarray(front_ranges)
      self.back_ranges=N.asarray(back_ranges)

   def __call__(self,edges,event_class=-1):
      """Return a mask suitable for the passed set of edges and event class.
         edges is interpreted as a set of left edges and a final right edge."""
      
      if event_class == 0: return self.__mask__(edges,self.front_ranges)
      if event_class == 1: return self.__mask__(edges,self.back_ranges)
      return N.append(self.__mask__(edges,self.front_ranges),self.__mask__(edges,self.back_ranges))
      
         
   def __mask__(self,edges,ranges):
      init_mask = N.array([False] * len(edges))
      for i in xrange(len(edges)):
         if N.any(N.logical_and(edges[i]>=ranges[:,0],edges[i]<ranges[:,1])):
            init_mask[i]=True
      return N.logical_or(init_mask[:-1],init_mask[1:])

class FrontOnly(Mask):   
   def __init__(self,front_ranges=[[30,300000]]):
      self.front_ranges = N.asarray(front_ranges)
      self.back_ranges = N.asarray([[0,0]])

class BackOnly(Mask):   
   def __init__(self,back_ranges=[[30,300000]]):
      self.front_ranges = N.asarray([[0,0]])
      self.back_ranges = N.asarray(back_ranges)
      


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class Unfolder(object):
   """Unfold the estimated flux density."""
   def __init__(self,source):
      self.r_matrix=source.response.r_matrix()
      inverse=inv(self.r_matrix)
      s_and_e=source.signals
      n=s_and_e.shape[0]/2
      counts=s_and_e[:n,0]+s_and_e[n:,0] #Summed counts (front+back)
      errors=s_and_e[:n,1]**2+s_and_e[n:,1]**2 #Variances in quadrature
      
      self.unfolded=N.dot(inverse,counts)
      self.covariance=N.dot(inverse*errors,inverse) #Check this at some point
      self.err=N.diag(self.covariance)**0.5

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class SpectralFitter(object):
   """Parent class for all of the fitters.  Provides access to the models that have
      been fit."""
   
   def __init__(self):
      self.models=[]

   def __getitem__(self,key):
      """Return the Source object(s) corresponding to the integer/slice 'key'.
         Alternatively, search by the name associated with the psl."""
      if len(self.models)==0:
         print 'Fitter has no models!'
         return None
      t=type(key)
      if (t is IntType or t is SliceType):
         return self.models[key]
      elif t is StringType:
         for x in self.models:
            if x.name==key:
               return x
         print key+' not found in SourceList.'
      else:
         print 'Invalid key for SourceList.  Must provide integer, slice, or valid source name.'





#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#


class LSFitter( SpectralFitter ):
   """Fit spectrum with a completely general spectral model.  Fitting is via least squares.
      This is *not* motivated by a particular probability model, although it is equivalent 
      to maximum likelihood when the photon number is large *and* the likelihood for alpha
      is regular, the latter usually following from the former, but not always."""

   def __init__( self, source, printfit=True ):
      """Input:
               source -- a valid Source instance
               model  -- a valid object from the Models module
               mask -- a boolean array; true if energy band to be included in the fit
               printfit -- verbose output during fit
               event_class -- 0=front, 1=back, -1=front+back"""
      
      self.printfit=printfit
      self.source=source
      self.models=[]
      self.remask()

   def remask(self):
      """This is a klugey way to use the same fitter for multiple fits, allowing the mask
         or event class to change."""
      mask = self.source.global_data.mask()
      self.mask=N.logical_and( N.logical_and(self.source.alphas[:,0]>3e-5, self.source.photons>2), mask)

   def fit(self,model,**kwargs):
      """Perform the least squares fit.  Not intended to be called externally."""
      
      source,mask=self.source,self.mask
      signal,errors=source.signals[mask].transpose()
      response=source.response
      
      def residuals(parameters): #Simple function to feed the least squares routine
         model.p=parameters #Adjust model parameters
         expected=response(model=model)[mask] #Recalculate expected number
         return (expected-signal)/errors #Weighted residuals

      if signal.shape[0]>1:
         try:
            fit=leastsq(residuals,model.p,full_output=1)
         except:
            fit=[0]*5
      else: fit=[0]*5

      if fit[4]==1: #1 indicates convergence of the fit

         model.good_fit=True #Save parameters to model
         model.cov_matrix=fit[1] #Save covariance matrix
         model.p=fit[0]
         model.chi_sq=(residuals(fit[0])**2).sum() #Save chi-sq information
         model.dof=sum((1 for x in mask if x))
         if self.printfit:
            print '\nFit converged!  Function value at minimum: %.4f'%model.chi_sq
            print str(model)+'\n'

      else:
         model.good_fit=False
         if self.printfit:
            print 'Did not converge!'
      self.models+=[model] #Concatenate model to the master list


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
class PoissonLikelihood(object):
   """Parent class for Conditonal and Marginal Poisson likelihoods."""

   def remask(self):
      self.__init__(self.source)

   def chi_sq(self,model):
      #Calculate a "chi_sq" statistic - return a tuple with chi_sq, dof
      mask=N.logical_and(self.source.global_data.mask(),self.source.photons>5)
      signal=self.source.signals
      expected=self.source.response(model=model)
      return ( (((signal[:,0][mask]-expected[mask])/signal[:,1][mask])**2).sum(),sum((1 for x in mask if x)) )
      

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
class MarginalPoissonLikelihood(PoissonLikelihood):
   """Calculate the likelihood by marginalizing over the signal fraction provided by Pointlike."""

   def __init__(self,source,**kwargs):
      """Pre-calculate much of the integral to save processor time during fitting."""

      self.init()
      self.__dict__.update(kwargs)

      photons = source.photons
      mask = source.global_data.mask()
      if self.emulate_unbinned: mask = N.logical_and(mask,photons>0)
      alphas,sigmas = source.alphas.transpose()

      lm = N.logical_and(mask,photons<2)
      mm = N.logical_and(N.logical_and(mask,photons>1),photons<=self.photoncut)
      hm = N.logical_and(mask,photons>self.photoncut)

      lp = N.where(alphas[lm]>=0.5,photons[lm],0)
      mp = photons[mm]
      hp = photons[hm]

      save_vars = ['source','photons','lp','mp','hp','lm','mm','hm']
      for v in save_vars: exec 'self.%s = %s'%(v,v) in locals()
      

      #Process the mid-range photons
      alphas = alphas[mm]
      sigmas = sigmas[mm]
      slikes = source.slikes[mm]
      photons = photons[mm]

      if photons.shape[0]==0: return #This is a kluge that should be properly handled
      
      #Sample from likelihood and construct Simpson's rule factors for doing quick integration
      sp=self.sampling_points
      alpha_min=self.alpha_min #avoid singularities at alpha=0
      max_alpha_ll=N.array([x() for x in slikes]) #Max value of likelihood
      simps_weights=N.append( ([1.]+[4.,2.]*(sp/2))[:-1],1.)*(1./(3.*sp))
      self.points=N.array([N.linspace(alpha_min,1,sp+1) for i in xrange(len(photons))])
      self.wls=simps_weights*N.array([\
         N.exp( N.nan_to_num(N.fromiter( (-j(a)+max_alpha_ll[i] for a in self.points[i]) , float)) )\
         for i,j in enumerate(slikes)])
      self.wls=N.transpose((self.points[:,1]-self.points[:,0])*N.transpose(self.wls))
      #self.norm=1./(self.weighted_like_sample.sum(axis=1)) #Normalize integral over alpha
      self.wls=self.wls.transpose()
      self.points=self.points.transpose()

      self.backup_points = N.linspace(self.alpha_min,1.,100)
      self.delta = (1.-self.alpha_min)/100.

   def init(self):
      self.photoncut = 20000
      self.sampling_points = 300
      self.alpha_min = 1e-3
      self.maxll = 1e8
      self.min_sec_deriv = 1e-12
      self.emulate_unbinned = True

   
   def integ(self,alphas,*args):
      """Return the (log) integrand."""
      lambdas,photons,slikes = args
      return -N.asarray([slikes[i](alphas[i]) for i in xrange(len(alphas))]) + photons*N.log(lambdas/alphas) - lambdas/alphas

   def integ_single_e(self,alphas,*args):
      """Return the integrand."""
      return N.exp((-args[2](alphas) + args[1]*N.log(args[0]/alphas) - args[0]/alphas))

   def integ_single(self,alphas,*args):
      """Return the -(log) integrand."""
      return -(-args[2](alphas) + args[1]*N.log(args[0]/alphas) - args[0]/alphas)
      
   def i_prime(self,alphas,*args):
      """Return the first derivative of the (log) integrand for marginalization over alpha."""
      lambdas,photons,slikes = args
      first_derivs = -N.array([slikes[i].derivatives(alphas[i])[0] for i in xrange(len(alphas))])
      return first_derivs - photons/alphas + lambdas/alphas**2

   #def i_prime_single(self,alphas,*args):
   #   """Return the first derivative of the (log) integrand for marginalization over alpha."""
   #   
   #   return -args[2].derivatives(alphas)[0] - args[1]/alphas + args[0]/alphas**2

   def i_dprime(self,alphas,*args):
      """Return the second derivative of the (log) integrand for marginalization over alpha."""
      lambdas,photons,slikes = args
      second_derivs = -N.array([slikes[i].derivatives(alphas[i])[1] for i in xrange(len(alphas))])
      return second_derivs + photons/alphas**2 - 2*lambdas/alphas**3

   def i_dprime_jac(self,alphas,*args):
      """Return the second derivative of the (log) integrand for marginalization over alpha."""
      lambdas,photons,slikes = args
      second_derivs = -N.array([slikes[i].derivatives(alphas[i])[1] for i in xrange(len(alphas))])
      return N.eye(len(alphas))*(second_derivs + photons/alphas**2 - 2*lambdas/alphas**3)

   #def i_dprime_single(self,alphas,*args):
   #   """Return the second derivative of the (log) integrand for marginalization over alpha."""
   #   return -args[2].derivatives(alphas)[1] + args[1]/alphas**2 - 2*args[0]/alphas**3

   def pt_dprime(self,alphas,*args):
      """Return the second derivative of the (log) integrand for marginalization over alpha."""
      lambdas,photons,slikes = args
      second_derivs = -N.array([slikes[i].derivatives(alphas[i])[1] for i in xrange(len(alphas))])
      return second_derivs

   def __call__(self,parameters,*args):
      """Return the (negative, modified -- see below) log likelihood.  Includes normalization so
         can be used for LRT and such.""" 
      
      model = args[0]
      local_vars = ['wls','points','maxll','lp','mp','hp','lm','mm','hm']
      for v in local_vars: exec '%s = self.%s'%(v,v) in locals()
      

      model.p = parameters #Update the model
      if parameters[0]<0: return maxll*len(self.photons) #No negative fluxes!

      expected = self.source.response(model=model) #Expected number for source under the model
      le,he = expected[lm],expected[hm]

      #Do integral over low photon counts (that are not too low!)
      integral = (wls*N.nan_to_num(poisson.pmf(mp,expected[mm]/points))).sum(axis=0) #Integrate prob over alpha
      integral_mask = integral<=0.0

      #Find extremum for high-photon terms
      slikes = self.source.slikes[hm]
      #seeds = N.maximum(N.minimum(he/hp,0.95),0.05)
      seeds = N.maximum(N.minimum((((self.source.alphas).transpose()[0])[hm] + he/hp)/2.,0.95),0.5)
      #alphas_0 = fsolve(self.i_prime,seeds,fprime=self.i_dprime_jac,args=(he,hp,slikes),warning=True,factor=0.01)
      #alphas_0 = fsolve(self.i_prime,seeds,args=(he,hp,slikes),warning=False)
      #alphas_0 = N.minimum(N.maximum(alphas_0,self.alpha_min),1.)
      #print alphas_0[N.logical_or(alphas_0 < 0, alphas_0 > 1)]
      #alpha_mask = N.logical_or(alphas_0 < 0, alphas_0 > 1)
      #alphas_0 = 

      #alphas_0 = N.asarray([brentq(self.i_prime_single,1e-3,1,args=(he[i],hp[i],slikes[i])) for i in xrange(len(he))])
      #alphas_0 = N.asarray([fminbound(self.integ_single,1e-3,1,args=(he[i],hp[i],slikes[i])) for i in xrange(len(he))])
      #alphas_0 = N.asarray([newton(self.i_prime_single,he[i]/hp[i],fprime=self.i_dprime_single,args=(he[i],hp[i],slikes[i])) for i in xrange(len(he))])
      #alphas_0 = N.asarray([brent(self.integ_single,args=(he[i],hp[i],slikes[i]),brack=[1e-3,seeds[i],1.]) for i in xrange(len(he))])
      n = len(he)
      alphas_0 = [0]*n
      #high_int = [0]*n
      for i in xrange(n):
         try:
            #high_int[i] = quad(self.integ_single_e,self.alpha_min,1.,args=(he[i],hp[i],slikes[i]),full_output=1)[0]
            alphas_0[i] = brent(self.integ_single,args=(he[i],hp[i],slikes[i]),brack=[self.alpha_min,seeds[i],1.])
                          
         except:
            #print 'Being a little crafty about evaluating this integral.  Also, slow.'
            #Check for maximum at the endpoints
            high_vals = -self.integ_single(.99,he[i],hp[i],slikes[i]),-self.integ_single(1.,he[i],hp[i],slikes[i])
            low_vals = -self.integ_single(1e-5,he[i],hp[i],slikes[i]),-self.integ_single(2e-5,he[i],hp[i],slikes[i])
            if high_vals[1] > high_vals[0]:
               #print high_vals,'hey hi'
               alphas_0[i] = 1.
            
            
            elif low_vals[0] > low_vals[1]:
               #print low_vals,'hey low'
               alphas_0[i] = self.alpha_min
            else:
               #print he[i],hp[i]
               #print -self.integ_single(1e-20,he[i],hp[i],slikes[i]),-self.integ_single(seeds[i],he[i],hp[i],slikes[i]),high_vals[1]
               #punt -- brute force find maximum
               al = self.backup_points[N.argmax([-self.integ_single(a,he[i],hp[i],slikes[i]) for a in self.backup_points])]
               if al < 1 and al > self.alpha_min:
                  try:
                     alphas_0[i] = brent(self.integ_single,args=(he[i],hp[i],slikes[i]),brack=[al-self.delta,al,al+self.delta])
                     #print 're-seeding worked'
                  except:
                     alphas_0[i] = al
                     #print 'going with brute force'
               else:
                  #print al,'well, that was surprising'
                  alphas_0[i] = al
               """
               s = self.source
               for j,photon in enumerate(s.photons):
                  if photon == hp[i]:
                     print photon,j
                     break
                  
               from SourcePlot import visi
               import pylab as P
               P.clf()
               visi(s,j,he[i])
               from time import sleep
               sleep(5)
               """
      """
      alphas_0 = N.asarray(alphas_0)
      
      #print alphas_0
      d2 = N.abs(self.i_dprime(alphas_0,he,hp,slikes))
      m = d2<self.min_sec_deriv
      if N.any(m):
         print 'Warning: possible loss of numerical precision in calculating second derivative of integrand!'
         print 'Setting to a (hopefully) safe value; but, beware estimated errors.'
         print 'If you see this message, you may want to adjust photoncut higher.'
         pt = self.pt_dprime(alphas_0,he,hp,slikes)[m]
         f = self.integ
         a0 = alphas_0
         h = 0.001
         fd = ( (f(a0+h,he,hp,slikes) - 2*f(a0,he,hp,slikes) + f(a0-h,he,hp,slikes))/h**2 )[m]
         d2[m] = N.abs(fd)
         print hp[m],he[m],alphas_0[m],pt,fd
      htmult = 1#N.where(N.logical_or(alphas_0 == 1,alphas_0 == self.alpha_min),0.5,1.0)
      high_terms = (htmult*(self.integ(alphas_0,he,hp,slikes)  - 0.5*N.log(d2))).sum()
      """
      high_terms = 0

      #First term is a "penalty" to mimic to log(very small number)
      logl =  maxll*(integral_mask.sum())\
         - (N.log(integral[N.logical_not(integral_mask)])).sum() \
         - (lp*N.log(le) - le).sum() \
         - high_terms
      #print parameters,high_terms,logl
      return logl

      
               

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ConditionalPoissonLikelihood(PoissonLikelihood):
   """Calculate the *log* likelihood conditional on alpha, the output of Pointfit."""

   def __init__(self,source):
      """Set up variables and cut on poorly-fit bins."""
      mask=source.global_data.mask()
      alphas,sigmas = source.alphas.transpose()
      self.mask=N.logical_and(N.logical_and(alphas>3e-5,sigmas<0.5),mask)
      self.zero_mask=alphas<0
      self.photons=source.photons[self.mask]    
      self.alphas=alphas[self.mask]
      self.sigmas=sigmas[self.mask] #Added for GLM
      self.source=source
      #Factorial terms for normalization, calculate just once
      self.norm_terms=N.fromiter( (N.log(N.arange(2,x+1)).sum() for x in self.photons) ,float).sum()
 
   def __call__(self,parameters,*args):
      """Return the (negative) log likelihood.  Includes the factorial term gof tests."""
      model=args[0]
      model.p=parameters
      expected=self.source.response(model=model)[self.mask]
      return self.norm_terms + expected[self.zero_mask].sum() +\
         (expected/self.alphas-self.photons*N.log(expected/self.alphas)).sum()

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ExtendedPoissonLikelihood(PoissonLikelihood):
   """Use the exended likelihood functionality of pointlike to fit a source spectrum."""

   def __init__(self,source,**kwargs):
      """Set up variables and cut on poorly-fit bins."""
      self.init()
      self.__dict__.update(kwargs)

      mask=source.global_data.mask()
      self.zero_mask = N.logical_and(mask,source.photons == 0)
      self.mask = N.logical_and(mask,source.photons > 0)
      self.slikes = source.slikes[self.mask]      
      self.photons = source.photons[self.mask]
      self.bg = N.array([x.background() for x in self.slikes])
      self.source = source

   def init(self):
      self.emulate_unbinned = True
 
   def __call__(self,parameters,*args):
      """Return the (negative) log likelihood."""
      model=args[0]
      model.p=parameters
      expected=self.source.response(model=model) #already accounts for PSF correction factor
      mask = self.mask
      alphas = expected[mask] / (self.bg + expected[mask])
      spatial_likes = sum( (self.slikes[i](alphas[i]) for i in xrange(len(alphas))) )
      poiss_likes = (expected[mask] + self.bg - self.photons*N.log(expected[mask] + self.bg)).sum()
      #likes = sum( (self.slikes[i].logLikelihood(expected[self.mask][i]) for i in xrange(len(self.slikes))) )
      zterms = 0 if self.emulate_unbinned else expected[self.zero_mask].sum()
      #print spatial_likes,poiss_likes,zterms
      return spatial_likes + poiss_likes + zterms

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PoissonFitter(SpectralFitter):
   """Fit spectrum with a completely general spectral model.  Fitting is done to the posterior
      mode as calculated with a uniform prior on the signal fraction and the spectral parameters.
      This is a largely correct probabilistically, and in particular works well for low photon
      counts."""

   def __init__(self, source, method='Marginal', printfit=True, lsfirst=True):

      self.printfit,self.lsfirst,self.source = printfit,lsfirst,source
      self.models=[]
      self.maxll=1e8
      self.systematics =self.one_off = False
      exec(''.join(['self.l=',method,'PoissonLikelihood(source)']))

   def remask(self):
      """Change energy band mask or event classes under consideration."""
      self.l.remask()

   def hessian(self,p,mf,m):
      """Calculate the Hessian; p is the best-fit point, f is the minimizing function, m is the model."""
      delt=0.01
      hessian=N.empty([len(p),len(p)])
      for i in xrange(len(p)):
         for j in xrange(i,len(p)): #Second partials by finite difference
            
            xhyh=p.copy()
            xhyh[i]*=(1+delt)
            xhyh[j]*=(1+delt)
            #print p

            xhyl=p.copy()
            xhyl[i]*=(1+delt)
            xhyl[j]*=(1-delt)
            #print p

            xlyh=p.copy()
            xlyh[i]*=(1-delt)
            xlyh[j]*=(1+delt)
            #print p

            xlyl=p.copy()
            xlyl[i]*=(1-delt)
            xlyl[j]*=(1-delt)
            #print p
            #print mf(xhyh,m),mf(xhyl,m),mf(xlyh,m),mf(xlyl,m)
            hessian[i][j]=hessian[j][i]=(mf(xhyh,m)-mf(xhyl,m)-mf(xlyh,m)+mf(xlyl,m))/(p[i]*p[j]*4*delt**2)
      #print 1./N.diag(hessian)**0.5/p
      #print hessian
      return hessian

   def __ls__(self,model):

      p = [x for x in model.p]
      #self.source.fit(model = model.name, p = p, e0 = model.e0, method = 'LS', printfit = False)
      self.source.fit(model = model.name, method = 'LS', printfit = False, **model.__dict__)
      if self.source.LS.models[-1].good_fit:
         model.p=[x for x in self.source.LS.models[-1].p]
         
   def fit(self,model,**kwargs):
      
      self.__dict__.update(kwargs)
      x0 = [x for x in model.p] #Make a copy of x0
      param_vals = []
      
      repetitions = 200 if self.systematics else 1
      if self.systematics and self.printfit:
         print 'Calculating systematics with %d Monte Carlo trials.  Go fix a drink, for this will take a while.'%repetitions

      for i in xrange(repetitions):

         if i > 0:
            self.lsfirst = False #Use previous value for seed, for speed
            self.source.response.update(random = True)
            printfit = self.printfit
            self.printfit = False

         #Do a simplex search for maximum likelihood, possibly using least squares estimate as seed
         if self.l.photons.shape[0]>1:
      
            if self.lsfirst:
               if self.printfit: print 'Using least squares to find seed position.'
               self.__ls__(model)
            
            #Uncomment this statement for easier debugging
            fit=fmin(self.l,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=1000)

            try:
               fit=fmin(self.l,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=1000)
               warnflag=(fit[4]==1 or fit[4]==2) or fit[1]>self.maxll/10.
               if warnflag: raise Exception

            except:
            
               if not self.lsfirst:
                  if self.printfit: print 'Did not converge from seed, trying least squares to seed.'
                  model.p=x0
                  self.__ls__(model)

                  try:

                     fit=fmin(self.l,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=1000)
                     warnflag=(fit[4]==1 or fit[4]==2) or fit[1]>self.maxll/10.

                  except: warnflag=True
               else: warnflag = True
               
         else: warnflag=True


         if not warnflag: #Good fit (claimed, anyway!)
            
            try:
               model.cov_matrix=inv(self.hessian(fit[0],self.l,model))
               model.good_fit=True
               model.p=fit[0]
               model.logl=-fit[1] #Store the log likelihood at best fit
               model.chi_sq,model.dof=self.l.chi_sq(model)
               if self.printfit:
                  print '\nFit converged!  Function value at minimum: %.4f'%fit[1]
                  print str(model)+'\n'
            except:
                  print 'Hessian inversion failed!'

         else:
            if self.printfit: print 'Fit did not converge :('

         if i == 0:
            self.models += [model]
            if self.systematics:
               x0 = [x for x in model.p] #Use future value for seeding
               logl = model.logl
               chi_sq = model.chi_sq
               cov_matrix = model.cov_matrix.copy()
            
         else:
            param_vals += [[x for x in model.p]]

      if self.systematics:
         #Put values back in -- what a kluge
         self.models[-1].p = x0
         self.models[-1].logl = logl
         self.models[-1].cov_matrix = cov_matrix

         #Now, calculate systematics
         param_vals = N.array(param_vals)
         #print param_vals
         systematics=N.empty([len(x0),3])
         contain_68 = 0.68*repetitions
         bite = int(round((repetitions-contain_68)/2.))
         for i in xrange(len(x0)):
            the_vals = N.sort(param_vals[:,i])
            systematics[i,0] = 1 - the_vals[bite]/x0[i]
            systematics[i,1] = the_vals[-bite]/x0[i] - 1
            systematics[i,2] = the_vals.std()/x0[i]
         self.models[-1].systematics = systematics
         #self.models[-1].systematics = param_vals.std(0)
         self.printfit = printfit

      if self.one_off: del(self.l) #Release memory for likelihood object; calling fit again will error!

if __name__=='__main__':
   pass
   #Write some test routines!

def cubic_roots(a,b,c):
   """Find (real) roots for x^3 + ax^2 + bx +c = 0 for real co-efficients."""

   q = (a**2 - 3*b)/9.
   r = (2*a**3 - 9*a*b + 27*c)/54.

   three_rr = r**2 < q**3
   if three_rr:
      theta = N.arccos(r/q**(3./2))
      prelims = N.array([theta/3.,(theta+N.pi*2)/3.,(theta-N.pi*2)/3])
      return -2*q**0.5*N.cos(prelimes) - a/3.
   
   multi = 1 if r < 0 else -1 #-sgn(R)
   A = multi * ( N.abs(r) + (r**2-q**3)**0.5 )**(1./3)
   B = q/A if a != 0 else 0
   return A+B-a/3.

   #for reference, the complex roots are -0.5*(A+B)-a/3 +/- i*3**0.5/2*(A-B)
