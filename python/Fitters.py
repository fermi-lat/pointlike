import numpy as N
import math as M
import pointlike as pl
from Models import *
from types import *
from scipy.optimize import leastsq,fmin,brute,fmin_powell
from scipy.integrate import quad
from scipy.stats import poisson
from numpy.linalg import inv

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
      self.set_psf_c()

   def __getitem__(self,slice_or_int):
      return N.array(self.joint[slice_or_int])

   def __call__(self,infinity=True):
      """Return the left band edges and (optionally) rightmost edge of last band."""
      if infinity: return N.append(self.s_bands,self.r_bands[0])
      else: return N.array(self.s_bands)

   def set_psf_c(self,umax=50):
      gamma = 2.25
      pl.PointSourceLikelihood.setDefaultUmax(umax)
      self.psf_corrections = N.array([ 1-(1+umax/gamma)**(1-gamma) ]*len(self.s_bands))
   
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
      mask=self.mask
      signal=self.source.signals
      expected=self.source.response(model=model)
      return ( (((signal[:,0][mask]-expected[mask])/signal[:,1][mask])**2).sum(),sum((1 for x in mask if x)) )
      

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
class MarginalPoissonLikelihood(PoissonLikelihood):
   """Calculate the likelihood by marginalizing over the signal fraction provided by Pointlike."""

   def __init__(self,source,maxll=1e8):
      """Pre-calculate much of the integral to save processor time during fitting."""

      photons = source.photons
      mask = source.global_data.mask()
      zero_mask = N.logical_and(mask,photons==0)
      ones_mask = N.logical_and(mask,photons==1)
      #mask = N.logical_and(mask,photons>1)
      mask = N.logical_and(mask,photons>0)

      alphas,sigmas = source.alphas[mask].transpose()
      slikes = source.slikes[mask]
      
      self.photons,self.mask,self.source,self.maxll,self.zero_mask,self.ones_mask = \
         photons,mask,source,maxll,zero_mask,ones_mask
     
      if self.photons.shape[0]==0: return

      photons=photons[mask]
      if photons.shape[0]==0: return
      alpha_min=1e-2 #avoid singularities at alpha=0

      #Sample from likelihood and construct Simpson's rule factors for doing quick integration
      sampling_points=1000
      max_alpha_ll=N.array([x() for x in slikes]) #Max value of likelihood
      simps_weights=N.append( ([1.]+[4.,2.]*(sampling_points/2))[:-1],1.)*(1./(3.*sampling_points))
      self.points=N.array([N.linspace(alpha_min,1,sampling_points+1) for i in xrange(len(photons))])
      self.weighted_like_sample=simps_weights*N.array([\
         N.exp( N.nan_to_num(N.fromiter( (-j(a)+max_alpha_ll[i] for a in self.points[i]) , float)) )\
         for i,j in enumerate(slikes)])
      self.weighted_like_sample=N.transpose((self.points[:,1]-self.points[:,0])*N.transpose(self.weighted_like_sample))
      self.norm=1./(self.weighted_like_sample.sum(axis=1)) #Normalize integral over alpha
      self.weighted_like_sample=self.weighted_like_sample.transpose()
      self.points=self.points.transpose()
      self.alphas = alphas



   def __call__(self,parameters,*args):
      """Return the (negative, modified -- see below) log likelihood.  Includes normalization so
         can be used for LRT and such."""
      
      wls,photons,zero_photons,points,model,maxll,norm = \
         self.weighted_like_sample,self.photons[self.mask],self.photons[self.zero_mask],self.points,\
         args[0],self.maxll,self.norm
      model.p=parameters #Update the model
      if parameters[0]<0: return maxll*len(photons)
      expected=self.source.response(model=model) #Expected number for source under the model
      z_expected=expected[self.zero_mask].sum() #Poisson prob for 0 photons
      o_expected=expected[self.ones_mask]
      alpha_expected = expected[self.mask]/points #Divide by alpha to get total expected number
      integral = norm*(wls*N.nan_to_num(poisson.pmf(photons,alpha_expected))).sum(axis=0) #Integrate prob over alpha
      integral_mask = integral<=0.0
      #if N.any(integral_mask): print photons[integral_mask]

      #First term is a "penalty" to mimic to log(very small number)
      return maxll*integral_mask.sum() \
         - (N.log(integral[N.logical_not(integral_mask)])).sum() + z_expected\
         - (N.log(o_expected) - o_expected).sum()

      
               

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
      
      #Experimental GLM terms
      #return ((expected-self.photons*self.alphas*N.log(expected))/(1+10*self.sigmas**2)).sum()
      #return ((expected-self.photons*self.alphas*N.log(expected))/ \
      #   (self.alphas+self.photons*self.sigmas**2/self.alphas)).sum()

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
         for j in xrange(len(p)): #Second partials by finite difference
            
            xhyh=p.copy()
            xhyh[i]*=(1+delt)
            xhyh[j]*=(1+delt)

            xhyl=p.copy()
            xhyl[i]*=(1+delt)
            xhyl[j]*=(1-delt)

            xlyh=p.copy()
            xlyh[i]*=(1-delt)
            xlyh[j]*=(1+delt)

            xlyl=p.copy()
            xlyl[i]*=(1-delt)
            xlyl[j]*=(1-delt)

            hessian[i][j]=(mf(xhyh,m)-mf(xhyl,m)-mf(xlyh,m)+mf(xlyl,m))/(p[i]*p[j]*4*delt*delt)
      return hessian

   def __ls__(self,model):

      p = [x for x in model.p]
      self.source.fit(model = model.name, p = p, e0 = model.e0, method = 'LS', printfit = False)
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
    
            try:
               fit=fmin(self.l,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=1000)
               warnflag=(fit[4]==1 or fit[4]==2) or fit[1]>self.maxll/10.
               if warnflag: raise Exception

            except:
            
               if self.printfit: print 'Did not converge from seed, trying least squares to seed.'
               model.p=x0
               self.__ls__(model)

               try:

                  fit=fmin(self.l,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=1000)
                  warnflag=(fit[4]==1 or fit[4]==2) or fit[1]>self.maxll/10.

               except: warnflag=True
               
         else: warnflag=True


         if not warnflag: #Good fit (claimed, anyway!)
            
            try: model.cov_matrix=inv(self.hessian(fit[0],self.l,model))
            except: warnflag=True

         if not warnflag: #Good fit and successful inversion of Hessian

            model.good_fit=True
            model.p=fit[0]
            model.logl=-fit[1] #Store the log likelihood at best fit
            model.chi_sq,model.dof=self.l.chi_sq(model)
            if self.printfit:
               print '\nFit converged!  Function value at minimum: %.4f'%fit[1]
               print str(model)+'\n'

         else:
            if self.printfit: print 'Fit did not converge :('
            model.good_fit=False

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
