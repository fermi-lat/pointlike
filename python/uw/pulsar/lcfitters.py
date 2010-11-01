
"""
A module implementing an unbinned maximum likelihood fit of phases from some aperture
to a profile template.  The template is formed from an arbitrary number of components,
each represented by its own class.

LCPrimitives are combined to form a light curve (LCTemplate).  LCFitter then performs
a maximum likielihood fit to determine the light curve parameters.

LCFitter also allows fits to subsets of the phases for TOA calculation.

$Header:  $

author: M. Kerr <matthew.kerr@gmail.com>

"""

import numpy as np

from lcprimitives import *
from operator import mul
from scipy.special import gamma
from scipy.optimize import fmin,leastsq

SECSPERDAY = 86400.

#===============================================================================================#
class LCTemplate(object):
   """Manage a lightcurve template (collection of LCPrimitive objects).
   
      IMPORTANT: a constant background is assumed in the overall model, so there is
                 no need to furnish this separately.
   """

   def __init__(self,primitives=None,template=None):
      """primitives -- a list of LCPrimitive instances."""
      if template is not None:
         self.read_template(template)
      elif primitives is not None:
         self.primitives = primitives
      else:
         print 'No light curve components or template provided!'
         raise Exception
      self.shift_mode = np.any([p.shift_mode for p in self.primitives])

   def read_template(self,template):
      """template -- a filename with a (multi)gaussian profile."""
      self.primitives = []
      toks = [line.strip().split() for line in file(template) if len(line.strip()) > 0]
      t    = toks[0]
      if 'gauss' in t:
         self.__read_gaussian__(toks[1:])
      elif 'kernel' in t:
         self.__read_kde__(toks[1:])
      elif 'fourier' in t:
         self.__read_fourier__(toks[1:])
      else:
         print 'Template format not recognized!'
         raise Exception

   def __read_gaussian__(self,toks):
      
      for i,tok in enumerate(toks):
         if tok[0].startswith('phas'):
            g = LCGaussian()
            g.p[2] = float(tok[2])
            g.errors[2] = float(tok[4])
            self.primitives += [g]
         elif tok[0].startswith('fwhm'):
            g = self.primitives[-1]
            g.p[1] = float(tok[2])/g.fwhm_scale
            g.errors[1] = float(tok[4])/g.fwhm_scale
         elif tok[0].startswith('ampl'):
            g = self.primitives[-1]
            g.p[0] = float(tok[2])
            g.errors[0] = float(tok[4])

   def __read_kde__(self,toks):
      self.primitives += [LCKernelDensity(input_file=toks)]

   def __read_fourier__(self,toks):
      self.primitives += [LCEmpiricalFourier(input_file=toks)]      

   def set_parameters(self,p):
      start = 0
      for prim in self.primitives:
         n = int(prim.free.sum())
         prim.set_parameters(p[start:start+n])
         start += n

   def set_errors(self,errs):
      start = 0
      for prim in self.primitives:
         n = int(prim.free.sum())
         prim.errors = np.zeros_like(prim.p)
         prim.errors[prim.free] = errs[start:start+n]
         start += n

   def get_parameters(self):
      r = np.asarray([])
      for prim in self.primitives:
         r = np.append(r,prim.get_parameters())
      return r

   def set_overall_phase(self,ph):
      """Put the peak of the first component at phase ph."""
      if self.shift_mode:
         self.primitives[0].p[0] = ph
         return
      shift = ph - self.primitives[0].get_location()
      for prim in self.primitives:
         new_location = prim.get_location() + shift
         new_location = new_location - int(new_location) #wrap around if necessary
         prim.set_location(new_location)

   def get_location(self):
      return self.primitives[0].get_location()

   def norm(self):
      self.last_norm = sum( (prim.integrate() for prim in self.primitives) )
      return self.last_norm

   def __call__(self,phases,ignore_cache=False,suppress_bg=False):
      rval = np.zeros_like(phases)
      for prim in self.primitives:
         #rval += prim.cache_vals if (prim.cache and not ignore_cache) else prim(phases)
         rval += prim(phases)
      if suppress_bg: return rval/self.norm()
      else          : return (1.-self.norm()) + rval

   def gradient(self,phases):
      from collections import deque
      r = np.empty([len(self.get_parameters()),len(phases)])
      c = 0
      for prim in self.primitives:
         n = prim.free.sum()
         r[c:c+n,:] = prim.get_gradient(phases)
         c += n
      return r

   def __sort__(self):
      def cmp(p1,p2):
         if   p1.p[-1] <  p2.p[-1]: return -1
         elif p1.p[-1] == p2.p[-1]: return 0
         else: return 1
      self.primitives.sort(cmp=cmp)

   def __str__(self):
      self.__sort__()
      return '\n'+'\n\n'.join( ['P%d -- '%(i+1)+str(prim) for i,prim in enumerate(self.primitives)] ) + '\n'

   def prof_string(self,outputfile=None):
      """Return a string compatible with the format used by pygaussfit.
         Assume all primitives are gaussians."""
      rstrings = []
      dashes = '--------------------------------------'
      norm,errnorm = 0,0
      
      for nprim,prim in enumerate(self.primitives):
         phas = prim.get_location(error=True)
         fwhm = prim.get_width(error=True,fwhm=True)
         ampl = prim.get_norm(error=True)
         norm += ampl[0]
         errnorm += (ampl[1]**2)
         for st,va in zip(['phas','fwhm','ampl'],[phas,fwhm,ampl]):
            rstrings += ['%s%d = %.5f +/- %.5f'%(st,nprim+1,va[0],va[1])]
      const = 'const = %.5f +/- %.5f'%(1-norm,errnorm**0.5)
      rstring = [dashes] + [const] + rstrings + [dashes]
      if outputfile is not None:
         f = open(outputfile,'w')
         f.write('# gauss\n')
         for s in rstring:
            f.write(s+'\n')
      return '\n'.join(rstring)
       
   def random(self,n):
      """Return n random variables drawn from the distribution given
         by this light curve template.

         Implementation is by the (relatively) slow but distribution-agnostic
         accept-reject method.

         TODO: implement a multinomial approach for lcs with components
         that can implement MC generation "analytically".
      """
      if n < 1: return 0
      # find the peak of the light curve for envelope function
      seed   = np.argmax(self(np.linspace(0,1,101)))/100.
      mphase = fmin(lambda x: -self(x),[seed],full_output=1,disp=0)
      if abs(mphase[0]-seed) > 0.01:
         print 'Failed to find maximum of light curve! Aborting.'
         return
      M = -mphase[1]*1.1 # 10 percent fudge factor just in case
      
      # get n samples
      rvals    = np.empty(n)
      position = 0
      rfunc    = np.random.rand

      while True:
         cand_phases = rfunc(n)
         cand_phases = cand_phases[rfunc(n) < self(cand_phases)/M]
         ncands = len(cand_phases)
         if ncands == 0: continue
         rvals[position:position + ncands] = cand_phases[:n-position]
         position += ncands
         if position >= n:
            break
     
      return rvals

#===============================================================================================#

def get_gauss2(pulse_frac=1,x1=0.1,x2=0.55,ratio=1.5,width1=0.01,width2=0.02):
    """Return a two-gaussian template.  Convenience function."""
    n1,n2 = np.asarray([ratio,1.])*(pulse_frac/(1.+ratio))
    return LCTemplate(primitives=[LCGaussian(p=[n1,width1,x1]),
                                  LCGaussian(p=[n2,width2,x2])])

#===============================================================================================#

def get_gauss1(pulse_frac=1,x1=0.5,width1=0.01):
    """Return a one-gaussian template.  Convenience function."""
    return LCTemplate(primitives=[LCGaussian(p=[pulse_frac,width1,x1])])

#===============================================================================================#

class LCFitter(object):
   """Perform the maximum likelihood fit template to the unbinned phases."""

   def __init__(self,template,phases,times=1):
      """
         template -- an instance of LCTemplate or a file with a pref-fit Gaussian template
         phases   -- a list or array of phase values
      """
      if type(template) == type(""):
         self.template = LCTemplate(gaussian_template=template)
      else:
         self.template = template
      self.phases   = np.asarray(phases)
      self.times    = np.asarray(times) #times (MET) for phases

      self.__hist_setup__()

   def __hist_setup__(self):
      """ Setup binning for a quick chi-squared fit."""
      nbins = max(min( int( len(self.phases)/50. ) , 32),101)
      hist = np.histogram(self.phases,bins=np.linspace(0,1,nbins),new=True)
      x = ((hist[1][1:] + hist[1][:-1])/2.)[hist[0]>0]
      counts = hist[0][hist[0]>0]
      y    = counts / (x[1] - x[0]) / counts.sum()
      yerr = counts**0.5 / (x[1] - x[0]) / counts.sum()
      self.chistuff = x,y,yerr

   def loglikelihood(self,p,*args):
      if not self.template.shift_mode and np.any(p < 0):
         #guard against negative parameters
         #a better solution would be a log transform for the width/norm and
         #an arctangent transformation for the location -- later iteration!
         return 2e20
      args[0].set_parameters(p)
      return -np.log(self.template(self.phases)).sum() #negative log likelihood for minimization

   def gradient(self,p,*args):
      args[0].set_parameters(p); t = self.template
      return -(t.gradient(self.phases)/t(self.phases)).sum(axis=1)

   def chi(self,p,*args):
      x,y,yerr = self.chistuff
      if not self.template.shift_mode and np.any(p < 0):
         return 2e100*np.ones_like(x)/len(x)
      args[0].set_parameters(p)
      chi = (self.template(x,ignore_cache=True) - y)/y**0.5
      if self.template.last_norm > 1:
         return 2e100*np.ones_like(x)/len(x)
      else:
         return chi

   def toa_loglikelihood(self,p,*args):
      self.template.set_overall_phase(p[0])
      return -np.log(self.template(args[0])).sum()

   def quick_fit(self):
      f = leastsq(self.chi,self.template.get_parameters(),args=(self.template))
         
   def fit(self,quick_fit_first = False, estimate_errors = True):
      
      # an initial chi squared fit to find better seed values
      if quick_fit_first:
         self.quick_fit()

      f = fmin(self.loglikelihood,self.template.get_parameters(),args=(self.template,),disp=0,ftol=1e-6,full_output=True)
      self.fitval = f[0]
      self.ll  = -f[1]
      if estimate_errors:
         self.__errors__()
         self.template.set_errors(np.diag(self.cov_matrix)**0.5)
         return self.fitval,np.diag(self.cov_matrix)**0.5
      else:
         return self.fitval

   def fit_cg(self):
      from scipy.optimize import fmin_cg
      fit = fmin_cg(self.loglikelihood,self.template.get_parameters(),fprime=self.gradient,args=(self.template,),full_output=1,disp=1)
      return fit

   def __errors__(self):
      from numpy.linalg import inv
      h = hessian(self.template,self.loglikelihood)
      try:
         self.cov_matrix = inv(h)
      except:
         print 'Unable to invert hessian!'
         self.cov_matrix = np.zeros_like(h)

   def __toa_error__(self,val,*args):
      f  = self.toa_loglikelihood
      val= val[0]
      d  = 0.01
      d2 = (f( [(1+d)*val], *args) - 2*f([val], *args) + f( [(1-d)*val], *args))/(d*val)**2 #check
      return d2**-0.5
      

   def __str__(self):
      return '\nLog Likelihood for fit: %.2f\n'%(self.ll) + str(self.template)

   
   def find_toa(self,phases,pr=False):
      f = self.toa_loglikelihood
      fit   = fmin(f,[self.template.get_location()],args=(phases,),disp=0,ftol=1e-6,full_output=True)[0]
      err   = self.__toa_error__(fit,phases)
      return fit[0],err


   def write_template(self,outputfile='template.gauss'):
      s = self.template.prof_string(outputfile=outputfile)


   def plot(self,nbins=50,fignum=2):
      import pylab as P
      dom = np.linspace(0,1,200)
      cod = self.template(dom,ignore_cache=True)
      P.figure(fignum)
      P.hist(self.phases,bins=np.linspace(0,1,nbins+1),histtype='step',ec='red',normed=True,lw=1)
      P.plot(dom,cod,color='blue',lw=1)
      P.xlabel('Phase')
      P.grid(True)

#===============================================================================================#

class SpectralLCFitter(LCFitter):
   """ An extension of the light curve fitting to support additional information from spectroscopy."""

   def __init__(self,template,phases,ratios,times=1,suppress_bg=True):

      self.ratios = ratios
      self.suppress_bg = suppress_bg
      self.__class__.__base__.__init__(self,template,phases,times)

   def __call__(self):
      return self.template(self.phases,suppress_bg=self.suppress_bg)

   def loglikelihood(self,p,*args):
      if not self.template.shift_mode and np.any(p < 0): return 2e20
      args[0].set_parameters(p)
      t = self()
      #if np.any(t < 0): return 2e20
      return -np.log(1+self.ratios*t).sum() #negative log likelihood for minimization

   def gradient(self,p,*args):
      args[0].set_parameters(p)
      f = self()
      return (self.ratios*self.template.gradient(self.phases)/(1+self.ratios*self())).sum(axis=1)

   def fit_cg(self):
      from scipy.optimize import fmin_cg
      fit = fmin_cg(self.loglikelihood,self.template.get_parameters(),fprime=self.gradient,args=(self.template,),full_output=1,disp=1)
      return fit

   def test_statistic(self):
      ts = 2*(-np.log(1 + self.ratios).sum() - self.loglikelihood(self.template.get_parameters(),self.template))
      from scipy.stats import chi2
      return ts,chi2.sf(ts,len(self.template.get_parameters()))

#===============================================================================================#

class TaylorMapper(object):
   """Implement the mapping of time to phase.  Initial version is a simple
      Taylor expansion.  All units are assumed to be seconds.
      
      The parameters are F, Fdot, Fdotdot, etc."""
  
   def __init__(self,p,epoch,**kwargs):
      #self.init()
      self.__dict__.update(kwargs)
      self.p     = np.asarray(p)
      self.epoch = epoch
      self.free  = np.asarray([False]*len(p))

   def __call__(self,times,mod=True,apply_epoch=True):
      #times = (times/SECSPERDAY + self.mjdref - self.period_epoch)*SECSPERDAY
      if apply_epoch:
         times = times - self.epoch
      p = self.p
      rvals = times * p[0]
      if len(p) == 1:
         return np.mod(rvals,1) if mod else rvals
      for i in xrange(1,len(p)):
         rvals += p[i] * times**(i + 1) / gamma(i + 2)
      return np.mod(rvals,1) if mod else rvals

   def frequency(self,times,delta_t = 0.001):
      """Return the frequency at the provided times by calculating a numerical derivative."""

      # do this carefully to preserve maximum precision for the typical case
      times_int = times.astype(int)
      times_dec = times - times_int

      pe_int = int(self.epoch)
      pe_dec = self.epoch - pe_int

      times1  = (times_int - pe_int) + (times_dec - pe_dec)
      times2  = (times_int - pe_int) + (times_dec + delta_t - pe_dec)

      deltaphi = self(times2,apply_epoch=False) - self(times1,apply_epoch=False)
      deltaphi += (deltaphi < 0) # account for phase-wrapping

      return deltaphi/delta_t


   def get_parameters(self)   : return self.p[self.free]
   def set_parameters(self,p) : self.p[self.free ]

#===============================================================================================#

class FullFitter(object):
   """A meta-class to handle the simultaneous fitting of a light curve and a
      mapping of time to phase (period, period derivative, etc.)."""

   def __init__(self):
      pass

#===============================================================================================#

def hessian(m,mf,*args):
   """Calculate the Hessian; mf is the minimizing function, m is the model,args additional arguments for mf."""
   p = m.get_parameters().copy()
   delt = 0.01
   hessian=np.zeros([len(p),len(p)])
   for i in xrange(len(p)):
      for j in xrange(i,len(p)): #Second partials by finite difference; could be done analytically in a future revision
         
         xhyh,xhyl,xlyh,xlyl=p.copy(),p.copy(),p.copy(),p.copy()
         xdelt = delt if p[i] >= 0 else -delt
         ydelt = delt if p[j] >= 0 else -delt

         xhyh[i]*=(1+xdelt)
         xhyh[j]*=(1+ydelt)

         xhyl[i]*=(1+xdelt)
         xhyl[j]*=(1-ydelt)

         xlyh[i]*=(1-xdelt)
         xlyh[j]*=(1+ydelt)

         xlyl[i]*=(1-xdelt)
         xlyl[j]*=(1-ydelt)

         hessian[i][j]=hessian[j][i]=(mf(xhyh,m,*args)-mf(xhyl,m,*args)-mf(xlyh,m,*args)+mf(xlyl,m,*args))/\
                                       (p[i]*p[j]*4*delt**2)

   mf(p,m,*args) #call likelihood with original values; this resets model and any other values that might be used later
   return hessian

#===============================================================================================#