from os.path import expandvars

import numpy as N

from astopy.io import fits as PF
from types import FunctionType,MethodType
from skymaps import SkyDir
from fitstools import rad_extract

TWOPI = 2*N.pi

def phase_cut(eventfile,outputfile=None,phaseranges=[[0,1]],phase_col_name='PULSE_PHASE'):
    """Select phases within a set of intervals.
    
        outputfile  - set to change the default output (eventfile_PHASECUT.fits)
        phaseranges - a set of ranges on which to make inclusive cuts
        
        NB -- there was a problem with using the mask as below.  Seems to be resolved,
              but check for consistency in output FITS file to be safe."""

    from numpy import array
    ef = PF.open(expandvars(eventfile))
    ph = array(ef['EVENTS'].data.field(phase_col_name)).astype(float)
    mask = array([False]*len(ph))
    
    for r in phaseranges:
        for i,myph in enumerate(ph):
            if (r[0]<= myph) and (myph <= r[1]): mask[i]=True
            
    duty_cycle = sum( (x[1] - x[0] for x in phaseranges) )
    print ('Selecting %d / %d photons (duty cycle = %.2f)'%(mask.sum(),len(mask),duty_cycle))

    hdu = PF.new_table(ef['EVENTS'].columns,nrows=mask.sum())
    for i in range(len(ef['EVENTS'].columns)):
        hdu.data.field(i)[:]=ef['EVENTS'].data.field(i)[mask]
    ef['EVENTS'].data=hdu.data
    if outputfile:
        ef.writeto(outputfile,clobber=True)
    else:
        ef.writeto(eventfile.replace('.fits','_PHASECUT.fits'),clobber=True)
    ef.close()

def phase_ltcube(ltcube,outputfile,phase):
    """ Multiply the ltcube by the phase fraction. This is useful for studying
        the DC component in the off pulse window.

        N.B. phase can be a float between 0 and 1 or must be an 
        instance of uw.pulsar.phase_range.PhaseRange

        Also, the 'EXPOSURE' and 'WEIGHTED_EXPOUSRE' table's header
        is updated so PHASE=the weighted phase fraction

        Disclaimer: I have no prove that this is the correct way to weight
        an exposure cube, but it feels right to me (and has been tested
        fairly well) --Josh
    """
    from uw.pulsar.phase_range import PhaseRange
    ltcube = PF.open(expandvars(ltcube))
    for table in ['exposure', 'weighted_exposure']:
        cb=ltcube[table].data.field('cosbins')
        fraction=phase.phase_fraction if isinstance(phase,PhaseRange) else phase
        cb*=fraction
        ltcube[table].header.update('PHASE',fraction)

    ltcube.writeto(outputfile,clobber=True)
    ltcube.close()


def constant_count_histogram(phases,photons_per_bin=100):
   """Return a set of bins and rates for a 'constant count' lightcurve.
      phases -- a list of photon phases"""
   from collections import deque
   counter = 0
   ev = N.asarray(phases)
   ev.sort()

   edges = deque()
   edges.append(0)
   for i,j in enumerate(ev):
      if counter == photons_per_bin-1:
         edges.append(j)
         counter = 0
      else: counter += 1
   edges.append(1)
   edges = N.asarray(edges)
   delta_phi = edges[1:]-edges[:-1]
   rates = photons_per_bin/delta_phi
   rates[-1] = counter/delta_phi[-1]
   return [edges[:-1],rates,delta_phi] #Left edges, rates, and widths, appropriate for a bar chart

class PulsarLightCurve(object):

   def init(self):
      self.cookie_cutter_radius = None
      self.energy_range         = [100,2e5]
      self.phase_col_name       = 'PULSE_PHASE'
      self.psf_irf              = 'P6_v3_diff'
      self.max_radius           = [5.,5.] #front, back
      self.percentage           = 68.
      self.tyrel                = False
      self.damien               = False
      self.psr_cat_ap           = -1

   def __init__(self,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.set_radius_function()

   def set_radius_function(self):

      if self.damien:
         def r(e,event_class=0):
            return N.minimum((e/1000.)**-0.85,2.2)
      elif self.tyrel:
         def r(e,event_class=0):
            return N.where(event_class,4.9,2.8)*(e/100)**-0.75
      elif self.psr_cat_ap > 0:
         def r(e,event_class=0):
            return N.maximum(N.minimum(self.psr_cat_ap,0.8*(1000./e)**0.75),0.35)
      elif self.cookie_cutter_radius is not None:
         def r(e,event_class=0):
            return N.where(event_class,self.cookie_cutter_radius[1],self.cookie_cutter_radius[0])

      else:
         from psf import PSF         
         psf_obj = PSF(psf_irf = self.psf_irf)

         b = N.logspace(N.log10(self.energy_range[0]),N.log10(self.energy_range[1]))
         bin_cents = (b[:-1]*b[1:])**0.5
      
         radii = N.append([psf_obj.conf_region(e,event_class=0,percentage=self.percentage) for e in bin_cents],\
                 [psf_obj.conf_region(e,event_class=1) for e in bin_cents])

         def r(e,event_class=0):
            rads = N.zeros_like(e)
            for i in range(len(e)):
               for j in range(len(b)-1):
                  if b[j]<= e[i] and e[i] < b[j+1]:
                     rads[i] = radii[j*(1+event_class[i])]
                     break
                     
            rads = N.minimum(N.where(event_class,self.max_radius[1],self.max_radius[0]),rads)
            return rads
      self.radius_function = r

   def get_phases(self,event_files,skydir):

      self.event_files = event_files #what is this used for?   
      results = rad_extract(event_files,skydir,self.radius_function,return_cols=['PULSE_PHASE','TIME','ENERGY'])
      ens = results['ENERGY']
      mask = N.logical_and(ens >= self.energy_range[0], ens < self.energy_range[1])
      self.phases = results['PULSE_PHASE'][mask]
      self.times  = results['TIME'][mask]
      self.differences = results['DIFFERENCES'][mask]

   def get_npb(self):
      return max(int(round(len(self.phases)/30.)),20)


   def plot(self,fignum=2,show_trend=True,show_errorbars=False,outfile=None,num_per_bin=None,
            two_cycles=False,point_source=None,relative=True,**kwargs):
      if 'linecolor' not in kwargs: kwargs['linecolor'] = 'black'
      if 'errorcolor' not in kwargs: kwargs['errorcolor']='red'
      if 'fill' not in kwargs: kwargs['fill'] = False
      
      import pylab as P
      #if show_trend: P.figure(fignum,(6,10))
      if show_trend: P.figure(fignum, (10,6))
      else: P.figure(fignum)

      npb = num_per_bin or self.get_npb()

      ev_cop = self.phases.copy()
      times = (self.times-self.times.min())/(24*3600.) #in days
      phases,rates,delta_phi = constant_count_histogram(ev_cop,npb)
      rates_pre = rates.copy()
      if relative:
         rates = rates/rates[rates>0].min()

      #ax1 = P.axes([0.15,0.1,0.75,0.4]) if show_trend else P.axes()
      ax1 = P.subplot(121) if show_trend else P.axes()
      label = '%d cts/bin'%(npb)
      if kwargs['fill']:
         ax1.bar(phases,rates,N.append(phases[1:],1)-phases,alpha=1.0,fc='red',ec='red')
      ax1.step(N.append(phases,1),N.append(rates[0],rates),where='pre',color=kwargs['linecolor'],label=label)

      ax = ax1.axis()
      ax1.axis([0,1,0,ax[3]])
      
      if two_cycles:
         ax1.step(N.append(phases,1)+1,N.append(rates[0],rates),where='pre',color=kwargs['linecolor'])
         ax1.axis([0,2,0,ax[3]])
      ax1.set_ylabel('Relative Count Rate')
      ax1.set_xlabel('Phase')
      ax1.grid(b=True)
      ax1.legend(loc=0)

      if show_errorbars:
         tph = N.append(phases,1)
         delta_phi = tph[1:]-tph[:-1]
         if relative:
            errs = ((rates_pre*delta_phi)**0.5/rates_pre[rates_pre>0].min())/delta_phi
         else:
            errs = ((rates_pre*delta_phi)**0.5)/delta_phi
         tph = (tph[:-1]+tph[1:])/2.
         P.errorbar(tph,rates,yerr=errs,ls=' ',capsize=0,color=kwargs['errorcolor'])
         if two_cycles:
            P.errorbar(tph+1,rates,yerr=errs,ls=' ',capsize=0,color=kwargs['errorcolor'])

      if show_trend:
         ax2 = P.subplot(122)
         stride = max(1,int(len(ev_cop)/500.))
         stride = 1
         ax2.scatter(times[::stride],self.phases[::stride],marker='+',color='red',s=10)
         ax2.axis([0,times[::stride].max(),0,1])
         ax2.grid(b=True)
         ax2.set_xlabel('Days')
         ax2.set_ylabel('Phase',color='red')
         for tl in ax2.get_yticklabels():
            tl.set_color('red')

         #Do H-statistic trend
         ax3 = ax2.twinx()
         t0 = self.times.min()
         delta_t = (self.times.max() - t0)/10
         ts = N.asarray([t0 + i*delta_t for i in range(10+1)])
         hs = []
         for t in ts[1:]:
            hs += [h_statistic(self.phases[self.times <= t]) ]
         ts = (ts - ts[0])/(24.*3600.)
         #ax3.plot(hs,ts[1:],marker='o',color='k',label='$\mathrm{Observed}$',ls=' ',ms=5)
         ax3.plot(ts[1:],hs,marker='o',color='blue',ls='-',ms=6)
         #fit a trend line
         from scipy import polyval,polyfit
         pfit = polyfit(ts[1:],hs,1)
         #ax3.plot(polyval(pfit,ts[1:]),ts[1:],linestyle='--',color='blue',label='$\mathrm{Linear\ Trend}$',lw=2)
         ax3.plot(ts[1:],polyval(pfit,ts[1:]),linestyle='--',color='blue',lw=3)
         #ax3.set_ylabel('Integration Time (Days)')
         sig = r'>5\sigma' if hs[-1] > 50 else '%.1f\sigma'%(sig2sigma(h_sig(hs[-1])))
         if point_source is not None:
            wh,wsig = self.weighted_h(point_source)
            tit = '$\mathrm{Aperture\ H:\ %s,\ Weighted\ H:\ %.1f\sigma}$'%(sig,wsig) 
         else:
            tit = '$\mathrm{Aperture\ H:\ %s}$'%(sig) 
         ax3.axhline(hs[-1],linestyle='-',color='blue',label='$\mathrm{Aperture\ H=%.1f\ (%s)}$'%(hs[-1],sig),lw=2)
         ax3.set_ylabel('H Test Statistic',color='blue')
         ax = ax3.axis()
         ax3.axis([0,times.max(),ax[2],ax[3]])
         for tl in ax3.get_yticklabels():
            tl.set_color('blue')
         #ax3.legend(loc=0)
         ax3.set_title(tit)

      #P.xlabel('Phase')
      if outfile is not None: P.savefig(outfile)

   def h_trend(self,intervals=10,fignum=1,outfile=None,show=True,point_source=None):
      import pylab as P
      t0 = self.times.min()
      delta_t = (self.times.max() - t0)/intervals
      ts = N.asarray([t0 + i*delta_t for i in range(intervals+1)])
      hs = []
      for t in ts[1:]:
         hs += [h_statistic(self.phases[self.times <= t]) ]
      P.figure(fignum)
      ts = (ts - ts[0])/(24.*3600.)
      
      #fit a trend line
      from scipy import polyval,polyfit
      pfit = polyfit(ts[1:],hs,1)

      if point_source is not None:
         wh,wsig = self.weighted_h(point_source)
      
      P.plot(ts[1:],hs,marker='o',color='k',label='$\mathrm{Observed}$',ls=' ',ms=9)
      P.xlabel('Integration Time (Days)')
      P.ylabel('H Test Statistic')
      sig = r'>5\sigma' if hs[-1] > 50 else '%.1f\sigma'%(sig2sigma(h_sig(hs[-1])))
      P.axhline(hs[-1],linestyle='-',color='red',label='$\mathrm{Aperture\ H=%.1f\ (%s)}$'%(hs[-1],sig),lw=2)
      if point_source is not None:
         P.axhline(hs[-1],linestyle='',color='green',label='$\mathrm{Weighted\ Test:\ %.1f\sigma}$'%(wsig),lw=2)
      P.plot(ts[1:],polyval(pfit,ts[1:]),linestyle='--',color='blue',label='$\mathrm{Linear\ Trend}$',lw=2)
      P.legend(loc=0)
      P.grid(b=True)
      if outfile is not None: P.savefig(outfile)
      #if show: P.show()

   def weighted_h(self,point_source,mc_iterations=1000):
      pp = PhotonProbability(max_radius=2)
      v,diffs,phs,ens = pp.get_probs(self.event_files,point_source,return_cols=['PULSE_PHASE','ENERGY'])
      lam = 1./N.mean( [h_statistic(N.random.rand(len(v)),v) for x in range(mc_iterations)] )
      weighted_h = h_statistic(phs,v)
      sig = N.exp(-lam*weighted_h)
      h_equiv = sig2h(sig)
      sigma   = sig2sigma(sig)
      return h_equiv,sigma


def h_sig(h):
   """Convert the H-test statistic to a chance probability."""
   if h <= 23:
      return 0.9999755*N.exp(-0.39802*h)
   if h > 50:
      return 4e-8
   return 1.210597*N.exp(-0.45901*h + 0.0022900*h**2)

def sig2h(sig):
   if sig <= 4e-8: return 50.
   if sig >= h_sig(23):
      return N.log(sig/0.9999755)/(-0.39802)
   def invhsig(x,*args):
      return h_sig(x) - args[0]
   from scipy.optimize import fsolve
   return fsolve(invhsig,[30],args=(sig,))

def sig2sigma(sig):
   from scipy.special import erfc,erfcinv
   if sig > 1e-15: return erfcinv(sig)*2**0.5
   def inverfc(x,*args):
      return erfc(x/2**0.5)-args[0]
   from scipy.optimize import fsolve
   return fsolve(inverfc,[8],(sig,))

def sigma2sig(sigma):
   # this appears to handle up to machine precision with no problem
   from scipy.special import erfc
   return erfc(sigma/2**0.5)

def sigma_trials(sigma,trials):
   # correct a sigmal value for a trials factor
   if sigma < 20:
       p = sigma2sig(sigma)*trials
       if p >= 1: return 0
       return sig2sigma(p)
   else:
      # use an asymptotic expansion -- this needs to be checked!
      return (sigma**2 - 2*N.log(trials))**0.5


def z2m(phases,m=2):
   """ Return the Z_m^2 test for each harmonic up to the specified m."""

   phases = N.asarray(phases)*TWOPI #phase in radians

   n = len(phases)


   if n < 5e3:  #faster for 100s to 1000s of phases, but requires ~20x memory of alternative

      s = (N.cos(N.outer(N.arange(1,m+1),phases))).sum(axis=1)**2 +\
          (N.sin(N.outer(N.arange(1,m+1),phases))).sum(axis=1)**2

   else:

      s = (N.asarray([(N.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
          (N.asarray([(N.sin(k*phases)).sum() for k in range(1,m+1)]))**2

   return 2./n*N.cumsum(s)

def z2mw(phases,weights,m=2):
   """ Return the Z_m^2 test for each harmonic up to the specified m.

       The user provides a list of weights.  In the case that they are
       well-distributed or assumed to be fixed, the CLT applies and the
       statistic remains calibrated.  Nice!

       NB -- the phases must be uniformly distributed, i.e., have 0 mean
       and a variance of 0.5.  Then, the 2nd central moment is just
       0.5 * the expectation of the square of the weights.
    """

   phases = N.asarray(phases)*(2*N.pi) #phase in radians

   s = (N.asarray([(N.cos(k*phases)*weights).sum() for k in range(1,m+1)]))**2 +\
       (N.asarray([(N.sin(k*phases)*weights).sum() for k in range(1,m+1)]))**2

   return N.cumsum(s) / (0.5*(weights**2).sum())

def em_four(phases,m=2,weights=None):
   """ Return the empirical Fourier coefficients up to the mth harmonic."""
   
   phases = N.asarray(phases)*TWOPI #phase in radians

   n = len(phases) if weights is None else weights.sum()
   weights = 1. if weights is None else weights

   aks = (1./n)*N.asarray([(weights*N.cos(k*phases)).sum() for k in range(1,m+1)])
   bks = (1./n)*N.asarray([(weights*N.sin(k*phases)).sum() for k in range(1,m+1)])

   return aks,bks

def em_lc(coeffs,dom):

   dom = N.asarray(dom)*(2*N.pi)

   aks,bks = coeffs
   rval = N.ones_like(dom)
   for i in range(1,len(aks)+1):
      rval += 2*(aks[i-1]*N.cos(i*dom) + bks[i-1]*N.sin(i*dom))
   return rval

def weighted_z1(phases,weights):
   """ Return the Z_1^2 test with weights."""

   phases = N.asarray(phases)*(2*N.pi) #phase in radians

   a = ((N.cos(phases)*weights).sum())**2
   b = ((N.sin(phases)*weights).sum())**2

   return 2/weights.sum()*(a + b)


def h_statistic(phases,m=20):

   phases = N.asarray(phases)*(2*N.pi) #phase in radians
   n = len(phases)

   if n < 5e3:  #faster for 100s to 1000s of phases, but requires ~20x memory of alternative

      s = (N.cos(N.outer(N.arange(1,m+1),phases))).sum(axis=1)**2 +\
          (N.sin(N.outer(N.arange(1,m+1),phases))).sum(axis=1)**2

   else:

      s = (N.asarray([(N.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
          (N.asarray([(N.sin(k*phases)).sum() for k in range(1,m+1)]))**2

   return (2./n*N.cumsum(s) - 4*N.arange(0,m)).max()


def simple_h_statistic(phases):

   phases = N.asarray(phases)*(2*N.pi) #phase in radians

   n = len(phases)

   s = (N.asarray([(N.cos(k*phases)).sum() for k in range(1,21)]))**2 +\
       (N.asarray([(N.sin(k*phases)).sum() for k in range(1,21)]))**2

   return (2./n*N.cumsum(s) - 4*N.arange(0,20)).max()

def weighted_h_statistic(phases,weights,m=20):

   phases = N.asarray(phases)*(2*N.pi) #phase in radians

   n = len(phases)

   s = (N.asarray([(weights*N.cos(k*phases)).sum() for k in range(1,m+1)]))**2 +\
       (N.asarray([(weights*N.sin(k*phases)).sum() for k in range(1,m+1)]))**2

   return ((2./(weights**2).sum())*N.cumsum(s) - 4*N.arange(0,m)).max()
