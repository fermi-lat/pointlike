import numpy as np
from skymaps import SkyDir
from uw.utilities.fitstools import rad_extract
from uw.like.pypsf import OldPsf,PsfOverlap
from lcprimitives import LCGaussian,LCHarmonic
from lcfitters    import LCFitter,LCTemplate,SpectralLCFitter
from scipy.integrate import simps
from scipy.stats import chi2
from scipy.special import gamma
from cPickle import dump

TWOPI = 2*np.pi

class EdepAperture(object):

   def init(self):
      self.flaperture = np.radians(1)

   def __init__(self):
      pass

   def flat_aperture(en,ct):
      return self.flaperture

class LCLikelihood(object):
   """ Implement a likelihood test using all available information about a set of photons.
       Uses the spectral fit from an ROIAnalysis to calculate the source and background rates
       for each photon.  Then supports a search in phase space for 1- and 2-peak light curves
       to find a maximum likelihood.

       Twice the difference in log likelihood is distributed as chi2 with 2(5) degrees of
       freedom for a 1(2)-peak light curve.

       The idea is that by incorporating all of the available information, the background
       "rejection" is improved and the sensitivity to pulsation increased.

       Additionally, since the test can be made prescriptive (by, e.g., fixing the shape
       of the light curves test), the sensitivity can be directly related to the DC
       sensitivity.  (This could be a lie -- just thinking out loud.)

       THIS COULD ALSO BE USED IN A BLIND SEARCH ALGORITHM!!!

   """

   def init(self):
     
      self.phase_col_name  = 'PULSE_PHASE'
      self.scan_width      = 0.03 # width of gaussian to use in phase scan
      self.scan_resolution = 0.01  # steps in phase to take in phase scan
      self.period          = None
      self.source_number   = 0
      self.flaperture      = 1
      self.edepaperture    = None
      self.psrcataperture  = None
      self.prms            = 0
      self.eventclass      = 3
      self.cuts            = None # cuts on FT1 columns in radial extraction

      self.emin = 1e2
      self.emax = 1e4

   def __init__(self,roi,ft1files=None,**kwargs):

      self.init()
      self.__dict__.update(kwargs)
      self.roi = roi
      self.exp = roi.sa.exposure.exposure

      # if ft1file is None, try to get it from roi
      ft1files = ft1files if ft1files is not None else roi.sa.pixeldata.ft1files

      # in principle, might want to vary the location in a blind search, but the sensitivity
      # of the pulsed photons to position is much greater than the DC likelihood
      skydir = roi.psm.point_sources[self.source_number].skydir

      # instantiate a PSF using the IRF/livetime of the ROI object
      #psf = OldPsf(roi.sa.CALDB,irf=roi.sa.irf)
      psf = roi.sa.psf
      psf.set_weights(roi.sa.ltcube,skydir)
      self.psf = psf

      # one can specify a pulsar catalog-like aperture which takes precedent
      if self.psrcataperture is not None:
         def r(e,event_class=0):
            return np.maximum(np.minimum(self.psrcataperture,0.8*(1000./e)**0.75),0.35)
      # alternatively, one can pass a function taking as arguments energy and event conversion type
      elif self.edepaperture is not None:
         r = self.edepaperture
      # finally, the default is a cookie cutter aperture
      else:
         r = self.flaperture

      return_cols = ['TIME','CONVERSION_TYPE','EVENT_CLASS'] + ([self.phase_col_name] if self.period is None else [])
      edict = rad_extract(ft1files,skydir,r,return_cols=return_cols,cuts=self.cuts)
      mask = np.logical_and(edict['ENERGY'] > self.emin,edict['ENERGY'] < self.emax)
      mask = np.logical_and(edict['EVENT_CLASS'] >= self.eventclass, mask)
      for k,v in edict.iteritems(): edict[k] = v[mask]
      edict['CONVERSION_TYPE'] = edict['CONVERSION_TYPE'].astype(int)
      
      self.edict  = edict # save values for now for debugging; don't need in likelihood algorithm
      self.signals,self.backs = self.calc_ratios(edict,psf)
      self.weights,self.ratios = self.signals / (self.signals + self.backs), self.signals / self.backs
      self.phases = edict[self.phase_col_name] if self.period is None else np.mod((edict['TIME']-edict['TIME'].min())/self.period,1)
      self.times  = edict['TIME']

      self.null   = np.log(1+self.ratios).sum() # value of log likelihood under null (no pulsation) hypothesis

   def pickle(self,filename):
      self.roi = None
      self.exp = None
      self.psf = None
      dump(self,file(filename,'wb'))

   def calc_ratios(self,edict,psf):
      """ edict is a dictionary containing the event information.
          Calculate for each photon the ratio of the source to the background.
      """

      en = edict['ENERGY']
      ct = edict['CONVERSION_TYPE']
      sn = self.source_number

      psm = self.roi.psm
      bgm = self.roi.bgm

      dirs = map(SkyDir,edict['RA'],edict['DEC'])
      
      signals = np.empty(len(en))
      backs   = np.empty(len(en))
      #ratios  = np.empty(len(en))
      #weights = np.empty(len(en))

      for i,(e,c,d) in enumerate(zip(en,ct,dirs)):
         
         psfluxes  = np.asarray([m(e) for m in psm.models])
         psdiffs   = np.asarray([d.difference(ps.skydir) for ps in psm.point_sources])
         psrates   = psfluxes * psf(e,c,psdiffs,density=True)

         bgrates   = np.asarray([m.get_dmodel(c).value(d,e)*m.smodel(e) for m in bgm.bgmodels])

         signals[i] = psrates[sn]
         backs[i]   = (psrates.sum() - psrates[sn] + bgrates.sum())

      return signals,backs

   def calc_integrals(self,edict,psf):
      """ Calculate the integral over the ROI of the model(s) """

      class SimpleBand(object):
         def __init__(self,en,ct,lcl):
            self.e = en; self.ct = ct
            self.psf = psf.band_psf(self)
            self.radius_in_rad = np.radians(lcl.flaperture)

      m  = self.roi.psm.models[self.source_number]
      md = self.roi.psm.point_sources[self.source_number].skydir
      rd = self.roi.psm.roi_dir
      po = PsfOverlap()

      # energy sampling points
      pts = np.logspace(np.log10(self.emin),np.log10(self.emax), 17)
      flu = np.asarray([m(e) for e in pts])
      ovl = np.asarray([ [po(SimpleBand(en,ct,self),rd,md) for ct in [0,1]] for en in pts ])
      exp = np.asarray([ [self.exp[ct].value(md,en)        for ct in [0,1]] for en in pts ])
      flu = flu * (ovl*exp).sum(axis=1)

      self.source_integral = simps(flu,x=pts)

def phase_scan(self,verbose=True,max_peaks=3,null_test=False,harmonics=None):

   r   = self.ratios
   ph  = np.random.rand(len(self.phases)) if null_test else self.phases

   gaussian = harmonics is None

   lct  = LCTemplate([]) #constant function
   lct.shift_mode = True
   lcf = SpectralLCFitter(lct,ph,r,1.,suppress_bg=False)

   dom = np.arange(0,1,self.scan_resolution)
   ts  = np.empty_like(dom)

   def add_peak():

      num_peaks = len(lct.primitives) + 1

      norm = 0.5 if num_peaks==1 else lct.norm()/num_peaks

      if gaussian:
         g1  = LCGaussian(p = [norm,self.scan_width,0])
      else:
         g1  = LCHarmonic(order=num_peaks,p=[norm,0.])
      lct.primitives += [g1]

      if lct.norm () >= 1:
         n = lct.norm()
         for p in lct.primitives:
            p.p[0] /= n


      # do a grid search first to find best peak candidate
      for i,p in enumerate(dom):
         lct.primitives[-1].set_location(p)
         ts[i] = lcf.loglikelihood(lct.get_parameters(),lct)

      import pylab as P
      P.plot(dom,ts)

      sig = ts.max() - ts.min()
      if verbose: print 'Delta Log L from scan: %.2f'%(sig)

      g1.set_location(dom[np.argmin(ts)])
      if verbose: print 'Peak trial position found at %.2f'%(g1.get_location())
      g1.free[0] = not(num_peaks==1); g1.free[1] = True;
      if gaussian: g1.free[2]  = False
      
      #if num_peaks > 1:
      # tune new peak with others fixed    
      lcf.fit()
      if verbose:
         print '\nFinishing tuning peak:'
         print g1

      # finally, re-adjust all params
      for prim in lct.primitives:
         prim.free[:] = True
         #if num_peaks==1 and gaussian: prim.free[0] = False
      lcf.fit()
      if verbose:
         print '\nFull refit results:'
         print lct
      prim.free[:] = True
      #if g1.p[2] < -1 or g1.p[2] > 2 or g1.p[1] > 2:
      #   lcf.bad_fit = True
      #else:
      #   lcf.bad_fit = False
      #print g1.p
      return sig

   for i in xrange(max_peaks):
      sig = add_peak()
      #if sig < 50 and not (quick_fit and i==0): break

   if null_test:
      if lcf.bad_fit: return -1
      ts = lcf.test_statistic()[0]
      #from scipy.integrate import quad
      #print quad(lct,0,1,args=(False,True))
      #print ts
      return ts
   return lcf.test_statistic()[0] if null_test else lcf

def null_values(self,max_peaks=1,num_trials=10):
   nvals = np.zeros(num_trials)
   for i in xrange(num_trials):
      lcf = phase_scan(self,max_peaks=max_peaks,null_test=True,verbose=False)
      nvals[i] = lcf.test_statistic()[0]
      print nvals[i]
   return nvals

def get_wrapper(filename):
   phases,ratios = np.asarray([line.strip().split() for line in file(filename)]).astype(float).transpose()

   class Wrapper(object):
      def __init__(self,phases,ratios):
         self.ratios = ratios; self.phases = phases
         self.scan_width      = 0.04
         self.scan_resolution = 0.01

   return Wrapper(phases,ratios)

def light_curve(phases,weights=None,bins=25,ec='blue',ls='solid',label=None,axes=None,gimme=False):
   
   import pylab as P
   import numpy as N
   from collections import deque
   axes = P.gca() if axes is None else axes

   mc   = deque()
   bins = np.linspace(0,1,bins+1)
   n    = len(phases)

   for i in xrange(200):
      #cod = np.histogram(phases[np.argsort(np.random.rand(n))],bins=bins,weights=weights,new=True,normed=True)[0]
      cod = np.histogram(np.random.rand(n),bins=bins,weights=weights,new=True,normed=True)[0]
      mc.append(cod)

   err = np.concatenate(mc).std()
   print err

   cod,dom = np.histogram(phases,bins=bins,weights=weights,new=True,normed=True)
   dom = np.append(dom[:-1],dom+1)
   cod = np.append(cod,cod)
   x = np.zeros(len(dom) * 2 )
   y = np.zeros(len(dom) * 2 )
   x[0::2],x[1::2] = dom,dom
   y[1:-1:2],y[2::2] = cod,cod
   axes.fill(x,y,closed=False,fill=False,edgecolor=ec,ls=ls,label=label)

   dom = (dom[1:] + dom[:-1])/2.
   axes.errorbar(dom,cod,yerr=err,color=ec,capsize=0,ls=' ',marker=' ')

   if gimme: return dom,cod,err

def compplot(phases,nbinsmin=5,nbinsmax=101,compsig=None):
   import pylab as P
   from uw.utilities.phasetools import sig2sigma,z2m
   from scipy.stats import chi2

   l10 = np.log(10)

   ax1 = P.subplot(121)
   #ax1.set_yscale('log')
   ax2 = P.subplot(122, sharey=ax1)

   chidom = np.arange(nbinsmin,nbinsmax+1,5)
   chisig = np.empty_like(chidom).astype(float)
   for i,nbins in enumerate(chidom):
      h = np.histogram(phases,bins=np.linspace(0,1,nbins+1),new=True)[0]
      m = h.mean()
      ts=((h-m)**2/m).sum()
      chisig[i]= ts2sig(ts,nbins-1)

   z2dom   = np.arange(1,21)
   z2mvals = z2m(phases,m=20)
   z2sigs  = np.asarray([ts2sig(z2mvals[m],2*(m+1)) for m in xrange(len(z2mvals))])

   ax1.plot(chidom,chisig,marker='o',ms=5,color='green')
   ax1.set_title(mathrm('\chi^2 Test'))
   ax1.set_xlabel(mathrm('number of bins'))
   ax1.set_ylabel(mathrm('Significance (\sigma units)'))
   if compsig is not None:
      for ics,cs in enumerate(compsig):
         ax1.axhline(cs,color=['blue','green','red'][ics%3],label=mathrm('%d-peak Spectral'%(ics)),lw=2)
   #ax1.legend(loc='center')
   ax1.grid()
   
   #ax2.set_yscale('log')
   ax2.plot(z2dom,z2sigs,marker='o',ms=5,color='green')
   ax2.set_title(mathrm('Z^2_m and H-test'))
   #ax2.axhline(hsig,color='red',label=mathrm('H Test upper limit'),lw=2)
   ax2.axhline(z2sigs[1],color='green',ls='--',label=mathrm('Z^2_2 Test'))
   ax2.set_xlabel(mathrm('numer of harmonics'))
   if compsig is not None:
      for ics,cs in enumerate(compsig):
         ax2.axhline(cs,color=['blue','green','red'][ics%3],label=mathrm('%d-peak Spectral'%(ics)),lw=2)
   #ax2.legend(loc='center')
   ax2.grid()

def mathrm(string):
   return r'$\mathrm{'+string.replace(' ','\ ')+'}$'

def logchi2sf(test_statistic,dof):
   
   dof = np.asarray(dof).astype(float)
   ts  = np.asarray(test_statistic).astype(float)

   p = chi2.sf(ts,dof)
   try:
      mask = p>0
      p[mask]= np.log(p)
      nmask = np.logical_not(p)
      p[nmask] = -ts[nmask]/2 + (dof[nmask]/2-1)*np.log(dof[nmask]/2) - np.log(gamma(dof[nmask]/2))
      return p
   except:
      if p >0: return np.log(p)
      else: return -ts/2 + (dof/2-1)*np.log(dof/2) - np.log(gamma(dof/2))

   # underflow conditions - use asymptotic expansion
   return p

def ts2sig(ts,dof):

   if ts < 500:
      from uw.utilities.phasetools import sig2sigma
      from scipy.stats import chi2
      return sig2sigma(chi2.sf(ts,dof))
   
   from scipy.optimize import fsolve

   ts = float(ts)
   dof= float(dof)
   x0sq = ts - (dof - 2)*np.log(dof/2) + 2*np.log(gamma(dof/2))
   eta  = np.log(np.pi/2) - x0sq

   f    = lambda xsq: xsq + np.log(xsq) + eta

   return fsolve(f,x0sq)**0.5

# to do -- select peaks from the first grid scan using some "peak finding" algorithm
# based on prominence or what-have-you; say it should be maximum within ~delta 0.1 in phase
# and be > ~100 than the mean of the log likelihoods scan

def iter_peak_phase_scan(self,verbose=True,max_peaks=3):

   r   = self.ratios
   ph  = self.phases

   lct  = LCTemplate([]) #constant function
   lct.shift_mode = True #allow wrapping to negative positions
   lcf = SpectralLCFitter(lct,ph,r,1.,suppress_bg=False)

   dom = np.arange(0,1,self.scan_resolution)
   ts  = np.empty_like(dom)

   for i in xrange(max_peaks):

      num_peaks = len(lct.primitives) + 1

      norm = 0.5 #if num_peaks==1 else lct.norm()/num_peaks

      g1  = LCGaussian(p = [norm,self.scan_width,0])
      lct.primitives += [g1]

      if lct.norm () >= 1:
         n = lct.norm()
         for p in lct.primitives:
            p.p[0] /= n


      # do a grid search first to find best peak candidate
      for i,p in enumerate(dom):
         lct.primitives[-1].set_location(p)
         ts[i] = lcf.loglikelihood(lct.get_parameters(),lct)

      import pylab as P
      P.plot(dom,ts)

      sig = ts.max() - ts.mean()
      if verbose: print 'Delta Log L from scan: %.2f'%(sig)

      g1.set_location(dom[np.argmin(ts)])
      if verbose: print 'Peak trial position found at %.2f'%(g1.get_location())

      g1.free[:] = False
      g1.free[0] = True # fit norm only
      # tune new peak with others fixed    
      lcf.fit()
      if verbose:
         print '\nFinishing tuning peak:'
         print g1

      # code to check whether to break/remove last peak
      
   if verbose: print 'Performing final joint fit...'
   if lct.norm () >= 1:
      n = lct.norm()
      for p in lct.primitives:
         p.p[0] /= n

   for p in lct.primitives: p.free[:] = True
   lcf.fit()
   if verbose: print lcf
   return lcf

def twop_scan(self,verbose=True,max_peaks=1,free_norm=False):

   r   = self.ratios
   ph  = self.phases

   lct  = LCTemplate([]) #constant function
   lct.shift_mode = True #allow wrapping to negative positions
   lcf = SpectralLCFitter(lct,ph,r,1.,suppress_bg=False)

   dom = np.arange(0,1,self.scan_resolution)
   ts  = np.empty_like(dom)

   for i in xrange(max_peaks):

      num_peaks = len(lct.primitives) + 1

      norm = np.arcsin(2*(1-lct.norm()) - 1)

      g1  = LCGaussian(p = [norm,self.scan_width,0])
      lct.primitives += [g1]


      # do a grid search first to find best peak candidate
      for i,p in enumerate(dom):
         lct.primitives[-1].set_location(p)
         ts[i] = lcf.loglikelihood(lct.get_parameters(),lct)

      #import pylab as P
      #P.plot(dom,ts)

      sig = ts.max() - ts.mean()
      if verbose: print 'Delta Log L from scan: %.2f'%(sig)

      g1.set_location(dom[np.argmin(ts)])
      if verbose: print 'Peak trial position found at %.2f'%(g1.get_location())

     
   if verbose: print 'Performing final joint fit...'

   lct.primitives[0].free[0] = free_norm
   lct.primitives[0].free[1:] = True
   lcf.fit()
   if verbose: print lcf
   return lcf

def iter_peak_phase_scan2(self,verbose=True,max_peaks=3):

   r   = self.ratios
   ph  = self.phases

   lct  = LCTemplate([]) #constant function
   lct.shift_mode = True #allow wrapping to negative positions
   lcf = SpectralLCFitter(lct,ph,r,1.,suppress_bg=False)

   #dom = np.arange(0,1,self.scan_resolution)
   domw = np.asarray([0.025,0.05,0.075,0.10,0.125,0.15])
   domp = np.arange(0,1,0.02)
   ts   = np.empty(len(domw)*len(domp))

   for i in xrange(max_peaks):

      #norm = np.arcsin(2*0.5 -1) if i == 0 else np.arcsin(2*lct.norm() - 1)
      #norm = np.arcsin(2*(1-lct.norm()) - 1)
      #norm = np.arcsin(2./(i+1) - 1)
      norm  = 1./(i+1)

      g1  = LCGaussian(p = [norm,self.scan_width,0],prms=self.prms)
      lct.primitives += [g1]
      for prim in lct.primitives:
         prim.p[0] = norm

      # do a grid search first to find best peak candidate
      n = len(domw)
      for i,p in enumerate(domp):
         for j,w in enumerate(domw):
            g1.set_location(p)
            g1.p[1] = w
            ts[i*n + j] = lcf.loglikelihood(lct.get_parameters(),lct)
            #print p,w,ts[i*n + j]

      #import pylab as P
      #P.plot(dom,ts)

      #sig = ts.max() - ts.mean()
      #if verbose: print 'Delta Log L from scan: %.2f'%(sig)

      a = np.argmin(ts)
      g1.set_location(domp[a / n])
      g1.p[1] = domw[a % n]
      if verbose:
         print 'Peak posit found at %.3f'%(g1.get_location())
         print 'Peak width found at %.3f'%(g1.p[1])

      #if max_peaks > 1:
      #   g1.free[:] = False
      #   g1.free[0] = True # fit norm only
      #   # tune new peak with others fixed    
      #   lcf.fit()
      #   if verbose:
      #      print '\nFinishing tuning peak:'
      #      print g1
      #   #g1.free[:] = False

      # code to check whether to break/remove last peak
      
   if verbose: print 'Performing final joint fit...'
   if verbose: print lct
   for p in lct.primitives:
      p.free[:] = True
      p.free[0] = False
   lcf.fit(estimate_errors=False)
   #if max_peaks == 1:
   #   p.free[0] = False
   #   lcf.fit()
   #   p.free[0] = True
   #lcf.fit()
   if verbose: print lcf
   return lcf

def recalc_signals(self):
    edict = self.edict
    en  = edict['ENERGY']
    ct  = edict['CONVERSION_TYPE']
    psm = self.roi.psm
    psf = self.psf
    sd  = psm.point_sources[self.source_number].skydir

    dirs  = map(SkyDir,edict['RA'],edict['DEC'])
    diffs = np.asarray([sd.difference(d) for d in dirs])

    psf_vals = np.asarray([psf(e,c,diff,density=True)[0] for e,c,diff in zip(en,ct,diffs)])
    self.signals = psm.models[self.source_number](en) * psf_vals
    self.weights,self.ratios = self.signals / (self.signals + self.backs), self.signals / self.backs