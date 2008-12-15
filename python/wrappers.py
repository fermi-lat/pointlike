"""A module for spectral fitting based on the (spatial) likelihood provided by pointlike.
   The code is structured as wrappers for the C++ pointlike classes, with functionality
   pushed as close to the "band" level as possible.  E.g., each SimpleLikelihoodWrapper
   represents a single energy band and event class and knows how to calculate its likelihood,
   signal, flux, expected counts, background, etc.  PointSourceLikelihoodWrapper makes the
   transition from the band level to energy space.
   
   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/wrappers.py,v 1.6 2008/11/21 17:58:05 kerrm Exp $

   author: Matthew Kerr
   """

import numpy as N

class Wrapper(object):
    """Parent class for a wrapper.  Inherited classes are supposed to get a list of kwargs that may override
       the defaults established by init.  Class-specified initialization may be done with the setup method.
    """

    def init(self): pass
    def setup(self): pass

    def __init__(self,*args,**kwargs):
        self.init()
        for key,value in kwargs.items():
            if key in self.__dict__: self.__dict__[key] = value
        self.args = args
        self.setup()

class ExposureWrapper(Wrapper):
    """A wrapper class to create and wrap an Exposure object.         

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  CALDB       [None] If not specified, will use environment variable
  irf         ['P6_v1_diff'] Used for effective area
  =========   =======================================================
    """

    def init(self):
        self.CALDB = None
        self.quiet = False
        self.irf   ='P6_v1_diff'

    def setup(self):
        
        from skymaps import Exposure, EffectiveArea
        if self.CALDB is not None: EffectiveArea.set_CALDB(self.CALDB)
        lt = self.args[0].lt
        inst = ['front', 'back']
        self.ea  = [EffectiveArea(self.irf+'_'+x) for x in inst]
        if not self.quiet: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        self.exposure = [Exposure(lt,ea) for ea in self.ea]

    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)

    def __call__(self, event_class): return self.exposure[event_class]

class BackgroundWrapper(Wrapper):
    """A wrapper class to create and wrap a Background object.         

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  isotropic   [1.5e-5,2.1] Sreekumar EGRET values for an isotropic power law diffuse
  galactic_scale [1.] Scale factor for galprop model
  =========   =======================================================
    """

    def init(self):
        self.isotropic      = [1.5e-5,2.1]
        self.galactic_scale = 1.

    def setup(self):
        from skymaps import Background, DiffuseFunction, CompositeSkySpectrum, IsotropicPowerLaw, HealpixDiffuseFunc
        diffusefile,exposure = self.args
        import pyfits
        q = pyfits.open(diffusefile, memmap=1)
        if q[1].name == 'SKYMAP': # first extension name: is it healpix?
            self.galactic_diffuse = HealpixDiffuseFunc(self.diffusefile)
        else:
            self.galactic_diffuse = DiffuseFunction(diffusefile)
        self.isotropic_diffuse = IsotropicPowerLaw(self.isotropic[0],self.isotropic[1])
        self.diffuse = CompositeSkySpectrum(self.galactic_diffuse,self.galactic_scale);
        self.diffuse.add(self.isotropic_diffuse,1.)
        self.background = Background(self.diffuse, exposure(0), exposure(1)) # array interface does not work

    def __call__(self): return self.background

class SimpleLikelihoodWrapper(Wrapper):
   """
        call signature::

  slw = SimpleLikelihoodWrapper(sl,sky_dir,exposure, **kwargs)

Create a wrapper for SimpleLikelihood and spectral fitting.  

    sl: an instance of SimpleLikelihood
    exposure: an instance of ExposureWrapper
    sky_dir: an instance of SkyDir giving the position of the source.

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  n_simps     [8] sampling points in log space for Simpson's rule                                                                 
  emulate_unbinned [True] if True, ignore bins with zero counts                                                                           
  extended_likelihood  [False] if True, use the *Python* version of extended likelihood                                                                                                       
  prior_cut   [1] estimated value of signal fraction beyond which to use a power law
              Bayesian prior on the signal fraction
  =========   =======================================================
   """

   def init(self):
      self.n_simps = 8
      self.emulate_unbinned = False
      self.extended_likelihood = False
      self.prior_cut = 1.00 

   def setup(self):
      sl,exposure,sky_dir = self.args
      self.sl,self.sky_dir,self.exposure = sl,sky_dir,exposure
      self.emin,self.emax = sl.band().emin(),sl.band().emax()
      
      #Pieces for calculating the expected counts under a given model
      self.sampling_points  = N.logspace(N.log10(self.emin),N.log10(self.emax),self.n_simps+1)
      self.exposure_points  = N.fromiter( (exposure.value(sky_dir,e,sl.band().event_class()) for e in self.sampling_points) , float)
      self.simpsons_weights = N.log(self.sampling_points[-1] / self.sampling_points[0]) * \
                              N.asarray([1.] + ([4.,2.]*(self.n_simps/2))[:-1] + [1.])/(3.*self.n_simps)
      self.psf_correction   = 1 - (1+sl.umax()/sl.gamma())**(1-sl.gamma()) #fraction of signal in ROI

      if not (self.sl.extended_likelihood() or self.extended_likelihood): self.marg_setup()

   def marg_setup(self):
      #Pieces for calculating marginalization integral; needed for speed!
      n_points = min(max((int(100*(self.sl.photons()/20.)**0.5) >> 1) << 1,100),1000) #force n_points to be even
      self.points = N.linspace(0.001,1,n_points+1)
      self.weights = N.asarray([1.] + ([4.,2.]*(n_points/2))[:-1] + [1.])/(3.*n_points)
      self.point_likes = N.exp( N.nan_to_num(N.fromiter( (-self.sl(a)+self.sl() for a in self.points) , float)) )
      self.prior = 1.
      if self.sl.photons() > 0: #assign Bayesian prior
         ahat = 1. - self.sl.background()/self.sl.photons()
         if ahat > self.prior_cut:
            m = max((2*ahat - .1)/(1.-ahat),15)
            self.prior = (m + 1)*self.points**m

   def __call__(self,sky_dir,model=None):
      """Emulate PSF SkyFunction behavior, but with ability to use spectral model."""
      if self.sl.photons()==0: return 0 #kluge
      expected = (self.expected(model)/self.sl.signal()/self.psf_correction if model is not None else 1.)/self.psf_correction
      return expected * self.sl(sky_dir)

   def expected(self,model):
      """Return expected counts in the ROI for a source model.      
         Integration is by Simpson's rule in log space."""
      return self.psf_correction*\
            (self.simpsons_weights*self.sampling_points*model(self.sampling_points)*self.exposure_points).sum()

   def avg_exposure(self):
      """Return the logarithmic average of the exposure."""
      return self.expected(model=lambda e: 1./e)/N.log(self.emax/self.emin)/self.psf_correction

   def logLikelihood(self,expected):
      """Return the (negative) log likelihood.  If extended likelihood is enabled, use that.
         Otherwise, return the alpha marginalization."""
      
      if self.sl.photons() == 0: return 0 if self.emulate_unbinned else expected

      if self.extended_likelihood: #this should only be used if SimpleLikelihood::extended_likelihood==False!
         alpha = float(expected / (expected + self.sl.background()))
         pointlike,background = self.sl(alpha),self.sl.background()
         poisson = expected + background - self.sl.photons()*N.log(expected + background)
         return pointlike + poisson

      if self.sl.extended_likelihood():
         return self.sl.logLikelihood(float(expected))
      
      from scipy.stats import poisson
      poiss_likes = N.nan_to_num(poisson.pmf(self.sl.photons(),expected/self.points))
      integral = (self.weights*self.point_likes*poiss_likes*self.prior).sum()
      if integral <= 0.: return 1e6
      return -N.log(integral)

   def signal(self):
      """Model-independent estimation of the signal and error in this band."""
      if self.sl.photons()==0: return (0,3) #Feldman-Cousins upper limit (put in actual value!)
      try:
         from scipy.optimize import fmin
         delta,l = 0.01,self.logLikelihood
         signal  = fmin(l,max(1,self.sl.signal()),disp=0)[0]
         if signal < 0:
            #use Feldman-Cousins upper limit?
            signal,signal_err = 0,3
         else:
            signal_err = ((l( signal*(1+delta) ) + l( signal*(1-delta) ) - 2*l(signal))/(delta*signal)**2)**-0.5
      except: #revert to cruder estimate
         print 'Warning: signal estimation failed at energy %d MeV, event class %d'%(int(self.energy()),self.sl.band().event_class())
         signal = self.sl.signal()
         signal_err = signal*(1./self.sl.photons() + (self.sl.sigma_alpha()/self.sl.alpha())**2)**0.5
      return (signal,signal_err)

   def flux(self,e_weight=0,cgs=True,mev=True):
      """Return the differential flux estimate for this band.         

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  e_weight    [2] factor of energy by which to multiply dN/dE
  cgs         [True] boolean; True to use ergs for abscissa energy; False uses ordinate units
  mev         [True] boolean; True to use MeV for ordinate units; False to use GeV
  =========   =======================================================
      """
      
      units = (1.60218e-6 if mev else 1.60218e-9)**(e_weight-1) if cgs else (1. if mev else 10**(-3*(e_weight-1)))
      multi = units*self.energy()**e_weight/self.avg_exposure()/(self.emax-self.emin)/self.psf_correction
      signal,error = self.signal()
      return (multi*signal,multi*error)

   def energy(self):
      """Return energy 'associated' with this band, taken to mean where the flux should optimally be evaluated."""
      #return (self.sl.band().emax() + self.sl.band().emin())/2.
      return (self.emin*self.emax)**0.5

class PointSourceLikelihoodWrapper(Wrapper):
   """
        call signature::

  pslw = PointSourceLikelihoodWrapper(psl,exposure,**kwargs)

Create a wrapper for PointSourceLikelihood and spectral fitting.  

    psl: an instance of PointSourceLikelihood to wrap
    exposure: an instance of ExposureWrapper

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  emin       [200] minimum energy bin to consider for spectral fitting                                                                 
  emax       [None] maximum energy for fitting; if None, no limit                                                                         
  quiet      [False] if True, suppress printed output                                                                                                       
  =========   =======================================================
   """

   def init(self):
      self.emin = 200
      self.emax = None
      self.quiet = False

   def setup(self):
      self.count = 0 #for iter interface
      self.psl,self.exposure = self.args
      self.update()
      self.bin_centers = N.sort(list(set([x.energy() for x in self])))
      ens = N.asarray([ [x.emin,x.emax] for x in self]).flatten()
      self.emin,self.emax = ens.min(),ens.max() #replace defaults with observed

   def update(self):
      self.psl.maximize()
      check = self.BandCheck(self.emin,self.emax)
      #for sl in self.psl:
      #   print 'Energy: %d Event_class: %d Accepted: %s'%(sl.band().emin(),sl.band().event_class(),'True' if check(sl) else 'False')

      self.sl_wrappers = [SimpleLikelihoodWrapper(sl,self.exposure,self.psl.dir()) for sl in self.psl if check(sl)]

   def __iter__(self): return self

   def next(self):
      if self.count >= len(self.sl_wrappers): self.count = 0; raise StopIteration
      else: self.count+=1; return self.sl_wrappers[self.count-1]

   def __len__(self): return len(self.sl_wrappers)
   
   def display(self,sky_dir,mode=4,model=None):
      """Return various spatial values at sky_dir(s).  All are photon densities.
         mode == 0: point source expectation
         mode == 1: observed values
         mode == 2: background prediction
         mode == 3: total prediction
         mode == 4: residuals
         
         If a model is provided, this model is used to calculate the normalization for
         the point source; otherwise, the spatial likelihood and observed counts are used."""

      options = ['sum( (slw(sky_dir,model) for slw in self) )',
                 'sum( (slw.sl.band()(sky_dir)/slw.sl.band().pixelArea() for slw in self) )',
                 'sum( (slw.sl.background_function()(sky_dir) for slw in self) )',
                 'self.display(sky_dir,mode=0,model=model) + self.display(sky_dir,mode=2,model=model)',
                 'self.display(sky_dir,mode=1,model=model) - self.display(sky_dir,mode=3,model=model)']
      exec 'result = %s'%options[mode] in locals()
      return result

   def spectrum(self,e_weight=0,cgs=True,mev=True):
      """Return the differential flux as estimated from the band likelihoods.
The signal and error are estimated from the maximum value and curvature of the band likelihood.
The spectrum is averaged over front and back estimates, and the errors added in quadrature.
         

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  e_weight    [2] factor of energy by which to multiply dN/dE
  cgs         [True] boolean; True to use ergs for abscissa energy; False uses ordinate units
  mev         [True] boolean; True to use MeV for ordinate units; False to use GeV
  =========   =======================================================
      """

      keys,d = [str(e) for e in self.bin_centers.astype(int)],dict()
      for key in keys: d[key] = N.zeros(3).astype(float) #flux, error, num
      for slw in self:
         key = str(int(slw.energy()))
         f,e = slw.flux(e_weight=e_weight,cgs=cgs,mev=mev)
         d[key] += N.asarray([f,e**2,1.])
      results = N.asarray([ N.asarray([d[key][0],d[key][1]**0.5])/d[key][2] for key in keys])
      return (self.bin_centers if mev else self.bin_centers/1000.,results[:,0],results[:,1])

   def counts(self,model=None,background=False):
      """Return the counts predicted by a spectral model for this band.      

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  model       [None] return observed counts if None; else, return model integrated over exposure
  background  [True] if True, return the background counts for this band
  =========   =======================================================
      """

      keys,d = [str(e) for e in self.bin_centers.astype(int)],dict()
      for key in keys: d[key] = 0.
      for slw in self:
         key = str(int(slw.energy()))
         d[key] += slw.sl.background() if background else \
                   (slw.sl.photons() if model is None else slw.expected(model))
      return (self.bin_centers,N.asarray([float(d[key]) for key in keys]))

   def TS(self):
      """Return TS summed over bands."""
      return sum( (slw.sl.TS() for slw in self) )

   def enable_extended_likelihood(self,val=True):
      for slw in self: slw.extended_likelihood = val

   def set_prior_cut(self,cut=0.75):
      """Enable the Bayesian prior by overriding the default cut of 1."""
      for slw in self: slw.prior_cut = cut; slw.marg_setup()

   class BandCheck(object):
      """One-off class to determine whether a band should be included based on energy cuts."""

      def __init__(self,emin,emax):
         self.femin,self.bemin,self.emax = emin,max(300,emin),emax #temporary
      def __call__(self,sl):
         emin = self.femin if sl.band().event_class() == 0 else self.bemin
         return round(sl.band().emin()) >= round(emin) and (self.emax is None or round(sl.band().emax()) <= round(self.emax))

