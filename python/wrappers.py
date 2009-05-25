"""A module for spectral fitting based on the (spatial) likelihood provided by pointlike.
   The code is structured as wrappers for the C++ pointlike classes, with functionality
   pushed as close to the "band" level as possible.  E.g., each SimpleLikelihoodWrapper
   represents a single energy band and event class and knows how to calculate its likelihood,
   signal, flux, expected counts, background, etc.  PointSourceLikelihoodWrapper makes the
   transition from the band level to energy space.
   
   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/wrappers.py,v 1.9 2009/05/25 17:18:46 kerrm Exp $

   author: Matthew Kerr
   """

import numpy as N

class Singleton(object):

   __instances = {}

   def __init__(self,constructor,*args,**kwargs):
      inst = Singleton._Singleton__instances
      key  = str(constructor)
      if key not in inst.keys():
         inst[key] = constructor(*args,**kwargs)

   def __call__(self,constructor):
      inst = Singleton._Singleton__instances
      key  = str(constructor)
      if key in inst: return inst[key]

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
  irf         ['P6_v3_diff'] Used for effective area
  =========   =======================================================
    """

    def init(self):
        self.CALDB = None
        self.quiet = False
        self.exp_irf   ='P6_v3_diff'

    def setup(self):
        
        from skymaps import Exposure, EffectiveArea
        if self.CALDB is not None: EffectiveArea.set_CALDB(self.CALDB)
        else:
            from os import environ
            EffectiveArea.set_CALDB(environ['CALDB'])
        lt = self.args[0].lt
        inst = ['front', 'back']
        self.ea  = [EffectiveArea(self.exp_irf+'_'+x) for x in inst]
        if not self.quiet: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        self.exposure = [Exposure(lt,ea) for ea in self.ea]

    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)

    def __call__(self, event_class = 0): return self.exposure[event_class]

class BackgroundWrapper(Wrapper):
    """A wrapper class to create and wrap a Background object.         

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  isotropic   [1.5e-5,2.1] Sreekumar EGRET values for an isotropic power law diffuse
  galactic_scale [1.] Scale factor for galprop model
  iso_file    [None] an optional MapCube representing the isotropic diffuse
  iso_scale   [1.] scale for the optional MapCube; useful for phase selection
  =========   =======================================================
    """

    def init(self):
        self.isotropic      = [1.5e-5,2.1]
        self.galactic_scale = 1.
        self.iso_file       = None
        self.iso_scale      = 1.

    def setup(self):
        print 'Opening backgrounds...'
        from skymaps import Background, DiffuseFunction, CompositeSkySpectrum, IsotropicPowerLaw, HealpixDiffuseFunc
        diffusefile,exposure = self.args
        import pyfits
        q = pyfits.open(diffusefile, memmap=1)
        if q[1].name == 'SKYMAP': # first extension name: is it healpix?
            self.unscaled_galactic_diffuse_singl = Singleton(HealpixDiffuseFunc,diffusefile,HealpixDiffuseFunc.DIFFDENSITY)
            self.unscaled_galactic_diffuse = self.unscaled_galactic_diffuse_singl(HealpixDiffuseFunc)
        else:
            self.unscaled_galactic_diffuse_singl = Singleton(DiffuseFunction,diffusefile)
            self.unscaled_galactic_diffuse = self.unscaled_galactic_diffuse_singl(DiffuseFunction)
        q.close()
        if self.galactic_scale == 1.:
            self.galactic_diffuse = self.unscaled_galactic_diffuse
        else:
            self.galactic_diffuse  = CompositeSkySpectrum(self.unscaled_galactic_diffuse,self.galactic_scale)
        if self.iso_file is None:
           self.isotropic_diffuse = IsotropicPowerLaw(self.isotropic[0],self.isotropic[1])
        else:
            self.unscaled_isotropic_diffuse = DiffuseFunction(self.iso_file)
            if self.iso_scale == 1.:
                self.isotropic_diffuse = self.unscaled_isotropic_diffuse
            else:
                self.isotropic_diffuse = CompositeSkySpectrum(self.unscaled_isotropic_diffuse,self.iso_scale)
        self.diffuse = CompositeSkySpectrum(self.galactic_diffuse,1.);
        self.diffuse.add(self.isotropic_diffuse,1.)
        self.background = Background(self.diffuse, exposure(0), exposure(1)) # array interface does not work
        print 'Finished with the backgrounds...'

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
      self.empty = False
      
      #Pieces for calculating the expected counts under a given model
      self.sampling_points = sampling_points = N.logspace(N.log10(self.emin),N.log10(self.emax),self.n_simps+1)
      exposure_points      = N.fromiter( (exposure.value(sky_dir,e,sl.band().event_class()) for e in sampling_points) , float)
      simpsons_weights     = N.log(sampling_points[-1] / sampling_points[0]) * \
                             N.asarray([1.] + ([4.,2.]*(self.n_simps/2))[:-1] + [1.])/(3.*self.n_simps)
      
      self.pre_product    = sampling_points*exposure_points*simpsons_weights
      self.psf_correction = 1 - (1+sl.umax()/sl.gamma())**(1-sl.gamma()) #fraction of signal in ROI

      if not (self.sl.extended_likelihood() or self.extended_likelihood): self.__marg_setup__()

   def __marg_setup__(self):
      #Pieces for calculating marginalization integral; needed for speed!
      n_points = min(max((int(100*(self.sl.photons()/20.)**0.5) >> 1) << 1,100),1000) #force n_points to be even
      self.points = N.logspace(-3,0,n_points+1) if self.sl.alpha() < 0.1 else N.linspace(1e-3,1,n_points+1)
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

   def expected(self,model,normalize=True):
      """Return expected counts in the ROI for a source model.      
         Integration is by Simpson's rule in log space."""
      return (self.psf_correction if normalize else 1.)*\
             (self.pre_product*model(self.sampling_points)).sum()

   def avg_exposure(self):
      """Return the logarithmic average of the exposure."""
      return self.expected(model=lambda e: 1./e,normalize=False)/N.log(self.emax/self.emin)

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
      from scipy.integrate import simps
      poiss_likes = N.nan_to_num(poisson.pmf(self.sl.photons(),expected/self.points))
      integral = simps(self.point_likes*poiss_likes*self.prior,x=self.points)
      return -N.log(integral) if integral > 0 else 1e6

   def signal_ll(self,log_expected): return self.logLikelihood(N.exp(log_expected))

   def __signal_error__(self,signal,delta=0.01):
      """Signal is the maximum likelihood estimator."""
      logsignal,l = N.log(signal),self.signal_ll
      hessian = (l( logsignal*(1+delta) ) + l( logsignal*(1-delta) ) - 2*l(logsignal))/(delta*logsignal)**2
      if hessian >= 0:
         return signal*(N.exp(hessian**-0.5)-N.exp(-hessian**-0.5))/2.
      else: return signal*0.99

   def TS(self): return self.sl.TS()

   def signal(self):
      
      if self.sl.photons()==0: return (0,3) #Feldman-Cousins upper limit (put in actual value!)

      try: #maximize likelihood for expected source counts; use log transformation
         from scipy.optimize import fmin
         delta,l = 0.01,self.signal_ll
         seed = N.log(self.sl.signal()) if self.sl.signal() > 0 else 0.
         signal  = N.exp(fmin(l,seed,disp=0)[0])         
         if (signal < 1e-2 ): return 0,3 #Feldman-Cousins for signal consistent with 0
         signal_err = self.__signal_error__(signal)

      except: #revert to estimate from shape analysis
         print 'Warning: signal estimation failed at energy %d MeV, event class %d'%(int(self.energy()),self.sl.band().event_class())
         signal = self.sl.signal()
         signal_err = ( (self.sl.photons()*self.sl.sigma_alpha())**2 + self.sl.alpha()**2*self.sl.photons() )**0.5
      
      return (signal/self.psf_correction,signal_err/self.psf_correction)
   
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
      multi = units*self.energy()**e_weight/self.avg_exposure()/(self.emax-self.emin)
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
      self.emax = 5e5
      self.quiet = False
      self.back_multi = 1.

   def setup(self):
      self.count = 0 #for iter interface
      self.psl,self.exposure = self.args
      self.update()
            

   def update(self):
      self.psl.maximize()
      check = self.BandCheck(self.emin,self.emax,self.back_multi)
      self.sl_wrappers = [SimpleLikelihoodWrapper(sl,self.exposure,self.psl.dir(),**self.__dict__) for sl in self.psl if check(sl)]
      
      self.bin_centers = N.sort(list(set([x.energy() for x in self])))
      self.bin_edges = N.sort(list(set([x.emin for x in self] + [x.emax for x in self])))

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
      return sum( (slw.TS() for slw in self) )

   def enable_extended_likelihood(self,val=True):
      for slw in self: slw.extended_likelihood = val

   def set_prior_cut(self,cut=0.75):
      """Enable the Bayesian prior by overriding the default cut of 1."""
      for slw in self: slw.prior_cut = cut; slw.__marg_setup__()

   class BandCheck(object):
      """One-off class to determine whether a band should be included based on energy cuts."""

      def __init__(self,emin,emax,back_multi=1.):
         self.femin,self.bemin,self.emax = emin,emin*back_multi,emax #temporary
      def __call__(self,sl):
         emin = self.femin if sl.band().event_class() == 0 else self.bemin
         return round(sl.band().emin()) >= round(emin) and round(sl.band().emax()) <= round(self.emax)

