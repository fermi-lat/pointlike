"""A module for spectral fitting based on the (spatial) likelihood provided by pointlike.
   The code is structured as wrappers for the C++ pointlike classes, with functionality
   pushed as close to the "band" level as possible.  E.g., each SimpleLikelihoodWrapper
   represents a single energy band and event class and knows how to calculate its likelihood,
   signal, flux, expected counts, background, etc.  PointSourceLikelihoodWrapper makes the
   transition from the band level to energy space.
   
   $Header$
   """

import numpy as N

class SimpleLikelihoodWrapper(object):
   """Wrap up SimpleLikelihood as the basic object for spectral fitting."""

   def init(self):
      """Establish defaults as class members."""
      self.n_simps = 8 #sampling points in log space for Simpson's rule
      self.emulate_unbinned = True
      self.extended_likelihood = False

   def __init__(self,sl,exposure,sky_dir,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
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

   def __call__(self,sky_dir,model=None):
      """Emulate PSF SkyFunction behavior, but with ability to use spectral model."""
      if model is None: return self.sl(sky_dir)
      expected = self.expected(model)
      s,g = self.sl.sigma(),self.sl.gamma()
      u = 0.5*(sky_dir.difference(self.sky_dir)/self.sl.sigma())**2
      return expected * (1-1./g)*(1+u/g)**-g / (2*3.1416*s**2)/self.psf_correction

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
      integral = (self.weights*self.point_likes*poiss_likes).sum()
      if integral <= 0.: return 1e6
      return -N.log(integral)

   def signal(self):
      """Model-independent estimation of the signal and error in this band."""
      if self.sl.photons()==0: return (0,3) #Feldman-Cousins upper limit (put in actual value!)
      try:
         from scipy.optimize import fmin
         delta,l = 0.01,self.logLikelihood
         signal  = fmin(l,max(1,self.sl.signal()),disp=0)[0]
         if signal < 0: raise Exception
         signal_err = ((l( signal*(1+delta) ) + l( signal*(1-delta) ) - 2*l(signal))/(delta*signal)**2)**-0.5
      except: #revert to cruder estimate
         print 'Warning: signal estimation failed at energy %d MeV'%(int(self.energy()))
         signal = self.sl.signal()
         signal_err = signal*(1./self.sl.photons() + (self.sl.sigma_alpha()/self.sl.alpha())**2)**0.5
      return (signal,signal_err)

   def flux(self,e_weight=0,cgs=True):
      """Return the differential flux, multiplied by energy**e_weight."""
      
      units = (1.60218e-6)**(e_weight-1) if cgs else 1 #convert MeV to ergs
      multi = units*self.energy()**e_weight/self.avg_exposure()/(self.emax-self.emin)/self.psf_correction
      signal,error = self.signal()
      return (multi*signal,multi*error)

   def energy(self):
      """Return energy 'associated' with this band, taken to mean where the flux should optimally be evaluated."""
      #return (self.sl.band().emax() + self.sl.band().emin())/2.
      return (self.emin*self.emax)**0.5

class PointSourceLikelihoodWrapper(object):

   def init(self):
      self.emin = 100
      self.emax = None
      self.quiet = False
      self.count = 0 #for iter interface

   def __init__(self,psl,exposure,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.psl,self.exposure = psl,exposure
      self.update()
      self.bin_centers = N.sort(list(set([x.energy() for x in self])))
      ens = N.asarray([ [x.emin,x.emax] for x in self]).flatten()
      self.emin,self.emax = ens.min(),ens.max() #replace defaults with observed

   def __iter__(self): return self

   def next(self):
      if self.count >= len(self.sl_wrappers): self.count = 0; raise StopIteration
      else: self.count+=1; return self.sl_wrappers[self.count-1]

   def update(self):
      self.psl.maximize()
      check = BandCheck(self.emin,self.emax)
      self.sl_wrappers = [SimpleLikelihoodWrapper(sl,self.exposure,self.psl.dir()) for sl in self.psl if check(sl)]
   
   def display(self,sky_dir,mode=4,model=None):
      """Return various spatial values at sky_dir(s).  All are photon densities.
         mode == 0: point source expectation
         mode == 1: observed values
         mode == 2: background prediction
         mode == 3: total prediction
         mode == 4: residuals"""
      options = ['sum( (slw(sky_dir,model) for slw in self) )',
                 'sum( (slw.sl.band()(sky_dir)/slw.sl.band().pixelArea() for slw in self) )',
                 'sum( (slw.sl.background_function()(sky_dir) for slw in self) )',
                 'self.display(sky_dir,mode=0,model=model) + self.display(sky_dir,mode=2,model=model)',
                 'self.display(sky_dir,mode=1,model=model) - self.display(sky_dir,mode=3,model=model)']
      exec 'result = %s'%options[mode] in locals()
      return result

   def spectrum(self,e_weight=0,cgs=True):
      """Return the differential flux as estimated from the band likelihoods; error estimated from curvature.
         The spectrum is averaged over front and back estimates, errors added in quadrature.
         
         Keyword arguments:
         e_weight    factor of energy by which to multiply dN/dE
         cgs         boolean; True to use ergs for energy, False for MeV"""
      keys,d = [str(e) for e in self.bin_centers.astype(int)],dict()
      for key in keys: d[key] = N.zeros(3).astype(float) #flux, error, num
      for slw in self:
         key = str(int(slw.energy()))
         f,e = slw.flux(e_weight=e_weight,cgs=cgs)
         d[key] += N.asarray([f,e**2,1.])
      results = N.asarray([ N.asarray([d[key][0],d[key][1]**0.5])/d[key][2] for key in keys])
      return (self.bin_centers,results[:,0],results[:,1])

   def counts(self,model=None,background=False):
      """
         background == True                   -> return background counts
         model == None && background == False -> return observed counts
         model != None && background == False -> return expected counts under model."""
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

class BandCheck(object):

   def __init__(self,emin,emax):
      self.femin,self.bemin,self.emax = emin,300,emax #temporary
   def __call__(self,sl):
      emin = self.femin if sl.band().event_class() == 0 else self.bemin
      return sl.band().emin() >= emin and (self.emax is None or sl.band().emax() <= self.emax)

class SpectralModelFitter(object):

   @staticmethod
   def least_squares(pslw,model,quiet=False):
      
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