"""
Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
   on the PSF.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_analysis.py,v 1.13 2010/05/18 22:25:36 kerrm Exp $

author: Matthew Kerr
"""

import numpy as N
import math
from roi_managers import *
from roi_bands import *
from roi_plotting import *
from roi_localize import *
from pointspec_helpers import PointSource

from specfitter import SpectralModelFitter,mycov

from collections import deque,defaultdict
from cPickle import dump

from scipy.optimize import fmin,fmin_powell,fmin_bfgs
from scipy.stats.distributions import chi2
from numpy.linalg import inv

EULER_CONST  = N.exp(1)
LOG_JACOBIAN = 1./N.log10(EULER_CONST)

###====================================================================================================###

class ROIAnalysis(object):

   def init(self):

      self.fit_emin = [125,125] #independent energy ranges for front and back
      self.fit_emax = [1e5 + 100,1e5 + 100] #0th position for event class 0
      self.quiet   = False
      self.verbose = False

      self.catalog_aperture = -1 # pulsar catalog analysis only -- deprecate
      self.phase_factor = 1.

   def __init__(self,roi_dir,ps_manager,ds_manager,spectral_analysis,**kwargs):
      """ roi_dir    -- the center of the ROI
          ps_manager -- an instance of ROIPointSourceManager
          ds_manager -- an instance of ROIDiffuseManager

          Optional Keyword Arguments
          ==========================
          fit_emin -- a two-element list giving front/back minimum energy
          fit_emax -- a two-element list giving front/back maximum energy
          phase_factor -- a correction that can be made if the data has
                          undergone a phase selection -- between 0 and 1
      """
      self.init()
      self.__dict__.update(**kwargs)
      self.roi_dir = roi_dir
      self.psm  = ps_manager
      self.bgm  = self.dsm = ds_manager

      self.sa   = spectral_analysis
      self.logl = None
      self.prev_logl = None
      self.__setup_bands__()
      self.bin_centers = N.sort(list(set([b.e for b in self.bands])))
      self.bin_edges   = N.sort(list(set([b.emin for b in self.bands] + [b.emax for b in self.bands])))

      self.param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
      self.param_vals  = N.concatenate([m.p for m in self.psm.models] + [m.p for m in self.bgm.models])
      self.psm.cache(self.bands)
      self.psm.update_counts(self.bands)

      #### TEMPORARY ###
      from roi_managers import ROIDiffuseManager
      if isinstance(ds_manager,ROIDiffuseManager):
         self.gradient = self.gradient_new
      else:
         self.gradient = self.gradient_old

   def __setup_bands__(self):

      self.bands = deque()
      for band in self.sa.pixeldata.dmap:
         evcl = band.event_class() & 1 # protect high bits

         if band.emin() >= self.fit_emin[evcl] and band.emax() < self.fit_emax[evcl]:
            self.bands.append(ROIBand(band,self.sa,self.roi_dir,catalog_aperture=self.catalog_aperture))

      self.bands = N.asarray(self.bands)

      self.psm.setup_initial_counts(self.bands)
      self.bgm.setup_initial_counts(self.bands)

      for band in self.bands: band.phase_factor = self.phase_factor

   def setup_energy_bands(self):

      groupings = dict()
      for bc in self.bin_centers:
         groupings[bc] = [band for band in self.bands if band.e==bc]

      self.energy_bands = [ROIEnergyBand(groupings[bc]) for bc in self.bin_centers]


   def logLikelihood(self,parameters,*args):

      bands = self.bands
      pf    = self.phase_factor

      if N.any(N.isnan(parameters)):
         # pretty ridiculous that this check must be made, but fitter passes NaNs...
         return 1e6
         # not sure if should "set parameters" in this case

      self.set_parameters(parameters)
      self.bgm.update_counts(bands)
      self.psm.update_counts(bands)

      ll = 0

      for b in bands:

         ll +=  (
                   #integral terms for ROI (go in positive)
                   (b.bg_all_counts + b.ps_all_counts)*pf

                   -

                   #pixelized terms (go in negative) -- note, no need for phase factor here
                   (b.pix_counts *
                       N.log(  b.bg_all_pix_counts + b.ps_all_pix_counts )
                   ).sum() if b.has_pixels else 0.
                )

      #print ll,parameters
      return 1e6 if N.isnan(ll) else ll

   def gradient_old(self,parameters,*args):
      """ Implement the gradient of the log likelihood wrt the model parameters."""

      bands    = self.bands
      pf       = self.phase_factor  # can just multiply the aperture term
      models   = self.psm.models
      bgmodels = self.bgm.models

      # sanity check -- for efficiency, the gradient should be called with the same params as the log likelihood
      if not N.allclose(parameters,self.parameters(),rtol=0):
         self.set_parameters(parameters)
         self.bgm.update_counts(bands)
         self.psm.update_counts(bands)

      gradient = N.zeros_like(parameters)

      # say byebye to abstraction and hello to disgusting code!

      # the positions in the list of point sources with free parameters
      indices = N.arange(len(models))[N.asarray([N.any(m.free) for m in models])]
      nparams = N.asarray([model.free.sum() for model in models])

      for b in bands:
         cp = 0
         if b.has_pixels:
            #tot_pix = b.bg_all_pix_counts + b.ps_all_pix_counts
            pix_weights = b.pix_counts / (b.bg_all_pix_counts + b.ps_all_pix_counts)
         else:
            pixterm = 0

         # do the backgrounds -- essentially, integrate the spectral gradient for each data pixel
         # and for the aperture
         for ind,model in enumerate(bgmodels):
            if not N.any(model.free):
               continue
            pts = model.gradient(b.bg_points)
            if len(pts.shape) == 1: pts = N.asarray([pts])
            for j in xrange(len(model.p)):
               if not model.free[j]: continue
               apterm = pf*(b.ap_evals[ind,:] * pts[j,:]).sum()
               if b.has_pixels:
                  pixterm = (pix_weights*(b.pi_evals[:,ind,:] * pts[j,:]).sum(axis=1)).sum()
               gradient[cp] += apterm - pixterm
               cp += 1

         # do the point sources -- call the band method, which will integrate the gradient
         # over the exposure -- NB -- inconsistency in that no correction is being made for
         # the position of the point sources...
         for ind,model in zip(indices,models[indices]):
            grad   = b.gradient(model)[model.free]*b.er[ind] # correct for exposure
            np     = nparams[ind]
            apterm = pf*b.overlaps[ind]
            if b.has_pixels:
               pixterm = (pix_weights*b.ps_pix_counts[:,ind]).sum()
            gradient[cp:cp+np] += grad * (apterm - pixterm)
            cp += np

      gradient *= 10**parameters * 1./N.log10(N.exp(1.)) # Jacobian
      return gradient

   def gradient_new(self,parameters,*args):
        """ Implement the gradient of the log likelihood wrt the model parameters."""

        bands    = self.bands
        pf       = self.phase_factor  # can just multiply the aperture term
        models   = self.psm.models

        # sanity check -- for efficiency, the gradient should be called with the same params as the log likelihood
        if not N.allclose(parameters,self.parameters(),rtol=0):
            self.set_parameters(parameters)
            self.bgm.update_counts(bands)
            self.psm.update_counts(bands)

        # do the point sources
        indices  = N.arange(len(models))[N.asarray([N.any(m.free) for m in models])]
        nparams  = N.asarray([model.free.sum() for model in models])
        gradient = N.zeros(nparams.sum())

        for b in bands:
            cp = 0
            if b.has_pixels:
                b.pix_weights = pix_weights = b.pix_counts / (b.bg_all_pix_counts + b.ps_all_pix_counts)
            else:
                pixterm = 0

            for ind,model in zip(indices,models[indices]):
                grad   = b.gradient(model)[model.free]*b.er[ind] # correct for exposure
                np     = nparams[ind]
                apterm = pf*b.overlaps[ind]
                if b.has_pixels:
                    pixterm = (pix_weights*b.ps_pix_counts[:,ind]).sum()
                gradient[cp:cp+np] += grad * (apterm - pixterm)
                cp += np

        # add in diffuse components
        gradient  = N.append(self.bgm.gradient(bands,pf),gradient)
        
        # transform into log space and return
        return gradient * 10**parameters * LOG_JACOBIAN
        

   def parameters(self):
      """Merge parameters from background and point sources."""
      return N.asarray(self.bgm.parameters()+self.psm.parameters())

   def get_parameters(self):
      """Support for hessian calculation in specfitter module."""
      return self.parameters()

   def get_free_errors(self):
      """Return the diagonal elements of the covariance matrix -- useful for step sizes in minimization, if known."""
      return N.asarray(self.bgm.get_free_errors() + self.psm.get_free_errors())

   def set_parameters(self,parameters):
      """Support for hessian calculation in specfitter module."""
      self.bgm.set_parameters(parameters,current_position=0)
      self.psm.set_parameters(parameters,current_position=len(self.bgm.parameters()))
      self.fit_parameters = parameters

   def fit_background(self):
      old_psm_frees = []
      for m in self.psm.models:
         old_psm_frees.append(m.free.copy())
         #m.free = N.asarray([False]*len(m.free))
         m.free[:] = False
      self.fit(fit_bg_first = False,estimate_errors=False)
      for n,nm in enumerate(self.psm.models):
         nm.free[:] = old_psm_frees[n]

   def __pre_fit__(self):

      #cache frozen values
      param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
      param_vals  = N.concatenate([m.p for m in self.psm.models] + [m.p for m in self.bgm.models])

      if len(param_state)  != len(self.param_state) or \
         N.any(param_state != self.param_state) or \
         N.any(param_vals  != self.param_vals):

         self.psm.cache(self.bands)

         ### NOTA BENE
         #self.bgm.cache() # remove if don't adopt with new paradigm
         #############

         self.param_state = param_state
         self.param_vals  = param_vals

   def fit(self,method='simplex', tolerance = 0.01, save_values = True, do_background=True,
                fit_bg_first = False, estimate_errors=True, error_for_steps=False,
                use_gradient = False, gtol = 1e-1):
      """Maximize likelihood and estimate errors.

         method    -- ['powell'] fitter; 'powell' or 'simplex'
         tolerance -- (approximate) absolute tolerance of log likelihood value
      """

      if fit_bg_first:
         self.fit_background()

      self.__pre_fit__()

      if not self.quiet: print '.....performing likelihood maximization...',
      if method == 'minuit':
         from uw.utilities.minuit import Minuit
         temp_params = self.parameters()
         npars = self.parameters().shape[0]
         param_names = ['p%i'%i for i in xrange(npars)]
         
         if use_gradient:
            gradient       = self.gradient
            force_gradient = 1
         else:
            gradient       = None
            force_gradient = 0

         if error_for_steps:
            steps = self.get_free_errors()
            steps[steps<1e-6] = 0.04 # for models without error estimates, put in the defaults
            steps[steps > 1]  = 1    # probably don't want to step more than 100%...
            m = Minuit(self.logLikelihood,temp_params,up=.5,maxcalls=20000,tolerance=tolerance,printMode=-self.quiet,param_names=param_names,steps=steps)
         else:
            m = Minuit(self.logLikelihood,temp_params,up=.5,maxcalls=20000,tolerance=tolerance,printMode=-self.quiet,param_names=param_names)

         params,fval = m.minimize()

         if save_values:
            if estimate_errors == True:
                self.__set_error_minuit(m,'HESSE')
            self.logLikelihood(params) # reset values to the ones found by minimization step
            self.prev_logl = self.logl if self.logl is not None else -fval
            self.logl = -fval
         #Saving this reference seems to cause a memory leak.
         #self._minuit = m
         return -fval
      else:
         ll_0 = self.logLikelihood(self.parameters())
         if use_gradient:
            f = self._save_bfgs = fmin_bfgs(self.logLikelihood,self.parameters(),self.gradient,full_output=1,maxiter=500,gtol=gtol,disp=0)
         else:
            minimizer  = fmin_powell if method == 'powell' else fmin
            f = minimizer(self.logLikelihood,self.parameters(),full_output=1,
                          maxiter=10000,maxfun=20000,ftol=0.01/abs(ll_0), disp=0 if self.quiet else 1)
         if not self.quiet: print 'Function value at minimum: %.8g'%f[1]
         if save_values:
            self.set_parameters(f[0])
            if estimate_errors: self.__set_error__(do_background,use_gradient)
            self.prev_logl = self.logl if self.logl is not None else -f[1]
            self.logl = -f[1]

      ## check for error conditions here
      #   if not self.quiet: print 'good fit!'
      #   return -f[1]

   def __set_error__(self,do_background=True,use_gradient=False):

      n = len(self.bgm.parameters())
      if use_gradient:
         hessian = mycov(self.gradient,self.parameters(),full_output=True)[1]
      else:
         hessian = SpectralModelFitter.hessian(self,self.logLikelihood)[0] #does Hessian for free parameters
      success = False
      # TODO -- check the return code

      try:
         if not do_background: raise Exception # what is this about?
         if not self.quiet: print 'Attempting to invert full hessian...'
         self.cov_matrix = cov_matrix = inv(hessian)
         if N.any(N.isnan(cov_matrix)):
            if not self.quiet: print 'Found NaN in covariance matrix!'
            raise Exception
         self.bgm.set_covariance_matrix(cov_matrix,current_position=0)
         self.psm.set_covariance_matrix(cov_matrix,current_position=n)
         success = True
      except:
         if len(self.psm.parameters()) > 0:
            if not self.quiet: print 'Skipping full Hessian inversion, trying point source parameter subset...'
            try:
               self.cov_matrix = cov_matrix = inv(hessian[n:,n:])
               if N.any(N.isnan(cov_matrix)):
                  if not self.quiet: print 'Found NaN in covariance matrix!'
                  raise Exception
               self.psm.set_covariance_matrix(cov_matrix,current_position=0)
               success = True
            except:
               if not self.quiet: print 'Error in calculating and inverting hessian.'
         else:
            np = len(self.get_parameters())
            self.cov_matrix = N.zeros([np,np])

      return success

   def __set_error_minuit(self,m,method='HESSE'):
      """Compute errors for minuit fit."""

      #Not sure yet if there will be problems with including the backgrounds.
      self.cov_matrix = m.errors(method=method)
      self.bgm.set_covariance_matrix(self.cov_matrix,current_position = 0)
      self.psm.set_covariance_matrix(self.cov_matrix,current_position = len(self.bgm.parameters()))

   def __str__(self):
      bg_header  = '======== BACKGROUND FITS =============='
      ps_header  = '======== POINT SOURCE FITS ============'
      if (self.logl is not None) and (self.prev_logl is not None):
         ll_string  = 'Log likelihood change: %.2f'%(self.logl - self.prev_logl)
      else:
         ll_string  = ''
      return '\n\n'.join([ps_header,self.psm.__str__(),bg_header,self.bgm.__str__(),ll_string])

   def TS(self,quick=True,which=0,method='simplex'):
      """Calculate the significance of the central point source.

         quick -- if set True, just calculate likelihood with source flux set to 0
                  if set False, do a full refit of all other free sources

         which -- the index of source to calculate -- default to central."""

      if quick:
         self.zero_ps(which)
         ll_0 = self.logLikelihood(self.get_parameters())
         self.unzero_ps(which)
         ll_1 = self.logLikelihood(self.get_parameters())
         return 2*(ll_0 - ll_1)

      save_params = self.parameters().copy() # save free parameters
      self.zero_ps(which)
      ll_0 = self.fit(save_values = False,method=method)
      print self
      self.unzero_ps(which)
      self.set_parameters(save_params) # reset free parameters
      self.__pre_fit__() # restore caching
      ll = -self.logLikelihood(save_params)
      return -2*(ll_0 - ll)

   def localize(self,which=0, tolerance=1e-3,update=False, verbose=False, bandfits=False, seedpos=None):
      """Localize a source using an elliptic approximation to the likelihood surface.

         which     -- index of point source; default to central
                      ***if localizing non-central, ensure ROI is large enough!***
         tolerance -- maximum difference in degrees between two successive best fit positions
         update    -- if True, update localization internally, i.e., recalculate point source contribution
         bandfits  -- if True, use a band-by-band (model independent) spectral fit; otherwise, use broabband fit
         seedpos   -- if set, use this position instead of the source position

         return fit position
      """
      rl = ROILocalizer(self,which=which,bandfits=bandfits,tolerance=tolerance,update=update,verbose=verbose)
      if seedpos is not None:
         rl.sd = seedpos  # override 
      return rl.localize()

   def upper_limit(self,which = 0,confidence = .95,e_weight = 0,cgs = False):
       """Compute an upper limit on the flux of a source.

          which      -- index of point source; default to central
          confidence -- confidence level for the upper limit, default to 95%

        The flux returned is an upper limit on the integral flux for the model
        above 100 MeV.

        N.B.-This still needs some work.  The fmin call will sometimes get lost
        and return an absurdly low limit (~10^-23 ph/cm^2/s, e.g.)
       """
       delta_logl = chi2.ppf(2*confidence-1,1)/2.
       params = self.parameters().copy()
       self.psm.models[which].p[0]  = -20
       zp = self.logLikelihood(self.parameters())

       def f(norm):
           self.psm.models[which].p[0] = N.log10(norm)
           ll = self.logLikelihood(self.parameters())
           return abs(ll - zp - delta_logl)

       limit = fmin(f,N.array([10**-6]),disp=0)[0]
       self.psm.models[which].p[0] = N.log10(limit)
       uflux = self.psm.models[which].i_flux(e_weight = e_weight,cgs = cgs)
       self.set_parameters(params)
       return uflux

   def printSpectrum(self,sources=None):
      """Print total counts and estimated signal in each band for a list of sources.

      Sources can be specified as PointSource objects, source names, or integers
      to be interpreted as indices for the list of point sources in the roi. If
      only one source is desired, it needn't be specified as a list. If no sources
      are specified, all sources with free fit parameters will be used."""
      if sources is None:
         sources = [s for s in self.psm.point_sources if N.any(s.model.free)]
      elif type(sources) != type([]):
         sources = [sources]
      if sources == []: return # No point sources in ROI
      bad_sources = []
      for i,s in enumerate(sources):
         if type(s) == PointSource:
            if not s in self.psm.point_sources:
               print 'Source not found in source list:\n%s\n'%s
               bad_sources += [s]
         elif type(s) == int:
            try:
               sources[i] = self.psm.point_sources[s]
            except IndexError:
               print 'No source #%i. Only %i source(s) specified.'\
                     %(s,len(self.psm.point_sources))
               bad_sources += [s]
         elif type(s) == type(''):
            names = [ps.name for ps in self.psm.point_sources]
            try:
               sources[i] = self.psm.point_sources[names.index(s)]
            except ValueError:
               print 'No source named %s'%s
               bad_sources += [s]
         else:
            print 'Unrecognized source specification:', s
            bad_sources += [s]
      sources = set([s for s in sources if not s in bad_sources])
      indices = [list(self.psm.point_sources).index(s) for s in sources]
      self.setup_energy_bands()

      fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic']\
                +[' '*15+'Signal']*len(sources)
      outstring = 'Spectra of sources in ROI about %s at ra = %.2f, dec = %.2f\n'\
                    %(self.psm.point_sources[0].name, self.roi_dir.ra(), self.roi_dir.dec())
      outstring += ' '*54+'  '.join(['%21s'%s.name for s in sources])+'\n'
      outstring += '  '.join(fields)+'\n'
      print outstring
      for eb in self.energy_bands:
         print eb.spectralString(which=indices)


   def save_fit(self,outfile,additional_data=None):
      """Save the spectral models (and locations) for all point sources and diffuse models.

         This saves the need to refit.  A future iteration should actually save all of the
         pixel predictions to avoid lengthy recalculation, too.

         additional_data: an optional dictionary with keys to add to output; note that
                          all entries should be serializable!"""

      d = defaultdict(list)
      for ps in self.psm.point_sources:
         m = ps.model
         m.ra  = ps.skydir.ra()
         m.dec = ps.skydir.dec()
         m.source_name = ps.name
         d['point_sources'].append(m)
      for bg in self.bgm.models:
         d['backgrounds'].append(bg)
      try:
         d['localization'] = [self.ldir.ra(),self.ldir.dec(),self.lsigma,self.qform.par]
      except:
         print 'No localization to save.'
      if additional_data is not None:
         try:    d.update(additional_data)
         except: print 'Warning! Could not merge requested keys into output dictionary.'
      f = open(outfile,'w')
      dump(d,f)
      f.close()

   def __call__(self,v):

      pass #make this a TS map? negative -- spatialLikelihood does it, essentially

   def add_ps(self,ps):
      """Add a new PointSource object to the model."""
      self.psm.add_ps(ps,self.bands)

   def del_ps(self,which):
      """Remove the PointSource at position given by which from the model."""
      return self.psm.del_ps(which,self.bands)

   def zero_ps(self,which):
      """Set the flux of point source given by which to 0."""
      return self.psm.zero_ps(which,self.bands)

   def unzero_ps(self,which):
      """Restore a previously-zeroed flux."""
      self.psm.unzero_ps(which,self.bands)

   def modify_loc(self,skydir,which):
      """Move point source given by which to new location given by skydir."""
      rl = ROILocalizer(self,which=which,update=True)
      rl.spatialLikelihood(skydir)
      self.psm.point_sources[which].skydir = skydir
      
   def print_summary(self, sdir=None, galactic=False, maxdist=5, title=''):
       """ formatted table point sources positions and parameter in the ROI
       Parameters
       ----------
           sdir : SkyDir, optional, default None for center
               for center: default will be the first source
           galactic : bool, optional, default False
              set true for l,b
           maxdist : float, optional, default 5
              radius in degrees

       """
       if sdir is None: sdir = self.psm.point_sources[0].skydir
       print '\n\t Nearby sources within %.1f degrees %s' % (maxdist,title)
       colstring = 'name dist ra dec flux8 index'
       if galactic: colstring =colstring.replace('ra dec', 'l b')
       colnames = tuple(colstring.split())
       n = len(colnames)-1
       print ('%-20s'+n*'%10s')% colnames
       for ps in self.psm.point_sources:
           dist=math.degrees(sdir.difference(ps.skydir))
           if maxdist and dist>maxdist:  break
           loc = (ps.skydir.l(),ps.skydir.b()) if galactic else (ps.skydir.ra(),ps.skydir.dec())
           fmt = '%-20s'+(n-2)*'%10.3f'+' '+2*'%9.2f%1s'
           freeflag = [ '*' if f else ' ' for f in ps.model.free]
           values = ((ps.name, dist) +loc
                   +( ps.model.fast_iflux()/1e-8, freeflag[0], 10**ps.model.p[1], freeflag[1]))
           print fmt % values

