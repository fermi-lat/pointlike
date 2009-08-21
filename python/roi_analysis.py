"""
Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
   on the PSF.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/roi_analysis.py,v 1.19 2009/08/20 22:37:15 wallacee Exp $

author: Matthew Kerr
"""

import numpy as N
from roi_modules import *
from roi_plotting import *
from psmanager import *

 

###====================================================================================================###

class ROIAnalysis(object):

   def init(self):

      self.fit_emin = [99,99] #independent energy ranges for front and back
      self.fit_emax = [1e5 + 100,1e5 + 100] #0th position for event class 0
      self.quiet = False

      self.catalog_aperture = -1

      self.phase_factor = 1.
  
   def __init__(self,ps_manager,bg_manager,spectral_analysis,**kwargs):
      self.init()
      self.__dict__.update(**kwargs)
      self.psm  = ps_manager
      self.bgm  = bg_manager

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

   def __setup_bands__(self):
      
      from collections import deque
      self.bands = deque()
      for band in self.sa.pixeldata.dmap:

         if band.emin() >= self.fit_emin[band.event_class()] and band.emax() < self.fit_emax[band.event_class()]:
            self.bands.append(ROIBand(band,self.sa,catalog_aperture=self.catalog_aperture))

      self.bands = N.asarray(self.bands)

      self.psm.setup_initial_counts(self.bands)
      self.bgm.setup_initial_counts(self.bands)

      for band in self.bands: band.phase_factor = self.phase_factor

   def setup_energy_bands(self):

      from collections import defaultdict

      groupings = defaultdict(list)
      for bc in self.bin_centers:
         for band in self.bands:
            if band.e == bc:
               groupings[bc].append(band)

      self.energy_bands = [ROIEnergyBand(groupings[bc]) for bc in self.bin_centers]

   def add_ps(self,ps):
      self.psm.add_ps(ps,self.bands)

   def logLikelihood(self,parameters,*args):

      bands = self.bands

      self.set_parameters(parameters)
      self.bgm.update_counts(bands)
      self.psm.update_counts(bands)

      ll = 0

      for b in bands:

         ll +=  ( 
                   #integral terms for ROI (go in positive)
                   (b.bg_all_counts + b.ps_all_counts)*self.phase_factor

                   -

                   #pixelized terms (go in negative)
                   (b.pix_counts *
                       N.log(  (b.bg_all_pix_counts + b.ps_all_pix_counts)*self.phase_factor )
                   ).sum() if b.has_pixels else 0.
                )

      return 1e6 if N.isnan(ll) else ll

   def spatialLikelihood(self,skydir,update=False,which=0):
      """Calculate log likelihood as a function of position a point source.
      
         which   -- index of point source; default to central
                    ***if localizing non-central, ensure ROI is large enough!***
      """
      
      ro = ROIOverlap()
      rd = self.sa.roi_dir
      ll = 0

      from skymaps import PsfSkyFunction

      for i,band in enumerate(self.bands):
   
         sigma,gamma,en,exp,pa = band.s,band.g,band.e,band.exp,band.b.pixelArea()
         exposure_ratio        = exp.value(skydir,en)/exp.value(rd,en) #needs fix? -- does not obey which?
         psf                   = PsfSkyFunction(skydir,gamma,sigma)
         
         overlap               = ro(band,rd,skydir) * exposure_ratio * band.solid_angle / band.solid_angle_p
         ool                   = band.overlaps[which]
         psc                   = band.ps_counts[which]

         tot_term              = (band.bg_all_counts + band.ps_all_counts + psc * (overlap - ool)) * self.phase_factor

         if band.has_pixels:

            ps_pix_counts = N.asarray(psf.wsdl_vector_value(band.wsdl))*((pa/(2*N.pi*sigma**2))*exposure_ratio)

            pix_term = (band.pix_counts * N.log
                           ( self.phase_factor * 
                              (
                              band.bg_all_pix_counts + 
                              band.ps_all_pix_counts + 
                              psc*(ps_pix_counts - band.ps_pix_counts[:,which])
                              )
                           )
                       ).sum()

         else:
            pix_term = 0

         ll += tot_term - pix_term
         if N.isnan(ll):
            raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f, band %d' %(skydir.ra(),skydir.dec(),i))

         if update:
            band.overlaps[which] = overlap 
            band.ps_all_counts  += psc * (overlap - ool)
            if band.has_pixels:
               band.ps_all_pix_counts      += psc * (ps_pix_counts - band.ps_pix_counts[:,which])
               band.ps_pix_counts[:,which]  = ps_pix_counts

      if N.isnan(ll):
         raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f' %(skydir.ra(),skydir.dec()))
      return ll


   def parameters(self):
      """Merge parameters from background and point sources."""
      return N.asarray(self.bgm.parameters()+self.psm.parameters())

   def get_parameters(self):
      """Support for hessian calculation in specfitter module."""
      return self.parameters()

   def set_parameters(self,parameters):
      """Support for hessian calculation in specfitter module."""
      self.bgm.set_parameters(parameters,current_position=0)
      self.psm.set_parameters(parameters,current_position=len(self.bgm.parameters()))
      self.fit_parameters = parameters

   def fit_background(self):
      old_psm_frees = []
      for m in self.psm.models:
         old_psm_frees.append(m.free.copy())
         m.free = N.asarray([False]*len(m.free))
      self.fit(fit_bg_first = False)
      for n,nm in enumerate(self.psm.models):
         nm.free = old_psm_frees[n]

   def __pre_fit__(self):
      
      #cache frozen values
      param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
      param_vals  = N.concatenate([m.p for m in self.psm.models] + [m.p for m in self.bgm.models])
      
      if len(param_state)  != len(self.param_state) or \
         N.any(param_state != self.param_state) or \
         N.any(param_vals  != self.param_vals):
         
         self.psm.cache(self.bands)
         self.bgm.cache()
         self.param_state = param_state
         self.param_vals  = param_vals
   
   def fit(self,method='simplex', tolerance = 1e-8, save_values = True, do_background=True, fit_bg_first = False):
      """Maximize likelihood and estimate errors.

         method -- ['powell'] fitter; 'powell' or 'simplex'
      """

      if fit_bg_first:
         self.fit_background()

      self.__pre_fit__()

      if not self.quiet: print 'Performing likelihood maximization...'
      from scipy.optimize import fmin,fmin_powell
      minimizer  = fmin_powell if method == 'powell' else fmin
      ll_0 = self.logLikelihood(self.parameters())
      f = minimizer(self.logLikelihood,self.parameters(),full_output=1,\
                    maxiter=10000,maxfun=20000,ftol=0.01/abs(ll_0) )
      if not self.quiet: print 'Function value at minimum: %.8g'%f[1]
      if save_values:
         self.set_parameters(f[0])
         self.__set_error__(do_background)
         self.prev_logl = self.logl if self.logl is not None else -f[1]
         self.logl = -f[1]
      return -f[1] 

   def __set_error__(self,do_background=True):
      from specfitter import SpectralModelFitter
      from numpy.linalg import inv
      n = len(self.bgm.parameters())
      hessian = SpectralModelFitter.hessian(self,self.logLikelihood) #does Hessian for free parameters

      try:
         if not do_background: raise Exception
         if not self.quiet: print 'Attempting to invert full hessian...'
         cov_matrix = inv(hessian)
         self.bgm.set_covariance_matrix(cov_matrix,current_position=0)
         self.psm.set_covariance_matrix(cov_matrix,current_position=n)
      except:
         if not self.quiet: print 'Skipping full Hessian inversion, trying point source parameter subset...'
         try:
            cov_matrix = inv(hessian[n:,n:])
            self.psm.set_covariance_matrix(cov_matrix,current_position=0)
         except:
            if not self.quiet: print 'Error in calculating and inverting hessian.'

   def __str__(self):
      bg_header  = '======== BACKGROUND FITS =============='
      ps_header  = '======== POINT SOURCE FITS ============'
      if (self.logl is not None) and (self.prev_logl is not None):
         ll_string  = 'Log likelihood change: %.2f'%(self.logl - self.prev_logl)
      else:
         ll_string  = ''
      return '\n\n'.join([ps_header,self.psm.__str__(),bg_header,self.bgm.__str__(),ll_string])
         
   def TS(self,quick=True,which=0):
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
      ll_0 = self.fit(save_values = False)
      print self
      self.unzero_ps(which)
      self.set_parameters(save_params) # reset free parameters
      self.__pre_fit__() # restore caching
      ll = -self.logLikelihood(save_params)
      return -2*(ll_0 - ll)

   def localize(self,which=0, tolerance=1e-3,update=False, verbose=False):
      """Localize a source using an elliptic approximation to the likelihood surface.

         which     -- index of point source; default to central 
                      ***if localizing non-central, ensure ROI is large enough!***
         tolerance -- maximum difference in degrees between two successive best fit positions
         update    -- if True, update localization internally, i.e., recalculate point source contribution

         return fit position
      """
      import quadform
      rl = ROILocalizer(self,which=which)
      l  = quadform.Localize(rl,verbose = verbose)
      ld = SkyDir(l.dir.ra(),l.dir.dec())
      ps = self.psm.point_sources[which]
      ll_0 = self.spatialLikelihood(self.psm.point_sources[which].skydir,update=False,which=which)

      if not self.quiet:
         fmt ='Localizing source %s, tolerance=%.1e...\n\t'+7*'%10s'
         tup = (ps.name, tolerance,)+tuple('moved delta ra    dec   a    b  qual'.split())
         print fmt % tup
         print ('\t'+4*'%10.4f')% (0,0,ps.skydir.ra(), ps.skydir.dec())
         diff = l.dir.difference(ps.skydir)*180/N.pi
         print ('\t'+7*'%10.4f')% (diff,diff, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
            
      for i in xrange(5):
         try:
            l.fit(update=True)
         except:
            l.recenter()
            if not self.quiet: print 'trying a recenter...'
            continue
         diff = l.dir.difference(ld)*180/N.pi
         delt = l.dir.difference(ps.skydir)*180/N.pi
         if not self.quiet: print ('\t'+7*'%10.4f')% (diff, delt, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
         if diff < tolerance:
            break
         ld = SkyDir(l.dir.ra(),l.dir.dec())

      if update:
         self.psm.point_sources[which].skydir = l.dir
         
      ll_1 = self.spatialLikelihood(l.dir,update=update,which=which)
      if not self.quiet: print 'TS change: %.2f'%(2*(ll_0 - ll_1))

      self.qform   = l
      self.ldir    = l.dir
      self.lsigma  = l.sigma
      self.rl      = rl
      self.delta_loc_logl = (ll_0 - ll_1)
      return l.dir

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
      indices = [self.psm.point_sources.index(s) for s in sources]
      self.setup_energy_bands()

      fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic']\
                +[' '*15+'Signal']*len(sources)
      outstring = 'Spectra of sources in ROI about %s at ra = %f, dec = %f\n'\
                    %(self.psm.point_sources[0].name, self.sa.roi_dir.ra(), self.sa.roi_dir.dec())
      outstring += ' '*54+'  '.join(['%21s'%s.name for s in sources])+'\n'
      outstring += '  '.join(fields)+'\n'
      print outstring
      for eb in self.energy_bands:
         print eb.spectralString(which=indices)


   def save_fit(self,outfile):
      """Save the spectral models (and locations) for all point sources and diffuse models.
      
         This saves the need to refit.  A future iteration should actually save all of the
         pixel predictions to avoid lengthy recalculation, too."""

      from cPickle import dump
      from collections import defaultdict
      d = defaultdict(list)
      for ps in self.psm.point_sources:
         m = ps.model
         m.ra  = ps.skydir.ra()
         m.dec = ps.skydir.dec()
         m.source_name = ps.name
         d['point_sources'].append(m)
      for bg in self.bgm.models:
         d['backgrounds'].append(bg)
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
      self.spatialLikelihood(skydir,update=True,which=which)