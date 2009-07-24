"""
Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
   on the PSF.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/roi_analysis.py,v 1.10 2009/07/22 03:17:22 burnett Exp $

author: Matthew Kerr
"""

import numpy as N
from roi_modules import *
from roi_plotting import *
from psmanager import *

 
###====================================================================================================###

class ROIAnalysis(object):

   def init(self):

      self.fit_emin = [100,100] #independent energy ranges for front and back
      self.fit_emax = [5e5,5e5] #0th position for event class 0
      self.threading = False
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
      self.setup_bands()

      self.bin_centers = N.sort(list(set([b.e for b in self.bands])))
      self.bin_edges   = N.sort(list(set([b.emin for b in self.bands] + [b.emax for b in self.bands])))

      self.param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
      self.psm.cache(self.bands)
      self.psm.update_counts(self.bands)

   def setup_bands(self):
      
      from collections import deque
      self.bands = deque()
      for band in self.sa.pixeldata.dmap:

         if band.emin() >= self.fit_emin[band.event_class()] and band.emax() < self.fit_emax[band.event_class()]:
            self.bands.append(ROIBand(band,self.sa,catalog_aperture=self.catalog_aperture))

      self.bands = N.asarray(self.bands)

      self.psm.setup_initial_counts(self.bands)
      self.bgm.setup_initial_counts(self.bands)


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

      
   def bandLikelihood(self,parameters,*args):

      #it would probably be better to parameterize this with an actual flux...
      # (or, at least, a log scale! haha...)
      scale = parameters[0]
      if scale < 0: return 1e6
      band  = args[0]
      which = args[1]

      psc = band.ps_counts[which]

      tot_term = band.bg_all_counts + band.ps_all_counts + (psc*band.overlaps[which])*(scale-1)

      pix_term = (band.pix_counts * N.log(
                  band.bg_all_pix_counts + band.ps_all_pix_counts + (psc*band.ps_pix_counts[:,which]*(scale-1))
                 )).sum() if band.has_pixels else 0.

      return tot_term - pix_term
   
   # note -- not adapted for phase factor yet!
   def spatialLikelihood(self,skydir,which=0):
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
         exposure_ratio        = exp.value(skydir,en)/exp.value(rd,en)
         psf                   = PsfSkyFunction(skydir,gamma,sigma)
         
         overlap               = ro(band,rd,skydir) * exposure_ratio * band.solid_angle / band.solid_angle_p
         ool                   = band.overlaps[which]
         psc                   = band.ps_counts[which]

         tot_term              = band.bg_all_counts + band.ps_all_counts + psc * (overlap - ool)

         if band.has_pixels:

            ps_pix_counts = N.asarray(psf.wsdl_vector_value(band.wsdl))*((pa/(2*N.pi*sigma**2))*exposure_ratio)

            pix_term = (band.pix_counts * N.log
                           (
                           band.bg_all_pix_counts + 
                           band.ps_all_pix_counts + 
                           psc*(ps_pix_counts - band.ps_pix_counts[:,which])
                           )
                       ).sum()

         else:
            pix_term = 0

         ll += tot_term - pix_term

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
   
   def fit(self,method='simplex', tolerance = 1e-8, save_values = True, do_background=True):
      """Maximize likelihood and estimate errors.

         method -- ['powell'] fitter; 'powell' or 'simplex'
      """
      #cache frozen values
      param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
      if not N.all(param_state == self.param_state):
         self.psm.cache(self.bands)
         self.bgm.cache()
         self.param_state = param_state

      print 'Performing likelihood maximization...'
      from scipy.optimize import fmin,fmin_powell
      minimizer  = fmin_powell if method == 'powell' else fmin
      f = minimizer(self.logLikelihood,self.parameters(),full_output=1,\
                    maxiter=10000,maxfun=20000,ftol=tolerance,xtol=tolerance)
      if not self.quiet: print 'Function value at minimum: %.8g'%f[1]
      if save_values:
         self.set_parameters(f[0])
         self.__set_error__(do_background)
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
      return '\n\n'.join([ps_header,self.psm.__str__(),bg_header,self.bgm.__str__()])
         
   def TS(self):
      """Calculate the significance of the central point source."""

      save_params = self.parameters().copy() #save parameters
      m = self.psm.models[0]
      m.p[0] = -200 #effectively 0 flux
      save_free = N.asarray(m.free).copy()
      self.logLikelihood(self.parameters()) #update counts before freezing
      for i in xrange(len(m.free)): m.free[i] = False #freeze all parameters
      alt_ll = self.fit(save_values = False)
      for i in xrange(len(m.free)): m.free[i] = save_free[i] #unfreeze appropriate
      self.psm.cache(self.bands)
      ll = -self.logLikelihood(save_params) #reset predicted counts
      return -2*(alt_ll - ll)

   def localize(self,tolerance=1e-4,update=False,which=0):
      """Localize a source using an elliptic approximation to the likelihood surface.

         tolerance -- maximum difference in degrees between two successive best fit positions
         update    -- if True, update localization internally, i.e., recalculate point source contribution
                      [NOT IMPLEMENTED]
         which     -- index of point source; default to central [NOT IMPLEMENTED]
                      ***if localizing non-central, ensure ROI is large enough!***
      """
      import quadform
      rl = ROILocalizer(self)
      l  = quadform.Localize(rl,verbose = False)
      ld = SkyDir(l.dir.ra(),l.dir.dec())
      for i in xrange(5):
         l.fit(update=True)
         diff = l.dir.difference(ld)*180/N.pi
         if not self.quiet: print 'Difference from previous fit: %.5f deg'%(diff)
         if diff < tolerance:
            if not self.quiet: print 'Converged!'
            break
         ld = SkyDir(l.dir.ra(),l.dir.dec())

      self.qform   = l
      self.ldir    = l.dir
      self.lsigma  = l.sigma
      self.rl      = rl

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
      spectra = dict.fromkeys(sources)
      for s in sources:
         spectra[s] = band_spectra(self,source=self.psm.point_sources.index(s))
      iso,gal,src,obs = counts(self)[1:5]
      outstring = 'Spectra of sources in ROI about %s at ra = %f, dec = %f\n'\
                    %(self.psm.point_sources[0].name, self.sa.roi_dir.ra(), self.sa.roi_dir.dec())
      fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic']\
                +[' '*10+'Signal']*len(sources)
      outstring += ' '*54+'  '.join(['%18s'%s.name for s in sources])+'\n'
      outstring += '  '.join(fields)+'\n'
      for i,band in enumerate(zip(self.bands[::2],self.bands[1::2])):
         values = (band[0].emin, band[0].radius_in_rad*180/N.pi,
                   band[1].radius_in_rad*180/N.pi,obs[i],gal[i],iso[i])
         for s in sources:
            values+=(spectra[s][1][i],.5*(spectra[s][3][i]-spectra[s][2][i]))          
         string = '  '.join(['%6i','%6.2f','%6.2f','%7i','%8i','%9i']+
                            ['%8.1f +/-%6.1f']*len(sources))%values
         outstring += string+'\n'
      print outstring


   def __call__(self,v):
      
      pass #make this a TS map? negative -- spatialLikelihood does it, essentially