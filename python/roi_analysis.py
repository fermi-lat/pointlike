#TODO
#     - spatial residuals
#     - (hard) energy dispersion
#     - PSF convolution of background

"""Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
   on the PSF.  The heavy lifting is done in pointlike, where the spatial distributions are calculated
   initially.  Fitting of spectral models is accomplished to reasonably accuracy by scaling these
   initial estimates band-by-band, based on the revised spectral models, with each iteration."""

import numpy as N
from roi_modules import *
from roi_plotting import *
from psmanager import *

 
###====================================================================================================###

class ROIAnalysis(object):

   def init(self):
      self.fit_emin=0
      self.fit_emax=1e6
      self.small_scale = None

      self.localization_dir = None
      self.localization_err = None

      self.spatial = False
      self.logl = None
  
   def __init__(self,ps_manager,bg_manager,spectral_analysis,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.ps_manager,self.bg_manager,self.spectral_analysis = ps_manager,bg_manager,spectral_analysis
      self.setup_fitter()

   def setup_fitter(self):
      print 'Making and maximizing primary PointSourceLikelihood...'
      sa,ps_manager = self.spectral_analysis,self.ps_manager
      f = sa.fitter(ps_manager.name(),ps_manager.ROI_dir())
      self.fitter = f

      from skymaps import WeightedSkyDirList
      print 'Getting pixels...'
      for slw in f.pslw:
         
         slw.radius_in_rad  = (2*slw.sl.umax())**0.5*slw.sl.sigma()
         slw.use_empty      = False #slw.npix < 500
         slw.wsdl           = WeightedSkyDirList(slw.sl.band(),ps_manager.ROI_dir(),slw.radius_in_rad,slw.use_empty) # get empties
         slw.pix_counts     = N.asarray([x.weight() for x in slw.wsdl]) #observed counts
         slw.npix           = slw.wsdl.total_pix() #slw.sl.band().total_pix(ps_manager.ROI_dir(),slw.radius_in_rad)
         slw.solid_angle    = slw.npix*slw.sl.band().pixelArea() #ragged edge
         slw.solid_angle_p  = 2*N.pi*(1-N.cos(slw.radius_in_rad)) #solid angle for a pure circle
         
         #Simpson's rule factors for quadrature of background
         #note Jacobean!
         slw.simps_pts      = N.logspace(N.log10(slw.emin),N.log10(slw.emax),5)
         slw.simps_vec      = (slw.simps_pts*N.asarray([1,4,2,4,1]))/(3.*4.)*N.log(slw.emax/slw.emin)

      self.bg_manager.setup_initial_counts(f.pslw)
      self.ps_manager.setup_initial_counts(f.pslw)

      """
      #Make total counts for spatial likelihood (background frozen)
      for slw in f.pslw:

         #sum up totals
         slw.total_bg = slw.gal_counts + slw.iso_counts
         if len(self.ps_manager.point_sources) > 1:
            total_ps = (slw.ps_counts[1:]*slw.overlaps[1:]).sum()
         else:
            total_ps = 0
         slw.total_bg += total_ps

         #sum up pixels
         slw.total_bg_pix = slw.gal_pix_counts + slw.iso_pix_counts
         if len(self.ps_manager.point_sources) > 1 and slw.has_pixels:
            total_ps_pix = (slw.ps_pix_counts[:,1:]*slw.ps_counts[1:]).sum(axis=1)
         else:
            total_ps_pix = 0
         slw.total_bg_pix += total_ps_pix
         
         if slw.has_pixels:
            slw.qs = (slw.ps_pix_counts[:,0]/slw.total_bg_pix)*(slw.total_bg/slw.overlaps[0])
         else:
            slw.qs = 0
      """
   """
   def logLikelihood(self,parameters,*args):

      f = self.fitter

      self.set_parameters(parameters)
      self.bg_manager.update_counts(f.pslw)
      self.ps_manager.update_counts(f.pslw)

      m = self.ps_manager.mask

      tot_term  =  sum( (
                   slw.gal_counts + slw.iso_counts + \
                   (slw.ps_counts[m]*slw.overlaps[m]).sum() + slw.frozen_total_counts\
                   for slw in f.pslw if
                   (slw.energy() > self.fit_emin and slw.energy() < self.fit_emax) ) )

      pix_term   =  sum( ( (slw.pix_counts * N.log(
                    slw.gal_pix_counts + slw.iso_pix_counts +\
                    (slw.ps_pix_counts.transpose()[m].transpose()*slw.ps_counts[m]).sum(axis=1) +\
                    slw.frozen_pix_counts)).sum()
                    for slw in f.pslw if (slw.has_pixels
                    and slw.energy() > self.fit_emin and slw.energy() < self.fit_emax) ) )     

      r = tot_term - pix_term
      return 1e6 if N.isnan(r) else r
   """

   def logLikelihood(self,parameters,*args):

      f = self.fitter
      m = self.ps_manager.mask

      self.set_parameters(parameters)
      self.bg_manager.update_counts(f.pslw)
      self.ps_manager.update_counts(f.pslw)

      r = 0

      for s in f.pslw:
         if s.energy() < self.fit_emin or s.energy() > self.fit_emax: continue

         r +=  ( 
                #integral terms for ROI (go in positive)
                s.gal_counts + 
                s.iso_counts + 
                (s.ps_counts[m]*s.overlaps[m]).sum() + 
                s.frozen_total_counts -

                #pixelized terms (go in negative)
                (s.pix_counts * N.log(
                    s.gal_pix_counts + 
                    s.iso_pix_counts +
                    (s.unfrozen_pix_counts*s.ps_counts[m]).sum(axis=1) +
                    s.frozen_pix_counts)
                ).sum() if s.has_pixels else 0.
               )

      """
      r = sum( (
                #integral terms for ROI (go in positive)
                s.gal_counts + 
                s.iso_counts + 
                (s.ps_counts[m]*s.overlaps[m]).sum() + 
                s.frozen_total_counts -

                #pixelized terms (go in negative)
                (s.pix_counts * N.log(
                    s.gal_pix_counts + 
                    s.iso_pix_counts +
                    (s.unfrozen_pix_counts*s.ps_counts[m]).sum(axis=1) +
                    s.frozen_pix_counts)
                ).sum() if s.has pixels else 0

                for s in f.pslw if
                  s.energy() >= self.fit_emin and
                  s.energy() <= self.fit_emax
               )
             )
      """
      return 1e6 if N.isnan(r) else r

   def spatialLikelihood(self,parameters,*args):

      f = self.fitter

      self.set_parameters(parameters)
      self.ps_manager.update_counts(f.pslw)
      self.bg_manager.update_counts(f.pslw)

      logl = sum( (
                    (slw.pix_counts * N.log(1 + (slw.ps_counts[0]/(slw.ps_counts[0]+slw.total_bg))*(slw.qs - 1.)  ) ).sum()
                  for slw in f.pslw if (slw.has_pixels and slw.energy() > self.fit_emin and slw.energy() < self.fit_emax) )
                )
      
      
      #logl = sum( (
      #              (slw.pix_counts * N.log( (slw.ps_counts[0]/(slw.ps_counts[0]+slw.total_bg))*slw.qs + slw.total_bg_pix/slw.total_bg ) ).sum()
      #            for slw in f.pslw if (slw.has_pixels and slw.energy() > self.fit_emin and slw.energy() < self.fit_emax) )
      #          )

      return 1e6 if N.isnan(logl) else -logl
      
   def bandLikelihood(self,parameters,*args):

      scale = parameters[0]
      if scale < 0: return 1e6
      slw   = args[0]
      which = args[1]

      old_counts = slw.ps_counts[which]

      slw.ps_counts[which] *= scale

      pix_term   =  (slw.pix_counts * N.log(
                    slw.gal_pix_counts + slw.iso_pix_counts +\
                    (slw.ps_pix_counts*slw.ps_counts).sum(axis=1))).sum()\
                    if slw.has_pixels else 0.
                    

      total_term = slw.gal_counts + slw.iso_counts + \
                   (slw.ps_counts*slw.overlaps).sum()

      slw.ps_counts[which] = old_counts

      return total_term - pix_term
     

   def parameters(self):
      """Merge parameters from background and point sources."""
      return N.asarray(self.bg_manager.parameters()+self.ps_manager.parameters())

   def get_parameters(self):
      """Support for hessian calculation in specfitter module."""
      return self.parameters()

   def set_parameters(self,parameters):
      """Support for hessian calculation in specfitter module."""
      self.bg_manager.set_parameters(parameters,current_position=0)
      self.ps_manager.set_parameters(parameters,current_position=len(self.bg_manager.parameters()))
      self.fit_parameters = parameters
   
   def fit(self,method='simplex', tolerance = 1e-8, save_values = True, do_background=True, spatial=False):
      """Maximize likelihood and estimate errors.

         method -- ['powell'] fitter; 'powell' or 'simplex'
      """
      #cache frozen values
      self.ps_manager.cache(self.fitter.pslw)

      print 'Performing likelihood maximization...'
      from scipy.optimize import fmin,fmin_powell
      minimizer  = fmin_powell if method == 'powell' else fmin
      likelihood = self.spatialLikelihood if spatial else self.logLikelihood 
      f = minimizer(likelihood,self.parameters(),full_output=1,maxiter=10000,maxfun=20000,ftol=tolerance,xtol=tolerance)
      print 'Function value at minimum: %.8g'%f[1]
      if save_values:
         self.set_parameters(f[0])
         self.__set_error__(do_background)
         self.logl = -f[1]
      return -f[1] 

   def __set_error__(self,do_background=True, spatial=False):
      from specfitter import SpectralModelFitter
      from numpy.linalg import inv
      n = len(self.bg_manager.parameters())
      hessian = SpectralModelFitter.hessian(self,self.logLikelihood if not spatial else self.spatialLikelihood) #does Hessian for free parameters

      try:
         if not do_background: raise Exception
         print 'Attempting to invert full hessian...'
         cov_matrix = inv(hessian)
         self.bg_manager.set_covariance_matrix(cov_matrix,current_position=0)
         self.ps_manager.set_covariance_matrix(cov_matrix,current_position=n)
      except:
         print 'Skipping full Hessian inversion, trying point source parameter subset...'
         try:
            cov_matrix = inv(hessian[n:,n:])
            self.ps_manager.set_covariance_matrix(cov_matrix,current_position=0)
         except:
            print 'Error in calculating and inverting hessian.'

   def __str__(self):
      bg_header  = '======== BACKGROUND FITS =============='
      ps_header  = '======== POINT SOURCE FITS ============'
      return '\n\n'.join([ps_header,self.ps_manager.__str__(),bg_header,self.bg_manager.__str__()])
         
   def TS(self):
      """Calculate the significance of the central point source."""

      save_params = self.parameters().copy() #save parameters
      m = self.ps_manager.models[0]
      m.p[0] = -200 #effectively 0 flux
      save_free = N.asarray(m.free).copy()
      self.logLikelihood(self.parameters()) #update counts before freezing
      for i in xrange(len(m.free)): m.free[i] = False #freeze all parameters
      alt_ll = self.fit(save_values = False)
      for i in xrange(len(m.free)): m.free[i] = save_free[i] #unfreeze appropriate
      self.ps_manager.cache(self.fitter.pslw)
      ll = -self.logLikelihood(save_params) #reset predicted counts
      return -2*(alt_ll - ll)

   def localize(self):
      self.localization_err = self.spectral_analysis.fitter.psl.localize()
      self.localization_dir = self.spectral_analysis.fitter.psl.dir()


   def __call__(self,v):
      """Implement SkyFunction."""

      #PS manager part

      #BG manager part
      pass