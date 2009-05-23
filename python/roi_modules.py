#TODO
#     - (less easy) plotting - residuals sortof implemented
#     - spatial residuals
#     - (hard) energy dispersion

"""Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
   on the PSF.  The heavy lifting is done in pointlike, where the spatial distributions are calculated
   initially.  Fitting of spectral models is accomplished to reasonably accuracy by scaling these
   initial estimates band-by-band, based on the revised spectral models, with each iteration."""

import numpy as N
from pointspec import SpectralAnalysis
from skymaps import SkyDir,WeightedSkyDirList,SkyIntegrator
from Models import *
from psf import PSF

class ROIOverlap(object):
   """Routines to calculate how much of the emission of a point source falls onto an ROI."""

   def init(self):
      self.quadrature_tol = 1e-5

   def __init__(self,**kwargs):
      self.init()
      self.__dict__.update(kwargs)


   def __call__(self,slw,roi_dir,ps_dir,max_radius = None,fixed_radius=None):
      """Return an array of fractional overlap for a point source at location skydir.
         Note radius arguments are in radians."""

      from skymaps import PsfSkyFunction
      from scipy.integrate import quad
      from math import cos,sin
      
      sigma,gamma,roi_u = slw.sl.sigma(),slw.sl.gamma(),slw.sl.umax()
      roi_rad           = sigma*(roi_u*2)**0.5
      
      if max_radius is not None:
         roi_rad = min(max_radius,roi_rad)
         roi_u   = 0.5*(roi_rad/sigma)**2

      if fixed_radius is not None:
         roi_rad = fixed_radius
         roi_u   = 0.5*(roi_rad/sigma)**2
      
      offset = roi_dir.difference(ps_dir)

      if offset < 1e-5:
         return 1-(1+roi_u/gamma)**(1-gamma) #point source in center of ROI, symmetry

      #Testing code -- compare numeric to analytic results below      
      #fpsf    = PsfSkyFunction(ps_dir,gamma,sigma)
      #avg     = fpsf.average(roi_dir,roi_rad,self.solid_angle_epsilon)
      #solid   = N.pi*roi_rad**2
      #jac     = 2*N.pi*sigma**2
      #overlap_num = avg*solid/jac

      if offset < roi_rad:

         def interior(x):

            c  = cos(x)
            eff_rad = ( roi_rad**2 + offset**2*(c**2 - 1) )**0.5 - offset*c
            u = 0.5 * (eff_rad / sigma)**2

            return 1.- (1. + u/gamma )**(1-gamma)

         overlap = quad(interior,0,N.pi,epsabs=self.quadrature_tol,epsrel=self.quadrature_tol)[0]/N.pi

      else:

         def exterior(x):

            c    = cos(x)
            s    = (1-c**2)**0.5

            r2 = ( roi_rad**2 - (offset*s)**2 )**0.5
            de = offset*c - r2

            u1 = 0.5 * (de/sigma)**2
            u2 = 0.5 * ((de + 2*r2)/sigma)**2

            return (1. + u1/gamma )**(1-gamma) - (1. + u2/gamma )**(1-gamma)

         limit = N.arcsin(roi_rad / offset)
         overlap = quad(exterior,0,limit,epsabs=self.quadrature_tol,epsrel=self.quadrature_tol)[0]/N.pi

      return overlap

      
###====================================================================================================###

class ROIModelManager(object):
   """Parent class for point source manager and background manager.  Provides universal
      methods to manager parameters and error matricies."""

   def parameters(self):
      #if len(self.models[self.free])==0: return []
      return list(N.concatenate([m.get_parameters() for m in self.models]))

   def set_parameters(self,parameters,current_position=0):
      for m in self.models:
         cp,np = current_position,current_position+len(m.get_parameters())
         m.set_parameters(parameters[cp:np])
         current_position += np-cp

   def set_covariance_matrix(self,cov_matrix,current_position=0):
      for m in self.models:
         cp,np = current_position,current_position+len(m.get_parameters())
         m.set_cov_matrix(cov_matrix[cp:np,cp:np])
         current_position += np-cp
      

###====================================================================================================###

class ROISourceManager(ROIModelManager):
   """Manage all point sources."""
   
   def __init__(self,point_sources):
      self.point_sources = point_sources
      self.models = N.asarray([point_source.model for point_source in point_sources])
      self.mask   = N.asarray([True]*len(self.models))

   def __len__(self): return len(self.point_sources)

   def __str__(self):
      return '\n\n'.join(['%s fitted with %s\n'%(ps.name,ps.model.pretty_name)+ps.model.__str__()\
                          for ps in self.point_sources])

   def ROI_dir(self): return self.point_sources[0].skydir

   def name(self): return self.point_sources[0].name

   def setup_initial_counts(self,pslw):

      print 'Setting up point sources...'

      #update the exposure to do event class properly (if gonna go, go all the way...)

      exp     = pslw.exposure().value
      roi_dir = self.ROI_dir()
      overlap = ROIOverlap()
      from skymaps import PsfFunction,PsfSkyFunction

      for i,slw in enumerate(pslw):

         sigma,gamma,band,en  = slw.sl.sigma(),slw.sl.gamma(),slw.sl.band(),slw.energy()
         exposure_ratios      = N.asarray([exp(ps.skydir,en)/exp(roi_dir,en) for ps in self.point_sources]) #make a first-order correction for exposure variation
         #psf                  = PsfFunction(gamma)

         slw.has_pixels       = slw.sl.photons() >0
         #slw.ps_pix_counts    = N.asarray([[psf(ps.skydir,d,sigma)*exposure_ratios[nps] for nps,ps in enumerate(self.point_sources) ] for d in slw.wsdl])\
         #                       * (band.pixelArea()/(2*N.pi*sigma**2))
         slw.ps_pix_counts    = N.empty([len(slw.wsdl),len(self.point_sources)])
         for nps,ps in enumerate(self.point_sources):
            psf = PsfSkyFunction(ps.skydir,gamma,sigma)
            slw.ps_pix_counts[:,nps] = psf.wsdl_vector_value(slw.wsdl)
         slw.ps_pix_counts   *= (band.pixelArea()/(2*N.pi*sigma**2)*exposure_ratios)
         #slw.ps_pix_counts    = slw.ps_pix_counts.transpose()

         if slw.use_empty:
            slw.overlaps      = slw.ps_pix_counts.sum(axis=0)
         else:
            slw.overlaps      = N.asarray([overlap(slw,roi_dir,ps.skydir) for ps in self.point_sources])
            slw.overlaps     *= (slw.solid_angle/slw.solid_angle_p) #small correction for ragged edge
         slw.overlaps        *= exposure_ratios
         slw.ps_counts        = N.asarray([slw.expected(model,normalize=False) for model in self.models])
         slw.src_scale        = 1.

         slw.frozen_pix_counts = 0
         slw.frozen_total_counts = 0
         slw.unfrozen_pix_counts = 0

      
   def update_counts(self,pslw):
      #updated to save computation for frozen models - is the cost prohibitive for full varying?
      for i,model in enumerate(self.models):
         if self.mask[i]:
            for slw in pslw: slw.ps_counts[i] = slw.expected(model,normalize=False)

   def cache(self,pslw):
      """Cache values for models with no degrees of freedom.  Helps a lot near galactic plane."""
      self.mask = m = N.asarray([N.any(model.free) for model in self.models])
      nm = N.logical_not(self.mask)
      s = nm.sum()
      for slw in pslw.sl_wrappers:

         if slw.has_pixels:
            slw.unfrozen_pix_counts = slw.ps_pix_counts.transpose()[m].transpose()
         else:
            slw.unfrozen_pix_counts = 0

         if s > 0:
            slw.frozen_total_counts = (slw.ps_counts[nm]*slw.overlaps[nm]).sum()
            if not slw.has_pixels:
               slw.frozen_pix_counts = 0
               continue
            if s > 1:
               slw.frozen_pix_counts   = (slw.ps_pix_counts.transpose()[nm].transpose()*slw.ps_counts[nm]).sum(axis=1)
            else:
               slw.frozen_pix_counts   =  slw.ps_pix_counts.transpose()[nm].transpose()*slw.ps_counts[nm]
         else:
            slw.frozen_total_counts = 0
            slw.frozen_pix_counts   = 0
         

   def update_localization(self,pslw,newdir):
      
      ro  = ROIOverlap()
      rd  = self.ROI_dir()
      
      exp            = pslw.exposure().value

      from skymaps import PsfFunction,PsfSkyFunction

      for i,slw in enumerate(pslw):
   
         sigma,gamma,pa,en  = slw.sl.sigma(),slw.sl.gamma(),slw.sl.band().pixelArea(),slw.energy()
         exposure_ratio     = exp(newdir,en)/exp(rd,en)
         #psf                = PsfFunction(gamma)
         psf                = PsfSkyFunction(newdir,gamma,sigma)
         
         slw.overlaps[0]    = ro(slw,rd,newdir) * exposure_ratio
         if slw.has_pixels:
            #slw.ps_pix_counts[:,0] = N.asarray([psf(newdir,d,sigma) for d in slw.wsdl])*\
            #                         (pa/(2*N.pi*sigma**2)*exposure_ratio)
            slw.ps_pix_counts[:,0] = N.asarray(psf.wsdl_vector_value(slw.wsdl))*(pa/(2*N.pi*sigma**2)*exposure_ratio)
            slw.unfrozen_pix_counts[:,0] = slw.ps_pix_counts[:,0]

         
###====================================================================================================###

class ROIBackgroundManager(ROIModelManager):
   """Manager background models.  Two spatial dependences are supported via pointlike -- a cube-based
      image for the galactic diffuse, and an isotropic model for extragalactic+residual.

      The isotropic model can be any spectral model, and the galactic model can be scaled by any
      spectral model."""

   def init(self):
      self.use_isotropic = True
      self.use_galactic  = True

   def __init__(self,models,spectral_analysis,roi_dir,**kwargs):
      """Take two spectral models for the galactic and isotropic backgrounds, and use the
         exposure and default parameters provided in an instance of SpectralAnalysis."""
      self.init()
      self.__dict__.update(**kwargs)
      self.models,self.spectral_analysis,self.roi_dir = N.asarray(models),spectral_analysis,roi_dir
      self.gal_model = models[0]
      self.iso_model = models[1] if (self.use_galactic and self.use_isotropic) else models[0]

   def __str__(self): return '\n\n'.join([model.__str__() for model in self.models])

   def setup_initial_counts(self,pslw):
      """Take a PointSourceLikelihoodWrapper object which has already had initial setup done by
         ROIAnalysis and calculate the initial background counts in the ROI and each pixel."""
      
      for slw in pslw: #initialize to values appropriate for zero background
         slw.gal_pix_counts = slw.gal_counts = 0
         slw.iso_pix_counts = slw.iso_counts = 0
      
      if self.use_galactic:  self.__setup__(pslw,bg_type='galactic')
      if self.use_isotropic: self.__setup__(pslw,bg_type='isotropic')

   def __setup__(self,pslw,bg_type):
      """A joint internal method to setup the background type specified by bg_type ('galactic' or 'isotropic')."""

      from skymaps import Background,SkyIntegrator
      from Models import PowerLawFlux,Constant

      sh = bg_type[:3] #'gal' or 'iso'
      sa = self.spectral_analysis
      bg = Background(sa.background.__dict__[bg_type+'_diffuse'],sa.exposure.exposure[0],sa.exposure.exposure[1])

      SkyIntegrator.set_tolerance(0.02)
      
      print 'Calculating %s background counts...'%(bg_type)

      for slw in pslw:
         band                            = slw.sl.band()
         model                           = self.__dict__[sh+'_model']
         slw.__dict__[sh+'_ang_e_vals']  = N.empty(5)
         slw.__dict__[sh+'_pix_e_vals']  = N.empty([len(slw.wsdl),5])
         
         for n,e in enumerate(slw.simps_pts):
            bg.setEnergy(e)
            bg.set_event_class(band.event_class())
            slw.__dict__[sh+'_ang_e_vals'][n]   = SkyIntegrator.ss_average(bg,self.roi_dir,slw.radius_in_rad)
            #slw.__dict__[sh+'_pix_e_vals'][:,n] = [bg.value(d,e) for d in slw.wsdl]
            slw.__dict__[sh+'_pix_e_vals'][:,n] = N.asarray(bg.wsdl_vector_value(slw.wsdl))
         
         slw.__dict__[sh+'_ang_e_vals'] *= slw.solid_angle
         slw.__dict__[sh+'_pix_e_vals'] *= band.pixelArea()

         v = model(slw.simps_pts) * slw.simps_vec
         slw.__dict__[sh+'_counts']     = (slw.__dict__[sh+'_ang_e_vals'] * v).sum()
         slw.__dict__[sh+'_pix_counts'] = (slw.__dict__[sh+'_pix_e_vals'] * v).sum(axis=1)
          
   
   def update_counts(self,pslw):
      if self.use_galactic and N.any(self.gal_model.free):
            for slw in pslw:
               v = self.gal_model(slw.simps_pts) * slw.simps_vec
               slw.gal_counts     = (slw.gal_ang_e_vals * v).sum()
               slw.gal_pix_counts = (slw.gal_pix_e_vals * v).sum(axis=1)
      if self.use_isotropic and N.any(self.iso_model.free):
            for slw in pslw:
               v = self.iso_model(slw.simps_pts) * slw.simps_vec
               slw.iso_counts     = (slw.iso_ang_e_vals * v).sum()
               slw.iso_pix_counts = (slw.iso_pix_e_vals * v).sum(axis=1)

###====================================================================================================###

class ROILocalizer(object):

   def __init__(self,roi):
      self.roi = roi
      self.rd  = roi.ps_manager.ROI_dir()
            
   def TSmap(self,skydir):

      self.roi.ps_manager.update_localization(self.roi.fitter.pslw,skydir)
      return (-2)*self.roi.logLikelihood(self.roi.get_parameters())
      #self.roi.ps_manager.update_localization(self.roi.fitter.pslw,self.roi.ps_manager.ROI_dir())

   def dir(self):
      return self.rd

   def errorCircle(self):
      return 0.05 #initial guess

   def __call__(self,v):

      from skymaps import SkyDir,Hep3Vector
      return self.TSmap(SkyDir(Hep3Vector(v[0],v[1],v[2])))