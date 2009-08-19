"""
Provides modules for managing point sources and backgrounds for an ROI likelihood analysis.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/roi_modules.py,v 1.10 2009/08/12 01:20:14 kerrm Exp $

author: Matthew Kerr
"""

import numpy as N
from pointspec import SpectralAnalysis
from skymaps import SkyDir,WeightedSkyDirList,SkyIntegrator
from Models import *
from psf import PSF
#from threading import Thread

class ROIOverlap(object):
   """Routines to calculate how much of the emission of a point source falls onto an ROI."""

   def init(self):
      self.quadrature_tol = 1e-4

   def __init__(self,**kwargs):
      self.init()
      self.__dict__.update(kwargs)


   def __call__(self,band,roi_dir,ps_dir,radius_in_rad=None):
      """Return an array of fractional overlap for a point source at location skydir.
         Note radius arguments are in radians."""

      from skymaps import PsfSkyFunction
      from scipy.integrate import quad
      from math import cos,sin
      
      sigma,gamma,roi_rad = band.s,band.g,radius_in_rad or band.radius_in_rad

      offset = roi_dir.difference(ps_dir)

      if offset < 1e-5:
         roi_u = 0.5 * (roi_rad / sigma)**2
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

         overlap = quad(interior,0,N.pi,epsabs=self.quadrature_tol)[0]/N.pi

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
         overlap = quad(exterior,0,limit,epsabs=self.quadrature_tol)[0]/N.pi

      return overlap

      
###====================================================================================================###

class ROIModelManager(object):
   """Parent class for point source manager and background manager.  Provides universal
      methods to manager parameters and error matrices."""

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

class ROIPointSourceManager(ROIModelManager):
   """Manage all point sources."""
   
   def __init__(self,point_sources,roi_dir, quiet=False):
      self.roi_dir = roi_dir
      self.point_sources = point_sources
      self.models = N.asarray([point_source.model for point_source in point_sources])
      self.mask   = N.asarray([True]*len(self.models))
      self.quiet = quiet

   def __len__(self): return len(self.point_sources)

   def __str__(self,verbose=False):
      mask = [True]*len(self.point_sources) if verbose else [N.any(p.model.free) for p in self.point_sources]
      return '\n\n'.join(['%s fitted with %s\n'%(ps.name,ps.model.pretty_name)+ps.model.__str__()\
                          for i,ps in enumerate(self.point_sources) if mask[i]])

   def ROI_dir(self): return self.roi_dir

   def name(self): return self.point_sources[0].name

   def setup_initial_counts(self,bands):

      if not self.quiet: print 'Setting up point sources...'

      roi_dir = self.roi_dir
      overlap = ROIOverlap()
      from skymaps import PsfSkyFunction

      for i,band in enumerate(bands):

         sigma,gamma,en,exp,pa = band.s,band.g,band.e,band.exp.value,band.b.pixelArea()

         #make a first-order correction for exposure variation
         denom = exp(roi_dir,en)
         if denom==0:
            raise Exception('ROIPointSourceManager: exposure is zero for  energy %f'%en)
         exposure_ratios       = N.asarray([exp(ps.skydir,en)/denom for ps in self.point_sources])

         #unnormalized PSF evaluated at each pixel, for each point source
         band.ps_pix_counts    = N.empty([len(band.wsdl),len(self.point_sources)])
         for nps,ps in enumerate(self.point_sources):
            psf = PsfSkyFunction(ps.skydir,gamma,sigma)
            band.ps_pix_counts[:,nps] = psf.wsdl_vector_value(band.wsdl)
         band.ps_pix_counts   *= ( (pa/(2*N.pi*sigma**2)) *exposure_ratios)

         #fraction of PSF contained within the ROI
         band.overlaps      = N.asarray([overlap(band,roi_dir,ps.skydir) for ps in self.point_sources])
         band.overlaps     *= (band.solid_angle/band.solid_angle_p) #small correction for ragged edge
         band.overlaps     *= exposure_ratios
         
         #initial model prediction
         band.ps_counts        = N.asarray([band.expected(model) for model in self.models])

   
   def reload_data(self,bands):

      roi_dir = self.roi_dir
      from skymaps import PsfSkyFunction

      for i,band in enumerate(bands):

         sigma,gamma,en,exp,pa = band.s,band.g,band.e,band.exp.value,band.b.pixelArea()

         #make a first-order correction for exposure variation
         exposure_ratios       = N.asarray([exp(ps.skydir,en)/exp(roi_dir,en) for ps in self.point_sources])

         #unnormalized PSF evaluated at each pixel, for each point source
         band.ps_pix_counts    = N.empty([len(band.wsdl),len(self.point_sources)])
         for nps,ps in enumerate(self.point_sources):
            psf = PsfSkyFunction(ps.skydir,gamma,sigma)
            band.ps_pix_counts[:,nps] = psf.wsdl_vector_value(band.wsdl)
         band.ps_pix_counts   *= ( (pa/(2*N.pi*sigma**2)) *exposure_ratios)


   def update_counts(self,bands):
      """Update models with free parameters."""
      ma = self.mask
      for band in bands:
         band.ps_counts[ma] = [band.expected(m) for m in self.models[ma]]
         band.ps_all_counts = (band.overlaps[ma]*band.ps_counts[ma]).sum() + band.frozen_total_counts
         if band.has_pixels:
            band.ps_all_pix_counts = \
                  (band.unfrozen_pix_counts * band.ps_counts[ma]).sum(axis=1) + band.frozen_pix_counts

      #for i,model in enumerate(self.models):
      #   if self.mask[i]:
      #      for band in bands: band.ps_counts[i] = band.expected(model)


   def cache(self,bands):
      """Cache values for models with no degrees of freedom.  Helps a lot near galactic plane."""
      self.mask = m = N.asarray([N.any(model.free) for model in self.models])
      nm = N.logical_not(self.mask)
      s = nm.sum()

      for band in bands:
         band.ps_counts = N.asarray([band.expected(model) for model in self.models])

      for band in bands:

         if band.has_pixels:
            band.unfrozen_pix_counts = band.ps_pix_counts.transpose()[m].transpose()
         else:
            band.unfrozen_pix_counts = 0

         if s > 0:
            band.frozen_total_counts = (band.ps_counts[nm]*band.overlaps[nm]).sum()
            if not band.has_pixels:
               band.frozen_pix_counts = 0
               continue
            if s > 1:
               band.frozen_pix_counts   = (band.ps_pix_counts.transpose()[nm].transpose()*band.ps_counts[nm]).sum(axis=1)
            else:
               band.frozen_pix_counts   =  band.ps_pix_counts.transpose()[nm].transpose()*band.ps_counts[nm]
         else:
            band.frozen_total_counts = 0
            band.frozen_pix_counts   = 0

           
   def add_cutoff(self,which):
      """Replace a power law point source with one with a cutoff.  Useful for catalog pulsars."""

      from Models import ExpCutoff
      e = ExpCutoff()
      m = self.models[which]

      e.p[0]      = N.log10(m(1000.))
      e.p[1]      = m.p[1]
      e.p[2]      = N.log10(5e5)
      e.free[:2]  = m.free[:2]
      e.free[2]   = False

      self.models[which] = e
      self.point_sources[which].model = e

   def add_ps(self, ps, bands):
      """Add a new PointSource object to the model and re-calculate the point source
         contribution."""

      self.point_sources = N.append(self.point_sources,ps)
      self.models        = N.append(self.models,ps.model)
      self.setup_initial_counts(bands) # not most efficient, but fewest loc!

   def del_ps(self, which, bands):
      ops = self.point_sources[which]
      self.point_sources = N.delete(self.point_sources,which)
      self.models        = N.delete(self.models,which)
      self.setup_initial_counts(bands) # ditto
      return ops

   def zero_ps(self, which, bands):
      m = self.models[which]
      m.old_flux = m.p[0]
      m.p[0] = -100
      m.old_free = m.free.copy()
      m.free[:] = False
      self.cache(bands)

   def unzero_ps(self, which, bands):
      m = self.models[which]
      try:
         m.p[0] = m.old_flux
         m.free = m.old_free.copy()
         self.cache(bands)     
      except:
         print 'Source indicated was not zeroed in the first place!'



###====================================================================================================###

class ROIBackgroundModel(object):
   """A wrapper to neatly associate spectral scaling models with their SkySpectrum instances."""

   def __init__(self,diffuse_model,scaling_model,name):
      """diffuse_model --- a SkySpectrum object describing a diffuse background; can
                           be a tuple of there is a separate version for front & back

         scaling_model --- a spectral model from the Models module; only one of these

         name          --- model name for pretty printing
      """

      self.dmodel = diffuse_model
      self.smodel = scaling_model
      self.name   = name

   def get_dmodel(self,event_class=0):
      if len(N.ravel(self.dmodel)) == 1: return self.dmodel
      else:
         return self.dmodel[event_class]
      
   def __str__(self):
      return '%s scaled with %s\n'%(self.name,self.smodel.pretty_name)+self.smodel.__str__()



         
###====================================================================================================###

class ROIBackgroundManager(ROIModelManager):
   """Manage.  The input is a set of diffuse models (SkySpectrum instances) and a matching set
      of position-independent energy spectral models with which to scale the diffuse models, e.g.,
      an overall scale or a power law."""

   def init(self):
      self.nsimps = 4

      # remove this?
      from Models import PowerLaw
      self.gal_model = PowerLaw(p=[1,1],free=[True,False],index_offset=1)
      self.iso_model = PowerLaw(p=[1,1],free=[True,False],index_offset=1)
      self.quiet     = False

   def __init__(self,spectral_analysis,models=None,**kwargs):
      """."""
      self.init()
      self.__dict__.update(**kwargs)
      self.sa       = sa = spectral_analysis

      if models is None: #default background model; note indices are fixed
         if sa.fb_maps:
            from skymaps import DiffuseFunction,IsotropicSpectrum
            from wrappers import Singleton
            self.gfsingl   = Singleton(DiffuseFunction,'front_diffuse',sa.ae.galactic_front)
            self.gbsingl   = Singleton(DiffuseFunction,'back_diffuse', sa.ae.galactic_back)
            self.gal_front = self.gfsingl('front_diffuse')
            self.gal_back  = self.gbsingl('back_diffuse')

            gal_model  = ROIBackgroundModel([self.gal_front,self.gal_back],
                                            self.gal_model, 'Galactic Diffuse')

            self.iso_front = IsotropicSpectrum(sa.ae.isotropic_front)
            self.iso_back  = IsotropicSpectrum(sa.ae.isotropic_back)

            iso_model  = ROIBackgroundModel([self.iso_front,self.iso_back],
                                            self.iso_model, 'Isotropic Diffuse')

         else:
            gal_model  = ROIBackgroundModel(sa.background.galactic_diffuse,
                                            self.gal_model, 'GalacticDiffuse')

            iso_model  = ROIBackgroundModel(sa.background.isotropic_diffuse,
                                            self.iso_model, 'Isotropic Diffuse')

         models     = [gal_model,iso_model]

      self.bgmodels = models
      self.models   = N.asarray([bgm.smodel for bgm in self.bgmodels])
      for model in self.models:
         model.background = True

      self.roi_dir  = sa.roi_dir
      
   def __str__(self): return '\n\n'.join([model.__str__() for model in self.bgmodels])

   def setup_initial_counts(self,bands):
      """Evaluate initial values of background models; these will be scaled in likelihood maximization."""

      if not self.quiet: print 'Calculating initial background counts...'

      from skymaps import Background,SkyIntegrator

      SkyIntegrator.set_tolerance(0.02)
      exp = self.sa.exposure.exposure
      ns  = self.nsimps
      nm  = len(self.models)
      rd  = self.roi_dir
        
      front_bgs = [Background(model.get_dmodel(0),exp[0]) for model in self.bgmodels]
      back_bgs  = [Background(model.get_dmodel(1),exp[1]) for model in self.bgmodels]

      for nband,band in enumerate(bands):

         #for bg in bgs: bg.set_event_class(band.ec)
         bgs = back_bgs if band.ec else front_bgs
         
         band.bg_points = sp = N.logspace(N.log10(band.emin),N.log10(band.emax),ns+1)
         band.bg_vector = sp * (N.log(sp[-1]/sp[0])/(3.*ns)) * \
                          N.asarray([1.] + ([4.,2.]*(ns/2))[:-1] + [1.])

         #figure out best way to handle no pixel cases...
         band.ap_evals  = N.empty([nm,ns + 1])      
         band.pi_evals  = N.empty([len(band.wsdl),nm,ns + 1]) if band.has_pixels else 0

                  
         for ne,e in enumerate(band.bg_points):
            for nbg,bg in enumerate(bgs):
               bg.setEnergy(e)
               band.ap_evals[nbg,ne]      = SkyIntegrator.ss_average(bg,rd,band.radius_in_rad)
               if band.has_pixels:
                  band.pi_evals[:,nbg,ne] = N.asarray(bg.wsdl_vector_value(band.wsdl))

         band.ap_evals *= (band.solid_angle   * band.bg_vector)
         band.pi_evals *= (band.b.pixelArea() * band.bg_vector)

         band.mo_evals = N.empty([nm,ns + 1])
         for n,m in enumerate(self.models):
            band.mo_evals[n,:] = m(band.bg_points)

         #discrepancy between initial calculation and updated -- can they be made consistent/better?   
         band.bg_counts     = (band.ap_evals * band.mo_evals).sum(axis = 1)
         band.bg_all_counts = band.bg_counts.sum()
         if band.has_pixels:
            band.bg_pix_counts = (band.pi_evals * band.mo_evals).sum(axis = 2) if band.has_pixels else 0
            band.bg_all_pix_counts = band.bg_pix_counts.sum(axis=1)

      self.cache() # check that this doesn't cause problems -- it shouldn't


   def cache(self):

      self.free_mask  = N.asarray([N.any(m.free) for m in self.models])
      self.edep_mask  = N.asarray([N.any(m.free[1:]) for m in self.models])
      self.init_norms = N.asarray([m.p[0] for m in self.models])
          
   """
   class updateCountsThread(Thread):

      def __init__(self,bands,psm):
         Thread.__init__(self)
         self.bands = bands
         self.psm   = psm

      def run(self):
         self.psm.update_counts(self.bands)
   """

   def update_counts(self,bands):

      em,fm = self.edep_mask,self.free_mask

      if not N.any(fm): return
            
      for nm,m in enumerate(self.models):
         if not fm[nm]: continue
         if em[nm]:
            for band in bands:
               pts = m(band.bg_points)
               band.bg_counts[nm] = (band.ap_evals[nm,:]   * pts).sum()
               if band.has_pixels:
                  band.bg_pix_counts[:,nm] = (band.pi_evals[:,nm,:] * pts).sum(axis=1)               
                  
         else:
            ratio = 10**(m.p[0]-self.init_norms[nm])
            self.init_norms[nm] = m.p[0]
            for band in bands: 
               band.bg_counts[nm] *= ratio
               if band.has_pixels:
                  band.bg_pix_counts[:,nm] *= ratio               

      for band in bands:
         band.bg_all_counts = band.bg_counts.sum() #inefficient!
         if band.has_pixels: band.bg_all_pix_counts   = band.bg_pix_counts.sum(axis=1)

   """

   def update_counts(self,bands):
      
      #cache the check for free? also, could be done more elegantly with a mask
      #also, want to implement a check for models that can be done with scaling only!
      
      for nm,m in enumerate(self.models):
         if N.any(m.free):
            for band in bands:
               pts = m(band.bg_points)
               band.bg_counts[nm] = (band.ap_evals[nm,:]   * pts).sum()
               band.bg_all_counts = band.bg_counts.sum() #inefficient!
               if band.has_pixels:
                  band.bg_pix_counts[:,nm] = (band.pi_evals[:,nm,:] * pts).sum(axis=1)               
                  band.bg_all_pix_counts   = band.bg_pix_counts.sum(axis=1)

   """

   def reload_data(self,bands):

      from skymaps import Background

      exp = self.sa.exposure.exposure
      ns  = self.nsimps
      nm  = len(self.models)
      rd  = self.roi_dir

      front_bgs = [Background(model.get_dmodel(0),exp[0]) for model in self.bgmodels]
      back_bgs  = [Background(model.get_dmodel(1),exp[1]) for model in self.bgmodels]

      for nband,band in enumerate(bands):

         if not band.has_pixels:
            band.pi_evals = 0
            band.bg_pix_counts = 0
            continue

         #for bg in bgs: bg.set_event_class(band.ec)
         bgs = back_bgs if band.ec else front_bgs
         
         band.pi_evals  = N.empty([len(band.wsdl),nm,ns + 1])
                  
         for ne,e in enumerate(band.bg_points):
            for nbg,bg in enumerate(bgs):
               bg.setEnergy(e)
               band.pi_evals[:,nbg,ne] = N.asarray(bg.wsdl_vector_value(band.wsdl))

         #there really needs to be a reconciliation between generated and updating
         band.mo_evals = N.empty([nm,ns + 1])
         for n,m in enumerate(self.models):
            band.mo_evals[n,:] = m(band.bg_points)
                  
         band.bg_counts     = (band.ap_evals * band.mo_evals).sum(axis = 1)
         band.bg_all_counts = band.bg_counts.sum()
         band.pi_evals *= (band.b.pixelArea() * band.bg_vector)
         band.bg_pix_counts = (band.pi_evals * band.mo_evals).sum(axis = 2)
         band.bg_all_pix_counts = band.bg_pix_counts.sum(axis=1)



###====================================================================================================###

class ROIBand(object):
   """Wrap a Band object, and provide additional functionality for likelihood."""

   
   def init(self):

      self.umax      = 50

      self.nsp_simps = 8
      self.nbg_simps = 4

      self.catalog_aperture = -1


   def __init__(self,band,spectral_analysis,**kwargs):

      self.init()
      self.__dict__.update(**kwargs)
      self.b   = band
      self.s   = band.sigma()
      self.g   = band.gamma()
      self.emin, self.emax = band.emin(),band.emax()
      self.e   = (self.emin*self.emax)**0.5
      self.sa  = spectral_analysis
      self.sd  = self.sa.roi_dir
      self.ec  = band.event_class()
      self.exp = self.sa.exposure.exposure[self.ec]

      self.__setup_data__()
      self.__setup_sp_simps__()


   def __setup_data__(self):
      """Get all pixels within the ROI in this band."""

      from skymaps import WeightedSkyDirList
      mi,ma              = N.asarray([self.sa.minROI,self.sa.maxROI])*(N.pi/180.)
      self.radius_in_rad = max(min((2*self.umax)**0.5*self.s,ma),mi)
      ###### begin PSR cat code
      if self.catalog_aperture > 0:
         th = 0.8*(self.e/1000.)**-0.75
         if th > self.catalog_aperture: th = self.catalog_aperture
         if th < 0.35: th = 0.35
         self.radius_in_rad = th * N.pi / 180
      ###### end PSR cat code
      self.wsdl          = WeightedSkyDirList(self.b,self.sd,self.radius_in_rad,False)
      self.pix_counts    = N.asarray([x.weight() for x in self.wsdl]) if len(self.wsdl) else 0.
      self.photons       = self.wsdl.counts()
      self.has_pixels    = self.photons > 0
      self.npix          = self.wsdl.total_pix()
      self.solid_angle   = self.npix*self.b.pixelArea() #ragged edge
      self.solid_angle_p = 2*N.pi*(1-N.cos(self.radius_in_rad)) #solid angle for a pure circle
             

   def __setup_sp_simps__(self):
      """Cache factors for quickly evaluating the counts under a given spectral model."""

      from pointlike import DoubleVector
      self.sp_points = sp = N.logspace(N.log10(self.emin),N.log10(self.emax),self.nsp_simps+1)
      exp_points     = N.asarray(self.exp.vector_value(self.sd,DoubleVector(sp)))
      simps_weights  = (N.log(sp[-1]/sp[0])/(3.*self.nsp_simps)) * \
                       N.asarray([1.] + ([4.,2.]*(self.nsp_simps/2))[:-1] + [1.])
      self.sp_vector = sp * exp_points * simps_weights


   def reload_data(self,band):
      """For Monte Carlo, when everything stays same but realization of data."""
      self.b = band
      self.__setup_data__()

   def expected(self,model):
      """Integrated the passed spectral model over the exposure and return expected counts."""
      
      return (model(self.sp_points)*self.sp_vector).sum()

   def bandLikelihood(self, parameters, *args):
      """Implement a model independent likelihood for the number of counts of a particular source.
         Other sources (diffuse, neighboring point sources) are treated as fixed."""

      new_counts = parameters[0]
      which = args[0] if len(args) > 0 else 0
      band = self

      old_counts = band.ps_counts[which]

      tot_term = (self.bg_all_counts + self.ps_all_counts + self.overlaps[which]*(new_counts - old_counts))*self.phase_factor

      pix_term = (
                  self.pix_counts * 
                     N.log(
                        self.phase_factor*(
                           self.bg_all_pix_counts + self.ps_all_pix_counts + self.ps_pix_counts[:,which]*(new_counts - old_counts)
                        )
                     )
                 ).sum() if self.has_pixels else 0.

      return tot_term - pix_term

   # note no phase factor!
   def loglikelihood(self,tot_only=False,pix_only=False,tl=False):

      tot = self.bg_all_counts + self.ps_all_counts
      
      if (tot_only):# or (not self.has_pixels):
         return (self.photons*N.log(tot) - tot) if tl else tot #non extended likelihood
      
      if self.has_pixels:
         pix = (self.pix_counts * N.log(self.bg_all_pix_counts + self.ps_all_pix_counts)).sum()
      else:
         pix = 0

      if pix_only: return pix #log likelihood for pixels only

      else: return tot - pix #-log likelihood

###====================================================================================================###

class ROIEnergyBand(object):
   """Wrap 1 or 2 ROIBand objects corresponding to the same energy level 
      but different conversion classes.  Implement a likelihood as a
      function of energy."""

   def __init__(self,bands,emin=None,emax=None):

      self.bands = bands
      self.emin = self.bands[0].emin if emin is None else emin
      self.emax = self.bands[0].emax if emax is None else emax

   def bandLikelihood(self,parameters,*args):
      m = args[0]
      m.set_parameters(parameters)
      return sum( (b.bandLikelihood([b.expected(m)],*args[1:]) for b in self.bands) )

   def bandFit(self,which=0):
      """Fit a model-independent flux to a point source."""

      self.m = PowerLaw(free=[True,False],e0=(self.emin*self.emax)**0.5) # fix index to 2
      from scipy.optimize import fmin,fsolve
      f = self.bandLikelihood

      self.fit = fmin(f,self.m.get_parameters(),disp=0,full_output=1,args=(self.m,which))

      def upper_limit():

         flux_copy = self.m.p[0]
         zp        = self.bandLikelihood(N.asarray([-20]),self.m,which)
         def f95(parameters):
            return abs(self.bandLikelihood(parameters,self.m,which) - zp - 1.92)
         
         # for some reason, can't get fsolve to work here.  good ol' fmin to the rescue
         self.uflux = 10**fmin(f95,N.asarray([-11.75]),disp=0)[0]
         self.lflux = None
         self.flux  = None

         self.m.p[0] = flux_copy

      # if flux below a certain level, set an upper limit
      if self.m.p[0] < -20:
         upper_limit()

      else:
         from specfitter import SpectralModelFitter
         from numpy.linalg import inv      
         hessian = SpectralModelFitter.hessian(self.m,self.bandLikelihood,which) #does Hessian for free parameters
         self.m.set_cov_matrix(inv(hessian))

         e = self.m.statistical(absolute=True,two_sided=True)
         self.flux  = e[0][0]
         self.uflux = self.flux + e[1][0]
         self.lflux = self.flux - e[2][0]

         if self.uflux / self.lflux > 1e3:
            upper_limit()

###====================================================================================================###

class NSMap(object):
   """Attempt to add and fit a new source to an existing model and display the
      results as a TS map."""

   def __init__(self,roi):
      self.roi = roi
      self.phase_factor = roi.phase_factor
      self.ro = ROIOverlap()
      self.mo = PowerLaw(free=[True,False])
      self.ll = roi.logLikelihood(roi.parameters())

   def __call__(self,v):
      from skymaps import SkyDir,Hep3Vector,PsfSkyFunction
      skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
      #skydir = v
      ro = self.ro
      rd = self.roi.sa.roi_dir

      self.mo.p[0] = -11.5

      # calculate new spatial structure for source
      for i,band in enumerate(self.roi.bands):

         sigma,gamma,en,exp,pa = band.s,band.g,band.e,band.exp,band.b.pixelArea()
         exposure_ratio        = exp.value(skydir,en)/exp.value(rd,en) #needs fix? -- does not obey which?
         psf                   = PsfSkyFunction(skydir,gamma,sigma)
         
         band.ns_overlap       = ro(band,rd,skydir) * exposure_ratio * band.solid_angle / band.solid_angle_p

         if band.has_pixels:
            band.ns_ps_pix_counts = N.asarray(psf.wsdl_vector_value(band.wsdl))*((pa/(2*N.pi*sigma**2))*exposure_ratio)

         band.base_counts = band.expected(self.mo)
      from scipy.optimize import fmin
      f = fmin(self.logLikelihood,N.asarray([-11.5]),full_output=1,disp=0)
      return 2*(self.ll - f[1])

   def logLikelihood(self,parameters,*args):

      bands = self.roi.bands
      ll    = 0
      mo    = self.mo
      mo.set_parameters(parameters)

      for b in bands:

         new_ps_counts = b.base_counts * 10**(mo.p[0] + 11.5)


         ll +=  ( 
                   #integral terms for ROI (go in positive)
                   (b.bg_all_counts + b.ps_all_counts + new_ps_counts*b.ns_overlap)*self.phase_factor

                   -

                   #pixelized terms (go in negative)
                   (b.pix_counts *
                       N.log(  (b.bg_all_pix_counts + b.ps_all_pix_counts + new_ps_counts*b.ns_ps_pix_counts)*self.phase_factor )
                   ).sum() if b.has_pixels else 0.
                )

      return 1e6 if N.isnan(ll) else ll


###====================================================================================================###

class ROILocalizer(object):

   def __init__(self,roi,which=0):
      self.roi,self.which = roi, which
      self.rd  = roi.psm.point_sources[which].skydir #note -- not necessarily ROI center!
      self.tsref=0
      self.tsref = self.TSmap(self.rd)
            
   def TSmap(self,skydir):
      return (-2)*self.roi.spatialLikelihood(skydir,which=self.which)-self.tsref

   def dir(self):
      return self.rd

   def errorCircle(self):
      return 0.05 #initial guess

   def __call__(self,v):
      from skymaps import SkyDir,Hep3Vector
      return self.TSmap(SkyDir(Hep3Vector(v[0],v[1],v[2])))