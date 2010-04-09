"""
Provides modules for managing point sources and backgrounds for an ROI likelihood analysis.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_managers.py,v 1.3 2010/03/11 19:23:29 kerrm Exp $

author: Matthew Kerr
"""

import numpy as N
from skymaps import SkyDir,SkyIntegrator,Background
from pointlike import DoubleVector
from Models import *
from pypsf import *
from roi_bands import *

      
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

   def get_free_errors(self):
      return list(N.concatenate([m.get_free_errors() for m in self.models]))
      

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
      return '\n\n'.join(['%d -- %s fitted with %s\n'%(i,ps.name,ps.model.pretty_name)+ps.model.__str__()\
                          for i,ps in enumerate(self.point_sources) if mask[i]])

   def ROI_dir(self): return self.roi_dir

   def name(self): return self.point_sources[0].name

   def setup_initial_counts(self,bands):

      if not self.quiet: print '.....setting up point sources (%d in ROI)...'%(len(self.point_sources)),

      roi_dir = self.roi_dir
      overlap = PsfOverlap()
      dv      = DoubleVector()

      for i,band in enumerate(bands):

         en,exp,pa = band.e,band.exp.value,band.b.pixelArea()

         # make a first-order correction for exposure variation
         denom = exp(roi_dir,en)
         if denom==0:
            raise Exception('ROIPointSourceManager: exposure is zero for  energy %f'%en)
         band.er = er       = N.asarray([exp(ps.skydir,en)/denom for ps in self.point_sources])

         #unnormalized PSF evaluated at each pixel, for each point source
         band.ps_pix_counts = N.empty([len(band.wsdl),len(self.point_sources)])
         for nps,ps in enumerate(self.point_sources):
            #diffs = N.asarray([ps.skydir.difference(w) for w in band.wsdl])
            band.wsdl.arclength(ps.skydir,dv)
            band.ps_pix_counts[:,nps] = band.psf(N.fromiter(dv,dtype=float),density=True)
         band.ps_pix_counts*= pa

         #fraction of PSF contained within the ROI
         band.overlaps      = N.asarray([overlap(band,roi_dir,ps.skydir) for ps in self.point_sources])
         band.overlaps     *= (band.solid_angle/band.solid_angle_p) #small correction for ragged edge
         
         #initial model prediction
         band.ps_counts     = N.asarray([band.expected(model) for model in self.models])*er
   
      if not self.quiet: print 'done!'
   
   def reload_data(self,bands):

      for i,band in enumerate(bands):

         # only need to re-calculate the models for new data pixels
         band.ps_pix_counts    = N.empty([len(band.wsdl),len(self.point_sources)])
         for nps,ps in enumerate(self.point_sources):
            #diffs = N.asarray([ps.skydir.difference(w) for w in band.wsdl])
            band.wsdl.arclength(ps.skydir,dv)
            band.ps_pix_counts[:,nps] = band.psf(N.fromiter(dv,dtype=float),density=True)
         band.ps_pix_counts   *= b.pixelArea()


   def update_counts(self,bands):
      """Update models with free parameters."""
      ma = self.mask
      for band in bands:
         new_counts = N.asarray([band.expected(m) for m in self.models[ma]])*band.er[ma]
         band.ps_counts[ma] = new_counts
         band.ps_all_counts = (band.overlaps[ma]*new_counts).sum() + band.frozen_total_counts
         if band.has_pixels:
            band.ps_all_pix_counts = \
               (band.unfrozen_pix_counts * new_counts).sum(axis=1) + band.frozen_pix_counts

   def cache(self,bands):
      """Cache values for models with no degrees of freedom.  Helps a lot near galactic plane."""
      self.mask = m = N.asarray([N.any(model.free) for model in self.models])
      nm = ~self.mask
      s  = nm.sum()
      
      # note on implementation:
      # frozen_pix_counts is always a 1-vector with dimension equal to the number of data pixels in a band
      # unfrozen_pix_counts is always an (npix,n_free_sources) matrix, even in the case where there are
      # no free sources, where it degenerates to a (npix,0)-dim matrix.
      # This structure is desired as the total counts that go into the likelihood, as determined in
      # update_counts, performs a sum over the second (source) axis of unfrozen_pix_counts, giving
      # an (npix,) vector that can be safely added to frozen_pix_counts.

      for band in bands:
         band.ps_counts = N.asarray([band.expected(model) for model in self.models]) * band.er

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
            band.frozen_pix_counts = (band.ps_pix_counts.transpose()[nm].transpose()*band.ps_counts[nm]).sum(axis=1)

         else:
            band.frozen_total_counts = 0
            band.frozen_pix_counts   = 0

           
   def add_cutoff(self,which):
      """Replace a power law point source with one with a cutoff.  Useful for catalog pulsars."""

      e = ExpCutoff()
      m = self.models[which]

      e.p[0]      = N.log10(m(1000.))
      e.p[1]      = m.p[1]
      e.p[2]      = N.log10(5e5)
      e.free[:2]  = m.free[:2]
      e.free[2]   = False
      e.e0        = m.e0

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

   def __init__(self,diffuse_model,scaling_model=None,name=None):
      """diffuse_model --- a SkySpectrum object describing a diffuse background; can
                           be a tuple of there is a separate version for front & back

         scaling_model --- [None] a spectral model from the Models module
                           if None, default to a simple scaling of the normalization

         name          --- [None] model name for pretty printing
      """

      self.dmodel = diffuse_model
      if scaling_model is None:
         self.smodel = Constant(p=[1],free=[True])
      else:
         self.smodel = scaling_model
      self.name   = name if name is not None else "Diffuse Background"

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
      #self.gal_model = PowerLaw(p=[1,1],free=[True,False],index_offset=1)
      #self.iso_model = PowerLaw(p=[1,1],free=[True,False],index_offset=1)
      self.quiet     = False

   def __init__(self,spectral_analysis,models,skydir,**kwargs):
      """."""
      self.init()
      self.__dict__.update(**kwargs)
      self.sa       = sa = spectral_analysis


      self.bgmodels = models
      self.models   = N.asarray([bgm.smodel for bgm in self.bgmodels])
      for model in self.models:
         model.background = True

      self.roi_dir  = skydir
      
   def __str__(self): return '\n\n'.join([model.__str__() for model in self.bgmodels])

   def setup_initial_counts(self,bands):
      """Evaluate initial values of background models; these will be scaled in likelihood maximization."""

      if not self.quiet:
         print '.....setting up diffuse backgrounds...'
         for bg in self.bgmodels:
            print '..........using %s'%(bg.name)

      SkyIntegrator.set_tolerance(0.02)
      exp = self.sa.exposure.exposure
      ns  = self.nsimps
      nm  = len(self.models)
      rd  = self.roi_dir
        
      front_bgs = [Background(model.get_dmodel(0),exp[0]) for model in self.bgmodels]
      back_bgs  = [Background(model.get_dmodel(1),exp[1]) for model in self.bgmodels]

      for nband,band in enumerate(bands):

         #for bg in bgs: bg.set_event_class(band.ec)
         bgs = back_bgs if band.ct else front_bgs
         
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

         # temporary code to try to address NaN problem
         band.backup_bg_counts     = band.bg_counts.copy()
         if band.has_pixels:
            band.backup_bg_pix_counts = band.bg_pix_counts.copy()
         self.init_norms = N.asarray([m.p[0] for m in self.models])

      self.cache() # check that this doesn't cause problems -- it shouldn't
      if not self.quiet: print 'done!'

   def cache(self):

      self.free_mask  = N.asarray([N.any(m.free) for m in self.models])
      self.edep_mask  = N.asarray([N.any(m.free[1:]) for m in self.models]) # first parameter is norm
      #self.init_norms = N.asarray([m.p[0] for m in self.models])

   def update_counts(self,bands):

      em = self.edep_mask; fm = self.free_mask

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
            if ratio == 0 or N.isnan(ratio): print ratio,m.p[0]
            #self.init_norms[nm] = m.p[0]
            for band in bands: 
               #band.bg_counts[nm] *= ratio
               band.bg_counts[nm] = ratio * band.backup_bg_counts[nm]
               if band.has_pixels:
                  #band.bg_pix_counts[:,nm] *= ratio
                  band.bg_pix_counts[:,nm] = ratio * band.backup_bg_pix_counts[:,nm]

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
