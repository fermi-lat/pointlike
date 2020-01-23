"""
Provides classes for managing point sources and backgrounds for an ROI likelihood analysis.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_managers.py,v 1.35 2011/11/10 04:05:39 lande Exp $

author: Matthew Kerr
"""

import sys
import numpy as N
from skymaps import SkyDir,SkyIntegrator,Background
from roi_diffuse import ROIDiffuseModel,DiffuseSource
from Models import *
from pypsf import *
from roi_bands import *
from abc import abstractmethod

class ROIModelManager(object):
    """Parent class for point source manager and background manager.  Provides universal
        methods to manager parameters and error matrices."""

    def parameters(self):
        if len(self.models)==0: return []
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

    @abstractmethod
    def add_source(self, model, bands): pass

    @abstractmethod
    def del_source(self, which, bands): pass

    @abstractmethod
    def zero_source(self, which, bands): pass

    @abstractmethod
    def unzero_source(self, which, bands): pass
        
###======================================================================###

class ROIPointSourceManager(ROIModelManager):
    """Manage all point sources."""
    
    def __init__(self,point_sources,roi_dir, quiet=False):
        self.roi_dir = roi_dir
        self.point_sources = N.asarray(point_sources)
        self.names = N.asarray([point_source.name for point_source in point_sources])
        self.models = N.asarray([point_source.model for point_source in point_sources])
        self.mask    = N.asarray([True]*len(self.models),bool)
        self.quiet = quiet

    def __len__(self): return len(self.point_sources)

    def __str__(self,verbose=False):
        mask = [True]*len(self.point_sources) if verbose else [N.any(p.model.free) for p in self.point_sources]
        return '\n\n'.join(['%d -- %s fitted with %s\n'%(i,ps.name,ps.model.full_name())+ps.model.__str__()\
                                  for i,ps in enumerate(self.point_sources) if mask[i]])

    def ROI_dir(self): return self.roi_dir

    def name(self): return self.point_sources[0].name

    def setup_initial_counts(self,bands):

        if not self.quiet: print ('.....setting up point sources (%d in ROI)...'%(len(self.point_sources)),)

        roi_dir = self.roi_dir
        overlap = PsfOverlap()

        for i,band in enumerate(bands):

            en,exp = band.e,band.exp.value
            rvals  = N.empty(len(band.wsdl),dtype=float)
            cpsf   = band.psf.cpsf

            # make a first-order correction for exposure variation
            denom = exp(roi_dir,en)
            if denom==0:
                raise Exception('ROIPointSourceManager: exposure is zero for  energy %f'%en)
            band.er = er = N.asarray([exp(ps.skydir,en)/denom for ps in self.point_sources])

            #unnormalized PSF evaluated at each pixel, for each point source
            band.ps_pix_counts = N.empty([len(band.wsdl),len(self.point_sources)])
            for nps,ps in enumerate(self.point_sources):
                cpsf.wsdl_val(rvals,ps.skydir,band.wsdl)
                band.ps_pix_counts[:,nps] = rvals
            band.ps_pix_counts*= band.b.pixelArea()

            # note this will use a numerical integration if the ragged edge is impt.
            band.overlaps      = N.asarray([overlap(band,roi_dir,ps.skydir) for ps in self.point_sources])

            #initial model prediction
            band.ps_counts     = N.asarray([band.expected(model) for model in self.models])*er
    
        if not self.quiet: print ('done!')
    
    def reload_data(self,bands):
        """ Recalculate the PS model for new data.  NB -- only appropriate for Monte Carlo
            since no exposure changes are factored in."""
        for i,band in enumerate(bands):
            cpsf = band.psf.cpsf
            rvals = N.empty(len(band.wsdl))
            # only need to re-calculate the models for new data pixels
            for nps,ps in enumerate(self.point_sources):
                cpsf.wsdl_val(rvals,ps.skydir,band.wsdl)
                band.ps_pix_counts[:,nps] = rvals
            band.ps_pix_counts    *= b.pixelArea()

    def update_counts(self,bands):
        """Update models with free parameters."""
        ma = self.mask
        have_ps = len(self.point_sources) > 0
        have_fs = ma.sum() > 0
        for band in bands:
            new_counts = N.asarray([band.expected(m) for m in self.models[ma]])*band.er[ma]
            band.ps_counts[ma] = new_counts
            band.ps_all_counts = (band.overlaps[ma]*new_counts).sum() + band.frozen_total_counts
            if band.has_pixels:
                if have_ps:
                    if have_fs:
                         band.ps_all_pix_counts = \
                             (band.unfrozen_pix_counts * new_counts).sum(axis=1) + band.frozen_pix_counts
                    else:
                         band.ps_all_pix_counts = band.frozen_pix_counts
                else:
                    band.ps_all_pix_counts = 0

    def cache(self,bands):
        """Cache values for models with no degrees of freedom.  Helps a lot near galactic plane."""
        self.mask = m = N.asarray([N.any(model.free) for model in self.models],bool)
        nm = ~self.mask
        nm_sum  = nm.sum()
        m_sum    = m.sum()
        
        # note on implementation:
        # frozen_pix_counts is always a 1-vector with dimension equal to the number of data pixels in a band
        # unfrozen_pix_counts is always an (npix,n_free_sources) matrix, even in the case where there are
        # no free sources, where it degenerates to a (npix,0)-dim matrix.
        # This structure is desired as the total counts that go into the likelihood, as determined in
        # update_counts, performs a sum over the second (source) axis of unfrozen_pix_counts, giving
        # an (npix,) vector that can be safely added to frozen_pix_counts.

        for band in bands:

            band.ps_counts = N.asarray([band.expected(model) for model in self.models]) * band.er

            if (m_sum > 0) and band.has_pixels:
                band.unfrozen_pix_counts = band.ps_pix_counts.transpose()[m].transpose()
            else:
                band.unfrozen_pix_counts = 0

            if nm_sum > 0:
                band.frozen_total_counts = (band.ps_counts[nm]*band.overlaps[nm]).sum()                    
                if band.has_pixels:
                     band.frozen_pix_counts = (band.ps_pix_counts.transpose()[nm].transpose()*band.ps_counts[nm]).sum(axis=1)
                else:
                     band.frozen_pix_counts = 0

            else:
                band.frozen_total_counts = 0
                band.frozen_pix_counts    = 0

              
    # comment out since so dependent on representaion
    #def add_cutoff(self,which, free_cutoff=False):
    #    """Replace a power law point source with one with a cutoff.  Useful for catalog pulsars."""

  # #     e = ExpCutoff()
    #    m = self.models[which]

  # #     e.p[0]        = N.log10(m(1000.))
    #    e.p[1]        = m.p[1]
    #    e.p[2]        = N.log10(5e5)
    #    e.free[:2]  = m.free[:2]
    #    e.free[2]    = free_cutoff
    #    e.e0          = m.e0

  # #     self.models[which] = e
    #    self.point_sources[which].model = e

    def add_source(self, ps, bands, nosort=False):
        """Add a new PointSource object to the model and re-calculate the point source
           contribution. When adding the source, resort the point sources
           based upon their distance from the center. 
           
           Use the flag nosort to not sort the point sources after adding the new one
           (for backwards compatability). """

        if nosort:
            sort_index=N.arange(len(self.point_sources)+1)
        else:
            sd = [ i.skydir for i in self.point_sources ] + [ ps.skydir ]
            sort_index=N.argsort([self.roi_dir.difference(i) for i in sd])

        self.point_sources = N.append(self.point_sources,ps)[sort_index]
        self.models        = N.append(self.models,ps.model)[sort_index]
        self.names         = N.append(self.names,ps.name)[sort_index]
        self.mask          = N.append(self.mask,N.any(ps.model.free))[sort_index]
        self.setup_initial_counts(bands) # not most efficient, but fewest loc!
        self.cache(bands)

    def del_source(self, which, bands):
        ops = self.point_sources[which]
        self.point_sources = N.delete(self.point_sources,which)
        self.models        = N.delete(self.models,which)
        self.names         = N.delete(self.names,which)
        self.mask          = N.delete(self.mask,which)

        self.setup_initial_counts(bands) # ditto
        self.cache(bands)
        return ops

    def zero_source(self, which, bands):
        m = self.models[which]
        m.zero()
        self.cache(bands)

    def unzero_source(self, which, bands):
        m = self.models[which]
        try:
            m.unzero()
            self.cache(bands)      
        except:
            print ('Source %d indicated was not zeroed in the first place!' %which)

###======================================================================###

class ROIDiffuseManager(ROIModelManager):
    """ Manage a set of ROIDiffuseModels as they interact with the likelihood."""

    def init(self):
        self.quiet  = False

    def __init__(self,models,roi_dir,**kwargs):
        """ models -- a list of ROIDiffuseModels
            skydir -- the center of the ROI
        """
        self.init()
        self.__dict__.update(**kwargs)

        self.bgmodels = models
        self.names = N.asarray([bgm.name for bgm in self.bgmodels])
        self.models   = N.asarray([bgm.smodel for bgm in self.bgmodels])
        self.diffuse_sources = N.asarray([bgm.diffuse_source for bgm in self.bgmodels])

        self.roi_dir  = roi_dir
      
    def __str__(self): return '\n\n'.join([model.__str__() for model in self.bgmodels])

    def setup_initial_counts(self,bands):
        """ Cause the diffuse models to perform any calculations needed
            to evaluate the diffuse contributions to the ROI. This
            includes evaluating the models in position and energy.
            Typically, this will only be called once, but extended sources
            with shape parameters should make use of it whenever the value
            of a shape parameter changes."""

        nm  = len(self.models)
        rd  = self.roi_dir

        # setup arrays to hold model counts
        for band in bands:
            band.bg_counts = N.empty(nm)
            band.bg_pix_counts = N.empty([len(band.wsdl),nm]) if band.has_pixels else 0
        # initialize models and get inital counts
        if not self.quiet:
            print ('.....setting up diffuse/extended backgrounds for %d bands...'%len(bands))
        for ibg,bg in enumerate(self.bgmodels):
            if not self.quiet:
                print ('.......... %s'%(bg.name) ,)
                sys.stdout.flush()
            bg.initialize_counts(bands)
            bg.update_counts(bands,ibg)


    def update_counts(self,bands):
        """ Update the counts in the bands objects to reflect the current
            value of the spectral scaling model."""

        for ibg,bg in enumerate(self.bgmodels):
            bg.update_counts(bands,ibg)

        for band in bands:
            band.bg_all_counts = band.bg_counts.sum()
            if band.has_pixels:
                band.bg_all_pix_counts = band.bg_pix_counts.sum(axis=1)

    def gradient(self,bands):
        """ Calculate the gradient for the free spectral parameters.
            Requires input from *all* sources, so there is necessarily
            some loss of abstraction.  The band objects *must* have a
            member, pix_weights, giving the ratio of the observed counts
            to the total expected counts from all models for the pixel.
        """

        return N.concatenate([bg.gradient(bands,ibg)
                              for ibg,bg in enumerate(self.bgmodels)])

    def add_source(self, model, bands):
        """Add a new diffuse object to the model and re-calculate the diffuse source
            contribution. Model must be an roi_diffuse.ROIDiffuseModel object. """
        if not isinstance(model,ROIDiffuseModel):
            raise Exception("Can only add to ROI an ROIDiffuseModel object.")

        self.bgmodels        = N.append(self.bgmodels,model)
        self.models          = N.append(self.models,model.smodel)
        self.names           = N.append(self.names,model.name)
        self.diffuse_sources = N.append(self.diffuse_sources,model.diffuse_source)

        # Add in this one source to the band, instead of recalculating
        # all other background models.
        index = len(self.models)-1
        for band in bands:
            band.bg_counts = N.append(band.bg_counts,N.empty(1))
            band.bg_pix_counts = N.append(band.bg_pix_counts, N.empty((len(band.wsdl),1)),axis=1) if band.has_pixels else 0

        self.bgmodels[index].initialize_counts(bands)
        self.update_counts(bands)

    def del_source(self, which, bands):
        """ which must be an index to the desired diffuse source 
            to delete. """
        ops = self.diffuse_sources[which]
        self.bgmodels        = N.delete(self.bgmodels,which)
        self.models          = N.delete(self.models,which)
        self.names           = N.delete(self.names,which)
        self.diffuse_sources = N.delete(self.diffuse_sources,which)

        for band in bands:
            band.bg_counts = N.delete(band.bg_counts,which)
            band.bg_pix_counts = N.delete(band.bg_pix_counts, which, axis=1)

        self.update_counts(bands)
        return ops

    def zero_source(self, which, bands):
        m = self.models[which]
        m.zero()
        self.update_counts(bands)

    def unzero_source(self, which, bands):
        m = self.models[which]
        try:
            m.unzero()
            self.update_counts(bands)
        except:
            print ('Source indicated was not zeroed in the first place!')

###======================================================================###

class ROIBackgroundManager(ROIModelManager):
    """Manage.  The input is a set of diffuse models (SkySpectrum instances) and a matching set
        of position-independent energy spectral models with which to scale the diffuse models, e.g.,
        an overall scale or a power law."""

    def init(self):
        self.nsimps = 4
        self.quiet  = False
        self.on_the_fly = False

    def __init__(self,spectral_analysis,models,skydir,**kwargs):
        """."""
        self.init()
        self.__dict__.update(**kwargs)
        self.sa         = sa = spectral_analysis


        self.bgmodels = models
        self.models    = N.asarray([bgm.smodel for bgm in self.bgmodels])
        for model in self.models:
            model.background = True

        self.roi_dir  = skydir
        
    def __str__(self): return '\n\n'.join([model.__str__() for model in self.bgmodels])

    def setup_initial_counts(self,bands):
        """Evaluate initial values of background models; these will be scaled in likelihood maximization."""

        if not self.quiet:
            print ('.....setting up diffuse backgrounds...')
            for bg in self.bgmodels:
                print ('..........using %s'%(bg.name))


        ns  = self.nsimps
        nm  = len(self.models)
        rd  = self.roi_dir

        constructor = ROIBackgroundCountsOTF if self.on_the_fly else ROIBackgroundCounts
        self.mybgs = bgs = [constructor(self,bgm) for bgm in self.bgmodels]

        for nband,band in enumerate(bands):
            
            band.bg_points = sp = N.logspace(N.log10(band.emin),N.log10(band.emax),ns+1)
            band.bg_vector = sp * (N.log(sp[-1]/sp[0])/(3.*ns)) * \
                                  N.asarray([1.] + ([4.,2.]*(ns/2))[:-1] + [1.])

            #figure out best way to handle no pixel cases...
            band.ap_evals  = N.empty([nm,ns + 1])        
            band.pi_evals  = N.empty([len(band.wsdl),nm,ns + 1]) if band.has_pixels else 0

                        
            for ne,e in enumerate(band.bg_points):
                for nbg,bg in enumerate(bgs):
                    bg.set_state(e,band.ct)
                    band.ap_evals[nbg,ne] = bg.ap_value(rd,band.radius_in_rad)
                    if band.has_pixels:
                        band.pi_evals[:,nbg,ne] = bg.pix_value(band.wsdl)

            band.ap_evals *= (band.solid_angle    * band.bg_vector)
            band.pi_evals *= (band.b.pixelArea() * band.bg_vector)

            band.mo_evals = N.empty([nm,ns + 1])
            for n,m in enumerate(self.models):
                band.mo_evals[n,:] = m(band.bg_points)

            #discrepancy between initial calculation and updated -- can they be made consistent/better?    
            band.bg_counts      = (band.ap_evals * band.mo_evals).sum(axis = 1)
            band.bg_all_counts = band.bg_counts.sum()
            if band.has_pixels:
                band.bg_pix_counts = (band.pi_evals * band.mo_evals).sum(axis = 2) if band.has_pixels else 0
                band.bg_all_pix_counts = band.bg_pix_counts.sum(axis=1)

            # temporary code to try to address NaN problem
            band.backup_bg_counts      = band.bg_counts.copy()
            if band.has_pixels:
                band.backup_bg_pix_counts = band.bg_pix_counts.copy()
            self.init_norms = N.asarray([m[0] for m in self.models])

        self.cache() # check that this doesn't cause problems -- it shouldn't
        if not self.quiet: print ('done!')

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
                    band.bg_counts[nm] = (band.ap_evals[nm,:]    * pts).sum()
                    if band.has_pixels:
                        band.bg_pix_counts[:,nm] = (band.pi_evals[:,nm,:] * pts).sum(axis=1)                    
                        
            else:
                ratio = (m[0]/self.init_norms[nm])
                if ratio == 0 or N.isnan(ratio): print (ratio,m[0])
                #self.init_norms[nm] = m.p[0]
                for band in bands: 
                    #band.bg_counts[nm] *= ratio
                    band.bg_counts[nm] = ratio * band.backup_bg_counts[nm]
                    if band.has_pixels:
                        #band.bg_pix_counts[:,nm] *= ratio
                        band.bg_pix_counts[:,nm] = ratio * band.backup_bg_pix_counts[:,nm]

        for band in bands:
            band.bg_all_counts = band.bg_counts.sum() #inefficient!
            if band.has_pixels: band.bg_all_pix_counts    = band.bg_pix_counts.sum(axis=1)


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
                        
            band.bg_counts      = (band.ap_evals * band.mo_evals).sum(axis = 1)
            band.bg_all_counts = band.bg_counts.sum()
            band.pi_evals *= (band.b.pixelArea() * band.bg_vector)
            band.bg_pix_counts = (band.pi_evals * band.mo_evals).sum(axis = 2)
            band.bg_all_pix_counts = band.bg_pix_counts.sum(axis=1)
