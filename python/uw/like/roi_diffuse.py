"""
Provides classes to encapsulate and manipulate diffuse sources.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_diffuse.py,v 1.1 2010/05/18 22:19:28 kerrm Exp $

author: Matthew Kerr
"""
import numpy as N
from uw.utilities.convolution import BackgroundConvolution
from skymaps import SkyIntegrator,Background

class SmallBand(object):
    """ A little holder."""
    pass

###=========================================================================###
class DiffuseSource(object):
    """ Associate a spatial model with a spectral scaling model."""
    __counter = 0

    def __init__(self,diffuse_model,scaling_model,name=None):

        self.dmodel = diffuse_model
        self.smodel = scaling_model

        if name is None:
            self.name = 'Diffuse Source %d'%(DiffuseSource.__counter)
            DiffuseSource.__counter += 1
        else: self.name = name

        if not hasattr(self.dmodel,'__len__'):
            self.dmodel = [self.dmodel]
   
    def __str__(self): return self.name


###=========================================================================###
class ROIDiffuseModel(object):
    """ Associate a SkySpectrum with a spectral model.
    
        Provide the interface that ROIDiffuseModel et al. should satisfy."""

    def __init__(self,spectral_analysis,diffuse_source,roi_dir,name=None,*args,**kwargs):

        self.sa = spectral_analysis
        self.roi_dir = roi_dir
        self.dmodel  = diffuse_source.dmodel
        self.smodel  = diffuse_source.smodel
        self.name    = diffuse_source.name

        self.init()
        self.__dict__.update(kwargs)
        self.setup()
   
    def init(self):  pass
    def setup(self): pass
    
    def __str__(self):
        return '%s scaled with %s\n'%(self.name,self.smodel.pretty_name)+self.smodel.__str__()

    def get_dmodel(self,event_class=0):
        return self.dmodel[event_class if (len(self.dmodel) > 1) else 0]

    def initialize_counts(self,bands,roi_dir=None):
        """ This method is responsible for establishing a state from
            which the model counts can be calculated.  E.g., evaluating
            the model on a series of energy subplanes and storing the
            results for a later Simpson's rule.

            The user may override the ROI direction that was used to
            setup the objects.  E.g., if the user is localizing an
            extended source.
        """
        raise NotImplementedError,'Classes must implement this method!'

    def update_counts(self,bands,model_index):
        """ This method *must* set the following members of each band:
            band.bg_counts[model_index] -- the total counts expected for
                    the model in the aperture
            band.bg_pix_counts[:,model_index] -- if the band has pixels,
                    an npixel vector with the expected counts from the
                    model for each data pixel
        """
        raise NotImplementedError,'Classes must implement this method!'

    def gradient(self,bands,model_index,phase_factor=1):
        """ This method should return the gradient with respect to the
            free parameters of the spectral model.
        """
        raise NotImplementedError,'Classes must implement this method!'


###=========================================================================###

class ROIDiffuseModel_OTF(ROIDiffuseModel):
    """ Use an on-the-fly numerical scheme to convolve the model
        with the PSF.
        
        Use a Simpson's rule integration over the exposure.
        
        To use a different convolution scheme but maintain the Simpson's
        rule exposure integral, child methods may override __ap_value 
        and __pix_value."""

    def init(self):
        self.pixelsize = 0.5
        self.npix      = 51
        self.nsimps    = 4   # consider making a function?

    def setup(self):
        exp = self.sa.exposure.exposure; psf = self.sa.psf
        if len(self.dmodel) == 1:
            self.bg  = [Background(self.dmodel[0],exp[0],exp[1])]
            self.bgc = [BackgroundConvolution(self.roi_dir,self.bg[0],psf,
                        npix=self.npix,pixelsize=self.pixelsize)]
        else:
            self.bg  = map(Background,self.dmodel,exp)
            self.bgc = [BackgroundConvolution(self.roi_dir,bg,psf,
                        npix=self.npix,pixelsize=self.pixelsize) for bg in self.bg]

    def set_state(self,energy,conversion_type):
        self.active_bgc = self.bgc[conversion_type if (len(self.bgc) > 1) else 0]
        self.active_bgc.do_convolution(energy,conversion_type)

    def initalize_counts(self,bands,roi_dir=None):
        ns = self.nsimps; rd = self.roi_dir if roi_dir is None else roi_dir
        self.bands = [SmallBand() for i in xrange(len(bands))]

        for myband,band in zip(self.bands,bands):
            myband.bg_points = sp = N.logspace(N.log10(band.emin),N.log10(band.emax),ns+1)
            myband.bg_vector = sp * (N.log(sp[-1]/sp[0])/(3.*ns)) * \
                                     N.asarray([1.] + ([4.,2.]*(ns/2))[:-1] + [1.])

            #figure out best way to handle no pixel cases...
            myband.ap_evals  = N.empty(ns + 1)      
            myband.pi_evals  = N.empty([len(band.wsdl),ns + 1]) if band.has_pixels else 0

                  
            for ne,e in enumerate(myband.bg_points):
                self.set_state(e,band.ct)
                myband.ap_evals[ne] = self.__ap_value(rd,band.radius_in_rad)
                if band.has_pixels:
                    myband.pi_evals[:,ne] = self.__pix_value(band.wsdl)

            myband.ap_evals *= (band.solid_angle   * myband.bg_vector)
            myband.pi_evals *= (band.b.pixelArea() * myband.bg_vector)

            # calculate integral counts for later
            myband.mo_evals  = self.smodel(myband.bg_points)
            myband.ap_counts = (myband.ap_evals * myband.mo_evals).sum()
            if band.has_pixels:
                myband.pi_counts = (myband.pi_evals * myband.mo_evals).sum(axis=1)

        self.init_norm = self.smodel.p[0]

    def update_counts(self,bands,model_index):

        mi = model_index
        sm = self.smodel
        if not N.any(sm.free): return

        # model is not energy-dependent
        if not N.any(sm.free[1:]):
            ratio = 10**(sm.p[0]-self.init_norm)
            for myband,band in zip(self.bands,bands): 
               band.bg_counts[mi] = ratio * myband.ap_counts
               if band.has_pixels:
                  band.bg_pix_counts[:,mi] = ratio * myband.pi_counts

        else:
            for myband,band in zip(self.bands,bands):
                pts = sm(myband.bg_points)
                band.bg_counts[mi] = (myband.ap_evals * pts).sum()
                if band.has_pixels:
                    band.bg_pix_counts[:,mi] = (myband.pi_evals * pts).sum(axis=1)               
                  

    def __ap_value(self,center,radius):
        """ Return the integral of the model over the aperture for the
            current state (energy/conversion type).
            Convolution could/should be done.
        """
        return self.active_bgc.ap_average(radius)

    def __pix_value(self,pixlist):
        """ Return the model evaluated at each data pixel for the current
            state (energy/conversion type).
            Convolution could/should be done.
        """
        return self.active_bgc(pixlist,self.active_bgc.cvals)

    def gradient(self,bands,model_index,phase_factor=1):
        sm  = self.smodel
        np  = len(sm.p)
        nfp = sm.free.sum()

        # special case -- no free parameters
        if nfp == 0: return []

        # special case -- only normalization free
        if (nfp == 1) and sm.free[0]:
            apterm  = sum( (b.bg_counts[model_index] for b in bands) )
            pixterm = sum( ( (b.bg_pix_counts[:,model_index]*b.pix_weights).sum() for b in bands if b.has_pixels) )
            return [(phase_factor*apterm - pixterm)/10**sm.p[0]]

        # general case -- this is a little gross, improve if time
        gradient = [0]*nfp
        for myband,band in zip(self.bands,bands):
            pts = sm.gradient(myband.bg_points)
            if nfp == 1: pts = N.asarray([pts])
            cp = 0
            for j in xrange(np):
                if not sm.free[j]: continue
                apterm = phase_factor*(myband.ap_evals * pts[j,:]).sum()
                if band.has_pixels:
                    pixterm = (band.pix_weights*(myband.pi_evals * pts[j,:]).sum(axis=1)).sum()
                gradient[cp] += apterm - pixterm
                cp += 1
        return gradient


###====================================================================================================###
class ROIDiffuseModel_PC(ROIDiffuseModel_OTF):
    """ The diffuse model is assumed to be pre-convolved.  This class then
        manages the exposure integral and model evaluation."""

    def init(self):
        tolerance = 0.02

    def setup(self):
        SkyIntegrator.set_tolerance(self.tolerance)
        exp = self.sa.exposure.exposure; psf = self.sa.psf
        if len(self.dmodel) == 1:
            self.bgs  = [Background(self.dmodel[0],exp[0],exp[1])]
        else:
            self.bgs  = map(Background,self.dmodel,exp)

    def set_state(self,energy,conversion_type):
        self.active_bg = self.bgs[conversion_type if (len(self.bgs) > 1) else 0]
        self.active_bg.setEnergy(energy)

    def __ap_value(self,center,radius):
        return SkyIntegrator.ss_average(self.active_bg,center,radius)

    def __pix_value(self,pixlist):
        return N.asarray(self.active_bg.wsdl_vector_value(pixlist))


