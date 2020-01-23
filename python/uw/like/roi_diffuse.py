"""
Provides classes to encapsulate and manipulate diffuse sources.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_diffuse.py,v 1.32 2012/09/13 06:16:44 kerrm Exp $

author: Matthew Kerr
"""
import sys
import numpy as np
from uw.utilities.convolution import BackgroundConvolution
from uw.utilities import keyword_options
from skymaps import SkyIntegrator,Background,IsotropicSpectrum,IsotropicConstant,IsotropicPowerLaw,DiffuseFunction
import copy
import collections

class SmallBand(object):
    """ A little holder."""
    pass

class DiffuseSource(object):
    """ Associate a spatial model with a spectral scaling model."""
    __counter = 0

    def __init__(self,diffuse_model,scaling_model,name=None):

        self.dmodel = diffuse_model
        self.smodel = scaling_model

        self.smodel.background = True

        if name is None:
            self.name = 'Diffuse Source %d'%(DiffuseSource.__counter)
            DiffuseSource.__counter += 1
        else: self.name = name

        if not isinstance(self.dmodel,collections.Iterable):
            self.dmodel = [self.dmodel]
   
    def __str__(self): return '\n'.join((self.name,'\t'+self.dmodel.__str__(),
            '\t'+self.smodel.__str__()))

    def __getstate__(self):
        """ this is a poor man's pickeling. Convert all the IsotropicSpectrum and DiffuseFunction
            objects into pickelable filenames. This implementation should be more robust. 
            
            You should be able to pickle Diffuse sources:

                >>> import pickle
                >>> import os
                >>> from uw.utilities import path
                >>> from uw.like.Models import Constant, FileFunction
                >>> from skymaps import IsotropicPowerLaw

            Pickle isotropic:

                >>> iso = path.expand('$GLAST_EXT/diffuseModels/v2r0p1/isotrop_2year_P76_source_v1.txt')
                >>> ds = DiffuseSource(IsotropicSpectrum(iso),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))

            After pickling the object, the original object should be uncahgned:

                >>> print (type(ds.dmodel[0]))
                <class 'skymaps.IsotropicSpectrum'>
            
            Pickle galactic diffuse:

                >>> gal = path.expand('$GLAST_EXT/diffuseModels/v2r0p1/ring_2year_P76_v0.fits')
                >>> ds = DiffuseSource(DiffuseFunction(gal),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))
            
            Pickle Isotropic PowerLaw:

                >>> ds = DiffuseSource(IsotropicPowerLaw(),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))

            Pickle Isotropic Constant:

                >>> ds = DiffuseSource(IsotropicConstant(),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))
            
        """
        d=copy.copy(self.__dict__)

        def convert_spectrum(spectrum):
            if isinstance(spectrum,IsotropicSpectrum) or isinstance(spectrum,DiffuseFunction):
                return (type(spectrum),spectrum.name())
            elif isinstance(spectrum,IsotropicPowerLaw):
                return (type(spectrum),spectrum.flux(),spectrum.index())
            elif isinstance(spectrum,IsotropicConstant):
                return (type(spectrum),spectrum.constant())
            else:
                # unrecognized type
                return spectrum

        d['dmodel'] = map(convert_spectrum,d['dmodel'])
        return d

    def __setstate__(self,state):
        """ recreate the dmodel. """
        self.__dict__ = state
        def unconvert_spectrum(spectrum):
            if type(spectrum) == tuple:
                return spectrum[0](*spectrum[1:])
            else:
                return spectrum

        self.dmodel = map(unconvert_spectrum,self.dmodel)

    def copy(self):
        """ Make a copy of a diffuse source. 
        
            First, create the DS

                >>> from uw.like.Models import Constant
                >>> from skymaps import IsotropicPowerLaw
                >>> ds = DiffuseSource(IsotropicConstant(),Constant())

            Nowe, we can copy it:

                >>> ds_copy = ds.copy()
                >>> print (type(ds.dmodel[0]))
                <class 'skymaps.IsotropicConstant'>
                >>> print (type(ds_copy.dmodel[0]))
                <class 'skymaps.IsotropicConstant'>
        """
        return DiffuseSource(
            name = self.name,
            diffuse_model = self.dmodel,
            scaling_model = self.smodel.copy())

class ROIDiffuseModel(object):
    """ Associate a SkySpectrum with a spectral model.
    
        Provide the interface that ROIDiffuseModel et al. should satisfy."""

    defaults = (
        ('quiet', False, 'Set True to quiet'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,spectral_analysis,diffuse_source,roi_dir,name=None,**kwargs):

        self.sa = spectral_analysis
        self.roi_dir = roi_dir
        self.diffuse_source = diffuse_source
        self.dmodel         = diffuse_source.dmodel
        self.smodel         = diffuse_source.smodel
        self.name           = diffuse_source.name

        keyword_options.process(self, kwargs)
        self.setup()
   
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
        raise NotImplementedError('Classes must implement this method!')

    def update_counts(self,bands,model_index):
        """ This method *must* set the following members of each band:
            band.bg_counts[model_index] -- the total counts expected for
                    the model in the aperture
            band.bg_pix_counts[:,model_index] -- if the band has pixels,
                    an npixel vector with the expected counts from the
                    model for each data pixel
        """
        raise NotImplementedError('Classes must implement this method!')

    def gradient(self,bands,model_index):
        """ This method should return the gradient with respect to the
            free parameters of the spectral model.
        """
        raise NotImplementedError('Classes must implement this method!')


class ROIDiffuseModel_OTF(ROIDiffuseModel):
    """ Use an on-the-fly numerical scheme to convolve the model
        with the PSF.
        
        Use a Simpson's rule integration over the exposure.
        
        To use a different convolution scheme but maintain the Simpson's
        rule exposure integral, child methods may override _ap_value 
        and _pix_value."""

    defaults = ROIDiffuseModel.defaults + (
        ('pixelsize',0.25,'Pixel size for convolution grid'),
        ('npix',101,"Note -- can be overridden at the band level"),
        ('nsimps',4,"Note -- some energies use a multiple of this"),
        ('r_multi',1.0,"Multiple of r95 to set max dimension of grid"),
        ('r_max',20,"An absolute maximum (half)-size of grid (deg)"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        super(ROIDiffuseModel_OTF,self).__init__(*args,**kwargs)

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

    def set_state(self,energy,conversion_type,band,**kwargs):
        self.active_bgc = self.bgc[conversion_type if (len(self.bgc) > 1) else 0]
        multi = 1 + 0.01*(energy==band.emin) -0.01*(energy==band.emax)
        r95 = self.sa.psf.inverse_integral(energy*multi,conversion_type,95)
        rad = self.r_multi*r95 + np.degrees(band.radius_in_rad)
        rad = max(min(self.r_max,rad),np.degrees(band.radius_in_rad)+2.5)
        npix = int(round(2*rad/self.pixelsize))
        npix += (npix%2 == 0)
        self.active_bgc.setup_grid(npix,self.pixelsize)
        #self.active_bgc.do_convolution(energy,conversion_type,override_en=band.e)
        self.active_bgc.do_convolution(energy,conversion_type)

    @staticmethod 
    def sub_energy_binning(band,nsimps):
        """ This code is in a static method so it can be used consistently with 
            the model predicted counts code (in roi_image.py). """

        # use a higher nsimps at low energy where effective area is jagged
        ns = (2 if band.emin<200 else 1)*nsimps
        if ns > 0:
            bg_points = sp = np.logspace(np.log10(band.emin),np.log10(band.emax),ns+1)
            bg_vector = sp * (np.log(sp[-1]/sp[0])/(3.*ns)) * \
                                     np.asarray([1.] + ([4.,2.]*(ns/2))[:-1] + [1.])
        else:
            bg_points = sp = np.asarray([band.e])
            bg_vector = sp * np.log(band.emax/band.emin)
        return ns,bg_points,bg_vector

    def initialize_counts(self,bands,roi_dir=None):
        rd = self.roi_dir if roi_dir is None else roi_dir
        self.bands = [SmallBand() for i in xrange(len(bands))]

        for iband,(myband,band) in enumerate(zip(self.bands,bands)):
            if not self.quiet: 
                status_string = '...convolving band %2d/%2d'%(iband+1,len(self.bands))
                print (status_string,;sys.stdout.flush())

            ns,myband.bg_points,myband.bg_vector = ROIDiffuseModel_OTF.sub_energy_binning(band,self.nsimps)

            #figure out best way to handle no pixel cases...
            myband.ap_evals  = np.empty(ns + 1)      
            myband.pi_evals  = np.empty([len(band.wsdl),ns + 1]) if band.has_pixels else 0
                  
            for ne,e in enumerate(myband.bg_points):
                self.set_state(e,band.ct,band)
                myband.ap_evals[ne] = self._ap_value(rd,band.radius_in_rad)
                if band.has_pixels:
                    myband.pi_evals[:,ne] = self._pix_value(band.wsdl)

            myband.ap_evals *= (band.solid_angle   * myband.bg_vector)
            myband.pi_evals *= (band.b.pixelArea() * myband.bg_vector)

            # calculate integral counts for later
            myband.mo_evals  = self.smodel(myband.bg_points)
            myband.ap_counts = (myband.ap_evals * myband.mo_evals).sum()
            if band.has_pixels:
                myband.pi_counts = (myband.pi_evals * myband.mo_evals).sum(axis=1)
            if not self.quiet: 
                print ('\b'*(2+len(status_string)),;sys.stdout.flush())

        self.init_p = self.smodel.get_all_parameters(internal=True)
        self.init_norm = self.smodel[0]
        self.prev_p = self.smodel.get_all_parameters(internal=True) +1e-5 # kluge
        if not self.quiet: print()
        
    def update_counts(self,bands,model_index):

        mi = model_index
        sm = self.smodel
        #if N.all(self.prev_p == sm.p): return
        #self.prev_p[:] = sm.p
        smp = sm.get_all_parameters(internal=True)
        if np.all(self.prev_p == smp): return
        self.prev_p[:] = smp
        

        # counts can just be scaled from initial integral
        smp = sm.get_all_parameters(internal=True)
        if np.all(smp[1:] == self.init_p[1:]):
            ratio = sm[0]/self.init_norm
            for myband,band in zip(self.bands,bands): 
               band.bg_counts[mi] = ratio * myband.ap_counts
               if band.has_pixels:
                  band.bg_pix_counts[:,mi] = ratio * myband.pi_counts

        # update requires new integral over energy
        else:
            for myband,band in zip(self.bands,bands):
                pts = sm(myband.bg_points)
                band.bg_counts[mi] = (myband.ap_evals * pts).sum()
                if band.has_pixels:
                    band.bg_pix_counts[:,mi] = (myband.pi_evals * pts).sum(axis=1)               
                  
    def _ap_value(self,center,radius):
        """ Return the integral of the model over the aperture for the
            current state (energy/conversion type).
            Convolution could/should be done.
        """
        return self.active_bgc.ap_average(radius)

    def _pix_value(self,pixlist):
        """ Return the model evaluated at each data pixel for the current
            state (energy/conversion type).
            Convolution could/should be done.
        """
        return self.active_bgc(pixlist,self.active_bgc.cvals)

    def gradient(self,bands,model_index):
        sm  = self.smodel
        npar  = len(sm.get_all_parameters())
        nfp = sm.free.sum()

        # special case -- no free parameters
        if nfp == 0: return []

        # general case -- this is a little gross, improve if time
        gradient = [0]*nfp
        for myband,band in zip(self.bands,bands):
            pts = sm.gradient(myband.bg_points)
            cp = 0
            for j in xrange(npar):
                if not sm.free[j]: continue
                apterm = band.phase_factor*(myband.ap_evals * pts[j,:]).sum()
                if band.has_pixels:
                    pixterm = (band.pix_weights*(myband.pi_evals * pts[j,:]).sum(axis=1)).sum()
                else:
                    pixterm = 0
                gradient[cp] += apterm - pixterm
                cp += 1
        return gradient


class ROIDiffuseModel_PC(ROIDiffuseModel_OTF):
    """ The diffuse model is assumed to be pre-convolved.  This class then
        manages the exposure integral and model evaluation."""

    defaults = ROIDiffuseModel.defaults + (
        ('tolerance',0.02,'SkyIntegrator tolerance'),
        ('nsimps',4,'Number of subenergies to evalulate the simpson integral over'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        super(ROIDiffuseModel_OTF,self).__init__(*args,**kwargs)

    def setup(self):
        SkyIntegrator.set_tolerance(self.tolerance)
        exp = self.sa.exposure.exposure; psf = self.sa.psf
        if len(self.dmodel) == 1:
            self.bgs  = [Background(self.dmodel[0],exp[0],exp[1])]
        else:
            self.bgs  = map(Background,self.dmodel,exp)

    def set_state(self,energy,conversion_type,band,**kwargs):
    #def set_state(self,energy,conversion_type,**kwargs):
        self.active_bg = self.bgs[conversion_type if (len(self.bgs) > 1) else 0]
        self.active_bg.setEnergy(energy)

    def _ap_value(self,center,radius):
        return SkyIntegrator.ss_average(self.active_bg,center,radius)

    def _pix_value(self,pixlist):
        return np.asarray(self.active_bg.wsdl_vector_value(pixlist))


if __name__ == "__main__":
    import doctest
    doctest.testmod()
