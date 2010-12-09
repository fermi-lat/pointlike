"""
This code provides obejcts to generate counts and model counts images with
two dimensional projections. The code interfaces with the ROIAnalysis
object for making the calculations, skymaps.SkyImage object for storing
the data, and the image.ZEA object for plotting.  The high level object
roi_plotting.ROIDisplay can use to access these objects form a high
level plotting interface.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_image.py,v 1.1 2010/12/05 09:53:32 lande Exp $

author: Joshua Lande
"""
from skymaps import SkyImage,SkyDir,PythonUtilities
import numpy as N

from uw.like.roi_diffuse import ROIDiffuseModel_OTF
from uw.like.roi_extended import ROIExtendedModel
from uw.utilities import keyword_options
from uw.utilities.fitstools import get_fields
from uw.utilities.image import ZEA


import numpy

class ROIImage(object):
    """ This object is suitable for creating a SkyImage object
        and filling it with some physically meaningful
        quantity gotten from an ROIAnalysis object. 
        The acutal work is done by subclasses. """

    defaults = ZEA.defaults 
    
    defaults += (
            ('center',    None,  'Center of image'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):
        keyword_options.process(self, kwargs)
        
        self.roi = roi

        # by default, use get energy range and image center from roi.
        if self.center is None: self.center=self.roi.roi_dir

        # set up, then create a SkyImage object to perform the projection
        # to a grid and manage an image
        self.skyimage = SkyImage(self.center, self.fitsfile, self.pixelsize, 
                                 self.size, 1, self.proj, self.galactic, False)

        self.fill()

        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()
        self.image=N.array(self.skyimage.image()).reshape((self.ny, self.nx))

    def get_ZEA(self,axes=None,nticks=None):
        """ axes and nticks can be created by this object's constructor, but are
            more logically specified here. If they are not specified, get values from
            initial object creation. """
        # get out of the object all parameters which should be passed to ZEA.
        zea_dict = dict((d[0],self.__dict__[d[0]]) for d in ZEA.defaults if hasattr(d,'__iter__'))
        if axes is not None: zea_dict['axes']=axes
        if nticks is not None: zea_dict['nticks']=nticks

        zea=ZEA(self.center,**zea_dict)
        zea.skyimage = self.skyimage
        zea.image = self.image

        # The old one gets removed by python's garbage collector (when zea.skyimage is replaced).
        zea.projector = zea.skyimage.projector() 

        zea.vmin,zea.vmax = zea.skyimage.minimum(), zea.skyimage.maximum()
        return zea


class CountsImage(ROIImage):

    defaults = ROIImage.defaults + (
        ('mc_src_id', None, ''),
        ('cuts',      None, '')
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):

        self.ft1files=roi.sa.pixeldata.ft1files

        super(CountsImage,self).__init__(roi,**kwargs)

    def process_filedata(self,fields):

        emin = self.roi.bin_edges[0]
        emax = self.roi.bin_edges[-1]
        conv_type = self.roi.sa.conv_type

        base_cuts = ['ENERGY > %s'%(emin),'ENERGY < %s'%(emax),'ZENITH_ANGLE < 105']
        if conv_type >= 0:        base_cuts += ['EVENT_CLASS == %d'%(conv_type)]
        if self.mc_src_id is not None: base_cuts += ['MC_SRC_ID == %d'%(self.mc_src_id)]
        cuts = base_cuts if self.cuts is None else self.cuts + base_cuts

        data = get_fields(self.ft1files,fields,cuts)

        return data

    def fill(self):
        print 'Currently, this function is not correctly applying Theta/GTI/Zenith cuts.'
        data = self.process_filedata(['RA','DEC'])

        for i in xrange(len(data['RA'])):
            photon_dir=SkyDir(float(data['RA'][i]),float(data['DEC'][i]))
            self.skyimage.addPoint(photon_dir)


class ModelImage(ROIImage):
    """ This ROIImage subclass fills the sky image with the model
        predicted counts for a fermi sky model described by an ROIAnalysis
        object.

        This code is forced to deal with the fact that model intensity
        can vary significantly across a spatial pixel. The rest of the
        pointlike code can avoid this whole issue by scaling the healpix
        pixel size with the PSF to ensure that pixels are always small
        compared to the instrument's intrisic resolution. But since
        model predicted counts maps can be generated of arbitary pixel size,
        this issue must directly be dealt with.

        The solution that this code uses to deal with this issue is to
        simply sample from a grid finer by an integer number of pixels
        in each dimensions.  After calculating the model predictions,
        the nearby blocks of model predictions are averaged to downsample
        to create the model predictions.  This formulation assumes that
        each of the subpixels has the same solid angle and so it only
        suitable for relativly small images where pixels have equal area.
        For that reason, it is advised to use the ZEA projection.

        For point and extended sources, the characteristic scale with
        which the convolution must be small compared to is the PSF. So
        the formula for determining the factor is

        factor = ceil(pixelsize/r10)

        Where pxielsize is the plotting pixel size and r10 is the 10%
        containment radius of the PSF.

        For background sources, the characteristic scale is not the PSF
        but the convolution grid pixelsize. So the formula for determining
        the factor is instead
        
        factor = ceil(pixelsize/(conv_pixelsize/4))

        Where conv_pixelsize is the size of the convolution grid's pixels.
        
        For background sources, this algorithm is generally efficiency
        since we except the background to vary on this smaller scale all
        across the image.  But for point and (small) extended sources,
        this algorithm is generally very poor because it requires
        calculating the PSF (or PDF) at many points where the value
        is very close to 0. A better algorithm would be an adaptive
        quadrature integration algorithm which evaluate the integral in
        each pixel, then did a more accurate integral and iterated until
        the integral converged. This would avoid having to evaluate the
        model predictions for a source very finely far from the source.
        On the other hand, adding this feature (presumably to C++
        for optimization) would be very costly, and this code runs
        fast enough...

        """

    defaults = ROIImage.defaults
    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        if kwargs.has_key('proj') and kwargs['proj'] != 'ZEA':
            print "Warning, it is strongly advised to use the 'ZEA projection when creating model counts maps."

        super(ModelImage,self).__init__(*args,**kwargs)

    def fill(self):
        self.wsdl = self.skyimage.get_wsdl()

        self.solid_angle = N.radians(self.pixelsize)**2

        model_counts = N.zeros(len(self.wsdl),dtype=float)

        model_counts += self.all_point_sources_counts()
        model_counts += self.all_diffuse_sources_counts()

        PythonUtilities.set_wsdl_weights(model_counts,self.wsdl)
        
        self.skyimage.set_wsdl(self.wsdl)

    @staticmethod
    def downsample(myarr,factor):
        """
        Code taken from http://code.google.com/p/agpy/source/browse/trunk/agpy/downsample.py

        Downsample a 2D array by averaging over *factor* pixels in each axis.
        Crops upper edge if the shape is not a multiple of factor.

        This code is pure numpy and should be fast.
        """
        xs,ys = myarr.shape
        crarr = myarr[:xs-(xs % int(factor)),:ys-(ys % int(factor))]
        dsarr = numpy.concatenate([[crarr[i::factor,j::factor] 
            for i in range(factor)] 
            for j in range(factor)]).mean(axis=0)
        return dsarr

    def bigger_wsdl(self,band,compare=None):
        """ Want to sample on a grid that is comparable in size (or
            smaller than) 10% of the psf.
            to ensure we get a reasonable. """
        if compare is None:
            r10=band.psf.inverse_integral_on_axis(0.10)
            compare=r10
        
        self.factor = int(N.ceil(self.pixelsize/compare))

        if self.factor == 1:
            return self.wsdl
        else:
            # hold onto this thing since it is needed by downsample_model
            self.fine_skyimage = SkyImage(self.center, self.fitsfile, float(self.pixelsize)/self.factor,
                                     self.size, 1, self.proj, self.galactic, False)
            wsdl = self.fine_skyimage.get_wsdl() 
            return wsdl

    def downsample_model(self,rvals):
        if self.factor==1:
            return rvals
        else:
            rvals = rvals.reshape((self.fine_skyimage.naxis2(), self.fine_skyimage.naxis1()))
            rvals = ModelImage.downsample(rvals,self.factor).flatten()
            return rvals

    def all_point_sources_counts(self):
        """ Calculate the point source contributions. """
        roi=self.roi
        bands=roi.bands
        point_sources = roi.psm.point_sources
        point_counts = N.zeros(len(self.wsdl),dtype=float)

        for band in bands:
            cpsf = band.psf.cpsf

            # generate a list of skydirs on a finer grid.
            wsdl = self.bigger_wsdl(band)

            rvals  = N.empty(len(wsdl),dtype=float)

            for nps,ps in enumerate(point_sources):
                # evaluate the PSF at the center of each pixel
                cpsf.wsdl_val(rvals,ps.skydir,wsdl)

                # average the finer grid back to original resolution.
                temp = self.downsample_model(rvals)

                temp *= self.solid_angle #multiply by pixel solid angle
                temp *= band.ps_counts[nps] # scale by total expected counts
                point_counts += temp

        return point_counts

    def extended_source_counts(self,extended_model):

        rd = self.roi.roi_dir 

        es = extended_model.extended_source
        sm = es.smodel

        extended_counts = N.zeros(len(self.wsdl),dtype=float)

        for band in self.roi.bands:
            extended_model.set_state(band)
            exposure=band.exp.value

            er = exposure(es.spatial_model.center,extended_model.current_energy)/exposure(rd,extended_model.current_energy)
            es_counts = band.expected(sm)*er

            wsdl = self.bigger_wsdl(band)

            es_pix_counts = extended_model._pix_value(wsdl)*self.solid_angle

            es_pix_counts= self.downsample_model(es_pix_counts)

            bg_pix_counts = es_pix_counts * es_counts

            extended_counts += bg_pix_counts

        return extended_counts

    def otf_source_counts(self,bg):
        roi=self.roi
        bands=roi.bands

        mo=bg.smodel

        background_counts = N.zeros(len(self.wsdl),dtype=float)

        for band in bands:

            # use a higher nsimps at low energy where effective area is jagged
            ns = (2 if band.emin<200 else 1)*bg.nsimps
            bg_points = sp = N.logspace(N.log10(band.emin),N.log10(band.emax),ns+1)
            bg_vector = sp * (N.log(sp[-1]/sp[0])/(3.*ns)) * \
                                     N.asarray([1.] + ([4.,2.]*(ns/2))[:-1] + [1.])

            pi_evals  = N.empty([len(self.wsdl),ns + 1])

            wsdl = self.bigger_wsdl(band,compare=bg.pixelsize/4.0)

            for ne,e in enumerate(bg_points):
                bg.set_state(e,band.ct,band)
                temp = self.downsample_model(bg._pix_value(wsdl))
                pi_evals[:,ne] = temp

            pi_evals *= (self.solid_angle * bg_vector)
            mo_evals  = mo(bg_points)
            pi_counts = (pi_evals * mo_evals).sum(axis=1)

            background_counts += pi_counts

        return background_counts

    def diffuse_source_counts(self,bg):
        if isinstance(bg,ROIDiffuseModel_OTF):
            return self.otf_source_counts(bg)
        elif isinstance(bg,ROIExtendedModel):
            return self.extended_source_counts(bg)
        else:
            raise Exception("Unable to calculate model predictions for diffuse source %s", bg.name)

    def all_diffuse_sources_counts(self):
        """ Calculate the diffuse source contributions. """
        return sum(self.diffuse_source_counts(bg)
                   for bg in self.roi.dsm.bgmodels)
