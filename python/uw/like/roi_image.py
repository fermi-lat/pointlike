"""
This code provides obejcts to generate counts and model counts images with
two dimensional projections. The code interfaces with the ROIAnalysis
object for making the calculations, skymaps.SkyImage object for storing
the data, and the image.ZEA object for plotting.  The high level object
roi_plotting.ROIDisplay can use to access these objects form a high
level plotting interface.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_image.py,v 1.41 2012/07/12 19:28:44 lande Exp $

author: Joshua Lande
"""
from skymaps import SkyImage,SkyDir,PythonUtilities,Band,WeightedSkyDirList
import numpy as np
import astropy.io.fits as pyfits
import scipy
import scipy.ndimage


from . pypsf import PretendBand
from . roi_diffuse import ROIDiffuseModel_OTF
from . roi_extended import ROIExtendedModel,ROIExtendedModelAnalytic
from . SpatialModels import RadiallySymmetricModel
from . pointspec_helpers import get_default_diffuse_mapper
from . roi_tsmap import TSCalc,TSCalcPySkyFunction
from uw.utilities import keyword_options
from uw.utilities.fitstools import get_fields
from uw.utilities.decorators import memoize
from pypsf import PsfOverlap
import collections
from abc import abstractmethod
import numbers



class ROIImage(object):
    """ This object is suitable for creating a SkyImage object
        and filling it with some physically meaningful
        quantity gotten from an ROIAnalysis object. 
        The acutal work is done by subclasses. """

    defaults = (
        ('size',     2,     'size of image in degrees'), 
        ('pixelsize',0.1,   'size, in degrees, of pixels'), 
        ('galactic', False, 'galactic or equatorial coordinates'), 
        ('proj',     'ZEA', 'projection name: can change if desired'),
        ('center',    None, 'Center of image. If None, use roi center.'),
        ('conv_type',   -1, 'Conversion type'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):
        """ Note, unlike ZEA, can support non-square images. To specify a nonsquare
            image, set the size parameter to a lenght two tuple:
                
                size=(10,5) # dx=10 degrees, dy=5 degrees. 
        """
        keyword_options.process(self, kwargs)

        if self.size < self.pixelsize:
            raise Exception("Can only create images with >=1 pixel in them.")
        
        self.roi = roi

        self.selected_bands = tuple(self.roi.bands if self.conv_type < 0 else \
            [ band for band in self.roi.bands if band.ct == self.conv_type ])

        # by default, use get energy range and image center from roi.
        if self.center is None: self.center=self.roi.roi_dir

        # set up, then create a SkyImage object to perform the projection
        # to a grid and manage an image
        if not isinstance(self.size,collections.Iterable):

            # make sure size and pixelsize are commensurate (helpful for
            # various downsampling code later).
            self.size = int(self.size/self.pixelsize + 0.01)*self.pixelsize
            self.skyimage = SkyImage(self.center, '', self.pixelsize, 
                                     self.size, 1, self.proj, self.galactic, False)
        else:
            self.skyimage = SkyImage(self.center, '', self.pixelsize, 
                                     float(self.size[0]), 1, self.proj, self.galactic, False, float(self.size[1]))

        self.fill()

        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()
        self.image=ROIImage.skyimage2numpy(self.skyimage)

    @staticmethod
    def skyimage2numpy(skyimage):
        nx, ny = skyimage.naxis1(), skyimage.naxis2()
        image=np.array(skyimage.image()).reshape((ny, nx))
        return image

    @abstractmethod
    def fill(self): pass

    def get_ZEA(self,axes=None,nticks=None):
        """ axes and nticks can be created by this object's constructor, but are
            more logically specified here. If they are not specified, get values from
            initial object creation. """
        # get out of the object all parameters which should be passed to ZEA.

        if hasattr(self.size,'__iter__'):
            raise Exception("Can only create ZEA object for square objects.")

        zea_dict = dict((d[0],self.__dict__[d[0]]) for d in ZEA.defaults if hasattr(d,'__iter__') and \
                hasattr(self,d[0]))
        if axes is not None: zea_dict['axes']=axes
        if nticks is not None: zea_dict['nticks']=nticks

        from uw.utilities.image import ZEA
        zea=ZEA(self.center,**zea_dict)
        zea.skyimage = self.skyimage
        # recalculate, in case the sky image has changed
        zea.image = ROIImage.skyimage2numpy(self.skyimage)

        # The old one gets removed by python's garbage collector (when zea.skyimage is replaced).
        zea.projector = zea.skyimage.projector() 

        zea.vmin,zea.vmax = zea.skyimage.minimum(), zea.skyimage.maximum()
        return zea

    def get_pyfits(self):
        """ Create and return a pyfits object that corresponds to the ROIImage object. 
            The fits file created is supposed to be consistent with the internal
            representation that SkyImage/SkyProj uses. """

        if self.galactic: 
            ctype1="GLON-%s" % self.proj
            ctype2="GLAT-%s" % self.proj
            # for some reason, SkyDir(0,0,SkyDir.GALACTIC).l() = 360
            crval1,crval2=self.center.l() % 360,self.center.b()
        else:
            ctype1="RA-%s" % self.proj
            ctype2="DEC-%s" % self.proj
            crval1,crval2=self.center.ra(),self.center.dec()

        cdelt1,cdelt2=-self.pixelsize,self.pixelsize

        # from SkyImage.cxx like 92:
        #   "center pixel; WCS convention is that center of a pixel is a half-integer"
        crpix1,crpix2=(self.skyimage.naxis1()+1)/2.0,(self.skyimage.naxis2()+1)/2.0

        values = [
            ["TELESCOP", "GLAST"],
            ["INSTRUME", "LAT"],
            ["DATE-OBS", ""],
            ["DATE-END", ""],
            ["EQUINOX", 2000.0, "Equinox of RA & DEC specifications"],
            ["CTYPE1", ctype1, "[RA|GLON]---%%%, %%% represents the projection method such as AIT"],
            ["CRPIX1", crpix1, "Reference pixel"],
            ["CRVAL1", crval1, "RA or GLON at the reference pixel"],
            ["CDELT1", cdelt1, "X-axis incr per pixel of physical coord at position of ref pixel(deg)"],
            ["CTYPE2", ctype2, "[DEC|GLAT]---%%%, %%% represents the projection method such as AIT"],
            ["CRPIX2", crpix2, "Reference pixel"],
            ["CRVAL2", crval2, "DEC or GLAT at the reference pixel"],
            ["CDELT2", cdelt2, "Y-axis incr per pixel of physical coord at position of ref pixel(deg)"],
            ["CROTA2",  0, "Image rotation (deg)"],
        ]
        for i in values: 
            if len(i)>2 and len(i[2])>47: i[2]=i[2][0:47]

        cards = [ pyfits.Card(*i) for i in values]

        header=pyfits.Header(cards=cards)

        hdu=pyfits.PrimaryHDU(data=self.image, header=header)
        fits = pyfits.HDUList([hdu])

        return fits

class ROITSMapImage(ROIImage):
    """ Subclass of ROIImage representing a residual TS map. """

    defaults = ROIImage.defaults + TSCalc.defaults

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        super(ROITSMapImage,self).__init__(*args,**kwargs)

    def fill(self):

        tscalc = TSCalc(self.roi,**keyword_options.defaults_to_kwargs(self,TSCalc))
        temp=TSCalcPySkyFunction(tscalc)
        self.skyimage.fill(temp.get_pyskyfun())

class CountsImage(ROIImage):
    """ This ROIImage subclass fills the sky image with the observed Fermi counts. """


    @staticmethod
    @memoize
    def process_filedata(roi,selected_bands):
        """ The radius parameter will apply a radius cut. 
            Assume that all bands are contiguous (hope this is always true!) """

        radius=roi.sa.maxROI

        emin = min(b.emin for b in selected_bands)
        emax = max(b.emax for b in selected_bands)

        ft1files=roi.sa.pixeldata.ft1files

        cuts = ['ENERGY > %s'% emin,
                'ENERGY < %s'% emax,
                'ZENITH_ANGLE < %s' % roi.sa.pixeldata.zenithcut,
                'THETA < %s' % roi.sa.pixeldata.thetacut,
                'EVENT_CLASS >= %s' % roi.sa.pixeldata.event_class]

        data = get_fields(ft1files,['RA','DEC','TIME','ENERGY','CONVERSION_TYPE'],cuts)
        # convert into skydirs
        skydirs = [ SkyDir(float(data['RA'][i]),float(data['DEC'][i])) for i in xrange(len(data['RA']))]

        # apply the same gti cut used to read in the initial WSDL.
        gti=roi.sa.pixeldata.gti
        good_dirs = []

        front_bins = [b for b in selected_bands if b.ct == 0]
        front_emin = min(b.emin for b in front_bins) if len(front_bins)>0 else None
        front_emax = max(b.emax for b in front_bins) if len(front_bins)>0 else None

        back_bins = [b for b in selected_bands if b.ct == 1]
        back_emin = min(b.emin for b in back_bins) if len(back_bins)>0 else None
        back_emax = max(b.emax for b in back_bins) if len(back_bins)>0 else None

        good_photons = []
        for skydir,time,energy,ct in zip(skydirs,data['TIME'],data['ENERGY'],data['CONVERSION_TYPE']):
            if gti.accept(time) and np.degrees(skydir.difference(roi.roi_dir)) < radius:
                if ct == 0 and \
                   front_emin is not None and front_emax is not None and \
                   energy > front_emin and energy < front_emax:
                    good_photons.append(skydir)
                if ct == 1 and \
                   back_emin is not None and back_emax is not None and \
                   energy > back_emin and energy < back_emax:
                    good_photons.append(skydir)

        return good_photons

    def fill(self):
        dirs = CountsImage.process_filedata(self.roi,self.selected_bands)

        for photon_dir in dirs:
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

    defaults = ROIImage.defaults + (
            ('override_point_sources', None, """ If either is specified, use override_point_sources these
                                                 and override_diffuse_sources to generate the image instead
                                                 of the sources in the ROI."""),
            ('override_diffuse_sources', None, 'Same as override_point_sources'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        if kwargs.has_key('proj') and kwargs['proj'] != 'ZEA':
            print "Warning, it is strongly advised to use the 'ZEA projection when creating model counts maps."

        super(ModelImage,self).__init__(*args,**kwargs)

    def fill(self):
        self.wsdl = self.skyimage.get_wsdl()

        self.solid_angle = np.radians(self.pixelsize)**2

        model_counts = np.zeros(len(self.wsdl),dtype=float)

        model_counts += self.all_point_source_counts()
        model_counts += self.all_diffuse_sources_counts()
        model_counts *= self.roi.phase_factor # don't forget about the phase factor!
        #NB -- this will need to be fixed if want to account for bracketing IRFs

        PythonUtilities.set_wsdl_weights(model_counts,self.wsdl)
        
        self.skyimage.set_wsdl(self.wsdl)

    @staticmethod
    def downsample(myarr,factor):
        """
        Code taken from http://code.google.com/p/agpy/source/browse/trunk/agpy/downsample.py

        Downsample a 1D or 2D array by averaging over *factor* pixels in each axis.
        Crops upper edge if the shape is not a multiple of factor.

        This code is pure numpy and should be fast.
        """
        assert isinstance(factor,numbers.Integral)
        assert len(myarr.shape) <= 2

        if len(myarr.shape) == 1:
            xs = myarr.shape[0]
            assert xs % factor == 0
            dsarr = np.concatenate([[myarr[i::factor]]
                                   for i in range(factor)]).mean(axis=0)
            return dsarr

        elif len(myarr.shape) == 2:
            xs,ys = myarr.shape
            assert xs % factor == 0 and ys % factor == 0
            dsarr = np.concatenate([[myarr[i::factor,j::factor] 
                for i in range(factor)] 
                for j in range(factor)]).mean(axis=0)
            return dsarr

    def bigger_wsdl(self,band,compare=None):
        """ Want to sample on a grid that is comparable in size (or
            smaller than) 10% of the psf to ensure we get a reasonable
            sampling of the grid. """
        if compare is None:
            r10=band.psf.inverse_integral_on_axis(0.10)
            compare=r10
        
        self.factor = int(np.ceil(self.pixelsize/compare))

        if self.factor == 1:
            return self.wsdl
        else:
            # hold onto this thing since it is needed by downsample_model
            if not hasattr(self.size,'__iter__'):
                self.fine_skyimage = SkyImage(self.center, '', float(self.pixelsize)/self.factor,
                                         self.size, 1, self.proj, self.galactic, False)
            else:
                self.fine_skyimage = SkyImage(self.center, '', float(self.pixelsize)/self.factor,
                                         float(self.size[0]), 1, self.proj, self.galactic, False, float(self.size[1]))

            wsdl = self.fine_skyimage.get_wsdl() 
            return wsdl

    def downsample_model(self,rvals):
        if self.factor==1:
            return rvals
        else:
            rvals = rvals.reshape((self.fine_skyimage.naxis2(), self.fine_skyimage.naxis1()))
            rvals = ModelImage.downsample(rvals,self.factor).flatten()
            return rvals

    @staticmethod
    def get_point_sources(roi,override_point_sources,override_diffuse_sources):
        if override_point_sources is None and override_diffuse_sources is None:
            return roi.psm.point_sources 
        if override_point_sources is None:
            return []
        elif not isinstance(override_point_sources,collections.Iterable):
            return [override_point_sources]
        else:
            return override_point_sources

    def all_point_source_counts(self):
        """ Calculate the point source contributions. """
        point_sources=ModelImage.get_point_sources(self.roi,self.override_point_sources,self.override_diffuse_sources)
        if len(point_sources)==0: return 0

        point_counts = np.zeros(len(self.wsdl),dtype=float)

        for band in self.selected_bands:
            cpsf = band.psf.cpsf

            # generate a list of skydirs on a finer grid.
            wsdl = self.bigger_wsdl(band)

            rvals  = np.empty(len(wsdl),dtype=float)

            for nps,ps in enumerate(point_sources):
                # evaluate the PSF at the center of each pixel
                cpsf.wsdl_val(rvals,ps.skydir,wsdl)

                # average the finer grid back to original resolution.
                temp = self.downsample_model(rvals)

                temp *= self.solid_angle #multiply by pixel solid angle
                temp *= band.expected(ps.model) # scale by total expected counts
                point_counts += temp

        return point_counts

    def extended_source_counts(self,extended_model):

        rd = self.roi.roi_dir 

        es = extended_model.extended_source
        sm = es.smodel

        extended_counts = np.zeros(len(self.wsdl),dtype=float)

        for band in self.selected_bands:
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

        mo=bg.smodel

        background_counts = np.zeros(len(self.wsdl),dtype=float)

        for band in self.selected_bands:

            ns,bg_points,bg_vector = ROIDiffuseModel_OTF.sub_energy_binning(band,bg.nsimps)

            pi_evals  = np.empty([len(self.wsdl),ns + 1])

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

    @staticmethod
    def get_diffuse_sources(roi,override_point_sources,override_diffuse_sources):

        if override_point_sources is None and override_diffuse_sources is None:
            return roi.dsm.bgmodels
        else:
            mapper=get_default_diffuse_mapper(roi.sa,roi.roi_dir,roi.quiet)
            if override_diffuse_sources is None:
                return []
            elif not isinstance(override_diffuse_sources,collections.Iterable):
                return [mapper(override_diffuse_sources)]
            else:
                return [mapper(ds) for ds in override_diffuse_sources]

    def all_diffuse_sources_counts(self):
        """ Calculate the diffuse source contributions. """
        
        bgmodels=ModelImage.get_diffuse_sources(self.roi,self.override_point_sources,self.override_diffuse_sources)

        return sum(self.diffuse_source_counts(bg)
                   for bg in bgmodels)


class ResidualImage(ROIImage):
    """ Has the same arguments at both ModelImage and Counts image but
        is a ROIImage object representing the difference between the
        source Counts and the ROI Model. """

    # get unique items in the Model + Counts defaults
    defaults = tuple(set(ModelImage.defaults).union(CountsImage.defaults))

    def __init__(self,roi,**kwargs):
        super(ResidualImage,self).__init__(roi,**kwargs)

        self.model=ModelImage(self.roi,**keyword_options.defaults_to_kwargs(self,ModelImage))
        self.counts=CountsImage(self.roi,**keyword_options.defaults_to_kwargs(self,CountsImage))

        self.image = self.counts.image - self.model.image
        SmoothedImage.add_to_skyimage(self.skyimage,self.image)


class RadialImage(object):
    """ This object is similar to ROIImage but performs a radial
        integral around around a given direction in the sky. """


    defaults = (
            ('center',       None,            'Center of image'),
            ('size',            2, 'Size of image (in degrees)'), 
            ('pixelsize', 0.00625, """ size of each image pixel. This is a little misleading because the
                                       size of each pixel varies, but regardless this is used to determine
                                       the total number of pixels with npix=size/pixelsize """),
            ('npix',         None, """ If specified, use this value instead of pixelsize. """),
            ('conv_type',    None,            'Conversion type'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):
        keyword_options.process(self, kwargs)

        if self.npix is None:
            self.npix = float(self.size)/self.pixelsize
        
        self.roi = roi

        self.selected_bands = tuple(self.roi.bands if self.conv_type < 0 else \
            [ band for band in self.roi.bands if band.ct == self.conv_type ])

        # by default, use get energy range and image center from roi.
        if self.center is None: self.center=self.roi.roi_dir

        # bins in theta^2
        self.bin_edges_deg = np.linspace(0.0,self.size**2,self.npix+1)
        self.bin_centers_deg = (self.bin_edges_deg[1:] + self.bin_edges_deg[:-1])/2.0
        
        # two factors of radians b/c theta^2
        self.bin_edges_rad = np.radians(np.radians(self.bin_edges_deg))
        self.bin_centers_rad = np.radians(np.radians(self.bin_centers_deg))

        # the lower and upper agle for each bin.
        self.theta_pairs_rad = zip(np.sqrt(self.bin_edges_rad[:-1]),
                                   np.sqrt(self.bin_edges_rad[1:]))


        self.fill()

    @abstractmethod
    def fill(self):
        """ This should fill up self.image appropriatly."""
        pass


class RadialCounts(RadialImage):
    """ This subclass of RadialImage calculates the counts within
        each radial bin. """

    defaults = RadialImage.defaults

    def fill(self):
        dirs = CountsImage.process_filedata(self.roi,self.selected_bands)
        diffs = [self.center.difference(i) for i in dirs]
        self.image=np.histogram(diffs,bins=np.sqrt(self.bin_edges_rad))[0]


class RadialSource(RadialImage):
    """ Subclass of RadialImage where returns the
        model predicted counts integrated radially
        for an extended source. Unlike RadialModel,
        the PSF is not convolved with the RadialSource
        (useful for comparing the extended source shape
        with the PSF). """

    defaults = RadialImage.defaults + (
            ('extended_source', None, 'A spatial moel'),
    )

    def fill(self):
        if self.extended_source is None:
            raise Exception("RadialSource must be given a spatial_model.")

        es=self.extended_source
        sm = es.model

        if not isinstance(es.spatial_model,RadiallySymmetricModel):
            raise Exception("spatial_model must be an instance of RadiallySymmetricModel.")

        total_counts = sum(band.expected(sm) for band in self.roi.bands)
        solid_angle = RadialModel.solid_angle_cone(np.radians(self.size))/self.npix

        fraction = es.spatial_model.at_r(np.sqrt(self.bin_centers_rad))
        fraction*=solid_angle

        self.image = total_counts*fraction



class RadialModel(RadialImage):

    defaults = RadialImage.defaults + (
            ('override_point_sources', None, """ If either is specified, use override_point_sources these
                                                 and override_diffuse_sources to generate the image instead
                                                 of the sources in the ROI."""),
            ('override_diffuse_sources', None, 'Same as override_point_sources'),
    )

    def fill(self):

        # Create fake bands just big enough to enclose the radial model.
        # This will speed up the convolution. 
        # The fake band will cause the diffuse sources be be convolved
        # in a smaller area, which will thus make them more accurate,
        # and create better plots.
        self.smaller_bands = []
        for band in self.selected_bands:
            rad=self.center.difference(self.roi.roi_dir) + np.radians(self.size)
            pb = PretendBand(band.e,band.ct, psf=band.psf, radius_in_rad=rad,
                             sd=self.center, emin=band.emin, emax=band.emax)
            self.smaller_bands.append(pb)

        self.image = np.zeros_like(self.bin_centers_rad)
        self.image += self.all_point_source_counts()
        self.image += self.all_diffuse_sources_counts()
        self.image *= self.roi.phase_factor # don't forget about the phase factor!

    def all_point_source_counts(self):
        """ Calculate the point source contributions. """
        point_sources=ModelImage.get_point_sources(self.roi,self.override_point_sources,self.override_diffuse_sources)
        if len(point_sources)==0: return 0

        point_counts = np.zeros_like(self.bin_centers_rad)

        overlap = PsfOverlap()

        for i,(theta_min,theta_max) in enumerate(self.theta_pairs_rad):

            model_counts=0

            for j,ps in enumerate(point_sources):
                for band in self.selected_bands:

                    # this code requires a redundant call to overlap. Improve if time.
                    # Note that ragged_edge is not appropriate here because our counts
                    # are always sumed in the exact range.
                    fraction=overlap(band,self.center,ps.skydir,radius_in_rad=theta_max, ragged_edge=np.inf) - \
                             overlap(band,self.center,ps.skydir,radius_in_rad=theta_min, ragged_edge=np.inf)

                    model_counts += band.expected(ps.model)*fraction

            point_counts[i] = model_counts

        return point_counts

    def extended_source_counts(self,extended_model):
        if type(extended_model) not in [ROIExtendedModel,ROIExtendedModelAnalytic]:
            raise Exception("Unknown extended model.")

        roi=self.roi
        sm = extended_model.extended_source.model

        extended_counts = np.zeros_like(self.bin_centers_rad)

        for band,smaller_band in zip(self.selected_bands,self.smaller_bands):

            extended_model.set_state(smaller_band)

            if type(extended_model) == ROIExtendedModel:

                nside = RadialModel.get_nside(self.size,self.npix)

                temp_band = Band(nside)
                wsdl = WeightedSkyDirList(temp_band,self.center,np.radians(self.size),True)
                vals=extended_model._pix_value(wsdl)

                rvals=np.empty(len(wsdl),dtype=float)
                PythonUtilities.arclength(rvals,wsdl,self.center)

                # get average value in each ring by averaging values.
                fraction = np.histogram(rvals,weights=vals,bins=np.sqrt(self.bin_edges_rad))[0]/\
                           np.histogram(rvals,bins=np.sqrt(self.bin_edges_rad))[0]

                # multiply intensities by solid angle in ring
                fraction *= RadialModel.solid_angle_cone(np.radians(self.size))/self.npix

            elif type(extended_model) == ROIExtendedModelAnalytic:

                fraction = np.empty_like(self.bin_centers_rad)

                for i,(theta_min,theta_max) in enumerate(self.theta_pairs_rad):

                    fraction[i]=extended_model._overlaps(self.center,band,theta_max) - \
                             extended_model._overlaps(self.center,band,theta_min)

            # total counts from source * fraction of PDF in ring = model predictions in each ring.
            extended_counts += band.expected(sm)*fraction

        return extended_counts

    @staticmethod
    def solid_angle_cone(radius_in_radians):
        return 2*np.pi*(1-np.cos(radius_in_radians))

    @staticmethod
    def get_nside(size,npix,num_points_per_ring=200):
        """ Solid angle of each healpix pixel is 4pi/(12*ns^2)
            Solid angel of each ring is pi*(size)^2/npix
            Want size of each ring > num_points_per_ring*size of each healpix (so that
            we get 20 pixels to sample each ring).

            Solving for ns, the size of the healpixels required, we get the required formula
        """
        total_solid_angle=4*np.pi
        image_solid_angle=RadialModel.solid_angle_cone(np.radians(size))

        solid_angle_per_ring=image_solid_angle/npix

        nside = int(np.ceil(np.sqrt(num_points_per_ring*total_solid_angle/(12*solid_angle_per_ring))))
        return nside

    def otf_source_counts(self,bg):

        roi=self.roi

        mo=bg.smodel

        background_counts = np.zeros_like(self.bin_centers_rad)

        for band,smaller_band in zip(self.selected_bands,self.smaller_bands):

            ns,bg_points,bg_vector = ROIDiffuseModel_OTF.sub_energy_binning(band,bg.nsimps)

            nside = RadialModel.get_nside(self.size,self.npix)

            temp_band = Band(nside)
            wsdl = WeightedSkyDirList(temp_band,self.center,np.radians(self.size),True)

            ap_evals = np.empty([len(self.bin_centers_rad),len(bg_points)])

            for ne,e in enumerate(bg_points):

                bg.set_state(e,band.ct,smaller_band)

                rvals=np.empty(len(wsdl),dtype=float)
                PythonUtilities.arclength(rvals,wsdl,self.center)
                vals=bg._pix_value(wsdl)

                # get average value in each ring by averaging values.
                ap_evals[:,ne] = np.histogram(rvals,weights=vals,bins=np.sqrt(self.bin_edges_rad))[0]/\
                                 np.histogram(rvals,bins=np.sqrt(self.bin_edges_rad))[0]

            # multiply intensities by solid angle in ring
            ap_evals *= RadialModel.solid_angle_cone(np.radians(self.size))/self.npix

            ap_evals          *= bg_vector
            mo_evals           = mo(bg_points)
            background_counts += (ap_evals * mo_evals).sum(axis=1)

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
        bgmodels=ModelImage.get_diffuse_sources(self.roi,self.override_point_sources,self.override_diffuse_sources)

        return sum(self.diffuse_source_counts(bg)
                   for bg in bgmodels)


class SmoothedImage(ROIImage):

    """ Represnts a ROIImage objects that is smoothed. This
        is an abstract base class, but subclasses represent
        smoothed versions of particular ROIImage objects.

        SmoothedCounts -> CountsImage
        SmoothedModel -> ModelImage
        SmoothedResidual -> ResidualImage

        Note, after smoothing, the image is normalized to represent
        
    """
    
    smoothed_options = (
        ('kerneltype',      'tophat', 'Type of kernel to use'),
        ('kernel_rad',          0.25, 'Sum counts within radius degrees.'),
        ('per_solid_angle',    False, 'If true, after smoothing divide by solid angle (counts/[deg]^2'),
    )

    defaults = ROIImage.defaults + smoothed_options



    @staticmethod
    def convolve_array(image,width):
        """ Summ all this pixels within a given pixel width
        
            code taken from 
            
            http://code.google.com/p/agpy/source/browse/trunk/agpy/convolve.py """

        return out

    @staticmethod
    def add_to_skyimage(skyimage,image):
        """ Take in a two dimensional numpy array representing the
            fits data and put it into the skyimage. """
        temp=image.reshape(image.shape[0]*image.shape[1])

        wsdl = skyimage.get_wsdl()
        PythonUtilities.set_wsdl_weights(temp,wsdl)
        skyimage.set_wsdl(wsdl)

    @staticmethod
    def convolve_skyimage(skyimage,kernel):
        """ Take in a skyimage object and sum each pixel with all the
            pixels within a given width. """
        image=ROIImage.skyimage2numpy(skyimage)

        convolved=scipy.ndimage.filters.convolve(image, kernel)

        SmoothedImage.add_to_skyimage(skyimage,convolved)

        return skyimage

    @staticmethod
    def get_kernel(kerneltype,kernel_rad,pixelsize):
        """ returns the kernel. Code inspired by
                http://code.google.com/p/agpy/source/browse/trunk/agpy/convolve.py
        """

        width = int(np.ceil(kernel_rad/pixelsize))

        if kerneltype == 'tophat':
            kernelsize=4*width

            kernel = np.zeros([kernelsize,kernelsize],dtype=float)

            xx,yy = np.indices(kernel.shape)
            rr = np.sqrt((xx-kernel.shape[0]/2.)**2+(yy-kernel.shape[1]/2.)**2)
            kernel[rr<=width] = 1

        elif kerneltype == 'gaussian':
            kernelsize = 8*width

            kernel = np.zeros([kernelsize,kernelsize],dtype=float)
            xx,yy = np.indices(kernel.shape)
            rr = np.sqrt((xx-kernel.shape[0]/2.)**2+(yy-kernel.shape[1]/2.)**2)
            kernel = np.exp(-(rr**2)/(2*width**2))

            # Gaussian kernels should be normalized so they don't change overall normalization
            kernel /= kernel.sum()

        else:
            raise Exception("...")

        return kernel


    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):

        super(SmoothedImage,self).__init__(*args,**kwargs)

        self.kernel=self.get_kernel(self.kerneltype,self.kernel_rad,self.pixelsize)

        self.kernelsize=self.kernel.shape[0]

        # Make an image, bigger then the desired one, by twice
        # the smooth radius in each direction.
        self.pass_dict=keyword_options.defaults_to_kwargs(self,self.object)
        self.pass_dict['size']=self.size+self.kernelsize*self.pixelsize
        self.smoothed=self.object(self.roi,**self.pass_dict)


        self.smoothed.skyimage=SmoothedImage.convolve_skyimage(self.smoothed.skyimage,self.kernel)
        self.smoothed.image=ROIImage.skyimage2numpy(self.smoothed.skyimage)

        # now, shrink down smoothed image and replace the current skyimage with it

        self.image = self.smoothed.image[self.kernelsize/2:-self.kernelsize/2,self.kernelsize/2:-self.kernelsize/2]

        if self.per_solid_angle:
            # convert from counts to counts per square degree
            self.image = self.image/self.pixelsize**2

        SmoothedImage.add_to_skyimage(self.skyimage,self.image)

class SmoothedCounts(SmoothedImage):

    defaults = CountsImage.defaults + SmoothedImage.smoothed_options

    object=CountsImage

class SmoothedModel(SmoothedImage):

    defaults = ModelImage.defaults + SmoothedImage.smoothed_options

    object=ModelImage

class SmoothedResidual(SmoothedImage):

    defaults = ResidualImage.defaults + SmoothedImage.smoothed_options

    object=ResidualImage
