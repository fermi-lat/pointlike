"""
This code provides obejcts to generate counts and model counts images with
two dimensional projections. The code interfaces with the ROIAnalysis
object for making the calculations, skymaps.SkyImage object for storing
the data, and the image.ZEA object for plotting.  The high level object
roi_plotting.ROIDisplay can use to access these objects form a high
level plotting interface.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_image.py,v 1.16 2011/04/05 22:46:06 lande Exp $

author: Joshua Lande
"""
from skymaps import SkyImage,SkyDir,PythonUtilities,Band,WeightedSkyDirList
import numpy as N
import pyfits
import scipy
import scipy.ndimage


from . roi_diffuse import ROIDiffuseModel_OTF
from . roi_extended import ROIExtendedModel,ROIExtendedModelAnalytic
from . pointspec_helpers import get_default_diffuse_mapper
from uw.utilities import keyword_options
from uw.utilities.fitstools import get_fields
from uw.utilities.image import ZEA
from pypsf import PsfOverlap
import collections
from abc import abstractmethod
import numbers

def memoize(function):
    """ This decorator allows painless caching of a function.
        
        From http://programmingzen.com/2009/05/18/memoization-in-ruby-and-python/,

        It should probably be placed somewhere more generally accessible. """
    cache = {}
    def decorated_function(*args):
        try:
            return cache[args]
        except KeyError:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function

def defaults_to_kwargs(obj,defaults):
    """ Take in a defaults list (used by keyword_options) and an object 
        which recognizes the keyword_options. A dictionary is
        returned with each of the defaults pointing to the value
        found in the object. This is useful for recreating 
        the object. """
    return dict([[i[0],getattr(obj,i[0])] for i in defaults if isinstance(i,list) or isinstance(i,tuple)])


class ROIImage(object):
    """ This object is suitable for creating a SkyImage object
        and filling it with some physically meaningful
        quantity gotten from an ROIAnalysis object. 
        The acutal work is done by subclasses. """

    defaults = ZEA.defaults + (
        ('center',    None, 'Center of image'),
        ('conv_type',   -1, 'Conversion type'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):
        """ Note, unlike ZEA, can support non-square images. To specify a nonsquare
            image, set the size parameter to a lenght two tuple:
                
                size=(10,5) # dx=10 degrees, dy=5 degrees. 
        """
        keyword_options.process(self, kwargs)
        
        self.roi = roi

        # by default, use get energy range and image center from roi.
        if self.center is None: self.center=self.roi.roi_dir

        # set up, then create a SkyImage object to perform the projection
        # to a grid and manage an image
        if not isinstance(self.size,collections.Iterable):
            self.skyimage = SkyImage(self.center, self.fitsfile, self.pixelsize, 
                                     self.size, 1, self.proj, self.galactic, False)
        else:
            self.skyimage = SkyImage(self.center, self.fitsfile, self.pixelsize, 
                                     float(self.size[0]), 1, self.proj, self.galactic, False, float(self.size[1]))

        self.fill()

        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()
        self.image=ROIImage.skyimage2numpy(self.skyimage)

    @staticmethod
    def skyimage2numpy(skyimage):
        nx, ny = skyimage.naxis1(), skyimage.naxis2()
        image=N.array(skyimage.image()).reshape((ny, nx))
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

        zea_dict = dict((d[0],self.__dict__[d[0]]) for d in ZEA.defaults if hasattr(d,'__iter__'))
        if axes is not None: zea_dict['axes']=axes
        if nticks is not None: zea_dict['nticks']=nticks

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
            crval1,crval2=self.center.l(),self.center.b()
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

        cards = [ pyfits.Card(*i) for i in values]

        header=pyfits.Header(cards=cards)

        hdu=pyfits.PrimaryHDU(data=self.image, header=header)
        fits = pyfits.HDUList([hdu])

        return fits

class CountsImage(ROIImage):
    """ This ROIImage subclass fills the sky image with the observed Fermi counts. """

    defaults = ROIImage.defaults

    @staticmethod
    @memoize
    def process_filedata(roi,conv_type=None,extra_cuts=None,radius=None):
        """ The radius parameter will apply a radius cut. """

        if radius is None: radius=roi.sa.maxROI

        emin = roi.bin_edges[0]
        emax = roi.bin_edges[-1]
        conv_type = roi.sa.conv_type if conv_type==None else conv_type

        ft1files=roi.sa.pixeldata.ft1files

        base_cuts = ['ENERGY > %s'% emin,
                     'ENERGY < %s'% emax,
                     'ZENITH_ANGLE < %s' % roi.sa.pixeldata.zenithcut,
                     'THETA < %s' % roi.sa.pixeldata.thetacut,
                     'EVENT_CLASS >= %s' % roi.sa.pixeldata.event_class]
        if conv_type >= 0:        base_cuts += ['CONVERSION_TYPE == %d'%(conv_type)]
        cuts = base_cuts if extra_cuts is None else extra_cuts + base_cuts

        data = get_fields(ft1files,['RA','DEC','time'],cuts)
        # convert into skydirs
        skydirs = [ SkyDir(float(data['RA'][i]),float(data['DEC'][i])) for i in xrange(len(data['RA']))]

        # apply the same gti cut used to read in the initial WSDL.
        gti=roi.sa.pixeldata.gti
        skydirs = [ skydir for skydir,time in zip(skydirs,data['time']) if gti.accept(time)]

        # apply ROI radial cut
        skydirs = [ skydir for skydir in skydirs if N.degrees(skydir.difference(roi.roi_dir)) < radius ]

        return skydirs

    def fill(self):
        dirs = CountsImage.process_filedata(self.roi,self.conv_type)

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
        self.selected_bands = self.roi.bands if self.conv_type < 0 else \
            [ band for band in self.roi.bands if band.ct == self.conv_type ]

        self.wsdl = self.skyimage.get_wsdl()

        self.solid_angle = N.radians(self.pixelsize)**2

        model_counts = N.zeros(len(self.wsdl),dtype=float)

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
            dsarr = N.concatenate([[myarr[i::factor]]
                                   for i in range(factor)]).mean(axis=0)
            return dsarr

        elif len(myarr.shape) == 2:
            xs,ys = myarr.shape
            assert xs % factor == 0 and ys % factor == 0
            dsarr = N.concatenate([[myarr[i::factor,j::factor] 
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

    def all_point_source_counts(self):
        """ Calculate the point source contributions. """
        roi=self.roi
        if self.override_point_sources is None and self.override_diffuse_sources is None:
            point_sources = roi.psm.point_sources 
        else:
            point_sources = self.override_point_sources if self.override_point_sources is not None else []

        point_counts = N.zeros(len(self.wsdl),dtype=float)

        for band in self.selected_bands:
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

        background_counts = N.zeros(len(self.wsdl),dtype=float)

        for band in self.selected_bands:

            ns,bg_points,bg_vector = ROIDiffuseModel_OTF.sub_energy_binning(band,bg.nsimps)

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

        if self.override_point_sources is None and self.override_diffuse_sources is None:
            bgmodels = self.roi.dsm.bgmodels
        else:
            mapper=get_default_diffuse_mapper(self.roi.sa,self.roi.roi_dir)
            if self.override_diffuse_sources is None:
                bgmodels = []
            else:
                bgmodels = [mapper(ds) for ds in self.override_diffuse_sources] 

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

        self.model=ModelImage(self.roi,**defaults_to_kwargs(self,ModelImage.defaults))
        self.counts=CountsImage(self.roi,**defaults_to_kwargs(self,CountsImage.defaults))

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
            ('conv_type',      -1,            'Conversion type'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):
        keyword_options.process(self, kwargs)

        self.npix = float(self.size)/self.pixelsize
        
        self.roi = roi

        # by default, use get energy range and image center from roi.
        if self.center is None: self.center=self.roi.roi_dir

        # bins in theta^2
        self.bin_edges_deg = N.linspace(0.0,self.size**2,self.npix+1)
        self.bin_centers_deg = (self.bin_edges_deg[1:] + self.bin_edges_deg[:-1])/2.0
        
        # two factors of radians b/c theta^2
        self.bin_edges_rad = N.radians(N.radians(self.bin_edges_deg))
        self.bin_centers_rad = N.radians(N.radians(self.bin_centers_deg))

        # the lower and upper agle for each bin.
        self.theta_pairs_rad = zip(N.sqrt(self.bin_edges_rad[:-1]),
                                   N.sqrt(self.bin_edges_rad[1:]))


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
        dirs = CountsImage.process_filedata(self.roi,self.conv_type)
        diffs = [self.center.difference(i) for i in dirs]
        self.image=N.histogram(diffs,bins=N.sqrt(self.bin_edges_rad))[0]

class RadialModel(RadialImage):

    defaults = RadialImage.defaults

    def fill(self):

        self.selected_bands = self.roi.bands if self.conv_type < 0 else \
            [ band for band in self.roi.bands if band.ct == self.conv_type ]

        self.image = N.zeros_like(self.bin_centers_rad)
        self.image += self.all_point_source_counts()
        self.image += self.all_diffuse_sources_counts()
        self.image *= self.roi.phase_factor # don't forget about the phase factor!

    def all_point_source_counts(self):
        """ Calculate the point source contributions. """
        roi=self.roi
        point_sources = roi.psm.point_sources

        point_counts = N.zeros_like(self.bin_centers_rad)

        overlap = PsfOverlap()

        for i,(theta_min,theta_max) in enumerate(self.theta_pairs_rad):

            model_counts=0

            for j,ps in enumerate(point_sources):
                for band in self.selected_bands:

                    # this code requires a redundant call to overlap. Improve if time.
                    fraction=overlap(band,self.center,ps.skydir,radius_in_rad=theta_max) - \
                             overlap(band,self.center,ps.skydir,radius_in_rad=theta_min)

                    model_counts += band.ps_counts[j]*fraction

            point_counts[i] = model_counts

        return point_counts

    def extended_source_counts(self,extended_model):
        if type(extended_model) not in [ROIExtendedModel,ROIExtendedModelAnalytic]:
            raise Exception("Unknown extended model.")

        roi=self.roi

        extended_counts = N.zeros_like(self.bin_centers_rad)

        # get the corresponding mybands.
        mybands=[extended_model.bands[N.where(roi.bands==band)[0][0]] for band in self.selected_bands]

        for band,myband in zip(self.selected_bands,mybands):
            extended_model.set_state(band)

            if type(extended_model) == ROIExtendedModel:

                nside = RadialModel.get_nside(self.size,self.npix)

                temp_band = Band(nside)
                wsdl = WeightedSkyDirList(temp_band,self.center,N.radians(self.size),True)
                vals=extended_model._pix_value(wsdl)

                rvals=N.empty(len(wsdl),dtype=float)
                PythonUtilities.arclength(rvals,wsdl,self.center)

                # get average value in each ring by averaging values.
                fraction = N.histogram(rvals,weights=vals,bins=N.sqrt(self.bin_edges_rad))[0]/\
                           N.histogram(rvals,bins=N.sqrt(self.bin_edges_rad))[0]

                # multiply intensities by solid angle in ring
                fraction *= RadialModel.solid_angle_cone(N.radians(self.size))/self.npix

            elif type(extended_model) == ROIExtendedModelAnalytic:

                fraction = N.empty_like(self.bin_centers_rad)

                for i,(theta_min,theta_max) in enumerate(self.theta_pairs_rad):

                    fraction[i]=extended_model._overlaps(self.center,band,theta_max) - \
                             extended_model._overlaps(self.center,band,theta_min)

            # total counts * fraction = model predictions in each ring.
            extended_counts += myband.es_counts*fraction

        return extended_counts

    @staticmethod
    def solid_angle_cone(radius_in_radians):
        return 2*N.pi*(1-N.cos(radius_in_radians))

    @staticmethod
    def get_nside(size,npix,num_points_per_ring=1000):
        """ Solid angle of each healpix pixel is 4pi/(12*ns^2)
            Solid angel of each ring is pi*(size)^2/npix
            Want size of each ring > num_points_per_ring*size of each healpix (so that
            we get 20 pixels to sample each ring).

            Solving for ns, the size of the healpixels required, we get the required formula
        """
        total_solid_angle=4*N.pi
        image_solid_angle=RadialModel.solid_angle_cone(N.radians(size))

        solid_angle_per_ring=image_solid_angle/npix

        nside = int(N.ceil(N.sqrt(num_points_per_ring*total_solid_angle/(12*solid_angle_per_ring))))
        return nside

    def otf_source_counts(self,bg):

        roi=self.roi

        mo=bg.smodel

        background_counts = N.zeros_like(self.bin_centers_rad)

        for band in self.selected_bands:

            ns,bg_points,bg_vector = ROIDiffuseModel_OTF.sub_energy_binning(band,bg.nsimps)

            nside = RadialModel.get_nside(self.size,self.npix)

            temp_band = Band(nside)
            wsdl = WeightedSkyDirList(temp_band,self.center,N.radians(self.size),True)

            ap_evals = N.empty([len(self.bin_centers_rad),len(bg_points)])

            for ne,e in enumerate(bg_points):

                bg.set_state(e,band.ct,band)

                rvals=N.empty(len(wsdl),dtype=float)
                PythonUtilities.arclength(rvals,wsdl,self.center)
                vals=bg._pix_value(wsdl)

                # get average value in each ring by averaging values.
                ap_evals[:,ne] = N.histogram(rvals,weights=vals,bins=N.sqrt(self.bin_edges_rad))[0]/\
                                 N.histogram(rvals,bins=N.sqrt(self.bin_edges_rad))[0]

            # multiply intensities by solid angle in ring
            ap_evals *= RadialModel.solid_angle_cone(N.radians(self.size))/self.npix

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
        return sum(self.diffuse_source_counts(bg)
                   for bg in self.roi.dsm.bgmodels)


class SmoothedImage(ROIImage):
    
    smoothed_options = (
        ('kerneltype', 'tophat', 'Type of kernel to use'),
        ('kernel_rad',   0.25, 'Sum counts within radius degrees.'),
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

        width = int(N.ceil(kernel_rad/pixelsize))

        if kerneltype == 'tophat':
            kernelsize=4*width

            kernel = N.zeros([kernelsize,kernelsize],dtype=float)

            xx,yy = N.indices(kernel.shape)
            rr = N.sqrt((xx-kernel.shape[0]/2.)**2+(yy-kernel.shape[1]/2.)**2)
            kernel[rr<width] = 1

        elif kerneltype == 'gaussian':
            kernelsize = 8*width

            kernel = N.zeros([kernelsize,kernelsize],dtype=float)
            xx,yy = N.indices(kernel.shape)
            rr = N.sqrt((xx-kernel.shape[0]/2.)**2+(yy-kernel.shape[1]/2.)**2)
            kernel = N.exp(-(rr**2)/(2*width**2)) / (width**2 * (2*N.pi))


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
        self.pass_dict=defaults_to_kwargs(self,self.object.defaults)
        self.pass_dict['size']=self.size+self.kernelsize*self.pixelsize
        self.smoothed=self.object(self.roi,**self.pass_dict)


        self.smoothed.skyimage=SmoothedImage.convolve_skyimage(self.smoothed.skyimage,self.kernel)
        self.smoothed.image=ROIImage.skyimage2numpy(self.smoothed.skyimage)

        # now, shrink down smoothed image and replace the current skyimage with it

        self.image = self.smoothed.image[self.kernelsize/2:-self.kernelsize/2,self.kernelsize/2:-self.kernelsize/2]
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
