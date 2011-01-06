"""
This code provides obejcts to generate counts and model counts images with
two dimensional projections. The code interfaces with the ROIAnalysis
object for making the calculations, skymaps.SkyImage object for storing
the data, and the image.ZEA object for plotting.  The high level object
roi_plotting.ROIDisplay can use to access these objects form a high
level plotting interface.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_image.py,v 1.3 2010/12/12 07:05:39 lande Exp $

author: Joshua Lande
"""
from skymaps import SkyImage,SkyDir,PythonUtilities,Band,WeightedSkyDirList
import numpy as N

from roi_diffuse import ROIDiffuseModel_OTF
from roi_extended import ROIExtendedModel,ROIExtendedModelAnalytic
from uw.utilities import keyword_options
from uw.utilities.fitstools import get_fields
from uw.utilities.image import ZEA
from pypsf import PsfOverlap
import collections

class ROIImage(object):
    """ This object is suitable for creating a SkyImage object
        and filling it with some physically meaningful
        quantity gotten from an ROIAnalysis object. 
        The acutal work is done by subclasses. """

    defaults = ZEA.defaults + (
            ('center',    None,  'Center of image'),
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
        self.image=N.array(self.skyimage.image()).reshape((self.ny, self.nx))

    def fill():
        raise NotImplementedError("Subclasses should implement this!")

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
        zea.image = self.image

        # The old one gets removed by python's garbage collector (when zea.skyimage is replaced).
        zea.projector = zea.skyimage.projector() 

        zea.vmin,zea.vmax = zea.skyimage.minimum(), zea.skyimage.maximum()
        return zea


class CountsImage(ROIImage):
    """ This ROIImage subclass fills the sky image with the observed Fermi counts. """

    defaults = ROIImage.defaults + (
        ('mc_src_id', None, ''),
        ('cuts',      None, '')
    )

    @staticmethod
    def process_filedata(roi,mc_src_id,extra_cuts):

        emin = roi.bin_edges[0]
        emax = roi.bin_edges[-1]
        conv_type = roi.sa.conv_type

        ft1files=roi.sa.pixeldata.ft1files

        base_cuts = ['ENERGY > %s'% emin,
                     'ENERGY < %s'% emax,
                     'ZENITH_ANGLE < %s' % roi.sa.pixeldata.zenithcut,
                     'THETA < %s' % roi.sa.pixeldata.thetacut,
                     'EVENT_CLASS == %s' % roi.sa.pixeldata.event_class]
        if conv_type >= 0:        base_cuts += ['EVENT_CLASS == %d'%(conv_type)]
        if mc_src_id is not None: base_cuts += ['MC_SRC_ID == %d'%(mc_src_id)]
        cuts = base_cuts if extra_cuts is None else extra_cuts + base_cuts

        data = get_fields(ft1files,['RA','DEC','time'],cuts)
        # convert into skydirs
        skydirs = [ SkyDir(float(data['RA'][i]),float(data['DEC'][i])) for i in xrange(len(data['RA']))]

        # apply the same gti cut used to read in the initial WSDL.
        gti=roi.sa.pixeldata.gti
        skydirs = [ skydir for skydir,time in zip(skydirs,data['time']) if gti.accept(time)]
        return skydirs

    def fill(self):
        dirs = CountsImage.process_filedata(self.roi,self.mc_src_id,self.cuts)

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
        dsarr = N.concatenate([[crarr[i::factor,j::factor] 
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
                self.fine_skyimage = SkyImage(self.center, self.fitsfile, float(self.pixelsize)/self.factor,
                                         self.size, 1, self.proj, self.galactic, False)
            else:
                self.fine_skyimage = SkyImage(self.center, self.fitsfile, float(self.pixelsize)/self.factor,
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
        return sum(self.diffuse_source_counts(bg)
                   for bg in self.roi.dsm.bgmodels)



class RadialImage(object):
    """ This object is similar to ROIImage but performs a radial
        integral around around a given direction in the sky. """


    defaults = (
            ('center',       None, 'Center of image'),
            ('size',            2, 'Size of image (in degrees)'), 
            ('pixelsize', 0.00625, """ size of each image pixel. This is a little misleading because the
                                       size of each pixel varies, but regardless this is used to determine
                                       the total number of pixels with npix=size/pixelsize """),
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

    def fill():
        """ This should fill up self.image appropriatly."""
        raise NotImplementedError("Subclasses should implement this!")


class RadialCounts(RadialImage):
    """ This subclass of RadialImage calculates the counts within
        each radial bin. """

    defaults = RadialImage.defaults + (
        ('mc_src_id', None, ''),
        ('cuts',      None, '')
    )

    def fill(self):
        dirs = CountsImage.process_filedata(self.roi,self.mc_src_id,self.cuts)
        diffs = [self.center.difference(i) for i in dirs]
        self.image=N.histogram(diffs,bins=N.sqrt(self.bin_edges_rad))[0]

class RadialModel(RadialImage):

    defaults = RadialImage.defaults

    def fill(self):

        self.image=N.zeros_like(self.bin_centers_rad)

        self.image += self.all_point_sources_counts()
        self.image += self.all_diffuse_sources_counts()

    def all_point_sources_counts(self):
        """ Calculate the point source contributions. """
        roi=self.roi
        bands=roi.bands
        point_sources = roi.psm.point_sources

        point_counts = N.zeros_like(self.bin_centers_rad)

        overlap = PsfOverlap()

        for i,(theta_min,theta_max) in enumerate(self.theta_pairs_rad):

            model_counts=0

            for j,ps in enumerate(point_sources):
                for band in bands:

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
        bands=roi.bands

        extended_counts = N.zeros_like(self.bin_centers_rad)

        for band,myband in zip(bands,extended_model.bands):
            extended_model.set_state(band)

            if type(extended_model) == ROIExtendedModel:

                num_points_per_ring=20
                nside = int(N.ceil(N.sqrt(num_points_per_ring*(4./12.)*(self.npix/N.radians(self.size)**2))))

                temp_band = Band(nside)
                wsdl = WeightedSkyDirList(temp_band,self.center,N.radians(self.size),True)
                vals=extended_model._pix_value(wsdl)

                rvals=N.empty(len(wsdl),dtype=float)
                PythonUtilities.arclength(rvals,wsdl,self.center)

                # get average value in each ring by averaging values.
                fraction = N.histogram(rvals,weights=vals,bins=N.sqrt(self.bin_edges_rad))[0]/\
                           N.histogram(rvals,bins=N.sqrt(self.bin_edges_rad))[0]

                # multiply intensities by solid angle in ring
                fraction *= 2*N.pi*(1-N.cos(N.radians(self.size)))/self.npix

            elif type(extended_model) == ROIExtendedModelAnalytic:

                fraction = N.empty_like(self.bin_centers_rad)

                for i,(theta_min,theta_max) in enumerate(self.theta_pairs_rad):

                    fraction[i]=extended_model._overlaps(self.center,band,theta_max) - \
                             extended_model._overlaps(self.center,band,theta_min)

            # total counts * fraction = model predictions in each ring.
            extended_counts += myband.es_counts*fraction

        return extended_counts


    def otf_source_counts(self,bg):
        """
            Solid angle of each healpix pixel is 4pi/(12*ns^2)
            Solid angel of each ring is 4*(size)^2/npix
            Want size of each ring > 20*size of each healpix (so that
            we get 20 pixels to sample each ring).

            Solving for ns, the size of the healpixels required, we get

            ns=ceil(20*(4/12)*npix/size^2) 
            
            with where size is in radians. 
        """

        roi=self.roi
        bands=roi.bands

        mo=bg.smodel

        background_counts = N.zeros_like(self.bin_centers_rad)

        for band in bands:

            ns,bg_points,bg_vector = ROIDiffuseModel_OTF.sub_energy_binning(band,bg.nsimps)

            num_points_per_ring=20
            nside = int(N.ceil(N.sqrt(num_points_per_ring*(4./12.)*(self.npix/N.radians(self.size)**2))))

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
                ap_evals[:,ne] *= 2*N.pi*(1-N.cos(N.radians(self.size)))/self.npix

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
