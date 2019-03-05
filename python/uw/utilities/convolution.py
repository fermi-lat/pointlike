"""Module to support on-the-fly convolution of a mapcube for use in spectral fitting.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/convolution.py,v 1.47 2013/01/25 00:53:35 lande Exp $

authors: M. Kerr, J. Lande

"""
from skymaps import SkyDir,BaseWeightedSkyDirList,Hep3Vector,SkyIntegrator,PySkyFunction,Background,PythonUtilities
from pointlike import DoubleVector
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad,romberg,cumtrapz
from scipy.special import hyp2f1
from uw.like.pypsf import BandCALDBPsf,PretendBand
from uw.like.SpatialModels import SpatialMap
from numpy.fft import fftshift,ifft2,fft2
import keyword_options

### TODO -- find a way to evaluate the PSF on a finer grid ###

class Grid(object):
    """ Define a grid with vertices approximately separated by a uniform arclength.

        Although the routine is agnostic, Galactic coordinates are used
        throughout.  Any lon/lat argument is understood to be GLON/GLAT.

        The default grid is established on the equator, at the azimuth
        specified by GLON of the center, and is then rotated to the 
        appropriate latitude.  Since a CAR projection is used, this rotation
        minimizes the distortion and non-distance-preserving properties of the
        projection.

        The convention for the coordinates is that the origin is in the upper
        lefthand corner.  The function is evaluated at pixel *corners*.  This puts the
        center of the image in the center of the (unique) center pixel.  The convention
        makes interpolation easier, i.e., there are no edge effects.

        Note -- if npix is a power of 2, the FFT algorithm will operate most quickly.
        However, with an even number of pixels, the function is not evaluated at the
        exact center of the grid, which can be deleterious for evaluating the PSF at
        high energies.  That is to say...

        *******************************************************************
        ***ALWAYS USE AN ODD NUMBER OF PIXELS IF CONVOLVING WITH THE PSF***
        *******************************************************************

        Note -- longitude increases "right to left".

        Note -- to display, imshow(Grid.vals.transpose()[::-1]) works -- a little cumbersome.
    """

    def __init__(self,center,npix=4*25+1,pixelsize=1./4,bounds_error=True,fill_value=0):
        """ center -- a SkyDir giving the center of the grid
            npix   -- the number of pixels along a grid side
            pixelsize -- the length of the side of a pixel (deg)
        """
        self.set_center(center)
        self.setup_grid(npix,pixelsize)

        self.bounds_error=bounds_error
        self.fill_value=fill_value

    def set_center(self,center):
        self.center = center
        self.clon,self.clat = center.l(),0
        self.rot_axis = SkyDir(self.clon+90,0).dir() # Hep3Vector

    def setup_grid(self,npix=4*25+1,pixelsize=1./4):

        self.npix,self.pixelsize = npix,pixelsize
        self.delta_lon = (npix-1)*pixelsize
        self.delta_lat = (npix-1)*pixelsize
        self.lon0 = (self.clon + self.delta_lon/2.)%360
        self.lat0 = self.clat - self.delta_lat/2.
        self.wrap = self.lon0 <= self.delta_lon # True if origin (0/360) is in frame

        mpix  = (float(npix)-1)/2.
        v     = (np.arange(0,npix) - mpix)**2
        self.dists = np.asarray([v+v[i] for i in xrange(npix)])**0.5*(np.radians(self.pixelsize))

        self.lons = self.lon0 - np.linspace(0,1,self.npix)*self.delta_lon
        self.lons = DoubleVector(np.where(self.lons < 0,self.lons+360,self.lons)) # if origin in frame
        self.lats = DoubleVector(np.linspace(0,1,self.npix)*self.delta_lat + self.lat0)

    def rot_lon_lat(self,lon,lat):
        """ Rotate the specified lon/lat to the equatorial grid. """
        v1 = DoubleVector(); v2 = DoubleVector()
        Background.rot_grid(v1,v2,DoubleVector([lon]),DoubleVector([lat]),self.center)
        return v1[0],v2[0]

    def in_grid(self,lon,lat,skydir=None):
        """ check whether the provided lon/lat is within the grid. """
        if skydir is not None:
            lon = skydir.l(); lat = skydir.b()
        lon,lat = self.rot_lon_lat(lon,lat)
        delta_lon = abs(lon - self.clon)
        delta_lon = min(360-delta_lon,delta_lon)
        lon_ok = delta_lon < self.delta_lon/2.
        lat_ok = abs(lat - self.clat) < self.delta_lat/2.

        return lat_ok and lon_ok

    def pix(self,skydir):
        """ return the pixel corresponding to the argument skydir. """
        lon,lat = self.rot_lon_lat(skydir.l(),skydir.b())
        if self.wrap:
            lon = np.where(lon > 180, lon - 360 , lon)
        npix = self.npix - 1
        x = (self.lon0 - lon)/self.delta_lon*npix
        y = (lat - self.lat0)/self.delta_lat*npix
        return x,y
    
    def __call__(self,skydir,v):
        """ Using v, an array of values that has been evaluated on the Grid (e.g., by fill),
            find the value(s) corresponding to skydir (can be single or a list) using
            bilinear interpolation.
            
            The skydir(s) are rotated onto the equatorial grid."""
        
        if hasattr(skydir,'EQUATORIAL'): skydir = [skydir]
        rvals = np.empty(len(skydir)*2,dtype=float)
        PythonUtilities.rot_grid(rvals,skydir,self.center)
        lon = rvals[::2]; lat = rvals[1::2]
        if self.wrap:
            # adopt negative longitudes for the nonce; fine for calculating differences
            lon = np.where(lon > 180, lon - 360 , lon)
        npix = self.npix - 1
        x = (self.lon0 - lon) / (self.delta_lon / npix)
        y = (lat - self.lat0) / (self.delta_lat / npix)
        #if np.any( (x<0) || (x > npix) || (y < 0) || (y > npix) ):
        #    return np.nan
        xlo,ylo = np.floor(x+1e-6).astype(int),np.floor(y+1e-6).astype(int)
        xhi,yhi = xlo+1,ylo+1
        dx = np.maximum(0,x - xlo)
        dy = np.maximum(0,y - ylo)
        v  = np.asarray(v)
        if self.bounds_error: 
            return v[xlo,ylo]*(1-dx)*(1-dy) + v[xlo,yhi]*(1-dx)*dy + v[xhi,ylo]*dx*(1-dy) + v[xhi,yhi]*dx*dy
        else:
            return np.where((xlo<0) | (ylo<0) | (xhi>npix) | (yhi>npix),self.fill_value,
                           v[np.clip(xlo,0,npix),np.clip(ylo,0,npix)]*(1-dx)*(1-dy) + \
                           v[np.clip(xlo,0,npix),np.clip(yhi,0,npix)]*(1-dx)*dy + \
                           v[np.clip(xhi,0,npix),np.clip(ylo,0,npix)]*dx*(1-dy) + \
                           v[np.clip(xhi,0,npix),np.clip(yhi,0,npix)]*dx*dy)

    def fill(self,skyfun):
        """ Evaluate skyfun along the internal grid and return the resulting array.
        
            In this process, the internal grid is transformed to the center of
            the ROI, and the skyfun evaluated for each of the resulting
            positions.
            
            The rotation is rather fast, likely no need to pre-compute.
        """
        v = np.empty(self.npix*self.npix)
        PythonUtilities.val_grid(v,self.lons,self.lats,self.center,skyfun)
        return v.reshape([self.npix,self.npix])
        

class BackgroundConvolution(Grid):

    def __init__(self,center,bg,psf,*args,**kwargs):
        """ center -- a SkyDir giving the center of the grid on which to convolve bg
            bg     -- an instance of skymaps::Background
            psf    -- and instance of CALDBPsf

            Additional kwargs are passed to Grid.
        """

        self.bg,self.psf = bg,psf
        super(BackgroundConvolution,self).__init__(center,**kwargs)
        if (self.npix%2 != 1):
            print """WARNING! You are attempting to convolve a map on a
                     grid with an even number of pixels.  This is NOT
                     recommended as it will smear the map inappropriately
                     at high energies."""

    def do_convolution(self,energy,conversion_type,override_skyfun=None,override_vals=None,
                       override_en=None):
        """ Perform a convolution at the specified energy, for the specified
            conversion type.  The values are stored internally as "cvals".
        """
        if override_skyfun is None and override_vals is None:
            ## here is where the conversion type is used to define which exposure to use
            ## isotropic might have only one 
            ct = conversion_type if conversion_type<self.bg.exposures() else 0
            self.bg.set_skyfun(ct, energy) 
            self.bg_vals = self.fill(self.bg)
        elif override_vals is None:
            self.bg_vals = self.fill(override_skyfun)
        else:
            self.bg_vals = override_vals
        self.bg_vals[np.isnan(self.bg_vals)] = 0

        pb = PretendBand(energy,conversion_type)
        bpsf = BandCALDBPsf(self.psf,pb,override_en=override_en,adjust_mean=False)
        self.psf_fill(bpsf)
        self.convolve()


    def psf_fill(self,psf):
        """ Evaluate a band psf over the grid."""
        psf_vals = psf(self.dists,density=True).reshape([self.npix,self.npix])
        self.psf_vals = psf_vals / psf_vals.sum()
        #self.psf_vals = psf_vals*np.radians(self.pixelsize)**2

    def convolve(self):
        """ Perform the convolution with the current values of the bg
            and psf evaluated over the grid."""
        fft_kernel = fft2(self.psf_vals)
        fft_image  = fft2(self.bg_vals)
        self.cvals = c = np.real(fftshift(ifft2(fft_kernel*fft_image)))

        # swap the 0th component into the proper place
        new = np.empty_like(c)
        new[0,:]  = c[-1,:]
        new[1:,:] = c[:-1,:]
        c[:,0]  = new[:,-1]
        c[:,1:] = new[:,:-1]

    def ap_average(self,radius,convolved=True):
        """ Estimate the average of the background over a radial aperture by
            averaging all (convolved) bg pixels in the interior."""
        v = self.cvals if convolved else self.bg_vals
        return (v[self.dists <= radius]).mean()

    def show(self):
        """ A quick triple plot for a sanity check on the convolution."""
        import pylab as P
        marker = float(self.npix)/2
        P.subplot(131)
        P.imshow(self.psf_vals,interpolation='nearest')
        norm = P.normalize(np.log10(self.bg_vals.min()),np.log10(self.bg_vals.max()))
        P.subplot(132)
        P.imshow(np.log10(self.bg_vals).transpose()[::-1],norm=norm,interpolation='nearest')
        P.axvline(marker,color='k')
        P.axhline(marker,color='k')
        P.subplot(133)
        P.imshow(np.log10(self.cvals).transpose()[::-1],norm=norm,interpolation='nearest')
        P.axvline(marker,color='k')
        P.axhline(marker,color='k')

#===============================================================================================#

class ExtendedSourceConvolution(BackgroundConvolution):
    """ This object is suitable for spatial models which are supposed to be 
        normalized to 1. It is also intended for spatial models where the entire 
        convolution radius is inside of the total roi.
        
        Note that this implementation assumes that there is no exposure variation
        across the source and ap_average must be mulitplied by the exposure
        outside of this function."""

    defaults =[
        ['pixelsize', 0.025, 'convolution resolution'],
        ['npix'   ,   101,  'Initial value gets reset automatically by do_convolution.'],
        ['r_multi',   2.0,  'multiple of r95 to set max dimension of grid'],
        ['r_max',     20,   'an absolute maximum (half)-size of grid (deg)'],
        ]
    @staticmethod
    def set_pixelsize(psize):
        was = ExtendedSourceConvolution.defaults[0][1]
        ExtendedSourceConvolution.defaults[0][1] =psize
        return was
    @keyword_options.decorate(defaults)
    def __init__(self,extended_source,psf, **kwargs):
        """ extended_source : object containing a spatial_model object
            psf : PSF object used for convolution
        """

        keyword_options.process(self, kwargs)
        self.extended_source = extended_source

        # Pass in none for the 
        super(ExtendedSourceConvolution,self).__init__(self.extended_source.spatial_model.center,None,psf,
                npix=self.npix,pixelsize=self.pixelsize,
                bounds_error=False,fill_value=0)

    def convolve(self,*args,**kwargs):
        super(ExtendedSourceConvolution,self).convolve(*args,**kwargs)
        self.cvals /= self.cvals.sum()*np.radians(self.pixelsize)**2

    def overlap(self,roi_center,roi_radius):
        """ Calculate the fraction of PDF contained within ROI (which
            is defined by the input radius and center). This is the extended
            source analog of the pypsf.PsfOverlap object. It is convenient
            to return just the overlap fraction since it is known that
            the entire spatial part is normalized.

            Note that this formula assumes the entire extended source is
            within the grid, which should be enforced by a reasonable
            choice of r_max and r_multi.
            
            roi_radius is in radians."""
        
        # If grid is entirly inside of roi, overlap is 1
        if self.center.difference(roi_center) + \
                np.radians(self.pixelsize)*self.npix < roi_radius:
            return 1.0

        x,y = self.pix(roi_center)
        dx =((np.arange(0,self.npix)-x)**2).reshape(self.npix,1)
        dy =((np.arange(0,self.npix)-y)**2)
        d  = np.sqrt(dx+dy)*np.radians(self.pixelsize)
        return (self.cvals[d <= roi_radius]).sum()*np.radians(self.pixelsize)**2

    def do_convolution(self,band):

        # recenter convolution grid to (possibly) new extended source center.
        self.set_center(self.extended_source.spatial_model.center)

        # Use the 'optimal' energy (calculated by the ADJUST_MEAN flag) if it exists.
        energy = band.psf.eopt if band.psf.__dict__.has_key('eopt') else band.e

        # Extended sources might be much bigger then the ROI if the convolution is only
        # being redone in a small area (for example, to make a radial profile.
        edge_distance=np.degrees(band.sd.difference(self.extended_source.spatial_model.center) + band.radius_in_rad)
        extended_src_edge=self.extended_source.spatial_model.effective_edge()
        edge=min(edge_distance,extended_src_edge)

        r95 = self.psf.inverse_integral(energy,band.ct,95)
        rad = self.r_multi*r95 + edge
        rad = max(min(self.r_max,rad),edge+2.5)
        self.npix = int(round(2*rad/self.pixelsize))
        self.npix += (self.npix%2 == 0)

        self.setup_grid(self.npix,self.pixelsize)

        if hasattr(self.extended_source.spatial_model,'fill_grid'):
            # Generally, extended soruces should implment the fill_grid function.
            vals=self.extended_source.spatial_model.fill_grid(self,energy)
            super(ExtendedSourceConvolution,self).do_convolution(energy,band.ct,override_vals=vals)

        elif isinstance(self.extended_source.spatial_model,SpatialMap):
            # For the SpatialMap object, there is already a C++ implmentation of the SkyFunction,
            # so we can quickly fill from it.
            super(ExtendedSourceConvolution,self).do_convolution(energy,band.ct,
                                   override_skyfun=self.extended_source.spatial_model.skyfun)
        else:
            # This part of the code wraps the spatial_model object as a PySkySpectrum object 
            # and is very poorly optimized. Hopefully it will only exist as a fallback 
            # because extended sources will either
            # (a) be radially symmetric
            # (b) Immediatly fill the convolution grid in python
            # (c) Have a C++ implementation of the sky function interface so that 
            #     it can be filled in python.
            super(ExtendedSourceConvolution,self).do_convolution(energy,band.ct)


#===============================================================================================#

class AnalyticConvolution(object):
    """ Calculates the convolution of the psf with a radially symmetric spatial_model. 
        This object has a similar interface to BackgroundConvolution.         

        Note that this implementation assumes that there is no exposure variation
        across the source.        
        """

    defaults = (
        ['num_points',     220, 'Number of points to calculate the PDF at. Interpolation is done in between.'],
        ['roi_pad', 1, 'Number of degress away from ROI to perform convolution (helps with caching'],
    )
    @staticmethod
    def set_points(psize):
        was = AnalyticConvolution.defaults[0][1]
        AnalyticConvolution.defaults[0][1] =psize
        return was

    @keyword_options.decorate(defaults)
    def __init__(self,extended_source,psf,**kwargs):
        """ spatial_model  a SpatialModel object
            psf            a PSF object.
        """

        keyword_options.process(self, kwargs)

        self.extended_source = extended_source
        self.psf=psf


    def _convolve(self,rlist,g,s,energy):
        """ Function to calculate the pdf at the givin radii points
            for a given psf sigma and gamma. (and a radially symmetric
            spatial model self.extended_source.spatial_model

            For each givin radii, we must integrate to to the 
            maximum source radius

        integrate to the edge of the source, which is decided on a
        source by source basis (so that sharp edge sources can define
        automatically. The default is a some constant times r68 since
        the intgral includes at term which is the PDF, the integral
        will presumably contribute very littel further away then this.
        """
        int_max = 0.5*(np.radians(self.extended_source.spatial_model.effective_edge(energy))/s)**2

        # u value corresponding to the given r.
        ulist=0.5*(self.rlist/s)**2

        # pdf values for each r (in a given energy bin)
        pdf=np.empty_like(self.rlist)

        integrand = lambda v,u: self.extended_source.spatial_model.at_r(np.sqrt(2*v)*s)\
                                          *((g-1)/g)*(g/(g+u+v))**g\
                                          *hyp2f1(g/2.,(1+g)/2.,1.,4.*u*v/(g+u+v)**2)
        for i,u in enumerate(ulist):
            pdf[i]=quad(integrand,0,int_max,args=(u,),epsrel=1e-3,
                        epsabs=1e-3,full_output=True)[0]
            # Sometimes this algorithm is not robust. In that case,
            # try to do the integral more accuratly
            if np.isnan(pdf[i]) or np.isinf(pdf[i]) or pdf[i]<0:

                # try integrating to infinity
                pdf[i]=quad(integrand,0,np.inf,args=(u,),epsrel=1e-3,
                            epsabs=1e-3,full_output=True)[0]

                if np.isnan(pdf[i]) or np.isinf(pdf[i]) or pdf[i]<0:

                    # try to use romberg to do the integral.
                    pdf[i]=romberg(integrand,0,int_max,args=(u,),divmax=20)
                                   
        return pdf
    
    def _get_pdf(self,energy,conversion_type,band):
        """ This function has to calculate self.rlist and self.pdf
            and is abstracted from the rest of the do_convolution
            function so that it can be overloaded for caching
            by AnalyticConvolutionCache. """

        pb = PretendBand(energy,conversion_type)
        self.bpsf = BandCALDBPsf(self.psf,pb)

        self.edge_distance=band.sd.difference(self.extended_source.spatial_model.center) + \
                band.radius_in_rad + self.roi_pad

        self.rmax=self.edge_distance

        self.rlist=np.linspace(0,np.sqrt(self.rmax),self.num_points)**2
        self.pdf=np.zeros_like(self.rlist)

        if self.bpsf.newstyle:

            nclist,ntlist,gclist,gtlist,\
                    sclist,stlist,wlist = self.bpsf.par

            # apparently, in the new BandCALDBPsf implementation,
            # nc and nt are no longer fractions from 0 to 1 but also
            # have the density term in it. Too keep compatable
            # with the extended soruce code, undo the density
            # part of nc and nt.
            nclist*=2*np.pi*sclist**2 
            ntlist*=2*np.pi*stlist**2 

            smax = np.append(sclist,stlist).max()

            # For P7SOURCE_V6, there is no theta dependence to the
            # psf so we can just do the convolution once!
            if np.allclose(nclist,1) and np.allclose(ntlist,0) and \
               np.allclose(gclist,gclist[0]) & np.allclose(sclist,sclist[0]):

                self.pdf += self._convolve(self.rlist,gclist[0],sclist[0],energy)

            else:
                for nc,gc,sc,nt,gt,st,w in zip(nclist,gclist,sclist,\
                                               ntlist,gtlist,stlist,wlist):
                    self.pdf += w*(nc*self._convolve(self.rlist,gc,sc,energy)+
                                   nt*self._convolve(self.rlist,gt,st,energy))
        else:
            raise Exception("I am not at all sure that the extended source code works with the old-style PSF any more.")

            # g = gamma, s = sigma, w = weight
            glist,slist,wlist = self.bpsf.par
            smax=slist.max()

            for g,s,w in zip(glist,slist,wlist):
                self.pdf += w*self._convolve(self.rlist,g,s,energy)

        # for some reason, I incorrectly got a negative pdf value for r=0 and for
        # especially large gaussian MC sources. Not sure how the integral of the 
        # hypergeometric function could do this, but it is best to simply remove 
        # it from the list since we know they are unphysical.
        bad = np.isnan(self.pdf)|np.isinf(self.pdf)|(self.pdf<0)
        if np.any(bad):
            message='WARNING! Bad values found in PDF. Removing them from interpolation.' % sum(bad)
            if np.any(np.isnan(self.pdf)): message += ' (%d NaN values)' % sum(np.isnan(self.pdf))
            if np.any(np.isinf(self.pdf)): message += ' (%d Inf values)' % sum(np.isinf(self.pdf))
            if np.any(self.pdf<0):        message +=' (%d negative values)' % sum(self.pdf<0)
            raise Exception(message)

    def do_convolution(self,band):
        """ Generate points uniformly in r^2 ~ u, to ensure there are more 
            points near the center, where the PDF is bigger and changing
            rapidly. Then, do a cubic spline interpolation of the log of
            the PDF values. Why like this? 
            
            The PSF is (1-1/g)(1+u/g)^-g, so

            log(PSF) = log(1-1/g) -g*log(1+u/g)

            So a plot of log(PSF) vs r^2 is close to being linear. In
            the limit of u/g being small (near the peak of the PSF where
            the convolution is most important), we have log(1+u/g) ~ u/g
            so 
            
            log(PSF) ~ A + B*u

            which is why it makes sense to do an interpolation of log(PSF)
            vs u which should be close to a straight line and interpolate 
            with little error. 
            
            Of course, this argument isn't really correct because we are
            interpolationt he PDF not the PSF, but I assume convolving
            the PSF with an extended shape like a gaussian would also
            create something that looks like a gaussian. 
            """
        # Use the 'optimal' energy (calculated by the ADJUST_MEAN flag) if it exists.
        energy = band.psf.eopt if band.psf.__dict__.has_key('eopt') else band.e

        # do most of the work in this subfunction, which allows for caching
        # by the AnalyticConvolutionCache subclass.
        self._get_pdf(energy,band.ct,band)

        # Assume pdf is 0 outside of the bound, which is reasonable if 
        # rmax is big enough also, interpolate the log of the intensity, which 
        # should for the types of shapes we are calculating be much more linear.
        self.log_pdf=np.log(self.pdf)
        self.interp=interp1d(self.rlist,self.log_pdf,kind=3,bounds_error=False,fill_value=self.log_pdf[0])

        self.val=lambda x: np.exp(self.interp(x))*(x<=self.rlist[-1])

        # Calculate the integral of the pdf. Then do an interploation of it.
        # This is used for the integrate function, which is supposed to
        # be analogous to the BandCALDBPsf.integral function.
        self.int_pdf=cumtrapz(x=self.rlist,y=2*np.pi*self.rlist*self.pdf)
        self.int_rlist=self.rlist[1:]
        self.int_interp=interp1d(self.int_rlist,self.int_pdf,
                                        bounds_error=False,fill_value=1)
        # fill 0 for smaller radii. Fill 1 above
        self.int_val=lambda x:self.int_interp(x)*(x>=self.int_rlist[0])

    def __call__(self,skydir):
        """ This funciton is analogous to the BandCALDBPsf.__call__ function
            except that it always returns the density (probability per unit
            area). Also, it is different in that it takes in a skydir or WSDL 
            instead of a radial distance. """
        if isinstance(skydir,BaseWeightedSkyDirList):
            difference = np.empty(len(skydir),dtype=float)
            PythonUtilities.arclength(difference,skydir,self.extended_source.spatial_model.center)
            return self.val(difference)
        elif type(skydir)== np.ndarray:
            return self.val(skydir)
        elif type(skydir)==list and len(skydir)==3:
            skydir = SkyDir(Hep3Vector(skydir[0],skydir[1],skydir[2]))

        elif type(skydir)==SkyDir:
            return float(self.val(skydir.difference(self.extended_source.spatial_model.center)))
        else:
            raise Exception("Unknown input to AnalyticConvolution.__call__()")

    def integral(self,dmax,dmin=0):
        """ This funciton is analogous to the BandCALDBPsf.integral function
            Note that when dmax=0 (with dmin=0), the integral equals 0. When
            dmax->infty (with dmin=0), the integral aproaches 1. """
        return float(self.int_val(dmax)-self.int_val(dmin))

class AnalyticConvolutionCache(AnalyticConvolution):
    """ A child of the AnalyticConvolution object which should
        perform the exact same cacluation. 
        
        This object implements a caching mechanism so that when the
        same spatial parameters (except for the first 2) are passed  
        as the last time, the previous value is
        returned instead of being recalculated. """

    defaults = AnalyticConvolution.defaults

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):

        super(AnalyticConvolutionCache,self).__init__(*args,**kwargs)

        self.last_p         = {}
        self.last_pdf       = {}
        self.last_rlist     = {}
        self.last_center    = {}

    def _get_pdf(self,energy,conversion_type,band):

        sm=self.extended_source.spatial_model

        if self.last_p.has_key(band) and \
           np.all(self.last_p[band] == sm.p[2:]) and \
           np.degrees(sm.center.difference(self.last_center[band])) < self.roi_pad:

            self.rlist   = self.last_rlist[band]
            self.pdf     = self.last_pdf[band]

        else:
            super(AnalyticConvolutionCache,self)._get_pdf(
                    energy,conversion_type,band)

            self.last_p[band]=sm.p[2:].copy()
            self.last_rlist[band]=self.rlist
            self.last_pdf[band]=self.pdf
            self.last_center[band]=sm.center

"""
a little sanity check class
from skymaps import *
from uw.utilities.image import ZEA
class GridPySkyFun(object):

    def __init__(self,grid,sf):
        self.grid = grid
        self.grid_vals = self.grid.fill(sf)
        self.sf = sf

    def __call__(self,v):
        sd = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        return self.grid(sd,self.grid_vals)

    def get_pyskyfun(self):
        return PySkyFunction(self)

    def do_zea(self,size=10):

        z = ZEA(self.grid.center,size=size,
                pixelsize=self.grid.pixelsize/2,galactic=True)
        z.fill(self.get_pyskyfun())
        z.set_axes()
        z.imshow()
        #z.axes.imshow(z.image,origin='lower')
        self.z = z
"""
