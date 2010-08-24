"""Module to support on-the-fly convolution of a mapcube for use in spectral fitting.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/convolution.py,v 1.17 2010/08/17 02:05:51 lande Exp $

authors: M. Kerr, J. Lande

"""
from skymaps import SkyDir,WeightedSkyDirList,Hep3Vector,SkyIntegrator,PySkyFunction,Background,PythonUtilities
from pointlike import DoubleVector
import numpy as N
from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.special import hyp2f1
from uw.like.pypsf import BandCALDBPsf,PretendBand
from numpy.fft import fftshift,ifft2,fft2

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
        self.center = center
        self.clon,self.clat = center.l(),0
        self.rot_axis = SkyDir(self.clon+90,0).dir() # Hep3Vector
        self.setup_grid(npix,pixelsize)

        self.bounds_error=bounds_error
        self.fill_value=fill_value

    def setup_grid(self,npix=4*25+1,pixelsize=1./4):

        self.npix,self.pixelsize = npix,pixelsize
        self.delta_lon = (npix-1)*pixelsize
        self.delta_lat = (npix-1)*pixelsize
        self.lon0 = (self.clon + self.delta_lon/2.)%360
        self.lat0 = self.clat - self.delta_lat/2.
        self.wrap = self.lon0 <= self.delta_lon # True if origin (0/360) is in frame

        mpix  = (float(npix)-1)/2.
        v     = (N.arange(0,npix) - mpix)**2
        self.dists = N.asarray([v+v[i] for i in xrange(npix)])**0.5*(N.radians(self.pixelsize))

        self.lons = self.lon0 - N.linspace(0,1,self.npix)*self.delta_lon
        self.lons = DoubleVector(N.where(self.lons < 0,self.lons+360,self.lons)) # if origin in frame
        self.lats = DoubleVector(N.linspace(0,1,self.npix)*self.delta_lat + self.lat0)

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
            lon = N.where(lon > 180, lon - 360 , lon)
        np = self.npix - 1
        x = (self.lon0 - lon)/self.delta_lon*np
        y = (lat - self.lat0)/self.delta_lat*np
        return x,y
    
    def __call__(self,skydir,v):
        """ Using v, an array of values that has been evaluated on the Grid (e.g., by fill),
            find the value(s) corresponding to skydir (can be single or a list) using
            bilinear interpolation.
            
            The skydir(s) are rotated onto the equatorial grid."""
        
        if hasattr(skydir,'EQUATORIAL'): skydir = [skydir]
        """ # this block is DEPRECATED -- remove if C++ changes successful for everyone
        lon = N.asarray(map(SkyDir.l,skydir))
        lat = N.asarray(map(SkyDir.b,skydir))
        rlon = DoubleVector() ; rlat = DoubleVector()
        Background.rot_grid(rlon,rlat,lon,lat,self.center)
        lon = N.asarray(rlon) ; lat = N.asarray(rlat)
        """ # end DEPRECATED
        rvals = N.empty(len(skydir)*2,dtype=float)
        PythonUtilities.rot_grid(rvals,skydir,self.center)
        lon = rvals[::2]; lat = rvals[1::2]
        if self.wrap:
            # adopt negative longitudes for the nonce; fine for calculating differences
            lon = N.where(lon > 180, lon - 360 , lon)
        np = self.npix - 1
        x = (self.lon0 - lon) / (self.delta_lon / np)
        y = (lat - self.lat0) / (self.delta_lat / np)
        #if N.any( (x<0) || (x > np) || (y < 0) || (y > np) ):
        #    return N.nan
        xlo,ylo = N.floor(x+1e-6).astype(int),N.floor(y+1e-6).astype(int)
        xhi,yhi = xlo+1,ylo+1
        dx = N.maximum(0,x - xlo)
        dy = N.maximum(0,y - ylo)
        v  = N.asarray(v)
        if self.bounds_error: 
            return v[xlo,ylo]*(1-dx)*(1-dy) + v[xlo,yhi]*(1-dx)*dy + v[xhi,ylo]*dx*(1-dy) + v[xhi,yhi]*dx*dy
        else:
            return N.where((xlo<0) | (ylo<0) | (xhi>np) | (yhi>np),self.fill_value,
                           v[N.clip(xlo,0,np),N.clip(ylo,0,np)]*(1-dx)*(1-dy) + \
                           v[N.clip(xlo,0,np),N.clip(yhi,0,np)]*(1-dx)*dy + \
                           v[N.clip(xhi,0,np),N.clip(ylo,0,np)]*dx*(1-dy) + \
                           v[N.clip(xhi,0,np),N.clip(yhi,0,np)]*dx*dy)

    def fill(self,skyfun):
        """ Evaluate skyfun along the internal grid and return the resulting array.
        
            In this process, the internal grid is transformed to the center of
            the ROI, and the skyfun evaluated for each of the resulting
            positions.
            
            The rotation is rather fast, likely no need to pre-compute.
        """
        """ # this block is DEPRECATED -- remove if C++ changes successful for everyone
        v = DoubleVector()
        Background.val_grid(v,self.lons,self.lats,self.center,skyfun)
        return N.resize(v,[self.npix,self.npix])
        """ # end DEPRECATED
        v = N.empty(self.npix*self.npix)
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

    def do_convolution(self,energy,conversion_type,override_skyfun=None,override_en=None):
        """ Perform a convolution at the specified energy, for the specified
            conversion type.  The values are stored internally as "cvals".
        """
        if override_skyfun is None:
            self.bg.set_skyfun(conversion_type,energy)
            self.bg_vals = self.fill(self.bg)
        else:
            self.bg_vals = self.fill(override_skyfun)
        pb = PretendBand(energy,conversion_type)
        bpsf = BandCALDBPsf(self.psf,pb,override_en=override_en,adjust_mean=False)
        self.psf_fill(bpsf)
        self.convolve()


    def psf_fill(self,psf):
        """ Evaluate a band psf over the grid."""
        psf_vals = psf(self.dists,density=True).reshape([self.npix,self.npix])
        self.psf_vals = psf_vals / psf_vals.sum()
        #self.psf_vals = psf_vals*N.radians(self.pixelsize)**2

    def convolve(self):
        """ Perform the convolution with the current values of the bg
            and psf evaluated over the grid."""
        fft_kernel = fft2(self.psf_vals)
        fft_image  = fft2(self.bg_vals)
        self.cvals = c = N.real(fftshift(ifft2(fft_kernel*fft_image)))

        # swap the 0th component into the proper place
        new = N.empty_like(c)
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
        norm = P.normalize(N.log10(self.bg_vals.min()),N.log10(self.bg_vals.max()))
        P.subplot(132)
        P.imshow(N.log10(self.bg_vals).transpose()[::-1],norm=norm,interpolation='nearest')
        P.axvline(marker,color='k')
        P.axhline(marker,color='k')
        P.subplot(133)
        P.imshow(N.log10(self.cvals).transpose()[::-1],norm=norm,interpolation='nearest')
        P.axvline(marker,color='k')
        P.axhline(marker,color='k')

#===============================================================================================#

class BackgroundConvolutionNorm(BackgroundConvolution):
    """ This object is suitable for a spatial model which are supposed to be normalized to 1.
        It is also intended for spatial models where the entire convolution radius is
        inside of the total roi.
        
        Note that this implementation assumes that there is no exposure variation
        across the source and ap_average must be mulitplied by the exposure
        outside of this function."""

    def convolve(self,*args,**kwargs):
        super(BackgroundConvolutionNorm,self).convolve(*args,**kwargs)
        self.cvals /= self.cvals.sum()*N.radians(self.pixelsize)**2

#===============================================================================================#

class AnalyticConvolution(object):
    """ Calculates the convolution of the psf with a radially symmetric spatial_model. 
        This object has a similar interface to BackgroundConvolution.         

        Note that this implementation assumes that there is no exposure variation
        across the source and ap_average must be mulitplied by the exposure
        outside of this function. """

    def init(self):
        self.num_points     = 200

        # number of points to use in the intergral
        self.num_int_points = 200
        self.skyfun         = PySkyFunction(self)
        self.fitpsf         = False

    def __init__(self,spatial_model,psf,**kwargs):
        """
Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  num_points     [200]   Number of points to calculate the PDF at.
                         interpolation is done in between.
  num_int_points [200]   Number of points to use in evaluating the 
                         integral.
  fitpsf         [False] Use emperical fits to the psf instead of
                         weighting over cos theta.
  =========     =======================================================
        """

        self.init()
        self.__dict__.update(**kwargs)


        self.spatial_model=spatial_model
        self.psf=psf


    def _get_pdf(self,rlist,g,s,energy):
        """ Function to calculate the pdf at the givin radii points
            for a given psf sigma and gamma. (and a radially symmetric
            spatial model self.spatial_model

            For each givin radii, we must integrate to to the 
            maximum source radius

        integrate to the edge of the source, which is decided on a
        source by source basis (so that sharp edge sources can define
        automatically. The default is a some constant times r68 since
        the intgral includes at term which is the PDF, the integral
        will presumably contribute very littel further away then this.
        """
        self.int_max = 0.5*(self.spatial_model.effective_edge(energy)/s)**2

        # u value corresponding to the given r.
        ulist=0.5*(self.rlist/s)**2

        vlist=N.linspace(0,self.int_max,self.num_int_points)

        v,u=N.meshgrid(vlist,ulist)

        integrand = self.spatial_model.at_r(N.sqrt(2*v)*s,energy)\
                                          *((g-1)/g)*(g/(g+u+v))**g\
                                          *hyp2f1(g/2.0,(1+g)/2.0,1.0,4.0*u*v/(g+u+v)**2)
        pdf = simps(integrand,v)
        return pdf


    def do_convolution(self,energy,conversion_type,band):
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
            create something that looks like a gaussian. """
        pb = PretendBand(energy,conversion_type)
        self.bpsf = BandCALDBPsf(self.psf,pb)

        if self.bpsf.newstyle:
            nclist,ntlist,gclist,gtlist,\
                    sclist,stlist,wlist = self.bpsf.par

            smax = N.append(sclist,stlist).max()

        else:
            # g = gamma, s = sigma, w = weight
            glist,slist,wlist = self.bpsf.par
            smax=slist.max()

        # Distance to furthest pixel is the farthest one ever needs to calculate PDF to.
        # Note that this generally interfeers with ap_average. Will fix 
        # eventually.
        self.edge_distance=band.sd.difference(self.spatial_model.center) + \
                band.radius_in_rad

        self.rmax=1.1*self.edge_distance

        self.rlist=N.linspace(0,N.sqrt(self.rmax),self.num_points)**2


        # pdf is the probability per unit area at a given radius.
        self.pdf=N.zeros_like(self.rlist)

        if self.fitpsf:
            # For new & old style psf, fit a single king function to the data.
            self.pdf = self._get_pdf(self.rlist,band.fit_gamma,band.fit_sigma,energy)
        else:
            if self.bpsf.newstyle:
                for nc,gc,sc,nt,gt,st,w in zip(nclist,gclist,sclist,\
                                               ntlist,gtlist,stlist,wlist):
                    self.pdf += w*(nc*self._get_pdf(self.rlist,gc,sc,energy)+
                                   nt*self._get_pdf(self.rlist,gt,st,energy))
            else:
                for g,s,w in zip(glist,slist,wlist):
                    self.pdf += w*self._get_pdf(self.rlist,g,s,energy)

        # Assume pdf is 0 outside of the bound, which is reasonable if 
        # rmax is big enough also, interpolate the log of the intensity, which 
        # should for the types of shapes we are calculating be much more linear.
        self.log_pdf=N.log(self.pdf)
        self.interp=interp1d(self.rlist,self.log_pdf,kind=3,bounds_error=False,fill_value=self.log_pdf[-1])

        self.val=lambda x: N.exp(self.interp(x))*(x<=self.rlist[-1])

    def ap_average(self,center,radius):
        # This function needs to be fixed to really integrate for when the extended source isn't all in.
        solid_angle=2*N.pi*(1-N.cos(radius))
        return 1/solid_angle

    def __call__(self,skydir):
        if type(skydir) == WeightedSkyDirList:
            dv = DoubleVector()
            skydir.arclength(self.spatial_model.center,dv)
            difference = N.fromiter(dv,dtype=float)
            return self.val(difference)

        elif type(skydir)==list and len(skydir)==3:
            skydir = SkyDir(Hep3Vector(skydir[0],skydir[1],skydir[2]))

        elif type(skydir)==SkyDir:
            return float(self.val(skydir.difference(self.spatial_model.center)))
        else:
            raise Exception("Unknown input to AnalyticConvolution.__call__()")


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
