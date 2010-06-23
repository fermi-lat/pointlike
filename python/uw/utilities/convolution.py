"""Module to support on-the-fly convolution of a mapcube for use in spectral fitting.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/utilities/convolution.py,v 1.7 2010/06/22 21:55:43 lande Exp $

authors: M. Kerr, J. Lande

"""
from skymaps import SkyDir,WeightedSkyDirList,Hep3Vector,SkyIntegrator,PySkyFunction,Background
from pointlike import DoubleVector
import numpy as N
from scipy.interpolate import interp2d,interp1d
from scipy.integrate import quad,Inf
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

    def __init__(self,center,npix=4*25+1,pixelsize=1./4):
        """ center -- a SkyDir giving the center of the grid
            npix   -- the number of pixels along a grid side
            pixelsize -- the length of the side of a pixel (deg)
        """

        self.npix,self.pixelsize = npix,pixelsize
        
        l,b    = center.l(),center.b()
        self.center = center
        
        self.clon,self.clat = l,0
        self.delta_lon = (npix-1)*pixelsize
        self.delta_lat = (npix-1)*pixelsize
        self.lon0 = (self.clon + self.delta_lon/2.)%360
        self.lat0 = self.clat - self.delta_lat/2.
        self.wrap = self.lon0 <= self.delta_lon # True if origin (0/360) is in frame

        npix  = self.npix
        mpix  = (float(self.npix)-1)/2.
        v     = (N.arange(0,npix) - mpix)**2
        self.dists = N.asarray([v+v[i] for i in xrange(npix)])**0.5*(N.radians(self.pixelsize))

        self.lons = self.lon0 - N.linspace(0,1,self.npix)*self.delta_lon
        self.lons = DoubleVector(N.where(self.lons < 0,self.lons+360,self.lons)) # if origin in frame
        self.lats = DoubleVector(N.linspace(0,1,self.npix)*self.delta_lat + self.lat0)

        self.rot_axis = SkyDir(l+90,0).dir() # Hep3Vector

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
            find the value(s) corresponding to skydir (can be single ora list) using
            bilinear interpolation.
            
            The skydir(s) are rotated onto the equatorial grid."""
        
        if hasattr(skydir,'EQUATORIAL'): skydir = [skydir]
        lon = N.asarray(map(SkyDir.l,skydir))
        lat = N.asarray(map(SkyDir.b,skydir))
        rlon = DoubleVector() ; rlat = DoubleVector()
        Background.rot_grid(rlon,rlat,lon,lat,self.center)
        lon = N.asarray(rlon) ; lat = N.asarray(rlat)
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
        return v[xlo,ylo]*(1-dx)*(1-dy) + v[xlo,yhi]*(1-dx)*dy + v[xhi,ylo]*dx*(1-dy) + v[xhi,yhi]*dx*dy

    def fill(self,skyfun):
        """ Evaluate skyfun along the internal grid and return the resulting array.
        
            In this process, the internal grid is transformed to the center of
            the ROI, and the skyfun evaluated for each of the resulting
            positions.
            
            ***THIS SHOULD BE PROFILED AT SOME POINT TO SEE IF WE SHOULD
            ***PRE-CALCULATE THE TRANSFORMATION.
        """
        v = DoubleVector()
        Background.val_grid(v,self.lons,self.lats,self.center,skyfun)
        return N.resize(v,[self.npix,self.npix])

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

    def do_convolution(self,energy,conversion_type,override_skyfun=None):
        """ Perform a convolution at the specified energy, for the specified
            conversion type.  The values are stored internally as "cvals".
        """
        if override_skyfun is None:
            self.bg.setEnergy(energy)
            self.bg.set_event_class(conversion_type)
            self.bg_vals = self.fill(self.bg)
        else:
            self.bg_vals = self.fill(override_skyfun)
        pb = PretendBand(energy,conversion_type)
        bpsf = BandCALDBPsf(self.psf,pb)
        self.psf_fill(bpsf)
        self.convolve()


    def psf_fill(self,psf):
        """ Evaluate a band psf over the grid."""
        psf_vals = psf(self.dists).reshape([self.npix,self.npix])
        self.psf_vals = psf_vals / psf_vals.sum()

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

    def convolve(self,*args,**kwargs):
        super(BackgroundConvolutionNorm,self).convolve(*args,**kwargs)
        self.cvals /= self.cvals.sum()

#===============================================================================================#

class AnalyticConvolution(object):
    """ Calculates the convolution of the psf with a radially symmetric spatial_model. 
        Has a very similar interface to BackgroundConvolution. Maybe one day they will
        derive from the same base convolution object. For now they are different enough."""

    def init(self):
        self.weight_before  =  False
        self.tolerance      =  1e-3
        self.num_points     =  50
        self.skyfun         = PySkyFunction(self)

    def __init__(self,spatial_model,psf,**kwargs):
        """
Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  weight_before [False] Weight sigmas & gammas before convolving.
  tolerance     [1e-4]  Tolerance set when integrating hypergeometirc
                        function to calcualte PDF.
  num_points    [100]   Number of points to calculate the PDF at.
                        interpolation is done in between.
  =========   =======================================================
        """

        self.init()
        self.__dict__.update(**kwargs)


        self.spatial_model=spatial_model
        self.psf=psf


    def _get_pdf(self,rlist,g,s):
        # integrate until you get to 10x the r68() of the source.
        # since the intgral includes at term which is the PDF,
        # the integral will presumably contribute very littel
        # further away then this.
        self.int_max = .5*(10*self.spatial_model.r68()/s)**2

        # u value corresponding to the given r.
        ulist=0.5*(self.rlist/s)**2

        # pdf values for each r (in a given energy bin)
        pdf =N.empty_like(self.rlist)

        for i,u in enumerate(ulist):
            integrand = lambda v: self.spatial_model.at_r(N.sqrt(2*v)*s)\
                                              *((g-1)/g)*(g/(g+u+v))**g\
                                              *hyp2f1(g/2.,(1+g)/2.,1.,4.*u*v/(g+u+v)**2)

            pdf[i]=quad(integrand,0,self.int_max,epsabs=self.tolerance,full_output=True)[0]

        return pdf 

    def do_convolution(self,energy,conversion_type):
        pb = PretendBand(energy,conversion_type)
        self.bpsf = BandCALDBPsf(self.psf,pb)

        if self.bpsf.newstyle:
            nclist,ntlist,gclist,gtlist,\
                    sclist,stlist,wlist = self.bpsf.par

            self.rmax=5*(N.append(sclist,stlist).max() +
                          self.spatial_model.r68())

            self.rlist=N.linspace(0,self.rmax,self.num_points)

            # pdf is the probability per unit area at a given radius.
            self.pdf=N.zeros_like(self.rlist)

            if self.weight_before:
                nc=(nclist*wlist).sum()
                nt=(ntlist*wlist).sum()
                gc=(gclist*wlist).sum()
                gt=(gtlist*wlist).sum()
                sc=(sclist*wlist).sum()
                st=(stlist*wlist).sum()

                self.pdf = nc*self._get_pdf(self.rlist,gc,sc)+\
                        nt*self._get_pdf(self.rlist,gt,st)
            else:
                for nc,gc,sc,nt,gt,st,w in zip(nclist,gclist,sclist,\
                                               ntlist,gtlist,stlist,wlist):
                    self.pdf += w*(nc*self._get_pdf(self.rlist,gc,sc)+
                                   nt*self._get_pdf(self.rlist,gt,st))

        else:
            # g = gamma, s = sigma, w = weight
            glist,slist,wlist = self.bpsf.par

            self.rmax=5*(slist.max()+self.spatial_model.r68())

            self.rlist=N.linspace(0,self.rmax,self.num_points)

            # pdf is the probability per unit area at a given radius.
            self.pdf=N.zeros_like(self.rlist)

            if self.weight_before:
                # weight the convolved shape by gamma, sigma, and weight
                g=(glist*wlist).sum()
                s=(slist*wlist).sum()
                self.pdf = self._get_pdf(self.rlist,g,s)
            else:
                # weight the convolved shape by gamma, sigma, and weight
                for g,s,w in zip(glist,slist,wlist):
                    self.pdf += w*self._get_pdf(self.rlist,g,s)

        # for some reason, I incorrectly got a negative pdf value for r=0 and for
        # especially large gaussian MC sources. Not sure how the integral of the 
        # hypergeometric function could do this, but it is best to simply remove 
        # it from the list since we know they are unphysical.
        bad = N.isnan(self.pdf)|N.isinf(self.pdf)|(self.pdf<0)

        if N.any(bad):
            print 'WARNING! Bad values found in PDF. Removing them from interpolation.' % sum(bad),
            if N.any(N.isnan(self.pdf)): print ' (%d nan values)' % sum(N.isnan(self.pdf)),
            if N.any(N.isinf(self.pdf)): print ' (%d inf values)' % sum(N.isinf(self.pdf)),
            if N.any(self.pdf<0):        print ' (%d negative values)' % sum(self.pdf<0),
            print
            self.rlist=self.rlist[~bad]
            self.pdf=self.pdf[~bad]

        # Assume pdf is 0 outside of the bound, which is reasonable if 
        # rmax is big enough also, interpolate the log of the intensity, which 
        # should for the types of shapes we are calculating be much more linear.
        self.log_pdf=N.log(self.pdf)
        self.interp=interp1d(self.rlist,self.log_pdf,kind='cubic',bounds_error=False,fill_value=0)

        # Just in case the first entry is a NaN, set values below equal to the first non-nan/inf
        # pdf value.
        self.val=lambda x: N.exp(self.interp(x))*(x>=self.rlist[0])+self.pdf[0]*(x<self.rlist[0])

    def ap_average(self,center,radius):
        solid_angle=2*N.pi*(1-N.cos(radius))
        return 1/solid_angle

    def __call__(self,skydir,not_needed=None):
        if type(skydir)==WeightedSkyDirList:
            dv = DoubleVector()
            skydir.arclength(self.spatial_model.center,dv)
            difference = N.fromiter(dv,dtype=float)
            return self.val(difference)

        if type(skydir)==list and len(skydir)==3:
            skydir = SkyDir(Hep3Vector(skydir[0],skydir[1],skydir[2]))

        if type(skydir)==SkyDir:
            return float(self.val(skydir.difference(self.spatial_model.center)))
        elif type(skydir)== WeightedSkyDirList:
            dv = DoubleVector()
            skydir.arclength(self.spatial_model.center,dv)
            difference = N.fromiter(dv,dtype=float)
            return self.val(difference)
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
