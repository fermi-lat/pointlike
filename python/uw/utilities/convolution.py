"""Module to support on-the-fly convolution of a mapcube for use in spectral fitting.
$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/utilities/convolution.py,v 1.3 2010/05/24 17:32:10 wallacee Exp $
"""
from skymaps import SkyDir,WeightedSkyDirList,Hep3Vector,SkyIntegrator,PySkyFunction
from pointlike import DoubleVector
import numpy as N
from scipy.interpolate import interp2d,interp1d
from scipy.integrate import quad,Inf
from scipy.special import hyp2f1
from uw.like.pypsf import BandCALDBPsf,PretendBand
from numpy.fft import fftshift,ifft2,fft2

class Grid(object):
    """ Define a grid with vertices separated by a uniform arclength, i.e., flat.

        The default coordinate system is Galactic.  When the center is placed too
        close to the poll, algorithm will switch to celestial.

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
        usegal = self.usegal = abs(b) < 70
        
        self.clon,self.clat = (center.l(),center.b()) if usegal else (center.ra(),center.dec())
        self.delta_lon = (npix-1)*pixelsize/N.cos(N.radians(self.clat))
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

    def in_grid(self,lon,lat,skydir=None):
        """ check whether the provided lon/lat is within the grid. """
        if skydir is not None:
            lon,lat = (skydir.l(),skydir.b()) if self.usegal else (skydir.ra(),skydir.dec())
        delta_lon = abs(lon - self.clon)
        delta_lon = min(360-delta_lon,delta_lon)
        lon_ok = delta_lon < self.delta_lon/2.
        lat_ok = abs(lat - self.clat) < self.delta_lat/2.

        return lat_ok and lon_ok

    def pix(self,skydir):
        lon,lat = (skydir.l(),skydir.b()) if self.usegal else (skydir.ra(),skydir.dec())
        if self.wrap:
            lon = N.where(lon > 180, lon - 360 , lon)
        np = self.npix - 1
        x = (self.lon0 - lon)/self.delta_lon*np
        y = (lat - self.lat0)/self.delta_lat*np
        return x,y
    
    def __call__(self,skydir,v):
        """ Using v, an array of values that has been evaluated on the Grid (e.g., by fill),
            find the value(s) corresponding to skydir (can be single ora list) using
            bilinear interpolation."""
        if hasattr(skydir,'EQUATORIAL'):
            skydir = [skydir]
            #lon,lat = (skydir.l(),skydir.b()) if self.usegal else (skydir.ra(),skydir.dec())
        #else:
        lon = N.asarray([sd.l() if self.usegal else sd.ra() for sd in skydir])
        lat = N.asarray([sd.b() if self.usegal else sd.dec() for sd in skydir])
        if self.wrap:
            # adopt negative longitudes for the nonce; fine for calculating differences
            lon = N.where(lon > 180, lon - 360 , lon)
        np = self.npix - 1
        x = (self.lon0 - lon)/self.delta_lon*np
        y = (lat - self.lat0)/self.delta_lat*np
        #if N.any( (x<0) || (x > np) || (y < 0) || (y > np) ):
        #    return N.nan
        xlo,ylo = N.floor(x+1e-6).astype(int),N.floor(y+1e-6).astype(int)
        xhi,yhi = xlo+1,ylo+1
        dx = N.maximum(0,x - xlo)
        dy = N.maximum(0,y - ylo)
        #print x,xlo,xhi,dx
        #print y,ylo,yhi,dy
        v = N.asarray(v)
        return v[xlo,ylo]*(1-dx)*(1-dy) + v[xlo,yhi]*(1-dx)*dy + v[xhi,ylo]*dx*(1-dy) + v[xhi,yhi]*dx*dy

    def fill(self,skyfun):
        """Evaluate skyfun along the internal grid and return the resulting array."""
        lons,lats = self.lons,self.lats
        v = N.empty([len(lons),len(lats)])
        s = [SkyDir.GALACTIC if self.usegal else SkyDir.EQUATORIAL] * len(lons)
        for ilon,lon in enumerate(lons):
            sds = map(SkyDir,[lon]*len(lons),lats,s)
            v[ilon,:] = [skyfun(sd) for sd in sds]
        return v

class BackgroundConvolution(Grid):

    def __init__(self,center,bg,psf,*args,**kwargs):
        """ center -- a SkyDir giving the center of the grid on which to convolve bg
            bg     -- an instance of skymaps::Background
            psf    -- and instance of CALDBPsf

            Additional kwargs are passed to Grid.
        """

        self.bg,self.psf = bg,psf
        super(BackgroundConvolution,self).__init__(center,**kwargs)

    def do_convolution(self,energy,conversion_type,override_skyfun=None):
        """ Perform a convolution at the specified energy, for the specified
            conversion type.  The values are stored internally as "cvals".
        """
        if override_skyfun is None:
            self.bg.setEnergy(energy)
            self.bg.set_event_class(conversion_type)
            self.bg_fill()
        else:
            self.bg_vals = self.fill(override_skyfun)
        pb = PretendBand(energy,conversion_type)
        bpsf = BandCALDBPsf(self.psf,pb)
        self.psf_fill(bpsf)
        self.convolve()

    def bg_fill(self):
        """ Evaluate the internal bg member over the grid."""
        lons,lats = self.lons,self.lats
        rval = DoubleVector()
        self.bg.grid_values(rval,lons,lats,SkyDir.GALACTIC if self.usegal else SkyDir.EQUATORIAL)
        self.bg_vals = N.resize(rval,[self.npix,self.npix])

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
        derive from a base convolution object."""

    def init(self):
        self.fast       = False                 # Weight sigmas & gammas before convolving.
        self.num_points = 50                    # Number of points between 0 and rmax.
        self.tolerance  = 0.02                  # Tolerance set when integrating the function
        self.skyfun     = PySkyFunction(self)

    def __init__(self,spatial_model,psf,**kwargs):
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

            pdf[i]=quad(integrand,0,self.int_max,epsabs=1e-3,full_output=True)[0]

        return pdf 

    def do_convolution(self,energy,conversion_type):
        pb = PretendBand(energy,conversion_type)
        self.bpsf = BandCALDBPsf(self.psf,pb)

        if self.bpsf.newstyle:
            # nc,nt,gc,gt,sc,st,w = band.bpsf.par
            raise Exception("AnalyticPDF Not implemented for newstyle psf.")

        else:
            # g = gamma, s = sigma, w = weight
            glist,slist,wlist = self.bpsf.par

            self.rmax=10*(slist.max()+self.spatial_model.r68())

            self.rlist=N.linspace(0,self.rmax,self.num_points)

            # pdf is the probability per unit area at a given radius.
            self.pdf=N.zeros_like(self.rlist)

            if self.fast:
                # weight the convolved shape by gamma, sigma, and weight
                g=(glist*wlist).sum()
                s=(slist*wlist).sum()
                self.pdf = self._get_pdf(self.rlist,g,s)
            else:
                # weight the convolved shape by gamma, sigma, and weight
                for g,s,w in zip(glist,slist,wlist):
                    self.pdf += w*self._get_pdf(self.rlist,g,s)

            # Assume pdf is 0 outside of the bound,
            # reasonable if rmax is big enough
            self.interp=interp1d(self.rlist,self.pdf,kind='cubic',bounds_error=False,fill_value=0)

    def ap_average(self,center,radius):
        solid_angle=2*N.pi*(1-N.cos(radius))
        if center.difference(self.spatial_model.center)+self.rmax < radius:
            # the source is properly normalized, so there is only 
            # a need to calculate the contribution if the PDF spills
            # outside of the ROI.
            return 1/solid_angle
        else:
            SkyIntegrator.set_tolerance(self.tolerance)
            return SkyIntegrator.average(self.skyfun,center,radius)

    def __call__(self,skydir,not_needed=None):
        if type(skydir)==WeightedSkyDirList:
            dv = DoubleVector()
            skydir.arclength(self.spatial_model.center,dv)
            difference = N.fromiter(dv,dtype=float)
            return self.interp(difference)

        if type(skydir)==list and len(skydir)==3:
            skydir = SkyDir(Hep3Vector(skydir[0],skydir[1],skydir[2]))

        if type(skydir)==SkyDir:
            return float(self.interp(skydir.difference(self.spatial_model.center)))

        raise Exception("Unknown input to AnalyticConvolution.__call__()")
