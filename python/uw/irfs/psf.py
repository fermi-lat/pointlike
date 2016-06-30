"""Module providing handling of the LAT point spread function.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/irfs/psf.py,v 1.2 2016/06/29 18:05:34 wallacee Exp $
Author: Eric Wallace
"""
__version__='$Revision: 1.2 $'

import os

import numpy as np
from astropy.io import fits
from scipy import integrate

from uw.utilities import keyword_options
from . import caldb, IrfError

class PSF(object):
    """Object representing the LAT PSF."""
    
    def __init__(self,filename,exposure=None,rpsf_extension='RPSF',psf_scaling_extension='PSF_SCALING'):
        self._load_data(filename,rpsf_extension,psf_scaling_extension)
        self.set_weights(exposure)

    def _load_data(self,filename,rpsf,scaling):
        rpsf = fits.getdata(filename,rpsf)[0]
        psf_scaling = fits.getdata(filename,scaling)[0]
        self.ebins = np.vstack([rpsf.field('ENERG_LO'),rpsf.field('ENERG_HI')]).T
        self.cthetabins = np.vstack([rpsf.field('CTHETA_LO'),rpsf.field('CTHETA_HI')]).T
        self.scale_parameters = psf_scaling.field('PSFSCALE')
        def _normalize():
            #Normalize parameters. Logic copied from like.pypsf.
            #I'm still confused about some of the conventions -EEW
            sf = self.scale_function(self.ebins.prod(axis=1)**.5)
            normc = self.psf_base_integral(np.pi/2,self.score*sf,self.gcore)
            normt = self.psf_base_integral(np.pi/2,self.score*sf,self.gtail)
            # NB leave scale factor out here so we can adjust norm to
            # a particular energy (effectively cancelling in integral)
            norm = (2*np.pi*(normc*self.score**2+
                             normt*self.ntail*self.stail**2))**-1
            self.ncore = norm # adjust NCORE
            self.ntail *= norm # set to NTAIL*NCORE

        if 'SCORE' in rpsf.array.dtype.names:
            self.ncore = rpsf.field('NCORE')
            self.score = rpsf.field('SCORE')
            self.gcore = rpsf.field('GCORE')
            self.ntail = rpsf.field('NTAIL')
            self.stail = rpsf.field('STAIL')
            self.gtail = rpsf.field('GTAIL')
            _normalize()
        else:
            #Old style (P6)
            #TODO: Check that this ncore gives the correct normalization
            self.score = self.stail = rpsf.field('SIGMA')
            self.gcore = rpsf.field('GCORE')
            self.gtail = rpsf.field('GTAIL')
            self.ncore = np.ones_like(rpsf.field('SIGMA'))
            self.ntail = np.zeros_like(rpsf.field('SIGMA'))
    
    def __getitem__(self,mask):
        """Return PSF parameters for a given energy and cos(theta) selection"""
        return np.concatenate([getattr(self,p)[mask][None]
                    for p in ('ncore','ntail','score','stail','gcore','gtail')])
        
    def set_weights(self,exposure=None):
        """Set weights to use for exposure-weighted averages."""
        self.weights = None
        #if exposure is None:
        #    self.weights = None
        #else:
        #    self.weights = 

    def scale_function(self,e):
        """Compute the PSF scale factor for energy `e`.
        
        Parameter
        ---------

        e 
            Energy to evaluate the scale factor at, in MeV. May be a scalar
            or a numpy array.  

        Returns
        -------

        sp 
            The scale factor at the requested energy, defined as
            :math:`\sqrt{\left(c_0(\frac{E}{100 MeV})^{-\beta}\right)^2+c_1^2}`.
            The return type is the same as that of `e`.
        """

        c0,c1,beta = self.scale_parameters
        
        return np.sqrt( (c0*(np.asarray(e)/100.)**-beta)**2 + c1**2)

    def __call__(self,delta,e):
        """Compute the PSF density at angular deviation delta and energy e.

        If exposure weights have been set by the `set_weights` method, the
        calculated density will be a weighted average over inclination angle.
        Otherwise, only the on-axis value will be calculated.

        Parameters
        ----------
        delta : array_like
            Angular deviation(s) at which to evaluate the PSF, in radians.
        e : array_like
            Energy (or energies) at which to evaluate the PSF, in MeV.
        ctheta : array_like, optional
            Cosine of the inclination angle, theta, for which the PSF is
            to be evaluated. The default is None, indicating that the 
            computed density should be averaged over the inclination angle.
            
            If exposure information has been provided through the 
            `set_exposure` method, the average over theta will be weighted
            by the exposure.

        Returns
        -------
        density: float or array
            The PSF density at angular distance `delta` and energy `e`. If
            either parameter is an array, the return value is an array of
            shape (n_delta, n_e). 
        
        """
        scale = self.scale_function(e)

        mask = np.fmin(np.searchsorted(self.ebins[:,1],e),self.ebins.shape[0]-1)
        if self.weights is None:
            nc,nt,sc,st,gc,gt = self[-1,mask]
            kc,kt = [self.psf_base(delta,s*scale,g)
                 for s,g in zip((sc,st),(gc,gt))]
            return (nc*kc+nt*kt)/scale**2
        else:
            nc,nt,sc,st,gc,gt = self[:,mask]
            kc,kt = [self.psf_base(delta,s*scale,g)
                 for s,g in zip((sc,st),(gc,gt))]
            return (self.weights[:,mask]*(nc*kc+nt*kt)/scale**2).sum(axis=-2)


    def psf_base(self,delta,sigma,gamma):
        """Evaluate the King function at angular deviation delta.

        Parameters
        ----------

        delta : array_like
            The angular deviation in radians at which the function is to be 
            evaluated. May be a scalar or a numpy array.
        sigma, gamma : array_like
            The parameters of the King function. May be scalars or arrays of
            the same size.

        Returns
        -------

        psf_base : float or array
            If `delta`, `sigma`, and `gamma` are all scalars, a scalar is
            returned. Otherwise, the return value is an array of dimension
            len(`delta`) by len(`sigma`).

        """
        return_scalar = np.all([np.isscalar(x) for x in (delta,sigma,gamma)])
        d,s,g = (np.asarray(x) for x in (delta,sigma,gamma))
        if s.shape!=g.shape:
            raise ValueError('Arrays for sigma and gamma must have the same shape')
        u = (.5*np.outer(d,1/s)**2).reshape(d.shape+s.shape)
        k = (1-1/g)*(1+u/g)**-g
        if return_scalar:
            return k.item()
        else:
            return k

    def psf_base_integral(self,dmax,sigma,gamma,dmin=0):
        """Integral of the PSF base function; g = gamma, s = sigma (scaled),
           delta = deviation in radians."""
        return_scalar = np.all([np.isscalar(x) for x in (dmax,sigma,gamma,dmin)])
        dmax,s,g,dmin = (np.asarray(x) for x in (dmax,sigma,gamma,dmin))
        if s.shape!=g.shape:
            raise ValueError('Arrays for sigma and gamma must have the same shape')
        if (dmin>0) and (dmax.shape!=dmin.shape):
            raise ValueError('Arrays for dmin and dmax must have the same shape')
        u0 = (.5*np.outer(dmin,1/s)**2).reshape(dmin.shape+s.shape)
        u1 = (.5*np.outer(dmax,1/s)**2).reshape(dmax.shape+s.shape)
        i = (1+u0/g)**(1-g)-(1+u1/g)**(1-g)
        if return_scalar:
            return i.item()
        else:
            return i

    def integral(self,dmax,e,dmin=0):
        """Integral of the PSF from dmin to dmax at energy e."""
        scale = self.scale_function(e)
        mask = np.fmin(np.searchsorted(self.ebins[:,1],e),self.ebins.shape[0]-1)
        if self.weights is None:
            nc,nt,sc,st,gc,gt = self[-1,mask]
            icore = np.pi*2*sc**2*nc*self.psf_base_integral(dmax, sc*scale, gc, dmin)
            itail = np.pi*2*st**2*nt*self.psf_base_integral(dmax, st*scale, gt, dmin)
            return (icore+itail)
        else:
            nc,nt,sc,st,gc,gt = self[:,mask]
            icore = np.pi*2*sc**2*nc*self.psf_base_integral(dmax, sc*scale, gc, dmin)
            itail = np.pi*2*st**2*nt*self.psf_base_integral(dmax, st*scale, gt, dmin)
            return (self.weights[:,mask]*(icore+itail)).sum(axis=-2)
    
        
    def band_psf(self,energy):
        return BandPSF(self,energy)
    

class BandPSF(PSF):
    """Representation of the PSF for a specific energy band."""

    def __init__(self,psf,energy):
        if hasattr(energy,'__iter__'):
            raise ValueError('BandPSF can only be defined for a single energy.') 
        self.energy = energy

        ind = min(np.searchsorted(psf.ebins[:,1],energy),psf.ebins.shape[0]-1)
        self._scale = psf.scale_function(energy)
        for p in ('ncore','ntail','score','stail','gcore','gtail'):
            setattr(self,p,getattr(psf,p)[:,ind])
        if psf.weights is not None:
            self.weights = psf.weights[:,ind]
        else:
            self.weights = None

    def scale_function(self):
        return self._scale

    def set_weights(self,exposure=None):
        if exposure is not None:
            PSF.set_weights(self,exposure)
            emask = min(np.searchsorted(psf.ebins[:,1],self.energy),psf.ebins.shape[0]-1)
            self.weights = self.weights[:,emask]
        else:
            self.weights = None
        
    def __call__(self,delta):
        scale = self.scale_function()
        if self.weights is None:
            nc,nt,sc,st,gc,gt = self[-1]
            kc,kt = [self.psf_base(delta,s*scale,g)
                 for s,g in zip((sc,st),(gc,gt))]
            return (nc*kc+nt*kt)/scale**2
        else:
            nc,nt,sc,st,gc,gt = self[:,]
            kc,kt = [self.psf_base(delta,s*scale,g)
                 for s,g in zip((sc,st),(gc,gt))]
            return (self.weights[:,mask]*(nc*kc+nt*kt)/scale**2).sum(axis=-2)
    
    def integral(self,dmax,dmin=0):
        """Integral of the PSF from dmin to dmax at energy e."""
        scale = self.scale_function()
        if self.weights is None:
            nc,nt,sc,st,gc,gt = self[-1]
            icore = np.pi*2*sc**2*nc*self.psf_base_integral(dmax, sc*scale, gc, dmin)
            itail = np.pi*2*st**2*nt*self.psf_base_integral(dmax, st*scale, gt, dmin)
            return (icore+itail)
        else:
            nc,nt,sc,st,gc,gt = self[:]
            icore = np.pi*2*sc**2*nc*self.psf_base_integral(dmax, sc*scale, gc, dmin)
            itail = np.pi*2*st**2*nt*self.psf_base_integral(dmax, st*scale, gt, dmin)
            return (self.weights[:,mask]*(icore+itail)).sum(axis=-2)

    def overlap(self, roi_dir, radius, skydir):
        """Calculate the fractional PSF overlap with a circle."""
        #NOTE: radius in degrees currently. Seems preferable for it to be radians for
        #consistency with other PSF code, but would require changes to clients
        radius = np.radians(radius)
        if hasattr(skydir,'__iter__'):
            scalar = False
            offset = np.asarray([roi_dir.difference(sd) for sd in skydir])
        else:
            scalar = True
            offset = np.asarray([roi_dir.difference(skydir)])
        ret = np.zeros(offset.shape)
        interior = offset<=radius
        lim = np.arcsin(radius/offset)
        ret[interior] = np.array([integrate.quad(self._interior_integrand(o,radius),0,np.pi)[0]/np.pi for o in offset[interior]])
        ret[~interior] = np.array([integrate.quad(self._exterior_integrand(o,radius),0,l)[0]/np.pi for o,l in zip(offset[~interior],lim[~interior])])
        if scalar:
            return ret.item()
        else:
            return ret

    def _interior_integrand(self,offset,radius):
        def integrand(theta):
            ctheta = np.cos(theta)
            rmax = (radius**2+offset**2*(ctheta**2-1))**.5 - ctheta*offset
            return self.integral(rmax)
        return integrand

    def _exterior_integrand(self,offset,radius):
        def integrand(theta):
            ctheta = np.cos(theta)
            x = (radius**2+offset**2*(ctheta**2-1))**.5
            return self.integral(offset*ctheta+x,offset*ctheta-x)
        return integrand

