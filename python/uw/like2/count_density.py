"""
Create observed count density files

$Header$

author:  Toby Burnett
"""
import numpy as np
import pandas as pd
from astropy.io import fits
import healpy

from .pub import healpix_map #make this local in future
from . import (diffuse, convolution )


class CountDensity(object):
    
    def __init__(self, filename, irfman, event_type, newfilename=None):
        """
        """
        print 'loading filename {}'.format(filename)
        self.hpcube= diffuse.HealpixCube(filename)
        self.hpcube.load()
        self.irfman = irfman
        self.event_type= event_type

        if newfilename is not None:
            self.replace_data()
            self.writeto(newfilename)
 
    def convert_layer(self, index):
        cube = self.hpcube
        energy = cube.energies[index]
        print 'converting layer for energy {:.0f}'.format(energy)
        cube.setEnergy(energy)
        iem_map = cube.column(energy)

        exp_map = exposure_map(self.irfman, 
            self.event_type, energy, nside=cube.nside )
        
        psf = self.irfman.psf(self.event_type, energy)
        
        count_density_map=convolve_healpix(iem_map * exp_map, psf, 5*psf.r68)
        return count_density_map

    def convert_all(self):
        return np.array([self.convert_layer(i) for i in range(len(self.hpcube.energies))])

    def replace_data(self):
        tt = self.convert_all()
        data =self.hpcube.hdulist[1].data
        for i in range(len(data)):
            data[i] = (tt[:,i],) #sort of convert to tuple with one element

    def writeto(self, newfilename):
        hdus = self.hpcube.hdulist
        header =hdus[0].header 
        header.add_comment('converted to convolved, integrated count density')
        header.add_comment('Input file {}'.format(hdus.filename()))
        header.add_comment('live-time cube {}'.format(self.irfman.ltcube))
        header.add_comment('Event type {}'.format(self.event_type))
        print 'Writing to file {}'.format(newfilename)
        hdus.writeto(newfilename)
        
def convolve_healpix(input_map,  func, thetamax, quiet=True ):
    """
    Convolve a HEALPix map with a function

    input_map : array of float
        a HEALPix array, RING indexing, nside a power of 2
    func : function of theta to convolve with
    thetamax : float
        estimate of largest theta for function

    Returns: the convolved map
    """
    nside = int(np.sqrt(len(input_map)/12))
    assert 12*nside**2 == len(input_map),'Bad length'

    alm = healpy.map2alm(input_map);
    lmax = healpy.Alm.getlmax(len(alm))
    if lmax < 0:
        raise TypeError('Wrong alm size for the given '
                        'mmax (len(alms[%d]) = %d).'%(ialm, len(alm)))
    ell = np.arange(lmax + 1.)
    if not quiet:
        print 'Creating spherical harmonic content, lmax={}'.format(lmax)
    fact = convolution.SphericalHarmonicContent(func,lmax, thetamax, quiet=quiet)(ell)

    healpy.almxfl(alm, fact, inplace=True)
    return healpy.alm2map(alm, nside=nside, verbose=False)

def exposure_map(irfman, event_type, energy, nside=128):
    """
    Return an exposure map for the given event type and entry
    """
    exp = irfman.exposure(event_type, energy)
    return healpix_map.HPskyfun('exposure', exp, nside=nside).getcol()

    