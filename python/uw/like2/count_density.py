"""
Create observed count density files

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/count_density.py,v 1.1 2017/11/17 22:50:36 burnett Exp $

author:  Toby Burnett
"""
import numpy as np
import pandas as pd
from astropy.io import fits
import os, healpy, argparse

from .pub import healpix_map #make this local in future
from . import (diffuse, convolution, configuration )


class CountDensity(object):
    """Manage computation of expected count density from a diffuse map.
    Includes exposure and convolution with PSF
    """
    
    def __init__(self, filename, irfman, event_type, energies=None, newfilename=None, overwrite=True):
        """
        filename : string
            name of a file in HEALPix format
        irfman :
        event_type : int
            0 or 1 for Front, Back

        energies : list of float | None
            if None, use the layers in the diffuse 
        """
        print 'loading filename {}, convolving with event type {}'.format(filename, event_type)
        self.hpcube= diffuse.HealpixCube(filename)
        self.hpcube.load()
        self.nside = self.hpcube.nside
        self.irfman = irfman
        self.event_type= event_type
        self.energies = self.hpcube.energies if energies is None else energies

        if newfilename is not None:
            self.replace_skymaps()
            self.replace_energies()
            self.writeto(newfilename, overwrite=overwrite)
 
    def convert_layer(self, index):
        cube = self.hpcube
        energy = self.energies[index]
        print 'converting layer for energy {:.0f}'.format(energy)
        cube.setEnergy(energy)
        iem_map = cube.column(energy) #note will interpolate

        exp_map = exposure_map(self.irfman, 
            self.event_type, energy, self.nside )
        
        psf = self.irfman.psf(self.event_type, energy)
        
        count_density_map=convolve_healpix(iem_map * exp_map, psf, 5*psf.r68)
        return count_density_map

    def replace_skymaps(self):
        tt = np.array([self.convert_layer(i) for i in range(len(self.energies))])
        newcols=[fits.Column(name='Bin{:d}'.format(i), format='E', array=tt[i]) for i,t in enumerate(tt)]
        skymap_hdu = fits.BinTableHDU.from_columns(newcols)

        skymap_hdu.header.update(
            LASTPIX = 12*self.nside**2-1,                                                  
            EXTNAME = 'SKYMAP2 ',                                                            
            PIXTYPE = 'HEALPIX ' ,                                                           
            ORDERING= 'RING    ' ,                                                           
            NSIDE   = self.nside,                                                  
            INDXSCHM= 'IMPLICIT',                                                            
            OBJECT  = 'FullSky ',                                                            
            NBRBINS = len(self.energies),                                                  
            COORDSYS= 'GAL     ',
        )
        self.hpcube.hdulist[1] = skymap_hdu

    def replace_energies(self):
        energy_hdu = fits.BinTableHDU.from_columns([fits.Column(name='energies', format='D',
             unit='MeV', array=self.energies)], name='ENERGIES')
        self.hpcube.hdulist[2]=energy_hdu

    def writeto(self, newfilename, overwrite=True):
        hdus = self.hpcube.hdulist
        header =hdus[0].header 
        header.add_comment('converted to convolved, integrated count density')
        header.add_comment('Input file {}'.format(hdus.filename()))
        header.add_comment('live-time cube {}'.format(self.irfman.ltcube))
        header.add_comment('Event type {}'.format(self.event_type))
        print 'Writing to file {}'.format(newfilename)
        hdus.writeto(newfilename, overwrite=overwrite)
        
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


def main(factor=1.0, suffix=''):
    """
    factor: float
        adjust energies
    """
    config= configuration.Configuration('.', quiet=True, postpone=True)
    iem_file = config.diffuse['ring']['filename']
    print 'IEM file: {}'.format(iem_file)
    print 'Livetime cube: {}'.format(config.irfs.ltcube)
    bf = fits.open(config.dataset.binfile)
    bands = bf['BANDS']
    emin,emax = bands.data['E_MIN'], bands.data['E_MAX']
    egmean = np.sqrt(emin*emax)
    if factor!=1.0: egmean *=factor
    energies = sorted(list(set(egmean*1e-3)))

    fname = iem_file.split('/')[-1]; 
    for i,fb in enumerate(config.event_type_names):
        
        newfile = os.path.expandvars(
            '$FERMI/diffuse/convolved/'+fname.replace('.fits','_8years_{}{}.fits'.format(suffix,fb) ))
        cd = CountDensity(iem_file, config.irfs, event_type=i, energies=energies, newfilename=newfile)

if __name__=='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description=""" Create count density files, in FITS HEALPix format, for the specified diffuse model
            uses configuration in current folder to obtain the IRF and exposure info
    """)
    #parser.add_argument('stream', nargs='*', default=None, help='optional Stream number')
    args = parser.parse_args()
    main()
     