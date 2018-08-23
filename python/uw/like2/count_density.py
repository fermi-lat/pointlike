"""
Create observed count density files

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/count_density.py,v 1.1 2017/11/17 22:50:36 burnett Exp $

author:  Toby Burnett
"""
import numpy as np
import pandas as pd
from astropy.io import fits
import os, healpy, argparse

from uw.like2.pub import healpix_map #make this local in future
from uw.like2 import (diffuse, convolution, configuration )

class Simpson(object):
    def __init__(self, cube):
        """cube: functor, of energy, returning HEALPix map of sky"""
        self.cube = cube
    def __call__(self, energy):
        # filer the function for Nan values, replace with zeros
        r=self.cube(energy)
        bad = np.isnan(r)
        if np.any(bad):
            print 'Warning: {} Nan pixels at {:.0f} MeV'.format(np.sum(bad),energy)
            r[bad]=0.
        return r
    def average(self, a, b, n=1):
        """return logarithmic average of f over energies from a to b, using Simpson rule index n"""
        if n==1: return self(np.sqrt(a*b))
        x = np.logspace(np.log10(a),np.log10(b),n+1)
        w = np.ones(n+1)/float(n); w[0]=w[-1]=0.5/float(n)
        r = w[0]*self(x[0])
        for i in range(1,n+1):
            r += w[i] * self(x[i])
        return r

class CountDensity(object):
    """Manage computation of expected count density from a diffuse map.
    Includes exposure and convolution with PSF
    (But not, for now, energy dispersion)
    """    
    def __init__(self, filename, irfman, event_type, 
            energies=None, 
            ebins=None, simpson_index=[], 
            psf_energy_factor=1.0,
            newfilename=None, overwrite=True, adjuster=None):
        """
        filename : string
            name of a file in HEALPix format
        irfman : irf/IrfMan object for exposure, PSF
        event_type : int
            0 or 1 for Front, Back
        energies : list of float | None
            if None, use the layers in the diffuse 
        ebins : list of float | None
            if specified, assume to be energy bin edges (one more than bins)
            energies will then be the geometric means
        simpson : list of int | None
            Simpson's rule index to use for energy integration if ebins specified
            Use 1 beyond range of list
            Note that "1" means evaluate at logE bin center
        psf_energy_factor : float, defalut 1.0
            Adjust the energy at which the PSF is evaluated
        adjuster : None | function object
            Multiply the resulting map by adjuster(energy), perhaps for dispersion correction            
        """
        print 'loading filename {}\nConvolving with event type {}, psf factor {}'.format(filename, event_type, psf_energy_factor)
        print 'Will write to {}'.format(newfilename) if newfilename is not None else 'No output file specified'
        
        self.hpcube= diffuse.HealpixCube(filename)
        self.hpcube.load()
        self.nside = self.hpcube.nside
        self.irfman = irfman
        self.event_type= event_type
        if ebins is not None:
            #set up for Simpson integration
            self.simpson = Simpson(self)
            self.simpson_index = simpson_index
            self.emin, self.emax = ebins[:-1], ebins[1:]
            self.energies = np.sqrt(self.emin * self.emax)
        else:
            self.simpson=None
            self.energies = self.hpcube.energies if energies is None else energies
        self.last_energy=0; self.last_map=None
        print 'Will process {} energies from {:.0f} to {:.0f}'.format(
                    len(self.energies), self.energies[0],self.energies[-1])
        self.psf_energy_factor=psf_energy_factor
        self.adjuster = adjuster
        if self.adjuster is not None:
            print 'IEM adjuster specified, {}'.format(self.adjuster.__class__)

        if newfilename is not None:
            self.replace_skymaps()
            self.replace_energies()
            self.writeto(newfilename, overwrite=overwrite)
 
    def __call__(self, energy):
        print 'Convolving for energy {:.0f}'.format(energy)
        if energy==self.last_energy: return self.last_map
        self.last_energy=energy
        self.hpcube.setEnergy(energy)
        iem_map = self.hpcube.column(energy) #note will interpolate
        exp_map = exposure_map(self.irfman, 
            self.event_type, energy, self.nside )
        psf = self.irfman.psf(self.event_type, energy* self.psf_energy_factor)
        count_density_map=convolve_healpix(iem_map * exp_map, psf, 5*psf.r68, lmax_limit=500)
        if self.adjuster is not None and energy<10000.:
            count_density_map *= self.adjuster(energy, self.event_type)
        self.last_map=count_density_map
        return count_density_map

    def convert_layer(self, bin_index):
        if self.simpson is None or bin_index>len(self.simpson_index)-1:
            return self(self.energies[bin_index])
        n = self.simpson_index[bin_index]
        print 'Simpson with n={}'.format(n)
        return self.simpson.average(self.emin[bin_index], self.emax[bin_index], n)
        # cube = self.hpcube
        # energy = self.energies[index]
        # print 'Map for energy {:.0f}'.format(energy)
        # cube.setEnergy(energy)
        # iem_map = cube.column(energy) #note will interpolate

        # exp_map = exposure_map(self.irfman, 
        #     self.event_type, energy, self.nside )
        
        # psf = self.irfman.psf(self.event_type, energy* self.psf_energy_factor)
        
        # count_density_map=convolve_healpix(iem_map * exp_map, psf, 5*psf.r68, lmax_limit=500)
        # if self.adjuster is not None and energy<10000.:
        #     count_density_map *= self.adjuster(energy, self.event_type)
        # return count_density_map

    def replace_skymaps(self):
        tt = np.array([self.convert_layer(i) for i in range(len(self.energies))])
        newcols=[fits.Column(name='Bin{:d}'.format(i), format='E', array=tt[i]) for i,t in enumerate(tt)]
        skymap_hdu = fits.BinTableHDU.from_columns(newcols)

        skymap_hdu.header.update(
            LASTPIX = 12*self.nside**2-1,                                                  
            EXTNAME = 'SKYMAP  ',                                                            
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
        
def convolve_healpix(input_map,  func, thetamax, quiet=True, lmax_limit=None ):
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
    if lmax_limit is not None and lmax>lmax_limit:
        lmax=lmax_limit
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


def main(factor=1.0,energies=None, adjuster=None, simpson_index=[4,4,2,2,1], ebins=None):
    """

    """
    config= configuration.Configuration('.', quiet=True, postpone=True)
    iem_file = config.diffuse['ring']['iemfile'] #original file
    filename = config.diffuse['ring']['filename'] #will be the convolved files
    print 'IEM file: {}'.format(iem_file)
    print 'Livetime cube: {}'.format(config.irfs.ltcube)

    if ebins is None and energies is None:
        # get energies from data
        bf = fits.open(config.dataset.binfile)
        bands = bf['BANDS']
        emin,emax = bands.data['E_MIN'], bands.data['E_MAX']
        energies = sorted(list(set(np.sqrt(emin*emax) * 1e-3) ))#from keV
    
    factor = config.diffuse['ring'].get('psf_energy_factor', 1.0)

    for i,fb in enumerate(config.event_type_names):
        
        newfile = os.path.expandvars('$FERMI/diffuse/'+filename.replace('*',fb) )
        cd = CountDensity(iem_file, config.irfs, event_type=i, 
            psf_energy_factor=factor, adjuster=adjuster,
            energies=energies, newfilename=newfile, ebins=ebins, simpson_index=simpson_index)

if __name__=='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description=""" Create count density files, in FITS HEALPix format, for the specified diffuse model
            uses configuration in current folder to obtain the IRF and exposure info
    """)
    #parser.add_argument('stream', nargs='*', default=None, help='optional Stream number')
    args = parser.parse_args()
    main()
     