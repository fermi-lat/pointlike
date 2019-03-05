"""

Module providing handling of the LAT effective area.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/irfs/effective_area.py,v 1.1 2016/06/22 17:02:51 wallacee Exp $
Author: Eric Wallace

"""

__version__ = "$Revision: 1.1 $"

import os

import numpy as np
from astropy.io import fits
from scipy import interpolate

import skymaps

#To get something working for now, just use C++ implementation

class EffectiveArea(skymaps.EffectiveArea):
    """Placeholder interface to skymaps::EffectiveArea"""
    def __init__(self,filename,
                 aeff_extension="EFFECTIVE_AREA",
                 eff_params_extension="EFFICIENCY_PARAMETERS"):
        self.eff_params_extension = eff_params_extension
        skymaps.EffectiveArea.__init__(self,'',str(filename),str(aeff_extension))


"""
#TODO: 
#      Implement interpolation
#      Implement livetime correction

class EffectiveArea(object):

    def __init__(self,filename, 
                 aeff_extension="EFFECTIVE_AREA",
                 aeff_params_extension = "EFFICENCY_PARAMETERS",
                 thetacut = 66.4):
        hdus = fits.open(filename)
        aeff = hdus[aeff_extension]
        self.aeff = aeff.data['EFFAREA'][0]
        self.efficiency_pars = hdus[eff_params_extension].data[0]
        self.event_type = aeff.header['DETNAM'].lower()
        self.ebins = np.vstack([aeff.data['ENERG_LO'][0],aeff.data['ENERG_HI'][0]]).T
        self.cthetabins = np.vstack([aeff.data['CTHETA_LO'][0],aeff.data['CTHETA_HI'][0]]).T
        #mask = self.cthetabins[:,0]>=cthetamin
        #self.cthetabins = self.cthetabins[mask]
        #self.aeff = aeff['EFFAREA'][mask]
        #NB: Slightly different interpolation than used in skymaps.EffectiveArea
        self._interp = interpolate.RectBivariateSpline(self.cthetabins.mean(axis=1),
                                                       np.log10(self.ebins).mean(axis=1),
                                                       self.aeff)

    

    def value(self,energy=1000.,ctheta=1.0,interpolate=True):
        if interpolate:
            return self.interpolate(energy,ctheta)*1e4
        else:
            mask = (np.fmin(np.searchsorted(self.cthetabins[:,1],ctheta),self.cthetabins.shape[0]-1),
                np.fmin(np.searchsorted(self.ebins[:,1],energy),self.ebins.shape[0]-1))
            return self.aeff[mask]*1e4
    
    def interpolate(self,energy,ctheta):
        return self._interp(ctheta,np.log10(energy),grid=False)

    def __call__(self,energy=1000.,ctheta=1.0):
        return self.value(energy,ctheta)

    def get_livetime_factors(self,energy):
        pass

"""


        

