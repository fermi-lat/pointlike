"""A module providing utilities for combining data with IRFs

Author: E. Wallace, M. Kerr
"""

__version__="$Revision: 1.3 $"
#$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pointspec2.py,v 1.3 2017/08/23 23:40:24 burnett Exp $

import os
import warnings
from math import cos,radians

import pointlike
import skymaps
from uw.data import dataman
from uw.like import pycaldb,pypsf
from uw.utilities import keyword_options
from uw.like.pointspec import SpectralAnalysis as OldSpectralAnalysis


class ExposureManager(object):
    """A class to handle the combination of effective area and livetime."""

    defaults = (('verbose',False,'If true, print a little more information'),)

    @keyword_options.decorate(defaults)
    def __init__(self,dm,caldb,**kwargs):
        """Initialize an ExposureManager
        
        Arguments:
            dm: a DataManager object containing the LivetimeCube(s) to be used.
            caldb: a CALDBManager specifying the CALDB information
        """
        keyword_options.process(self,kwargs)
        self.data_manager = dm
        self.caldb = caldb
        skymaps.EffectiveArea.set_CALDB(caldb.CALDB)
        skymaps.Exposure.set_cutoff(cos(radians(
                                    dm.dataspec.theta_cut['upper'])))
        inst = ['front', 'back']
        aeff_files = caldb.get_aeff()
        ok = [os.path.exists(file) for file in aeff_files]
        if not all(ok):
            raise Exception('one of CALDB aeff files not found: %s' %aeff_files)
        ### THB: Adjust this very late in game, all this pretty obsolete, for pass 8 CALDB format
        self.ea  = [skymaps.EffectiveArea('', filename, 'EFFECTIVE AREA_'+fb) for filename, fb in zip(aeff_files,['FRONT','BACK'])]
        if self.verbose: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        if dm.dataspec.use_weighted_livetime:
            self.exposure = [skymaps.Exposure(dm.lt,dm.weighted_lt,ea) for ea in self.ea]
        else:
            self.exposure = [skymaps.Exposure(dm.lt,ea) for ea in self.ea]

    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)

class SpectralAnalysis(OldSpectralAnalysis):
    """ Minimal inheritance from SpectralAnalysis to use new data interface."""

    # revised keywords to remove data specification from SpectralAnalysis interface
    defaults = (
        'keywords controlling data binning and livetime calculation',
        ('roi_dir',None,'aperture center; if None, assume all-sky analysis'),
        ('use_weighted_livetime',False,'Use the weighted livetime'),
        'keywords controlling instrument response',
        ('irf',None,'Which IRF to use'),
        ('psf_irf',None,'specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!'),
        ('CALDB',None,'override the CALDB specified by $CALDB.'),
        ('custom_irf_dir',None,'override the CUSTOM_IRF_DIR specified by the env. variable'),
        ('keywords controlling spectral analysis'),
        ('maxROI',10,'maximum ROI for analysis; note ROI aperture is energy-dependent = min(maxROI,r95(e,conversion_type))'),
        ('minROI',5,'minimum ROI for analysis; note ROI aperture is energy-dependent = max(minROI,r95(e,conversion_type))'),
        'miscelleneous keywords',
        ('quiet',False,'Set True to suppress (some) output'),
        ('verbose',False,'More output'),
    )

    def __init__(self, data_specification, **kwargs):
        """
        Create a new spectral analysis object.

        data_specification: an instance of DataSpec with links to the 
                            FT1/FT2, and/or binned data / livetime cube needed
                            for analysis (see docstring for that class). 
                            
                            Instance may be given as a pickle file.
        """
        if not isinstance(data_specification,dataman.DataSpec):
            if os.path.exists(data_specification):
                from cPickle import load
                try:
                    data_specification = load(file(data_specification))
                except UnpicklingError:
                    print 'Invalid pickle file for DataSpecification.'
                    raise Exception
        self.ae = self.dataspec = data_specification
        keyword_options.process(self, kwargs)

        #pixeldata name for backward compatibility
        self.dataman=self.pixeldata = self.dataspec()
        self.CALDBManager = pycaldb.CALDBManager(irf=self.irf,psf_irf=self.psf_irf,
            CALDB=self.CALDB,custom_irf_dir=self.custom_irf_dir)

        self.exposure  = ExposureManager(self.dataman,self.CALDBManager,verbose=self.verbose)
        self.psf = pypsf.CALDBPsf(self.CALDBManager)

    @property
    def emin(self):
        """Alias for backward compatibility"""
        warnings.warn(DeprecationWarning('EMIN SHOULD BE A PROPERTY OF DATASPEC'))
        return 1e1
    @property
    def emax(self):
        """Alias for backward compatibility"""
        warnings.warn(DeprecationWarning('EMIN SHOULD BE A PROPERTY OF DATASPEC'))
        return 1e6
    @property
    def conv_type(self):
        """Alias for backward compatibility"""
        warnings.warn(DeprecationWarning('CONV_TYPE SHOULD BE SPECIFICED AT FIT TIME'))
        return -1

    def set_psf_weights(self,skydir):
        """ Set the PSF to a new position.  Weights by livetime."""
        self.psf.set_weights(self.dataspec.ltcube,skydir)
