"""Module for managing instrument response functions.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/irfs/irfman.py,v 1.2 2016/06/27 23:06:35 wallacee Exp $
Author: Eric Wallace
"""
__version__="$Revision: 1.2 $"

import numpy as np

from . import (caldb, psf, effective_area, exposure, IrfError)


class IrfManager(object):
    """An object to manage loading sets of IRFs."""
    event_type_names = ('front','back', 'psf0','psf1','psf2','psf3','edisp0',
                         'edisp1','edisp2','edisp3')
    event_type_partitions = dict(fb = (0,1),
                                 psf = (2,3,4,5),
                                 edisp = (6,7,8,9))

    def __init__(self, dataset, irf_dir="$CALDB"):
        self.caldb = caldb.CALDB(irf_dir)
        irfname = dataset.irf
        irfname_parts = irfname.split('_')
        self.event_class = irfname_parts.pop(1)
        self.irf_version = '_'.join(irfname_parts)
        if dataset.psf_event_types:
            ets = 'psf'
        else:
            ets = 'fb'
        self.event_types = self._parse_event_type(ets)
        if dataset.legacy:
            #No DSS keywords, assume dataset has thetacut member
            #TODO: Fix dataman.DataSpec to provide a sensible default
            #      or find something to reference that works for both cases
            self.cthetamin = np.cos(np.radians(dataset.thetacut))
        else:
            self.cthetamin = np.cos(np.radians(dataset.theta_cut.get_bounds()[1]))
        self.dataset = dataset
        self._load_irfs()

    def _load_irfs(self):
        psf_info = self.caldb('psf',version=self.irf_version,
                                    event_class = self.event_class,
                                    event_type = self.event_types)
        aeff_info = self.caldb('aeff',version=self.irf_version,
                                      event_class = self.event_class,
                                      event_type = self.event_types)
        self._aeff = {et:effective_area.EffectiveArea(d['filename'],
                                     aeff_extension=d['extensions']['EFF_AREA'],
                                     eff_params_extension=d['extensions']['EFFICIENCY_PARS'])
                            for et,d in aeff_info.items()}
        self._exposure = {et:exposure.Exposure(self.dataset.lt,aeff,cthetamin=self.cthetamin)
                                for et,aeff in self._aeff.items()}
        self._psf = {et:psf.PSF(d['filename'],
                                     rpsf_extension=d['extensions']['RPSF'],
                                     psf_scaling_extension=d['extensions']['PSF_SCALING'],
                                     exposure = self._exposure[et])
                            for et,d in psf_info.items()}

        
    def psf(self,event_type,energy):
        """Return a BandPSF for the given energy and event_type."""
        et = self._parse_event_type(event_type)
        return self._psf[et].band_psf(energy)

    def exposure(self,event_type,energy):
        """Return a BandExposure obect for the given energy and event_type."""
        et = self._parse_event_type(event_type)
        return self._exposure[et].band_exposure(energy)

    def _parse_event_type(self,event_type):
        """Find the event type or types for a given event_type selection""" 
        if hasattr(event_type,'__iter__'):
            return [self._parse_event_type(et) for et in event_type]

        try:
            event_type = event_type.lower()
        except AttributeError: #Not a string, assume index
            if event_type in range(10):
                return event_type
            else:
                raise IrfError("Invalid event type index {}. Should be 0-9.".format(event_type))

        try:
            return self.event_type_names.index(event_type)
        except ValueError as exc:
            try: return self.event_type_partitions[event_type]
            except KeyError:
                raise exc
