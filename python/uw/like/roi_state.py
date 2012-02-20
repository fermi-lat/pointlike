import copy

import numpy as np
from uw.like.roi_managers import ROIPointSourceManager, ROIDiffuseManager
from uw.like.pointspec_helpers import get_default_diffuse_mapper
from uw.utilities.keyword_options import defaults_to_kwargs

class PointlikeState(object):
    """Save the parameter state of a pyLikelihood object and provide a
        method to restore everything."""

    def _copy_ps(self,ps):
        return [ i.copy() for i in ps ]

    def _copy_bg(self,bg):
        roi = self.roi
        mapper=get_default_diffuse_mapper(roi.sa, roi.roi_dir, roi.quiet)
        return [ mapper(i.diffuse_source.copy()) for i in bg ]

    def __init__(self, roi):
        self.roi = roi
        self.point_sources = self._copy_ps(roi.psm.point_sources)
        self.bgmodels = self._copy_bg(roi.dsm.bgmodels)

        from uw.like.roi_analysis import ROIAnalysis
        self.roi_kwargs=defaults_to_kwargs(roi,ROIAnalysis)

    def _restore_spectra(self, roi=None):
        """ Restore the spectral of all sources in the region. """

        if roi is None: roi = self.roi

        map = dict()
        map = {ps.name:ps.model for ps in self.point_sources}
        map.update({bgm.name:bgm.smodel for bgm in self.bgmodels})

        for name in np.append(roi.psm.names, roi.dsm.names):
            roi.modify(which=name, model=map[name].copy(), keep_old_flux=False)


    def _restore_everything(self, roi=None, **kwargs):
        """ Restore the spectral + spatial information of all sources in the region. """

        if roi is None: roi = self.roi

        roi_dir=roi.roi_dir
        quiet=roi.quiet

        psm = ROIPointSourceManager(self._copy_ps(self.point_sources), roi_dir, quiet=quiet)
        dsm = ROIDiffuseManager(self._copy_bg(self.bgmodels), roi_dir, quiet=quiet)

        sa=roi.sa

        k = copy.deepcopy(self.roi_kwargs)
        k.update(kwargs)

        roi.__init__(roi_dir,psm,dsm,sa, **k)

    def restore(self, roi=None, just_spectra=False, **kwargs):
        if just_spectra:
            self._restore_spectra(roi, **kwargs)
        else:
            self._restore_everything(roi, **kwargs)


