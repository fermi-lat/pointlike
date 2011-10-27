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

    def restore(self, roi=None):

        if roi is None: roi = self.roi

        roi_dir=roi.roi_dir
        quiet=roi.quiet

        psm = ROIPointSourceManager(self._copy_ps(self.point_sources), roi_dir, quiet=quiet)
        dsm = ROIDiffuseManager(self._copy_bg(self.bgmodels), roi_dir, quiet=quiet)

        sa=roi.sa

        from uw.like.roi_analysis import ROIAnalysis
        kwargs=defaults_to_kwargs(roi,ROIAnalysis)

        roi.__init__(roi_dir,psm,dsm,sa, **kwargs)
