import numpy as np
from uw.like.roi_managers import ROIPointSourceManager, ROIDiffuseManager
from uw.like.pointspec_helpers import get_default_diffuse_mapper

class PointlikeState(object):
    """Save the parameter state of a pyLikelihood object and provide a
        method to restore everything."""
    def __init__(self, roi):
        self.roi = roi
        self.point_sources = [ i.copy() for i in roi.psm.point_sources ]

        mapper=get_default_diffuse_mapper(roi.sa, roi.roi_dir, roi.quiet)
        self.bgmodels = [ mapper(i.copy()) for i in roi.dsm.diffuse_sources ]

    def restore(self):

        roi = self.roi
        roi_dir=roi.roi_dir
        quiet=roi.quiet

        psm = ROIPointSourceManager(self.point_sources, roi_dir, quiet=quiet)
        dsm = ROIDiffuseManager(self.bgmodels, roi_dir, quiet=quiet)

        sa=roi.sa
        roi.__init__(roi_dir,psm,dsm,sa)
