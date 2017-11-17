"""
Load a skymodel fror a dataframe (under development)
$Header$
"""
import os, glob
import cPickle as pickle
import numpy as np
import skymaps
from pointlike import IntVector

from . import (roimodel, sources, diffuse, extended)

class ROImodelFromDict(roimodel.ROImodel):
    """Load from a dictionary
    """

    def load_sources(self, roi_index, source_radius=12):
        fn = glob.glob('skymodel*.pkl')[-1]
        self.skymodel = pickle.load(open('skymodel.pkl'))

    def __repr__(self):
        return '%s.%s : %d global, %d local, %d total sources for ROI %d' \
            % (self.__module__, self.__class__.__name__, self.global_count,self.local_count,  
            len(self), self.index)


