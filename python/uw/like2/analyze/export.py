"""
Export processing
$Header$

"""
import os

from uw.like2.pipeline import diagnostic_plots as dp
from .. import to_xml
from .. import to_fits
import numpy as np
import pandas as pd
import pylab as plt

class Export(dp.SourceInfo):
    """Manage, and document an export step
    """
    def setup(self, **kw):
        super(Export, self).setup(**kw)
        self.plotdir = 'export'
        print 'setup Exporting...'
        
        pass
        
    def all_plots(self):
        pass
