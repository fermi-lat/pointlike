"""
Analyze the contents of tsmap HEALPix tables

$Header$

"""

import os, glob, time
from astropy.io import fits as pyfits
import numpy as np
import pylab as plt
import pandas as pd

from matplotlib.colors import LogNorm
from uw.like2.pub import healpix_map
from skymaps import Band, SkyDir
from . import analysis_base, sourceinfo
from ..pipeline import check_ts, ts_clusters
from analysis_base import html_table, FloatFormat

class TStables(sourceinfo.SourceInfo):
    """Process a set of TS tables for source finding
    <br>
    %(tsmap_analysis)s
    """
    def setup(self, **kw):
        super(TStables, self).setup(**kw)
        self.plotfolder='sourcefinding'

    def all_plots(self):
        self.runfigures([]);    