"""
Configuration information

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/config.py,v 1.5 2015/07/24 17:56:02 burnett Exp $

"""

import os
import numpy as np
import pylab as plt
import pandas as pd

from . import analysis_base, _html


class Configuration(analysis_base.AnalysisBase):
    """Configuration information
    """

    def setup(self,**kw):
        self.plotfolder='config'
        
    def make_config(self):
        """Summary log files from the processing of this model
        %(html)s
        """
        self.html =''
        for filename in ('config.txt', '../config.txt', 'dataset.txt', 'converge.txt', 'summary_log.txt'):
            if not os.path.exists(filename): continue
            self.html += '<h4>%s</h4>\n<pre>%s</pre>' % (filename, open(filename).read())
        return None

    def all_plots(self):
        self.runfigures([self.make_config,])
