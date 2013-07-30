"""
Create a Pivot collection

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/sourceinfo.py,v 1.7 2013/07/21 15:17:52 burnett Exp $

"""
import os, pickle
import numpy as np
import pylab as plt
import pandas as pd

from uw.utilities import makepivot
from . import sourceinfo
from . diagnostics import FloatFormat
from . _html import HTMLindex

class Collection(sourceinfo.SourceInfo):
    """Catalog presented as a Pivot collection
    All the sources with TS>10 are shown.
    """
    def setup(self, **kw):
        super(Collection, self).setup(**kw)
        self.plotfolder='collection'
        
    def make_collection(self):
        """Pivot table with all sources
        Inputs are a table of sources, a <a href="%(csv_file)s?upload=True">csv file</a>
        and a <a href="%(sedfig_zip)?upload=true">zip file</a> containing all SED plots.
        
        %(full_pivot_link)s
        """
        self.csv_file = 'sources_%s.csv' % self.skymodel
        self.sedfig_zip = 'sedfig.zip'
        assert os.path.exists(csv_file), 'Did not find the file %s'% self.csv_file
        try:
            pc =makepivot.MakeCollection('all sources %s' % self.skymodel, 'sedfig', self.csv_file)
            self.full_pivot_link = """\
                    <p>These can be examined with a 
                    <a href="%s/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
                    which requires Silverlight."""  % (makepivot.http_host, pc.cId)
        except Exception, msg: 
            print "**** Failed to make pivot table, perhaps need to run sedinfo first: %s" % msg

    def all_plots(self):
        self.runfigures([self.make_collection,])