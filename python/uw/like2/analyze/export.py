"""
Export processing
$Header$

"""
import os, glob

from uw.like2.pipeline import diagnostic_plots as dp
from .. import to_xml
from .. import to_fits
import numpy as np
import pandas as pd
import pylab as plt

class Export(dp.SourceInfo):
    """Manage, and document an export step
    <p>Expect that findpeak has been run to update the sources file
    <br>Source files found: %(sourcecsv)s
    """
    def setup(self, **kw):
        super(Export, self).setup(**kw)
        self.plotfolder = 'export'
        self.sourcecsv = sorted(glob.glob('source*.csv'))
        self.sourcelist=pd.read_csv(self.sourcecsv[-1], index_col=0)
    
    def log(self):
        """Log of analysis stream
        <pre>%(logstream)s</pre>"""
        self.startlog()
        to_xml.main()
        t = to_fits.MakeCat(self.sourcelist)
        outfile = '_'.join(os.path.abspath('.').split('/')[-2:])+'.fits'
        t(outfile)
        self.logstream=self.stoplog()
         
    def all_plots(self):
        self.runfigures([self.log,])
