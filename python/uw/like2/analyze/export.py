"""
Export processing
$Header#

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
    <p>Expect that <a href="../peak_finder/index.html?skipDecoration">findpeak</a> has been run to update the sources file
    <br>Source files found: %(sourcecsv)s
    <p>Generates XML and FITS files from csv source list.
    """
    def setup(self, **kw):
        super(Export, self).setup(**kw)
        self.plotfolder = 'export'
        self.sourcecsv = sorted(glob.glob('source*.csv'))
        self.sourcelist=pd.read_csv(self.sourcecsv[-1], index_col=0)
    
    def analysis(self):
        """Log of analysis stream
        <pre>%(logstream)s</pre>"""
        self.startlog()
        print 'Running to_xml...'
        to_xml.main()
        
        print 'Running to_fits...'
        t = to_fits.MakeCat(self.sourcelist)
        outfile = '_'.join(os.path.abspath('.').split('/')[-2:])+'.fits'
        t(outfile)
        self.logstream=self.stoplog()
         
    def files(self):
        """Links to output files
        <ul>
         <li>FITS <a href="../../%(fits)s">%(fits)s</a></li>
         <li>XML  <a href="../../%(xml)s">%(xml)s</ax></li>
        </ul>
        
        """
        self.fits=glob.glob('*fits')[0]
        self.xml = glob.glob('*.xml')[0]
        
    def all_plots(self):
        self.runfigures([self.analysis,self.files,])
