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
import pyfits

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
        print 'Running "to_xml"...'
        to_xml.main()
        
        print '\nRunning "to_fits"...'
        t = to_fits.MakeCat(self.sourcelist)
        self.fits_file = '_'.join(os.path.abspath('.').split('/')[-2:])+'.fits'
        t(self.fits_file )
        self.logstream=self.stoplog()
         
    def files(self):
        """Links to output files
        <ul>
         <li>FITS <a href="../../%(fits_file)s">%(fits_file)s</a></li>
         <li>XML  <a href="../../%(xml)s">%(xml)s</ax></li>
        </ul>
        
        """
        self.xml = glob.glob('*.xml')[0]
        
    def fits_summary(self):
        """FITS file summary
        Read back the FITS file, display column information
        %(fits_summary_table)s
        """
        t = pyfits.open(self.fits_file)[1].data
        df = pd.DataFrame(t)
        summary = dp.html_table(df.describe().T, float_format=dp.FloatFormat(3), href=False)
        self.fits_summary_table = summary.replace('%', '%%')
        
    def all_plots(self):
        self.runfigures([self.analysis,self.fits_summary, self.files,])
