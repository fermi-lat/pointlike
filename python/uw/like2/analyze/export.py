"""
Export processing
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/export.py,v 1.5 2013/06/18 22:08:19 burnett Exp $

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
    <p>Generates XML and FITS files from the file %(sourcecsv)s.
    <p>Expect that <a href="../peak_finder/index.html?skipDecoration">findpeak</a> has been run to update the sources file. In this case, there are a set of sources for which the standard localization analyis failed, and have had a localization moments analysis. For these sources, bit 4 of flags is set to indicate that the error ellipse parameters were replaced. The LocaliationQuality is not changed, however. 
    The error ellipse for any of these sources is thus only an approximate representation of the uncertainty, and the user is urged to examine the corresponding TS maps: png and FITS versions can be found in <a href="../../tsmap_fail/">this folder</a>.
    
    
    <p>A final selection on the error ellipse size, and systematic corrections are made here:
    <table border="1">
      <tr><td class="index">Cut on initial semi-major axis (deg)</td><td>%(error_box_cut)s</td> </tr>
      <tr><td class="index">Multiplicative systematic factor</td><td>%(error_box_factor)s</td> </tr>
      <tr><td class="index">Fixed systematic, added in quadrature (deg)</td><td>%(error_box_add)g</td> </tr>
    </table>
    The multiplicative factor is examined <a href="../associations/index.html?skipDecoration">here.</a>
    """
    def setup(self, **kw):
        super(Export, self).setup(**kw)
        self.plotfolder = 'export'
        self.sourcecsv = sorted(glob.glob('source*.csv'))[-1]
        self.sourcelist=pd.read_csv(self.sourcecsv, index_col=0)
        self.error_box_factor = 1.05
        self.error_box_add = 5e-3
        self.error_box_cut = 0.25
        self.cuts = '(sources.ts>10)*(sources.a<%.2f)+pd.isnull(sources.locqual)' %self.error_box_cut
    
    def analysis(self):
        """Analysis log
        <pre>%(logstream)s</pre>"""
        self.startlog()
        print 'Running "to_xml"...'
        to_xml.main(cuts=self.cuts)
        
        print '\nRunning "to_fits"...'
        self.fits_file = '_'.join(os.path.abspath('.').split('/')[-2:])+'.fits'
        to_fits.main(self.fits_file,  cuts=self.cuts,
                     localization_systematic = (self.error_box_factor, self.error_box_add))

        self.logstream=self.stoplog()
         
    def files(self):
        """Links to output files
        <ul>
         <li>FITS <a href="../../%(fits_file)s">%(fits_file)s</a></li>
         <li>XML  <a href="../../%(xml)s">%(xml)s</a></li>
        </ul>
        
        """
        self.xml = glob.glob('*.xml')[0]
        
    def fits_summary(self):
        """FITS file summary
        Read back the FITS file, display numerical column information.
        %(fits_summary_table)s
        * flux13 and unc_flux13 are not in the FITS file, but set to [Unc_]Flux_Density * 1e13 for numerical display.
        """
        t = pyfits.open(self.fits_file)[1].data
        df = pd.DataFrame(t)
        df['flux13*'] = df['Flux_Density']*1e13
        df['unc_flux13*'] = df['Unc_Flux_Density']*1e13
        summary = dp.html_table(df.describe().T, float_format=dp.FloatFormat(3), href=False)
        self.fits_summary_table = summary.replace('%', '%%')
        print 'Check: %s' % df
        
    def all_plots(self):
        self.runfigures([self.analysis,self.fits_summary, self.files,])
