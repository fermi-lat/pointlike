"""
Export processing
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/export.py,v 1.16 2016/10/28 20:48:14 burnett Exp $

"""
import os, glob

from . import sourceinfo 
from .. import to_xml
from .. import to_fits
from . analysis_base import FloatFormat, html_table
import numpy as np
import pandas as pd
import pylab as plt
from astropy.io import fits as pyfits

class Export(sourceinfo.SourceInfo):
    """Export to XML and FITS
    <p>Generates XML and FITS files from the file %(sourcecsv)s.
    <p>Expect that <a href="../peak_finder/index.html?skipDecoration">findpeak</a> has been run to update 
    the sources file. In this case, there are a set of sources for which the standard localization analyis failed, 
    and have had a localization moments analysis. For these sources, bit 4 of flags is set to indicate that the 
    error ellipse parameters were replaced. The LocaliationQuality is not changed, however. 
    The error ellipse for any of these sources is thus only an approximate representation of the uncertainty, 
    and the user is urged to examine the corresponding TS maps: png and FITS versions can be found in 
    <a href="../../tsmap_fail/">this folder</a>.
    
    <p>A final selection on the error ellipse size, and systematic corrections are made here:
    <table border="1">
      <tr><td class="index">Cut on initial semi-major axis (deg)</td><td>%(error_box_cut)s</td> </tr>
      <tr><td class="index">Multiplicative systematic factor</td><td>%(error_box_factor)s</td> </tr>
      <tr><td class="index">Fixed systematic, added in quadrature to r95 (deg)</td><td>%(error_box_add)g</td> </tr>
    </table>
    The multiplicative factor is examined <a href="../associations/index.html?skipDecoration">here.</a>
    """
    def setup(self, **kw):
        super(Export, self).setup(**kw)
        self.plotfolder = 'export'
        self.sourcecsv = sorted(glob.glob('source*.csv'))[-1]
        self.sourcelist=pd.read_csv(self.sourcecsv, index_col=0)
        systematic = self.config['localization_systematics']\
             if 'localization_systematics' in self.config.keys() else (1.1, 0.3)
        self.error_box_factor = systematic[0]
        self.error_box_add = systematic[1]/60.
        self.error_box_cut = 0.5
        self.cuts = '(sources.ts>10) & (sources.a<%.2f) | pd.isnull(sources.locqual)' %self.error_box_cut

        name_root = '_'.join(os.path.abspath('.').split('/')[-2:])
        versions = glob.glob(name_root+'*.fits')
        if len(versions)==0:
            self.fits_file = name_root+'.fits'
        else:
            last_version = os.path.splitext(versions[-1])[0].split('_')[-1];
            if last_version==self.skymodel:
                next_version='1'
            else:
                next_version = str(int(last_version)+1)
            self.fits_file = name_root+'_'+next_version+'.fits'
        print ('Will write to {}'.format(self.fits_file))

    def analysis(self, fits_only=True): # for now
        """Analysis log
        <pre>%(logstream)s</pre>"""
        self.startlog()

        
        print ('\nRunning "to_fits"...')
        #self.fits_file = '_'.join(os.path.abspath('.').split('/')[-2:])+'.fits'
        to_fits.main(self.fits_file,  cuts=self.cuts,
                     localization_systematic = (self.error_box_factor, self.error_box_add)
                     )
        self.xml_file=''
        if not fits_only:
            print ('Running "to_xml"...')
            self.xml_file = self.fits_file.replace('.fits', '.xml' )
            to_xml.main(filename=[self.xml_file], cuts=self.cuts)

        self.logstream=self.stoplog()
         
    def files(self):
        """Links to output file(s)
        <ul>
         <li>FITS <a href="../../%(fits_file)s?download=true">%(fits_file)s</a></li>
         %(xml_link)s
        </ul>
        
        """
        if os.path.exists(self.xml_file):
            self.xml_link='<li>XML  <a href="../../%(xml_file)s?download=true">%(xml_file)s</a></li>'.format(self.xml_file)
        else:
            self.xml_link=''
        # <a href="../../%(xml)s?download=true">%(xml)s</a></li>
        #self.xml = glob.glob('*.xml')[0]
        
    def fits_summary(self):
        """FITS file summary
        Read back the FITS file, display numerical column information.
        %(fits_summary_table)s
        * flux13 and unc_flux13 are not in the FITS file, but set to [Unc_]Flux_Density * 1e13 for numerical display.
        """
        t = pyfits.open(self.fits_file)[1].data
        # remove columns that have multiple dimensions
        for j in range(3):#??? why
            for i,col in enumerate(t.columns):
                if len(col.array.shape)>1:
                    t.columns.del_col(i)

        tt=pyfits.BinTableHDU.from_columns(t.columns)
        df = pd.DataFrame(tt.data)
        df['flux13*'] = df['Flux_Density']*1e13
        df['unc_flux13*'] = df['Unc_Flux_Density']*1e13
        summary = html_table(df.describe().T, float_format=FloatFormat(3),
                heading='', href=False, maxlines=50)
        self.fits_summary_table = summary.replace('%', '%%')
        # creates error??
        #print ('Check: %s' % df)
        
    def all_plots(self):
        self.runfigures([self.analysis,self.fits_summary, self.files,])
