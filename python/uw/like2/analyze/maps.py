"""
Sky maps of various types

$Header$

"""
import os, sys, glob, pyfits
import pandas as pd
import numpy as np
import pylab as plt
import Image
# workaround for PIL not being nice?
# symptom: accessinit-hash-collision-3-both-1-and-1
#import PIL.Image
#sys.modules['Image']=PIL.Image

from skymaps import SkyDir
from uw.utilities import colormaps

from . import analysis_base
from . analysis_base import FloatFormat, html_table
from ..pub import display_map

def make_thumbnail(outfile):
        im = Image.open(outfile)
        im.thumbnail((200,100))
        im.save(outfile.replace('.png', '_thumbnail.png'))

class Maps(analysis_base.AnalysisBase):
    """Displays of various HEALPix-based sky maps
    <br> (Needs some text explaining how they are generated, specific definitions.)
    """
    def setup(self, **kw):
        self.plotfolder='maps'
    
    def ait_plots(self,   **kwargs):
        """Images from HEALPix FITS files
        %(healpix_plots)s
        """
        show_kw_dict=dict( 
            kde=dict(nocolorbar=True, scale='log',vmin=4.5,vmax=7.5, cmap=colormaps.sls),
            ts = dict(nocolorbar=False, vmin=10, vmax=25),
            galactic = dict(nocolorbar=True, scale='log', cmap=colormaps.sls),
            counts = dict(scale='log'),
            )
        dpi = kwargs.pop('dpi', 120)
        infits = glob.glob('hptables*.fits')
        if  len(infits)==0:
            raise Exception('No hptables found')
        hp="<p>Click on any of the thumbnails below to see an expanded version. All are AIT projections."
        for filename in infits:
            try:
                t = pyfits.open(filename)[1].data
            except Exception, msg:
                raise Exception('Failed to open fits file "%s": %s'% (filename, msg))
            nside = int(np.sqrt(len(t)/12.))
            
            hp += """<p>File <a href="../../%(filename)s?download=true">"%(filename)s"</a>:
                nside=%(nside)s, %(cols)d columns."""% dict(filename=filename, nside=nside, cols=len(t.dtype.names))
            for table in t.dtype.names:
                print 'processing table %s, length=%d' % (table, len(table))
                outfile = self._check_exist(table+'_ait.png')
                if outfile is not None: 
                    dm = display_map.DisplayMap(t.field(table))
                    show_kw=show_kw_dict.get(table, {})
                    dm.fill_ait(show_kw=show_kw, **kwargs)
                    #plt.title( '%s for %s' % (field, self.outdir))
                    plt.savefig(outfile, bbox_inches='tight', dpi=dpi)
                    make_thumbnail(outfile)
                    print 'wrote %s image and thumbnail' % outfile
                else: 
                    print 'using existing %s and thumbnail' % table
                hp+="""\n <p><b>%(table)s</b> <br><a href="%(path)s_ait.png"> 
                        <img alt="%(path)s_ait.png"  
                        src="%(path)s_ait_thumbnail.png" /></a> <br/>""" % dict(table=table, path=table)
            
        self.healpix_plots = hp + '\n'
        
    def _check_exist(self, filename, overwrite=False):
        """ return full filename, remove if exists"""
        fn = os.path.join(self.plotfolder, filename)
        if os.path.exists(fn):
            if overwrite: os.remove(fn)
            else: return None
        return fn


        
    def all_plots(self):
        self.runfigures([self.ait_plots,])
        pass