"""
Sky maps of various types

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/maps.py,v 1.3 2013/09/07 11:51:08 burnett Exp $

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
    
    def ait_plots(self, pattern='hptables*.fits',  **kwargs):
        show_kw_dict=dict( 
            kde=dict(nocolorbar=True, scale='log',vmin=4.5,vmax=7.5, cmap=colormaps.sls),
            ts = dict(nocolorbar=False, vmin=10, vmax=25),
            galactic = dict(nocolorbar=True, scale='log', cmap=colormaps.sls),
            counts = dict(scale='log'),
            )
        show_kw_default = dict(vmin=kwargs.pop('vmin',None), vmax=kwargs.pop('vmax',None))
        dpi = kwargs.pop('dpi', 120)
        infits = glob.glob(pattern)
        if  len(infits)==0:
            raise Exception('No files match pattern %s' %pattern)
        hp="<p>Click on any of the thumbnails below to see an expanded version. All are AIT projections."
        for filename in infits:
            try:
                t = pyfits.open(filename)[1].data
            except Exception, msg:
                raise Exception('Failed to open fits file "%s": %s'% (filename, msg))
            nside = int(np.sqrt(len(t)/12.))
            
            hrow = '\n<table><tr>'
            mrow = '</tr>\n<tr>'
            for table in t.dtype.names:
                print 'processing table %s, length=%d' % (table, len(table))
                outfile = self._check_exist(table+'_ait.png')
                if outfile is not None: 
                    dm = display_map.DisplayMap(t.field(table))
                    show_kw=show_kw_dict.get(table, show_kw_default)
                    dm.fill_ait(show_kw=show_kw, **kwargs)
                    plt.savefig(outfile, bbox_inches='tight', dpi=dpi)
                    make_thumbnail(outfile)
                    print 'wrote %s image and thumbnail' % outfile
                else: 
                    print 'using existing %s and thumbnail' % table
                hrow += '<td>%s</td>' % table 
                mrow += """\n<td><a href="%(path)s_ait.png"> 
                        <img alt="%(path)s_ait.png"  
                        src="%(path)s_ait_thumbnail.png" /></a> </td>""" % dict(path=table)
                #hp+="""\n <p><b>%(table)s</b> <br><a href="%(path)s_ait.png"> 
                #        <img alt="%(path)s_ait.png"  
                #        src="%(path)s_ait_thumbnail.png" /></a> <br/>""" % dict(table=table, path=table)
            
        hp += + hrow + mrow +'</tr>\n</table>'
        hp += """<p>These were projected from the nside=%(nside)s HEALPix FITS file <a href="../../%(filename)s?download=true">"%(filename)s"</a>:
           containing %(cols)d images.
           It can be examined with <a href="http://aladin.u-strasbg.fr">Aladin</a>"""% dict(filename=filename, nside=nside, cols=len(t.dtype.names))
        return hp
        
    def hptables(self):
        """Images from HEALPix FITS files
        %(healpix_plots)s
        """
        self.healpix_plots = self.ait_plots('hptables*.fits')
        
    def diffuse_corrections(self, vmin=0.9, vmax=1.1):
        """Images from the diffuse correction
        <p>The maps show the locations of the applied diffuse correction, for the first four bands, from 100 MeV to 1 GeV.
        %(diffuse_corr)s
        """
        self.diffuse_corr = self.ait_plots('diffuse_corr.fits', vmin=vmin, vmax=vmax)
    
    def _check_exist(self, filename, overwrite=False):
        """ return full filename, remove if exists"""
        fn = os.path.join(self.plotfolder, filename)
        if os.path.exists(fn):
            if overwrite: os.remove(fn)
            else: return None
        return fn


        
    def all_plots(self):
        self.runfigures([self.hptables, self.diffuse_corrections,])
        pass