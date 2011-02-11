"""
Utilities for managing Healpix arrays
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pub/healpix_map.py,v 1.1 2011/01/24 22:03:44 burnett Exp $
"""
import os,glob,pickle
import pylab as plt
import numpy as np
import pyfits
from skymaps import Band, SkyDir, Hep3Vector, PySkyFunction
import skymaps
from uw.utilities import image
from uw.like import pycaldb
from .. import maps 

class HParray(object):
    """ base class, implement a HEALPix array, provide AIT plot
    """
    def __init__(self, name, vec): 
        self.name=name
        self.vec = vec
        self.nside = int(np.sqrt(len(vec)/12))
        assert len(self.vec)==12*self.nside**2, 'length of %s not consistent with HEALPix' % self.name
    def __getitem__(self, index):
        return self.vec[index]
    def __len__(self): return 12*self.nside**2
    def getcol(self, type=np.float32): return np.asarray(self.vec, type)
    def skyfun(self, skydir):
        return self[Band(self.nside).index(skydir)]
    def plot(self, title='', axes=None, fignum=30, ait_kw={}, **kwargs):
        """ make an AIT skyplot of a HEALpix array
        crec : array
            must be sorted according to the HEALpix index
        title : string
            set the figure title
        ait_kw : dict
            to set kwargs for image.AIT, perhaps pixelsize
        
        Other args passed to imshow
        """
        band = Band(self.nside)
        def skyplotfun(v):
            skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
            index = band.index(skydir)
            return self[index]
        if axes is None:
            plt.close(fignum)
            fig = plt.figure(fignum, figsize=(12,6))
        ait=image.AIT(PySkyFunction(skyplotfun) ,axes=axes, **ait_kw)
        ait.imshow(title=title, **kwargs)
        return ait
        
class HPtables(HParray):
    """ assemble an array from tables in a set of ROIs
    """
    def __init__(self, tname, outdir, nside=512, roi_nside=12):
        """ combine the tables generarated at each ROI
            nside must be consistent with the sub division for tables 
        """
        self.name = tname
        self.nside = nside
        files =glob.glob(os.path.join(outdir, '%s_table'%tname,'*.pickle'))
        nf = len(files)
        assert nf>0, 'no pickle files found in %s' % os.path.join(outdir, '%s_table'%tname)
        if nf<1728: print 'warning: missing %d files' % (1728-nf)
        files.sort()

        self.vec = np.zeros(12*nside**2)
        self.vec.fill(np.nan)
        pklist = [pickle.load(open(f)) for f in files]
        i12 = [int(f[-11:-7]) for f in files]
        index_table = maps.make_index_table(roi_nside, nside)
        for index, pk in zip(i12,pklist):
            indeces = index_table[index]
            for i,v in enumerate(pk):
                self.vec[indeces[i]]=v
    
class HPskyfun(HParray):
    """generate from a sky function
    """
    def __init__(self, name, skyfun, nside):
        """ skyfun : function of position
        """
        self.name = name
        self.skyfun = skyfun
        self.nside=nside
        self.dirfun = Band(self.nside).dir
    def __getitem__(self, index):
        return self.skyfun(self.dirfun(index))
        
    def getcol(self, type=np.float32):
        return np.asarray([self[index] for index in xrange(12*self.nside**2)],type)

class HPresample(HParray):
    """ resample from another HParray object with different nside
    """

    def __init__(self, other, nside=256):
        self.name=other.name
        self.vec = other.vec
        self.nside =nside
        self._nside = int(np.sqrt(len(self.vec)/12) ) # internal 
        self.dirfun = Band(self.nside).dir # external direction from index
        self._indexfun = Band(self._nside).index #internal index from direction

    def __getitem__(self, index):
        return self.vec[self._indexfun(self.dirfun(index))]
    def getcol(self, type=np.float32):
        return np.asarray([self[index] for index in xrange(12*self.nside**2)],type)

class HEALPixFITS(list):

    def __init__(self, cols, nside=None):
        """ initialize
        """
        assert len(cols)>0, 'no HParray object list'
        self.nside = nside if nside is not None else cols[0].nside
        for col in cols:
            self.append( col if col.nside==self.nside else HPresample(col, self.nside) )
            print 'appended column %s' %col.name
            
    def write(self, outfile, unit=None):
        makecol = lambda v: pyfits.Column(name=v.name, format='E', unit=unit, array=v.getcol(np.float32))
        cols = map(makecol, self)
        nside = self.nside
        cards = [pyfits.Card(*pars) for pars in [ 
                ('PIXTYPE',  'HEALPIX',          'Pixel algorithm',),
                ('ORDERING', 'RING',             'Ordering scheme'),
                ('NSIDE' ,    nside,        'Resolution Parameter'),
                ('NPIX',     12*nside**2,   '# of pixels'),
                ('FIRSTPIX',  0,                 'First pixel (0 based)'),
                ('LASTPIX',  12*nside**2-1, 'Last pixel (0 based)')]]
        table = pyfits.new_table(cols, header=pyfits.Header(cards))
        table.name = 'healpix' 
        hdus =  [pyfits.PrimaryHDU(header=None),  #primary
                 table,      # this table
                ]
        if os.path.exists(outfile):
            os.remove(outfile)
        pyfits.HDUList(hdus).writeto(outfile)
        print '\nwrote FITS file to %s' % outfile
   
 
def tables(outdir='uw40', names=('kde', 'ts')):
    """ create list of HParray guys from ROI tables
    """
    return [HPtables(name, outdir) for name in names]
    
    
def diffusefun(gal ='ring_24month_P74_v1.fits', energy=1000, nside=512):
    """ an HParray object wrapping the galactic diffuse """
    f = os.path.expandvars(os.path.join('$FERMI', 'diffuse', gal))
    t = skymaps.DiffuseFunction(f)
    def skyplotfun(skydir):
        return t.value(skydir, energy)
    return HPskyfun('galactic', skyplotfun, nside=nside)
 
def make_fits(outdir='uw40', outfile='../../pivot/24M7_uw40/aladin512.fits'):
    t = tables()
    d = diffusefun()
    HEALPixFITS(t+[d]).write(outfile) 
 
 
# this crashes now, no idea  
#def exposure(ltcube = 'P7_V4_SOURCE/24M7_lt.fits',
#        irf='P7SOURCE_V4PSF', energy=1000, nside=384):
#    """
#    return and HParray object wrapping the exposure
#    """
#    caldb = pycaldb.CALDBManager(irf)
#    skymaps.Exposure.set_cutoff(0.4)
#    skymaps.EffectiveArea.set_CALDB(caldb.CALDB)
#
#    f = os.path.expandvars(os.path.join('$FERMI','data', ltcube))
#    lc = skymaps.LivetimeCube(f) 
#    #ea  = [skymaps.EffectiveArea('',f) for f in caldb.get_aeff()]
#    aeff = map(lambda f: skymaps.EffectiveArea('',f), caldb.get_aeff())
#    #exposure = [skymaps.Exposure(lc,x) for x in ea]
#    exposure = map(lambda x: skymaps.Exposure(lc,x), aeff)
#    def skyplotfun(skydir):
#        return sum([exp.value(skydir, energy) for exp in exposure])
#    return skyplotfun #HPskyfun('exposure', skyplotfun, nside=nside)
#

    
