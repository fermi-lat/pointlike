"""
Utilities for managing Healpix arrays
$Header$
"""
import os,glob,pickle
import pylab as plt
import numpy as np
from skymaps import Band, SkyDir, Hep3Vector, PySkyFunction, DiffuseFunction
from uw.utilities import image

def load_tables(name, outdir, nside=12, subnside=12*32):
    """
    
    """
    files =glob.glob(os.path.join(outdir, '%s_table'%name,'*.pickle'))
    nf = len(files)
    assert nf>0, 'no pickle files found in %s' % os.path.join(outdir, '%s_table'%name)
    if nf<1728: print 'warning: missing %d files' % (1728-nf)
    files.sort()

    r = np.zeros(12*subnside**2)
    r.fill(np.nan)
    pklist = [pickle.load(open(f)) for f in files]
    i12 = [int(f[-11:-7]) for f in files]
    index_table = make_index_table(nside, subnside)
    for index, pk in zip(i12,pklist):
        indeces = index_table[index]
        for i,v in enumerate(pk):
            r[indeces[i]]=v
            
    return r

def dump_table(skyfun, tablename, nside=384):
    """ create a table from the skyfun
        To adapt a FITS image file, set skyfun = DiffuseFunction(filename)
    """
    bdir = Band(nside).dir
    table = np.asarray([skyfun(bdir(index)) for index in xrange(12*nside**2)],float)
    pickle.dump(table, open(tablename, 'wb'))

def make_index_table(nside, subnside):
    band, subband = Band(nside), Band(subnside)
    npix, nsubpix = 12*nside**2, 12*subnside**2
    t=np.array([band.index(subband.dir(i)) for i in xrange(nsubpix)])
    a = np.arange(nsubpix)
    index_table = [a[t==i] for i in xrange(npix)]
    return index_table


def skyplot(crec, title='', axes=None, fignum=30, ait_kw={}, **kwargs):
    """ make an AIT skyplot of a HEALpix array
    crec : array
        must be sorted according to the HEALpix index
    title : string
        set the figure title
    ait_kw : dict
        to set kwargs for image.AIT, perhaps pixelsize
    
    Other args passed to imshow
    """
    n = len(crec)
    nside = int(np.sqrt(n/12))
    assert n==12*nside**2, 'wrong length to be healpix array'
    band = Band(nside)
    def skyplotfun(v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        index = band.index(skydir)
        return crec[index]
    if axes is None:
        plt.close(fignum)
        fig = plt.figure(fignum, figsize=(12,6))
    ait=image.AIT(PySkyFunction(skyplotfun) ,axes=axes, **ait_kw)
    ait.imshow(title=title, **kwargs)
    return ait
