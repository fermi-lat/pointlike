"""

Process the HEALpix ROI information into a pivot collection
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pub/roi_pivot.py,v 1.5 2011/07/22 13:55:16 burnett Exp $
"""
import os, glob, pickle
import numpy as np
from . import pivot
from skymaps import Band

def hp12_index(skydir):
    return Band(12).index(skydir)
def hp12_dir(index):
    return Band(12).dir(index)
def hpindex(sdir):
    return Band(12).index(sdir)
def hpdir(index):
    return Band(12).dir(index)
    
def hpname(index):
    return 'HP12_%04d'%index
    
def table_array_max(name, outdir):
    files =glob.glob(os.path.join(outdir, '%s_table'%name,'*.pickle'))
    nf = len(files)
    assert nf>0, 'no pickle files found in %s' % os.path.join(outdir, '%s_table'%name)
    if nf<1728: print ('warning: missing %d files' % (1728-nf))
    files.sort()
    maxes = map(lambda f:pickle.load(open(f)).max(), files)
    return maxes

    
def make_pivot(z, outdir,
        pivot_dir, pivot_name, 
        pivot_file='rois.cxml',
        dzc = 'dzc.xml',
        ):
    """make the Pivot collection file (this first will make a DeepZoom collection (dzc) if the folder does not exist 
    """
    names = z.name
    assert os.path.exists(pivot_dir), 'you need to create the images'
    if pivot_name is None: pivot_name = '%s - ts maps'%outdir
    
    # run the basic Pivot creation class, but ignore stuff assuming these are sources
    p = pivot.MultiPivot(z, [dzc], pivot_dir, pivot_name, fill_default_facets=False)

    # add additional Facets here
    indices = range(len(names))
    sdir = map( Band(12).dir, indices)
    p.add_facet('index', 'Number', 'F0', indices)
    p.add_facet('ra',  'Number', 'F1', [s.ra() for s in sdir])
    p.add_facet('dec', 'Number', 'F1', [s.dec() for s in sdir])
    b = np.array([s.b() for s in sdir])
    p.add_facet('glat','Number', 'F1', b)
    p.add_facet('glon','Number', 'F1', [s.l() for s in sdir])
    p.add_facet('High Latitude', 'String', 'C', np.abs(b)>10)
    for cname in 'galnorm galindex isonorm limbnorm loglike chisq'.split():
        p.add_facet(cname, 'Number', 'F2', z.field(cname))
        
    # get the TS map arrays for maximum if there
    try:
        tsmax= table_array_max('ts', outdir)
        p.add_facet('maximum TS','Number', 'F1', tsmax)
    except: pass

    related = [[['sources', 'sources.cxml#ROI_num=EQ.%d'%(index)]] for index in range(len(names))]
    p.add_related(related)
    p.write(pivot_file)
  

