"""

Process the HEALpix ROI information into a pivot collection
$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/pipeline/roi_pivot.py,v 1.5 2011/01/01 15:50:05 burnett Exp $
"""
import os, glob, pickle
import numpy as np
from uw.thb_roi import pivot
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
def table_array(name, outdir):
    files =glob.glob(os.path.join(outdir, '%s_table'%name,'*.pickle'))
    nf = len(files)
    assert nf>0, 'no pickle files found in %s' % os.path.join(outdir, '%s_table'%name)
    if nf<1728: print 'warning: missing %d files' % (1728-nf)
    files.sort()

    pklist = np.array([pickle.load(open(f)) for f in files])
    return pklist

    
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
    for cname in 'galnorm galindex isonorm loglike chisq'.split():
        p.add_facet(cname, 'Number', 'F2', z.field(cname))
        
    # get the TS map arrays
    pks = table_array('ts', outdir)
    maxes = np.array([pks[i,:].max() for i in indices])
    p.add_facet('maximum TS','Number', 'F1', maxes)

    related = [[['sources', 'sources.cxml#healpix_12=EQ.%d'%(index)]] for index in range(len(names))]
    p.add_related(related)
    p.write(pivot_file)
  

