"""
Manage creation of a source Pivot collection from a set of sources
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pub/source_pivot.py,v 1.2 2011/02/11 21:27:34 burnett Exp $

"""

import os, pickle
import numpy as np
from . import roi_pivot, pivot
from skymaps import SkyDir, Band 

def hpindex(sdir):
    return Band(12).index(sdir)
def hpdir(index):
    return Band(12).dir(index)
    
def hpname(index):
    return 'HP12_%04d'%index
 
   
def make_pivot(z, outdir,
        pivot_dir, pivot_name, 
        pivot_file='sources.cxml',
        dzc = 'dzc.xml',
        variability=False,):
    """make the Pivot collection file for source
    """
    names = z.name
    assert os.path.exists(pivot_dir), 'pivot directory %s does not exist' %s
    p = pivot.Pivot(z, pivot_dir, pivot_name, dzc)
    # add additional Facets here
    
    try:  p.add_facet('beta', 'Number', 'F3',  p.limit(z.beta, 0,2.5, nan=0))
    except:
        print 'no beta found in the source rec array'
        pass
    p.add_facet('SpectrumType', 'String', 'C', z.modelname)
    
    def source_name(name):
        try:
            return {'1F':'1FGL', 'PG':'PGW', 'MR':'MRF', 'UW':'UW', 'MS':'MST',
                    'SE':'SEED', '24':'24M',
                    'Cy':'bin', 'LS':'bin','PS':'PSR','gc':'CG source'}[name[:2]] 
        except:
            if name in ['IC443','W28','HESS J1825-137','W44','W51C','MSH 1552','Vela X','LMC','SMC']:
                return 'Extended'
            else: return name[:2]
        
    p.add_facet('source', 'String', 'C', map(source_name, p.rec.name))
    ts = z.ts
    ts[ts<0]=0
    p.add_facet('fit quality', 'Number', 'F1', p.limit(z.band_ts-ts,0,100))
    good = (ts>10) * (p.rec.qual<25)
    try:    p.add_facet('10 GeV band TS', 'Number', 'F1', z.bts10)
    except: pass
    # out of memory??
    #p.add_facet('cutoff energy', 'Number', 'F1', p.limit(z.cutoff, 0, 1e5, nan=-1))
    #p.add_facet('good', 'String', 'C', good)
    if variability:
        p.add_facet('variability prob', 'Number', 'F3', p.limit(z.cprob, 0,1, nan=-1))
        p.add_facet('variability delta', 'Number', 'F3', p.limit(z.delta, -5, 5, nan=-6))
        p.add_facet('variability index', 'Number', 'F3', p.limit(z.varindex, 0, 1e3, nan=-1))
        
    # ROI-related fields    
    sdir = map(SkyDir, p.rec.ra, p.rec.dec)
    hpindices = map(roi_pivot.hpindex, sdir)
    hpdirs = map(roi_pivot.hpdir, hpindices)
    p.add_facet('ROI_num', 'Number', 'I', hpindices, filter='true')
    #roi_dist = [np.degrees(dir.difference(roidir)) for dir, roidir in zip(sdir, hpdirs)]
    roi_dist = map( lambda a,b: np.degrees(a.difference(b)), sdir, hpdirs)
    p.add_facet('ROI_dist', 'Number', 'F2', roi_dist, filter='true')
    # and the related collection
    related = [[[roi_pivot.hpname(index), 'rois.cxml#index=EQ.%d'%(index)]] for index in hpindices]
    p.add_related(related)

    p.write(pivot_file)
    
    