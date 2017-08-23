"""
Manage creation of a source Pivot collection from a set of sources
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pub/source_pivot.py,v 1.2 2012/08/13 19:52:23 burnett Exp $

"""

import os, pickle
import astropy.io.fits as pyfits
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
        source_names = None,
        variability=None,):
    """make the Pivot collection file for source
        source_names : string or None
            if a string, assume a FITS catalog file and check for NickName -> Source_Name
        variability: None, or dictionary of dictionaries
            keys: source name, varindex and cprob
    """
    names = z.name
    assert os.path.exists(pivot_dir), 'pivot directory %s does not exist' %s
    p = pivot.Pivot(z, pivot_dir, pivot_name, dzc)
    # add additional Facets here
   
    if source_names is not None:
        table = pyfits.open(source_names)[1].data
        sdict=dict(zip(list(table.NickName) ,list(table.Source_Name)))
        # note that the NickName column has blacks compressed
        sourcename = [sdict.get(n.replace(' ',''), '(none)') for n in z.name]
        p.add_facet('Source_Name', 'String', 'C', sourcename) 
        print 'added Source_Name facet from FITS file %s' % source_names
 
    try:  p.add_facet('beta', 'Number', 'F3',  p.limit(z.beta, 0,2.5, nan=0))
    except:
        print 'no beta found in the source rec array'
        pass
    p.add_facet('SpectrumType', 'String', 'C', z.modelname)
    
    def source_name(name, extended):
        if extended: return 'Extended'
        try:
            return {'1F':'1FGL', 'PG':'PGW', 'MR':'MRF', 'UW':'UW', 'MS':'MST',
                    'SE':'SEED', '24':'24M',
                    'Cy':'bin', 'LS':'bin','PS':'PSR','gc':'CG source'}[name[:2]] 
        except:
            return name[:2]
        
    p.add_facet('source', 'String', 'C', map(source_name, p.rec.name, z.extended))
    p.add_facet('extended', 'String','C', z.extended)
    ts = z.ts
    ts[ts<0]=0
    p.add_facet('fit quality', 'Number', 'F1', p.limit(z.band_ts-ts,0,100))
    good = (ts>10) * (p.rec.qual<25)
    try:    p.add_facet('10 GeV band TS', 'Number', 'F1', z.bts10)
    except: pass
    # out of memory??
    #p.add_facet('cutoff energy', 'Number', 'F1', p.limit(z.cutoff, 0, 1e5, nan=-1))
    #p.add_facet('good', 'String', 'C', good)
    if variability is not None:
        v = variability
        varindex = np.array([v[name]['varindex'] for name in names])
        cprob = np.array([v[name]['cprob'] for name in names])
        p.add_facet('variability prob',  'Number',  'F3', p.limit(cprob, 0,1, nan=-1))
        #?p.add_facet('variability delta', 'Number', 'F3', p.limit(v.delta, -5, 5, nan=-6))
        p.add_facet('variability index', 'Number', 'F3', p.limit(varindex, 0, 1e3, nan=-1))
        
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
    
    