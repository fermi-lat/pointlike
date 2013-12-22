"""
Main entry for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipe.py,v 1.42 2013/12/18 15:12:16 burnett Exp $
"""
import os, types, glob, pickle
import numpy as np
from uw.like2 import tools

def roirec(version, nside=12):
    roi_files = glob.glob('uw%02d/pickle/*.pickle'%version)
    roi_files.sort()
    n = 12*nside**2
    if len(roi_files)<n:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(n) if x not in t]
        print 'misssing roi files: %s' % missing
        raise Exception('misssing roi files: %s' % missing)
        
    recarray = tools.RecArray('name chisq loglike'.split())
    for fname in roi_files:
        p = pickle.load(open(fname))
        counts = p['counts']
        obs,mod = counts['observed'], counts['total']
        chisq =  ((obs-mod)**2/mod).sum()

        recarray.append(p['name'], chisq, p['logl'])
    return recarray()

def check_missing_files(folder):
    roi_files = sorted(glob.glob(os.path.join(folder, '*.pickle')))
    n = 1728
    missing = []
    if len(roi_files)<n:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(n) if x not in t]
        if len(missing)<10:
            print '\tmisssing roi files: %s' % missing
        else:
            print '\tMissing %d roi files: %s...' %( len(missing), missing[:10])
    return missing    



def roirec(outdir, check=False):
    roi_files = sorted(glob.glob(os.path.join(outdir, 'pickle/*.pickle')))
    n = 1728
    if len(roi_files)<n:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(n) if x not in t]
        if len(missing)<10:
            print '\tmisssing roi files: %s' % missing
        else:
            print '\tMissing %d roi files: %s...' %( len(missing), missing[:10])
        if check: return missing    
    if check: return[]#    return None # raise Exception('misssing roi files: %s' % missing)
        
    recarray = tools.RecArray('name chisq loglike prevlike niter'.split())
    bad = []
    for fname in roi_files:
        p = pickle.load(open(fname))
        if 'counts' not in p.keys():
            bad.append(fname)
            continue
        counts = p['counts']
        obs,mod = counts['observed'], counts['total']
        chisq =  ((obs-mod)**2/mod).sum()
        history = p.get('history', None)
        if history is not None:
            # new format
            logl = history[-1]['logl']
            histlen = len(history)
            prevlike = history[-2]['logl'] if histlen>1 else logl
        else:
            logl_list = p.get('prev_logl', [0])
            histlen = len( logl_list )
            logl = p.get('logl', 0)
            prevlike = logl_list[-1]
        recarray.append(p['name'], chisq, logl, prevlike, histlen)
    if len(bad)>0:
        print 'no fit info in file(s) %s' % bad if len(bad)<10 else (str(bad[:10])+'...')
    return recarray()

def check_converge(month, tol=10, add_neighbors=True, log=None):
    """ check for convergence, ROI that have been updated
    month: int or string
        if int, intrepret as a month, else a folder
        
    """
    from pointlike import IntVector
    from skymaps import Band
    outdir = 'month%02d'%month if type(month)==types.IntType else month
    #print '%s:' %outdir
    r = roirec(outdir)
    if r is None: return
    diff = r.loglike-r.prevlike
    dmin,dmax = diff.min(), diff.max()
    rmin,rmax = list(diff).index(dmin), list(diff).index(dmax)
    changed = set(np.arange(1728)[np.abs(diff)>tol])
    print >>log, '\tpass %d:  %d changed > %d, min, max: %d(#%d) %d(#%d)' % (max(r.niter),len(changed), tol, dmin,rmin,dmax,rmax),
    if not add_neighbors: return list(changed)
    nbrs = set()
    b12 = Band(12)
    for x in changed: 
        v = IntVector()
        b12.findNeighbors(int(x),v) # int is tricky
        for n in v:
            nbrs.add(n)
    q =list(changed.union( nbrs))
    print >>log, ' (total %d)' %len(q)
    if log is not None: log.flush()
    return q
    
       

class Create(object):
    pass
    
class Update(object): pass
class PulsarLimitTables(object):pass
class Finish(object): pass
class PulsarDetection(object):pass
class Tables(object): pass