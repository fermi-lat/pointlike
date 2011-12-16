"""
Check the residual TS maps for clusters
$Header$

"""

import os,pickle
import pyfits
from skymaps import Band, SkyDir
import numpy as np
import pylab as plt
from pointlike import IntVector

band = Band(256)
def sdir(index):
    return band.dir(int(index))

class TSdata(object):
    def __init__(self, outdir=None, filename='aladin256.fits', fieldname='ts'):
        if outdir is None:
            outdir = 'uw%02d' % int(open('version.txt').read())
        if not os.path.exists(os.path.join(outdir, filename)):
            project= os.path.split(os.getcwd())[-1]
            outdir = os.path.join('../..', 'pivot', '%s_%s' %(project, outdir))
        print 'using outdir %s' % outdir
        self.rts = pyfits.open(os.path.join(outdir,filename))[1].data.field(fieldname)
        self.glat= np.array([sdir(i).b() for i in range(len(self.rts))])
        self.glon= np.array([sdir(i).l() for i in range(len(self.rts))])
    def select(self, ts_min=0, b_min=0):
        cut = (self.rts>ts_min)* (abs(self.glat)>b_min)
        return self.rts[cut]
    def indices(self, ts_min=0, b_min=0):
        cut = (self.rts>ts_min)* (abs(self.glat)>b_min)
        return np.arange(len(self.rts))[cut]
        
def rts_hist(tsdata, fignum=1, bcut=10, bins=np.linspace(0,100,51)):
    plt.close(fignum);
    plt.figure(fignum)
    plt.hist(tsdata.select(),           bins, log=True, label='all')
    plt.hist(tsdata.select(b_min=bcut), bins, log=True, label='|b|>%d'%bcut)
    plt.legend(); plt.grid(True)
    plt.xlabel('TS')

def neighbors(i,j):
    iv = IntVector()
    n = band.findNeighbors(int(i),iv)
    return j in iv

def grow(indeces):
    i = 0
    cluster =list(indeces[0:1])
    remainder = indeces[1:]
    while i<len(cluster):
        newremain = []
        for j in remainder:
            if neighbors(cluster[i],j): cluster.append(j)
            else: newremain.append(j)
        i += 1
        remainder = newremain
    return cluster, remainder
        
def cluster(indices):
    ret = []
    rem = indices
    while len(rem)>0:
        clu, rem = grow(rem)
        ret.append(clu)
    return ret

class Cluster(object):
    def __init__(self, rts):
        self.rts = rts
    
    def group(self, clu):
        """ clu: list of pixel indices"""
        ts = np.array([self.rts[i] for i in clu])
        #print clu, ts
        dirs = [sdir(i) for i in clu]
        ra   = np.array([s.ra() for s in dirs])
        #print ra
        dec  = np.array([s.dec() for s in dirs])
        #print dec
        wra = sum(ts*ra)/sum(ts)
        wdec = sum(ts*dec)/sum(ts)
        self.ts = ts.max()
        self.sdir = SkyDir(float(wra), float(wdec))
        #print self.sdir
        
def make_seeds(tsdata, rcut=10, bcut=5, out=None, rec=None, seedroot='SEED-46'):
    """
    tsdata: object created by TSdata
    """
    indices  = tsdata.indices(rcut,bcut)
    clusters = cluster(indices)
    cl = Cluster(tsdata.rts)
    if out is not None: print >> out,  '# Region file format: DS9 version 4.0 global color=green'
    if rec is not None:
        print >> rec, 'name\t ra\t dec\t ts\t size\tl\tb'
    for i,x in enumerate(clusters):
        cl.group(x)
        if out is not None: 
            print >>out,'fk5; point(%8.3f, %8.3f) # point=cross text={%d:%d %.1f}'%\
                ( cl.sdir.ra(), cl.sdir.dec(),i, len(x), cl.ts)
        if rec is not None:
            print >>rec, '%s-%02d \t%8.3f \t%8.3f\t %8.1f\t%8d\t%8.3f \t%8.3f ' %\
                (seedroot, i,cl.sdir.ra(), cl.sdir.dec(),  cl.ts, len(x), cl.sdir.l(),cl.sdir.b())
        
    if rec is not None: rec.close()
    if out is not None: out.close()
    

    
