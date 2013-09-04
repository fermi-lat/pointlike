"""
Check the residual TS maps for clusters
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/check_ts.py,v 1.2 2013/01/20 14:06:04 burnett Exp $

"""

import os,pickle, argparse
import pyfits
from skymaps import Band, SkyDir
import numpy as np
import pylab as plt
from pointlike import IntVector
nside=512
band = Band(nside)
def sdir(index):
    return band.dir(int(index))

class TSdata(object):
    def __init__(self, outdir, filename, fieldname='ts'):
        assert os.path.exists(os.path.join(outdir, filename)), 'outdir, %s, not found' % outdir
        self.rts = pyfits.open(os.path.join(outdir,filename))[1].data.field(fieldname)
        assert len(self.rts)==12*(nside)**2, 'wrong nside in file %s: expect %d' %(filename, nside)
        self.glat= np.array([sdir(i).b() for i in range(len(self.rts))])
        self.glon= np.array([sdir(i).l() for i in range(len(self.rts))])
        print 'read %s ok' %filename
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
        
def make_seeds(tsdata, rcut=10, bcut=5, out=None, rec=None, seedroot='SEED-46', minsize=2):
    """
    tsdata: object created by TSdata
    rec: open file to write tab-delimited file to
    """
    indices  = tsdata.indices(rcut,bcut)
    clusters = cluster(indices)
    cl = Cluster(tsdata.rts)
    if out is not None: print >> out,  '# Region file format: DS9 version 4.0 global color=green'
    if rec is not None:
        print >> rec, 'name\tra\tdec\tts\tsize\tl\tb'
    for i,x in enumerate(clusters):
        if len(x)<minsize: continue
        cl.group(x)
        if out is not None: 
            print >>out,'fk5; point(%8.3f, %8.3f) # point=cross text={%d:%d %.1f}'%\
                ( cl.sdir.ra(), cl.sdir.dec(),i, len(x), cl.ts)
        if rec is not None:
            print >>rec, '%s-%03d\t%8.3f \t%8.3f\t %8.1f\t%8d\t%8.3f \t%8.3f ' %\
                (seedroot, i,cl.sdir.ra(), cl.sdir.dec(),  cl.ts, len(x), cl.sdir.l(),cl.sdir.b())
        
    if rec is not None: rec.close()
    if out is not None: out.close()
 
def main(args):
    print args
    global nside
    nside=args.nside
    assert os.path.exists(args.files[0]), 'did not find file %s'%args.files[0]
    tsdata = TSdata('.', args.files[0], args.tsfield)
    rec = open(args.files[1], 'w')
    make_seeds(tsdata, rcut=args.tsmin, bcut=float(args.bmin), rec=rec, seedroot=args.seedroot, minsize=args.minsize)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='examine a TS map file for clusters, create a file of source candidates')
    parser.add_argument('files', nargs=2, help='input FITS file and output text file')
    parser.add_argument('--tsmin',    help='minimum TS for cluster', default=10)
    parser.add_argument('--bmin',     help='minimum |b|',            default=5)
    parser.add_argument('--minsize',  help='minimum cluster size',   default=2)
    parser.add_argument('--seedroot', help='root for seed names',    default='SEED')
    parser.add_argument('--nside',    help='nside for healpix map',  default=512)
    parser.add_argument('--tsfield',  help='name of field with TS data',  default='ts')
    

    args = parser.parse_args()
    main(args)


    
