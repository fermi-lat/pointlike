"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_analysis.py,v 1.84 2011/04/20 00:36:30 lande Exp $

handle binning of data for extraction of TOAs

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""
import numpy as np
from stats import hm,hmw
import sys

class TOABinner(object):
    """ Implement an interface to calculate bin edges for a span of LAT
        observations.  Note the default implementation is uniform in time."""

    def __init__(self,ntoa,phase_data):
        self.counter = 0
        self.ntoa = ntoa
        self.data = phase_data
        self.make_edges()

    def __iter__(self):
        self.counter = 0
        return self

    def next(self):
        if self.counter >= self.ntoa: raise StopIteration
        idx,self.counter = self.counter,self.counter+1
        return self.starts[idx],self.stops[idx]

    def make_edges(self):
        bin_edges = np.linspace(self.data.mjd_start,self.data.mjd_stop,self.ntoa+1)
        self.starts = bin_edges[:-1]
        self.stops  = bin_edges[1:]
            
class UniformSigBinner(TOABinner):
    """ Attempt to select bin edges such that each bin has approximately the same pulsed
        significance as measured by the H test.

        TODO -- would be better to use a Z^2_m test to reduce fluctuations
        TODO -- a better approach in general is to pick a goal for the H (Z^2)
            test and pick the number of bins accordingly, rather than vice
            versa as current implementation.
    """
    
    def make_edges(self):
        phases,weights,mjds = self.data.ph,self.data.weights,self.data.mjds
        ntoa = self.ntoa
        def get_hs():
            if weights is None:    
                return np.asarray([hm(phases[i0[i][0]:i0[i][1]]) for i in xrange(ntoa)])
            else:
                return np.asarray([hmw(phases[i0[i][0]:i0[i][1]],weights[i0[i][0]:i0[i][1]]) for i in xrange(ntoa)])
        e0 = np.linspace(self.data.mjd_start,self.data.mjd_stop,self.ntoa+1)
        idx = np.arange(len(phases))
        i0 = []
        for i in xrange(ntoa):
            indices = idx[(mjds>e0[i])&(mjds<e0[i+1])]
            i0 += [ [indices[0],indices[-1]] ]
        print [x[1]-x[0] for x in i0]
        init_hs = get_hs()
        print 'Initial set of H-test values with std dev. %.2f'%(init_hs.std())
        print init_hs
        shift = max(1,int(float(len(phases))/ntoa*0.02))
        def do_shift(my_shift,idx):
            if idx == 0:
                i0[idx][1] += 2*my_shift
                i0[idx+1][0] += 2*my_shift
            elif idx == ntoa-1:
                i0[idx][0] -= 2*my_shift
                i0[idx-1][1] -= 2*my_shift
            else:
                i0[idx][0] -= my_shift
                i0[idx][1] += my_shift
                i0[idx-1][1] -= my_shift
                i0[idx+1][0] += my_shift
        goal = 0.1*init_hs.mean()
        old_std = 0
        old_std2 = 0
        def check_indices():
            for idx in i0:
                if idx[1]-idx[0] <= 10: return False 
            return True

        for j in xrange(300):
            hs = get_hs()
            std = hs.std()
            if std < goal:
                break
            #print std;sys.stdout.flush()
            if (std == old_std) or (std == old_std2):
                print 'Detected loop'
                shift -= 1
                if shift==0: shift = 2
                shift = max(1,shift)
            old_std2,old_std = old_std,std
            m = hs.mean()
            for i,h in enumerate(hs):
                my_shift = shift if (h < m) else -shift
                do_shift(my_shift,i)
                if not check_indices(): do_shift(-my_shift,i)
                ph = phases[i0[i][0]:i0[i][1]]
                nh = hm(phases[i0[i][0]:i0[i][1]])
                if ((h < m) and ((nh+1) < h)) or ((h >= m) and ((nh-1) > h)):
                    #print i,h,nh,m
                    do_shift(-my_shift,i)
        print 'Final set of H-test values with std. dev %.2f'%(hs.std())
        print hs

        # construct MJD bounds
        edges = np.empty(ntoa+1)
        edges[0] = self.data.mjd_start
        edges[-1] = self.data.mjd_stop
        for i in xrange(ntoa-1):
            edges[i+1] = (mjds[i0[i][1]] + mjds[i0[i+1][0]])/2
        #return edges
        print edges
        self.starts = edges[:-1]; self.stops = edges[1:]
