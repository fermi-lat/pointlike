"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/edf.py,v 1.1 2011/04/27 18:32:03 kerrm Exp $

a module for implementing EDF (empirical distribution function) TOA determination; heavily beta

Authors: Matthew Kerr <matthew.kerr@gmail.com>
"""

import numpy as np
import pylab as pl
#from lcfitters import get_gauss2,get_gauss1
from scipy.interpolate import interp1d
from scipy.optimize import fmin

#lct = get_gauss2(pulse_frac=0.1)
#lct = get_gauss1(pulse_frac=0.1,width1=0.02)
#lct = get_gauss2(x1=0.25,x2=0.75,pulse_frac=0.1,width1=0.02,width2=0.02,ratio=1)
#lct = get_gauss2(x1=0.3,x2=0.7,pulse_frac=0.3,width1=0.15,width2=0.08,ratio=1)
#ph1 = lct.random(1000)
#ph2 = lct.random(1000)

def comp_edfs2(edf1,edf2):
    """ Calculate a Kuiper test statistic for two edfs."""
    idx1 = np.arange(1,edf1.n+1,dtype=float)/edf1.n
    idx2 = np.arange(1,edf2.n+1,dtype=float)/edf2.n
    s = np.argsort(np.append(edf1.ph,edf2.ph))
    idx = np.append(idx1,idx2)
    mask = np.append([True]*edf1.n,[False]*edf2.n)
    idx = idx[s]
    mask = mask[s]
    e1 = 0; e2 = 0
    all_e1s = np.empty_like(idx)
    all_e2s = np.empty_like(idx)
    for i in range(len(s)):
        if mask[i]:  e1 = idx[i]
        else:        e2 = idx[i]
        all_e1s[i] = e1
        all_e2s[i] = e2
    return np.max(all_e1s-all_e2s) + np.max(all_e2s-all_e1s)

def timeit(e1,e2,n=10):
    for i in range(n):
        comp_edfs2(e1,e2) 
    

class EDF(object):

    def __init__(self,ph,comp_n=100000):
        self.ph = np.sort(ph)
        self.n = len(self.ph)
        #self.comp_vals = self(np.linspace(0,1,comp_n))

    def __call__(self,phi):
        if not hasattr(phi,'__len__'):
            phi = [phi]
        return np.searchsorted(self.ph,phi).astype(float)/self.n

    def comp(self,e):
        # Kuiper test
        return comp_edfs2(self,e)
        #return np.max(self.comp_vals-e.comp_vals) + \
        #       np.max(e.comp_vals-self.comp_vals)

    def random(self,n):
        # NB -- this needs to be checked/optimized
        rvals = np.random.rand(n)
        y = np.linspace(0,1,self.n+1)
        #idx = np.searchsorted(y,rvals)
        idx = (rvals*self.n+1).astype(int)
        ph = np.append(self.ph,1+self.ph[0])
        rvals = self.n*((y[idx]-rvals)*ph[idx-1] + (rvals-y[idx-1])*ph[idx])
        return np.mod(rvals,1)
        
        
def make_cum_plot(ph,fignum=3):
    s = np.argsort(ph)
    num = np.linspace(0,1,len(ph)+1)
    num = (num[:-1]+num[1:])/2
    pl.figure(fignum);pl.clf()
    pl.plot(ph[s],num,marker='.')
    pl.plot([0,1],[0,1],color='k',ls='--')
    pl.xlabel('Phase')
    pl.ylabel('Empirical Distribution Function')


def find_alignment(e1,ph2,ngrid=101,hint=0,hint_range=.50):

    def pop_grid(grid):
        rvals = np.empty(len(grid),dtype=float)
        for ix,x in enumerate(grid):
            ph = np.mod(ph2+x,1)
            e2 = EDF(ph)
            rvals[ix] = e1.comp(e2)
        return rvals

    dom = np.linspace(hint-hint_range,hint+hint_range,ngrid)
    dcopy = pop_grid(dom)
    hint = dom[np.argmin(pop_grid(dom))]
    hint_range = hint_range/(float(ngrid-1)/2)
    dom = np.linspace(hint-hint_range,hint+hint_range,ngrid)
    hint = dom[np.argmin(pop_grid(dom))]
    return hint,dcopy
    
def estimate_diff(e1,ph2,mc=100,hint=0,hint_range=0.50):
    diff = find_alignment(e1,ph2,hint=hint,hint_range=hint_range)
    n = len(ph2)
    diffs = np.empty(mc)
    for i in range(mc):
        #diffs[i] = find_alignment(e1,e1.random(n),hint=hint)
        diffs[i] = find_alignment(e1,lct.random(n),hint=hint)
    return diff,diffs

        
