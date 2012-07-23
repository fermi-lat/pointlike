"""
A module for establishing fractional upper limits on modulation for a
light curve *with calibrated photon weights*.  Uses include upper
limits on pulsars and on orbital modulation eclipsing MSP binaries.

The basic idea is to proceed by Monte Carlo.  For each iteration, psuedo-
random numbers are drawn from a uniform distribution and compared to the
weights to determine if a given photon is "from the source" or "from the
background".

This is accomplished with the usual LC* machinery, which also requires
a template, which the user should supply in the form of, e.g., a typical
pulsar light curve, a notch, sinusoidal modulation, etc.

Then, the amplitude of the modulated component is gradually increased
until the weighted H-test (or other statistic) surprasses some threshold
for some fraction, e.g. 2/3, of the simulated population.  This gives
the fractional detection threshold for orbital (or any) modulation.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/modtest.py,v 1.1 2012/03/14 00:46:25 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>
"""
import numpy as np
from stats import hmw
import lcprimitives
import lctemplate

def get_stat_dist(lct,weights,stat,n=100):
    """ Use the provided template, photon weights (CALIBRATED) and
        statistic (TS) to generate a distribution of the TS at the
        specified amplitude."""

    stats = np.empty(n)
    for i in xrange(n):
        phases = lct.random(len(weights),weights=weights)
        stats[i] = stat(phases,weights)
    return np.sort(stats)

def fill_grid(lct,weights,stat=hmw,gmin=1e-5,gmax=1-1e-5,ngrid=20,
              fraction=[0.05,0.5,0.95],n=100):
    """ Find the distribution of a pulsation test statistic over a grid of
        modulation normalizations.
        
        lct        -- instance of LCTemplate; the light curve shape
        weights    -- probability a photon comes from the pulsar/point source
        stat       -- [hmw] the pulsation test statistic to use
        gmin       -- [1e-5] minimum modulation to use
        gmax       -- [1-1e-5] maximum modulation to use
        fraction   -- positions in the cdf of the MC population
        n          -- [100] size of the Monte Carlo sample to use

        Returns: (grid,vals)
            grid   -- the modulation amplitudes tested
            vals   -- the specified [fraction] contours of the distribution
    """
    
    index = np.round(n*np.asarray(fraction)).astype(int)

    grid = np.linspace(1e-5,1-1e-5,ngrid)
    vals = np.empty([len(index),len(grid)])

    for i,g in enumerate(grid):
        lct.norms.set_total(g)
        s = get_stat_dist(lct,weights,stat,n)
        vals[:,i] = s[index]

    return grid,vals

def find_threshold(grid,vals,threshold,index=0,degree=2):
    """ Use polynomial fits to solve for the modulation fraction giving
        an appropriate confidence interval.

        grid/vals -- results of fill_grid
        threshold -- the stat value for which to compute the interval
        index     -- index into vals
    """
    v = vals[index,:]-threshold
    if not np.any(v>0):
        print 'No modulation value satisfies threshold.'
        return None
    p = np.polyfit(grid,v,degree)
    if degree==2:
        a2,a1,a0 = p
        s1 = (-a1 + (a1**2-4*a0*a2)**0.5)/(2*a2)
        s2 = (-a1 - (a1**2-4*a0*a2)**0.5)/(2*a2)
    dom = np.linspace(grid[0],grid[-1],1001)
    cod = np.polyval(p,dom)
    idx = np.searchsorted(cod,0)
    return dom[idx:idx+2].sum()/2

def test_sinusoid(weights,order=1):
    """ Calculate limits on (by default) sinusoidal modulation.  Higher
        harmonics are supported by the [order] kwarg."""
    p0 = lcprimitives.LCHarmonic(order=order)
    lct = lctemplate.LCTemplate([p0],[1])
    grid,vals = fill_grid(lct,weights)
    return grid,vals

def test_notch(weights,width=0.05):
    """ Test for a dropout of gamma-ray emission with a notch shape."""
    p0 = lcprimitives.LCTopHat(p=[1-width,0])
    lct = lctemplate.LCTemplate([p0],[1])
    grid,vals = fill_grid(lct,weights)
    return grid,vals

def test_pulsar(weights):
    """ Look for a canonical 2peak+bridge pulsar light curve."""
    x2 = 0.55 
    lor = True
    bfrac = 0.1
    skew = False
    lct = lctemplate.get_gauss2(pulse_frac=0.25,bridge_frac=bfrac,lorentzian=lor,skew=skew,x1=0.2,x2=x2,ratio=1.1,width1=0.02,width2=0.02)
    grid,vals = fill_grid(lct,weights)
    return grid,vals

import pylab as pl
import pyfits
f = pyfits.open('/edata/single_sources/j1124m3653/gamma_products/j1124m3653-ft1_gtselect_gtmktime_r2.fits')
weights = f[1].data.field('WEIGHT')
f.close()

g1,v1 = test_sinusoid(weights)
g2,v2 = test_notch(weights)
g3,v3 = test_pulsar(weights)

# display results and print threshold
for v in [v1,v2,v3]:
    pl.errorbar(x=g1,y=v[1],yerr=[v[1]-v[0],v[2]-v[1]])
    print find_threshold(g1,v,15)

