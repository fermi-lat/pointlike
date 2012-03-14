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

$Header: $

author: M. Kerr <matthew.kerr@gmail.com>
"""
import numpy as np
from stats import hmw
import lcprimitives
import lctemplate

def do_sims(lct,weights,stat,n=100):
    """ Use the provided template, photon weights (CALIBRATED) and
        statistic (TS) to generate a distribution of the TS at the
        specified amplitude."""

    stats = np.empty(n)
    for i in xrange(n):
        phases = lct.random(len(weights),weights=weights)
        stats[i] = stat(phases,weights)
    return np.sort(stats)

def find_threshold(lct,weights,stat=hmw,threshold=15,fraction=[0.05,0.5,0.95],n=100):
    """ Find the modulation threshold of stat by increasing the
        modulation amplitude until [fraction] of the simulated
        population of [n] surpass [threshold].  This is done by
        sampling the TS distribution at a series of modulation
        amplitudes and interpolating to [threshold]. """
    
    index = np.round(n*np.asarray(fraction)).astype(int)

    grid = np.linspace(1e-5,1-1e-5,20)
    vals = np.empty([len(grid),len(index)])

    for i,g in enumerate(grid):
        lct.norms.set_total(g)
        s = do_sims(lct,weights,stat,n)
        #v = np.searchsorted(s,threshold)
        #vals[i] = s[index]
        vals[i,:] = s[index]
        #if v < n:
            #vals[i] = s[v]
        #else:
            #vals[i] = s[-1] # "lower limit"

    return grid,vals

def test_sinusoid(weights,order=1):
    """ Calculate limits on (by default) sinusoidal modulation.  Higher
        harmonics are supported by the [order] kwarg."""
    p0 = lcprimitives.LCHarmonic(order=order)
    lct = lctemplate.LCTemplate([p0],[1])
    grid,vals = find_threshold(lct,weights)
    return grid,vals

def test_notch(weights,width=0.05):
    """ Test for a dropout of gamma-ray emission with a notch shape."""
    p0 = lcprimitives.LCTopHat(p=[1-width,0])
    lct = lctemplate.LCTemplate([p0],[1])
    grid,vals = find_threshold(lct,weights)
    return grid,vals

import pylab as pl

#weights = np.random.rand(500)
import pyfits
f = pyfits.open('/edata/single_sources/j1124m3653/gamma_products/j1124m3653-ft1_gtselect_gtmktime_r2.fits')
weights = f[1].data.field('WEIGHT')
f.close()
"""
x2 = 0.55 
lor = True
bfrac = 0.1
skew = False
lct = lctemplate.get_gauss2(pulse_frac=0.25,bridge_frac=bfrac,lorentzian=lor,skew=skew,x1=0.2,x2=x2,ratio=1.1,width1=0.02,width2=0.02)

grid,vals = find_threshold(lct,weights,hmw,threshold=50,n=100)
y = vals[:,1]
pl.errorbar(x=grid,y=vals[:,1],yerr=[y-vals[:,0],vals[:,2]-y])
"""
