import numpy as N
import pyfits as PF
from types import ListType,FunctionType

try:
   import uw.pointlike
except: pass

import pointlike as pl
from pointlike import SkyDir

def phase_circ(eventfiles,center=None,radius=6,phaseranges=[[0,1]],\
   erange=None,event_class=-1,ctbclasslevel=None,drift_check=False):
    """Construct a histogram of events interior to a circular region."""

    if not type(radius) is FunctionType:
        rval=radius
        radius=lambda e: e/e*rval

    if not type(eventfiles) is ListType: eventfiles=[eventfiles]

    center=center or SkyDir(128.83646,-45.17658)

    ef = [PF.open(e,memmap=1) for e in eventfiles]
    ra = [N.asarray(e['EVENTS'].data.field('RA')).astype(float) for e in ef]
    dec = [N.asarray(e['EVENTS'].data.field('DEC')).astype(float) for e in ef]
    ph = [N.asarray(e['EVENTS'].data.field('PULSE_PHASE')).astype(float) for e in ef]
    en = [N.asarray(e['EVENTS'].data.field('ENERGY')).astype(float) for e in ef]
    ec = [N.asarray(e['EVENTS'].data.field('EVENT_CLASS')).astype(int) for e in ef]
    if ctbclasslevel is not None:
      ctb = [N.asarray(e['EVENTS'].data.field('CTBCLASSLEVEL')).astype(int) for e in ef]
    if drift_check:
      dc = [N.asarray(e['EVENTS'].data.field('TIME')).astype(float) for e in ef]
    for e in ef: e.close()

    radii = [radius(e) for e in en]

    #Fast, quick and dirty cut
    mask = [N.abs(dec[i]-center.dec())<=radii[i] for i in xrange(len(dec))]  
    ra = [ra[i][mask[i]] for i in xrange(len(ra))]
    dec = [dec[i][mask[i]] for i in xrange(len(dec))]
    ph = [ph[i][mask[i]] for i in xrange(len(ph))]
    en = [en[i][mask[i]] for i in xrange(len(en))]
    ec = [ec[i][mask[i]] for i in xrange(len(ec))]
    if ctbclasslevel is not None:
      ctb = [ctb[i][mask[i]] for i in xrange(len(ctb))]
    if drift_check:
      dc = [dc[i][mask[i]] for i in xrange(len(dc))]
    radii = [radii[i][mask[i]] for i in xrange(len(radii))]

    #Slow, accurate cut
    scale=180/N.pi
    mask = [ N.array( [center.difference(SkyDir(ra[i][j],dec[i][j]))*scale<=radii[i][j]  \
                for j in xrange(len(ra[i])) ]) for i in xrange(len(ra)) ]
    ra = [ra[i][mask[i]] for i in xrange(len(ra))]
    dec = [dec[i][mask[i]] for i in xrange(len(dec))]
    ph = N.concatenate([ph[i][mask[i]] for i in xrange(len(ph))])
    en = N.concatenate([en[i][mask[i]] for i in xrange(len(en))])
    ec = N.concatenate([ec[i][mask[i]] for i in xrange(len(ec))])
    if ctbclasslevel is not None:
      ctb = N.concatenate([ctb[i][mask[i]] for i in xrange(len(ctb))])
    if drift_check:
      dc = N.concatenate([dc[i][mask[i]] for i in xrange(len(dc))])
    mask=N.array([False]*len(ph))
    for r in phaseranges:
        for i in xrange(len(ph)):
            if r[0]<= ph[i] and ph[i] < r[1]: mask[i]=True
    if erange is not None:
        mask = N.logical_and( N.logical_and(en>=erange[0],en<erange[1]), mask)
    if event_class >= 0:
        mask = N.logical_and(ec==event_class,mask)
    if ctbclasslevel is not None:
        mask = N.logical_and(ctb>=ctbclasslevel,mask)
    
    if drift_check:
        return [ph[mask],dc[mask]] #Return phase and time
    return ph[mask] #Just return phase

def phase_cut(eventfile,outputfile=None,phaseranges=[[0,1]]):
    """Select phases within a set of intervals.
    
        outputfile - set to change the default output (eventfile_PHASECUT.fits)
        phaseranges - a set of ranges"""
    from numarray import array as narray
    ef=PF.open(eventfile)
    ph=N.asarray(ef['EVENTS'].data.field('PULSE_PHASE')).astype(float)
    mask=N.array([False]*len(ph))
    for r in phaseranges:
        for i in xrange(len(ph)):
            if r[0]<= ph[i] and ph[i] < r[1]: mask[i]=True                    
    hdu=PF.new_table(ef['EVENTS'].columns,nrows=int(mask.sum()))
    mask=narray(mask)
    for i in xrange(len(ef['EVENTS'].columns)):
        hdu.data.field(i)[:]=ef['EVENTS'].data.field(i)[mask]
    ef['EVENTS'].data=hdu.data
    if outputfile:
        ef.writeto(outputfile,clobber=True)
    else:
        ef.writeto(eventfile.replace('.fits','_PHASECUT.fits'),clobber=True)
    ef.close()

def constant_count_histogram(phases,photons_per_bin=100):
   """Return a set of bins and rates for a 'constant count' lightcurve.
      phases -- a list of photon phases"""
   from collections import deque
   counter = 0
   ev = N.asarray(phases)
   ev.sort()

   edges = deque()
   edges.append(0)
   for i,j in enumerate(ev):
      if counter == photons_per_bin-1:
         edges.append(j)
         counter = 0
      else: counter += 1
   edges.append(1)
   edges = N.asarray(edges)
   delta_phi = edges[1:]-edges[:-1]
   rates = photons_per_bin/delta_phi
   rates[-1] = counter/delta_phi[-1]
   return [edges[:-1],rates,delta_phi] #Edges, rates, and widths, appropriate for a bar chart


#--------------------------------------------------------
    
if __name__=='__main__':
   a=Phase('ft1_first_diff_v2.fits')
   a.write()
   from pylab import hist
   hist(a.phaseHist(), bins = 20)
    

        
"""Same, but with periods.  This is harder, don't know why I bothered.
if pd==0:
    phases = [t/p for t in self.times] #Exact to 0th
elif pdd==0:
    phases = [1/pd*N.log(1.+pd/p*t) for t in self.times] #Exact to 1st
else:
    d=(4*pdd*p-pd**2)**0.5 #Discriminant
    phases = [2/d*(N.arctan( (pd+2*pdd*t)/d ) - N.arctan( pd/d )) for t in self.times] #Exact to 2nd
"""