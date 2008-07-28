import numpy as N
import pyfits as PF
from types import ListType,FunctionType

try:
   import uw.pointlike
except: pass

import pointlike as pl
from pointlike import SkyDir


def phase_circ(eventfiles,ra=128.83646,dec=-45.17658,radius=6,phaseranges=[[0,1]]):
    """Construct a histogram of events interior to a circular region."""

    if not type(radius) is FunctionType:
        rval=radius
        radius=lambda e: e/e*rval

    if not type(eventfiles) is ListType: eventfiles=[eventfiles]

    center=SkyDir(ra,dec)

    ef = [PF.open(e) for e in eventfiles]
    ra = [N.asarray(e['EVENTS'].data.field('RA')).astype(float) for e in ef]
    dec = [N.asarray(e['EVENTS'].data.field('DEC')).astype(float) for e in ef]
    ph = [N.asarray(e['EVENTS'].data.field('PULSE_PHASE')).astype(float) for e in ef]
    en = [N.asarray(e['EVENTS'].data.field('ENERGY')).astype(float) for e in ef]
    for e in ef: e.close()

    radii = [radius(e) for e in en]

    #Fast, quick and dirty cut
    mask = [N.abs(dec[i]-center.dec())<=radii[i] for i in xrange(len(dec))]  
    ra = [ra[i][mask[i]] for i in xrange(len(ra))]
    dec = [dec[i][mask[i]] for i in xrange(len(dec))]
    ph = [ph[i][mask[i]] for i in xrange(len(ph))]


    #Slow, accurate cut
    scale=180/N.pi
    mask = [ N.array( [center.difference(SkyDir(ra[i][j],dec[i][j]))*scale<=radii[i][j]  \
                for j in xrange(len(ra[i])) ]) for i in xrange(len(ra)) ]
    ra = [ra[i][mask[i]] for i in xrange(len(ra))]
    dec = [dec[i][mask[i]] for i in xrange(len(dec))]
    ph = N.squeeze([ph[i][mask[i]] for i in xrange(len(ph))])
    mask=N.array([False]*len(ph))
    for r in phaseranges:
        for i in xrange(len(ph)):
            if r[0]<= ph[i] and ph[i] < r[1]: mask[i]=True                    
    return ph[mask]

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