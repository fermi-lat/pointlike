"""
Manage plotting of the band energy flux

"""
from uw.like import roi_bands
from uw.utilities import makerec
import numpy as np
import pylab as plt

    
class BandFlux(object):

    def __init__(self, roi, which=0, merge=True):
        """ extracted from roi_plotting
            roi:   has a list of ROIEnergyBand objects (only dependence)
            which: 
            merge: flag to merge adjacent upper and lower bands with upper limits only
        """
        self.which = which
        roi.setup_energy_bands()
        self.bands = roi.energy_bands
        centers = np.array([(b.emin*b.emax)**0.5  for b in self.bands])

        for eb in self.bands:
            eb.bandFit(which=which)
            eb.merged = False

        if merge:
            # count in from high end to find how many high energy bins have only upper limits
            for nhi,neb in enumerate(self.bands[::-1]):
                if neb.flux is not None: break

            # count in from low end to find how many low energy bins have only upper limits
            for nlo,neb in enumerate(self.bands):
                if neb.flux is not None: break
        else:
            nlo = nhi = 0
        
        nt = len(self.bands)
        eblo = [self.merge(0,nlo)]  if nlo>0 else []
        ebhi = [self.merge(nt-nhi, nt)] if nhi>0 else []

        # the truncated list of bands
        bands = eblo + self.bands[nlo:nt-nhi] + ebhi

        # Save info for plots or print in a recarray
        rec = makerec.RecArray('elow ehigh flux lflux uflux'.split())
        for eb in bands:

            xlo,xhi = eb.emin,eb.emax
            fac = xlo*xhi 
            bc  = (fac)**0.5
            hw  = (xhi-xlo)/2.
            if eb.merged: hw /= eb.num
            
            if eb.flux is not None: 
                rec.append(xlo, xhi, fac*eb.flux, fac*eb.lflux, fac*eb.uflux)
            else:
                rec.append(xlo, xhi, 0, 0, fac*eb.uflux)
            
        self.rec = rec()
    
    def merge(self, n1, n2):
        """merge the set of bands into a new ROIEnergyBand
            n1, n2 -- band indices to mege
        """
        hbands = []
        ebands = self.bands[n1:n2]
        for eb in ebands:
            for band in eb.bands:
                hbands.append(band)
        ebmerged = roi_bands.ROIEnergyBand(hbands, ebands[0].emin, ebands[-1].emax)
        ebmerged.bandFit(which=self.which) # a new fit
        ebmerged.merged = True
        ebmerged.num = len(ebands)
        return ebmerged

        
    
    def __str__(self):
        n  = len(self.rec.dtype)-2
        return ((2*'%10s'+' '+n*'%-12s'+'\n') % self.rec.dtype.names)\
             +'\n'.join( [(2*'%10.0f'+n*'%12.3e') % tuple(row) for row in self.rec])
        
    def plot(self, fignum=5, axes= None, axis=None, **kwargs):
        
        if axes is None:
            plt.figure(fignum); plt.clf()
            axes = plt.gca()
        axes.set_xscale('log')
        axes.set_yscale('log')
        axes.grid(True)
        
        if 'color' not in kwargs:
            color=kwargs['color'] = 'k'
        
        for r in self.rec:
            xl, xh = r.elow, r.ehigh
            bc = (xl*xh)**0.5
            if r.flux >0:
                axes.plot([xl,xh], [r.flux, r.flux],  **kwargs)
                axes.plot([bc,bc], [r.lflux,r.uflux], **kwargs)
            else:
                x,y = bc, r.uflux
                axes.plot([xl,xh], [y,y] , **kwargs) # bar at upper limit
                # plot arrow 0.6 long by 0.4 wide, triangular head (in log coords)
                axes.plot([x, x,     x*1.2, x,     x/1.2, x],
                          [y, y*0.6, y*0.6, y*0.4, y*0.6, y*0.6], **kwargs)
 
        if axis is not None:
            axes.axis(axis)
            axes.set_autoscale_on(False)
        return axes                        
                                