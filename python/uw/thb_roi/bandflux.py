"""
Manage plotting of the band energy flux

"""
from uw.like import roi_bands
from uw.utilities import makerec
import numpy as np
import pylab as plt

    
class BandFlux(object):

    def __init__(self, roi, which=0, merge=True, scale_factor=1.0):
        """ extracted from roi_plotting
            roi:   has a list of ROIEnergyBand objects (only dependence)
            which: 
            merge: flag to merge adjacent upper and lower bands with upper limits only
            scale_factor [1.0] used to scale flux units
        """
        self.which = which
        roi.setup_energy_bands()
        self.bands = roi.energy_bands
        self.scale_factor = scale_factor 
        centers = np.array([(b.emin*b.emax)**0.5  for b in self.bands])

        for eb in self.bands:
            eb.bandFit(which=which)
            eb.merged = False

        
        if merge:
            # this function decides if there is a measurement, meaning the lower flux is finite
            checkbound = lambda x: x.lflux> 1e-15 
            # count in from high end to find how many high energy bins have no measurements
            for nhi,neb in enumerate(self.bands[::-1]):
                if checkbound(neb): break

            # count in from low end to find how many low energy bins have no measurements
            for nlo,neb in enumerate(self.bands):
                if checkbound(neb): break
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
            fac = xlo*xhi * scale_factor # note the scale is here, perhaps to ergs
            bc  = (xlo*xhi)**0.5
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
        ebmerged = roi_bands.ROIEnergyBand(hbands)
        ebmerged.bandFit(which=self.which) # a new fit
        ebmerged.merged = True
        ebmerged.num = len(ebands)
        return ebmerged

    def __call__(self, emin, emax):
        """ call: return total ROIEnergyBand in given range
        
        """
        hbands = []
        num = 0
        for eb in self.bands:
            if eb.emin>emin/1.01 and eb.emax<emax*1.01:
                num +=1
                for band in eb.bands:
                    hbands.append(band)
        ebmerged = roi_bands.ROIEnergyBand(hbands)
        ebmerged.bandFit(which=self.which) # a new fit
        ebmerged.merged = True
        ebmerged.num = num
        return ebmerged
        
    
    def __str__(self):
        n  = len(self.rec.dtype)-2
        return ((2*'%10s'+' '+n*'%-12s'+'\n') % self.rec.dtype.names)\
             +'\n'.join( [(2*'%10.0f'+n*'%12.3e') % tuple(row) for row in self.rec])
        
    def plot_data(self, axes, **kwargs):
        
        if 'color' not in kwargs:
            kwargs['color'] = 'k'
        
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
 
                      
    def plot_model(self, axes, m, dom,  butterfly, **kwargs):
        """ 
            m: the model, implements MOdels.Model
            dom: the domain, a set of points
        """
        stat = m.statistical()
        err = stat[0]*stat[1]
        energy_flux_factor = self.scale_factor
        axes.errorbar([m.e0], [energy_flux_factor*stat[0][0]*m.e0**2], fmt='or', 
                yerr=energy_flux_factor*err[0]*m.e0**2, elinewidth=2, markersize=8)

        axes.plot( dom, energy_flux_factor*m(dom)*dom**2, **kwargs)
        if butterfly:
            # 'butterfly' region
            dom_r = np.array([dom[-i-1] for i in range(len(dom))]) #crude reversal.
            a,gamma = stat[0]
            var = err**2
            # r is e/e0
            bfun = lambda r: r**-gamma * np.sqrt(var[0] + (a*np.log(r))**2 * var[1])
            upper = energy_flux_factor*(m(dom)  + bfun(dom/m.e0)  )*dom**2
            lower = energy_flux_factor*(m(dom_r)/(1 +bfun(dom_r/m.e0)/m(dom_r)))*dom_r**2
            ymin, ymax = plt.gca().get_ylim()
            lower[lower<ymin] = ymin
            upper[upper>ymax] = ymax
            t =axes.fill(np.hstack( [dom,   dom_r] ), 
                        np.hstack( [upper, lower] ), 'r')
            t[0].set_alpha(0.4)
               
     
    