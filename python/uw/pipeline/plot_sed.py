"""
Manage a SED plot
Two classes:
    SEDflux -- use the ROIAnalysis stuff to extract the flux per band
    PlotSED  -- from the flux per band, generate a SED plot
    
    To generate a plot from an roi object, given:
            roi: an ROIAnalysis object, 
            which: index of the source
    
        fluxes = SEDflux(roi, which).rec
        model = roi.psm.models[which]
        name  = roi.psm.names[which]
        PlotSED(fluxes)(model, name)
        
    Note that PlotSED is independent of the likelihood code, and the plot can be recreated
    from the fluxes, model and name

$Header$
"""
import os, sys
import numpy as np
import pylab as plt
from uw.utilities import makerec, image
from uw.like import roi_bands

class SEDflux(object):
    """ manage the band fluxes from a point source fit (extracted from roi_plotting)
    
    generates a recarray (member rec) with fields elow ehigh flux lflux uflux
    """

    def __init__(self, roi, which=0, merge=True, use_ergs=True):
        """ 
            roi:  ROIAnalysis object
                has a list of ROIEnergyBand objects (only dependence, except for extended stuff?)
            which: int or string
                index of the model (or a name? Josh code)
            merge: bool
                flag to merge adjacent upper and lower bands with upper limits only
            use_ergs : bool
                set False to revert to MeV units
                
        """
        self.which = which
        self.manager,self.index=roi.mapper(which)
        self.roi   = roi

        roi.setup_energy_bands()
        self.bands = roi.energy_bands
        scale_factor = 1.602e-6 if use_ergs else 1.0
        centers = np.array([(b.emin*b.emax)**0.5  for b in self.bands])

        for eb in self.bands:
            if self.manager == roi.psm:
                eb.bandFit(which=self.index)
            else:
                BandFitExtended(self.index,eb,self.roi).fit()

            eb.merged = False
        checkbound = lambda x: x.lflux> 1e-15 
        if merge:
            # this function decides if there is a measurement, meaning the lower flux is finite
            # count in from high end to find how many high energy bins have no measurements
            for nhi,neb in enumerate(self.bands[::-1]):
                if checkbound(neb): break

            # count in from low end to find how many low energy bins have no measurements
            for nlo,neb in enumerate(self.bands):
                if checkbound(neb): break
        else:
            nlo = nhi = 0
        nt = len(self.bands)
        if nlo != nt-1:
            eblo = [self.merge(0,nlo)]  if nlo>0 else []
            ebhi = [self.merge(nt-nhi, nt)] if nhi>0 else []

            # the truncated list of bands
            bands = eblo + self.bands[nlo:nt-nhi] + ebhi
        else:
            # Merge all bands
            bands = [self.merge(0,nt-1)]

        # Save info for plots or print in a recarray
        rec = makerec.RecArray('elow ehigh flux lflux uflux'.split())
        for eb in bands:

            xlo,xhi = eb.emin,eb.emax
            fac = xlo*xhi * scale_factor # note the scale is here, perhaps to ergs
            bc  = (xlo*xhi)**0.5
            hw  = (xhi-xlo)/2.
            if eb.merged: hw /= eb.num
            
            if checkbound(eb): 
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

        if self.manager == self.roi.psm:
            ebmerged.bandFit(which=self.index) # a new fit
        else:
            BandFitExtended(self.index,ebmerged,self.roi).fit()

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

        if self.manager == self.roi.psm:
            ebmerged.bandFit(which=self.which) # a new fit
        else:
            BandFitExtended(self.index,ebmerged,self.roi).fit()

        ebmerged.merged = True
        ebmerged.num = num
        return ebmerged
        
    
    def __str__(self):
        n  = len(self.rec.dtype)-2
        return ((2*'%10s'+' '+n*'%-12s'+'\n') % self.rec.dtype.names)\
             +'\n'.join( [(2*'%10.0f'+n*'%12.3e') % tuple(row) for row in self.rec])
        


class PlotSED(object):
    """
   
    
    """
    def __init__(self, fluxes, use_ergs=True):
        """ fluxes: a recarray prepared by BandFlux
        
         """
        
        self.rec=fluxes
        self.use_ergs = use_ergs
        self.scale_factor = 1.602e-6 if use_ergs else 1.0


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
            m: the model, implements Models.Model
            dom: the domain, a set of points
            butterfly: bool, 
            kwargs: pass to the line to plot
        """
        stat = m.statistical()
        err = stat[0]*stat[1]
        energy_flux_factor = self.scale_factor
        
        # show position of e0, possibly the pivot energy
        e0 = m.e0 if m.name!='LogParabola' else 10**m.p[-1]
        flux = m(e0); flux_unc = flux*stat[1][0]
        axes.errorbar([e0], 
                [energy_flux_factor*flux * e0**2], fmt='or', 
            yerr=energy_flux_factor*flux_unc * e0**2, elinewidth=2, markersize=8)

        #if m.name=='PowerLaw':
        #    axes.errorbar([m.e0], [energy_flux_factor*stat[0][0]*m.e0**2], fmt='or', 
        #        yerr=energy_flux_factor*err[0]*m.e0**2, elinewidth=2, markersize=8)
        ## plot the curve
        axes.plot( dom, energy_flux_factor*m(dom)*dom**2, **kwargs)
        #butterfly if powerlaw
        if butterfly and m.name=='PowerLaw':
            # 'butterfly' region
            dom_r = np.array([dom[-i-1] for i in range(len(dom))]) #crude reversal.
            a,gamma = stat[0][:2]
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
               
    def __call__(self, model, name,
                fignum=5, axes=None,
                axis=None, #(1e2,1e6,1e-8,1e-2),
                data_kwargs=dict(linewidth=2, color='k',),
                fit_kwargs =dict(lw=2,        color='r',),
                butterfly = True,
                use_ergs = True,
                outdir = None,
                ):
        """Plot the SED
        ========     ===================================================
        keyword      description
        ========     ===================================================
        model        spectral model object
        name         name of the source
        fignum       [5] if set, use (and clear) this figure. If None, use current Axes object
        axes         [None] If set use this Axes object
        axis         None, (1e2, 1e5, 1e-8, 1e-2) depending on use_ergs
        data_kwargs  a dict to pass to the data part of the display
        fit_kwargs   a dict to pass to the fit part of the display
        butterfly    [True] plot model with a butterfly outline
        use_ergs     [True] convert to ergs in the flux units and use GeV on the x-axis
        outdir       [None] if set, save sed into <outdir>/<source_name>_sed.png if outdir is a directory, save into filename=<outdir> if not.
        ========     ===================================================
        
        """
        energy_flux_unit = 'ergs' if self.use_ergs else 'MeV'
        energy_flux_factor = self.scale_factor
        # conversion 1.602E-19 * 1E6 eV/Mev * 1E7 erg/J * = 1.602E-6 erg/MeV
        oldlw = plt.rcParams['axes.linewidth']
        plt.rcParams['axes.linewidth'] = 2
        if axes is None: 
            fig=plt.figure(fignum, figsize=(4,4)); plt.clf()
            fig.add_axes((0.22,0.15,0.75,0.72))
            axes = plt.gca()
        axes.set_xscale('log')
        axes.set_yscale('log')
        if axis is None:
            axis = (1e2,1e6,1e-13,1e-8) if use_ergs else (1e2,1e6,1e-8,1e-2)
        axes.axis(axis)
        axes.grid(True)
        axes.set_autoscale_on(False)
       
        self.plot_data(axes, **data_kwargs)
        # and the model, perhaps with a butterfly
        #dom = np.logspace(np.log10(roi.fit_emin[0]), np.log10(roi.fit_emax[0]), 101)
        dom = np.logspace(np.log10(self.rec.elow[0]), np.log10(self.rec.ehigh[-1]), 101)
        self.plot_model(axes, model, dom, butterfly, **fit_kwargs)
        plt.rcParams['axes.linewidth'] = oldlw

        # the axis labels
        plt.ylabel(r'$\mathsf{Energy\ Flux\ (%s\ cm^{-2}\ s^{-1})}$' % energy_flux_unit)
        if use_ergs:
            plt.xlabel(r'$\mathsf{Energy\ (GeV)}$')
            axes.set_xticklabels(['','1','10','100', ''])
        else:
            plt.xlabel(r'$\mathsf{Energy\ (MeV)}$')
        plt.title(name)
        
        
        if outdir is not None: 
            if os.path.isdir(outdir):
                fname = name.replace(' ','_').replace('+','p')
                plt.savefig(os.path.join(outdir,'%s_sed.png'%fname))
            else :
                plt.savefig(outdir)
