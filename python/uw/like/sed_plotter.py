"""
Manage plotting of the band energy flux and model

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/sed_plotter.py,v 1.25 2012/09/12 21:49:56 lande Exp $

author: Toby Burnett <tburnett@uw.edu>

"""
import os,sys
from uw.like import roi_bands
from uw.like.roi_extended import BandFitExtended
from uw.like.Models import PowerLaw, LogParabola
from uw.utilities import makerec, image
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
        self.manager,self.index=roi.mapper(which)
        self.roi   = roi

        roi.setup_energy_bands()
        self.bands = roi.energy_bands
        self.scale_factor = scale_factor 
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
        
    def plot_data(self, axes, **kwargs):
        
        if 'color' not in kwargs:
            kwargs['color'] = 'k'
        printout = kwargs.pop('printout',False)
        
        def lsp(val,spaces=6):
            s = str(val)
            return ' '*(spaces-len(s))+s
        
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
            if printout:
                print '%s %s   %s   %s   %s'%(lsp(int(round(xl))),lsp(int(round(xh))),lsp('%.3g'%r.flux,8),lsp('%.3g'%r.lflux,8),lsp('%.3g'%r.uflux,8))
 
                      
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

        if hasattr(m,'e0'):
            # show position of e0 
            e0 = m.e0
            flux = m(e0); flux_unc = flux*stat[1][0]
            axes.errorbar([e0], 
                          [energy_flux_factor*flux * m.e0**2], fmt='or', 
                          yerr=energy_flux_factor*flux_unc * m.e0**2, elinewidth=2, markersize=8)

        # plot the curve
        axes.plot( dom, energy_flux_factor*m(dom)*dom**2, **kwargs)
        #butterfly if powerlaw
        if butterfly and isinstance(m,PowerLaw) or isinstance(m,LogParabola) and m[2]<=1e-3:
            # 'butterfly' region
            e0 = m.e0 if isinstance(m,PowerLaw) else m['e_break']
            dom_r = np.array([dom[-i-1] for i in range(len(dom))]) #crude reversal.
            a,gamma = stat[0][:2]
            var = err**2
            # r is e/e0
            bfun = lambda r: r**-gamma * np.sqrt(var[0] + (a*np.log(r))**2 * var[1])
            upper = energy_flux_factor*(m(dom)  + bfun(dom/e0)  )*dom**2
            lower = energy_flux_factor*(m(dom_r)/(1 +bfun(dom_r/e0)/m(dom_r)))*dom_r**2
            ymin, ymax = plt.gca().get_ylim()
            lower[lower<ymin] = ymin
            upper[upper>ymax] = ymax
            t =axes.fill(np.hstack( [dom,   dom_r] ), 
                        np.hstack( [upper, lower] ), 'r')
            t[0].set_alpha(0.4)
               
     
def plot_sed(roi, which=0, fignum=5, axes=None,
            axis=None, #(1e2,1e6,1e-8,1e-2),
            data_kwargs=dict(linewidth=2, color='k',),
            fit_kwargs =dict(lw=2,        color='r',),
            butterfly = True,
            use_ergs = False,
            energy_flux_unit = None,
            gev_scale = True,
            outdir = None,
            galmap = True,
            phase_corr=False,
            printout=False,
            title=None,
            merge=True,
            ):
    """Plot a SED
    ========     ===================================================
    keyword      description
    ========     ===================================================
    roi          a ROIAnalysis object
    which        [0] index of source to plot
    fignum       [5] if set, use (and clear) this figure. If None, 
                 use current Axes object
    axes         [None] If set use this Axes object
    axis         None, (80, 5e5, 1e-7, 1e-2) depending on use_ergs
    data_kwargs  a dict to pass to the data part of the display
    fit_kwargs   a dict to pass to the fit part of the display
    butterfly    [True] plot model with a butterfly outline
    use_ergs     [True] convert to ergs in the flux units (instead of MeV)
    energy_flux_unit    [None] If specified, one of 'erg', 'MeV', 'eV' or 'GeV': otherwise
                  set to 'erg' or 'MeV' based on use_ergs
    gev_scale    [True] use GeV instead of MeV units on x-axis
    outdir       [None] if set, save sed into 
                 <outdir>/<source_name>_sed.png if outdir is a 
                 directory, save into filename=<outdir> if not.
    galmap       [True] plot position on galactic map if set
    phase_corr   [False] multiply sed by phase_factor; appropriate 
                 for an on-pulse spectral analysis
    printout     [False] if True, print the sed points to stdout
    title        [None] Title for the plot, if specified. Otherwise, 
                 use source name
    merge        merge upper limits on edge.
    ========     ===================================================
    
    """
    self = roi # temp.
    if energy_flux_unit is None:
        energy_flux_unit = 'erg' if use_ergs else 'MeV'
    assert energy_flux_unit in ('erg', 'MeV', 'GeV', 'eV') , 'unrecognized energy flux unit'
    energy_flux_factor = dict(erg=1.602e-6, MeV=1, eV=1e6, GeV=1e-3)[energy_flux_unit]
    if phase_corr: energy_flux_factor *=roi.phase_factor
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
        axis = (80, 5e5, 1e-7*energy_flux_factor,1e-2*energy_flux_factor) 
    axes.axis(axis)
    axes.grid(True)
    axes.set_autoscale_on(False)
   
    #  create a BandFlux, and have it plot the band fluxes, merging adjacent limits at the ends
    bf = BandFlux(self, which=which, merge=merge, scale_factor= energy_flux_factor)
    data_kwargs['printout'] = printout
    bf.plot_data(axes, **data_kwargs)

    source = roi.get_source(which)
    model = source.model
    name  = source.name
    
    # and the model, perhaps with a butterfly
    dom = np.logspace(np.log10(roi.fit_emin[0]), np.log10(roi.fit_emax[0]), 101)
    bf.plot_model(axes, model, dom, butterfly, **fit_kwargs)
    plt.rcParams['axes.linewidth'] = oldlw

    # the axis labels
    plt.ylabel(r'$\mathsf{Energy\ Flux\ (%s\ cm^{-2}\ s^{-1})}$' % energy_flux_unit)
    def gevticklabel(x):
        if x<100 or x>1e5: return ''
        elif x==100: return '0.1'
        return '%d'% (x/1e3)
    if gev_scale:
        """ make it look nicer """
        axes.set_xticklabels(map(gevticklabel, axes.get_xticks()))
        axes.set_xlabel(r'$\mathsf{Energy\ (GeV)}$')
    else:
        axes.set_xlabel(r'$\mathsf{Energy\ (MeV)}$')

    plt.title(name if title is None else title)
    
    # a galactic map if requested
    if galmap: image.galactic_map(roi.roi_dir, color='lightblue', marker='o', markercolor='r')
    
    if outdir is not None: 
        if os.path.isdir(outdir):
            fname = name.replace(' ','_').replace('+','p')
            plt.savefig(os.path.join(outdir,'%s_sed.png'%fname))
        else :
            plt.savefig(outdir)

    return bf

