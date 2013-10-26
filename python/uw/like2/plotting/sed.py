"""
Manage a SED plot
    
    To generate a plot from a SourceFlux, given:
            sf an SourceFlux object, 
        Plot(sf)()

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/plotting/sed.py,v 1.16 2013/05/20 14:06:26 burnett Exp $
"""
import os, types
import numpy as np
import pylab as plt
from uw.utilities import image

def set_xlabels( axes, gev_scale) :   
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

     
class Plot(object):
    """
    
    """
    def __init__(self, source, energy_flux_unit='eV', gev_scale=True):
        """ source
        
        """
        self.name = source.name
        self.model = source.spectral_model
        self.rec = source.sedrec
        assert energy_flux_unit in ('erg', 'MeV', 'GeV', 'eV') , 'unrecognized energy flux unit'
        self.energy_flux_unit= energy_flux_unit
        self.scale_factor = dict(erg=1.602e-12, MeV=1e-6, eV=1., GeV=1e-9)[energy_flux_unit]
        self.gev_scale=gev_scale

    def plot_data(self, axes, **kwargs):
        ul_kwargs = kwargs.copy()
        ul_kwargs['color']='gray'
        if 'color' not in kwargs:
            kwargs['color'] = 'k'

        for r in self.rec:
            xl, xh = r.elow, r.ehigh
            fac = self.scale_factor 
            bc = (xl*xh)**0.5
            if r.flux >0:
                axes.plot([xl,xh], [r.flux*fac, r.flux*fac],  **kwargs)
                axes.plot([bc,bc], [r.lflux*fac,r.uflux*fac], **kwargs)
            else:
                x,y = bc, r.uflux*fac
                axes.plot([xl,xh], [y,y] , **ul_kwargs) # bar at upper limit
                # plot arrow 0.6 long by 0.4 wide, triangular head (in log coords)
                axes.plot([x, x,     x*1.2, x,     x/1.2, x],
                          [y, y*0.6, y*0.6, y*0.4, y*0.6, y*0.6], **ul_kwargs)
 
                      
    def plot_model(self,  m,  butterfly=False, **kwargs):
        """ 
            m: the model, implements Models.Model
            dom: the domain, a set of points
            butterfly: bool, 
            kwargs: pass to the line to plot
        """
        energy_flux_factor = self.scale_factor*1e6 # from MeV to eV
        dom = self.dom
        axes = self.axes
        
        ## plot the curve
        axes.plot( dom, energy_flux_factor*m(dom)*dom**2, **kwargs)
 
        # show position of e0, possibly the pivot energy
        if butterfly:
            try:
                self.plot_butterfly(m)
            except:
                print 'fail to plot butterfly for {}'.format(self.name)

    def plot_butterfly(self, m, ):
        energy_flux_factor = self.scale_factor*1e6 # from MeV to eV
        axes = self.axes
        dom = self.dom
        try:
            e0 = m.pivot_energy()
        except:
            print 'no butterfly plot for {}'.format(self.name)
            return
        flux = m(e0); 
        eflux = lambda e: energy_flux_factor * m(e) * e**2
        bfun  = lambda e: m.flux_relunc(e)

        axes.errorbar([e0], [eflux(e0)], yerr=[eflux(e0)*bfun(e0)], 
                    fmt='+r', elinewidth=2, markersize=8)
                
        dom_r = np.array([dom[-i-1] for i in range(len(dom))]) #crude reversal.
        upper = eflux(dom)  * (1 + bfun(dom)  ) 
        lower = eflux(dom_r) /(1 + bfun(dom_r) )
        ymin, ymax = axes.get_ylim()
        lower[lower<ymin] = ymin
        upper[upper>ymax] = ymax
        t =axes.fill(np.hstack( [dom,   dom_r] ), 
                    np.hstack( [upper, lower] ), 'r')
        t[0].set_alpha(0.4)
        
    def plot_residual(self, axes, model, dom, **kwargs):
        energy_flux_factor = self.scale_factor*1e6 # from MeV to eV
        energy = sqrt(self.rec.elow*self.rec.ehigh)
        mflux = model(energy)*energy**2*energy_flux_factor
        y = self.rec.flux/mflux
        yerr = self.rec.uflux/mflux-y, y-self.rec.lflux/mflux
        axes.errorbar(energy, y, yerr=yerr, fmt='o', label='UW data')
        setp(axes, xscale='log', xlim=(100, 10e3), ylim =(0.85, 1.15), 
                xlabel='Energy (Mev)')
        axes.grid()
        axes.axhline(1.0, color='k', lw=2)

    def __call__(self, model=None, name=None,
                fignum=5, axes=None,
                axis=None, #(1e2,1e6,1e-7,1e-2),
                data_kwargs=dict(linewidth=2, color='k',),
                fit_kwargs =dict(lw=2,        color='r',),
                butterfly = True,
                outdir = None,
                suffix = '_sed',
                galmap=None,
                annotate=None,
                grid=True,
                ):
        """Plot the SED
        ========     ===================================================
        keyword      description
        ========     ===================================================
        model        spectral model object
        name         name of the source
        fignum       [5] if set, use (and clear) this figure. If None, use current Axes object
        axes         [None] If set use this Axes object
        axis         None, (1e2, 1e5, 1e-8, 1e-2) depending on energy flux unit
        data_kwargs  a dict to pass to the data part of the display
        fit_kwargs   a dict to pass to the fit part of the display
        butterfly    [True] plot model with a butterfly outline
        outdir       [None] if set, save sed into <outdir>/<source_name>_sed.png if outdir is a directory, save into filename=<outdir> if noself.
        suffix       ['_sed'] Add to source name to form filename
        galmap       [None] if set to a SkyDir, create a little galactic map showing this position
        annotate     [None] if set, a tuple of (x, y, text), in axes coords
        grid         [True] Set False to turn off grid
        ========     ===================================================
        
        """
        if model is None: model=self.model
        if name is None: name = self.name
        energy_flux_factor = self.scale_factor
        # conversion 1.602E-19 * 1E6 eV/Mev * 1E7 erg/J * = 1.602E-6 erg/MeV
        oldlw = plt.rcParams['axes.linewidth']
        oldticksize = plt.rcParams['xtick.labelsize']
        plt.rcParams['axes.linewidth'] = 2
        plt.rcParams['xtick.labelsize']=plt.rcParams['ytick.labelsize']=10
        if axes is None: 
            fig=plt.figure(fignum, figsize=(4,4)); plt.clf()
            fig.add_axes((0.22,0.15,0.75,0.72))
            axes = plt.gca()
        self.axes = axes
        axes.set_xscale('log')
        axes.set_yscale('log')
        if axis is None:
            axis = (1e2,4e5, 0.2*self.scale_factor, 1e4*self.scale_factor) 
        axes.axis(axis)
        axes.grid(grid)
        axes.set_autoscale_on(False)
       
        self.plot_data(axes, **data_kwargs)
        # and the model, perhaps with a butterfly
        self.dom = np.logspace(np.log10(self.rec.elow[0]), np.log10(self.rec.ehigh[-1]), 26)
        self.plot_model( model, butterfly, **fit_kwargs)
        plt.rcParams['axes.linewidth'] = oldlw
        plt.rcParams['xtick.labelsize']=plt.rcParams['ytick.labelsize']=oldticksize

        # the axis labels (note reduced labelpad for y) 
        axes.set_ylabel(r'$\mathsf{Energy\ Flux\ (%s\ cm^{-2}\ s^{-1})}$' % self.energy_flux_unit, labelpad=0)
        axes.set_xlabel(r'$\mathsf{Energy\ (GeV)}$')
        if self.energy_flux_unit=='eV':
           axes.set_yticklabels(['', '1', '10', '100', r'$\mathdefault{10^{3}}$'])

        axes.set_title(name, size=14)
        set_xlabels(axes, self.gev_scale)
        # add a galactic map if requested
        if galmap is not None:
            image.galactic_map(galmap, axes=self.axes, color='lightblue', marker='s', markercolor='r', markersize=20)

        if annotate is not None:
            axes.text(annotate[0],annotate[1], annotate[2],transform=axes.transAxes, fontsize=8)
        if outdir is not None: 
            self.name=name
            self.savefig( outdir, suffix)
            

    def savefig(self,outdir, suffix=''):
        if os.path.isdir(outdir):
            fname = self.name.replace(' ','_').replace('+','p') + suffix
            outf = os.path.join(outdir,'%s.png'% fname)
            plt.savefig(outf)
            print 'saved sedfig to %s' %outf
        else :
            plt.savefig(outdir)
            print 'saved sedfig to %s' %outdir

def sed_table(roi, source_name=None):
    """
    """
    import pandas as pd
    si = roi.get_sed(source_name)
    return pd.DataFrame(dict(flux=si.flux.round(1), TS=si.ts.round(1), lflux=si.lflux.round(1),
            uflux=si.uflux.round(1), mflux=si.mflux.round(1), pull=si.pull.round(1)), 
                index=np.array(np.sqrt(si.elow*si.ehigh),int), columns='flux lflux uflux mflux TS pull'.split())
                

def stacked_plots(roi, source_name=None, outdir=None, fignum=6, **kwargs):
    """ 
    Make stacked plots
        
        roi : A ROIstat object
            Uses the name as a title unless title specified
        outdir : None or the name of a folder
            In the folder case, makes a file name from the ROI name
            
        Creates the two Axes objects, and returns them
    """
    plt.close(fignum)
    oldlw = plt.rcParams['axes.linewidth']
    plt.rcParams['axes.linewidth'] = 2
    fig, axes = plt.subplots(2,1, sharex=True, num=fignum, figsize=(4,5),dpi=100)
    fig.subplots_adjust(hspace=0)
    axes[0].tick_params(labelbottom='off')
    left, bottom, width, height = (0.15, 0.10, 0.75, 0.85)
    fraction = 0.8

    axes[0].set_position([left, bottom+(1-fraction)*height, width, fraction*height])
    axes[1].set_position([left, bottom, width, (1-fraction)*height])
    source = roi.get_source(source_name)
    if not hasattr(source, 'sedrec'):
        roi.get_sed(source_name)
    kw = dict(energy_flux_unit=kwargs.pop('energy_flux_unit', 'eV'))
    p = Plot(source, **kw)
    kwargs.update(outdir=None)
    suffix = kwargs.pop('suffix', '_sed')
    p(axes=axes[0],  **kwargs)
    axes[0].set_xlabel('') 

    energy = np.sqrt(p.rec.elow*p.rec.ehigh)

    pull = source.sedrec.pull
    axes[1].plot(energy, pull, 'ko')
    nhigh = sum(pull>3)
    if nhigh>0:  axes[1].plot(energy[pull>3], [3]*nhigh, '^r', markersize=10) 
    nlow = sum(pull<-3)
    if nlow >0:  axes[1].plot(energy[pull<-3], [-3]*nlow, 'vr', markersize=10) 
    axes[1].axhline(0, color='k')
    axes[1].grid()
    
    plt.rcParams['axes.linewidth'] = oldlw
    plt.setp(axes[1], xscale='log', ylabel='pull', ylim=(-3.5,3.5) )
    set_xlabels( axes[1], p.gev_scale )
    
    if outdir is not None:
        p.savefig( outdir, suffix)
    return axes

