"""
Plotting routines to display results of an ROI analysis.
Given an ROIAnalysis object roi:

     # Make an SED (somewhat duplicate functionality to sed_plotter.plot_sed)
     plot_spectra(roi)

     # Make a 2 dimensional map of the region
     ROIDisplay(roi).show()

     # Make a counts SED
     plot_counts(roi)

     # Plot the counts in a vertical slice
     ROISlice(roi).show()

     # Plot the integral counts in a radius.
     ROIRadialIntegral(roi).show()


$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_plotting.py,v 1.53 2011/06/22 00:31:48 lande Exp $

author: Matthew Kerr, Joshua Lande
"""
import exceptions
import pprint
import numpy as N
import copy
from . roi_bands import ROIEnergyBand
from . roi_image import ModelImage,CountsImage,RadialCounts,RadialModel,SmoothedCounts,SmoothedModel,SmoothedResidual
from . roi_extended import ExtendedSource
from . pointspec_helpers import PointSource
from . SpatialModels import SpatialMap, PseudoSpatialModel, RadiallySymmetricModel, \
        EllipticalSpatialModel, Disk, EllipticalDisk, Ring, EllipticalRing
from uw.utilities import colormaps
from uw.utilities.image import ZEA
from uw.utilities import region_writer 
from uw.utilities import keyword_options
from skymaps import SkyDir

from collections import deque

from scipy.stats import poisson,norm
from scipy.optimize import fmin,fsolve

import pylab as P
from matplotlib import rcParams,mpl,pyplot,ticker,font_manager,spines
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import FancyArrow,Circle,Ellipse

def band_spectra(r,source=0):
    
    en = N.asarray([band.e for band in r.bands])

    groupings = [deque() for x in xrange(len(r.bin_centers))]

    #group slw by energy
    for i,ei in enumerate(r.bin_centers):
        for band in r.bands:
            if band.e == ei:
                groupings[i].append(band)
        groupings[i] = list(groupings[i])

    scales     = N.zeros(len(groupings))
    scales_hi = N.zeros(len(groupings))
    scales_lo = N.zeros(len(groupings))

    m        = r.psm.point_sources[source].model
    ps      = N.asarray([sum([band.ps_counts[source] for band in g]) for g in groupings])
    exp     = N.asarray([sum([band.expected(m)/m.i_flux(band.emin,band.emax) for band in g]) for g in groupings])

    for i,gi in enumerate(groupings):
        obs = sum([band.photons for band in gi])
        """ #way to handle 0 count bins
        
        if obs == 0:
            bg = sum([slw.gal_exp*slw.gal_counts + slw.iso_exp*slw.iso_counts for slw in gi])
            so = sum([slw.ps_counts[source]*slw.overlaps[source] for slw in gi])
            eta = (1.92 - bg)*so
            print 'multi = %.2f'%(eta)
        """
        #define a joint likelihood for the energy band
        f = lambda scale: sum( [band.bandLikelihood(scale,source) for band in gi] )
        result = fmin(f,[1.],disp=0,full_output=1)
        if result[4] != 1 and result[4] != 2:
            scales[i] = result[0]
            f68 = lambda scale: -result[1] + f(scale) - 0.5 #68% confidence
            f95 = lambda scale: -result[1] + f(scale) - 1.92 #95% confidence
            
            if scales[i] < 0.05 or scales[i]*ps[i] < 0.1: #upper limit
                scales[i] = 0
                scales_lo[i] = 0
                scales_hi[i] = fsolve(f95,[max(obs/ps[i],scales[i])])

            else:
                scales_lo[i] = fsolve(f68,[0])
                scales_hi[i] = fsolve(f68,[scales[i]*1.1])

    cts     = scales*ps
    cts_hi = scales_hi*ps
    cts_lo = scales_lo*ps
    

    return r.bin_edges,cts,cts_lo,cts_hi,exp,r.psm.point_sources[source]

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def band_fluxes(r,which=0,axes=None,axis=None,outfile=None,emin=[0,0], **kwargs):

    plot_kwargs = {'color':'red', 'linewidth':1,}
    plot_kwargs.update(kwargs)
    if axes is None:
        ax = P.gca()
    else: ax = axes

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)

    r.setup_energy_bands(emin=emin)

    for eb in r.energy_bands:
        eb.bandFit(which=which)
        eb.merged = False

    # count in from high end to find how many high energy bins have only upper limits
    for n,neb in enumerate(r.energy_bands[::-1]):
        if neb.flux is not None: break
    # count in from low end to find how many low energy bins have only upper limits
    for nlo,neb in enumerate(r.energy_bands):
        if neb.flux is not None: break
    #reg_slice     = slice(0,len(r.bin_centers) - n,1)
    reg_slice     = slice(nlo,len(r.bin_centers) - n,1)
    merged_slice = slice(len(r.bin_centers)-n, len(r.bin_centers), 1)
    lo_merged_slice = slice(0,nlo,1)

    if n > 0:
        bands = []
        for eb in r.energy_bands[merged_slice]:
            for band in eb.bands:
                bands.append(band)
        ebhi = ROIEnergyBand(bands,r.energy_bands[merged_slice][0].emin,r.energy_bands[merged_slice][-1].emax)
        ebhi.bandFit(which=which)
        ebhi.merged = True
        ebhi.num = n
        ebhi = [ebhi]
    else:
        ebhi = []

    if nlo > 0:
        bands = []
        for eb in r.energy_bands[lo_merged_slice]:
            for band in eb.bands:
                bands.append(band)
        eblo = ROIEnergyBand(bands,r.energy_bands[lo_merged_slice][0].emin,r.energy_bands[lo_merged_slice][-1].emax)
        eblo.bandFit(which=which)
        eblo.merged = True
        eblo.num = nlo
        eblo = [eblo]
    else:
        eblo = []

    r.energy_bands = eblo + r.energy_bands[reg_slice] + ebhi

    for eb in r.energy_bands:

        xlo,xhi = eb.emin,eb.emax
        bc  = (xlo*xhi)**0.5
        fac = xlo*xhi 
        hw  = (xhi-xlo)/2.
        if eb.merged: hw /= eb.num
        
        if eb.flux is not None: 
            ax.plot([xlo,xhi],[fac*eb.flux,fac*eb.flux], **plot_kwargs)
            ax.plot([bc,bc]  ,[fac*eb.lflux,fac*eb.uflux],**plot_kwargs)
        else:
            ax.plot([xlo,xhi],[fac*eb.uflux,fac*eb.uflux],color='k')
            dy = -(fac*eb.uflux - 10**(N.log10(fac*eb.uflux) - 0.6))
            hl = 10**(N.log10(fac*eb.uflux) - 0.4) - 10**(N.log10(fac*eb.uflux) - 0.6)
            a = FancyArrow(bc,fac*eb.uflux,0,dy,width=(xhi-xlo)/100.,head_width=hw,
                                                      head_length=hl, length_includes_head=True,
                                                      facecolor='k',edgecolor='k',fill=False)
            ax.add_patch(a)

    if axis is not None: ax.axis(axis)
    if outfile is not None: P.savefig(outfile)


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def make_sed(r,which=0,axes=None,axis=None,plot_model=True, 
          data_kwargs=None, fit_kwargs=None,emin=[0,0]):
     
     if data_kwargs is None: data_kwargs={}
     if fit_kwargs is None: fit_kwargs={'color':'blue', 'linewidth':1}
     if axes is None: axes = P.gca()
     band_fluxes(r,which=which,axes=axes,emin=emin, **data_kwargs)
     axes.set_xlabel('Energy (MeV)')
     axes.set_ylabel('Energy Flux (MeV/cm2/s)')
     if plot_model:
          try:
                dom = N.logspace(N.log10(r.fit_emin[0]),N.log10(r.fit_emax[0]),51)
                cod = r.psm.models[which](dom)*dom**2
                axes.plot(dom,cod, **fit_kwargs)
          except exceptions.OverflowError:
                pass # failed, at least plot axes
          except:
                raise # anything else
     if axis is None:
          axes.axis([1e2,1e5,1e-10,1e-2])
     else: axes.axis(axis)
     axes.set_autoscale_on(False)
     axes.grid(True)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
    

def counts(r,integral=False):

    groupings = [deque() for x in xrange(len(r.bin_centers))]

    #group slw by energy
    for i,ei in enumerate(r.bin_centers):
        for band in r.bands:
            if band.e == ei:
                groupings[i].append(band)
        groupings[i] = list(groupings[i])

    #iso = N.asarray([ sum((band.bg_counts[1] for band in g)) for g in groupings]) * p
    #gal = N.asarray([ sum((band.bg_counts[0] for band in g)) for g in groupings]) * p
    dif = N.asarray([ N.asarray([band.phase_factor*band.bg_counts for band in g]).sum(axis=0) for g in groupings])
    obs = N.asarray([ sum((band.photons for band in g)) for g in groupings])
    src = N.asarray([ N.asarray([band.phase_factor*band.ps_counts*band.overlaps for band in g]).sum(axis=0) for g in groupings])
    
    if integral:
        for i in xrange(len(iso)):
            #iso[i] = iso[i:].sum()
            #gal[i] = gal[i:].sum()
            dif[i] = dif[i:].sum(axis=0)
            obs[i] = obs[i:].sum()
            src[i] = src[i:].sum(axis=0)

    return r.bin_edges,dif,src,obs,[b.name for b in r.bgm.bgmodels],[p.name for p in r.psm.point_sources]


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


def plot_counts(r,fignum=1,outfile=None,integral=False,max_label=10,merge_non_free=True):

    colors = ['blue','green','red','orange']

    en,dif,src,obs,bg_names,ps_names = counts(r,integral=integral)
    en = (en[1:]*en[:-1])**0.5

    # optionally merge all of the "frozen" sources for a cleaner plot
    if merge_non_free:
        free_mask = N.asarray([N.any(m.free) for m in r.psm.models])
        new_src    = N.zeros([len(en),free_mask.sum()+1])
        counter = 0
        for i in xrange(len(free_mask)):
            if free_mask[i]:
                new_src[:,counter] = src[:,i]
                counter += 1
            else:
                new_src[:,-1] += src[:,i]
        #return en,dif,src,obs,bg_names,ps_names,new_src
        src        = new_src
        ps_names = [ps_names[i] for i in xrange(len(ps_names)) if free_mask[i]]
        ps_names += ['Other Point Sources']
        
    P.figure(fignum,(14,8))
    P.clf()
    P.subplot(121)
    P.gca().set_xscale('log')
    P.gca().set_yscale('log')
    for i,name in enumerate(ps_names):
        label = name if (i < max_label) else '_nolegend_'
        P.loglog(en,src[:,i],linestyle='-',marker='',label=label)
    for i,name in enumerate(bg_names):
        if N.any(dif[:,i]==0): continue
        P.loglog(en,dif[:,i],linestyle='-',marker='',label=name)

    #tot = src.sum(axis=1)+iso+gal
    tot = src.sum(axis=1) + dif.sum(axis=1)
    P.loglog(en,tot,linestyle='steps-mid',color='black',label='Total Model')
    err = obs**0.5
    low_err = N.where(obs-err <= 0, 0.99*obs, err)
    P.errorbar(en,obs,yerr=[low_err,err],linestyle=' ',marker='x',label='Counts',color='black')
    ax = P.axis()
    P.axis([ax[0],ax[1],max(0.1,ax[2]),ax[3]])
    if integral:
        P.ylabel('Counts > E')
    else:
        P.ylabel('Counts per Bin')
    P.xlabel('Energy (MeV)')
    prop = font_manager.FontProperties(size='x-small')
    P.legend(loc=0,prop=prop)
    P.grid(b=True)

    P.subplot(122)
    P.gca().set_xscale('log')
    P.errorbar(en,(tot-obs)/(tot),yerr=tot**-0.5,linestyle=' ',marker='o',color='black')
    P.plot(N.linspace(P.axis()[0],P.axis()[1],100),[0]*100,color='black')
    ax = P.axis()
    ybound = min( 0.5, max( abs(ax[2]), abs(ax[3]) ))
    P.axis([ax[0],ax[1],-ybound,ybound])
    P.ylabel('Fractional Deviation (model-obs)/model')
    P.xlabel('Energy (MeV)')
    P.grid(b=True)

    if outfile is not None: P.savefig(outfile)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#



def plot_spectra(r, which=0, eweight=2,fignum=1,outfile=None,merge_bins=False,
                      return_vals=False,axis=None,use_legend=True,axes=None):

    colors = ['blue','green','red','orange']

    en,cts,ctslo,ctshi,exp,ps = band_spectra(r,which)

    
    if merge_bins:
        new_en     = []
        new_cts    = []
        new_ctslo = []
        new_ctshi = []
        new_exp    = []
        for i in range(len(cts))[::2]:
            new_en     += [en[i]]            
            new_cts    += [cts[i] + (cts[i+1] if i < (len(cts) - 1) else 0)]
            new_ctslo += [((cts[i]-ctslo[i])**2 + ( (cts[i+1]-ctslo[i+1])**2 if i < (len(ctslo) - 1) else 0))**0.5]
            new_ctshi += [((ctshi[i]-cts[i])**2 + ( (ctshi[i+1]-cts[i+1])**2 if i < (len(ctshi) - 1) else 0))**0.5]
            new_exp    += [(exp[i]*exp[i+(1 if i < (len(cts) - 1) else 0)])**0.5]
        new_en += [en[-1]]

        en     = N.asarray(new_en)
        cts    = N.asarray(new_cts)
        ctslo = cts - N.asarray(new_ctslo)
        ctshi = N.asarray(new_ctshi) + cts
        exp    = N.asarray(new_exp)
    
    bc = (en[1:]*en[:-1])**0.5

    if axes is None:
        P.figure(fignum)
        P.clf()
        axes = P.gca()

    fluxes     = bc**eweight*cts/exp/(en[1:]-en[:-1])
    fluxes_lo = bc**eweight*ctslo/exp/(en[1:]-en[:-1])
    fluxes_hi = bc**eweight*ctshi/exp/(en[1:]-en[:-1])
    P.gca().set_xscale('log')
    P.gca().set_yscale('log')
    if N.all(fluxes==0):
        axes.plot(bc[fluxes==0],fluxes_hi[fluxes==0],marker='v',ls=' ',mfc='black',mec='black')
    else:
        f    = fluxes[fluxes>0]
        flo = fluxes_lo[fluxes>0]
        flo[flo==0] = 1e-20#0.999*f
        fhi = fluxes_hi[fluxes>0]
        axes.errorbar(x=bc[fluxes>0],y=f,yerr=[f-flo,fhi-f],linestyle=' ',marker='o',mfc = 'white', mec = 'black',\
                              color='black',label='Band Fits',ms=10)
        axes.plot(bc[fluxes==0],fluxes_hi[fluxes==0],marker='v',ls=' ',mfc='black',mec='black')
        domain = N.logspace(N.log10(en[0]),N.log10(en[-1]),50)
        axes.plot(domain,domain**eweight*ps.model(domain),color='red',lw=2,label='Model')
        if use_legend: P.legend(loc=0,numpoints=1)
        axes.grid(b=True)
        axes.set_xlabel('Energy (MeV)')
        if axis is None:
            axes.axis([0.7*min(bc),1.3*max(bc),max(min(fluxes_lo[fluxes>0])*0.7,max(fluxes)/100.),max(fluxes_hi)*1.3])
        else:
            axes.axis(axis)

    if outfile is not None: P.savefig(outfile)
    if return_vals: return en,fluxes,fluxes_lo,fluxes_hi
        
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#



class ROIDisplay(object):
    """ Plotting a two dimensional map of the counts, model predicted,
        and residual counts in the ROI. Also plot a histogram of
        the weighted residuals and the p-values. """

    defaults = (
            ('figsize',        (7,6.5),                    'Size of the image'),
            ('fignum',          None,               'matplotlib figure number'),
            ('pixelsize',       0.25,               'size of each image pixel'),
            ('conv_type',         -1,                        'Conversion type'),
            ('size',              10,              'Size of the field of view'),
            ('nticks',             5,              'Number of axes tick marks'),
            ('label_sources',  False,               'Label sources duing plot'),
            ('galactic',        True,             'Coordinate system for plot'),
            ('countsfile',      None, 'Fits file to save the counts map data.'),
            ('modelfile',       None,  'Fits file to save the model map data.'),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        try:
            import pywcsgrid2
        except:
            raise Exception("You must install pywcsgrid2 to use this module.")
        
        self.roi = roi

        self.cm = CountsImage(self.roi,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)
        self.mm = ModelImage(self.roi,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        self.cm_pf, self.mm_pf = self.cm.get_pyfits(), self.mm.get_pyfits()

        self.cm_d, self.mm_d = self.cm_pf[0].data, self.mm_pf[0].data
        self.res_d=(self.cm_d-self.mm_d)/self.mm_d**0.5

        self.h=self.cm_pf[0].header

        if self.countsfile is not None: 
            self.cm_pf.writeto(self.countsfile,clobber=True)
        if self.modelfile is not None: 
            self.mm_pf.writeto(self.modelfile,clobber=True)

        # Use same scale for the counts and model map
        self.counts_max = N.ceil(max(N.max(self.cm.image),N.max(self.mm.image)))
        self.counts_min = max(N.floor(min(N.min(self.cm.image),N.min(self.mm.image))),1.0)

    def add_cbar(self,im,ax):
        from matplotlib.axes import Axes
        import pylab as P
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
        self.fig.add_axes(cax)
        cbar = P.colorbar(im, cax=cax)

    def model_plot(self):

        from uw.utilities import colormaps
        im=self.ax_model.imshow(self.mm_d, norm=self.counts_norm, 
                            cmap=colormaps.b,
                            **self.imshow_args)

        #self.grid[0].cax.colorbar(im)
        self.ax_model.set_title('Model Counts')
        self.ax_model.grid(linestyle='-')
        self.add_cbar(im,self.ax_model)
        
    def counts_plot(self):

        from uw.utilities import colormaps
        d=self.cm_d
        m=self.counts_min
        im=self.ax_counts.imshow(N.where(d>m,d,m), norm=self.counts_norm,
                            cmap=colormaps.b,
                            **self.imshow_args)
        self.ax_counts.set_title('Observed Counts')
        self.ax_counts.grid(linestyle='-')
        self.add_cbar(im,self.ax_counts)

    def resids_plot(self):

        im=self.ax_res.imshow(self.res_d,
                            norm=self.norm_res,
                            **self.imshow_args)
        self.ax_res.set_title('Weighted Residuals')
        self.ax_res.grid(linestyle='-')
        self.add_cbar(im,self.ax_res)

    def hist_plot(self):
        self.ax_pvals.set_title('P-Values')

        mc = self.mm_d.flatten()
        cc = self.cm_d.flatten()
        rc = self.res_d.flatten()

        pvals = 1-poisson.cdf(cc,mc)

        nb = 20
        av = float(len(pvals)) / nb
        self.ax_pvals.hist(pvals,bins=N.linspace(0,1,20),histtype='step')
        self.ax_pvals.axhline( av, color='red')  
        lo,hi = ppf( (50.-95/2)/100., av), ppf( (50. + 95/2)/100.,av)
        self.ax_pvals.axhspan( lo , hi , facecolor='red', alpha=0.3)

        from matplotlib.offsetbox import AnchoredText
        from matplotlib.font_manager import FontProperties
        from matplotlib.patheffects import withStroke

        font=dict(size='small');
        at=AnchoredText('$95\%$ Conf.', loc=1, prop=font,frameon=False)
        self.ax_pvals.add_artist(at)
        at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])

        self.ax_resplot.set_title('Weighted Residuals')

        self.ax_resplot.hist(rc, bins=N.linspace(-5,5,20), histtype='step',normed=True)
        b=N.linspace(-5,5,100)
        self.ax_resplot.plot(b,norm.pdf(b))
        self.ax_resplot.axvline(0,color='red')
        self.ax_resplot.grid(True)
        self.ax_resplot.set_xbound(lower=-5,upper=5)

    def show(self,filename=None):

        # taken from http://matplotlib.sourceforge.net/users/usetex.html
        from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
        import pywcsgrid2
        from matplotlib import mpl
        from matplotlib.gridspec import GridSpec,GridSpecFromSubplotSpec

        self.imshow_args = dict(interpolation='nearest', origin='lower')

        self.fig = P.figure(self.fignum,self.figsize)
        P.clf()

        self.counts_norm = mpl.colors.LogNorm(vmin=self.counts_min,vmax=self.counts_max)
        self.norm_res = mpl.colors.Normalize(vmin=-5,vmax=5)


        # first, divide the plot in 4
        gs = GridSpec(4, 4)
        gs.update(wspace=0.75, hspace=0.75) # wspace=0.75,
        self.ax_model = pywcsgrid2.subplot(gs[0:2, 0:2], header=self.h)
        self.ax_counts = pywcsgrid2.subplot(gs[0:2,2:4], header=self.h)
        self.ax_res = pywcsgrid2.subplot(gs[2:4, 0:2], header=self.h)

        # then divide the 4th space in two.
        self.ax_pvals = P.subplot(gs[2, 2:4])
        self.ax_resplot = P.subplot(gs[3, 2:4])

        self.model_plot()
        self.counts_plot()
        self.resids_plot()
        self.hist_plot()

        for ax in [self.ax_model, self.ax_counts, self.ax_res]:
            ROISignificance.plot_sources(self.roi,ax,self.h,label_sources=self.label_sources,
                                         show_extension=False, marker_scale=2, color='k')

        if filename is not None: P.savefig(filename)



#===============================================================================================#

def ppf(prob,mean):
    """Return the (approximate) Poisson percentage point function for given distribution.  Klugey."""
    if mean > 100: #normal approximation
        n = norm(mean,mean**0.5)
        return n.ppf(prob)        
    d = poisson(mean)
    prev = 0
    for i in xrange(1,200):        
        new = d.cdf(i)
        if new >= prob: break
        prev = new
    return (i-1)+(prob-prev)/(new-prev) #linear interpolation

#===============================================================================================#

def int2bin(n, count=24):
    """returns the binary of integer n, using count number of digits"""
    return "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])

#===============================================================================================#

class ROISlice(object):
    """ Create counts slice plot (integrating out one dimention and plotting
        the counts and model predicted counts in the other direction. """

    defaults = (
            ('which',           None,                       'Source to analyze'),
            ('figsize',        (7,6),                       'Size of the image'),
            ('fignum',          None,                'matplotlib figure number'),
            ('pixelsize',      0.125,                'size of each image pixel'),
            ('size',              10,               'Size of the field of view'),
            ('galactic',        True,              'Coordinate system for plot'),
            ('int_width',          2,            'Integration width for slice.'),
            ('conv_type',         -1,                         'Conversion type'),
            ('just_diffuse',    True, """Display the model predictions with 
                                               all point + extended sources 
                                             removed. The background is not 
                                                                     refit. """),
            ('aspoint',         True, """Display also the model predictions 
                                            for an extended source fit with 
                                           the point hypothesis. Only works 
                                          when which is an extended source. """),
            ('oversample_factor',  4, """ Calculate the model predictions 
                                          this many times more finely.  
                                          This will create a smoother 
                                          plot of model predictions. Set 
                                          to 1 if you want the model 
                                          predictions to 'look like' the 
                                          data.                             """),
            ('use_gradient',    True,            'Use gradient when refitting.'),
            ('title',           None,                      'Title for the plot'),
            ('black_and_white',False, """ If True, make the plot black and 
                                          white (better for printing and
                                          publishing)                       """),
    )

    @staticmethod
    def set_color_cycle():
        P.gca().set_color_cycle(['k', 'b', 'g', 'r', 'm', 'y', 'k'])


    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        if not type(self.oversample_factor) == int:
            raise Exception("oversample_factor must be an integer.")

        self.roi   = roi

        manager,index=self.roi.mapper(self.which)
        if manager == roi.psm:
            self.source=manager.point_sources[index]
            self.pretty_name = 'Point'
        else:
            self.source=manager.diffuse_sources[index]
            self.pretty_name = self.source.spatial_model.pretty_name

        self.center=self.source.skydir

        fig = P.figure(self.fignum,self.figsize)
        P.clf()

        self.get_counts()
        self.get_model()

    def get_counts(self):
        self.ci_x = CountsImage(self.roi,center=self.center,size=(self.size,self.int_width),
                                pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        self.counts_x=self.ci_x.image.sum(axis=0)
        self.counts_dx=ROISlice.range(self.counts_x,self.pixelsize,x_axis=True)

        if N.any(N.isnan(self.counts_x)):
            raise Exception("Error in calculating a slice, the x counts slice contains NaNs.")

        self.ci_y = CountsImage(self.roi,center=self.center,size=(self.int_width,self.size),
                                pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        self.counts_y=self.ci_y.image.sum(axis=1)
        self.counts_dy=ROISlice.range(self.counts_y,self.pixelsize,x_axis=False)

        if N.any(N.isnan(self.counts_y)):
            raise Exception("Error in calculating a slice, the y counts slice contains NaNs.")

    @staticmethod
    def cache_roi(roi):
        roi.old_parameters = roi.parameters().copy() # save free parameters

    @staticmethod
    def uncache_roi(roi):
        roi.set_parameters(roi.old_parameters) # reset free parameters
        roi.__update_state__()

    def get_model(self):
        model_pixelsize=float(self.pixelsize)/self.oversample_factor

        kwargs=dict(center=self.center,pixelsize=model_pixelsize,
                    galactic=self.galactic,conv_type=self.conv_type)

        self.names = []
        self.mi_x = []
        self.mi_y = []

        self.names.append(self.pretty_name)
        self.mi_x.append(ModelImage(self.roi,size=(self.size,self.int_width),**kwargs))
        self.mi_y.append(ModelImage(self.roi,size=(self.int_width,self.size),**kwargs))

        if self.aspoint and isinstance(self.source,ExtendedSource):
            # replace with a point source

            ROISlice.cache_roi(self.roi)

            es=self.roi.get_source(self.which)
            
            ps=PointSource(name=es.name,model=es.model.copy(),skydir=es.spatial_model.center,leave_parameters=True)
            self.roi.add_source(ps)

            # only zero it after making a copy of the spectral part!
            self.roi.zero_source(es)

            self.roi.fit(estimate_errors=False,use_gradient=self.use_gradient)

            self.names.append('Point')
            self.mi_x.append(ModelImage(self.roi,size=(self.size,self.int_width),**kwargs))
            self.mi_y.append(ModelImage(self.roi,size=(self.int_width,self.size),**kwargs))

            self.roi.del_source(ps)
            self.roi.unzero_source(es)

            ROISlice.uncache_roi(self.roi)

        if self.just_diffuse:
            # hide all point + extended sources.

            sources = list(self.roi.psm.point_sources) + \
                    [ i for i in self.roi.dsm.diffuse_sources if isinstance(i,ExtendedSource) ]
            # don't zero already zeroed sources
            sources = [ i for i in sources if i.model.getp(0,internal=True) != -100 ]

            ROISlice.cache_roi(self.roi)

            for source in sources: self.roi.zero_source(source)

            self.names.append('Diffuse')
            self.mi_x.append(ModelImage(self.roi,size=(self.size,self.int_width),**kwargs))
            self.mi_y.append(ModelImage(self.roi,size=(self.int_width,self.size),**kwargs))

            for source in sources: self.roi.unzero_source(source)

            ROISlice.uncache_roi(self.roi)

        self.models_x=[model.image.sum(axis=0)*self.oversample_factor for model in self.mi_x]
        self.models_y=[model.image.sum(axis=1)*self.oversample_factor for model in self.mi_y]

        self.model_dx = ROISlice.range(self.models_x[0],model_pixelsize,x_axis=True)
        self.model_dy = ROISlice.range(self.models_y[0],model_pixelsize,x_axis=False)

    @staticmethod
    def range(data,pixelsize,x_axis=False):
        if x_axis:
            return (len(data)/2.0-N.arange(len(data)))*pixelsize - pixelsize/2
        else:
            return (N.arange(len(data))-len(data)/2.0)*pixelsize + pixelsize/2

    @staticmethod
    def get_styles(black_and_white):
        return ['k-','k--','k.-'] if black_and_white else ['r-','b--','g-.']

    def plotx(self):

        P.subplot(211)

        ROISlice.set_color_cycle()

        styles=ROISlice.get_styles(self.black_and_white)[0:len(self.names)]
        for name,model,style in zip(self.names,self.models_x,styles):
            P.plot(self.model_dx,model,style,label=name)

        P.errorbar(self.counts_dx,self.counts_x,yerr=N.sqrt(self.counts_x),label='Counts', fmt='.')

        P.gca().set_xlim(self.counts_dx[0],self.counts_dx[-1])

        P.legend(loc='upper right',numpoints=1)

        P.gca().xaxis.set_major_formatter(FormatStrFormatter('$%g^\circ$'))

        P.xlabel(r'$\Delta l$' if self.galactic else r'$\Delta \mathrm{RA}$')
        P.ylabel('Counts')

    def ploty(self):

        P.subplot(212)

        ROISlice.set_color_cycle()

        styles=ROISlice.get_styles(self.black_and_white)[0:len(self.names)]
        for name,model,style in zip(self.names,self.models_y,styles):
            P.plot(self.model_dy,model,style,label=name)

        P.errorbar(self.counts_dy,self.counts_y,yerr=N.sqrt(self.counts_y),label='Counts', fmt='.')

        P.gca().xaxis.set_major_formatter(FormatStrFormatter('$%g^\circ$'))

        P.xlabel(r'$\Delta b$' if self.galactic else r'$\Delta \mathrm{Dec}$')
        P.ylabel('Counts')

        P.gca().set_xlim(self.counts_dy[0],self.counts_dy[-1])

        # only need to show legend once

    def save_data(self,datafile):
        """ Note, shrink model predicted counts to be same size as regular counts,
            for an easier file format. """

        x,y=['l','b'] if self.galactic else ['ra','dec']
            
        results_dict = {
            x : {
                'Counts': [self.counts_dx.tolist(), self.counts_x.tolist()]
            },
            y : {
                'Counts': [self.counts_dy.tolist(), self.counts_y.tolist()]
            }
        }

        for name,modelx,modely in zip(self.names,self.models_x,self.models_y):
            results_dict[x][name]=[self.model_dx.tolist(), modelx]
            results_dict[y][name]=[self.model_dy.tolist(), modely]

        file = open(datafile,'w')
        try:
            import yaml
            file.write(yaml.dump(results_dict))
        except:
            import pprint
            file.write(pprint.pformat(results_dict))

        file.close()

    def show(self,filename=None, datafile=None):

        self.plotx()
        self.ploty()

        if self.title is None:
            self.title = 'Counts Slice'
            if self.source is not None: self.title += ' for %s' % self.source.name

        P.suptitle(self.title)

        if datafile is not None: self.save_data(datafile)

        if filename is not None: P.savefig(filename)


class ROIRadialIntegral(object):
    """ Create a radial integral plot which integrates radially
        the counts and model predicted counts and bins uniformly in theta^2. """

    defaults = (
            ('which',           None,                       'Source to analyze'),
            ('figsize',        (7,6),                       'Size of the image'),
            ('size',               2,                'Size of image in degrees'), 
            ('pixelsize',     0.0625, """ size of each image pixel. This is a misleading because the
                                          size of each pixel varies, since the image is uniform in theta^2. 
                                          This value is used to determine the total number of pixels using
                                          the formula npix=size/pixelsize and represents something
                                          like an average pixelsize."""),
            ('npix',            None, """ If specified, use this value instead of pixelsize. """),
            ('fignum',          None, 'matplotlib figure number'),
            ('conv_type',         -1, 'Conversion type'),
            ('just_diffuse',    True, """ Display the model predictions with all point + extended
                                          sources removed. The background is not refit. """),
            ('aspoint',         True, """ Display also the model predictions for an extended source 
                                          fit with the point hypothesis. Only works when which is an
                                          extended source. """),
            ('oversample_factor',  4, """ Calculate the model predictions this many times more finely. 
                                          This will create a smoother plot of model predictions. Set 
                                          to 1 if you want the model predictions to 'look like' the 
                                          data."""),
            ('use_gradient',    True, """Use gradient when refitting."""),
            ('title',           None,   'Title for the plot'),
            ('black_and_white',False, """If True, make the plot black and white (better for 
                                         printing/publishing)"""),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        if not type(self.oversample_factor) == int:
            raise Exception("oversample_factor must be an integer.")


        self.roi   = roi

        manager,index=self.roi.mapper(self.which)
        if manager == roi.psm:
            self.source=manager.point_sources[index]
            self.pretty_name = 'Point'
        else:
            self.source=manager.diffuse_sources[index]
            self.pretty_name = self.source.spatial_model.pretty_name

        self.center=self.source.skydir

        fig = P.figure(self.fignum,self.figsize)
        P.clf()

        self.get_counts()
        self.get_model()

    def get_counts(self):
        self.ci = RadialCounts(self.roi,center=self.center,size=self.size,pixelsize=self.pixelsize,conv_type=self.conv_type, npix=self.npix)
        self.counts = self.ci.image

        self.theta_sqr_co=self.ci.bin_centers_deg

        if N.any(N.isnan(self.counts)):
            raise Exception("Error in calculating a radial integral plot, the counts contains NaNs.")


    def get_model(self):
        kwargs=dict(center=self.center,size=self.size,conv_type=self.conv_type, 
                    pixelsize=self.pixelsize/self.oversample_factor,
                    npix=self.npix*self.oversample_factor if self.npix is not None else None)

        # Use lists to preserve legend order
        self.mi,self.names = [],[]

        self.mi.append(RadialModel(self.roi,**kwargs))
        self.names.append(self.pretty_name)

        if self.aspoint and isinstance(self.source,ExtendedSource):

            ROISlice.cache_roi(self.roi)

            es=self.roi.get_source(self.which)

            ps=PointSource(name=es.name,model=es.model.copy(),skydir=es.spatial_model.center,leave_parameters=True)
            self.roi.add_source(ps)

            # only zero it after making a copy of the spectral part!
            self.roi.zero_source(es)

            self.roi.fit(estimate_errors=False,use_gradient=self.use_gradient)

            self.mi.append(RadialModel(self.roi,**kwargs))
            self.names.append('Point')

            self.roi.del_source(ps)
            self.roi.unzero_source(es)

            ROISlice.uncache_roi(self.roi)

        if self.just_diffuse:

            sources = list(self.roi.psm.point_sources) + \
                    [ i for i in self.roi.dsm.diffuse_sources if isinstance(i,ExtendedSource) ]
            # don't zero already zeroed sources
            sources = [ i for i in sources if i.model.getp(0,internal=True) != -100 ]

            ROISlice.cache_roi(self.roi)

            for source in sources: self.roi.zero_source(source)

            self.mi.append(RadialModel(self.roi,**kwargs))
            self.names.append('Diffuse')

            for source in sources: self.roi.unzero_source(source)

            ROISlice.uncache_roi(self.roi)

        self.theta_sqr_mo=self.mi[0].bin_centers_deg

        # scale the model to line up with the counts
        for i in self.mi: i.image*=self.oversample_factor
        self.models=[i.image for i in self.mi]

        for name,model in zip(self.names,self.models):
            if N.any(N.isnan(model)):
                raise Exception("Error in calculating a radial integral, model %s contains NaNs." % name)

    def save_data(self,datafile):

        results_dict = {}

        results_dict['Counts']=[ self.theta_sqr_co.tolist(), self.counts.tolist() ]
        for name,model in zip(self.names,self.models):
            results_dict[model]=[ self.theta_sqr_mo.tolist(), self.model.tolist() ]

        file = open(datafile,'w')
        try:
            import yaml
            file.write(yaml.dump(results_dict))
        except:
            import pprint
            file.write(pprint.pformat(results_dict))
        file.close()

    def show(self,filename=None,datafile=None):

        ROISlice.set_color_cycle()

        styles=ROISlice.get_styles(self.black_and_white)[0:len(self.names)]
        for name,model,style in zip(self.names,self.models,styles):
            P.plot(self.theta_sqr_mo,model,style,label=name)

        P.errorbar(self.theta_sqr_co,self.counts,yerr=N.sqrt(self.counts),label='Counts', fmt='.')

        P.legend(loc='upper right',numpoints=1)

        P.xlabel(r'$\Delta \theta^2 ([\mathrm{deg}]^2)$')
        P.ylabel('Counts')

        if self.title is None:
            self.title = 'Radially Integrated Counts'
            if self.source is not None: self.title += ' for %s' % self.source.name

        P.title(self.title)

        if datafile is not None: self.save_data(datafile)

        if filename is not None: P.savefig(filename)


class ROISignificance(object):
    """ Make a plot of the statistical significance (D-M)/sqrt(M)
        where D and M are the measured counts and the model predicted
        counts integrated within a circual aperature of radius . """

    defaults = (
            ('figsize',        (6,5),                        'Size of the image'),
            ('fignum',          None,                 'matplotlib figure number'),
            ('conv_type',         -1,                          'Conversion type'),
            ('size',               5,                'Size of the field of view'),
            ('galactic',        True,               'Coordinate system for plot'),
            ('kernel_rad',       0.25, 'Sum counts/model within radius degrees.'),
            ('label_sources',  False,                 'Label sources duing plot'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)
        
        self.roi = roi

        # Fit many pixels inside of the summing radius
        self.pixelsize=self.kernel_rad/10.0

        kwargs=dict(size=self.size,
                    pixelsize=self.pixelsize,
                    galactic=self.galactic,
                    conv_type=self.conv_type,
                    kerneltype='tophat',
                    kernel_rad=self.kernel_rad)

        self.counts=SmoothedCounts(self.roi,**kwargs)
        self.model=SmoothedModel(self.roi,**kwargs)

        self.significance = (self.counts.image - self.model.image)/N.sqrt(self.model.image)

        self.pyfits = self.counts.get_pyfits()
        self.pyfits[0].data = self.significance

    @staticmethod
    def plot_sources(roi, ax, header, color='black', show_sources=True, white_marker=True, marker_scale=4, 
            label_sources=True, show_extension=True, extension_color='white'):

        sources = roi.get_sources()

        if len(sources)<1: return

        ras = [source.skydir.ra() for source in sources]
        decs = [source.skydir.dec() for source in sources]
        
        # plot sources
        markersize=marker_scale*6
        if white_marker: ax["fk5"].plot(ras,decs,'w+',markersize=markersize)
        ax["fk5"].plot(ras,decs,'x',color=color,markersize=markersize)

        if label_sources: 
            from matplotlib.patheffects import withStroke
            names = [source.name for source in sources]
            for ra,dec,name in zip(ras,decs,names):
                myeffect = withStroke(foreground="w", linewidth=2)
                kwargs=dict(path_effects=[myeffect])
                ax["fk5"].annotate(name, (ra,dec), 
                        ha='center', va='top',
                        xytext=(0,markersize), textcoords='offset points',**kwargs)

        if show_extension:

            kwargs = dict(color=extension_color,fill=False)
            for source in roi.get_extended_sources():
                sm=source.spatial_model
                ra,dec=sm.center.ra(),sm.center.dec()

                if isinstance(sm,PseudoSpatialModel) or type(sm) == SpatialMap:
                    pass

                elif isinstance(sm,RadiallySymmetricModel):
                    if isinstance(sm,Disk) or isinstance(sm,Ring):
                        ax["fk5"].add_patch(Circle((ra,dec),sm.sigma,**kwargs))
                        if isinstance(sm,Ring):
                            ax["fk5"].add_patch(Circle((ra,dec),sm.frac*sm.sigma,**kwargs))
                    else:    
                        ax["fk5"].add_patch(Circle((ra,dec),sm.r68(),**kwargs))

                elif isinstance(sm,EllipticalSpatialModel):
                    # note ellipses in matplotlib have angle defiend from west
                    # instead of north, so we must rotate by 90 degrees. 
                    # In matplotlib, angles increase in wrong direction.
                    # Also, ellipses specify the total lenght in each
                    # direction (not lenght of semi-major/semi-minor axes),
                    # so we need to scale by a factor of 2.
                    if isinstance(sm,EllipticalDisk) or isinstance(sm,EllipticalRing):
                        sigma_x, sigma_y, theta = sm.sigma_x, sm.sigma_y, sm.theta
                        ax["fk5"].add_patch(Ellipse((ra,dec),2*sigma_x,2*sigma_y,90-theta,**kwargs))
                        if isinstance(sm,EllipticalRing):
                            frac=sm.frac
                            ax["fk5"].add_patch(Ellipse((ra,dec),2*frac*sigma_x,2*frac*sigma_y,90-theta,**kwargs))
                    else:    
                        a,b,c=sm.ellipse_68()
                        ax["fk5"].add_patch(Ellipse((ra,dec),2*a,2*b,90-c,**kwargs))
                else:
                    raise Exception("Unable to Plot Spatial Model %s" % type(sm))

    def show(self,filename=None):

        import pywcsgrid2
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        from matplotlib.axes import Axes

        self.fig = P.figure(self.fignum,self.figsize)
        P.clf()

        h, d = self.pyfits[0].header, self.pyfits[0].data

        ax = pywcsgrid2.subplot(111, header=h)

        # add colorbar axes
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
        self.fig.add_axes(cax)

        im = ax.imshow(d, origin="lower", interpolation='bilinear')
        cbar = P.colorbar(im, cax=cax)

        ax.set_title('Significance $(D-M)/\sqrt{M}$')
        
        ax.grid()

        ROISignificance.plot_sources(self.roi,ax,h,label_sources=self.label_sources)

        if filename is not None: P.savefig(filename)

class ROISmoothedSource(object):
    """ Make a smoothed residual plot to look at a particular source. 

        To do so, subtract all backgorund sources.
    
        This object requires pywcsgrid2. To get it, go to
        http://leejjoon.github.com/pywcsgrid2/users/overview.html """

    defaults = (
            ('which',             None,    'Draw the smoothed point version of this source.'),
            ('figsize',      (5.5,4.5),                                'Size of the image'),
            ('fignum',            None,                           'Matplotlib figure number'),
            ('conv_type',           -1,                                    'Conversion type'),
            ('size',                 3,                          'Size of the field of view'),
            ('galactic',          True,                         'Coordinate system for plot'),
            ('overlay_psf',       True, 'Add a smoothed reference PSF on top of the counts.'),
            ('label_psf',         True,                        'Add a label on the PSF box.'),
            ('psf_size',             1,                         'Size of the PSF insert box'), 
            ('psf_loc',              4,                       'Location to put the psf box.'), 
            ('show_sources',      True,                     'Put an x over all the sources.'),
            ('label_sources',    False,                           'Label sources duing plot'),
            ('kerneltype',  'gaussian',                'Type of kernel to smooth image with'),
            ('kernel_rad',         0.1,            'Sum counts/model within radius degrees.'),
            ('title',             None,                                 'Title for the plot'),
            ('show_extension',    True,                             'Overlay the extension.'),
            ('extension_color','white',                          'Color of extended sources'),
    )

    def get_residual(self,**kwargs):
        """ Allow the particular method for getting the residual image to be overloaded. """

        self.roi.zero_source(which=self.which)

        residual = SmoothedResidual(self.roi,**kwargs)

        self.roi.unzero_source(which=self.which)

        return residual

    @staticmethod
    def get_max_intensity(source,pyfits,roi):
        """ Return the maximum value in the pyfits file either 
            within the extended source's size or otherwise
            within the 68% containment radius of the PSF (at
            the lowest energy). """
        try:
            import pyregion

            if hasattr(source,'spatial_model'):
                # Get the maximum intensity value inside
                # the spatial model's extension
                extension_string='\n'.join(region_writer.unparse_extension(source.spatial_model,r68=True))
                reg = pyregion.parse(extension_string)
            else:
                # Get the maximum intensity inside
                emin=roi.bin_edges[0]
                ra,dec=source.skydir.ra(),source.skydir.dec()
                r68=roi.sa.psf.inverse_integral(emin,1,68)
                reg = pyregion.parse("fk5; circle(%.4f, %.4f, %.4f)" % (ra,dec,r68))

            mask = reg.get_mask(pyfits[0])

            return pyfits[0].data[mask].max()
        except:
            return pyfits[0].data.max()


    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        keyword_options.process(self, kwargs)

        self.which=which
        
        self.roi = roi

        self.cmap = colormaps.b

        # Fit many pixels inside of the summing radius
        self.pixelsize=self.kernel_rad/5.0

        self.source = roi.get_source(which)

        smoothed_kwargs=dict(size=self.size,
                    pixelsize=self.pixelsize,
                    galactic=self.galactic,
                    conv_type=self.conv_type,
                    center=self.source.skydir,
                    per_solid_angle=True,
                    kerneltype=self.kerneltype,
                    kernel_rad=self.kernel_rad)

        if self.overlay_psf:

            # convert it to a point source placed at the origin
            point_version=PointSource(name=self.source.name,
                                      skydir=self.source.skydir,
                                      model=self.source.model.copy())
        
        self.residual = self.get_residual(**smoothed_kwargs)

        self.residual_pyfits = self.residual.get_pyfits()

        if self.overlay_psf:
            # create an image of the PSF (for our model).
            # Shrink down the image of the psf
            psf_kwargs=copy.deepcopy(smoothed_kwargs)
            psf_kwargs['size']=self.psf_size
            self.psf_model=SmoothedModel(self.roi,
                    override_point_sources=[point_version],
                    **psf_kwargs)
            self.psf_pyfits = self.psf_model.get_pyfits()

    def show(self,filename=None):
        import pywcsgrid2
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.axes_grid import ImageGrid

        self.fig = P.figure(self.fignum,self.figsize)
        P.clf()

        h, d = self.residual_pyfits[0].header, self.residual_pyfits[0].data

        self.grid = grid = ImageGrid(self.fig,
            rect=(0.12,0.12,0.75,0.8),
            nrows_ncols = (1, 1),
            cbar_mode="single", cbar_pad="2%",
            cbar_location="right",
            axes_class=(pywcsgrid2.Axes, dict(header=h)))

        self.ax = ax = grid[0]

        self.max_intensity = ROISmoothedSource.get_max_intensity(self.source,self.residual_pyfits,self.roi)

        im=ax.imshow(d, origin="lower", cmap=self.cmap, vmin=0, vmax=self.max_intensity)

        cb_axes = grid.cbar_axes[0] # colorbar axes

        cbar = cb_axes.colorbar(im)
        cbar.ax.set_ylabel(r'$\mathrm{counts}/[\mathrm{deg}]^2$')

        ax.grid()

        if self.title is None: 
            self.title = 'Smoothed Counts'
            self.title += ' for %s' % self.source.name

        ax.set_title(self.title)

        if self.overlay_psf:

            # Normalize psf to have same maximum pixel scale
            # as residual image.
            self.psf_pyfits[0].data *= self.max_intensity/N.max(self.psf_pyfits[0].data)

            h_psf, d_psf = self.psf_pyfits[0].header, self.psf_pyfits[0].data
            axins = zoomed_inset_axes(ax, zoom=1, loc=self.psf_loc,
                              axes_class=pywcsgrid2.Axes,
                              axes_kwargs=dict(wcs=h_psf))

            # Note, match color maps with parent.
            axins.imshow(d_psf, cmap=im.cmap, origin="lower")
            axins.axis[:].toggle(all=False)
            axins.axis[:].line.set_color('white')

            if self.label_psf:
                axins.add_inner_title("PSF", loc=3)

        ROISignificance.plot_sources(self.roi,ax,h,
                show_sources=self.show_sources,label_sources=self.label_sources,
                show_extension=self.show_extension,extension_color=self.extension_color)

        if filename is not None: P.savefig(filename)

class ROISmoothedSources(ROISmoothedSource):
    """ Subclass ROISmoothedSource, but add only the diffuse emission to the background. """


    def get_residual(self,**kwargs):

        residual = SmoothedResidual(self.roi,
                override_diffuse_sources=[i for i in self.roi.dsm.diffuse_sources if not hasattr(i,'skydir')],
                **kwargs)

        return residual

class ROITSMapPlotter(object):
    """ Create a residual TS map plot and overlay
        on it all of the sources in the ROI."""
    
    defaults = (
        ('size',                        5,                          ),
        ('pixelsize',               0.125,                          ),
        ('galactic',                 True,                          ),
        ('figsize',               (5.5,4),'Size of figure in inches'),
        ('fignum',                   None,'Matplotlib figure number'),
        ('title',       'Residual TS Map',       'Title of the plot'),
        ('label_sources',           False,                          ),
        ('fitsfile',                 None,                          ), 
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):

        self.roi=roi

        keyword_options.process(self, kwargs)

        from uw.like.roi_image import ROITSMapImage

        self.image=ROITSMapImage(roi,
                center=self.roi.roi_dir,
                pixelsize=self.pixelsize,
                size=self.size,
                galactic=self.galactic,
        )

        self.pf=self.image.get_pyfits()

        if self.fitsfile is not None:
            self.pf.writeto(self.fitsfile,clobber=True)

    def show(self,filename=None):

        import pylab as P
        import pywcsgrid2
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        from matplotlib.axes import Axes
        from matplotlib import mpl
        from uw.utilities import colormaps

        self.fig =  fig = P.figure(self.fignum,self.figsize)
        P.clf()

        h, d = self.pf[0].header, self.pf[0].data

        ax = pywcsgrid2.subplot(111, header=h)

        # add colorbar axes
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
        fig.add_axes(cax)

        # since this is a residual tsmap, never let the scale go below 25.
        norm=mpl.colors.Normalize(vmin=0, vmax=25) if N.max(self.image.image) < 25 else None

        im = ax.imshow(d, interpolation='nearest', origin="lower", cmap=colormaps.b, norm=norm)
        cbar = P.colorbar(im, cax=cax)

        ax.set_title(self.title)

        ax.grid(color='w',linestyle='-')

        ROISignificance.plot_sources(self.roi,ax,h,label_sources=self.label_sources,color='black')

        if filename is not None: P.savefig(filename)


class ROISmoothedModel(object):
    """ Plot (on the left) the diffuse subtracted smoothed counts and
        (on the right) the diffuse subtrcted smoothed model predicted
        counts Useful to see if your model (qualitativly) looks like
        the right source. """

    defaults = (
            ('which',            None,                                   'Source to analyze'),
            ('figsize',          (8,4),                                  'Size of the image'),
            ('fignum',            None,                           'Matplotlib figure number'),
            ('conv_type',           -1,                                    'Conversion type'),
            ('size',                 3,                          'Size of the field of view'),
            ('galactic',          True,                         'Coordinate system for plot'),
            ('show_sources',      True,                     'Put an x over all the sources.'),
            ('label_sources',    False,                           'Label sources duing plot'),
            ('kerneltype',  'gaussian',                'Type of kernel to smooth image with'),
            ('kernel_rad',         0.1,            'Sum counts/model within radius degrees.'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi

        self.cmap = colormaps.b

        # Fit many pixels inside of the summing radius
        self.pixelsize=self.kernel_rad/5.0

        self.source = roi.get_source(self.which)

        kwargs=dict(size=self.size,
                    pixelsize=self.pixelsize,
                    galactic=self.galactic,
                    conv_type=self.conv_type,
                    center=self.source.skydir,
                    per_solid_angle=True,
                    kerneltype=self.kerneltype,
                    kernel_rad=self.kernel_rad)

        # Background subtracted counts
        self.counts = SmoothedResidual(self.roi,
                override_diffuse_sources=[i for i in self.roi.dsm.diffuse_sources if not hasattr(i,'skydir')],
                **kwargs)

        # Model counts for non-background sources.
        self.model = SmoothedModel(self.roi,
                override_point_sources=self.roi.psm.point_sources,
                override_diffuse_sources=[i for i in self.roi.dsm.diffuse_sources if hasattr(i,'skydir')],
                **kwargs)

        self.model_pyfits = self.model.get_pyfits()

        self.counts_pyfits = self.counts.get_pyfits()

    def plot_counts(self):
        ax = self.grid[0]
        h, d = self.counts_pyfits[0].header, self.counts_pyfits[0].data
        im=ax.imshow(d, origin="lower", cmap=self.cmap, vmin=0, vmax=self.max_intensity)
        ax.grid()

        ROISignificance.plot_sources(self.roi,ax,h,
                show_sources=self.show_sources,label_sources=self.label_sources,
                show_extension=False)

        ax.add_inner_title("Counts", loc=2)

    def plot_model(self):
        ax = self.grid[1]
        h, d = self.model_pyfits[0].header, self.model_pyfits[0].data
        im=ax.imshow(d, origin="lower", cmap=self.cmap, vmin=0, vmax=self.max_intensity)
        ax.grid()

        cb_axes = self.grid.cbar_axes[0]
        cbar = cb_axes.colorbar(im)
        cbar.ax.set_ylabel(r'$\mathrm{counts}/[\mathrm{deg}]^2$')

        ROISignificance.plot_sources(self.roi,ax,h,
                show_sources=self.show_sources,label_sources=self.label_sources,
                show_extension=False)

        ax.add_inner_title("Model", loc=2)

    def show(self,filename=None):
        import pywcsgrid2
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.axes_grid import ImageGrid

        self.fig = P.figure(self.fignum,self.figsize)
        P.clf()

        self.max_intensity = max(
                ROISmoothedSource.get_max_intensity(self.source,i,self.roi)
                for i in [self.counts_pyfits,self.model_pyfits])

        self.grid = grid = ImageGrid(self.fig, (1, 1, 1), nrows_ncols = (1, 2),
                         axes_pad=0.1, share_all=True,
                         cbar_mode="single", cbar_pad="2%",
                         cbar_location="right",
                         axes_class=(pywcsgrid2.Axes, dict(header=self.counts_pyfits[0].header)))

        self.plot_counts()
        self.plot_model()

        if filename is not None: P.savefig(filename)

