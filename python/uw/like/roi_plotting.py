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


$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_plotting.py,v 1.38 2011/04/24 03:04:18 lande Exp $

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
from . SpatialModels import Disk, Gaussian
from uw.utilities import colormaps
from uw.utilities.image import ZEA
from uw.utilities import region_writer 
from uw.utilities import keyword_options
from skymaps import SkyDir

from collections import deque

from scipy.stats import poisson,norm
from scipy.optimize import fmin,fsolve

import pylab as P
from matplotlib import rcParams,mpl,pyplot,ticker,font_manager
from matplotlib.patches import FancyArrow

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
        #print pslw.bin_centers[i]
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
        #if i < max_label:
        P.loglog(en,src[:,i],linestyle='-',marker='',label=name)
        #else:
        #    P.loglog(en,src[:,i],linestyle='-',marker='')
    for i,name in enumerate(bg_names):
        if N.any(dif[:,i]==0): continue
        P.loglog(en,dif[:,i],linestyle='-',marker='',label=name)
    #if not N.any(gal==0.):
    #    P.loglog(en,gal,linestyle='-',marker='',label='Galactic')

    #if not N.any(iso==0.):
    #    P.loglog(en,iso,linestyle='-',marker='',label='Isotropic')

    #tot = src.sum(axis=1)+iso+gal
    tot = src.sum(axis=1) + dif.sum(axis=1)
    P.loglog(en,tot,linestyle='steps-mid',label='Total Model',color='black')
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
    """ Manage the plotting of ROI info. The output is several two dimensional sky maps.
        These include counts, model counts, residual counts, and p-values. """

    defaults = (
            ('figsize',       (12,8),         'Size of the image'),
            ('fignum',             3,  'matplotlib figure number'),
            ('pixelsize',       0.25,  'size of each image pixel'),
            ('conv_type',         -1,           'Conversion type'),
            ('size',              10, 'Size of the field of view'),
            ('nticks',             5, 'Number of axes tick marks'),
            ('label_sources',  False,  'Label sources duing plot'),
            ('galactic',        True,'Coordinate system for plot'),
            ('countsfile',      '','File to save counts map'),
            ('modelfile',       '','File to save model map'),
    )

    @staticmethod
    def mathrm(st):
        return  r'$\mathrm{'+st.replace(' ','\/')+'}$'

    @staticmethod
    def matplotlib_format():
        rcParams['xtick.major.size']=10 #size in points
        rcParams['xtick.minor.size']=6
        rcParams['ytick.major.size']=10 #size in points
        rcParams['ytick.minor.size']=6
        rcParams['xtick.labelsize']=12
        rcParams['ytick.labelsize']=12
        rcParams['font.family']='serif'

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)
        
        self.roi = roi

        self.cm = CountsImage(self.roi,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)
        self.mm = ModelImage(self.roi,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        if self.countsfile is not '':
            fits=self.cm.get_pyfits()
            fits.writeto(self.countsfile,clobber=True)
        if self.modelfile is not '':
            fits=self.mm.get_pyfits()
            fits.writeto(self.modelfile,clobber=True)

        ROIDisplay.matplotlib_format()

        try:
            self.cmap_b   = colormaps.b
        except:
            self.cmap_b = mpl.cm.jet

        interactive=P.isinteractive(); P.ioff()
        fig = P.figure(self.fignum,self.figsize)
        P.clf()

        voff = -0.05
        imw = ( 0.55, 0.55, 0.55,  0.55,      0.25, 0.25) #0.38, 0.38)
        imh = ( 0.35, 0.35, 0.35,  0.35,      0.17, 0.17) # picture sizes
        imx = (-0.09, 0.23, 0.55, -0.09,      0.40, 0.40) #0.55, 0.55)
        imy = ( 0.60+voff, 0.60+voff, 0.60+voff,  0.15+voff,      0.32, 0.08)


        titles  = ['Modeled Counts', 'Observed Counts', 'Weighted Residuals', 'P-values','', '']
        xlabels = ['RA']*4 + ['p-values','weighted residuals'] 
        ylabels = ['Dec']*4 + ['Frequency','Frequency']

        for y in [titles,xlabels,ylabels]:
            for i in xrange(len(y)):
                y[i] = ROIDisplay.mathrm(y[i])

        for i,(t,xl,yl) in enumerate(zip(titles,xlabels,ylabels)):
            ax = self.__dict__['axes%d'%(i+1)] = fig.add_axes([imx[i],imy[i],imw[i],imh[i]])
            if i < 4: ax.set_aspect(1)
            ax.set_title(t)
            ax.set_xlabel(xl) 
            ax.set_ylabel(yl)
        
        # Use same scale for the counts and model map
        counts_max = N.ceil(max(N.max(self.cm.image),N.max(self.mm.image)))
        counts_min = max(N.floor(min(N.min(self.cm.image),N.min(self.mm.image))),1.0)

        self.norm1 = mpl.colors.Normalize(vmin=0,vmax=1)
        self.norm2 = mpl.colors.LogNorm(vmin=counts_min,vmax=counts_max,clip=True)
        self.norm3 = mpl.colors.Normalize(vmin=-5,vmax=5)

        #resid colorbar
        rcb_axes = fig.add_axes([imx[3]+imw[3]*0.72,imy[3],0.01,imh[3]])
        rcb        = mpl.colorbar.ColorbarBase(rcb_axes,
                                               norm=self.norm1,
                                               ticks = ticker.MaxNLocator(4),
                                               cmap  = self.cmap_b,
                                               orientation='vertical')

        #counts colorbar
        ccb_axes1 = fig.add_axes([imx[0]+imw[0]*0.72,imy[0],0.01,imh[0]])
        ccb1        = mpl.colorbar.ColorbarBase(ccb_axes1,
                                                norm=self.norm2,
                                                ticks = ticker.LogLocator(),
                                                cmap  = self.cmap_b,
                                                orientation='vertical')

        ccb_axes2 = fig.add_axes([imx[1]+imw[1]*0.72,imy[1],0.01,imh[1]])
        ccb2        = mpl.colorbar.ColorbarBase(ccb_axes2,
                                                norm=self.norm2,
                                                ticks = ticker.LogLocator(),
                                                cmap  = self.cmap_b,
                                                orientation='vertical')

        ccb_axes3 = fig.add_axes([imx[2]+imw[2]*0.72,imy[2],0.01,imh[2]])
        ccb3        = mpl.colorbar.ColorbarBase(ccb_axes3,
                                                norm=self.norm3,
                                                ticks = ticker.MaxNLocator(4),
                                                orientation='vertical')

        if interactive: P.ion()

    def model_plot(self):

        self.mm_zea = self.mm.get_ZEA(nticks=self.nticks, axes=self.axes1)
        self.axes1.imshow(self.mm_zea.image,origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm2)
        self.mm_zea.grid()
        self.mm_zea.scale_bar(color='white')
        self.plot_sources(self.mm_zea,mc='k')
        
    def counts_plot(self):

        self.cm_zea = self.cm.get_ZEA(nticks=self.nticks, axes=self.axes2)
        # note that we must explicitly clip bins with no counts to avoid confusing the colorbar with NaNs.
        self.axes2.imshow(self.cm_zea.image,origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm2)
        self.cm_zea.grid()
        self.cm_zea.scale_bar(color='white')
        self.plot_sources(self.cm_zea,mc='k')

    def resids_plot(self):

        self.resids_zea = ZEA(self.roi.roi_dir,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,axes=self.axes3)
        self.resids_zea.image = (self.cm_zea.image-self.mm_zea.image)/self.mm_zea.image**0.5
        self.axes3.imshow(self.resids_zea.image,origin='lower',interpolation='nearest',norm=self.norm3)
        self.resids_zea.grid()
        self.plot_sources(self.resids_zea,mc='k')
        
        pvals = 1 - poisson.cdf(self.cm_zea.image,self.mm_zea.image ) #0 problem?

        self.pvals_zea = ZEA(self.roi.roi_dir,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,axes=self.axes4)
        self.pvals_zea.image = pvals
        self.axes4.imshow(self.pvals_zea.image,origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm1)
        self.pvals_zea.grid()

        self.plot_sources(self.pvals_zea,mc='k')

    def hist_plot(self):

        mc = self.mm_zea.image.flatten()
        cc = self.cm_zea.image.flatten()

        pvals = 1-poisson.cdf(cc,mc)

        nb = 20
        av = float(len(pvals)) / nb
        self.axes5.hist(pvals,bins=N.linspace(0,1,20),histtype='step')
        self.axes5.axhline( av, color='red')
        lo,hi = ppf( (50.-95/2)/100., av), ppf( (50. + 95/2)/100.,av)
        self.axes5.axhspan( lo , hi , facecolor='red', alpha=0.3,label='95% Conf.')
        self.axes5.legend(loc='upper right')

        self.axes6.hist( (cc - mc)/mc**0.5, bins=N.linspace(-5,5,20), histtype='step',normed=True)
        b=N.linspace(-5,5,100)
        self.axes6.plot(b,norm.pdf(b))
        self.axes6.axvline(0,color='red')
        self.axes6.grid(True)
        self.axes6.set_xbound(lower=-5,upper=5)


    def show(self,to_screen=True,out_file=None):

        interactive=P.isinteractive(); P.ioff()

        t =self.label_sources
        self.model_plot()
        self.label_sources=False #only label the model plot
        self.counts_plot()
        self.resids_plot()
        self.hist_plot()
        self.label_sources = t

        if interactive: P.ion()

        if out_file is not None: P.savefig(out_file)
        if to_screen: P.show()

    def plot_sources(self, image, symbol='+',  fontsize=8, markersize=10, fontcolor='w', mc= 'green',**kwargs):
        """ Plot both point and diffuse sources. """
        nx = image.nx
        ny = image.ny

        src = [p for p in self.roi.get_sources() if hasattr(p,'skydir')]

        def allow(nx, ny, px, py, padx = 0.15, pady = 0.15):
            padx = padx * nx
            pady = pady * ny

            return (px > padx and px < nx - padx) and (py > pady and py < ny - pady)

        for p in src:
            x,y = image.pixel(p.skydir)
            if self.label_sources:
                padx = pady = 0.15
            else:
                padx = pady = 0.02
            if not allow(nx,ny,x,y,padx,pady): continue
            image.axes.plot([x], [y], symbol, markersize=markersize, mec = mc, mfc = mc,**kwargs)
            image.axes.plot([x], [y], 'x', markersize=markersize, mec = 'white', mfc = 'white',**kwargs)
            if self.label_sources:
                image.axes.text( x+nx/100., y+nx/100., p.name, fontsize=fontsize, color=fontcolor)




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
    """ Object to create counts slice plot. """

    defaults = (
            ('figsize',        (7,6), 'Size of the image'),
            ('fignum',             4, 'matplotlib figure number'),
            ('pixelsize',      0.125, 'size of each image pixel'),
            ('size',              10, 'Size of the field of view'),
            ('galactic',        True, 'Coordinate system for plot'),
            ('int_width',          2, 'Integration width for slice.'),
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
            ('smooth_model',    True, """Connect the model preditions with a line (instead of steps) 
                                         to make the model look smoother."""),
            ('use_gradient',    True, """Use gradient when refitting."""),
    )

    @staticmethod
    def set_color_cycle():
        P.gca().set_color_cycle(['k', 'b', 'g', 'r', 'm', 'y', 'k'])


    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        keyword_options.process(self, kwargs)

        if not type(self.oversample_factor) == int:
            raise Exception("oversample_factor must be an integer.")

        self.roi   = roi
        self.which = which

        manager,index=self.roi.mapper(self.which)
        if manager == roi.psm:
            self.source=manager.point_sources[index]
            self.pretty_name = 'Point'
        else:
            self.source=manager.diffuse_sources[index]
            self.pretty_name = self.source.spatial_model.pretty_name

        self.center=self.source.skydir

        ROIDisplay.matplotlib_format()

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

        self.mi_x = dict()
        self.mi_y = dict()

        self.mi_x[self.pretty_name]=ModelImage(self.roi,size=(self.size,self.int_width),**kwargs)
        self.mi_y[self.pretty_name]=ModelImage(self.roi,size=(self.int_width,self.size),**kwargs)

        if self.just_diffuse:
            # hide all point + extended sources.

            sources = list(self.roi.psm.point_sources) + \
                    [ i for i in self.roi.dsm.diffuse_sources if isinstance(i,ExtendedSource) ]
            # don't zero already zeroed sources
            sources = [ i for i in sources if i.model.getp(0,internal=True) != -100 ]

            ROISlice.cache_roi(self.roi)

            for source in sources: self.roi.zero_source(source)

            self.mi_x['Diffuse']=ModelImage(self.roi,size=(self.size,self.int_width),**kwargs)
            self.mi_y['Diffuse']=ModelImage(self.roi,size=(self.int_width,self.size),**kwargs)

            for source in sources: self.roi.unzero_source(source)

            ROISlice.uncache_roi(self.roi)

        if self.aspoint and isinstance(self.source,ExtendedSource):
            # replace with a point source

            ROISlice.cache_roi(self.roi)

            es=self.roi.get_source(self.which)
            
            ps=PointSource(name=es.name,model=es.model.copy(),skydir=es.spatial_model.center,leave_parameters=True)
            self.roi.add_source(ps)

            # only zero it after making a copy of the spectral part!
            self.roi.zero_source(es)

            self.roi.fit(estimate_errors=False,use_gradient=self.use_gradient)

            self.mi_x['Point']=ModelImage(self.roi,size=(self.size,self.int_width),**kwargs)
            self.mi_y['Point']=ModelImage(self.roi,size=(self.int_width,self.size),**kwargs)

            self.roi.del_source(ps)
            self.roi.unzero_source(es)

            ROISlice.uncache_roi(self.roi)

        self.models_x=dict([[name,model.image.sum(axis=0)*self.oversample_factor] for name,model in self.mi_x.items()])
        self.models_y=dict([[name,model.image.sum(axis=1)*self.oversample_factor] for name,model in self.mi_y.items()])

        self.model_dx = ROISlice.range(self.models_x[self.pretty_name],model_pixelsize,x_axis=True)
        self.model_dy = ROISlice.range(self.models_y[self.pretty_name],model_pixelsize,x_axis=False)

    @staticmethod
    def range(data,pixelsize,x_axis=False):
        if x_axis:
            return (len(data)/2.0-N.arange(len(data)))*pixelsize - pixelsize/2
        else:
            return (N.arange(len(data))-len(data)/2.0)*pixelsize + pixelsize/2

    def plot_dx(self):

        P.subplot(211)

        ROISlice.set_color_cycle()

        P.errorbar(self.counts_dx,self.counts_x,yerr=N.sqrt(self.counts_x),label='Counts', fmt='.')

        for name,model in self.models_x.items():
            P.plot(self.model_dx,model,drawstyle='steps' if not self.smooth_model else 'default',
                   label=name)

        P.gca().set_xlim(self.counts_dx[0],self.counts_dx[-1])

        P.legend(loc='upper right',numpoints=1)

        P.xlabel(ROIDisplay.mathrm('delta l' if self.galactic else 'delta ra'))
        P.ylabel(ROIDisplay.mathrm('Counts'))

    def plot_dy(self):

        P.subplot(212)

        ROISlice.set_color_cycle()

        P.errorbar(self.counts_dy,self.counts_y,yerr=N.sqrt(self.counts_y),label='Counts', fmt='.')

        for name,model in self.models_y.items():
            P.plot(self.model_dy,model,drawstyle='steps' if not self.smooth_model else 'default',
                   label=name)

        P.xlabel(ROIDisplay.mathrm('delta b' if self.galactic else 'delta dec'))
        P.ylabel(ROIDisplay.mathrm('Counts'))

        P.gca().set_xlim(self.counts_dy[0],self.counts_dy[-1])

        # only need to show legend once

    def save_data(self,data_file):
        """ Note, shrink model predicted counts to be same size as regular counts,
            for an easier file format. """

        x,y,dx,dy=['l','b','dl','db'] if self.galactic else ['ra','dec','dra','ddec']
            
        results_dict = {
            x : {
                dx: self.counts_dx.tolist(),
                'Counts': self.counts_x.tolist()
            },
            y : {
                dy: self.counts_dy.tolist(),
                'Counts': self.counts_y.tolist()
            }
        }

        of=self.oversample_factor
        for name,model in self.models_x.items():
            results_dict[x][name]=ModelImage.downsample(model,of).tolist()

        for name,model in self.models_y.items():
            results_dict[y][name]=ModelImage.downsample(model,of).tolist()

        file = open(data_file,'w')
        try:
            import yaml
            file.write(yaml.dump(results_dict))
        except:
            import pprint
            file.write(pprint.pformat(results_dict))

        file.close()

    def show(self,to_screen=True,out_file=None, data_file=None):

        self.plot_dx()
        self.plot_dy()

        title = 'Counts Slice'
        if self.source is not None: title += ' for %s' % self.source.name

        P.suptitle(title)

        if data_file is not None: self.save_data(data_file)

        if out_file is not None: P.savefig(out_file)
        if to_screen: P.show()


class ROIRadialIntegral(object):
    """ Object to create a radial integral plot, binned uniformly in theta^2. """

    defaults = (
            ('figsize',        (7,6), 'Size of the image'),
            ('size',               2, 'Size of image in degrees'), 
            ('pixelsize',     0.0625, """ size of each image pixel. This is a misleading because the
                                          size of each pixel varies, since the image is uniform in theta^2. 
                                          This value is used to determine the total number of pixels using
                                          the formula npix=size/pixelsize and represents something
                                          like an average pixelsize."""),
            ('fignum',             5, 'matplotlib figure number'),
            ('conv_type',         -1, 'Conversion type'),
            ('just_diffuse',    True, """ Display the model predictions with all point + extended
                                          sources removed. The background is not refit. """),
            ('aspoint',         True, """ Display also the model predictions for an extended source 
                                          fit with the point hypothesis. Only works when which is an
                                          extended source. """),
            ('smooth_model',    True, """Connect the model preditions with a line (instead of steps) 
                                         to make the model look smoother."""),
            ('use_gradient',    True, """Use gradient when refitting."""),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi   = roi
        self.which = which

        manager,index=self.roi.mapper(self.which)
        if manager == roi.psm:
            self.source=manager.point_sources[index]
            self.pretty_name = 'Point'
        else:
            self.source=manager.diffuse_sources[index]
            self.pretty_name = self.source.spatial_model.pretty_name

        self.center=self.source.skydir

        ROIDisplay.matplotlib_format()

        fig = P.figure(self.fignum,self.figsize)
        P.clf()

        self.get_counts()
        self.get_model()

    def get_counts(self):
        self.ci = RadialCounts(self.roi,center=self.center,size=self.size,pixelsize=self.pixelsize,conv_type=self.conv_type)
        self.counts = self.ci.image

        if N.any(N.isnan(self.counts)):
            raise Exception("Error in calculating a radial integral plot, the counts contains NaNs.")


    def get_model(self):
        kwargs=dict(center=self.center,size=self.size,pixelsize=self.pixelsize,conv_type=self.conv_type)

        self.mi = dict()

        self.mi[self.pretty_name]=RadialModel(self.roi,**kwargs)

        if self.just_diffuse:

            sources = list(self.roi.psm.point_sources) + \
                    [ i for i in self.roi.dsm.diffuse_sources if isinstance(i,ExtendedSource) ]
            # don't zero already zeroed sources
            sources = [ i for i in sources if i.model.p[0] != -100 ] 

            ROISlice.cache_roi(self.roi)

            for source in sources: self.roi.zero_source(source)

            self.mi['Diffuse']=RadialModel(self.roi,**kwargs)

            for source in sources: self.roi.unzero_source(source)

            ROISlice.uncache_roi(self.roi)

        if self.aspoint and isinstance(self.source,ExtendedSource):

            ROISlice.cache_roi(self.roi)

            es=self.roi.get_source(self.which)

            ps=PointSource(name=es.name,model=es.model.copy(),skydir=es.spatial_model.center,leave_parameters=True)
            self.roi.add_source(ps)

            # only zero it after making a copy of the spectral part!
            self.roi.zero_source(es)

            self.roi.fit(estimate_errors=False,use_gradient=self.use_gradient)

            self.mi['Point']=RadialModel(self.roi,**kwargs)

            self.roi.del_source(ps)
            self.roi.unzero_source(es)

            ROISlice.uncache_roi(self.roi)

        self.theta_sqr=self.ci.bin_centers_deg
        self.models=dict([[name,model.image] for name,model in self.mi.items()])

        for name,model in self.models.items():
            if N.any(N.isnan(model)):
                raise Exception("Error in calculating a radial integral, model %s contains NaNs." % name)

    def save_data(self,data_file):

        results_dict = dict(
            Radius=self.theta_sqr.tolist(),
            Counts=self.counts.tolist()
        )
        for name,model in self.models.items():
            results_dict[name]=model.tolist()

        file = open(data_file,'w')
        try:
            import yaml
            file.write(yaml.dump(results_dict))
        except:
            import pprint
            file.write(pprint.pformat(results_dict))
        file.close()

    def show(self,to_screen=True,out_file=None,data_file=None):

        ROISlice.set_color_cycle()

        P.errorbar(self.theta_sqr,self.counts,yerr=N.sqrt(self.counts),label='Counts', fmt='.')

        for name,model in self.models.items():
            P.plot(self.theta_sqr,model,drawstyle='steps' if not self.smooth_model else 'default',
                   label=name)

        P.legend(loc='upper right',numpoints=1)

        P.xlabel(ROIDisplay.mathrm(r'delta \theta^2 [deg^2]'))
        P.ylabel(ROIDisplay.mathrm('Counts'))

        title = 'Radial Integral Counts'
        if self.source is not None: title += ' for %s' % self.source.name
        P.title(title)

        if data_file is not None: self.save_data(data_file)

        if out_file is not None: P.savefig(out_file)
        if to_screen: P.show()


class ROISignificance(object):
    """ Make a plot of the statistical significance (D-M)/sqrt(M)
        where D and M are the measured counts and the model predicted
        counts integrated within a circual aperature of radius . """

    defaults = (
            ('figsize',        (5,5),         'Size of the image'),
            ('fignum',             3,  'matplotlib figure number'),
            ('conv_type',         -1,           'Conversion type'),
            ('size',               5, 'Size of the field of view'),
            ('galactic',        True,'Coordinate system for plot'),
            ('kernel_rad',       0.25, 'Sum counts/model within radius degrees.'),
            ('label_sources',  False,  'Label sources duing plot'),
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
    def plot_sources(roi, ax, header, label_sources, color='black'):
        """ Add to the pywcsgrid2 axes ax any sources in the ROI which
            have a center. Also, overlay the extended source shapes
            if there are any. Note, this function requires pyregion. 
            
            The code to do this was found at:
                http://leejjoon.github.com/matplotlib_astronomy_gallery/tutorial/w51c.html
                
        """

        import pyregion
        region_string = region_writer.get_region(roi,color=color,label_sources=label_sources, show_localization=False)
        reg = pyregion.parse(region_string).as_imagecoord(header)
        patch_list, artist_list = reg.get_mpl_patches_texts()

        for p in patch_list: ax.add_patch(p)
        
        for t in artist_list: 
            # make the markers bigger
            if hasattr(t,'set_markersize'): t.set_markersize(4*t.get_markersize())
            ax.add_artist(t)

    def show(self,to_screen=True,out_file=None):

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

        ROISignificance.plot_sources(self.roi,ax,h,self.label_sources)

        if out_file is not None: P.savefig(out_file)
        if to_screen: P.show()

class ROISmoothedSource(object):
    """ Make a smoothed residual plot to look at a paritcular source. 
    
        This object requires pywcsgrid2. To get it, go to
        http://leejjoon.github.com/pywcsgrid2/users/overview.html """

    defaults = (
            ('which',           None,    'Draw the smoothed point version of this source.'),
            ('figsize',        (5,5),                                  'Size of the image'),
            ('fignum',             3,                           'Matplotlib figure number'),
            ('conv_type',         -1,                                    'Conversion type'),
            ('size',               3,                          'Size of the field of view'),
            ('galactic',        True,                         'Coordinate system for plot'),
            ('overlay_psf',     True, 'Add a smoothed reference PSF on top of the counts.'),
            ('psf_size',           1,                         'Size of the PSF insert box'), 
            ('psf_loc',            4,                       'Location to put the psf box.'), 
            ('label_sources',  False,  'Label sources duing plot'),
            ('kerneltype','gaussian',  'Type of kernel to smooth image with'),
            ('kernel_rad',       0.1,            'Sum counts/model within radius degrees.'),
            ('title',           None,            'Title for the plot'),
    )

    def get_residual(self,**kwargs):
        """ Allow the particular method for getting the residual image to be overloaded. """

        self.roi.zero_source(which=self.which)

        residual = SmoothedResidual(self.roi,**kwargs)

        self.roi.unzero_source(which=self.which)

        return residual

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        if self.which is None: raise Exception("Parameter soruce must be specified.")
        
        self.roi = roi

        self.cmap = colormaps.b

        # Fit many pixels inside of the summing radius
        self.pixelsize=self.kernel_rad/10.0

        self.source = roi.get_source(self.which)

        smoothed_kwargs=dict(size=self.size,
                    pixelsize=self.pixelsize,
                    galactic=self.galactic,
                    conv_type=self.conv_type,
                    center=self.source.skydir,
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
            # Normalize psf to have same maximum pixel scale
            # as residual image.
            self.psf_pyfits[0].data *= \
                    N.max(self.residual_pyfits[0].data)/\
                    N.max(self.psf_pyfits[0].data)

    def show(self,to_screen=True,out_file=None):
        import pywcsgrid2
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
        from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

        self.fig = P.figure(self.fignum,self.figsize)
        P.clf()

        h, d = self.residual_pyfits[0].header, self.residual_pyfits[0].data

        grid = ImageGrid(self.fig, (1, 1, 1), nrows_ncols = (1, 1),
                         cbar_mode="single", cbar_pad="2%",
                         cbar_location="right",
                         axes_class=(pywcsgrid2.Axes, dict(header=h)))

        self.ax = ax = grid[0]

        im=ax.imshow(d, origin="lower", cmap=self.cmap)

        cb_axes = grid.cbar_axes[0] # colorbar axes

        cb_axes.colorbar(im)

        cb_axes.axis["right"].toggle(ticklabels=False)
        
        ax.grid()

        if self.title: ax.set_title(self.title)

        if self.overlay_psf:
            h_psf, d_psf = self.psf_pyfits[0].header, self.psf_pyfits[0].data
            axins = zoomed_inset_axes(ax, zoom=1, loc=self.psf_loc,
                              axes_class=pywcsgrid2.Axes,
                              axes_kwargs=dict(wcs=h_psf))

            # Note, match color maps with parent.
            im2=axins.imshow(d_psf, cmap=im.cmap,
                             origin="lower")
            axins.axis[:].toggle(all=False)

            axins.add_inner_title("PSF", loc=3)

        ROISignificance.plot_sources(self.roi,ax,h,self.label_sources)

        if out_file is not None: P.savefig(out_file)
        if to_screen: P.show()

class ROISmoothedSources(ROISmoothedSource):
    """ Subclass ROISmoothedSource, but add only the diffuse emission to the background. """


    def get_residual(self,**kwargs):

        residual = SmoothedResidual(self.roi,
                override_diffuse_sources=[i for i in self.roi.dsm.diffuse_sources if not hasattr(i,'skydir')],
                **kwargs)

        return residual
