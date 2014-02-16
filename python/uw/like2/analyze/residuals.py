"""
Residual plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/residuals.py,v 1.2 2014/02/13 04:11:16 burnett Exp $

"""

import pickle
import numpy as np
import pylab as plt
import pandas as pd
import matplotlib.gridspec as gridspec

from . import (roi_info,  analysis_base)
from .. import tools

class Residuals(roi_info.ROIinfo): 
    """Residual plots
    <br>Analysis of results from a measurement (using the uwpipeline task 'resid') of the nomalization likelihood
 functions performed for every source and every energy band. 
    """ 
    require='residuals.zip'
    def setup(self, **kwargs):
        super(Residuals, self).setup(**kwargs)

        self.plotfolder = 'residuals'
        files, self.pkls = analysis_base.load_pickles_from_zip('residuals.zip')
        self.energy = self.pkls[0]['ring']['all'].index
        self.sindec = np.sin(np.radians(np.asarray(self.df.dec,float)))
        self.singlat = np.sin(np.radians(np.array(self.df.glat, float)))

    def resid_array(self, source_name, column_name, event_type='all'):
        """ extract a component from the array of residuals
        
        source_name : str
            for now, the name of a global source
        column_name :
            one of the column names
        returns a 2-d array, shape (1728,14)
        
        """
        empty = [np.nan]*len(self.energy)
        return np.array([p[source_name][event_type][column_name] if source_name in p else empty for p in self.pkls])

    def norm_plot(self, name='isotrop', ax=None, ylim=(0.5,1.5)):
        """Isotropic Normalization vs Dec
        Only the isotropic component is allowed to vary; this is the resulting value.
        """
        lnorms =np.array([m[0] if m is not None else np.nan for m in self.diffuse_models(name)])
        high = np.abs(self.df.glat)>10
        if ax is None:
            fig,ax=plt.subplots(1,1, figsize=(10,5))
        else: fig=ax.figure
        ax.plot(self.sindec,lnorms.clip(*ylim),  '.r' , label='|b|<10')
        ax.plot(self.sindec[high],lnorms[high].clip(*ylim),  'og', label='|b|>10')
        plt.setp(ax, ylim=ylim, xlim=(-1,1), xlabel='sin(Dec)',
            ylabel='normalization factor', title='%s normalization vs sin(Dec)'%name)
        ax.grid(); ax.legend(prop=dict(size=10))
        return fig

    def pull_maps(self, source_name='ring', event_type='all', vmin=-5, vmax=5, bands=8):
        """ Maps of the residual pulls, or the normalized residuals for refitting this component alone.
        """
        nrows, ncols = ((bands+1)//4, 4 ) if bands>=4 else (1, bands)
        fig, axx = plt.subplots(nrows, ncols, figsize=(3+3*ncols,1+3*nrows), sharex=True, sharey=True)
        plt.subplots_adjust(right=0.9, hspace=0.15, wspace=0.1)
        resid = self.resid_array(source_name, 'pull', event_type=event_type)
        for ib,energy in enumerate(self.energy[:bands]):
            ax = axx.flatten()[ib]
            scat=self.basic_skyplot(ax, self.df.glon, self.singlat, resid[:,ib].clip(vmin,vmax),
                 title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, s=15, edgecolor='none', colorbar=False, labels=False)
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('normalized residual')
        fig.text(0.5,0.05, 'galactic longitude', ha='center')
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center')
        fig.text(0.5, 0.95, 'pulls for %s, %s events'%(source_name, event_type), ha='center', size=12)
        return fig
        
    def maxl_plots(self, source_name='isotrop', event_type='all', bands=8,bcut=5, ylim=(0.5,1.5)):
        """Refit normalizations per band. The peak values for the individual likelihood functions.
        """
        maxl = self.resid_array(source_name,'maxl', event_type=event_type)
        glat = np.array(self.df.glat,float)
        hilat = np.abs(glat)>bcut
        nrows, ncols = ((bands+1)//4, 4 ) if bands>=4 else (1, bands)
        fig, axx = plt.subplots(nrows, ncols, figsize=(3+3*ncols,0.5+4*nrows), sharex=True, sharey=True)
        plt.subplots_adjust(right=0.9, left=0.1, top=0.85, bottom=0.15,hspace=0.15, wspace=0.1)
        for k,ax in enumerate(axx.flatten()):
            scat=ax.scatter(self.sindec[hilat], maxl[:,k][hilat].clip(*ylim), c=abs(glat)[hilat], edgecolor='none');
            plt.setp(ax, ylim=ylim, xlim=(-1,1), title='%d MeV'%self.energy[k])
            ax.grid()
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('abs(|b|')
        fig.text(0.5,0.05, 'sin(Dec)', ha='center')
        fig.text(0.05, 0.5, 'refit normalization', rotation='vertical', va='center')
        fig.text(0.5, 0.95, 'normalizations for %s, %s events'%(source_name, event_type), ha='center', 
                 size=14)

        return fig

    def front_back_ridge(self):
        """front/back galactic residual ratio on ridge
        Ridge is within 10 degrees in latitude, 60 degrees in longitude.
        <br>Top three rows: histograms
        <br>Bottom plat: average
        """
        fb = [self.resid_array('ring', 'maxl', et) for et in ['front', 'back']]
        fbr = fb[0]/fb[1]; 
        ridge = (np.abs(self.df.glat)<10) & (np.abs(self.df.glon)<60); 
        fbrr = fbr[ridge]

        fig =plt.figure(figsize=(12,9))
        gs1 = gridspec.GridSpec(3,4, bottom=0.35)
        axx = np.array([[plt.subplot(gs1[i,j]) for i in range(3)] for j in range(4)])

        gs2 = gridspec.GridSpec(1,1, top=0.3, right=0.6)
        axb = plt.subplot(gs2[0,0])
        means = []
        for i,ax in enumerate(axx.flatten()):
            u = fbrr[:,i]
            uok = u[np.isfinite(u)]; means.append(uok.mean());
            ax.hist(u.clip(0.5,1.5), np.linspace(0.5, 1.5, 51));
            ax.axvline(1.0, color='b', ls='--')
            ax.text(0.1, 0.9, '%.0f MeV' % self.energy[i], size=10,  transform = ax.transAxes)
            
        ax=axb
        ax.semilogx(self.energy[:12], means, 'o-')
        ax.axhline(1.0, color='b', ls='--')
        ax.grid()
        plt.setp(ax, ylim=(0.85,1.15), xlabel='Energy [MeV]', ylabel='front/back ratio')
        return fig
        
    @tools.decorate_with(pull_maps, append=True)
    def pull_maps_ring(self):
        """Pull plots for galactic diffuse
        """
        return self.pull_maps('ring')
        
    @tools.decorate_with(pull_maps, append=True)
    def pull_maps_isotrop(self):
        """Pull plots for isotropic
        """
        return self.pull_maps('isotrop')
        
    @tools.decorate_with(pull_maps, append=True)
    def pull_maps_limb(self):
        """Limb """
        return self.pull_maps('limb', bands=2)
        
    @tools.decorate_with(maxl_plots, append=True)
    def maxl_plots_isotrop_back(self):
        """Max Likelihood for Isotropic back
        """
        return self.maxl_plots(event_type='back', bands=4)
        
    @tools.decorate_with(maxl_plots, append=True)
    def maxl_plots_isotrop_front(self):
        """Max Likelihood for Isotropic front
        """
        return self.maxl_plots(event_type='front', bands=4)
        
    @tools.decorate_with(maxl_plots, append=True)
    def maxl_plots_limb_back(self):
        """Max Likelihood for Limb back
        """
        return self.maxl_plots('limb', event_type='back', bands=2)
        
    @tools.decorate_with(maxl_plots, append=True)
    def maxl_plots_limb_front(self):
        """Max Likelihood for Limb front 
        """
        return self.maxl_plots('limb', event_type='front', bands=2)
    def all_plots(self):
       self.runfigures([
            self.pull_maps_ring, self.pull_maps_isotrop, 
            self.norm_plot,
            self.front_back_ridge,
            self.maxl_plots_isotrop_front, self.maxl_plots_isotrop_back,
            self.maxl_plots_limb_front, self.maxl_plots_limb_back,
            ])
