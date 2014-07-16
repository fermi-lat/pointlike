"""
Residual plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/residuals.py,v 1.6 2014/03/14 13:40:44 burnett Exp $

"""

import os, glob, pickle
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
        r= np.array([list(p[source_name][event_type][column_name]) if source_name in p else empty for p in self.pkls])
        assert r.shape==(1728,14), 'Failed shape requirement'
        return r
    def update_correction(self, vin, vout=None):
        """ update the Galactic Diffuse correction factor array with new residuals"""
        if vout is None: vout = vin+1
        # get current residual array, replace any nan's with 1.0
        t = self.resid_array('ring', 'maxl')
        ra = t[:,:8] # only upto 10 GeV
        ra[~np.isfinite(ra)]=1.0
        # read in input correction array
        infile = os.path.expandvars('$FERMI/diffuse/galactic_correction_v%d.csv'%vin)
        assert os.path.exists(infile), 'File %s not found' %infile
        cv_old = pd.read_csv(infile, index_col=0)
        # multiply input corrections by residuals ane write it out
        cv_new = ra * cv_old
        outfile = os.path.expandvars('$FERMI/diffuse/galactic_correction_v%d.csv'%vout)
        cv_new.to_csv(outfile)
        print 'wrote new diffuse correction file %s' % outfile
    
    def norm_plot(self, name='isotrop', ax=None, ylim=(0.5,1.5)):
        """Isotropic Normalization vs Dec
        Only the isotropic component is allowed to vary; this is the resulting value.
        """
        lnorms =np.array([m[0] if m is not None else np.nan for m in self.diffuse_models(name)])
        high = np.array(np.abs(self.df.glat)>10, bool)
        if ax is None:
            fig,ax=plt.subplots(1,1, figsize=(10,5))
        else: fig=ax.figure
        ax.plot(self.sindec,      lnorms.clip(*ylim),  '.r' , label='|b|<10')
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
        
    def maxl_plots(self, source_name='isotrop', event_type='all', bands=8,bcut=5, ylim=(0.5,1.5), nocolorbar=False):
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
        if nocolorbar: return fig
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
        <br>Bottom plot: average
        """
        fb = [self.resid_array('ring', 'maxl', et) for et in ['front', 'back']]
        fbr = fb[0]/fb[1]; 
        ridge = (np.abs(self.df.glat)<10) & (np.abs(self.df.glon)<60); 
        fbrr = fbr[ridge]

        fig =plt.figure(figsize=(12,9))
        gs1 = gridspec.GridSpec(3,4, bottom=0.35)
        axt = []; first=True
        for i in range(3):
            for j in range(4):
                axt.append(plt.subplot(gs1[i,j], sharex=None if first else axt[0] ))
                first=False
        
        gs2 = gridspec.GridSpec(1,1, top=0.3, right=0.6)
        axb = plt.subplot(gs2[0,0])
        means = []
        for i,ax in enumerate(axt):
            u = fbrr[:,i]
            uok = u[np.isfinite(u)]; means.append(uok.mean());
            ax.hist(u.clip(0.5,1.5), np.linspace(0.5, 1.5, 51));
            ax.axvline(1.0, color='b', ls='--')
            ax.text(0.1, 0.9, '%.0f MeV' % self.energy[i], size=10,  transform = ax.transAxes)
        plt.setp(axt[0] , xlim=(0.5,1.5))
   
        ax=axb
        ax.semilogx(self.energy[:12], means, 'o-')
        ax.axhline(1.0, color='b', ls='--')
        ax.grid()
        plt.setp(ax, ylim=(0.85,1.15), xlabel='Energy [MeV]', ylabel='front/back ratio')
        return fig
        
    def front_back_strong(self):
        """Front/back ratio for strongest sources
        """
        sources = glob.glob('sources_*.csv')
        assert len(sources)>0, 'no sources_*.csv files fouind'
        filename = sources[0]
        assert os.path.exists(filename), 'sources csv file not found'
        sdf = pd.read_csv(filename, index_col=0)
        t =sdf.sort_index(by='ts')
        strong_names = t.ts[-4:].index
        print 'selected strong names: %s' % list(strong_names)
        roinames = [sdf.ix[sn]['roiname'] for sn in strong_names]
        rois = map(lambda s: int(s[-4:]), roinames)
        
        rdict = {}; edict={}
        for roi,name in zip(rois, strong_names):
            t = self.pkls[roi][name]
            u =np.array([t[et]['flux'] for et in ('front','back')]); 
            du = np.array([ t[et]['uflux']-t[et]['flux'] for et in ('front', 'back')]) 
            rdict[name] = u[0,:8]/u[1,:8]
            edict[name] = np.sqrt((du[0,:8]/u[0,:8])**2 +
                                  (du[1,:8]/u[1,:8])**2 )

        fig, axx = plt.subplots(2,2, figsize=(10,10), sharex=True, sharey=True)
        for (key, values), errors, ax  in zip(rdict.items(),edict.values(), axx.flatten()):
            ax.errorbar(x=self.energy[:8],y=values,yerr=errors, fmt='o')
            ax.grid()
            ax.axhline(1.0, color='k', ls ='--')
            plt.setp(ax, xlabel='Energy [MeV]', xscale='log', title=key, ylim=(0.85,1.15));
        return fig
        
    def isotropic_hists(self, bmin=10, xlim=(0.4,1.6), etnames = ('front', 'back', 'all') ):
        """Isotropic normalization hists
        Maximum likelihood values for normalization for front (green) and back (red), all (black)
        """
        maxl = [self.resid_array('isotrop', 'maxl', et) for et in etnames ]
        gcut = abs(self.df.glat)>bmin; 
        fig,axx = plt.subplots(2,4, figsize=(12,8), sharex=True)
        for i, ax in enumerate(axx.flatten()):
            z = [maxl[j][gcut,i].clip(*xlim) for j in range(len(etnames))]
            ax.hist(z, np.linspace( *xlim), histtype='step', color=('g','r', 'k'), label=etnames);
            ax.text(0.1,0.9, '%d'%self.energy[i],  transform = ax.transAxes)
            ax.axvline(1.0, ls = '--')
        fig.text(0.5, 0.05, 'normalization factor', ha='center')
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
            self.front_back_ridge, self.front_back_strong,
            self.isotropic_hists,
            self.maxl_plots_isotrop_front, self.maxl_plots_isotrop_back,
            self.maxl_plots_limb_front, self.maxl_plots_limb_back,
            ])
