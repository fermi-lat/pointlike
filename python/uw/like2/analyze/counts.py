"""
Count plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/counts.py,v 1.1 2013/08/20 14:29:32 burnett Exp $

"""

import pickle
import numpy as np
import pylab as plt
import pandas as pd

from . import analysis_base
from . analysis_base import FloatFormat

class CountPlots(analysis_base.AnalysisBase): 
    """Count plots
    <br> Plots generated after each iteration, checking quality of counts histogram
    %(iteration_info)s
    """ 
    require='pickle.zip'
    def setup(self, **kwargs):
        self.plotfolder = 'counts'

         # get the basic pickles with the model
        files, pkls = self.load_pickles()
        self.pkls = pkls # for development
        assert len(pkls)==1728, 'expect to find 1728 pickled roi files'
        #sdirs = [r['skydir'] for r in pkls]
        #glon = np.array([r['skydir'].l() for r in pkls]); 
        #glon[glon>180]-=360
        #glat = np.array([r['skydir'].b() for r in pkls])
        #singlat = np.sin(np.radians(glat))
        #self.roinames=roinames = [p['name'] for p in pkls]
        def chisq10(counts):
            total, observed = counts['total'], counts['observed']
            return ((observed-total)**2/total)[:8].sum()
        def lat180(l): return l if l<180 else l-360
        rdict = dict()
        for r in pkls:
            rdict[r['name']]=dict(
                glon = lat180(r['skydir'].l()),
                glat = r['skydir'].b(),
                chisq = r['counts']['chisq'],
                chisq10= chisq10(r['counts']),
                last_diff= r['logl']-r['prev_logl'][-1] if 'prev_logl' in r else np.nan,
                n_iter = int(len(r['prev_logl'])+1 if 'prev_logl' in r else 0),
                )
        self.rois = pd.DataFrame(rdict).transpose()
        self.rois['singlat'] = np.sin(np.radians(self.rois.glat))
        # dict of dataframes with count info. columns are energies
        self.energy = pkls[0]['counts']['energies'] # extract list from first pickle
        counts = [p['counts'] for p in pkls]
        self.counts=dict()
        for key in ['observed', 'total']:
            self.counts[key]= pd.DataFrame([x[key] for x in counts], index=self.rois.index)
        try:
            self.add_model_info()
        except Exception, msg:
            print msg
        iters = self.rois.n_iter #np.array([ len(p['prev_logl']) for p  in self.pkls])
        logl = np.array([p['prev_logl']+[p['logl']] if 'prev_logl' in p else 0 for p  in self.pkls])
        def logsum(n):
            return sum([x[n] if n<len(x) else x[-1] for x in logl])
        k = int(iters.max())
        t=np.array(map(logsum, range(k)))
        def log_delta(n):
            return np.array([x[n]-x[n-1] if n<len(x) else 0 for x in logl])
        rois = np.histogram(self.rois.n_iter, range(k))[0]
        rois[rois==0]=1728; 
        config = eval(open('config.txt').read()) 
        ihist = str((t[1:]-t[:-1]).round(1))
        ihist = pd.DataFrame(dict( loglike=t,  rois=rois, delta_sum=list(t[1:]-t[:-1]), 
                                  delta_min= [log_delta(i).min() for i in range(1,k)],
                                  delta_max= [log_delta(i).max() for i in range(1,k)],
                    ), 
                columns='rois delta_min delta_max delta_sum'.split(),
                index=range(1,k))
        ihist.index.name='iteration'
        input_model=config['input_model']['path']
        self.iteration_info = """<p>Input model: <a href="../../%s/plots/index.html?skipDecoration">%s</a>
        <p>Minimum, maximum numbers of iterations: %d %d 
        <p>Iteration history: log likelihood change for each step: \n%s
        """ % (input_model,input_model, iters.min(), iters.max(), 
                ihist.T.to_html(float_format=FloatFormat(1)) )
    
    def add_model_info(self):
        for i,key in enumerate(['ring','isotrop', 'SunMoon', 'limb',]): # the expected order
            t = []
            for j,p in enumerate(self.pkls):
                if key in p['diffuse_names']:
                    y=p['counts']['models'][i]
                    assert y[0]==key, 'wrong key, roi %d: %s!=%s; list is %s'% (j,key, y[0], p['diffuse_names'])
                    t.append(y[1])
                else:
                    t.append(np.zeros(len(self.energy)))
            self.counts[key]= pd.DataFrame(t, index=self.rois.index)

    def counts_map(self):
        """ Sum, for E>100 Mev
        """
        obs = self.counts['observed']
        total = np.array([sum(x[1]) for x in obs.iterrows()])

    def residual(self, ib):
        """ residual DF array for energy band ib 
        """
        obs   = self.counts['observed'].transpose().ix[ib]
        model = self.counts['total'].transpose().ix[ib]
        resid = (obs-model)/np.sqrt(model)
        return resid
     
    def residual_hists(self):
        """ histograms of normalized residuals 
        subset for ridge (|b|<10, |l|<60) shown
        """
        fig,axx = plt.subplots(3,4, figsize=(12,12))
        ridge = ( np.abs(self.rois.glat)<10) * ( np.abs(self.rois.glon)<60 )

        for ib,ax in enumerate(axx.flatten()):
            resid = self.residual(ib)
            ax.hist(resid.clip(-5,5), np.linspace(-5,5,21))
            ax.hist(resid[ridge].clip(-5,5), np.linspace(-5,5,21))
            ax.set_title('%.0f MeV'% self.energy[ib], fontsize=10)
            ax.axvline(0, color='k')
            plt.setp(ax, xlim=(-5,5))
            ax.grid(True)
        return fig
    
    def residual_plot(self):
        """ plot of the average normalized residual
        """
        res = [self.residual(ib) for ib in range(len(self.energy))]
        means = [x.mean() for x in res]
        stds  = [x.std() for x in res]
        fig,ax= plt.subplots(figsize=(4,4))
        ax.errorbar(self.energy, y=means, yerr=stds, fmt='o')
        plt.setp(ax, xscale='log', xlabel='energy', ylabel='average residual')
        ax.set_title('count residuals')
        ax.grid()
        ax.axhline(0, color='k')
        
    def chisq_plots(self, use10=False, hsize=(1.0, 0.8, 1.5, 0.5), vmin=0, vmax=50, bcut=10):
        """ chi squared plots
        chi squared distribution
        """
        #fig, axs = plt.subplots( 1,2, figsize=(8,3))
        #plt.subplots_adjust(wspace=0.3)
        fig, axs = self.subplot_array( hsize, figsize=(11,5))
        chisq = self.rois.chisq if not use10 else self.rois.chisq10

        def chisky(ax):
            self.basic_skyplot(ax, self.rois.glon, self.rois.singlat, chisq, 
                s=60, vmin=vmin, vmax=vmax,  edgecolor='none', colorbar=True);
                
        def chihist(ax):
            bins = np.linspace(0, vmax, 26)
            lolat = np.abs(self.rois.glat)<bcut
            ax.hist(chisq.clip(0,vmax), bins, label='all: mean=%.1f'%chisq.mean())
            ax.hist(chisq.clip(0,vmax)[lolat], bins, color='red', label='|b|<%d (%.1f)'%(bcut, chisq[lolat].mean()))
            ax.legend(loc='upper right', prop=dict(size=10)) 
            plt.setp(ax, xlabel='chisq', xlim=(0,vmax))
            ax.grid(True)
            
        for f, ax in zip( (chihist, chisky), axs.flatten()): f(ax)
        return fig
        
    def residual_maps(self, vmin=-5, vmax=5):
        """ Maps of the residuals 
        """
        fig, axx = plt.subplots(3,4, figsize=(12,10), sharex=True, sharey=True)
        plt.subplots_adjust(right=0.9, hspace=0.15, wspace=0.1)
        for ib,energy in enumerate(self.energy[:12]):
            ax = axx.flatten()[ib]
            scat=self.basic_skyplot(ax, self.rois.glon, self.rois.singlat, self.residual(ib).clip(vmin,vmax),
                 title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, s=15, edgecolor='none', colorbar=False, labels=False)
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('normalized residual')
        fig.text(0.5,0.05, 'galactic longitude', ha='center')
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center')
        return fig
        
    def resid_vs_dec(self, ib=0, ax=None, ylim=(-8,8), labels=True):
        """ residual vs. declination angle
        """
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        r = self.residual(ib).clip(*ylim)
        ax.plot(self.rois.dec, r, '.', color='gray')
        galplane = np.abs(self.rois.glat)<5
        ax.plot(self.rois.dec[galplane], r[galplane], 'or', label='|b|<5')
        ax.axhline(0, color='k')
        plt.setp(ax, ylim=ylim, xlim=(-90,90))
        if labels: plt.setp(ax, xlabel='Dec', ylabel='normalized residual')
        ax.set_xticks((-90,-45,0,45,90))
        ax.set_title('%d MeV' % self.energy[ib], size=10)
        ax.legend(prop=dict(size=10))
        ax.grid()
        return fig
        
    def ridge_spectral_residuals(self, ax=None, glat=(-10,10), glon=(-60,60), nolabels=False):
        """ Spectral residuals along the Galactic ridge
       
        Designed to match, except for different ROI definitions, the Saclay standard plot.
        """
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        def cut( x, range):
            return (x>range[0])*(x<range[1])
        #ridge = ( np.abs(self.rois.glat)<10) * ( np.abs(self.rois.glon)<60 )
        ridge = cut(self.rois.glat, glat) * cut(self.rois.glon, glon)
        data =self.counts['observed'][ridge].sum()
        model = self.counts['total'][ridge].sum()
        x = self.energy
        y = data/model-1
        yerr = 1/np.sqrt(model) # needs correlation factor
        ax.errorbar(x, y, yerr=yerr, fmt='o')
        plt.setp(ax, xscale='log')
        if not nolabels:
            plt.setp(ax, xlabel='energy (MeV)', ylabel='(counts-model)/model')
        ax.axhline(0, color='gray')
        ax.grid()
        return fig

    def all_plots(self):
        self.runfigures([
            self.chisq_plots,
            self.residual_maps, 
            self.residual_plot, 
            self.residual_hists, 
            self.ridge_spectral_residuals,
            ])