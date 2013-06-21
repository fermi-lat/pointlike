"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/uwsourcecomparison.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import os
import numpy as np
import pylab as plt
import pandas as pd

from . import sourceinfo

class UWsourceComparison(sourceinfo.SourceInfo):
    """Comparision with another UW model: %(othermodel)s
    <br>Ratios are %(skymodel)s/%(othermodel)s.
    
    """
    def setup(self, othermodel='uw25'):
        super(UWsourceComparison,self).setup()
        self.plotfolder = 'comparison_%s' % othermodel

        otherfilename = '../%s/sources.pickle' %othermodel
        self.othermodel=othermodel
        assert os.path.exists(otherfilename), 'File %s not found' % otherfilename
        print 'loading %s' % otherfilename
        self.odf = pd.load(otherfilename)

    def compare(self, scat=True):
        """Ratios of values of various fit parameters
        """
        self.df['pindex_old']=self.odf.pindex
        self.df['ts_old'] = self.odf.ts
        self.df['eflux_old']=self.odf.eflux
        self.df['a_old'] = self.odf.a
        odf, df = self.odf, self.df
        plane = np.abs(df.glat)<5
        
        def plot_ratio(ax, y, cut, ylim, ylabel):
            ts = df.ts
            if scat:
                ax.semilogx(ts[cut], y.clip(*ylim)[cut], '.')
                ax.semilogx(ts[cut*plane], y.clip(*ylim)[cut*plane], '+r',label='|b|<5')
                plt.setp(ax, xlim=(10,1e4), ylim=ylim, ylabel=ylabel,)
            else:
                bins = np.logspace(1,4,13)
                x = np.sqrt(bins[:-1]*bins[1:])
                t = ts[cut]
                u = y.clip(*ylim)[cut]
                ybinned = np.array([ u[ (t >= bins[i])*(t < bins[i+1])] for i in range(len(bins)-1)])
                ymean = [t.mean() for t in ybinned]
                yerr = [t.std()/np.sqrt(len(t)) if len(t)>1 else 0 for t in ybinned] 
                ax.errorbar(x, ymean, yerr=yerr, fmt='o')
                sy = lambda y: 1+(y-1)/4.
                plt.setp(ax, xlim=(10,1e4), ylim=(sy(ylim[0]),sy(ylim[1])), ylabel=ylabel, xscale='log')
            ax.axhline(1, color='k')
            if scat: ax.legend(prop=dict(size=10))
            ax.grid()
            yhigh = y[cut*(df.ts>200)]
            print '%-6s %3d %5.3f +/- %5.3f ' % (ylabel, len(yhigh), yhigh.mean(),  np.sqrt(yhigh.std()/len(yhigh)))

        def plot_pindex(ax,  ylim=(0.9,1.1)):
            cut=df.beta<0.01
            y = df.pindex/df.pindex_old
            plot_ratio(ax, y, cut, ylim, 'index')
            ax.set_xlabel('TS')
        def plot_ts(ax, rdts=(0.5,1.5)):
            y = df.ts/(df.ts_old +0.1)
            cut=df.ts>10
            plot_ratio(ax, y, cut, rdts,  'TS')
        def plot_flux(ax, ylim=(0.5, 1.5)):
            y = df.eflux/df.eflux_old
            cut = df.ts>10
            plot_ratio(ax, y, cut, ylim, 'Eflux')
            
        def plot_semimajor(ax, ylim=(0.5,1.5)):
            y = df.a/df.a_old
            cut = df.ts>10
            plot_ratio(ax, y, cut, ylim, 'r95')
            
        fig, ax = plt.subplots(4,1, figsize=(12 if scat else 8,12), sharex=True)
        plt.subplots_adjust(hspace=0.05, left=0.1, bottom=0.1)
        for f, ax in zip((plot_ts, plot_flux, plot_pindex,plot_semimajor,), ax.flatten()): f(ax)
        fig.text(0.5, 0.05, 'TS', ha='center')
        return fig
    
    def compare_profile(self):
        """Ratios of values of various fit parameters: profile plots
        Same as the first plot, but showing means and errors for bins in TS.
        <br>Changes below TS=25 are due to a threshold selection effect, a bias towared higher TS for the Clean, as can be seen
        in the TS scatter plot above.
        """
        return self.compare(False)
    
    def band_compare(self):
        """Band flux ratios
        For each of the 12 energy bands from 100 MeV to 100 GeV, plot ratio of fits for each source in common
        """
        fnow=self.df.sedrec[0].flux
        fold=self.odf.sedrec[0].flux
        self.df['sedrec_old'] = self.odf.sedrec
        fnew = np.array([s.flux for s in self.df.sedrec])
        fold = np.array([s.flux if type(s)!=float else [np.nan]*14 for s in self.df.sedrec_old])
        energy = np.logspace(2.125, 5.125, 13) # 4/decade
        
        fig, axx = plt.subplots(3,4, figsize=(14,12), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.1, bottom=0.1)
        def plotone(ib, ax):
            ok = fold[:,ib]>1
            r =fnew[:,ib][ok]/fold[:,ib][ok]
            ax.plot(fold[:,ib][ok].clip(1,1e3), r, '.');
            plt.setp(ax, xscale='log', ylim=(0.5,1.5))
            ax.axhline(1.0, color='k')
            ax.text( 50, 1.4 ,'%d MeV'% energy[ib], fontsize=12)

        for ib,ax in enumerate(axx.flatten()): plotone(ib, ax)
        axx[0,0].set_xlim(1,1e3)
        fig.text(0.5, 0.05, 'Energy flux (eV/cm**2/s)', ha='center')
        
        return fig
        
    def quality_comparison(self):
        """FIt quality comparison
        Compare the spectral fit quality of the reference model with this one. All sources in common with TS>50, are shown.
        """
        fig, ax = plt.subplots(figsize=(5,5))
        lim=(0,30)
        cut = self.df.ts>50
        self.df['fitqual_old'] = self.odf.fitqual
        new, old = self.df.fitqual.clip(*lim)[cut], self.df.fitqual_old.clip(*lim)[cut]
        ax.plot(old, new ,  '.')
        ax.plot(lim, lim, color='r')
        plt.setp(ax, xlim=lim, ylim=lim, xlabel='old fit quality', ylabel='new fit quality')
        return fig
        
    def all_plots(self):
        self.runfigures([self.compare, self.compare_profile, self.band_compare, self.quality_comparison ])