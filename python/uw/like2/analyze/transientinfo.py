"""
Analysis plots of transient sources

$Header$

"""

import os, pickle, pyfits, glob
from uw.like2.analyze import (sourceinfo, associations,)
from pointlike import IntVector
from skymaps import Band, SkyDir

from collections import Counter
import numpy as np
import pylab as plt
import pandas as pd

from analysis_base import html_table, FloatFormat
from . import sourceinfo, localization, associations

class TransientInfo(sourceinfo.SourceInfo): 
    """Transient source properties 
    """
    require='pickle.zip'
    def setup(self, **kw):
        super(TransientInfo, self).setup(**kw)
        self.plotfolder='transients' #needed by superclass
       
        print 'transients:', sum(self.df.transient)
        self.df = self.df.ix[self.df.transient]
        assert len(self.df)>0, 'No transients found'
        self.loc=localization.Localization(df=self.df) # to make plots using code
        self.assoc = associations.Associations(df=self.df)
    
    def cumulative_ts(self):
        """ Cumulative test statistic TS
        
        A logN-logS plot, but using TS. Important thresholds at TS=10 and 25 are shown.
        """

        df = self.df
        df['prefix'] = [n[:2] for n in df.index]
        fig=super(TransientInfo, self).cumulative_ts(check_localized=False, 
                    ts = df.ts, label='all',
                    other_ts=[df.ts[df.prefix=='PG'],
                               df.ts[(df.prefix=='Sh') | (df.prefix=='S ')|(df.prefix=='TS')]] ,
                    other_label=['PGWave', 'TSmap'], 
                    legend=True,
                    );
        plt.setp(fig.axes[0], ylim=(1,1000), xlim=(9,100))

        return fig
    
   
    def association_summary(self):
        """Summary
        %(summary_html)s
        """
        self.assoc.summary()
        self.summary_html = self.assoc.summary_html
    
    def bzcat_study(self):
        """AGN analysis
        %(self.bzcat_html)s
        """
        fig = self.assoc.bzcat_study(tsmax=200, title='Cumulative bzcat associations')
        self.bzcat_html = self.assoc.bzcat_html
        return fig

    
    def all_plots(self):
        self.runfigures([
            self.census, self.cumulative_ts, self.non_psr_spectral_plots, self.fit_quality, self.pivot_vs_e0, 
            self.loc.localization, self.loc.locqual_hist, 
            self.association_summary, self.bzcat_study,
            ])
        

class ExtTransientInfo(TransientInfo):
    """ subclass invoked with a specific path
    """
    def __init__(self, model_path, quiet=True):
        self.curdir = os.getcwd()
        try:
            self.model_path=model_path
            os.chdir(model_path)
            if not quiet: print os.getcwd()
            self.skymodel = os.path.split(os.getcwd())[-1]
            self.setup(refresh=False, quiet=quiet)
            self.plotfolder=model_path+'/plots/'+self.plotfolder
        finally:
            os.chdir(self.curdir)
    
    def all_plots(self):
        os.chdir(self.model_path)
        try:
            super(ExtTransientInfo, self).all_plots()
        finally:
            os.chdir(self.curdir)
        
B12 = Band(12)
def neighbors(i):
    iv = IntVector()
    n = B12.findNeighbors(int(i),iv)
    return list(iv)

def legend(ax, **kw):
    leg = ax.legend(**kw)
    for box in leg.get_patches():
        box._height=0; box._y=0.5

class Analysis(object):
    """ a place to put code to make plots for all months
    """
    def __init__(self, last=72,  path='$FERMI/skymodels/P301_monthly'):
        os.chdir(os.path.expandvars(path))
        files =  sorted(glob.glob('month*/sources.pickle'));
        self.monthlist = [f.split('/')[0] for f in files];
        self.monthinfo=monthinfo=[]
        for month in self.monthlist[:last+1]:
            sinfo = sourceinfo.ExtSourceInfo(month, quiet=True)
            monthinfo.append(sinfo) 
            assoc = associations.ExtAssociations(month,quiet=True)
            for key in 'aprob acat aname aang adeltats'.split():
                sinfo.df[key] = assoc.df[key]
        # concatenate the transient sources
        dflist= []
        for i,month in enumerate(monthinfo):
            month.df['month'] = i+1
            month.df['has_assoc'] = [a is not None for a in month.df.associations]
            dflist.append( month.df[month.df.transient]['ra dec glat glon ts pindex eflux flux a locqual aprob acat adeltats has_assoc month'.split()])
        df = pd.concat(dflist) 
        df['skydir'] = map(SkyDir, df.ra, df.dec)
        df['roi'] = map( lambda s: Band(12).index(s), df.skydir)
        
        # select those with TS>10 and good localization
        good = np.logical_not((df.ts<10) | (df.locqual>8) | (df.a<0.01) | (df.a>0.25))
        self.df = df[good]
        print 'good sources: %d/%d' % (len(self.df), len(df))
        self.hilat = np.abs(self.df.glat)>10
        print 'High latitude (>10 deg)', sum(self.hilat)
        self.assoc = np.isfinite(self.df.aprob)
        print 'Associated', sum(self.assoc)
        
def monthly_info(self, sname):
    """Return monthly info for sources in all or most months
    """
    s6y = pd.read_csv('../P301_6years/uw972/sources_uw972.csv', index_col=0)
    assert sname in s6y.index, 'Source %s not found in 6-year list' % sname
    rec6y = s6y.ix[sname]
    eflux6y = rec6y.flux * rec6y.e0**2*1e6
    fluxunc67 = rec6y.flux_unc
    # compile monthly stuff
    months = dict()
    for i,month in enumerate(self.monthinfo):
        if sname not in month.df.index: continue
        rec = month.df.ix[sname]
        eflux = rec.flux * rec.e0**2*1e6
        relfluxunc = rec.flux_unc/rec.flux
        months[i+1]= dict(  
            eflux = eflux,
            relfluxunc = relfluxunc,
            pull = (1-eflux6y/eflux)/relfluxunc,
            a = rec.a,
            ts=rec.ts,
            locqual=rec.locqual,
            good= rec.ts>10 and rec.a>0.01 and rec.a<0.25 and rec.locqual<8, 
            )
    assert len(months)>10, 'source %s not in at least 10 monthly models' % sname
    monthly = pd.DataFrame(months).T
    good=monthly.good
    pull_rms = monthly[good].pull.std()
    pull_mean = monthly[good].pull.mean()                  
    ngood=sum(good)
    return dict( monthly=monthly, 
                 nmonths=len(monthly), 
                 ngood=ngood,
                 pull_rms=pull_rms,
                 mean=pull_mean,
                 ts = rec6y.ts,
                 eflux = eflux6y,
                 )


    
    
def plots(self):
    dfall = self.df
    cats =set(dfall.acat)
    agn = (dfall.adeltats<9) & [(t in 'agn bllac qso crates bzcat cgrabs'.split()) for t in dfall.acat]
    psr_lat = (dfall.adeltats<9) & (dfall.acat=='pulsar_lat')
    psr = (dfall.adeltats<9) & (dfall.acat=='pulsar_big')
    dfall['singlat'] =np.sin(np.radians(np.asarray(dfall.glat,float)))
    hist_kw = dict(lw=2, histtype='step', log=True)


    fig, axx = plt.subplots(1,2,figsize=( 10,6))
    def glat_fig(ax):
        bins = np.linspace(-1,1,41)
        ax.hist(dfall.singlat, bins, label='all', **hist_kw)
        ax.hist(dfall.singlat[dfall.locqual<5], bins,label='good loc', **hist_kw)
        ax.hist(dfall.singlat[agn], bins, label='AGN assoc', **hist_kw)
        ax.grid(True, alpha=0.5)
        legend(ax)
        plt.setp(ax, title='sin(glat) for transients', xlabel='sin(glat)');
    def dts_fig(ax):
        bins = np.linspace(-5,10,31) 
        ax.hist(dfall.adeltats.clip(-5,10),bins, label='associated', **hist_kw);
        ax.hist(dfall.adeltats[agn].clip(-5,10),bins, label='agn',color='red', **hist_kw);
        ax.hist(dfall.adeltats[psr_lat].clip(-5,10),bins, label='lat_psr', color='green',**hist_kw);
        ax.grid(True, alpha=0.5);
        legend(ax)
        plt.setp(ax, xlim=(-5,10), xlabel='Delta TS', ylim=(1,None), title='association TS')
    def locqual(ax):
        dfall['lq'] = np.asarray(dfall.locqual, float)
        bins = np.linspace(0,8,25) 
        ax.hist(dfall.lq[self.hilat].clip(0,8),bins, label='|b|>10', **hist_kw);
        ax.hist(dfall.lq[agn].clip(0,8),bins, label='agn',color='red', **hist_kw);
        ax.hist(dfall.lq[np.abs(dfall.glat<5)].clip(0,8),bins, label='|b|<5', color='green',**hist_kw);
        ax.grid(True, alpha=0.5);
        legend(ax)
        plt.setp(ax, xlim=(0,8), xlabel='localization quality', ylim=(1,None), title='loc qual')
    for f,ax in zip([dts_fig,], axx.flatten()): f(ax)
    return fig
        
def plots2(self):
    fig, axx = plt.subplots(2,3, figsize=(15,12))
    
    df = self.df[self.df.ts>10]
    df['lq'] = np.asarray(df.locqual, float)
    unassoc = np.logical_not(self.assoc)
    def make_hist(ax, var, bins, xlabel):
        
        hist_kw=dict(lw=2, histtype='step', log=True)
        #df['locqual'] = np.asarray(df.locqual, float)
        ax.hist(var, bins, label='all', **hist_kw);
        ax.hist(var[self.hilat], bins, color='orange', label='hilat', **hist_kw);
        ax.hist(var[self.assoc], bins, color='green', label='assoc', **hist_kw);
        ax.hist(var[unassoc], bins, color='red', label='unassoc', **hist_kw);
        leg = ax.legend()
        for box in leg.get_patches():
            box._height=0; box._y=0.5
        plt.setp(ax, ylim=(1,None), xlabel=xlabel, xlim=(bins[0], bins[-1]))
        ax.grid(True, alpha=0.5)

    def lq_hist(ax):
        bins=np.linspace(0,10,51)
        make_hist(ax, df.lq, bins, xlabel='localization quality')
        
    def pindex_hist(ax):
        bins=np.linspace(0,5,51)
        make_hist(ax, df.pindex.clip(0,5), bins, 'photon index')
    def eflux_hist(ax):
        bins = np.linspace(0,10,51)
        make_hist(ax, df.eflux, bins, 'energy flux (eV)')
    def sin_glat(ax):
        df['singlat'] = np.sin(np.radians(np.asarray(df.glat,float)))
        bins = np.linspace(-1,1,41)
        make_hist(ax, df.singlat, bins, 'sin(glat)')
    def ts(ax):
        bins = np.logspace(1,2,26)
        make_hist(ax, df.ts, bins, 'TS')
        plt.setp(ax, xscale='log')
        
        
    for f,ax in zip([lq_hist, pindex_hist, eflux_hist, sin_glat, ts,], axx.flatten()): f(ax)

def monthly(self):
    df=self.df
    fig, ax = plt.subplots(1,1,figsize=(8,5))
    hist_kw=dict(lw=2, histtype='step', log=True)
    bins=np.linspace(0.5,72.5,73)
    ax.hist(df.month, bins, label='all', **hist_kw);
    ax.hist(df.month[self.hilat], bins, label='hilat',color='orange', **hist_kw);
    ax.hist(df.month[self.assoc], bins, label='assoc',color='red', **hist_kw);
    plt.setp(ax, xlim=(0.5,72.5), xlabel='month', ylim=(10,None), title='Montly totals');
    legend(ax, loc='lower left');
    ax.grid(True, alpha=0.5)
    return fig

def pair_correlations(self, df):
    hpmask = [self.df.roi==n for n in range(1728)]
    def pair_diffs(k, df):
        r = df[hpmask[k]]
        rnb = pd.concat([df[hpmask[j]] for j in neighbors(k)])
        deltas = []
        inside = r.skydir
        outside = rnb.skydir
        for i,a in enumerate(inside):
            for j,b in enumerate(inside): #[i+1:]:
                if i==j: continue
                deltas.append(a.difference(b))
            for b in outside:
                deltas.append(a.difference(b))
        return np.degrees(deltas)

    return np.hstack([pair_diffs(k,df) for k in range(1728)])    

def cumulative_ts(self):
    ny = (len(self.monthinfo)+1)/6
    fig, axx = plt.subplots(ny,6, sharex=True, sharey=True, figsize=(16,5*ny))
    plt.subplots_adjust(hspace=0.3)
    for ax, info, month,  in zip(axx.flatten(), self.monthinfo, self.monthlist):
        df=info.df
        df['prefix'] = [n[:2] for n in df.index]
        df['fixed'] = [not x for x in df.transient]
        fig=info.cumulative_ts(ax=ax, check_localized=False, tscut=[],
                    ts = df.ts[df.fixed], label='6-year sources',
                    other_ts=[df.ts[df.prefix=='PG'],
                               df.ts[(df.prefix=='Sh') | (df.prefix=='S ')|(df.prefix=='TS')]] ,
                    other_label=['PGWave', 'TSmap'], 
                    legend=False,
                    );
        plt.setp(ax, title=month, ylim=(1,2000))
