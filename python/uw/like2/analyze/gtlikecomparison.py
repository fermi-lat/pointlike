"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/gtlikecomparison.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import os, pickle, glob
import numpy as np
import pylab as plt
import pandas as pd

from uw.utilities import makepivot
from . import sourcecomparison

class GtlikeComparison(sourcecomparison.SourceComparison):
    """ Comparison plots with a gtlike analysis
    """
    def setup(self, catpat='gll_psc4year*.fit', **kw):
        gllcats = sorted(glob.glob(os.path.expandvars(os.path.join('$FERMI','catalog', catpat))))
        assert len(gllcats)>0, 'No gtlike catalogs found'
        cat = gllcats[-1]
        super(GtlikeComparison, self).setup(cat=cat, catname=cat.split('_')[-1].split('.')[0], **kw)
        self.plotfolder = 'comparison_%s' % self.catname
        # make copy of the df with no blanks in names, for comparison, combination
        cname = [n.replace(' ','') for n in self.df.index]
        self.dfx = self.df.copy()
        self.dfx.index=cname
 
        #make a data frame from the analysis
        ff = glob.glob('gtlike/models/*.pickle')
        assert len(ff)>0, 'Expected to find pickle files in gtlike/models'
        tt = [pickle.load(open(f)) for f in ff]
        
        gtcat = dict()
        for t in tt:
            for s in t:
                gtcat[s['name']] = s
                #gtcat[s[0]]=s[1]
        self.gdf = pd.DataFrame(gtcat).T
        self.gdf.index = [n.replace(' ','') for n in self.gdf.name]
        print 'loaded analysis of gtlike fit models, found %d sources' % len(self.gdf)
        for col in self.gdf.columns:
            self.dfx[col] = self.gdf[col]
            
        df = self.dfx
        self.delta = df.ts-df.other_ts

        df['plane']= np.abs(df.glat)<5
        df['ts_gtlike']= self.cat.ts ## should be the TS from the catalog
        df['ts_delta'] = self.delta
        df['ts_gt'] = df.other_ts
        df['ts_pt'] = df.ts
        fixme = df[((self.delta>25)+(self.delta<-1))*(df.ts>10)]['name ts ts_gtlike glat plane fitqual ts_delta ts_gt ts_pt freebits beta roiname'.split()].sort_index(by='roiname')
        fixme.index = fixme.name
        fixme.index.name='name'
        fixme.to_csv('gtlike_mismatch.csv')
        print 'wrote %d entries to gtlike_mismatch.csv' % len(fixme)
        version = os.path.split(os.getcwd())[-1]
        pc=makepivot.MakeCollection('gtlike mismatch %s/v4'%version, 'gtlike/sed', 'gtlike_mismatch.csv')
        self.pivot_id=pc.cId

    def compare(self):
        """ Compare spectral quantities for sources common to both models
        
        """
        df = self.dfx
        fig, ax = plt.subplots(1,3, figsize=(14,4))
        plt.subplots_adjust(wspace=0.3)
        def plot1(ax):
            ax.hist(self.cat.ts.clip(0,1000), np.logspace(1,3,41))
            plt.setp(ax, xscale='log', ylim=(0,200), xlabel='gtlike TS')
            ax.grid()
            ax.axvline(25, color='g')
        def plot2(ax, lim=(10,1e4)):
            df['ts_gtlike'] = self.cat.ts
            ax.loglog(df.ts_gtlike, df.ts, '.')
            ax.plot(lim, lim, '--r')
            plt.setp(ax, xlabel='gtlike TS', ylabel='pointlike TS', xlim=lim,ylim=lim)
            ax.grid()
            ax.axvline(25, color='g')
        def plot_pivot(ax, xylim = (100,3e4)):
            psr = self.dfx.modelname=='PLSuperExpCutoff'
            self.dfx['pivot_gt'] = self.cat['pivot']
            ax.loglog(self.dfx.pivot_gt, self.dfx.e0, '.')
            ax.loglog(self.dfx.pivot_gt[psr], self.dfx.pivot_energy[psr], 'or', label='pulsars')
            plt.setp(ax, xlim=xylim, ylim=xylim, xlabel='gtlike pivot', ylabel='pointlike pivot')
            ax.plot(xylim, xylim, '--r')
            ax.grid(); ax.legend(loc='upper left', prop=dict(size=10))

        for f, ax in zip([plot1,plot2,plot_pivot,], ax): f(ax)
        return fig

    def delta_ts(self, dmax=25, dmin=-1):
        """ Delta TS
        Plots of the TS for the gtlike fit spectra detrmined with the pointlike analysis, compared with the pointlike value.<br>
        Outliers: %(over_ts)d with gtlike worse by 25; %(under_ts)d with pointlike worse by 1.
        <p> These can be examined with a 
        <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%(pivot_id)d">Pivot browser</a>,
        which requires Silverlight.  
       """
        df = self.dfx
        delta = self.delta
        x = np.array(delta, float).clip(dmin,dmax) # avoid histogram problem
        cut = (~np.isnan(x))*(df.ts>10)
        hilat = cut*(np.abs(df.glat)<5)
        self.under_ts = sum((delta<dmin)*cut)
        self.over_ts  = sum((delta>dmax)*cut)
        print 'under, over delta_ts: %d, %d' % (self.under_ts, self.over_ts)
        fig, ax = plt.subplots(1,2, figsize=(11,4))
        def plot1(ax):
            ax.plot(df.ts_pt[cut], delta[cut].clip(dmin,dmax), '.')
            ax.plot(df.ts_pt[hilat], delta[hilat].clip(dmin,dmax), 'or')
            plt.setp(ax, ylim=(dmin-1,dmax+1), xscale='log', xlim=(10,1e4), xlabel='pointlike TS', ylabel='TS diff')
            ax.grid(); #ax.legend(prop = dict(size=10))
        def plot2(ax):
            bins = np.linspace(dmin,dmax,1*(dmax-dmin+1))
            ax.hist( x[cut], bins)
            ax.hist(x[hilat], bins, color='red', label='|b|<5')
            plt.setp(ax, xlabel='TS diff', xlim=(dmin, dmax))
            ax.grid(); ax.legend(prop = dict(size=10))
        for f, ax in zip([plot1,plot2,], ax): f(ax)
        return fig
    
    def missing(self):
        """ Sources in skymodel not fit by gtlike        
        Examine sources in the model that were rejected by the gtlike analysis, mostly by the (gtlike) TS>25 requirement.
        """
        df = self.dfx
        fig, ax = plt.subplots(1,2, figsize=(10,4))
        ts = df.ts[pd.isnull(df.ts_gt)*(df.ts>10)]
        def plot1(ax):
            ax.hist(ts.clip(0,100), np.linspace(0,100,51))
            plt.setp(ax, xscale='linear', xlim=(0,100), xlabel='pointlike TS')
            ax.axvline(25, color='k')
            ax.grid()
        def plot2(ax):
            self.skyplot( ts, ax=ax, vmin=10, vmax=100, cbtext='pointlike TS')
        for f, ax in zip([plot1,plot2,], ax): f(ax)
        return fig

    def all_plots(self):
        """Results of comparison with glike version %(catname)s 
        """
        self.runfigures([self.compare, self.missing, self.delta_ts,  ])