"""
Comparison with a gtlike model

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/gtlikecomparison.py,v 1.11 2013/12/07 01:21:32 burnett Exp $

"""

import os, pickle, glob
import numpy as np
import pylab as plt
import matplotlib.gridspec as gridspec
import pandas as pd

from uw.utilities import makepivot
from uw.like import Models
from . import sourcecomparison
from . analysis_base import html_table, FloatFormat

class GtlikeComparison(sourcecomparison.SourceComparison):
    """Results of comparison with glike 
    Gtlike version:  %(catname)s Compare the two lists of sources and spectral parameters, assuming that the skymodel names 
    correspond to the "NickName" field in the gtlike-generated FITS file.
    Assume that a special run has been made to evaluate the likelihoods for the gtlike models.
    """
    def setup(self, catpat='gll_psc4year*.fit', **kw):
        gllcats = sorted(glob.glob(os.path.expandvars(os.path.join('$FERMI','catalog', catpat))))
        assert len(gllcats)>0, 'No gtlike catalogs found'
        cat = gllcats[-1]
        catname= os.path.split(cat)[-1]
        super(GtlikeComparison, self).setup(cat=cat, catname=catname ) #cat.split('_')[-1].split('.')[0], **kw)
        cat_version = (catname.split('_')[2]).lower()
        assert cat_version[0]=='v', 'expected version in 3rd token of %s' %catname
        self.plotfolder = 'comparison_%s' % cat_version
        
        # make copy of the df  index with no blanks in names, for comparison, combination
        cname = [n.replace(' ','') for n in self.df.index]
        gtname_dict= dict(zip(cname,self.df.index)) # look up original name, w/ blanks
        self.dfx = self.df.copy()
        self.dfx.index=cname
 
        #make a data frame from the analysis
        ff, tt = self.load_pickles('gtlike/models')
        
        gtcat = dict()
        for t in tt:
            for s in t:
                gtcat[s['name']] = s
                #gtcat[s[0]]=s[1]
        self.gdf = pd.DataFrame(gtcat).T
        self.gdf.index = [n.replace(' ','') for n in self.gdf.name]
        print 'loaded analysis of gtlike fit models, found %d sources' % len(self.gdf)
        # copy info from gtlike to columns of corresponding skymodel entries
        for col in self.gdf.columns:
            self.dfx[col] = self.gdf[col]
        df = self.dfx
        self.delta = df.ts-df.other_ts
        df['name3'] = self.cat.name3
        df['plane']= np.abs(df.glat)<5
        df['ts_gtlike']= self.cat.ts ## should be the TS from the catalog
        df['ts_delta'] = self.delta
        df['ts_gt'] = df.other_ts
        df['index_gt'] = [m['Index'] if isinstance(m, Models.Model) else np.nan for m in self.dfx.othermodel]
        df['ts_pt'] = df.ts
        df['no_gtlike'] = pd.isnull(df.ts_gtlike) * (~np.isinf(df.ts_gtlike))
        
        df['pivot_gt'] = self.cat['pivot']
        df['flux_gt'] = self.cat.flux
        df['modflux'] =[m(e) for (m, e ) in zip(df.model, df.pivot_gt)]
        df['flux_ratio']= df.modflux/df.flux_gt

        
        # finally, restore the blanks in the index for self.dfx
        df.index = [gtname_dict[n] for n in df.index]
        
    def check_contents(self):
        """Contents of the two lists
        %(content_differences)s
        """
        df = self.dfx
        gtnames = set(self.cat.index.values)
        ptnames = set(map(lambda s:s.replace(' ',''), df.index.values))
        added, lost = list(gtnames.difference(ptnames)), sorted(list(ptnames.difference(gtnames)))

        s  = html_table(self.cat.ix[added]['ra dec ts'.split()], 
             dict(ts='TS, gtlike TS', ra='RA,', dec='Dec,'), 
             heading = '\n<h4>Sources in gtlike not in pointlike</h4>',
             float_format=FloatFormat(2), maxlines=50)
             
        cut50 = (df.no_gtlike)*((df.ts>50)+df.psr*(df.ts>10)) 
        missing50 = df.ix[np.array(cut50,bool)]['ra dec ts fitqual locqual roiname'.split()]
        
        s += html_table(missing50,
            dict(ts='TS, pointlike TS', ra='RA,', dec='Dec,', fitqual=',Fit quality,', 
                locqual=',Localization quality; this is NaN for extended sources'),
            heading='\n<h4>Sources with pointlike TS>50 or LAT pulsars and TS>10, not in gtlike</h4>', 
            float_format=FloatFormat(2), maxlines=50)
        
        s +=  html_table(df.ix[np.isinf(df.ts_gtlike.values)] ['ra dec ts fitqual locqual roiname'.split()],
            dict(ts='TS, pointlike TS', ra='RA,', dec='Dec,', fitqual=',Fit quality,', 
                locqual=',Localization quality; this is NaN for extended sources'),
            heading = '\n<h4>Sources present, but not fit by gtlike</h4>',
            float_format=FloatFormat(2), maxlines=50)

        self.content_differences = s
        
    def compare_fits(self):
        """ Compare spectral quantities for sources common to both models
        <br><b>Top Left: </b> Gtlike TS distribution.
        <br><b>Top center:</b> Gtlike TS vs. pointlike TS, showing the high latitude subset.
        <br><b>Top right: </b> Comparison of the pivot energies.
        <br><b>Bottom: </b>Ratio of pointlike differential flux at gtlike pivot energy, vs ptlike flux.
        
        """
        df = self.dfx
        cat = self.cat
        lowlat = np.abs(df.glat)<5
        psr = df.modelname=='PLSuperExpCutoff'

        def plot1(ax):
            ax.hist(self.cat.ts.clip(0,1000), np.logspace(1,3,41))
            plt.setp(ax, xscale='log', ylim=(0,200), xlabel='gtlike TS')
            ax.grid()
            ax.axvline(25, color='g')
        def plot2(ax, lim=(10,1e4)):
            df['ts_gtlike'] = self.cat.ts
            ax.loglog(df.ts_gtlike, df.ts, '.')
            ax.plot(df.ts_gtlike[lowlat], df.ts[lowlat], '.r', label='|b|<5')
            ax.plot(lim, lim, '--g')
            plt.setp(ax, xlabel='gtlike TS', ylabel='pointlike TS', xlim=lim,ylim=lim)
            ax.grid()
            ax.legend(loc='upper left', prop=dict(size=10))
            ax.axvline(25, color='g')
        def plot_pivot(ax, xylim = (100,3e4)):
            self.dfx['pivot_gt'] = self.cat['pivot']
            ax.loglog(self.dfx.pivot_gt,      self.dfx.e0, '.')
            plt.setp(ax, xlim=xylim, ylim=xylim, xlabel='gtlike pivot', ylabel='pointlike pivot')
            ax.plot(xylim, xylim, '--g')
            ax.grid();
            
        def plot_flux_ratio(ax):
            df['pivot_gt'] = cat['pivot']
            df['flux_gt'] = cat.flux
            df['modflux'] =[m(e) for (m, e ) in zip(df.model, df.pivot_gt)]
            ax.loglog(df.modflux, df.modflux/df.flux_gt, '.')
            ax.loglog(df.modflux[lowlat], (df.modflux/df.flux_gt)[lowlat], '.r', label='|b|<5')
            plt.setp(ax, ylim=(0.5,2.0), xlim=(1e-15, 2e-9), xlabel='ptlike flux at gtlike pivot', ylabel='pt/gt flux ratio')
            ax.axhline(1, color='gray')
            for u in (0.95,1.05): 
                ax.axhline(u, ls='--', color='orange')
            ax.set_yticks((0.5, 0.8, 1.0, 1.25, 2.0))
            ax.set_yticklabels('0.5 0.8 1.0 1.25 2.0'.split())
            ax.legend()
            ax.grid()

        gs = gridspec.GridSpec(2,3)
        gs.update(wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(12,8))
        axx = [plt.subplot(gs[0,col]) for col in range(3)]
        axx.append(plt.subplot(gs[1,:]) )
        
        toplot = [plot1, plot2, plot_pivot, plot_flux_ratio]
        for f, ax in zip(toplot, axx):
            f(ax)
        return fig

    def galactic_distortion(self):
        """Galactic plane distortion
        Check for distortion of the spectral fits to moderate strength sources near the galactic plane.
        Sources were selected with flux to the range $10^{-13}$ to $10^{-11}$.
        </br>
        <b>Top row:</b> Histograms of the ratio of fluxes, spectral indeces, and the correlation.
        <br>
        <b>Bottom row:</b> Histogram of the TS difference, and that compared with the flux ratio for galactic sources.
        """
        df = self.dfx
        cut1 = (df.modflux>1e-13) & (df.modflux<1e-11) 
        cut2 = cut1 & (np.abs(df.glat)<10)
        cut2_label='|b|<10'
        
        def plot1(ax): 
            ax.hist((df.modflux/df.flux_gt)[cut1], np.linspace(0,2))
            ax.hist((df.modflux/df.flux_gt)[cut2], np.linspace(0,2), color='orange',label=cut2_label)
            plt.setp(ax, xlabel='pt/gt flux ratio')
            ax.legend(prop=dict(size=10)); ax.grid()
        
        def plot2(ax, space=np.linspace(0.75, 1.25)): 
            ax.hist((df.pindex/df.index_gt)[cut1], space)
            ax.hist((df.pindex/df.index_gt)[cut2], space, color='orange', label=cut2_label)
            plt.setp(ax, xlabel='pt/gt index ratio')
            ax.legend(prop=dict(size=10)); ax.grid()
        
        def plot3(ax):
            ax.plot((df.modflux/df.flux_gt)[cut2], (df.pindex/df.index_gt)[cut2], '.', color='blue')
            plt.setp(ax, xlim=(0,2), xlabel='flux ratio', ylim=(0.75, 1.25), ylabel='index ratio')
            ax.grid()
        
        tdlim=(-50,50); tsdiff_label='TS_pt-TS_gtlike'
        tsdiff =(df.ts_pt-df.ts_gtlike).clip(*tdlim) 
        tdspace=np.linspace(*tdlim)

        def plot4(ax):
            ax.hist(tsdiff[cut1], tdspace)
            ax.hist(tsdiff[cut2], tdspace, color='orange', label=cut2_label)
            ax.grid(); ax.legend(prop=dict(size=10))
            plt.setp(ax, xlabel=tsdiff_label, xlim=tdlim)
        
        def plot5(ax):
            ax.plot( (df.modflux/df.flux_gt)[cut2], tsdiff[cut2], '.')
            plt.setp(ax, xlim=(0.,2.0), xlabel='flux_ratio', ylim=tdlim, ylabel=tsdiff_label)
            ax.grid()
        
        fig, axx = plt.subplots(2,3, figsize=(12,8))
        plt.subplots_adjust(wspace=0.3)
        toplot = [plot1, plot2, plot3, plot4, plot5, None]
        for f, ax in zip(toplot, axx.flatten()):
            if f is not None: 
                f(ax)
            else: ax.set_visible(False)
        return fig
        
    def delta_ts(self, dmax=10, dmin=-1, pivotit=True):
        """ Delta TS
        Plots of the TS for the gtlike fit spectra determined with the pointlike analysis, 
        compared with the pointlike value.<br>
        Mismatches: %(over_ts)d with gtlike worse by %(dmax)d; %(under_ts)d with pointlike worse by %(dmin)d.<br>
        <br><b>Top left</b>": Scatter plot of $\Delta$ TS with the pointlike TS.
        <br><b>Top middle</b>: Histogram of $\Delta$ TS.
        <br><b>Top right</b>: Distribution with respect to galactic latitude, with the 
        subset with discrepancies shown.
        <br><b>Bottom left</b>
        <br><b>Bottom middle</b> For strong sources, this checks the possibilty that the origin of
        discrepancies for the gtlike fits applied to the pointlike data is a consequence of a exposure 
        bias. If so, there would be a correlation of the ratio of fluxes at the gtlike pivot 
        energy with $\Delta$ TS.
        <br>%(mismatch_table)s
        <br>%(pivot_info)s
       """
        df = self.dfx
        delta = df.ts_delta
        mismatch = (delta>dmax)+(delta<dmin)
        self.dmax = dmax; self.dmin=dmin
        df['logflux'] = np.log10(np.asarray(df.flux,float))
        fixme = df[mismatch]['name ts ts_gtlike glat plane fitqual ts_delta ts_gt ts_pt logflux flux_ratio freebits beta roiname'.split()].sort_index(by='roiname')
        fixme.index = fixme.name
        fixme.index.name='name'
        self.mismatch_table=html_table( fixme, columns={}, name=self.plotfolder+'/mismatch', 
            heading='Table of poor matches with delta_ts<%d or >%d'%(dmin,dmax), href=True,
            float_format=FloatFormat(2)) 
        fixme.to_csv('gtlike_mismatch.csv')
        print 'wrote %d entries to gtlike_mismatch.csv' % len(fixme)
        version = os.path.split(os.getcwd())[-1]
        if pivotit:
            try:
                pc=makepivot.MakeCollection('gtlike mismatch %s/%s'% (version, self.catname), 'gtlike/sed', 'gtlike_mismatch.csv', refresh=True)
                self.pivot_info="""<p> These can be examined with a 
            <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
            which requires Silverlight. """% pc.cId
            except Exception, msg:
                self.pivot_info = '<p>No pivot output; job failed %s' %msg
        else:
            self.pivot_info = '<p>(no pivot)'
        delta = self.delta
        x = np.array(delta, float).clip(dmin,dmax) # avoid histogram problem
        cut = (~np.isnan(x))#*(df.ts>10)
        hilat = np.array(cut*(np.abs(df.glat)<5),bool)
        self.under_ts = sum((delta<dmin)*cut)
        self.over_ts  = sum((delta>dmax)*cut)
        print 'under, over delta_ts: %d, %d' % (self.under_ts, self.over_ts)
        def plot1(ax):
            ax.plot(df.ts_pt[cut], delta[cut].clip(dmin,dmax), '.')
            ax.plot(df.ts_pt[hilat], delta[hilat].clip(dmin,dmax), '.r')
            plt.setp(ax, ylim=(dmin-1,dmax+1), xscale='log', xlim=(10,1e4), xlabel='pointlike TS', ylabel='TS diff')
            ax.grid(); #ax.legend(prop = dict(size=10))
        def plot2(ax):
            bins = np.linspace(dmin,dmax,4*(dmax-dmin)+1)
            ax.hist( x[cut], bins)
            ax.hist(x[hilat], bins, color='red', label='|b|<5')
            plt.setp(ax, xlabel='TS diff', xlim=(dmin, dmax))
            ax.grid(); ax.legend(prop = dict(size=10))
        def plot3(ax):
            bins = np.linspace(-1,1,51)
            singlat = np.sin(np.radians(np.array(df.glat, float)))
            ax.hist( singlat, bins)
            ax.hist( singlat[delta>10], bins, color='r', label='delta_ts>10')
            plt.setp(ax, xlabel='sin(glat)')
            ax.grid()
            ax.legend(prop = dict(size=10))

        def plot4(ax, dmax=25,cut=(~np.isnan(x)) & (df.ts>1e3)):
            bins = np.linspace(dmin,dmax,4*(dmax-dmin)+1)
            x = np.array(delta, float).clip(dmin,dmax) # avoid histogram problem
            ax.hist( x[cut], bins)
            ax.hist(x[hilat*cut], bins, color='red', label='|b|<5')
            plt.setp(ax, xlabel='TS diff', xlim=(dmin, dmax))
            ax.grid(); ax.legend(prop = dict(size=10))
        def plot5(ax, dmax=25):
            cut= df.logflux>-11
            ax.plot(df.ts_delta[cut], df.flux_ratio[cut],'.', label='flux>1e-11')
            plt.setp(ax, xlim=(-5, 25), ylim=(0.8, 1.2),xlabel='delta TS',ylabel='flux ratio')
            ax.legend(prop = dict(size=10))
            ax.axhline(1.0,color='gray')
        
        fig, ax = plt.subplots(2,3, figsize=(14,10))
        plt.subplots_adjust(left=0.1)
        for f, ax in zip([plot1,plot2,plot3, plot4, plot5,None], ax.flatten()): 
            if f is not None: f(ax)
            else: ax.set_visible(False)
        return fig
    
    def missing(self, tsmax=50):
        """ Sources in skymodel missing from gtlike       
        Examine pointlike TS and positions for sources in the model that were rejected by the gtlike analysis, mostly by the gtlike TS>25 requirement.
        """
        df = self.dfx
        fig, ax = plt.subplots(1,2, figsize=(10,4))
        plt.subplots_adjust(left=0.1)
        ts = df.ts[df.no_gtlike & (df.ts>10)]
        ts25=df.ts[df.no_gtlike & (df.ts>25)]
        def plot1(ax):
            ax.hist(ts.clip(0,tsmax), np.linspace(0,tsmax,51))
            plt.setp(ax, xscale='linear', xlim=(0,tsmax), xlabel='pointlike TS')
            ax.axvline(25, color='k')
            ax.grid()
        def plot2(ax):
            self.skyplot( ts25, ax=ax, vmin=25, vmax=tsmax, cbtext='pointlike TS')
        for f, ax in zip([plot1,plot2,], ax): f(ax)
        return fig

    def all_plots(self):
        self.runfigures([self.check_contents, self.missing, self.compare_fits,  self.delta_ts, self.galactic_distortion  ])