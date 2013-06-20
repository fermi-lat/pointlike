"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/seedcheck.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import numpy as np
import pylab as plt
import pandas as pd

from . import sourceinfo

class SeedCheck(sourceinfo.SourceInfo):
    require='seedcheck.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname='seedcheck'
        self.spectral_type='power law'
        self.load()
        
    def load(self):
        files, sources = self.load_pickles(self.seedname)
        sdict={}
        assoc={}
        for source in sources:
            name = source.name
            model = source.model
            pars = np.empty(4); pars.fill(np.nan)
            errs = np.empty(4); errs.fill(-2)
            free = np.zeros(4, bool)
            n = model.len()
            pars[:n] = model.parameters
            free[:n] = model.free
            try:
                diag = np.diag(model.get_cov_matrix())
                errs[:n] = [np.sqrt(x) if x>0 else -1 for x in diag[:n]]
                badfit = np.any(errs[model.free]<=0)
            except Exception, msg:
                print 'fail errors for %s:%s' % (name, msg)
                badfit = True
            has_adict = hasattr(source,'adict') and source.adict is not None
            sdict[name] = dict(
                ra =source.skydir.ra(), dec=source.skydir.dec(),
                ts=source.ts,
                delta_ts=source.ellipse[5] if hasattr(source, 'ellipse') else np.nan,
                r95 = 2.6*source.ellipse[2] if hasattr(source, 'ellipse') else np.nan,
                glat=source.skydir.b(), glon=source.skydir.l(),
                eflux=pars[0] * model.e0**2 *1e6,
                eflux_unc=errs[0] * model.e0**2 *1e6 if errs[0]>0 else np.nan,
                pindex = pars[1],
                pindex_unc = errs[1] if errs[1]>0 else np.nan,
                par2 = pars[2],
                par2_unc = errs[2] if errs[2]>0 else np.nan,
                e0 = model.e0,
                aprob = source.adict['prob'][0] if has_adict else 0,
                index = source.index,
                #gflux  = model.i_flux(), ## photon flux
                )
            assoc[name] = dict(
                acat = source.adict['cat'][0] if has_adict else None,
                aname= source.adict['name'][0] if has_adict else None,
                adelta_ts = source.adict['deltats'][0] if has_adict else None,
                aprob = source.adict['prob'][0] if has_adict else 0.,
                adict = source.adict if has_adict else None,
                )
        self.df = pd.DataFrame(pd.DataFrame(sdict).transpose(), 
            columns='ra dec glat glon ts delta_ts eflux eflux_unc pindex pindex_unc par2 par2_unc e0 r95 aprob index'.split() 
            )
        self.df.index.name='name'
        self.assoc = pd.DataFrame(assoc).transpose()
        self.assoc.index.name = 'name'
        # analyze associations, make summary
        acat=list(self.assoc.ix[self.assoc.aprob>0.8]['acat'].values)
        sa = list(set(acat))
        t = np.zeros(len(sa),int)
        for x in acat:
            t[sa.index(x)]+=1
        self.assoc_sum = zip(sa, t)
        #self.psr = [  s is not None and ('pulsar' in s or 'msp' in s ) for s in t] * associated
        #self.agn = [  s is not None and ('agn' in s or 'crates' in s or 'bzcat' in s or 'bllac' in s or 'qso' in s or 'cgrabs' in s) for s in t] * associated
    
    def select_candidates(self, tsmin):
        t=self.df[(self.df.ts>tsmin)*(self.df.r95<0.6)]['ra dec ts pindex r95 aprob'.split()]#.sort_index(by='ts')
        t['acat']=self.assoc.acat
        t['aname']=self.assoc.aname
        return t
    
    def seed_cumulative_ts(self, cut=None, label='all seeds'):
        """ Cumulative TS distribution for seeds 
        """
        v = self.df.ts
        if cut is not None: v=v[cut]
        fig = self.cumulative_ts(v, check_localized=False, label=label)
        ax = plt.gca()
        plt.setp(ax, ylim=(9,1000), xlim=(9,100))
        leg =ax.legend()
        pbox = leg.get_patches()[0]; pbox._height=0; pbox._y=5
        return fig
        
    def unassoc_seed_cumulative_ts(self):
        """ Cumulative TS distribution for seed sources that are not associated
        """
        return self.seed_cumulative_ts(cut=self.assoc.aprob<0.8, label='unassociated')
        
    def histo(self, ax, v, bins):
        ax.hist(v, bins)
        ax.hist(v[self.df.ts>10], bins, label='TS>10')
        ax.hist(v[self.df.ts>25], bins, label='TS>25')
        ax.legend(prop=dict(size=10))
        ax.grid()
    
    def localization(self):
        """ Localization results
        <br>Left: r95; right; delta TS
        """
        fig, ax = plt.subplots(1,2, figsize=(12,5))
        def r95(ax):
            v = 60.*self.df.r95; bins = np.linspace(0, 25,26)
            self.histo(ax, v[~pd.isnull(v)], bins)
            plt.setp(ax, xlabel='r95 (arcmin)')
        def delta_ts(ax):
            v = np.sqrt(list(self.df.delta_ts)); bins = np.linspace(0,10,26)
            self.histo(ax, v, bins)
            plt.setp(ax, xlabel='sqrt(delta_ts)', xlim=(0,10))

        for f, a in zip((r95, delta_ts), ax.flatten()):
            f(a)
        return fig

    def spectral_parameters(self, ax=None):
        """ Spectral fit parameters
        Flux vs. spectral index for %(spectral_type)s fit
        <br>histograms of sin(glat) and sqrt(delta_ts) for all, TS>10, and TS>25
        """
        fig, ax = plt.subplots(2,2, figsize=(12,12))
        good = self.df.ts>10
        super = self.df.ts>25
        def flux_index(ax):
            for cut, c,label in zip((good, super), ('.b', 'or'), ('TS>10', 'TS>25')):
                ax.plot(self.df.eflux[cut], self.df.pindex[cut], c, label=label)
            ax.grid()
            ax.legend(prop=dict(size=10))
            plt.setp(ax, ylim=(0.5,3.0), xlim=(0.1,10), xscale='log', ylabel='spectral index', xlabel='enegy flux (eV)')
        def singlat(ax):
            v = np.sin(np.radians(list(self.df.glat))); bins=np.linspace(-1,1,26)
            self.histo(ax, v, bins)
            plt.setp(ax, xlabel='sin(glat)')
        def skyplot(ax):
            glon = self.df.glon
            glon[glon>180]-=360
            ax.plot(glon, np.sin(np.radians(list(self.df.glat))), 'o')
            plt.setp(ax, xlim=(180,-180), xlabel='glon', ylabel='sin(glat)')
            #self.skyplot(self, self.df.ts, ax=ax)
        def index_vs_cutoff(ax):
            cutoff = self.df.par2
            for tsmin, marker in zip((10,25), ('.b', 'or')):
                cut = self.df.ts>tsmin
                ax.plot(cutoff[cut], self.df.pindex[cut],  marker, label='TS>%d'%tsmin)
            plt.setp(ax, ylabel='spectral index', xlabel='cutoff', ylim=(0.5,3.0), xlim=(0, 3000))
            ax.grid(); ax.legend(prop=dict(size=10))
        for f, a in zip((flux_index, index_vs_cutoff, singlat, skyplot), ax.flatten()):
            f(a)
            
        return fig

    def locations(self):
        """ Positions
        """
        return self.skyplot(self.df.ts)
        
    def all_plots(self):
        """ Results of analysis of seeds
        %(info)s
        """
        #t=self.df[(self.df.ts>6)*(-self.df.r95.isnull())]['ra dec ts pindex r95 aprob index'.split()]#.sort_index(by='ts')
        #t['acat']=self.assoc.acat
        #t['aname']=self.assoc.aname
        tsmin = 6
        t = self.select_candidates(tsmin)
        good_seeds = 'good_seeds.csv'
        t.to_csv(good_seeds)
        print 'wrote list that succeeded to %s' % good_seeds
        self.info = ' ' #self.df.describe().to_html()
        self.info += '<h3> Selected TS>%.0f and localized</h3>' %tsmin
        self.info += t.to_html()
        self.info += '<h3>Association summary:</h3>' #\n<pre>%s\n</pre>' %self.assoc_sum
        self.info += '<table border="1"><thead><tr><th>Catalog</th><th>Sources</th></tr></thead>\n<tbody>'
        for (c,n) in  self.assoc_sum:
            self.info += '<tr><td>%s</td><td>%5d</td></tr>' % (c,n)
        self.info += '</tbody></table>\n'
        self.runfigures([self.seed_cumulative_ts, self.spectral_parameters, self.localization,])