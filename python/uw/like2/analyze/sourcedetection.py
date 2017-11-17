"""
Seed finding plots

$Header$

"""

import os, pickle
import numpy as np
import pylab as plt
import pandas as pd

from . import analysis_base, sourceinfo
from ..import maps
from ..plotting import sed

# Can't do this now due to loop
# from ..pipeline import stagedict

class SourceDetection(sourceinfo.SourceInfo):
    """ Source detection analysis

        <p>Source detection is performed ...

    """
    def setup(self, **kwargs):
        #keys =stagedict.stagenames['sourcefinding']['pars']['table_keys']
        plt.style.use('seaborn-bright')
        modelname = os.getcwd().split('/')[-1]
        self.seedfile ='seeds/seeds_all.csv' 
        self.seeddf = pd.read_csv(self.seedfile, index_col=0)
        self.n_seeds = len(self.seeddf)
        super(SourceDetection, self).setup(**kwargs)
        self.plotfolder='sourcedetection'
        # get template info from maps
        self.keys=keys = 'hard soft ts tsp'.split()
        self.alias_dict=dict(hard='hard', soft='soft', ts='flat', tsp='pulsar')
        alias_names = map(lambda k: self.alias_dict[k], keys)        
        modelnames = [maps.table_info[key][1]['model'] for key in keys]; modelnames
        models = [eval('maps.'+modelname) for modelname in modelnames]
        class Source(object):
            def __init__(self, name, key, model):
                self.name=name; self.key=key; self.model= model
        self.templates = map(Source, alias_names, keys, models) 

        # add detection information to the seed dataframe
        df = self.df
        sdf = self.seeddf
        sdf['ok'] = [name in df.index for name in sdf.index]
        sdf['pindex'] = df.pindex
        sdf['ts_fit'] = df.ts
        sdf['associations']=df.associations
        aprob=[]; acat=[]
        for a in sdf.associations:
            if a is None or pd.isnull(a) :
                aprob.append(0)
                acat.append('')
            else:
                aprob.append(a['prob'][0])
                acat.append(a['cat'][0])
        sdf['aprob']=aprob
        sdf['acat'] = acat

    def template_info(self):
        """ Template information

        """
        fig,ax= plt.subplots(figsize=(4,4))
        for src in self.templates:
            sed.Plot(src)(axes=ax,name='', butterfly=False, fit_kwargs=dict(label=src.name, lw=2))
        ax.legend(loc='upper left', prop=dict(size=12, family='monospace'), );
        ax.set_ylim(1e-2, 20);
        ax.set_xlim(100,1e5)
        ax.set_title('Template SED plots')
        return fig

    def seed_plots(self,  bcut=5, subset=None, title=None):
        """ Seed plots

        Results of cluster analysis of the residual TS distribution. Analysis of %(n_seeds)d seeds from file 
        <a href="../../%(seedfile)s">%(seedfile)s</a>. 
        <br>Left: size of cluster, in 0.15 degree pixels
        <br>Center: maximum TS in the cluster
        <br>Right: distribution in sin(|b|), showing cut if any.
        """
        z = self.seeddf #self.seeds if subset is None else self.seeds[subset]
        keys = list(set(z.key)); 
        fig,axx= plt.subplots(1,3, figsize=(12,4))
        plt.subplots_adjust(left=0.1)
        bc = np.abs(z.b)<bcut
        histkw=dict(histtype='step', lw=2)
        def all_plot(ax, q, dom, label, log=False):
            for key in keys:
                ax.hist(q[z.key==key].clip(dom[0],dom[-1]),dom, label=self.alias_dict[key], log=log, **histkw)
            plt.setp(ax, xlabel=label, xlim=(None,dom[-1]))
            ax.grid()
            ax.legend(prop=dict(size=10, family='monospace'))
            if log: ax.set_ylim(ymin=0.9)
        all_plot(axx[0], z['size'], np.linspace(0.5,10.5,11), 'cluster size')
        all_plot(axx[1], z.ts.clip(10,25), np.linspace(10,25,16), 'TS', log=True)
        all_plot(axx[2], np.sin(np.radians(z.b)), np.linspace(-1,1,21), 'sin(b)')
        axx[2].axvline(0, color='k')
        fig.suptitle('{} seeds'.format( len(z)) )
        #     if title is None else title)
        fig.set_facecolor('white')
        return fig

    def detection_plots(self):
        """Detection plots

        Comparison of seed values with fits

        """
        sdf = self.seeddf
        ddf =sdf[sdf.ok & (sdf.ts_fit>10)]; 
        fig, axx =plt.subplots(2,2, figsize=(12,12))

        ax =axx[0,0]
        xlim = (1., 3.5)
        hkw = dict(bins=np.linspace(xlim[0],xlim[1],26), histtype='step', lw=2)
        for t in self.templates:
            if t.name=='pulsar': continue
            key =t.key
            ax.hist(sdf.pindex[sdf.ok & (sdf.key==key)].clip(*xlim), label=t.name, **hkw);
            ax.axvline(t.model['Index'], ls='--', color='grey')
        ax.legend(prop=dict(family='monospace'));
        ax.set_xlim(*xlim)
        ax.set_xlabel('fit photon index')
        ax.set_title('Power Law photon index')
        ax.grid(alpha=0.5);

        ax=axx[0,1]
        xlim = (0,20)
        ylim = (0,40)
        for t in self.templates:
            key =t.key
            sel =ddf.key==key
            ax.plot(ddf.ts[sel].clip(*xlim), ddf.ts_fit[sel].clip(*ylim), '.', label=t.name);
        ax.legend(prop=dict(size=12, family='monospace'));
        ax.plot([10,20],[10,20], '--', color='grey')
        plt.setp(ax, xlabel='seed TS', ylabel='fit TS', title="Test Statistic");
    
        ax = axx[1,0]
        sinb = np.sin(np.radians(ddf.b))
        hkw = dict(bins=np.linspace(-1,1,21), histtype='step', lw=2)
        for t in self.templates:
            key =t.key
            sel =ddf.key==key
            ax.hist(sinb[sel], label=t.name, **hkw);
        ax.legend(prop=dict(size=12, family='monospace'));
        ax.set_xlim(-1,1); ax.set_xlabel('sin(b)')
        
        axx[1,1].set_visible(False)
        return fig



    def all_plots(self):
        self.runfigures([self.template_info, 
            self.seed_plots,
            self.detection_plots,
        ])