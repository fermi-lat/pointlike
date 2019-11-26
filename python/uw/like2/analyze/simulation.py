"""
Simulation analysis


"""

import os, glob
# from collections import OrderedDict
# import astropy.io.fits as pyfits
from astropy.table import Table
import numpy as np
import pylab as plt
import pandas as pd
from scipy import stats
from . import (analysis_base, sourceinfo, fermi_catalog)

from skymaps import SkyDir

class Simulation(analysis_base.AnalysisBase): #sourceinfo.SourceInfo):
    """Analysis of this simulation including comparison with FL8Y or 4FGL

    %(info)s

    """
    def setup(self,  **kw):
        
        self.plotfolder='simulation'
        if not os.path.exists('plots/'+self.plotfolder):
            os.mkdir('plots/'+self.plotfolder)
        self.uwmodel= uwmodel = self.config['input_model']['path'].split('/')[-1]

        self.df = pd.read_csv('sources_{}.csv'.format(self.skymodel), index_col=0)
        self.dfa = dfa = pd.read_csv('../{0}/sources_{0}.csv'.format(uwmodel), index_col=0); 
        self.info =  '<br> UW model: {}: {}'.format(uwmodel,len(dfa))

        # this model
        self.sim = sim = self.skymodel 
        self.dfm = dfm = pd.read_csv('sources_{}.csv'.format(sim), index_col=0); 
        self.info += '<br> Simulated model: {}: {} sources'.format(sim, len(dfm))

        def set_glat(df):
            from skymaps import SkyDir
            sds = map(SkyDir, df.ra, df.dec);
            glat = map(lambda x: x.b(), sds)
            df['glat'] = glat
            df['hilat'] = (df.glat<-5) | (df.glat>5)
        set_glat(dfa) 
        set_glat(dfm)

        # tag MC seed sources
        mc_tag = kw.get('mc_tag', sim[-2:])
        def seedcheck(name):
            if name.find( mc_tag )>0:
                return 'MC seed'
            else: return 'fit source'

        grouped =dfm.groupby([seedcheck])
        try:
            self.mc_seed= grouped.get_group('MC seed')
            self.info += ' including {} spurious from seeds'.format(len(self.mc_seed))
        except:
            print 'No seeded sources found: searched for names containing text "{}"'.format(mc_tag)
            self.info += ' (no seeds found)'

        
        # load gll version of FL8Y or 4FGL
 
        fl8y_file = self.config.get('gllcat', None)
        self.cat_name = self.config.get('catname', 'FL8Y')
        assert fl8y_file is not None, '{} "gll" file not specified'.format(self.cat_name)
        if not fl8y_file.startswith('/'):
            fl8y_file = os.path.expandvars('$FERMI/catalog/'+fl8y_file)

        gtest = sorted(glob.glob(fl8y_file))
        assert len(gtest)>0, "no file found, looked for {}".format(fl8y_file) 
        file = gtest[0]

        df_fl8y= fermi_catalog.GLL_PSC2(file).df #Table.read(file, hdu=1).to_pandas()
        self.info += '<br> {} file: {} w/ {} sources'.format(self.cat_name, fl8y_file, len(df_fl8y))
        # df_fl8y.index=df_fl8y.NickName
        # del df_fl8y['NickName']
        dfa['in_fl8y']=[name in df_fl8y.index for name in dfa.index]
        # self.startlog()
        print 'Info:' + self.info.replace('<br>', '\n\t')
        # self.logstream= self.stoplog()


    def scat_ts_pindex(self, scat=True):
        """TS vs. Spectral Index

        """
        dfa, mc_seed = self.dfa, self.mc_seed

        def plotit(df, color, ax=None, title=None):
            if ax is None:
                fig, ax = plt.subplots(figsize=(8,8))
            else: fig=ax.figure
            lolat = np.abs(df.glat)<5
            if scat:
                ax.plot(df.ts, df.pindex.clip(0.5,3.5), '.', color='blue')
                ax.plot(df.ts[lolat], df.pindex[lolat].clip(0.5,3.5), '.', color='red', label='|b|<5')
                ax.axvline(25, ls=':', color='red')
                ax.set(xlim=(10,40), ylim=(0.5,3.5));
                ax.set(xlabel='TS', ylabel=r'$\Gamma$')
                ax.axhline(2.8, ls=':', color='red')
                ax.axhline(1.6, ls=':', color='red')

            else:
                ts_cut = (df.ts>16) & (df.ts<32)
                hkw = dict(bins=np.linspace(0.5,3.5,16), histtype='step', lw=2)
                ax.hist(df.pindex[ts_cut & ~lolat], label='lolat', color='green', **hkw)
                ax.hist(df.pindex[ts_cut & lolat], label='hilat', color='orange', **hkw)
                ax.legend()

                ax.axvline(2.8, ls=':', color='red')
                ax.axvline(1.5, ls=':', color='red')
            if title is not None: ax.set_title(title)

        fig,axx = plt.subplots(2,2, figsize=(10,10), sharex=True, sharey=False)
        axf = axx.flatten()
        plotit(dfa, 'green', ax=axf[0], title=self.uwmodel)
        plotit(dfa.query('in_fl8y==False'), 'grey', ax=axf[2], title=self.uwmodel+' not in {}'.format(self.cat_name))
        plotit(mc_seed, 'orange', ax=axf[3], title='MC seeds')
        plotit(dfa.query('in_fl8y==True'), 'blue', ax=axf[1], title=self.cat_name)

        return fig

    def purity_plots(self, cut=None, tsbins=(16, 60, 19), ylim=(50,105)):
        """Purity

        """
        mc_seed, dfa = self.mc_seed, self.df # compare with this oneself.dfa
    
        fig, axx =plt.subplots(3,2, figsize=(10,10), sharex=True, 
                gridspec_kw=dict(left=0.05, wspace=0.35, hspace=0.01))

        hkw=dict(bins=np.linspace(*tsbins), histtype='step',log=True, lw=2)

        for i, cut, label in zip(range(3), ['pindex<1.6', '3.5>pindex>2.8',None,], 
                    [r'$\mathsf{Hard\ (\Gamma<1.6)}$', r'$\mathsf{Soft\ (\Gamma>2.8)}$','All', ]):
            # get the MC group
            if cut is None:
                mc_ts=mc_seed.ts
                data_ts= self.df.ts #dfa.ts
            else:
                mc_ts=mc_seed.query(cut).ts
                data_ts=dfa.query(cut).ts

            ax=axx[i,0]
            ax.text(0.1,0.9, label, fontsize=16, transform=ax.transAxes)
            tslim=tsbins[:2]
            ax.hist(data_ts.clip(*tslim), label=self.uwmodel,color='green', **hkw);
            ax.hist(mc_ts.clip(*tslim),label='MC seed', color='orange', **hkw);
            ax.legend();
            ax.grid(alpha=0.5)
            ax.set(xlim=tslim, ylim=(0.8,None))
            if i==2: ax.set_xlabel('TS')
            ax.axvline(25, ls=':', color='red')

            ax=axx[i,1]
            bins=hkw['bins']
            mc = np.histogram(mc_ts,bins=bins)[0]
            data = np.array(np.histogram(data_ts, bins=bins)[0], float)
            purity = 1-(mc/data).clip(0,1)
            delta=bins[1]-bins[0]
            t=bins[:-1]+0.5*delta
            yerr=100.*np.sqrt(mc*(mc+data)/data**3)
            ax.errorbar(x=t,y=purity*100, xerr=delta/2,yerr=yerr, fmt='o', marker='o')
            ax.set( ylim=ylim)
            if i==2: ax.set_xlabel('TS')
            ax.set_ylabel(ylabel='purity [%]',fontsize=14)
 
            ax.axvline(25, ls=':', color='red')
            ax.axhline(100, ls=':', color='grey')
            ax.grid(alpha=0.5);
            
        #fig.suptitle(suptitle, fontsize=14)
        return fig

    def latitude_dependence(self, cut='pindex<3.4 and locqual<10 and delta_ts<4'):
        """Latitude dependence

        """
        mc_seed, dfa = self.mc_seed, self.dfa
        dfa_cut = dfa.query(cut)
        mc_cut = mc_seed.query(cut)
        titles = (self.uwmodel, self.uwmodel+ ' sources not in {}'.format(self.cat_name), 'MC seed sources')

        fig,axx = plt.subplots(len(titles),1, figsize=(8,8), sharex=True, sharey=True)
        hkw = dict(bins=np.logspace(1,3,26), log=True, histtype='step', lw=2)
        for df,ax,title in zip((dfa_cut, dfa_cut[~dfa_cut.in_fl8y], mc_cut), axx.flatten(),titles):
            ax.hist(df[df.hilat].ts.clip(10,1e3), label='high lat', **hkw);
            ax.hist(df[~df.hilat].ts.clip(10,1e3), label='low lat', **hkw)
            ax.set(xlabel='TS', xscale='log', ylim=(0.8,None), title=title)
            ax.grid(alpha=0.5)
            ax.legend();  
        return fig

    def spectral_check(self, ):
        """Compare spectral parameters to check uncertainties
        Normalized deviations of %(check_total)s sources with LogParabola fits, comparing the fits to the simulated data with actual values
        that were used to define the model.
        """
        a, b = self.dfa, self.dfm.copy()
        b['ts_a']=a.ts
        b['flux_a'] = a.flux
        b['dflux'] = (b.flux-b.flux_a)/b.flux_unc
        b['eflux100_a'] = a.eflux100
        b['deflux'] = (b.eflux100-b.eflux100_a)/b.eflux100_unc
        b['pindex_a'] = a.pindex
        b['gdelta'] = (b.pindex-b.pindex_a)/b.pindex_unc
        self.dfm = b # since copy

        fig,axx = plt.subplots(1,2, figsize=(10,5), sharey=True)
        hkw = dict(bins=np.linspace(-5,5,51), histtype='step', lw=2, density=True)

        cut =  (b.ts>50) & ~pd.isnull(b.deflux) & ~pd.isnull(b.gdelta) &\
                    (b.modelname=="LogParabola") & (b.pindex<3) & (b.pindex>0.5) &\
                    (b.e0>500) &(b.eflux100_unc>0) &(b.pindex_unc>0)
        self.check_total = sum(cut)
        for ax, title, val in zip(axx.flatten(), ['Energy Flux', 'Spectral index'], [b.deflux, b.gdelta]):    

            df=val[cut]
            ax.hist(df.clip(-5,5), label='mean {:5.2f}\nstd  {:5.2f}'.format(df.mean(),df.std()), **hkw);
            ax.grid(alpha=0.5); 
            x=np.linspace(-4,4)
            ax.plot(x, stats.norm.pdf(x), '--g' );
            ax.set(xlabel='normalized fit deviation', title=title, )
            ax.legend(loc='upper left',prop=dict(family='monospace'))
        fig.suptitle('Normalized devations of fit from model', fontsize=16);

        return fig

    def all_plots(self):
        self.runfigures([ 
            self.spectral_check,
            self.scat_ts_pindex,
            self.purity_plots,
            self.latitude_dependence,
            ])
