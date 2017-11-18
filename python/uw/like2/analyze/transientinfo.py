"""
Analysis plots of transient sources

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/transientinfo.py,v 1.6 2016/10/28 20:48:14 burnett Exp $

"""

import os, pickle,  glob
from astropy.io import fits as pyfits
from uw.like2.analyze import (sourceinfo, associations,)
from uw.like2.tools import DateStamp
from pointlike import IntVector
from skymaps import Band, SkyDir
from uw.like2.pub import healpix_map


import numpy as np
import pylab as plt
import pandas as pd

from analysis_base import (html_table, FloatFormat,)
from . import (sourceinfo, localization, associations,)

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
        other_ts=[]
        other_label=[]
        if sum(df.prefix=='PG')>0:
            other_ts.append(df.ts[df.prefix=='PG'])
            other_label.append('PGwave')
        tsmap_sel = (df.prefix=='Sh') | (df.prefix=='S ')|(df.prefix=='TS')
        if sum(tsmap_sel)>0:
            other_ts.append(df.ts[tsmap_sel])
            other_label.append('TSmap')
        fig=super(TransientInfo, self).cumulative_ts(check_localized=False, 
                    ts = df.ts, label='all',
                    other_ts=other_ts,
                    other_label=other_label, 
                    legend=True,
                    );
        plt.setp(fig.axes[0], ylim=(1,1000), xlim=(9,100))
        fig.set_facecolor('white')
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
        fig.set_facecolor('white')
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
        
#B12 = Band(12)
#def neighbors(i):
#    iv = IntVector()
#    n = B12.findNeighbors(int(i),iv)
#    return list(iv)

def legend(ax, **kw):
    leg = ax.legend(**kw)
    for box in leg.get_patches():
        box._height=0; box._y=0.5

def delta_ts_figure(df,  title='', ax=None, xlim=(-4,10), ymax=None, binsize=0.5):
    """ make a Delta TS figure using the dataframe assuming it has ts and deltats or delta_ts columns
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=(5,5))
    else: fig = ax.figure
    deltats= np.array(df.adeltats if 'adeltats' in df else df.deltats if 'deltats' in df else df.delta_ts, float)
    ts = np.array(df.ts, float)
    bins = np.linspace(xlim[0],xlim[-1],(xlim[-1]-xlim[0])/binsize+1 )
    hist_kw=dict(bins=bins, log=True, histtype='step', lw=2)
    ax.hist( deltats.clip(*xlim),  label='all', **hist_kw)
    ax.hist( deltats[ts<25].clip(*xlim), color='orange',
            label='TS<25',  **hist_kw)
    plt.setp(ax, xlabel='delta TS', xlim=xlim, ylim=(0.8,ymax), title=title)
    x = np.linspace(0,10,51)
    ax.plot(x, len(deltats)*binsize*0.5*np.exp(-x/2),ls='--', 
            color='r', label='expected')
    ax.grid(True, alpha=0.5); 
    ax.axvline(0, color='k')
    leg=ax.legend()
    for box in leg.get_patches():
        box._height=0; box._y=0.5
    fig.set_facecolor('white')
    return fig
    
class Analysis(object):
    """ a place to put code to make plots for all months
    """
    def __init__(self,  path=None, prefix='', quiet=False, 
                all_months_model='../P301_6years/uw972'):
        """ If prefix=='o', then it is obssim
        """
        if path is not None:
            os.chdir(os.path.expandvars(path))
        else:
            path=os.getcwd()
        files =  sorted(glob.glob(prefix+'month*/sources.pickle'));
        print 'Found {} "monthxx" folders in folder {}'.format(len(files), path)
        self.monthlist = [f.split('/')[0] for f in files];
        self.monthinfo=monthinfo=[]
        for month in self.monthlist: #[:last+1]:
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
            monthdf = month.df[month.df.transient]["""
                ra dec glat glon ts 
                modelname pindex eflux flux a b ang locqual aprob flux pindex index2 e0 flags
                flux_unc pindex_unc index2_unc cutoff cutoff_unc fitqual 
                acat adeltats aang has_assoc month""".split()]
            if prefix=='o':
                monthdf.index = [s.replace('TSxx','TS%02d' % i) for s in monthdf.index]
            dflist.append(monthdf)
        df = pd.concat(dflist) 
        df['skydir'] = skydirs = map(SkyDir, df.ra, df.dec)
        df['roi'] = map( lambda s: Band(12).index(s), skydirs)
        
        # add the photon density from the 6-year dataset
        filename=os.path.expandvars('$FERMI/skymodels/P301_6years/uw972/hptables_ts_kde_512.fits')
        kde6 = healpix_map.FromFITS(filename,'kde')
        kdepercent = np.percentile(kde6.vec, np.arange(0,100))
        kdevals = map(kde6, skydirs)
        df['kdepercentile'] = kdepercent.searchsorted(kdevals)
        
        # select those with TS>10 and good localization
        good = np.logical_not((df.ts<10) | (df.locqual>8) | (df.a<0.001) | (df.a>0.25))
        self.df = df[good]
        self.hilat = np.abs(self.df.glat)>10
        self.assoc = np.isfinite(self.df.aprob)
        if not quiet:
            print 'good sources: %d/%d' % (len(self.df), len(df))
            print 'High latitude (>10 deg)', sum(self.hilat)
            print 'Associated', sum(self.assoc)
        
        # get the 6-year dataframe for reference, and access to functions in SourceInfo
        self.sinfo = sourceinfo.ExtSourceInfo(all_months_model)
        self.df6y = self.sinfo.df 

    def all_6y_names(self):
        """ return list of all 6-year names found in the monthly models
        """
        return self.df6y[self.df6y.ts>10].index
        
    def psr_names(self):
        """ return list of unique PSR names found in the monthly models
        """
        names=[]
        for month in self.monthinfo:
            df = month.df
            for nm in list(df.index):
                if nm.startswith('PSR'):
                    names.append(nm)

        return list(set(names))

    def monthly_info(self, sname):
        """Return monthly info for sources in all or most months
        """
        s6y = self.df6y
        assert sname in s6y.index, 'Source %s not found in 6-year list' % sname
        rec6y = s6y.ix[sname]
        eflux6y = rec6y.flux * rec6y.e0**2*1e6
        fluxunc67 = rec6y.flux_unc
        ass = rec6y.associations
        if ass is not None:
            acat = ass['cat'][0]
            aprob = ass['prob'][0]
        else:
            aprob = 0
            acat=''
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
                good= rec.ts>10 and rec.a>0.001 and rec.a<0.25 and rec.locqual<8, 
                adeltats=rec.adeltats,
                aprob=aprob,
                acat=acat,
                )
        if len(months)==0:
            self.notdetected.append(sname)
            #print 'Source %s not detected in any month' % sname
            return None
        monthly = pd.DataFrame(months).T
        good=monthly.good
        pull_rms = monthly[good].pull.std()
        pull_mean = monthly[good].pull.mean()
        ratio = monthly[good].eflux.mean()/eflux6y
        adts_mean = monthly[good].adeltats.mean()
        adts_max  = monthly[good].adeltats.max()
        ngood=sum(good)
        tssum = sum(monthly.ts)
        tsmax = max(monthly.ts)
        return dict( monthly=monthly, 
                     nmonths=len(monthly),
                     tssum=tssum,
                     tsmax = tsmax,                        
                     ngood=ngood,
                     pull_rms=pull_rms,
                     mean=pull_mean,
                     flux_ratio=ratio,
                     ts = rec6y.ts,
                     eflux = eflux6y,
                     adts_mean=adts_mean,
                     adts_max = adts_max,
                     curvature=rec6y.curvature,
                     aprob=aprob,
                     acat=acat,
                     )

class ObssimAnalysis(Analysis):
    def __init__(self, path='$FERMI/skymodels/obssim_monthly'):
        super(ObssimAnalysis, self).__init__(path=path, prefix='o')

class SourceInfo(object):
    """Manage the all-year source information, adding monthly info
    """
    def __init__(self, transient_analysis=None, names=None, quiet=True):
        ta= transient_analysis if transient_analysis is not None else Analysis()
        self.ta=ta
        if names is None:
            names=ta.all_6y_names()
        print 'loading SourceInfo with %d sources' % len(names)
        info=dict()
        ta.notdetected=[]
        for s in names:
            if not quiet: print '.',
            month_dict = ta.monthly_info(s)
            info[s] = month_dict
        if len(ta.notdetected)>0:
            print '{} Sources where not detected in any month.'.format(len(ta.notdetected))
            #print '{}'.format(np.array(ta.notdetected))
        self.df= pd.DataFrame(info).T

    def __getitem__(self, sourcename):
        """access to a record by name"""
        return self.df.ix[sourcename]
        
    def make_simple_record(self, outfile='for_benoit.csv'):
        bigdict = dict()
        bigindex=0
        skipped = []
        for name, row in si.df.iterrows():
            #print name, row
            if row.monthly is None:
                skipped.append(name)
                continue
            for month, monthinfo in row.monthly.iterrows():
                if monthinfo.good:
                    #print month, monthinfo
                    bigdict[bigindex]=dict(name=name, month=month, ts=monthinfo.ts)
                    bigindex +=1
        bigdf =pd.DataFrame(bigdict).T
        print 'skipped {}'.format(skipped)
        bigdf.to_csv(outfile)

        
        
    
class PSRinfo(object):
    """ Manage PSR info for 6 years and the monthly models
    """
    def __init__(self, transient_analysis=None, names=None, quiet=True):
        ta= transient_analysis if transient_analysis is not None else Analysis()
        self.ta=ta
        if names is None:
            names=ta.psr_names()
        info=dict()
        for s in names:
            try:
                if not quiet: print '.',
                month_dict = ta.monthly_info(s)
                info[s] = month_dict
            except Exception,msg:
                if not quiet: print '\n',msg
        self.df= pd.DataFrame(info).T
        self.df['efficiency']=np.array(self.df.ngood/72.,float)

    def __getitem__(self, sourcename):
        """access to a record by name"""
        if sourcename.startswith('PSR '):
            return self.df.ix[sourcename]
        return self.df.ix['PSR '+sourcename]
        
    def psr_plots(self, sourcename, which=None, **kw):

        df = self.df.ix[sourcename]
        dfm=df.monthly
        hist_kw = dict(histtype='step', lw=2)
        good = dfm.good
        #print 'Good months: %d/%d' %(sum(good), len(good))
        if sum(good)==0: return None
        
        
        def histit(ax, name, bins, label=None, good_only=False, **kw):
            all = np.asarray(dfm[name],float).clip(bins[0],bins[-1])
            subset = np.asarray(dfm[good][name],float).clip(bins[0],bins[-1])
            if not good_only:
                ax.hist( all,  bins, color='red', **hist_kw)
            ax.hist( subset , bins, color='green', **hist_kw)
            ax.grid(True, alpha=0.5)
            ax.set_xlabel(name if label is None else label)
            if len(kw)>0: plt.setp(ax, **kw)

        def a_hist(ax):
            histit(ax, 'a' ,np.linspace(0,0.25,26),  'Major axis [deg]')

        def eflux_hist(ax, **kw):
            bins = kw.pop('bins', np.logspace(0,2,26) )
            histit(ax, 'eflux', bins, 'Energy Flux [eV]', xscale='log', **kw)
            ax.axvline(df.eflux, color='red')

        def ts_hist(ax):
            histit(ax, 'ts', np.logspace(1,2.5,21), 'TS' ,xscale='log')
            ax.axvline(df.ts/72, color='red')

        def lq_hist(ax):
            histit(ax, 'locqual', np.linspace(0,8,26),'Localization Quality')

        if which is not None:
            ax = kw.pop('ax', None)
            bins=kw.get('bins',None)
            label = kw.get('label', None)
            if ax is None:
                fig, ax = plt.subplots(1,1, figsize=(4,4))
            else:
                fig = ax.figure
            if which=='eflux':
                eflux_hist(ax, bins=bins, good_only=True)
            else:
                histit(ax, which, bins=bins, label=label, good_only=True)
            return fig
        
        fig, axx = plt.subplots(1,4, figsize=(12,4))
        for f,ax in zip([ts_hist, eflux_hist, a_hist, lq_hist], axx.flatten()): f(ax)
        plt.suptitle('Monthly analysis for %s' % sourcename)
        return fig 


    def good_months(self):
        psr_df = self.df
        fig, ax= plt.subplots(1,1,figsize=(5,5))
        ax.hist(psr_df.nmonths, np.linspace(0,72,73), histtype='step');
        ax.hist(psr_df.ngood, np.linspace(0,72,73), color='red', histtype='step');
        return fig
 
    def save_psr_plots(self, folder='psr_plots', psr_names=None):
        if not os.path.exists(folder):
            os.mkdir(folder)
        if psr_names is None: psr_names=self.ta.psr_names();
        for name in psr_names:
            print name,
            try:
                fig=self.psr_plots(name)
                if fig is not None:
                    fig.savefig('psr_plots/%s.png' %name.replace(' ','_') )
                    plt.close(fig)  
            except Exception, msg:
                print msg
    
    def pull_plot(self, xlim=(0.8,2.0), ylim=(0.2,2.1)):
        """
        
        """
        psr_df = self.df
        fig, ax= plt.subplots(1,1,figsize=(10,8))
        rms = np.asarray(psr_df.pull_rms, float)
        ngood = psr_df.ngood

        scat=ax.scatter(psr_df.flux_ratio.clip(*xlim), rms.clip(*ylim), s=50, 
            c=np.array(psr_df.ngood,float)*100/72., edgecolor='none' ) 
        plt.setp(ax, xlabel='flux ratio', ylabel='STDEV of pulls', 
            xlim=xlim, ylim=ylim, title='Pulsar monthly check of flux consistency')
        ax.axvline(1.0, color='k', ls='--'); ax.axhline(1.0,color='k', ls='--')
        ax.grid(True, alpha=0.5)
        cb = plt.colorbar(scat)
        cb.set_label('Efficiency')
        return fig
        
    def efficiency(self, value='eflux', bins=np.logspace(0,2,11), xlabel=None):
        """Make a plot of efficiency, by counting the number of months for a given pulsar,
        in groups of flux
        """
        df = self.df
        dom = bins
        f = df[value]
        if value=='ts': f/=72. # scale to expected TS
        cuts = [ (f<dom[i+1]) & (f>dom[i]) for i in range(len(dom)-1)]
        ng = df.ngood; 
        n = np.array([sum(c) for c in cuts],float)

        eff =[ng[c].mean()/72. for c in cuts]
        err = np.array([ng[c].std()/72. for c in cuts]) / np.sqrt(n)
        fig, ax = plt.subplots(figsize=(8,5))
        x = 0.5*(dom[1:]+dom[:-1])
        xerr = (x-dom[:-1], dom[1:]-x)
        ax.errorbar(x=x, xerr=xerr, y= eff, yerr=err, fmt='o')
        if xlabel is None:
            xlabel=r'$\mathsf{Energy\ Flux\ at\ pivot\ (eV\ cm^{-2}\ s^{-1})}$'
        plt.setp(ax, xscale='log', xlim=(bins[0],bins[-1]),
                  ylim=(0,1.1),
                title='Detection Efficiency')
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel('Efficiency', fontsize=14)
        ax.set_title('Pulsar monthly detection efficiency', fontsize=16)
        ax.grid(True, alpha=0.5)
        ax.axhline(1.0, color='k', ls=':')
        return fig 
        
    def association_ts(self, bins=np.linspace(0,5,21)):
        """Plot the mean of the association ts for all sources, and those with at least 50 efficiency
        """
        df = self.df
        fig, ax = plt.subplots(1,1, figsize=(5,5))
        hist_kw= dict(histtype='step', lw=2)
        ax.hist(df.adts_mean.clip(0,5), bins, label='all', **hist_kw );
        ax.hist(df.adts_mean[df.ngood>36].clip(0,5), bins,color='red',
                label='eff>0.5', **hist_kw );
        ax.grid(True, alpha=0.5)
        legend(ax)
        plt.setp(ax, xlabel='Mean of association TS', title='Association test');
        return fig

    def association_delta_ts(self, sname, ax=None, xlim=(-1,20)):
        """ histogram of the association delta ts for the source
        """
        if ax==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5))
        else: fig=ax.figure

        adts = self.df.ix[sname].monthly.adeltats
        ax.hist(np.asarray(adts,float).clip(*xlim), 
                np.linspace(xlim[0],xlim[1],22), histtype='step', lw=2);
        ax.grid(True, alpha=0.5)
        plt.setp(ax, xlabel='Association Delta TS', xlim=xlim, title=sname);
        return fig
    
class BZCATinfo(SourceInfo):

    def __init__(self, transient_analysis=None, quiet=True, path='$FERMI/skymodels/P301_monthly/month*'):
        if path.find('/month')>0:
            # get the 6-year bzcat summary, load sources 
            self.bz6=bz6 = pd.read_csv(os.path.expandvars(
                '$FERMI/skymodels/P301_6years/uw972/plots/associations/bzcat_summary.csv'),
                      index_col=0)
            bz6names = list(bz6.index)
            bz6snames = list(bz6.sname)
            super(BZCATinfo, self).__init__(transient_analysis, 
                                            names=bz6snames, quiet=quiet)
        else: 
            self.bz6=bz6=None
            
        # load each of the monthly
        import glob, pickle
        from IPython import display
        ff = sorted(glob.glob(
               os.path.expandvars(path)+'/plots/associations/bzcat_summary*.csv'));
        assert len(ff)>0, 'Falied to load any months!'
        self.bzm=bzm = [pd.read_csv(f, index_col=0) for f in ff]
        bznames=[] 
        for bz in bzm:
            bznames.append(list(bz.index))
        sets = map(set, bznames)
        
        # make a combined table of bzcat associations
        a = sets[0]
        for b in sets[1:]:
            a = a.union(b)
        bznames = sorted(a); 
        if not quiet:
            print 'Found %d unique names in %d months' % (len(bznames),len(ff)) 
        
        if bz6 is not None:
            # compare with the 6-year list
            insix =set(bznames).intersection(bz6names)
            print '%d monthly names are in the 6-year list of %d names'\
                % (len(insix), len(bznames))
        
        # make sumamry dataframe with monthly info for transient guys
        bb = dict()
        for bzname in bznames:
            tss = []
            msname = []
            monthnumber = []
            delta_ts = []
            for i,tlist in enumerate(bzm):
                if bzname not in tlist.index: continue
                j = list(tlist.index).index(bzname)
                tss.append( round(tlist.ts[j],1) )
                msname.append(tlist.sname[j] )
                monthnumber.append(i+1)
                delta_ts.append(round(tlist.deltats[j],1))
            bb[bzname] = dict(ts=np.array(tss),
                              msname=np.array(msname),
                             monthnumber=np.array(monthnumber),
                             delta_ts=np.array(delta_ts),
                             )
        df = pd.DataFrame(bb).T
        df.index.name='bzcat_name'
        df['tsmax'] = [max(ts) for ts in df.ts]
        df['months']= [sum(np.array(ts)>0) for ts in df.ts]
        df['type'] = [name[3] for name in df.index]
        if bz6 is not None:
            df['ts6'] = bz6.ts.round(1)
            df['y6name'] = bz6.sname
            df['ra']  = bz6.ra
            df['dec'] = bz6.dec
            df['roi'] = [-1 if np.isnan(ra) else Band(12).index(SkyDir(ra,dec))
                     for ra,dec in zip(df.ra,df.dec)]
            df['notin6y'] = [not isinstance(x, str) for x in df.y6name]; 
        else:
            df['notin6y'] = True
        self.bzdf = df
        
        # write out a file with new associations
        #df['notin6y'].to_csv('bzcat_new.csv')
        

    def bzcat_plots(self):
        si = self
        fig, axx = plt.subplots(2,2, figsize=(12,12))
        ax=axx[0,0]
        ax.hist(si.df.ngood, np.linspace(0,72,73), histtype='step')
        plt.setp(ax, xlabel='number of good months', 
                 xlim=(0,72), title='Efficiency')
        ax.grid(True, alpha=0.5)

        ax=axx[0,1]
        x = np.asarray(si.df.adts_mean, float)
        hkw = dict(bins=np.linspace(0,5,21), histtype='step',lw=2)
        ngood = np.asarray(si.df.ngood, float)
        ax.hist(x.clip(0,5), label='all', **hkw );
        ax.hist(x[ngood>18].clip(0,5),label='>18 good months', 
                color='red', **hkw);
        ax.hist(x[ngood>36].clip(0,5),label='>36 good months', 
                color='green', **hkw);
        legend(ax)
        plt.setp(ax, xlabel='mean(Delta TS)', title='Localization consistency')
        ax.grid(True, alpha=0.5)

        ax=axx[1,0]
        x = np.array(si.df.pull_rms, float).clip(0,5)
        hkw = dict(bins=np.linspace(0,5,21), histtype='step',lw=2, log=True)
        ax.hist(x,         label='all', **hkw)
        ax.hist(x[ngood>18], label='>18 good months', **hkw)
        ax.hist(x[ngood>36], label='>36 good months', **hkw)

        plt.setp(ax, xlabel='Pull RMS', ylim=(1,None), title='Variabity')
        ax.grid(True, alpha=0.5)
        legend(ax)

        ax=axx[1,1]
        x = np.array(si.df.flux_ratio, float).clip(0,5)
        hkw = dict(bins=np.linspace(0,5,21), histtype='step',lw=2, log=True)
        ax.hist(x,         label='all', **hkw)
        ax.hist(x[ngood>18], label='>18 good months', **hkw)
        ax.hist(x[ngood>36], label='>36 good months', **hkw)

        plt.setp(ax, xlabel='Flux ratio', ylim=(1,None), title='Variability')
        ax.grid(True, alpha=0.5)
        legend(ax);
        return fig

    def cumulative_ts(self):
        """
        """
        df = self.bzdf
        
        dfnew=df[df.notin6y] if self.bz6 is not None else df
        # make an integral logTS plot for the new 1-month guys

        fig, ax = plt.subplots(figsize=(6,6))
        tsmax=100
        bins = np.logspace(1,np.log10(tsmax),500)
        hist_kw = dict(bins=bins, cumulative=-1, lw=2, histtype='step', log=True)
        ax.axvline(25, color='k', ls='--', label='TS=25')

        #ax.hist(dfnew.tsmax[dfnew.months==1], dom, label='all', **hist_kw );

        for type,label in zip('QBGU', ['FSRQ', 'BL Lac', 'Galaxy', 'Unknown']):
            sel = (dfnew.type==type)
            if sum(sel)>0:
                try:
                    ax.hist(dfnew.tsmax[sel], label=label, **hist_kw)
                except: pass

        plt.setp(ax, xscale='log', xlabel='TS', ylim=(0.5,None), xlim=(10, tsmax),
                title='Integral TS distribution for new bzcat associations')
        ax.grid(True, alpha=0.8);

        legend(ax)
        return fig
        
    def monthly_detections_plot(self):
        assert False, 'Not working now'
        df=self.bzdf
        monthlist=range(1,73)
        nm =[sum( [ts[m]>0 for ts in df.ts]) for m in range(len(monthlist))]; 
        ngood=[sum( [((ts[m]>0) & b) for ts,b in zip(df.ts,df.notin6y)]) for m in range(len(monthlist))]; 
        #nm,ngood
        fig,ax = plt.subplots(figsize=(8,4))
        ax.plot(monthlist,[len(bz) for bz in self.bzm],
                'o', label='all associations')
        ax.plot(monthlist, ngood, 'sr', label='new associations')
        plt.setp(ax, xlim=(0.5, monthlist[-1]+0.5), 
                 xlabel='month number', ylim=(0,None),
                 ylabel='BZCAT sources found',
                title='BZCAT detections per month')
        ax.grid(True, alpha=0.5)
        ax.legend(loc='upper left');
        return fig
        
    def monthly_frequency_table(self):
        df = self.bzdf
        in6y = np.logical_not(df.notin6y)
        h1,b = np.histogram(df[in6y].months, range(1,15))
        h2 = np.histogram(df[df.notin6y].months, range(1,15))[0]
        return '<h3>Number of detected months</h3>'+\
            pd.DataFrame([h1,h2],index=['in 6y','new'],
                         columns=b[:-1]).to_html()
        
    def agn_type_table(self, name='monthly data'):
        df = self.bzdf
        dfnew = df[df.notin6y & (df.tsmax>10)]
        types = set(list(df.type))
        typedict = dict(B='BL Lac', G='Radio galaxy', Q='FSRQ', U='unknown')
        typenames= [typedict[t] for t in types] + ['Total']
        
        tc10 = [sum(dfnew.type==type) for type in types] + [len(dfnew)]
        tc25 = [sum(dfnew.type[dfnew.tsmax>25]==type) for type in types] + [sum(dfnew.tsmax>25)]
        tcdf = pd.DataFrame([tc10,tc25], index=['TS>10', 'TS>25',]).T
        tcdf.index=typenames
        return ('<h4>AGN type counts for {}</h4>'.format(name)+ tcdf.to_html())

class BZCAT(object):
    """Load the current catalog
        """
    def __init__(self, 
                 bzcat_file='/nfs/farm/g/glast/g/catalog/pointlike/fermi/catalog/srcid/cat/obj-blazar-bzcat.fits' 
                ):
        self.cat=pd.DataFrame(pyfits.open(bzcat_file)[1].data)
        # add useful things to the DataFrame
        #make a column without the question marks (should be another for that?
        self.cat['z'] = np.array(map(lambda x:x.replace('?',''), self.cat.Redshift), float)
        self.cat['xray']= np.array(self.cat['X-ray flux0.1-2.4 keV(1.e-12 cgs)'], float)
        self.catnames = list(self.cat['Source name'])
        self.cat['fsrq'] = np.array([n[3]=='Q' for n in self.catnames])
        from skymaps import SkyDir
        self.cat['glat'] = [
            SkyDir(s['RA (J2000.0)'], s['Dec (J2000.0)']).b()
                for n,s in self.cat.iterrows()]
        self.cat['hilat']= np.abs(self.cat.glat)>10

    def selection(self, colname, bznames):
        #add a column to tag the entries in the list bznames 
        self.cat[colname]= inbz = np.array([ name in bznames for name in self.catnames] ) 
    
    def tag(self, colname, bznames):
        #add a column to tag the entries in the list bznames 
        #assert np.all(np.array(bznames) in self.catnames)
        self.cat[colname]= inbz = np.array([ name in bznames for name in self.catnames] ) 
        print 'Tagged {} with {} BZCAT entries'.format(colname, len(bznames))

    def __getitem__(self, sourcename):
        assert sourcename in self.catnames, 'Name not found'
        return self.cat.ix[self.catnames.index(sourcename)]

    def describe(self):
        return self.cat.describe()

    def type_frequency(self):
        snames = self.cat['Source name']
        types = [n[3] for n in snames]
        typeset = set(types)
        trans = dict(B='BL Lac', Q='FSRQ', U='unknown', 
                     G='radio galaxy')

        type_frequency = dict(); 
        for x in typeset: type_frequency[trans[x]]=0
        for x in types:
            type_frequency[trans[x]] +=1
        df = pd.DataFrame(type_frequency, index=['frequency']).T
        df['fraction (%)']=(df.frequency/len(snames)*100).round(1)
        df.index.name='spectral type'
        return df

    def correlations(self, lacnames, bz6names, found_in_months):
        # unfinsished
        bzcatnames = list(self.cat['Source name'])
        inbz = np.array([ name in bz6names         for name in bzcatnames] ) 
        minbz =np.array([ name in found_in_months  for name in bzcatnames] )
        lacdfnames = list(lacdf['Source name'])
        
    def bzcat_hist(self, quantity, bins,
               label=None, title='', cut=None,):
        if label is None: label=quantity
        fig, ax = plt.subplots(figsize=(10,6))
        hist_kw = dict(bins=bins,  
                       histtype='stepfilled',alpha=0.2, lw=2, log=True)
        cat= self.cat
        x = cat[quantity]
        if cut is  None: 
            cut = np.ones(len(cat),bool)
        else: cut=cat[cut]
        clip = lambda x: np.array(x.clip(bins[0], bins[-1]),float)
        ax.hist(clip(x[cut]), label='all BZCAT', **hist_kw);
        ax.hist(clip(x[cut &cat.sixyr]), label='6-year sources',     hatch='\\', **hist_kw)
        ax.hist(clip(x[cut &cat.month]), label='monthly detections', hatch='/',  **hist_kw)
        ax.grid(True, alpha=0.5)
        plt.setp(ax, ylim=(0.6,None), xlabel=label, title=title)
        ax.legend(loc='best');


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
        
class CloseCut(object):
    """
    Manage selection of a set of sources, by removing those close to another set
    """
    def __init__(self, sources, reference):
        """ sources, reference : SkyDir lists
        """
        print 'Comparing {} sources with {} in reference set'.format(len(sources), len(reference))
        def closest(a, b):
            return np.degrees(min([a.difference(x) for x in b]))
        self.closediff = np.array([closest(a,reference) for a in sources])
    
    def plot(self, marker=0.5, bins=np.linspace(0,1,101), ax=None):
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(figsize=(6,6))
        else: fig=ax.figure
        ax.hist(self.closediff**2, bins, histtype='step', lw=2);
        ax.axvline(marker**2,color='red', label='{:.2f} deg'.format(marker))
        plt.setp(ax, xlabel='Distance**2 [deg**2]', title='Closest distance to reference source')
        ax.legend(); ax.grid(True, alpha=0.5)
        return fig
    
    def __call__(self, mindist=0.5):
        return self.closediff>mindist;
    
def plots2(self, cut=None, title='TS>10'):
    
    if cut is None: 
        cut = self.df.ts>10
    print 'Effect of cut {} for these plots: from {} to {}'.\
        format(title, len(self.df),sum(cut))
    df = self.df[cut]
    df['lq'] = np.asarray(df.locqual, float)
    unassoc = np.logical_not(self.assoc)
    lolat = np.logical_not(self.hilat)
    cuts = (lolat, self.hilat, self.assoc, unassoc)
    labels = 'lolat hilat assoc unassoc'.split()
    if cut is not None:
        cuts = [x & cut for x in cuts]
    print 'Number of sources with TS>25: {}'.format(sum(df.ts>25))
    
    #for c, name in zip(cuts, labels):
    #    print '{:10s}: {:6d}'.format(name, sum(c))
    print '         {:8>s} {:8>s}    sum'.format(labels[2], labels[3])
    m = np.zeros(4).reshape((2,2))
    sums = np.zeros(3)
    for i in (0,1): 
        for j in (2,3):
            m[i,j-2]= sum(cuts[i] & cuts[j])
    for i in (0,1):
        a,b =  m[i,0], m[i,1]
        sums+= np.array([a,b,a+b])
        print '{:6s}{:8.0f}{:8.0f}{:8.0f}'.format(labels[i],a,b,a+b)
    a,b = sum(cuts[0]), sum(cuts[1])
    print 'sum   {:8.0f}{:8.0f}{:8.0f}'.format(*sums)
    
    def make_hist(ax, var, bins, xlabel, legloc='upper right'):        
        hist_kw=dict(bins=bins, lw=2, histtype='step', log=True)
        #df['locqual'] = np.asarray(df.locqual, float)
        ax.hist(var, color='k', label='all', **hist_kw);
        #ax.hist(var[lolat],    color='orange',    label='lolat', **hist_kw);
        #ax.hist(var[self.hilat], color='orange', label='hilat', **hist_kw);
        ax.hist(var[self.assoc], color='green', label='assoc', **hist_kw);
        #ax.hist(var[unassoc],    color='red', label='unassoc', **hist_kw);
        leg = ax.legend(loc=legloc)
        for box in leg.get_patches():
            box._height=0; box._y=0.5
        plt.setp(ax, ylim=(1,None), xlabel=xlabel, xlim=(bins[0], bins[-1]))
        ax.set_xlabel(xlabel, fontsize=12)
        ax.grid(True, alpha=0.5)

    def lq_hist(ax):
        bins=np.linspace(0,10,51)
        make_hist(ax, df.lq, bins, xlabel='localization quality')
        
    def pindex_hist(ax):
        bins=np.linspace(0,4,41)
        make_hist(ax, df.pindex.clip(0,4), bins, 'photon index', legloc='upper left')
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
        ax.axvline(25, color='red', ls=':', lw=2, label='TS=25')
        legend(ax);
        plt.setp(ax, xscale='log')
    def kdpc(ax):
        bins= np.linspace(0,100, 51)
        make_hist(ax, df.kdepercentile, bins, 'Diffuse percentile', legloc='upper left')
        
    #fig, axx = plt.subplots(2,3, figsize=(15,12))            
    #for f,ax in zip([lq_hist, pindex_hist, eflux_hist, sin_glat, ts, kdpc,], 
    #                axx.flatten()): f(ax)
    #if title!='': fig.suptitle(title, size=16)
    fig, axx = plt.subplots(1,2, figsize=(10,5))            
    for f,ax in zip([pindex_hist,  ts, ], 
                    axx.flatten()): f(ax)
    if title!='': fig.suptitle(title, size=16)
    #return labels, cuts

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
    B12 = Band(12)
    def neighbors(i):
        iv = IntVector()
        n = B12.findNeighbors(int(i),iv)
        return list(iv)

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

def pair_correlation_plots(ta, cut=None):
    """ta: Analyis object
    """
    unassoc = np.logical_not(ta.assoc)
    if cut is None: 
        alldiffs =  pair_correlations(ta, ta.df)
        hilat_diffs = pair_correlations(ta, ta.df[ta.hilat])
        assoc_diffs =  pair_correlations(ta, ta.df[ta.assoc])
        unassoc_diffs = pair_correlations(ta, ta.df[unassoc])
    else:
        alldiffs =  pair_correlations(ta, ta.df[cut])
        hilat_diffs = pair_correlations(ta, ta.df[ta.hilat & cut])
        assoc_diffs =  pair_correlations(ta, ta.df[ta.assoc& cut])
        unassoc_diffs = pair_correlations(ta, ta.df[unassoc& cut])
    
    flg, axx = plt.subplots(1,2, figsize=(12,6))

    for bins, ax in zip(
            [np.linspace(0,10,51),np.linspace(0,0.15,31)], axx, ):
        ax.hist(alldiffs**2, bins, label='all',color='blue', lw=2, histtype='step');
        ax.hist(hilat_diffs**2, bins, label='|b|>10',color='orange', lw=2, histtype='step');
        try:
            ax.hist(assoc_diffs**2, bins, label='assoc',color='green', lw=2, histtype='step');
        except: pass
        ax.hist(unassoc_diffs**2, bins, label='unassoc',color='red', lw=2, histtype='step');
        ax.grid(True, alpha=0.5)
        plt.setp(ax, xlabel='distance squared [deg**2]', xlim=(0, bins[-1]))
        leg=ax.legend() 
        for box in leg.get_patches():
            box._height=0; box._y=0.5
    plt.suptitle('Two-source correlation', size=14);


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

def source_plots(self,bzcat, bzname):
    """ BZCATinfo object and a name, look up source info, start the model, 
    make plots"""
    from uw.like2 import process
    rec = self.bzdf.loc[bzname]
    sname = rec.msname[0] # get first if more than one
    month = rec.monthnumber[0]
    ts = rec.ts
    sel =bzcat.cat['Source name']==bzname
    z=float(bzcat.cat.z[sel])
    print 'redshift, source(s),  month(s), ts:', z, rec.msname, rec.monthnumber, ts
    os.chdir(
       os.path.expandvars('$FERMI/skymodels/P301_monthly/month%02d'%month))
    r=process.Process('.', quiet=True)
    # get the source info for this month, look up the source
    ss =pickle.load(open('sources.pickle'))
    r.setup_roi(int(ss.ix[sname].roiname[-4:]))
    fig, axx = plt.subplots(1,2, figsize=(10,5))
    r.plot_tsmap(sname, axes=axx[0]); r.plot_sed(axes=axx[1]);
    return fig
 
class ObssimDataComparison(object):
    def __init__(self, data_file='../P301_monthly/transients_cut.csv',
                        sim_file='obssim_transients_cut.csv'):
        self.data= pd.read_csv(data_file, index_col=0)
        self.data_months = self.data.month.max()
        print 'Read {} data sources from {}, {} months'.format(len(self.data), data_file, self.data_months)
        self.sim = pd.read_csv(sim_file, index_col=0)
        self.sim_months = len(set(self.sim.month))
        print 'Read {} simulated sources from {}, {} months'.format(len(self.sim),sim_file, self.sim_months)
        
    def __call__(self, val='ts', cut=None, xlabel='TS',bins= np.linspace(10,30,21),
                      ylim=(-2,None),xscale='linear'):
                      
        data_cut = (np.abs(self.data.glat)>10) & (self.data.ts>10)
        sim_cut =  (np.abs(self.sim.glat)>10)  & (self.sim.ts>10)
        print 'Initial selection: {} data and {} sim'.format(sum(data_cut), sum(sim_cut) )
        if cut is not None:
            print 'applying cut {}'.format(cut)
            data_cut= data_cut & self.data[cut]
            sim_cut = sim_cut & self.sim[cut]
        data_rate_sum = sum(data_cut)/float(self.data_months)
        sim_rate_sum = sum(sim_cut)/float(self.sim_months)
        data_vals = self.data[val][data_cut]
        sim_vals = self.sim[val][sim_cut]

        fig, ax = plt.subplots(figsize=(6,6))
        x = 0.5*(bins[:-1] + bins[1:])

        def whist( data, weight, err=None, **kw):
            x = 0.5*(bins[:-1] + bins[1:])
            h =ax.hist(data, weights=np.ones(len(data))*weight, bins=bins, 
                 alpha=0.4, lw=2, histtype='step', **kw)[0] 
            yerr = np.sqrt(h*weight) if err is None else err
            ax.errorbar(x, h, fmt='+' , yerr=yerr, capsize=0)
            return h, yerr
        data_rate, data_err =whist(data_vals, 1/float(self.data_months),
                                   color='b', label='data [{:.1f}]'.format(data_rate_sum))
        sim_rate, sim_err =whist(sim_vals, 1/float(self.sim_months), 
                                 color='r', label='sim [{:.1f}]'.format(sim_rate_sum))
        y = data_rate-sim_rate
        ax.hist( x, weights=y, bins=bins, color='g', alpha=0.2, 
            histtype='stepfilled', lw=2, 
            label='difference [{:.1f}]'.format(sum(y)) )
        ax.errorbar(x, y, yerr= np.sqrt(data_err**2+sim_err**2), fmt='+', color='g', lw=2)
        ax.legend()
        ax.grid(True, alpha=0.5)
        plt.setp(ax, xlabel=xlabel, xscale=xscale,ylim=ylim)
        ax.set_ylabel('Sources per month', fontsize=14)
        ax.set_title('Monthly Average Transient Sources vs {}'.format(val))

                      
