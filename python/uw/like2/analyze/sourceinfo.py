"""
Basic analyis of source spectra

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/sourceinfo.py,v 1.30 2016/04/27 02:22:06 burnett Exp $

"""

import os
import astropy.io.fits as pyfits
import cPickle as pickle
from collections import Counter
import numpy as np
import pylab as plt
import pandas as pd

from uw.utilities import makepivot
from . import analysis_base, _html
from .. import extended, configuration
from analysis_base import html_table, FloatFormat
from skymaps import SkyDir, Band

class SourceInfo(analysis_base.AnalysisBase): #diagnostics.Diagnostics):
    """Source spectral properties 
    <br>See <a href="../localization/index.html?skipDecoration"> localization </a> for localization plots.
    """
    require='pickle.zip'
    def setup(self, **kwargs):
        self.plotfolder='sources' #needed by superclass
        filename = 'sources.pickle'
        self.quiet = kwargs.pop('quiet', True)
        refresh = kwargs.pop('refresh', not os.path.exists(filename) or os.path.getmtime(filename)<os.path.getmtime('pickle.zip'))
        if refresh:
            files, pkls = self.load_pickles('pickle')
            assert len(files)==1728, 'Expected to find 1728 files'
            self.pkls = pkls # for debugging
            sdict= dict()
            try:
                get_cat3fgl = Cat_3fgl()
            except Exception, msg:
                print 'Could not load 3FGL: %s' % msg
                get_cat3fgl=None
            
            for pkl in pkls:
                roidir = pkl['skydir']
                for name, info in pkl['sources'].items():
                    model = info['model']
                    pars = np.empty(4); pars.fill(np.nan)
                    errs = np.empty(4); errs.fill(-2)
                    free = np.zeros(4, bool)
                    n = model.len()
                    pars[:n] = model.parameters
                    free[:len(model.free)] = model.free
                    try:
                        d = np.diag(model.get_cov_matrix()).copy()
                        d[d<0] =0
                        errs[:n] = np.sqrt(d)
                        errs[np.isnan(errs)]=-1
                        badfit = np.any(errs[model.free]<=0)
                        freebits= np.sum( int(b)*2**i for i,b in enumerate(model.free))
                        badbits = freebits &  np.sum( int(b)*2**i for i,b in enumerate(errs==0))
                        
                    except Exception, msg:
                        raise exception( 'fail errors for %s:%s' % (name, msg))
                        badfit = True
                        
                    try:
                        fitqual = round(sum(info['sedrec'].delta_ts),2)
                    except:
                        fitqual = np.nan
                    try:
                        fitndf  = sum(1-info['sedrec'].zero_fract)
                    except:
                        fitndf = 10
                    try:
                        npred = round(sum(info['sedrec'].npred),1)
                    except:
                        npred=np.nan
                    ellipse = info.get('ellipse', None)
                    moment  = info.get('moment', None)
                    has_moemnt = moment is not None
                    sdict[name] = info
                    pulsar = model.name.endswith('Cutoff')
                    betavalue = float(pars[2]) if not pulsar else np.nan
                    if pulsar: # factor to convert flux to prefactor
                        bvalue =1.0 if model.npar==3 else model['b']
                        prefactor = np.exp(-(model.e0/model['cutoff'])**bvalue)
                    else: prefactor = 1.0
                    sdict[name].update(
                        glat=info['skydir'].b(), glon=info['skydir'].l(),
                        npred=npred,
                        roiname=pkl['name'], 
                        pars= pars, errs=errs, free=free, badfit=badfit,
                        a = ellipse[2] if ellipse is not None else np.nan,
                        b = ellipse[3] if ellipse is not None else np.nan,
                        ang=ellipse[4] if ellipse is not None else np.nan,
                        locqual = round(ellipse[5],2) if ellipse is not None else np.nan,
                        delta_ts = ellipse[6] if ellipse is not None else np.nan,
                        ts_beta = info.get('ts_beta', None ),
                        moment=moment,
                        freebits= freebits,
                        badbits=badbits,
                        flux = prefactor * pars[0],
                        flux_unc = prefactor * errs[0],
                        pindex = pars[1],
                        pindex_unc = errs[1],
                        beta = betavalue,
                        beta_unc = errs[2] if not pulsar and pars[2]>0.002 else np.nan,
                        index2 = pars[3] if pulsar else pars[2],
                        index2_unc = errs[3] if pulsar and not np.isnan(pars[3]) else errs[2],
                        cutoff = pars[2] if pulsar else np.nan,
                        cutoff_unc = errs[2] if pulsar else np.nan,
                        e0 = model.e0,
                        modelname=model.name,
                        fitqual = fitqual,
                        fitndf  = fitndf,
                        eflux = prefactor*pars[0]*model.e0**2*1e6,
                        eflux100 = info.get('eflux', (np.nan,np.nan))[0],
                        eflux100_unc = info.get('eflux', (np.nan,np.nan))[1],
                        psr = pulsar,
                        cat3fgl = None if get_cat3fgl is None else get_cat3fgl(name),
                        transient= not info.get('fixed_spectrum', False) and not info['isextended'],
                        roi_dist= np.degrees(info['skydir'].difference(roidir)),
                        )
            df = pd.DataFrame(sdict).transpose()
            df.index.name='name'
            df['name'] = df.index
            ra = [x.ra() for x in df.skydir]
            dec = [x.dec() for x in df.skydir]
            df['ra'] = ra
            df['dec'] = dec
            self.df = df.sort_index(by='ra')
            self.df['hassed'] = np.array([self.df.ix[i]['sedrec'] is not None for i in range(len(self.df))])
            self.curvature(setup=True) # add curvature item
            self.df.to_pickle(filename)
            if not self.quiet:
                print 'saved %s' % filename

        else:
            if not self.quiet: print 'loading %s' % filename
            self.df = pd.read_pickle(filename)
        #self.df['flux']    = [v[0] for v in self.df.pars.values]
        #self.df['flux_unc']= [v[0] for v in self.df.errs.values]
        localized = ~np.array(pd.isnull(self.df.delta_ts))
        extended = np.array(self.df.isextended, bool)
        self.df['unloc'] = ~(localized | extended)
        self.df['poorloc'] = (self.df.a>0.2) | (self.df.locqual>8) | (self.df.delta_ts>2)
        self.df['flags'] = 0  #used to set bits below
        flags = self.df.flags
        pl = (self.df.poorloc | self.df.unloc) & (self.df.ts>10)
        flags[pl] += 8 ### bit 8
        #print '%d sources flagged (8) as poorly or not localized' % sum(pl)

 
        self.energy = np.sqrt( self.df.ix[0]['sedrec'].elow * self.df.ix[0]['sedrec'].ehigh )
            
    def skyplot(self, values, proj=None, ax=None, ecliptic=False, df=None,
                labels=True, title='', colorbar=True, cbtext='', **scatter_kw):
        """ 
        Make a sky plot of some quantity for a selected set of sources
        Parameters:
        ----------
            values: a DataFrame column, posibly a subset: 
                expect to have source name index to get position
            proj: None, or a function to map values to colors
            s : float
                size of dot to plot
                
            df : None | DataFrame
                if None, use self.df to look up glat and glon
        """
        assert hasattr(values, 'index'), 'skyplot: values arg must have index attribute'
        
        if df is None:
            df = self.df
        assert len(set(values.index).intersection(df.index))==len(values), 'skyplot: index values unknown'
        # generate arrays of glon and singlat using index 
        sd = df.ix[values.index, ['glat', 'glon']] # see page 101
        glon = sd.glon
        glon[glon>180]-=360
        singlat = np.sin(np.radians(list(sd.glat)))

        c = values if proj is None else map(proj, values)

        if ax is None:
            fig, ax = plt.subplots(figsize = (6,5))
        else: fig = ax.figure
        
        scatter_kw_default=dict(s=20, vmin=None, vmax=None, edgecolor='none')
        scatter_kw_default.update(scatter_kw)
        
        scat = self.basic_skyplot(ax, glon, singlat, c=c,
                title=title, ecliptic=ecliptic, colorbar=colorbar,cbtext=cbtext, **scatter_kw_default)
        return fig
        
    def fluxinfo(self, ib=0, cut=None):
        """ extract flux info for energy bin ib, return as a DataFrame
        """
        if cut is None: cut=self.df.ts>25
        s = self.df[cut]
        energy = self.energy[ib]
        fdata = np.array([s.ix[i]['sedrec'].flux[0] for i in range(len(s))])
        udata = np.array([s.ix[i]['sedrec'].uflux[0] for i in range(len(s))])
        ldata = np.array([s.ix[i]['sedrec'].lflux[0] for i in range(len(s))])
        fmodel = np.array([s.ix[i]['model'](energy)*energy**2*1e6 for i in range(len(s))])
        return pd.DataFrame(dict(fdata=fdata, udata=udata, ldata=ldata, fmodel=fmodel, 
                glat=s.glat, glon=s.glon, roiname=s.roiname),
            index=s.index).sort_index(by='roiname')

    def cumulative_ts(self, ts=None, tscut=(10,25), check_localized=True, 
            label=None,  other_ts=[], other_label=[], ax=None, legend=True):
        """ Cumulative test statistic TS
        
        A logN-logS plot, but using TS. Important thresholds at TS=10 and 25 are shown.
        """
        usets = self.df.ts if ts is None else ts
        df = self.df
        if ax is None:
            fig,ax = plt.subplots( figsize=(8,6))
        else: fig=ax.figure
        
        dom = np.logspace(np.log10(9),5,1601)
        ax.axvline(25, color='gray', lw=1, ls='--',label='TS=25')
        hist_kw=dict(cumulative=-1, lw=2,  histtype='step')
        ax.hist( usets ,dom,  color='k', label=label, **hist_kw)
        if len(other_ts)>0 :
            for ots, olab in zip(other_ts,other_label):
                ax.hist( ots, dom,  label=olab, **hist_kw)
        if check_localized:
            unloc = df.unloc
            ul = df[(unloc | df.poorloc) & (usets>tscut[0])] 
            n = len(ul)
            if n>10:
                ax.hist(ul.ts ,dom,  color='r', 
                    label='none or poor localization', **hist_kw)
                ax.text(12, n, 'none or poor localization (TS>%d) :%d'%(tscut[0],n), fontsize=12, color='r')
        plt.setp(ax,  ylabel='# sources with greater TS', xlabel='TS',
            xscale='log', yscale='log', xlim=(9, 1e4), ylim=(9,20000))
        ax.set_xticklabels([' ', ' ', '10', '100', '1000'])
        #ax.set_yticklabels(['', '10', '100', '1000'])
            
        # label the plot with number at given TS
        for t in tscut:
            n = sum(usets>t) 
            ax.plot([t,2*t], [n,n], '-k');
            ax.plot(t, n, 'og')
            ax.text(2*t, n, 'TS>%d: %d'%(t,n), fontsize=14, va='center')
                
        ax.grid(True, alpha=0.3)
        if (label is not None or other_label is not None) and legend: 
            leg =ax.legend()
            for patch in leg.get_patches():
                pbox = patch; pbox._height=0; pbox._y=5
        return fig

    def cumulative_counts(self):
        #assume rois.pickle available
        # extract the model counts for each source from it, and add to the DF
        if not hasattr(self,'roi_df'):
            self.roi_df = pickle.load(open('rois.pickle'))
        def counts(src_df,roi_df, name):
            """return fit counts for source name"""
            roiname = src_df.ix[name]['roiname']
            roi=roi_df.ix[roiname]
            names = list(roi['counts']['names'])
            i = names.index(name)
            if i<0: print name, i
            try: 
                t =np.sum(roi['counts']['models'][i][1])
                return t
            except:
                print 'fail to get counts for source', name, roiname, i
                return 0
        self.df['counts'] = [counts(self.df, self.roi_df, name) for  name in self.df.index]
        
    def non_psr_spectral_plots(self, index_min=1.0, index_max=3.5, beta_max=2.0):
        """ Plots showing spectral parameters for PowerLaw and LogParabola spectra
        Left: energy flux in eV/cm**2/s. This is the differential flux at the pivot energy
        <br> Center: the spectral index.
        <br> Right: the curvature index for the subset with log parabola fits.
        %(tail_check)s
        %(beta_check)s
        """
        fig, axx = plt.subplots( 1,3, figsize=(12,4))
        plt.subplots_adjust(wspace=0.2, left=0.05,bottom=0.15)

        t = self.df.ix[(self.df.ts>10)&(self.df.modelname=='LogParabola')]['ts flux pindex beta beta_unc freebits e0 roiname'.split()]
        t['eflux'] = t.flux * t.e0**2 * 1e6
        ax = axx[0]
        [ax.hist(t.eflux[t.ts>tscut].clip(4e-2,1e2), np.logspace(-2,2,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xscale='log', xlabel='energy flux', xlim=(4e-2,1e2)); ax.grid(); ax.legend(prop=dict(size=10))
        ax = axx[1]
        [ax.hist(t.pindex[t.ts>tscut].clip(index_min,index_max), np.linspace(index_min,index_max,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xlabel='spectral index'); ax.grid(); ax.legend(prop=dict(size=10))
        ax = axx[2]
        sel=(t.ts>tscut)&(t.beta>0.01)
        if sum(sel)>0:
            [ax.hist(t.beta[sel].clip(0,beta_max), np.linspace(0,beta_max,26), label='TS>%d' % tscut) for tscut in [10,25] ]
            plt.setp(ax, xlabel='beta'); ax.grid(); ax.legend(prop=dict(size=10))
        # get tails
        tail_cut = (t.eflux<5e-2) | ((t.pindex<index_min) | (t.pindex>index_max))& (t.beta==0) | (t.beta>beta_max) | (t.beta<0)
        
        if sum(tail_cut)>0:
            tails=t[tail_cut]['ts eflux pindex beta freebits roiname'.split()].sort_index(by='roiname')
            filename = 'non_pulsar_tails.html'
            html_file = self.plotfolder+'/%s' % filename
            #html = tails.sort_index(by='roiname').to_html(float_format=FloatFormat(2))
            self.tail_check = html_table(tails, name=self.plotfolder+'/pulsar_tails', 
                heading='<h4>Table of %d sources on tails</h4>'%len(tails),
                float_format=FloatFormat(2))
            #open(html_file,'w').write('<head>\n'+ _html.style + '</head>\n<body>'+ html+'\n</body>')
            #self.tail_check = '<p><a href="%s?skipDecoration">Table of %d sources on tails</a>: '% (filename, len(tails))
            self.tail_check += 'Criteria: require index between 1 and 3.5 for powerlaw, beta<2.0 for log parabola'
            
            # flag sources
            flags = self.df.flags
            tails = tails.index
            flags[tails] += 1 ### bit 1
            print '%d sources flagged (1) in tails of flux, index, or beta' % len(tails)
        else:
            self.tail_check ='<p>No sources on tails'

        # check errors, especially that beta is at least 2 sigma
        self.beta_check=''
        #beta_bad = (t.beta>0.001) & ((t.beta_unc==0) | (t.beta/t.beta_unc<2) | (t.freebits!=7))
        #if sum(beta_bad)>0:
        #    print '%d sources fail beta check' % sum(beta_bad)
        #    self.beta_check = html_table(t[beta_bad]['ts beta beta_unc freebits roiname'.split()], 
        #        name=self.plotfolder+'/beta_check',
        #        heading = '<h4>Table of %d sources failing beta 2-sigma check</h4>'%sum(beta_bad),
        #        float_format=FloatFormat(2))
            
        return fig
        
    def beta_check(self):
        """Check beta
        
        <p>Compare TS for power-law vs. log parabola
        
        Histograms of TS_beta, the TS difference between beta=0 and a best fit
        (beta_check_note)s
        """
        df=self.df
        check = np.array([np.isfinite(x) for x in self.df.ts_beta]);
        beta_check_note=''
        if sum(check)==0:
            print 'No ts_beta values set'
            self.beta_check_note='<p>No beta analysis was done'
            return
        LP= df.freebits==7
        PL = df.freebits==3
        fig,ax = plt.subplots(figsize=(6,6))
        xmax=25; dom=np.linspace(0,xmax,26)
        args= dict(log=True, histtype='step', lw=2)
        ax.hist(df.ts_beta[PL].clip(0,xmax), dom,  label='PowerLaw', **args)
        ax.hist(df.ts_beta[LP].clip(0,xmax),dom, color='orange', label='LogParabola', **args)
        plt.setp(ax, xlabel='TS_beta', ylim=(1,None))
        ax.grid(); ax.legend()
        return fig
    
    def pulsar_spectra(self, index_min=0.0, index_max=2.5, cutoff_max=8000):
        """ Distributions for sources fit with PLSuperExpCutoff spectral model, mostly LAT pulsars
        
        For each plot, the subset with a poor fit is shown.
        %(pulsar_tail_check)s
        %(pulsar_fixed)s
        %(pulsar_b)s
        """
        fig, axx = plt.subplots( 1,4, figsize=(14,4))
        plt.subplots_adjust(wspace=0.3, left=0.05,bottom=0.15)
        psrmodel = (self.df.ts>10) & (self.df.modelname=='PLSuperExpCutoff')
        t = self.df.ix[psrmodel]\
            ['ts flux pindex cutoff e0 index2 index2_unc roiname freebits fitqual'.split()]
        t['eflux'] = t.flux * t.e0**2 * 1e6
        badfit = t.fitqual>30

        def plot1(ax, efmin=1e-1,efmax=1e3):
            bins = np.logspace(np.log10(efmin),np.log10(efmax),26)
            vals = t.eflux.clip(efmin,efmax)
            ax.hist(vals, bins )
            if sum(badfit)>0:
                ax.hist(vals[badfit], bins, color='red', label='poor fit')
            plt.setp(ax, xscale='log', xlabel='energy flux', xlim=(efmin,efmax)); ax.grid(); 
            ax.legend(prop=dict(size=10))

        def plot2(ax):
            bins = np.linspace(index_min,index_max,26)
            vals = t.pindex.clip(index_min,index_max)
            ax.hist(vals, bins)
            if sum(badfit)>0:
                ax.hist(vals[badfit], bins, color='red', label='poor fit')
            plt.setp(ax, xlabel='spectral index'); ax.grid(); 
            ax.legend(prop=dict(size=10))
            
        def plot3(ax):
            bins = np.linspace(0,cutoff_max/1e3,26)
            vals = t.cutoff.clip(0,cutoff_max) /1e3
            ax.hist(vals, bins)
            if sum(badfit)>0:
                ax.hist(vals[badfit], bins, color='red', label='poor fit')
            plt.setp(ax, xlabel='cutoff energy (GeV)'); ax.grid()
            ax.legend(prop=dict(size=10))
            
        def plot4(ax, xlim=(0,cutoff_max)):
            xvals = t.cutoff.clip(*xlim) / 1e3
            yvals = t.pindex.clip(index_min,index_max)
            ax.plot(xvals, yvals, 'o')
            if sum(badfit)>0:
                ax.plot(xvals[badfit], yvals[badfit], 'or', label='poor fit')
            plt.setp(ax, xlabel='cutoff energy', ylabel='spectral index',
                xlim=(xlim[0]-0.1,1.03*xlim[1]/1e3), ylim=(index_min-0.1, index_max+0.1),
                )
            ax.grid(); 
            ax.legend(loc='lower right', prop=dict(size=10))

        for f,ax in zip((plot1,plot2,plot3,plot4,), axx.flatten()): f(ax)
        flags = self.df.flags
        tail_cut = (t.pindex<=index_min) | (t.pindex>index_max) | (t.cutoff>cutoff_max)
        tails = t.ix[tail_cut].index
        flags[tails] += 1 ### bit 1
        print '%d pulsar sources flagged (1) in tails of  index or cutoff' % sum(tail_cut)
        if sum(tail_cut)>0:
            tails=t[tail_cut]['ts eflux pindex cutoff freebits roiname'.split()]
            filename = 'pulsar_tails.html'
            html_file = self.plotfolder+'/%s' % filename
            #html = tails.sort_index(by='roiname').to_html(float_format=FloatFormat(2))
            html = html_table(tails.sort_index(by='roiname'), float_format=FloatFormat(2))
            open(html_file,'w').write('<head>\n'+ _html.style + '</head>\n<body>'+ html+'\n</body>')
            self.pulsar_tail_check = '<p><a href="%s?skipDecoration">Table of %d sources on tails</a>: '% (filename, len(tails))
            self.pulsar_tail_check += 'Criteria: require index between 0 and 2.5, cutoff<8 GeV'
        else:
            self.pulsar_tail_check ='<p>No sources on tails'

        # table of pulsars with b<1

        tt=t[t.index2<1]['ts fitqual pindex cutoff index2 index2_unc'.split()]
        tt['significance'] = (1-tt.index2)/tt.index2_unc
        self.pulsar_b = html_table(tt,
            name=self.plotfolder+'/pulsar_b',
            heading='<h4>Table of %d sources with b&lt;1</h4>' % len(tt),
            float_format=FloatFormat(2))
        print '%d pulsar sources with b<1' %len(tt)

        # table of fits with any fixed parame er other than b
        tt = t[((np.array(t.freebits,int)&7) != 7)]['ts fitqual pindex cutoff freebits roiname'.split()].sort_index(by='roiname')
        if len(tt)>0:
            print '%d pulsar-like sources with fixed parameters' %len(tt)
            self.pulsar_fixed= html_table(tt, name=self.plotfolder+'/pulsar_fixed', 
                heading='<h4>%d pulsar-like sources with fixed parameters</h4>' %len(tt),
                float_format=FloatFormat(2))
        else: self.pulsar_fixed=''
        return fig
    
    def ecliptic_hist(self, ax=None, title=''):
        ea = map(self.ecliptic_angle, self.df.skydir)
        fig, ax = self.get_figure(ax)
        ax.hist( ea, np.linspace(-90, 90, 91))
        plt.setp(ax, xlim=(-90,90), xlabel='ecliptic latitude')
        ax.set_xticks([-90, -45, 0, 45, 90])
        if title: ax.set_title(title)
        ax.grid(True)
        return fig
        
    def fit_quality(self, xlim=(0,50), ndf=10, tsbandcut=20, grid_flag=True, make_table=True, legend_flag=True):
        """ Spectral fit quality
        This is the difference between the TS from the fits in the individual energy bands, and that for the spectral fit.
        It should be distributed approximately as chi squared of at most 14-2 =12 degrees of freedom. 
        However, high energy bins usually do not contribute, so we compare with ndf=%(ndf)d.
        All sources with TS_bands>%(tsbandcut)d are shown.<br>
        <b>Left</b>: Power-law fits. Tails in this distribution perhaps could be improved by converting to log parabola. 
        <br><b>Center</b>: Log parabola fits.
        <br><b>Right</b>: Fits for the pulsars, showing high latitude subset.
        <br><br> Averages: %(fit_quality_average)s
        %(badfit_check)s
        %(poorfit_table)s

        """
        from scipy import stats
        fig, axx = plt.subplots(1,3, figsize=(12,6))
        plt.subplots_adjust(left=0.1)
        s = self.df
        psr = np.asarray(s.psr, bool)
        fq = np.array(s.fitqual, float)
        beta = s.beta
        logparabola = (~psr) & (beta>0.01)
        powerlaw = (~psr) & (beta.isnull() | (beta<0.01) )

        self.tsbandcut=tsbandcut
        cut = np.array((s.band_ts>tsbandcut) & (fq>0) , bool)
        
        dom = np.linspace(xlim[0],xlim[1],26)
        d = np.linspace(xlim[0],xlim[1],51); delta=dom[1]-dom[0]
        chi2 = lambda x: stats.chi2.pdf(x,ndf)
        fudge = 1.0 # to scale, not sure why
        hilat = np.abs(self.df.glat)>5
        self.average = [0]*4; i=0
        def tobool(a): return np.array(a, bool)
        for ax, label, cut_expr in zip(axx[:2], ('power law', 'log-normal',), ('powerlaw','logparabola')):
            mycut=tobool(cut&eval(cut_expr))
            count = sum(mycut)
            if count==0:
                print 'Not generating plot for %s' % label
                continue
            ax.hist(fq[mycut].clip(*xlim), dom, histtype='stepfilled',  label=label+' (%d)'%count)
            self.average[i] = fq[mycut].mean(); i+=1
            ax.plot(d, chi2(d)*count*delta/fudge, 'r', lw=2, label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
            ax.grid(grid_flag); ax.set_xlabel('fit quality')
            if legend_flag: ax.legend(prop=dict(size=10))
            else: ax.set_title(label)
            
        def right(ax, label='exponential cutoff', cut_expr='psr'):
            mycut= cut & (psr) #tobool(cut&eval(cut_expr))
            count = sum(mycut)
            if count==0: return
            if legend_flag:
                labels = [label+' (%d)' %count,label+' [|b|>5] (%d)' %sum(mycut*hilat),r'$\mathsf{\chi^2\ ndf=%d}$'%ndf]
            else:
                labels = ['all', '|b|>5', '_nolegend_']
                ax.set_title(label)
            ax.hist(fq[tobool(mycut)].clip(*xlim), dom, histtype='stepfilled', label=labels[0])
            ax.hist(fq[tobool(mycut&hilat)].clip(*xlim), dom, histtype='stepfilled', color='orange', 
                label=labels[1])
            self.average[i]   = fq[tobool(mycut&hilat)].mean()
            self.average[i+1] = fq[tobool(mycut&(~hilat))].mean()
            ax.plot(d, chi2(d)*count*delta/fudge, 'r', lw=2, label=labels[2])
            ax.grid(grid_flag);ax.set_xlabel('fit quality')
            ax.legend(loc='upper right', prop=dict(size=10))
        
        right(axx[2])
        self.df['badfit2'] =np.array(self.df.badfit.values, bool)
        t = self.df.ix[(self.df.badfit2) & (self.df.ts>10)].sort_index(by='roiname')
        print '%d sources with bad fits' %len(t)
        if len(t)>0:
            print '%d sources with missing errors' % len(t)
            self.badfit = t[['ts', 'freebits', 'badbits', 'pindex', 'beta', 'e0','roiname']]
            self.badfit_check = html_table(self.badfit, name=self.plotfolder+'/badfits', 
                heading='<h4>%d Sources with missing errors</h4>' % len(t), float_format=FloatFormat(1))
        else: self.badfit_check = '<p>All sources fit ok.'
        self.fit_quality_average =  ', '.join( map(lambda x,n :'%s: %.1f' %(n,x) ,
                            self.average, 'powerlaw logparabola expcutoff(hilat) expcutoff(lolat)'.split()) )
        self.ndf=ndf
        print 'fit quality averages:', self.fit_quality_average
        if make_table:
            # Make table of the poor fits
            s['pull0'] = np.array([x.pull[0] if x is not None else np.nan for x in s.sedrec ])
            t =s.ix[((s.fitqual>30) | (np.abs(s.pull0)>3)) & (s.ts>10) ]['ra dec glat fitqual pull0 ts modelname freebits index2 roiname'.split()].sort_index(by='roiname')
            if len(t)==0:
                self.poorfit_table= '<p>No poor fits found'
            else:
            
                self.poorfit_table  = html_table(t, name=self.plotfolder+'/poorfit', 
                        heading='<h4>Table of %d poor spectral fits</h4>'%len(t),
                        float_format=FloatFormat(2),
                        formatters=dict(ra=FloatFormat(3), dec=FloatFormat(3), ts=FloatFormat(0),index2=FloatFormat(3)))

                # flag sources that made it into the list
                self.df.flags[t.index] = np.asarray(self.df.flags[t.index],int) | 2
                print '%d sources flagged (2) as poor fits' %len(t)
        return fig
      
    def poor_fit_positions(self):
        """ Positions of poorly-fit sources
        Selection: fitqual>30 or |pull|>3
        """
        s = self.df
        s['pull0'] = np.array([x.pull[0] if x is not None else 0 for x in s.sedrec])
        poor = ( (s.fitqual>30) | (np.abs(s.pull0)>3)) & (s.ts>10) 
        return self.skyplot(s.fitqual[poor], vmin=30, vmax=100, cbtext='TS')
        
    def pivot_vs_e0(self, xylim=(100, 1e5)):
        """ pivot vs e0
        The reference energy, e0, is fixed except by a special run that iterates until the measured pivot energy, 
        which is the energy at which the differential flux uncertainty is minimum is the same. This plot checks that by measuring the pivot energy, and comparing it with the current reference. Note that e0 is required to be in the range 200 MeV to 100 GeV.
        """
        fig, ax = plt.subplots(figsize=(5,5))
        s = self.df
        cut = s.ts>10
        df = self.df
        offset = df.pivot_energy/df.e0-1.0
        df['offset'] = offset
        not_converged = cut & ( (abs(offset)>0.3) &(df.e0<99999) & (df.e0>201) |
                                 (df.e0>99999) & (df.pivot_energy<df.e0*0.95) 
                                | (df.e0<201) & (df.pivot_energy>df.e0*1.05)
                               )
        print 'Pivot needs fixing: %d sources' % sum(not_converged)
        self.pivotfix = tofix=df[not_converged]['ts pivot_energy e0 offset roiname'.split()].sort_index(by='roiname')

        ax.plot(s.e0[cut].clip(*xylim), s.pivot_energy[cut].clip(*xylim), '.')
        plt.setp(ax, xscale='log',xlabel='e0', xlim=xylim, 
                    ylabel='pivot', yscale='log', ylim=xylim)
        ax.set_title('compare calculated pivot with e0', fontsize=10)
        ax.grid()
        return fig
        
    def fitquality(self, grid=True):
        """Fit Quality
        
        Left: fit quality histogram; right fit quality vs. TS'

        """
        fig, axs = plt.subplots(1,2, figsize=(7,3))
        plt.subplots_adjust(wspace=0.35)
        s = self.df
        fitqual = s.band_ts-s.ts
        from scipy import stats
        ndf=12
        chi2 = lambda x: stats.chi2.pdf(x,ndf)
        d = np.linspace(0,100,51); delta=d[1]-d[0]
        ax =axs[0]
        ax.hist(fitqual, d, log=False);
        ax.hist(fitqual[s.ts>500], d, label='TS>500');
        ax.plot(d, chi2(d)*len(fitqual)*delta/1.6, 'r', label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
        plt.setp(ax, xlabel='fit qual', ylim=(0,500))
        ax.grid(grid); ax.legend(prop=dict(size=10))
        ax = axs[1]
        ax.plot(s.ts, fitqual, '.'); 
        plt.setp(ax, xscale='log', xlabel='TS', xlim=(10,1e5),
             ylabel='fit qual',ylim=(1,1e3),yscale='log')
        ax.grid(grid)
        return fig

    def flux_uncertainty(self):
        """ flux uncertainty compared with TS
        
        """
        fig, axx = plt.subplots(1,2, figsize=(10,4))
        plots=[]
        relflux_unc= self.df.flux_unc/self.df.flux
        ts = np.asarray(self.df.ts, float)
        ru= np.array(relflux_unc*100.,float)

        def plot1(ax):   
            dom = np.logspace(0,2,26)
            ax.hist(ru[ts>9], dom, label='%d sources'% sum(ts>9))
            for tsmin in (25,100,1000):
                ax.hist(ru[ts>tsmin], dom, label='TS<%d' % tsmin )
            plt.setp(ax, xscale='log', xlabel='relative flux uncertainty (%)', xlim=(1,100))
            ax.set_xticklabels([1,10,100])
            ax.grid()
            ax.legend(loc='upper left', prop=dict(size=10))
        plots.append(plot1)
            
        def plot2(ax):
            ax.plot(ts, ru*np.sqrt(ts)/100, '.')
            plt.setp(ax, xlabel='TS', xlim=(10,10000), xscale='log',
                 yscale='log',ylabel='ru*sqrt(ts)', ylim=(0.8,4))
            ax.plot([0.1,100], [0.1,100],'-g')
            ax.grid()
        plots.append(plot2)
            
        for plotf, ax in zip( (plots), axx.flatten(),):
            plotf(ax)
        return fig
        
    def spectral_fit_consistency_plots(self, energy=133., minflux=2.0, 
            title = 'low energy fit consistency',
            three_plots=False,
        ):
        """ Spectral fit consistency for the lowest energy bin
        
        These plots show the consistency of the lowest energy band with the spectrum
        defined by the full fit. <br>
        <b>Left</b>: distribution of the "pull" <br>
        <b>Right</b>: position in the sky of sources with |pull]>3 <br>
        """
        hassed = np.array([self.df.ix[i]['sedrec'] is not None for i in range(len(self.df))]) 
        nosed = (self.df.ts>10) & ~ hassed
        if sum(nosed)>0:
            print '+++Warning: %d TS>10 sources without sed info' % sum(nosed)
            print self.df[~hassed]['ts roiname'.split()][:min(20, sum(nosed))]
        cut= (self.df.ts>25) & hassed
        s = self.df[self.df.hassed] #[cut]
        
        fdata = np.array([s.ix[i]['sedrec'].flux[0] for i in range(len(s))])
        udata = np.array([s.ix[i]['sedrec'].uflux[0] for i in range(len(s))])
        ldata = np.array([s.ix[i]['sedrec'].lflux[0] for i in range(len(s))])
        pull = np.array([s.ix[i]['sedrec'].pull[0] for i in range(len(s))])
        fmodel = np.array([s.ix[i]['model'](energy)*energy**2*1e6 for i in range(len(s))])
        glat = np.array([x.b() for x in s.skydir])
        fluxcut = fmodel>minflux
        latcut  = abs(glat)>5.0
        hilat = fluxcut & (latcut)
        lolat = fluxcut & (~latcut)
        
        #lowebad = np.abs(pull)>3
        #s.flags[lowebad] += 4
        #print 'Tagged %d sources with lowebad bit (4)' % sum(lowebad)

        lowebad = np.asarray((np.abs(pull)>3) , bool)
        self.df.flags[lowebad] = np.array(self.df.flags[lowebad],int) | 4
        print 'Tagged %d sources with lowebad, abs(pull0)>3, bit (4)' % sum(lowebad)

        y = fdata/fmodel
        ylower, yupper =[(fdata-ldata)/fmodel,(udata-fdata)/fmodel]
        xhi,yhi,yerrhi = fmodel[hilat], y[hilat], [ylower[hilat],yupper[hilat]]
        xlo,ylo,yerrlo = fmodel[lolat], y[lolat], [ylower[lolat],yupper[lolat]]
        
        def error_bar(ax):
            ax.errorbar(x=xhi, y=yhi, yerr=yerrhi, fmt='og', label='%d hilat sources'%sum(hilat))
            ax.errorbar(x=xlo, y=ylo, yerr=yerrlo, fmt='or', label='%d lowlat sources'%sum(lolat))
            plt.setp(ax, xlabel=r'$\mathsf{model\ flux\ (eV\ cm^{-2} s^{-1}})$', xscale='log', 
                ylabel='data/model', ylim=(0,2.5), xlim=(minflux, 100) )
            #ax.set_xticks([2,5,10,20,50,100])
            ax.set_title( title, fontsize='medium')
            ax.legend(loc='upper left', prop=dict(size=10))
            ax.grid()  

        def hist(ax):
            hist_kw=dict(bins=np.linspace(-3,3,25), lw=2, histtype='step')
            q=pull.clip(-3,3) 
            
            ax.hist(q[hilat], color='g',  label='%d hilat sources'%sum(hilat),  **hist_kw)
            ax.hist(q[lolat], color='r',  label='%d lowlat sources'%sum(lolat), **hist_kw)
            ax.set_xlabel('pull')
            ax.axvline(0, color='k')
            ax.set_xlim((-3,3))
            ax.set_title( title, fontsize='medium')
            ax.legend(loc='upper right', prop=dict(size=10))
            ax.grid()  

        def skyplot(ax):
            pdf = pd.DataFrame(dict(pull=pull), index=s.index) # to atatch indx
            self.skyplot(pdf.pull[lowebad], ax=ax, vmin=-3, vmax=3, title=title, cbtext='pull')

        if three_plots:
            fig,ax = plt.subplots(1,3, figsize=(12,5))
            plt.subplots_adjust(wspace=0.3, left=0.05)
            for f, ax in zip( (hist, error_bar, skyplot), ax.flatten()):
                f(ax=ax)
        else:
            fig,ax = plt.subplots(1,2, figsize=(12,5))
            plt.subplots_adjust(wspace=0.3, left=0.05)
            for f, ax in zip( (hist, skyplot), ax.flatten()):
                f(ax=ax)
        return fig
        
    def census(self, primary_prefix='P86Y'): #'P7R4'):
        """Census
        
        %(census_html)s
        <p>
        In this table of prefixes, the columns are the number of sources with TS greater than the header value. 
        The row labels are the first four characters of the source name, except 'ext' means extended.
        %(suffix_html)s
        """
        df = self.df
        extended = np.asarray(df.isextended.values,bool)
        pointsource = ~extended

        def count(prefix, tsmin):
            if prefix=='ext':
                return sum(extended & (df.ts>tsmin))
            elif prefix=='total':
                return sum(df.ts>tsmin)
            names = df[pointsource & (df.ts>tsmin)]['name'].values    
            return sum([n.startswith(prefix) for n in names])
        if count('ext',0)>0:
            prefixes = list(set( n[:4] for n in df[pointsource]['name'])) +['ext', 'total']
        else:
            prefixes = list(set( n[:4] for n in df[pointsource]['name'])) +['total']
        
        census = dict()
        prefixes = sorted(prefixes)
        for x in (0, 5, 10, 25):
            census[x] = [count(prefix, x) for prefix in prefixes]
        self.census_data=pd.DataFrame(census, index=prefixes)
        self.census_html = '\n<h4>Prefixes</h4>\n'\
            +html_table(self.census_data, maxlines=20, href=False)
        
        # now check suffixes
        self.primary_prefix=primary_prefix
        suffixes=[s[-1] for s in self.df.index.values if s.startswith(primary_prefix) and s[-1]>'9']
        c=Counter(suffixes)
        scounts = lambda  r : int(sum([c[x] for x in r if x in c.keys()]))
        suffixranges = ('ABCDEF', 'GHIJKL', 'MNO', 'PQRSTUVW', 'XYZ')
        sdict = dict([(r[0]+'-'+r[-1], [scounts(r)]) for r in suffixranges])
        total=sum([ s[0] for s in sdict.values()])
        
        if total>0:
            self.suffix_html ="""\
            <br>This second table has the suffixes; for sources with the
            "%s" prefix. Each group, except 'X-Z', 
            represents a set of seeds added to the original
            list of sources. The last row, 'X-Z', were added by hand."""% primary_prefix\
            + '\n<h4>Suffixes</h4>\n'\
            + html_table(pd.DataFrame(sdict, index=['freq']).T, href=False)
        else:
            self.suffix_html = '\n<p>No suffixes found'

    def flag_proc(self):
        """ Flagged source summary:
        %(flagged_link)s
        """
        # Generate summary table for flagged sources (except 
        t =self.df[(self.df.flags& 7 >0) & (self.df.ts>10)]['ra dec ts fitqual pull0 eflux pindex beta cutoff index2 flags roiname'.split()]
        t.to_csv('flagged_sources.csv')
        print 'wrote %d sources to flagged_sources.csv' % len(t)
        
        num=[sum(self.df.flags & 2**b > 0) for b in range(4)] 
        flagtable=pd.DataFrame(dict(number=num, description=('tails','poor fits','low energy bad', 'poor localization') ))
        flagtable.index.name='bit'
        self.flagged_link = """\
        <p>A number of these sources have been flagged to indicate potential issues with the spectral fit. 
        The flag bits and number flagged as such are:
        %s<br>  """ % html_table(flagtable, href=False)
        try:
            pc =makepivot.MakeCollection('flagged sources %s' % os.path.split(os.getcwd())[-1], 'sedfig', 'flagged_sources.csv')
            self.flagged_link += """\
            <p>These can be examined with a 
            <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
            which requires Silverlight."""  % pc.cId
        except Exception, msg: 
            print "**** Failed to make pivot table, perhaps need to run sedinfo first: %s" % msg
        return None

    def flux_table(self, source_name):
        """ Return a DataFrame describing the sed info for the given source
        Columns, with energy fluxes in eV units, are:
            flux, lflux, uflux : measured flux, lower and upper uncertainty, or 0,0, 95% limit
            mflux : predicted flux for this bin
            TS : Test Statistic for the signal
            pull : signed square root of the 
        """
        si = self.df.ix[source_name]['sedrec']
        pull = np.sign(si.flux-si.mflux) * np.sqrt(si.delta_ts.clip(0,100))
        return pd.DataFrame(dict(flux=si.flux.round(1), TS=si.ts.round(1), lflux=si.lflux.round(1), 
            uflux=si.uflux.round(1), model=si.mflux.round(1), pull=pull.round(1) ),
                index=np.array(np.sqrt(si.elow*si.ehigh),int), columns='flux lflux uflux model TS pull'.split())
    
    def find_nearest_to(self, *pars):
        """Given the ra, dec, or a skydir, find nearest source, return name, distance"""
        from skymaps import SkyDir

        sd = SkyDir(*pars) if len(pars)==2 else pars[0]
        dists = [x.difference(sd) for x in self.df.skydir]
        t= min(dists)
        i = dists.index(t)
        return self.df.index[i], np.degrees(t)

    def curvature(self, setup=False):
        """Curvature
        
        Distribution of the curvature per source, equivalent to the beta parameter for a LogParabola spectral model.
        """
        if setup:
            #expect to be called initially
            self.df['curvature']= np.array([model.curvature() for model in self.df.model])
            return
        assert 'curvature' in self.df, 'Curvature not calculated'
        df = self.df
        psr = np.asarray([n.startswith('PSR') for n in df.index], bool)
        fig,ax = plt.subplots(figsize=(8,6))
        hkw = dict(bins=np.linspace(0,2,41), log=True, histtype='step', lw=2)
        ax.hist(df.curvature.clip(0,2), label='all sources', **hkw)
        ax.hist(df[df.psr].curvature.clip(0,2), label='EC model', **hkw)
        ax.hist(df[psr].curvature.clip(0,2), label='PSR souce', **hkw)
        plt.setp(ax, xlabel='Curvature', ylim=(0.5,None))
        ax.legend()
        ax.grid()
        return fig
    
    def roi_check(self):
        """Distance from ROI center
        The distribution of source positions from the center of the ROI.
        %(roi_check_html)s
        """
        df = self.df
        df['actual_roi'] = map(Band(12).index, df.skydir)
        df['roi'] = map( lambda n:int(n[-4:]), df.roiname)
        df['rname'] = np.array(['HP12_{:04d}'.format(x) for x in df.roi])
        df['roi_dist'] = map(lambda s,r: np.degrees(s.difference(Band(12).dir(r))),df.skydir,df.roi)
        check = df.roi!=df.actual_roi; 
        print 'Found {} sources in the wrong ROI'.format(sum(check))
        fig,ax = plt.subplots(figsize=(5,5))
        hist_kw = dict(bins=np.linspace(0,8,33), histtype ='stepfilled', log=True)
        ax.hist(df.roi_dist.clip(0,8), **hist_kw)
        if sum(check)>0:
            ax.hist(df.roi_dist[check].clip(0.8), color='red', label='wrong ROI', **hist_kw)
        ax.grid(True, alpha=0.5);
        ax.legend()
        plt.setp(ax, xlabel='Distance (deg)', ylim=(0.8,None))
        to_move = df[check]['roi actual_roi ts'.split()].sort_index(by='roi')
        to_move['roi_dist'] = df.roi_dist
        self.roi_check_html = ''
        if sum(check)>0:
            self.roi_check_html = html_table(to_move, name=self.plotfolder+'/outside_roi', 
                heading='<h4>%d Sources that are outside the HEALPix ROI boundary</h4>' %len(to_move),
                float_format=FloatFormat(2))
        return fig

    def extended_table(self):
        """Extended table
        Table of information for fits, including links to SEDs for all extended sources.
        %(extended_table_html)s
        """
        
        df = self.df
        ext = df.isextended; sum(ext)

        cols = 'ra dec ts fitqual pindex roiname'.split()
        extdf = pd.DataFrame(df[ext][cols]).sort();
        config = configuration.Configuration('.', quiet=True, postpone=True)
        ecat = extended.ExtendedCatalog(os.path.expandvars('$FERMI/catalog/')+config.extended)
        extdf['spatial_model'] = [ecat[name].dmodel.name for name in extdf.index]
        self.extended_table_html = html_table(extdf, name=self.plotfolder+'/extended_sources',                                        
                                        heading='<h4>Table of {} extended sources</h4>'.format(len(extdf)),
                                        float_format=FloatFormat(2))
        return None
        
        
    
    def all_plots(self):
        version = os.path.split(os.getcwd())[-1]
        plt.close('all')
        csvfile='sources_%s.csv' % version
        colstosave="""ra dec ts modelname freebits fitqual e0 flux flux_unc pindex pindex_unc index2 index2_unc
                 cutoff cutoff_unc eflux100 eflux100_unc locqual delta_ts a b ang flags roiname""".split()
        self.df.ix[self.df.ts>10][colstosave].to_csv(csvfile)
        print 'saved truncated csv version to "%s"' %csvfile
        
        self.runfigures([self.census, self.cumulative_ts, 
            self.fit_quality,self.spectral_fit_consistency_plots, self.poor_fit_positions,
            self.non_psr_spectral_plots, 
            #self.beta_check, 
            self.pulsar_spectra, self.curvature, self.pivot_vs_e0, self.roi_check, self.extended_table, ]
        )
    
class OldName(object):
    """convert  from post-3FGL to P86Y naming """
    def __init__(self, filename='../../P301_6years/uw965/rename.csv'):
        """load the conversion table, make 'newname' the index"""
        self.d = pd.read_csv(filename, index_col=-1)
    def __call__(self, newname):
        if newname not in self.d.index: return None
        return self.d.ix[newname]['name']
        
class Cat_3fgl(object):
    """get 3FGL info by oldname"""
    def __init__(self, catname='3FGL-v13r3_v6r9p1_3lacv12p1_v7.fits',):
        if catname[0]!='/':
            catname = os.path.expandvars('$FERMI/catalog/'+catname)
        assert os.path.exists(catname), 'Did not find file %s' %catname
        self.ft = ft = pyfits.open(catname)[1].data
        def nickfix(n):
            return n if n[:3]!='PSR' else 'PSR '+n[3:]
        self.index_3fgl = map(nickfix, [x.strip() for x in ft.NickName_3FGL]) #Source_Name 
        id_prob = [np.nan]*len(ft)
        try:
            id_prob = ft.ID_Probability_v6r9p1[:,0] ## should find that suffix
        except: 
            print 'warning: id_prob not set' 

        self.get_oldname= OldName()
        self.cat = pd.DataFrame(dict(name3=ft.Source_Name_3FGL_1, 
                nickname=map(nickfix, ft.NickName_3FGL), 
                ra=ft.RAJ2000,dec= ft.DEJ2000, 
                ts=ft.Test_Statistic, 
                #skydir=cat_skydirs,
                #glat=glat, glon=glon, 
                #pivot=ft.Pivot_Energy, flux=ft.Flux_Density, 
                #modelname=ft.SpectrumType, 
                eflux = ft.Energy_Flux100,
                id_prob=id_prob,
                a95=ft.Conf_95_SemiMajor, b95=ft.Conf_95_SemiMinor, ang95=ft.Conf_95_PosAng,
                flags=np.asarray(ft.Flags_3FGL, int),
                ), 
            columns = 'name3 nickname ra dec ts eflux a95 b95 ang95 id_prob flags'.split(), # this to order them
            index=self.index_3fgl, )

    def __call__(self, newname):
        oldname = self.get_oldname(newname)
        if oldname not in self.index_3fgl: return None
        return self.cat.ix[oldname]

class ExtSourceInfo(SourceInfo):
    """ subclass invoked with a specific path
    """
    def __init__(self, model_path, quiet=True):
        curdir = os.getcwd()
        try:
            os.chdir(model_path)
            if not quiet: print os.getcwd()
            self.setup(refresh=False, quiet=quiet)
            self.plotfolder=model_path+'/plots/'+self.plotfolder
        finally:
            os.chdir(curdir)


