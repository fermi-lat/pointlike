"""
Basic analyis of source spectra

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/sourceinfo.py,v 1.6 2013/07/12 13:37:17 burnett Exp $

"""

import os, pickle
from collections import Counter
import numpy as np
import pylab as plt
import pandas as pd

from uw.utilities import makepivot
from . import diagnostics
from . diagnostics import FloatFormat
from . _html import HTMLindex


class SourceInfo(diagnostics.Diagnostics):
    """ To be superclass for specific source plot stuff, creates or loads
        a DataFrame with all sources 
        """
    require='pickle.zip'
    def setup(self, **kwargs):
        self.plotfolder='sources' #needed by superclass
        filename = 'sources.pickle'
        refresh = kwargs.pop('refresh', not os.path.exists(filename) or os.path.getmtime(filename)<os.path.getmtime('pickle.zip'))
        if refresh:
            files, pkls = self.load_pickles('pickle')
            assert len(files)==1728, 'Expected to find 1728 files'
            sdict= dict()
            for pkl in pkls:
                for name, info in pkl['sources'].items():
                    model = info['model']
                    pars = np.empty(4); pars.fill(np.nan)
                    errs = np.empty(4); errs.fill(-2)
                    free = np.zeros(4, bool)
                    n = model.len()
                    pars[:n] = model.parameters
                    free[:n] = model.free
                    try:
                        d = np.diag(model.get_cov_matrix())
                        d[d<0] =0
                        errs[:n] = np.sqrt(d)
                        errs[np.isnan(errs)]=-1
                        badfit = np.any(errs[model.free]<=0)
                    except Exception, msg:
                        print 'fail errors for %s:%s' % (name, msg)
                        badfit = True
                    try:
                        fitqual = round(sum(info['sedrec'].delta_ts),2)
                    except:
                        fitqual = np.nan
                    ellipse = info.get('ellipse', None)
                    sdict[name] = info
                    pulsar = model.name.endswith('Cutoff')
                    betavalue = float(pars[2]) if not pulsar else np.nan
                    if pulsar: # factor to convert flux to prefactor
                        bvalue =1.0 if model.npar==3 else model['b']
                        prefactor = np.exp(-(model.e0/model['cutoff'])**bvalue)
                    else: prefactor = 1.0
                    sdict[name].update(
                        glat=info['skydir'].b(), glon=info['skydir'].l(),
                        roiname=pkl['name'], 
                        pars= pars, errs=errs, free=free, badfit=badfit,
                        a = ellipse[2] if ellipse is not None else np.nan,
                        b = ellipse[3] if ellipse is not None else np.nan,
                        ang=ellipse[4] if ellipse is not None else np.nan,
                        locqual = round(ellipse[5],2) if ellipse is not None else np.nan,
                        delta_ts = ellipse[6] if ellipse is not None else np.nan,
                        freebits= np.sum( int(b)*2**i for i,b in enumerate(model.free)),
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
                        eflux = prefactor*pars[0]*model.e0**2*1e6,
                        psr = pulsar,
                        )
            df = pd.DataFrame(sdict).transpose()
            df.index.name='name'
            df['name'] = df.index
            ra = [x.ra() for x in df.skydir]
            dec = [x.dec() for x in df.skydir]
            df['ra'] = ra
            df['dec'] = dec
            self.df = df.sort_index(by='ra')
            self.df.save(filename)
            print 'saved %s' % filename

        else:
            print 'loading %s' % filename
            self.df = pd.load(filename)
        #self.df['flux']    = [v[0] for v in self.df.pars.values]
        #self.df['flux_unc']= [v[0] for v in self.df.errs.values]
        localized = ~np.array(pd.isnull(self.df.delta_ts))
        extended = np.array(self.df.isextended, bool)
        self.df['unloc'] = ~(localized | extended)
        self.df['poorloc'] = (self.df.a>0.2) + (self.df.locqual>8) + (self.df.delta_ts>2)
        self.df['flags'] = 0  #used to set bits below
        flags = self.df.flags
        pl = (self.df.poorloc + self.df.unloc) * (self.df.ts>10)
        flags[pl] += 8 ### bit 8
        #print '%d sources flagged (8) as poorly or not localized' % sum(pl)

 
        self.energy = np.sqrt( self.df.ix[0]['sedrec'].elow * self.df.ix[0]['sedrec'].ehigh )
            
    def skyplot(self, values, proj=None, ax=None, ecliptic=False,
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
        """
        assert hasattr(values, 'index'), 'skyplot: values arg must have index attribute'
        
        # generate arrays of glon and singlat using index 
        sd = self.df.ix[values.index, ['glat', 'glon']] # see page 101
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

    def cumulative_ts(self, ts=None, other_ts=None, tscut=(10,25), check_localized=True, label=None, otherlabel=None):
        """ Cumulative test statistic TS
        
        A logN-logS plot, but using TS. Important thresholds at TS=10 and 25 are shown.
        """
        usets = self.df.ts if ts is None else ts
        df = self.df
        fig,ax = plt.subplots( figsize=(8,6))
        dom = np.logspace(np.log10(9),5,1601)
        ax.axvline(25, color='gray', lw=1)
        ax.hist( usets ,dom, cumulative=-1, lw=2, color='g', histtype='step',label=label)
        if other_ts is not None:
            ax.hist( other_ts ,dom, cumulative=-1, lw=2, color='b', histtype='step',label=otherlabel)
        if check_localized:
            unloc = df.unloc
            ul = df[(unloc+df.poorloc) * usets>tscut[0]] 
            n = len(ul)
            if n>10:
                ax.hist(ul.ts ,dom, cumulative=-1, lw=2, color='r', histtype='step',
                    label='none or poor localization')
                ax.text(12, n, 'none or poor localization (TS>%d) :%d'%(tscut[0],n), fontsize=12, color='r')
        plt.setp(ax,  ylabel='# sources with greater TS', xlabel='TS',
            xscale='log', yscale='log', xlim=(9, 1e4), ylim=(9,8000))
        ax.set_xticklabels([' ', '10', '100', '1000'])
        ax.set_yticklabels(['', '10', '100', '1000'])
            
        # label the plot with number at given TS
        for t in tscut:
            n = sum(usets>t) 
            ax.plot([t,2*t], [n,n], '-k');
            ax.plot(t, n, 'og')
            ax.text(2*t, n, 'TS>%d: %d'%(t,n), fontsize=14, va='center')
                
        ax.grid()
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
        """ Plots showing spectral parameters for PowerLaw and ExpCutoff spectra
        Left: energy flux in eV/cm**2/s. This is the differential flux at the pivot energy
        <br> Center: the spectral index.
        <br> Right: the curvature index for the subset with log parabola fits.
        %(tail_check)s
        %(beta_check)s
        """
        fig, axx = plt.subplots( 1,3, figsize=(12,4))
        plt.subplots_adjust(wspace=0.2, left=0.05,bottom=0.15)

        t = self.df.ix[(self.df.ts>10)*(self.df.modelname=='LogParabola')]['ts flux pindex beta beta_unc freebits e0 roiname'.split()]
        t['eflux'] = t.flux * t.e0**2 * 1e6
        ax = axx[0]
        [ax.hist(t.eflux[t.ts>tscut].clip(4e-2,1e2), np.logspace(-2,2,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xscale='log', xlabel='energy flux', xlim=(4e-2,1e2)); ax.grid(); ax.legend(prop=dict(size=10))
        ax = axx[1]
        [ax.hist(t.pindex[t.ts>tscut].clip(index_min,index_max), np.linspace(index_min,index_max,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xlabel='spectral index'); ax.grid(); ax.legend(prop=dict(size=10))
        ax = axx[2]
        [ax.hist(t.beta[(t.ts>tscut)*(t.beta>0.01)].clip(0,beta_max), np.linspace(0,beta_max,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xlabel='beta'); ax.grid(); ax.legend(prop=dict(size=10))
        # get tails
        tail_cut = (t.eflux<5e-2)+((t.pindex<index_min)+(t.pindex>index_max))*t.beta.isnull()+(t.beta>beta_max)
        
        if sum(tail_cut)>0:
            tails=t[tail_cut]['ts eflux pindex beta freebits roiname'.split()].sort_index(by='roiname')
            filename = 'non_pulsar_tails.html'
            html_file = self.plotfolder+'/%s' % filename
            #html = tails.sort_index(by='roiname').to_html(float_format=FloatFormat(2))
            html = diagnostics.html_table(tails, float_format=FloatFormat(2))
            open(html_file,'w').write('<head>\n'+ HTMLindex.style + '</head>\n<body>'+ html+'\n</body>')
            self.tail_check = '<p><a href="%s?skipDecoration">Table of %d sources on tails</a>: '% (filename, len(tails))
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
        beta_bad = (t.beta>0.001) * ((t.beta_unc==0) + (t.beta/t.beta_unc<2) + (t.freebits!=7))
        if sum(beta_bad)>0:
            print '%d sources fail beta check' % sum(beta_bad)
            self.beta_check ='<br>Sources failing beta 2-sigma significance check' +\
            diagnostics.html_table(t[beta_bad]['ts beta beta_unc freebits roiname'.split()], float_format=FloatFormat(2))
            
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
        psrmodel = (self.df.ts>10)*(self.df.modelname=='PLSuperExpCutoff')
        t = self.df.ix[psrmodel]\
            ['ts flux pindex cutoff e0 index2 index2_unc roiname freebits fitqual'.split()]
        t['eflux'] = t.flux * t.e0**2 * 1e6
        badfit = t.fitqual>30

        def plot1(ax, efmin=1e-1,efmax=1e3):
            bins = np.logspace(np.log10(efmin),np.log10(efmax),26)
            vals = t.eflux.clip(efmin,efmax)
            ax.hist(vals, bins )
            ax.hist(vals[badfit], bins, color='red', label='poor fit')
            plt.setp(ax, xscale='log', xlabel='energy flux', xlim=(efmin,efmax)); ax.grid(); 
            ax.legend(prop=dict(size=10))

        def plot2(ax):
            bins = np.linspace(index_min,index_max,26)
            vals = t.pindex.clip(index_min,index_max)
            ax.hist(vals, bins)
            ax.hist(vals[badfit], bins, color='red', label='poor fit')
            plt.setp(ax, xlabel='spectral index'); ax.grid(); 
            ax.legend(prop=dict(size=10))
            
        def plot3(ax):
            bins = np.linspace(0,cutoff_max/1e3,26)
            vals = t.cutoff.clip(0,cutoff_max) /1e3
            ax.hist(vals, bins)
            ax.hist(vals[badfit], bins, color='red', label='poor fit')
            plt.setp(ax, xlabel='cutoff energy (GeV)'); ax.grid()
            ax.legend(prop=dict(size=10))
            
        def plot4(ax, xlim=(0,cutoff_max)):
            xvals = t.cutoff.clip(*xlim) / 1e3
            yvals = t.pindex.clip(index_min,index_max)
            ax.plot(xvals, yvals, 'o')
            ax.plot(xvals[badfit], yvals[badfit], 'or', label='poor fit')
            plt.setp(ax, xlabel='cutoff energy', ylabel='spectral index',
                xlim=(xlim[0]-0.1,1.03*xlim[1]/1e3), ylim=(index_min-0.1, index_max+0.1),
                )
            ax.grid(); ax.legend(loc='lower right', prop=dict(size=10))

        for f,ax in zip((plot1,plot2,plot3,plot4,), axx.flatten()): f(ax)
        flags = self.df.flags
        tail_cut = (t.pindex<=index_min) + (t.pindex>index_max) + (t.cutoff>cutoff_max)
        tails = t.ix[tail_cut].index
        flags[tails] += 1 ### bit 1
        print '%d pulsar sources flagged (1) in tails of  index or cutoff' % sum(tail_cut)
        if sum(tail_cut)>0:
            tails=t[tail_cut]['ts eflux pindex cutoff freebits roiname'.split()]
            filename = 'pulsar_tails.html'
            html_file = self.plotfolder+'/%s' % filename
            #html = tails.sort_index(by='roiname').to_html(float_format=FloatFormat(2))
            html = diagnostics.html_table(tails.sort_index(by='roiname'), float_format=FloatFormat(2))
            open(html_file,'w').write('<head>\n'+ HTMLindex.style + '</head>\n<body>'+ html+'\n</body>')
            self.pulsar_tail_check = '<p><a href="%s?skipDecoration">Table of %d sources on tails</a>: '% (filename, len(tails))
            self.pulsar_tail_check += 'Criteria: require index between 0 and 2.5, cutoff<8 GeV'
        else:
            self.pulsar_tail_check ='<p>No sources on tails'

        # table of pulsars with b<1

        filename='pulsar_b.html'
        #t['significance']=np.where(t.index2_unc>0, (1-t.index2)/t.index2_unc, [np.nan]*len(t) )
        tt=t[t.index2<1]['ts fitqual pindex cutoff index2 index2_unc'.split()]
        tt['significance'] = (1-tt.index2)/tt.index2_unc
        html_file = self.plotfolder+'/%s' % filename
        html = diagnostics.html_table(tt,float_format=FloatFormat(2))
        open(html_file,'w').write('<head>\n'+ HTMLindex.style + '</head>\n<body>'+ html+'\n</body>')
        self.pulsar_b = '<p><a href="%s?skipDecoration">Table of %d sources with b&lt;1</a> '% (filename, len(tt))
        print '%d pulsar sources with b<1' %len(tt)

        # table of fits with any fixed parameter other than b
        tt = t[(t.freebits&7!=7)]['ts fitqual pindex cutoff freebits roiname'.split()].sort_index(by='roiname')
        if len(tt)>0:
            print '%d pulsar-like sources with fixed parameters' %len(tt)
            self.pulsar_fixed='<p>Sources with any fixed parameter other than b: %s' % diagnostics.html_table(tt, float_format=FloatFormat(2))
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
        
    def fit_quality(self, xlim=(0,50), ndf=10, tsbandcut=20):
        """ Spectral fit quality
        This is the difference between the TS from the fits in the individual energy bands, and that for the spectral fit.
        It should be distributed approximately as chi squared of at most 14-2 =12 degrees of freedom. 
        However, high energy bins usually do not contribute, so we compare with ndf=%(ndf)d.
        All sources with TS_bands>%(tsbandcut)d are shown.<br>
        Left: Power-law fits. Tails in this distribution perhaps could be improved by converting to log parabola. 
        <br>Center: Log parabola fits.
        <br>Right: Fits for the pulsars, showing high latitude subset.
        <br> Averages: %(fit_quality_average)s
        %(badfit_check)s
        %(poorfit_table)s

        """
        from scipy import stats
        fig, axx = plt.subplots(1,3, figsize=(12,6))
        plt.subplots_adjust(left=0.1)
        s = self.df
        psr = np.asarray(s.psr, bool)
        beta = s.beta
        logparabola = (~psr) * (beta>0.01)
        powerlaw = (~psr) * (beta.isnull() + (beta<0.01) )

        self.tsbandcut=tsbandcut
        cut=s.band_ts>tsbandcut
        
        dom = np.linspace(xlim[0],xlim[1],26)
        d = np.linspace(xlim[0],xlim[1],51); delta=dom[1]-dom[0]
        chi2 = lambda x: stats.chi2.pdf(x,ndf)
        fudge = 1.0 # to scale, not sure why
        hilat = np.abs(self.df.glat)>5
        self.average = [0]*4; i=0
        for ax, label in zip(axx[:2], ('powerlaw', 'logparabola',)):
            mycut=cut*eval(label)
            count = sum(mycut)
            ax.hist(s.fitqual[mycut].clip(*xlim), dom, label=label+' (%d)'%count)
            self.average[i]=s.fitqual[mycut].mean(); i+=1
            ax.plot(d, chi2(d)*count*delta/fudge, 'r', lw=2, label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
            ax.grid(); ax.set_xlabel('fit quality')
            ax.legend(prop=dict(size=10))
            
        def right(ax, label='PSR'):
            mycut = cut * (psr)
            count = sum(mycut)
            ax.hist(s.fitqual[mycut].clip(*xlim), dom, label=label+' (%d)' %count)
            ax.hist(s.fitqual[mycut*hilat].clip(*xlim), dom, label=label+' [|b|>5] (%d)' %sum(mycut*hilat))
            self.average[i]=s.fitqual[mycut*hilat].mean()
            self.average[i+1]=s.fitqual[mycut*(~hilat)].mean()
            ax.plot(d, chi2(d)*count*delta/fudge, 'r', lw=2, label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
            ax.grid();ax.set_xlabel('fit quality')
            ax.legend(loc='upper left', prop=dict(size=10))
        
        right(axx[2])
        self.df['badfit2'] =np.array(self.df.badfit.values, bool)
        t = self.df.ix[(self.df.badfit2)*(self.df.ts>10)].sort_index(by='roiname')
        print '%d sources with bad fits' %len(t)
        if len(t)>0:
            self.badfit = t[['ts', 'errs', 'roiname']]
            #self.badfit_check = '<h4>Sources with missing errors:</h4>'+self.badfit.to_html(float_format=FloatFormat(1))
            self.badfit_check = '<h4>Sources with missing errors:</h4>'+diagnostics.html_table(self.badfit, float_format=FloatFormat(1))
        else: self.badfit_check = '<p>All sources fit ok.'
        self.fit_quality_average =  ', '.join( map(lambda x,n :'%s: %.1f' %(n,x) ,
                            self.average, 'powerlaw logparabola expcutoff(hilat) expcutoff(lolat)'.split()) )
        self.ndf=ndf
        print 'fit quality averages:', self.fit_quality_average

        # Make tables (csv and html) of the poor fits
        s['pull0'] = np.array([x.pull[0] for x in s.sedrec])
        t =s.ix[((s.fitqual>30) | (np.abs(s.pull0)>3))*(s.ts>10) ]['ra dec glat fitqual pull0 ts modelname freebits index2 roiname'.split()].sort_index(by='roiname')
        poorfit_csv = 'poor_spectral_fits.csv'
        t.to_csv(poorfit_csv)
        bs =sorted(list(set(t.roiname)))
        print 'Wrote out list of poor fits to %s, %d with fitqual>30 or abs(pull0)>3, in %d ROIs' % (poorfit_csv, len(t), len(bs))
        # todo: make a function to do this nidcely
        poorfit_html = self.plotfolder+'/poorfits.html'
        #t_html = '<h3>Table of poorly-fit sources, model %s</h3>'%self.skymodel + t.to_html(float_format=FloatFormat(2),
        #        formatters=dict(ra=FloatFormat(3), dec=FloatFormat(3), ts=FloatFormat(0),index2=FloatFormat(3)))
        t_html = '<h3>Table of poorly-fit sources, model %s</h3>'%self.skymodel + diagnostics.html_table(t,float_format=FloatFormat(2),
                formatters=dict(ra=FloatFormat(3), dec=FloatFormat(3), ts=FloatFormat(0),index2=FloatFormat(3)))

        open(poorfit_html,'w').write('<head>\n'+ HTMLindex.style + '</head>\n<body>'+t_html+'\n</body>')
        self.poorfit_table = '<p> <a href="poorfits.html?skipDecoration"> Table of %d poor fits, with fitqual>30 or abs(pull0)>3</a>' % (  len(t) )
        # flag sources that made it into the list
        self.df.flags[t.index] |= 2
        print '%d sources flagged (2) as poor fits' %len(t)
        return fig
      
    def poor_fit_positions(self):
        """ Positions of poorly-fit sources
        Selection: fitqual>30 or |pull|>3
        """
        s = self.df
        s['pull0'] = np.array([x.pull[0] for x in s.sedrec])
        poor = ( (s.fitqual>30) | (np.abs(s.pull0)>3))*(s.ts>10) 
        return self.skyplot(s.fitqual[poor], vmin=30, vmax=100)
        
    
    def pivot_vs_e0(self, xylim=(100, 4e4)):
        """ pivot vs e0
        The reference energy, e0, is fixed except by a special run that iterates until the measured pivot energy, 
        which is the energy at which the differential flux uncertainty is minimum is the same. This plot checks that by measuring the pivot energy, and comparing it with the current reference. Note that e0 is required to be in the range 200 MeV to 20 GeV.
        """
        fig, ax = plt.subplots(figsize=(4,4))
        s = self.df
        cut = s.ts>10
        ax.plot(s.e0[cut].clip(*xylim), s.pivot_energy[cut].clip(*xylim), '.')
        plt.setp(ax, xscale='log',xlabel='e0', xlim=xylim, 
                    ylabel='pivot', yscale='log', ylim=xylim)
        ax.set_title('compare calculated pivot with e0', fontsize=10)
        ax.grid()
        return fig
        
    def fitquality(self):
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
        ax.grid(); ax.legend(prop=dict(size=10))
        ax = axs[1]
        ax.plot(s.ts, fitqual, '.'); 
        plt.setp(ax, xscale='log', xlabel='TS', xlim=(10,1e5),
             ylabel='fit qual',ylim=(1,1e3),yscale='log')
        ax.grid()
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
            title = 'low energy fit consistency'
        ):
        """ Spectral fit consistency for the lowest energy bin
        
        These plots show the consistency of the lowest energy band with the spectrum
        defined by the full fit. <br>
        Left: distribution of the "pull" <br>
        Center: data/model ratio with errors, vs. the model flux.<br>
        Right: position in the sky of flagged sources <br>
        """
        cut=self.df.ts>25
        s = self.df[cut]
        fdata = np.array([s.ix[i]['sedrec'].flux[0] for i in range(len(s))])
        udata = np.array([s.ix[i]['sedrec'].uflux[0] for i in range(len(s))])
        ldata = np.array([s.ix[i]['sedrec'].lflux[0] for i in range(len(s))])
        pull = np.array([s.ix[i]['sedrec'].pull[0] for i in range(len(s))])
        fmodel = np.array([s.ix[i]['model'](energy)*energy**2*1e6 for i in range(len(s))])
        glat = np.array([x.b() for x in s.skydir])
        fluxcut = fmodel>minflux
        latcut  = abs(glat)>5.0
        hilat = fluxcut*(latcut)
        lolat = fluxcut*(~latcut)
        
        lowebad = np.abs(pull)>3
        self.df.flags[lowebad] += 4
        print 'Tagged %d sources with lowebad bit (4)' % sum(lowebad)

        y = fdata/fmodel
        ylower, yupper =[(fdata-ldata)/fmodel,(udata-fdata)/fmodel]
        xhi,yhi,yerrhi = fmodel[hilat], y[hilat], [ylower[hilat],yupper[hilat]]
        xlo,ylo,yerrlo = fmodel[lolat], y[lolat], [ylower[lolat],yupper[lolat]]
        
        def error_bar(ax):
            ax.errorbar(x=xhi, y=yhi, yerr=yerrhi, fmt='og', label='%d hilat sources'%sum(hilat))
            ax.errorbar(x=xlo, y=ylo, yerr=yerrlo, fmt='or', label='%d lowlat sources'%sum(lolat))
            plt.setp(ax, xlabel=r'$\mathsf{model\ flux\ (eV\ cm^{-2} s^{-1}})$', xscale='log', 
                ylabel='data/model', ylim=(0,2.5), xlim=(minflux, 100) )
            ax.set_xticks([2,5,10,20,50,100])
            ax.set_title( title, fontsize='medium')
            ax.set_title( title, fontsize='medium')
            ax.legend(loc='upper left', prop=dict(size=10))
            ax.grid()  

        def hist(ax):
            dom = np.linspace(-3,3,26)
            hist_kw=dict(lw=2, histtype='step')
            q=pull.clip(-3,3) 
            
            ax.hist(q[hilat], dom, color='g',  label='%d hilat sources'%sum(hilat),  **hist_kw)
            ax.hist(q[lolat], dom, color='r',  label='%d lowlat sources'%sum(lolat), **hist_kw)
            ax.set_xlabel('pull')
            ax.axvline(0, color='k')
            ax.set_xlim((-3,3))
            ax.set_title( title, fontsize='medium')
            ax.legend(loc='upper left', prop=dict(size=10))
            ax.grid()  

        def skyplot(ax):
            pdf = pd.DataFrame(dict(pull=pull), index=s.index) # to atatch indx
            self.skyplot(pdf.pull[lowebad], ax=ax, vmin=-3, vmax=3, title=title, cbtext='pull')

        fig,ax = plt.subplots(1,3, figsize=(12,4))
        plt.subplots_adjust(wspace=0.3)
        for f, ax in zip( (hist, error_bar, skyplot), ax.flatten()):
            f(ax=ax)
        return fig
        
    def census(self, primary_prefix='P7R4'):
        """Census
        
        %(census_html)s
        <p>
        In the first table of prefixes, the columns are the number of sources with TS greater than the header value. 
        The row labels are the first four characters of the source name, except 'ext' means extended.
        <br>The second table has the suffixes; for sources with the "%(primary_prefix)s" prefix. Each group, except 'X-Z', 
        represents a set of seeds added to the original
        list of sources. The last row, 'X-Z', were added by hand.
        """
        df = self.df
        extended = np.asarray(df.isextended.values,bool)
        pointsource = ~extended

        def count(prefix, tsmin):
            if prefix=='ext':
                return sum(extended*(df.ts>tsmin))
            elif prefix=='total':
                return sum(df.ts>tsmin)
            names = df[pointsource*(df.ts>tsmin)]['name'].values    
            return sum([n.startswith(prefix) for n in names])
        prefixes = list(set( n[:4] for n in df[pointsource]['name'])) +['ext', 'total']
        census = dict()
        for x in (0, 10, 25):
            census[x] = [count(prefix, x) for prefix in prefixes]
        self.census_data=pd.DataFrame(census, index=prefixes)
        self.census_html = '\n<h4>Prefixes</h4>\n'+diagnostics.html_table(self.census_data)
        
        # now check suffixes
        self.primary_prefix=primary_prefix
        suffixes=[s[-1] for s in self.df.index.values if s.startswith(primary_prefix) and s[-1]>'9']
        c=Counter(suffixes)
        scounts = lambda  r : int(sum([c[x] for x in r if x in c.keys()]))
        suffixranges = ('ABCDEF', 'GHIJKL', 'MNO', 'PQRSTUVW', 'XYZ')
        sdict = dict([(r[0]+'-'+r[-1], [scounts(r)]) for r in suffixranges])
        self.census_html += '\n<h4>Suffixes</h4>\n'+diagnostics.html_table(pd.DataFrame(sdict, index=['freq']).T)

    def all_plots(self):
        """ Plots of source properties, from analysis of spectral fits. 
        See <a href="../localization/index.html?skipDecoration"> localization </a> for localization plots.
        """
        version = os.path.split(os.getcwd())[-1]
        plt.close('all')
        csvfile='sources_%s.csv' % version
        colstosave="""ra dec ts modelname freebits fitqual e0 flux flux_unc pindex pindex_unc index2 index2_unc
                 cutoff cutoff_unc locqual delta_ts a b ang flags roiname""".split()
        self.df.ix[self.df.ts>10][colstosave].to_csv(csvfile)
        print 'saved truncated csv version to "%s"' %csvfile
        
        self.runfigures([self.census, self.cumulative_ts, self.fit_quality,self.spectral_fit_consistency_plots, self.poor_fit_positions,
            self.non_psr_spectral_plots, self.pulsar_spectra, self.pivot_vs_e0, self.flag_proc, ]
        )

    def flag_proc(self):
        """ Flagged source summary:
        %(flagged_link)s
        """
        # Generate summary table for flagged sources
        t =self.df[(self.df.flags>0)*(self.df.ts>10)]['ra dec ts fitqual pull0 eflux pindex beta cutoff index2 flags roiname'.split()]
        t.to_csv('flagged_sources.csv')
        print 'wrote %d sources to flagged_sources.csv' % len(t)
        
        num=[sum(self.df.flags & 2**b > 0) for b in range(4)] 
        flagtable=pd.DataFrame(dict(number=num, description=('tails','poor fits','low energy bad', 'poor localization') ))
        flagtable.index.name='bit'
        self.flagged_link = """\
        <p>A number of these sources have been flagged to indicate potential issues. 
        The flag bits and number flagged as such are:
        %s<br>  """ % diagnostics.html_table(flagtable, href=False)
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

        