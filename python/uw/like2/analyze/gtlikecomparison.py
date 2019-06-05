"""
Comparison with a gtlike model

"""

import os, glob
import cPickle as pickle
import numpy as np
import pylab as plt
import matplotlib.gridspec as gridspec
import pandas as pd

from skymaps import SkyDir
#from uw.utilities import makepivot
from uw.like import Models
from . import (sourcecomparison, sourceinfo,fermi_catalog,)
from . analysis_base import html_table, FloatFormat


class FL8YComparison(sourceinfo.SourceInfo):
    """Comparison with 4FGL
            This analysis uses the %(cat)s catalog file %(fhl_file)s.
    <p>This is using the %(skymodel)s model, with the same 8-year data set, 
    with Source class events. There are some differences:
    <ul>
    <li>The zenith cut is 100 degrees, for all energies, while %(cat)s varies from 90 to 110.
    <li>It restricts theta<66.4 degrees, since the IRF is not reliable above this: about 3%% loss
    <li>It uses Front/Back event types, with Front only for E<316 MeV. This loses some localization resolution, but avoids the PSF3 effective area problem.
    <li>The diffuse models are not modified for each ROI.    </ul>
    """
    def setup(self, pattern='gll_pscP305uw8606', **kwargs):
        super(FL8YComparison, self).setup(**kwargs)
        self.cat='4FGL'
        self.plotfolder='{}_comparison'.format(self.cat)

        # make copy dataframe with compressed names
        self.old_index = self.df.index
        cindex = [n.replace(' ','') for n in self.df.index]
        self.df.index = cindex
        # add info on E>10 GeV
        systematic = self.config['localization_systematics']
        f95, quad = 2.45*systematic[0], systematic[1]/60. 
        self.df['r95'] = (f95**2*(self.df.a * self.df.b) + quad**2)** 0.5

        # get the catalog "gll" entries as a DataFrame and set corresponding values
        if pattern is None:
            pattern=self.config['gllcat']
        filename = glob.glob(os.path.expandvars('$FERMI/catalog/'+pattern))[-1]
        fcat = fermi_catalog.GLL_PSC2(filename)
        self.fhl_file = fcat.filename.split('/')[-1]
        self.gdf = gdf=  fcat.df
        gdf['uw_ts']    = self.df.ts
        gdf['uw_r95']   = self.df.r95
        gdf['uw_pindex']= self.df.pindex
        gdf['uw_eflux100']=self.df.eflux100

        # identify sources missing from FL8Y
        # 
        a = set(cindex)
        b = set(self.gdf.index); 
        lost = np.array(list(set(b.difference(a))))
        if len(lost)>10:
            print '{} {} sources not here:,{}...'.format(len(lost), self.cat, lost[:10])
        self.lost=lost # save for further analysis

        # add boolean for in FL8Y
        self.df['fl8y'] = [n in gdf.index for n in cindex]

    def source_info_plots(self, tt, tscut=100):

        sd = map (SkyDir, tt.ra, tt.dec)
        glon = np.array(map(lambda x: x.l(), sd))
        glon[glon>180]-=360
        glat = map(lambda x: x.b(), sd)
        singlat = np.sin(np.radians(glat))
        hights = tt.ts>tscut

        fig, axx = plt.subplots(1,2, figsize=(12,6))
        plt.subplots_adjust(wspace=0.25)

        ax = axx[0]
        ts=np.array(tt.ts, float)
        hkw=dict(bins=np.logspace(np.log10(20),3,41), log=True, lw=2)
        ax.hist(ts[~pd.isnull(ts)].clip(10,1e3),histtype='step', **hkw);
        ax.hist(ts[~pd.isnull(ts) & hights].clip(10,1e3), color='red',
          label='TS>{}: {}\nmax:{:.0f}'.format(tscut, sum(hights), tt.ts.max()),histtype='stepfilled', **hkw);
        ax.legend()
        ax.set(xscale='log', xlabel='TS', ylim=(0.9, None));  
        ax.axvline(25, color='green', ls='--')  
        
        ax = axx[1]
        self.basic_skyplot(ax, glon,singlat, 'blue',  s=10,  title='Locations')
        self.basic_skyplot(ax, glon[hights],singlat[hights],'red', s=30,  title='Locations')

        return fig

    def extra_source_info(self, **kwargs):
        """Info about additional sources
        
        TS and locations of the %(extra_count)s sources <b>not</b> in 4FGL.
        """
        gnames = set(self.gdf.index)
        pnames = set(self.df.index)
        gdiff = gnames.difference(pnames)
        pdiff = pnames.difference(gnames)
        self.extra_count = len(pdiff)
        extra = self.df.loc[list(pdiff)] 

        fig = self.source_info_plots(extra, **kwargs)

        return fig

    def lost_source_info(self, **kwargs):
        """Info on lost sources

        Left: TS values
        Right: locations, showing the high TS values
        
        <p> Link to a csv file containing a list of the %(number_lost)s sources that were lost:
        <a href="../../%(lost_sources)s?download=true">%(lost_sources)s</a>
        These are entries for which the NickName field does not have a corresponding source.
        """

        tt = self.gdf.query('~(uw_ts>0)')

        fig = self.source_info_plots(tt, **kwargs)
        
        lost_name = '{}/lost_sources.csv'.format(self.plotfolder)
        tt.index.name='name'
        self.number_lost = len(tt)
        tt['sname ra dec ts pindex eflux100 r95'.split()].to_csv(lost_name)
        print 'Wrote file "{}" with info on {} missing sources'.format(lost_name, self.number_lost)
        self.lost_sources = lost_name
   
        return fig

    def load_pickled_info(self, path='psc_check/info', debug=False):
        # if hasattr(self, 'ts_df'):
        #     return self.ts_df

        # get the TS and chisq values
        ff =sorted(glob.glob(path+'/*'))
        print 'read {} pickle files from {}'.format(len(ff), path)
        dd = map(lambda f:pickle.load(open(f)), ff)
        z = dict()
        gtmodel=dict()
        for roi,d in enumerate(dd):
            for a,b in d:
                try:
                    eflux_pt=a[1].i_flux(e_weight=1)*1e6;
                    eflux_gt=b[1].i_flux(e_weight=1)*1e6
                except Exception, msg:
                    print b[0],msg                    
                    eflux_pt=eflux_gt=np.nan
                z[b[0]] = dict(
                    ts_pt=a[2],        ts_gt=b[2], 
                    chisq_pt=a[3],     chisq_gt=b[3], 
                    eflux_pt=eflux_pt, eflux_gt=eflux_gt,
                    nickname=a[0], roi=roi,
                    index_pt=a[1][1],
                    )
                gtmodel[a[0]]=b[1] 

        q = self.ts_df=pd.DataFrame(z).T
        self.gtmodel =gtmodel
        if debug:
            return q
        # add positional info, using nickname field as a key into the model dataframe (which has compressed names)
        nicknames = map(lambda n:n.replace(' ',''), self.ts_df.nickname.values)

        # check for now missing nicknames
        indexset= set(self.df.index); 
        nicknameset = set(nicknames)
        missing_nicknames = list(nicknameset.difference(indexset))
        if len(missing_nicknames)>0:
            print 'Warning: following nicknames not in current model: {}'.format(np.array(missing_nicknames))
            nicknames = list(indexset.intersection(nicknameset))
            cnick = [n.replace(' ','') for n in nicknames]
            qv = [n.replace(' ','') for n in q.nickname.values];
            ok = np.array([name in cnick for name in qv], bool)
            q = self.ts_df = q[ok]
        if debug:
            return nicknames    
        sdir = self.df.loc[nicknames,'skydir'].values
        if debug: return sdir
        glon = np.array(map(lambda s:s.l(), sdir),float)
        glon[glon>180]-=360
        glat = map(lambda s:s.b(), sdir)
        singlat = np.sin(np.radians(glat))
        q['glon']=glon; q['glat']=glat

        # construct quality difference
        a,b = q.chisq_pt.values, q.chisq_gt.values
        for x in a,b:
            x[pd.isna(x)]=100
        delta = np.array( b-a, float  )#/ np.array(q.ts_pt,float)**power, float)
        q['delta']=delta

        # flux ratio
        q['eflux_ratio'] = q.eflux_pt/q.eflux_gt
        return q

    def comparison_plots(self):
        """Comparison plots for corresponding sources

        <br>Upper Left: Test Statistic Comparison; the UW value is for the full energy range, so is nearly always greater.
        <br>Center left: Localization radius comparison. The UW one is almost always better since it has more data
        <br>Center right: Spectral index comparison. This is always the index used for the fit spectral model
 
        """
        skymodel=self.skymodel
        df=self.gdf; dfuw=self.df
        psr = np.array([name.startswith('PSR') for name in df.index],bool); 
        dfok = df

        def cplot(ax, a,b, xlim, label, ylim=(0.,2.),xscale='log'):
            ax.semilogx(a.clip(*xlim), (b/a).clip(*ylim), '.b');
            ax.semilogx(a[psr].clip(*xlim), (b/a)[psr].clip(*ylim), '.r', label='PSR');
            ax.axhline(1.0, ls='--', color='g');
            ax.set( xlabel=label, ylabel='UW/gtlike ratio', xlim=xlim,
                ylim=ylim, xscale=xscale)
            ax.set_xlabel(label,fontsize=14)
            ax.grid(alpha=0.5)
            ax.legend()
            #ax.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
            #ax.set_yticklabels(['0.6','0.8','1.0', '1.2', '1.4', '1.6'])

        fig, axx = plt.subplots(4,1, figsize=(12,15))
        plt.subplots_adjust(left=0.05, top = 0.95, hspace=0.3    )

        cplot(axx[0], dfok.ts, dfok.uw_ts, (20,1e5), 'TS')
        cplot(axx[1], dfok.r95,dfok.uw_r95,(8e-3,0.5),'R95 [deg]')
        cplot(axx[2], dfok.eflux100, dfok.uw_eflux100/1.602e-6, (4e-7, 1e-3),
            'eflux100 [erg/(s cm^2)]')
        cplot(axx[3], dfok.pindex, dfok.uw_pindex, (1.0, 3.5),'pindex', xscale='linear')
        fig.suptitle('Comparison of values for common sources', fontsize=14);
        fig.set_facecolor('white')
        return fig

    def setup_quality(self):
        #takes some time
        self.info =  self.load_pickled_info()

    def quality_check_plots(self, ylim=(-5,25), tsmin=100):
        """Fit quality check
        Compare fit consistency values of this model with that for the %(cat)s model, that is, 
        as caculated with the pointlike implementation, 
        but using the %(cat)s spectra determined by gtlike.
        <br><b>Upper Left:</b> Scatter plot of the TS value difference, normalized by the square root of the uw value
        <br><b>Upper Right:</b> Histogram of the normalized TS difference, for TS_uw>100.
        <br><b>Lower Left: </b> Positions of sources in each tail
        <br><b>Lower Right: </b> Positions of sources in each tail, along gal. plane
        <b><h3>%(deltax2_positive)s</h3>
        <b><h3>%(deltax2_positive_psr)s</h3>
        <b><h3>%(deltax2_negative)s</h3>
        """
        if not hasattr(self, 'info'):
            self.setup_quality()
        q =   self.info

        delta_clip = q.delta.clip(*ylim)
        delta_label = '(chi2_uw - chi2_g)/sqrt(TS_uw)'

        # make a table of the outliers
        neg =(q.delta<=ylim[0]) & (q.ts_pt>tsmin)
        pos =(q.delta>=ylim[1]) & (q.ts_pt>tsmin)
        psr = np.array([name.startswith('PSR') for name in q.nickname])
  
        print 'Outliers (above TS={}): {} negative, {} positive'.format(tsmin, sum(neg), sum(pos))
        try:
            self.deltax2_positive=html_table(q[pos & ~psr].sort_values(by='delta', ascending=False), 
                name=self.plotfolder+'/deltax2_positive', 
                heading='<h4>pointlike better: (non-pulsars) {}</h4>'.format(sum(pos & ~psr)),
                href=True, href_pattern='psc_check/sed/%s*.jpg', )
            self.deltax2_positive_psr=html_table(q[pos & psr].sort_values(by='delta', ascending=False), 
                name=self.plotfolder+'/deltax2_positive_psr', 
                heading='<h4>pointlike better (pulsars): {}</h4>'.format(sum(pos & psr)),
                href=True, href_pattern='psc_check/sed/%s*.jpg', )
            if sum(neg)>0:
                self.deltax2_negative=html_table(q[neg].sort_values(by='delta', ascending=True), 
                    name=self.plotfolder+'/deltax2_negative', 
                    heading='<h4>gtlike better: {}</h4>'.format(sum(neg)),
                    href=True, href_pattern='psc_check/sed/%s*.jpg', )
            else:
                self.deltax2_negative = '<h4>No gtlike fits better than Delta TS = {}<.h4>'.format(np.abs(ylim[0]))
        except Exception, msg:
            print 'Failed to create tables of of outliers: "{}"'.format(msg)

        fig, axy = plt.subplots(2,2, figsize=(15,10))
        plt.subplots_adjust(wspace=0.25)
        axx = axy.flatten()

        ax = axx[0]  # a)
        ridge = np.array((abs(q.glon)<60.) & (abs(q.glat)<5), bool)
        ax.semilogx(q.ts_pt.clip(10, 1e5), delta_clip, '.b')
        ax.semilogx(q[ridge].ts_pt.clip(10, 1e5), delta_clip[ridge], '.r', label='ridge')
        ax.axhline(0, color='orange')
        ax.set(ylabel=delta_label)
        ax.legend(loc='lower right')

        ax = axx[1]  # b)
        hkw = dict(bins= np.linspace(ylim[0],ylim[1],36), histtype='step', lw=2, log=False)
        delta_clip_ts = delta_clip[q.ts_pt>100]
        ax.hist(delta_clip_ts,**hkw);
        hkw.update(histtype='stepfilled')
        ax.hist(delta_clip_ts[delta_clip_ts<=ylim[0]], color='green', **hkw)
        ax.hist(delta_clip_ts[delta_clip_ts>=ylim[1]], color='red', **hkw)
        ax.axvline(0, color='orange')
        ax.set_xlabel(delta_label)
        
        ax = axx[2]  # c)
        cut = (q.delta>=10) | (q.delta<=-2)
        singlat = np.sin(np.radians(q.glat))
        self.basic_skyplot(ax, q.glon[cut], singlat[cut],
            delta_clip[cut], s=20, cmap=plt.get_cmap('coolwarm'));

        ax = axx[3]  # d)
        self.basic_skyplot(ax, q.glon[cut], singlat[cut],
            delta_clip[cut], s=20, cmap=plt.get_cmap('coolwarm'), aspect=5*180.);
        ax.set(ylim = (-0.1,0.1))
        
        return fig

    def correlation_plots(self, df=None, dfuw=None):
        """correlation plots
        We assume that the UW sources correponding to 3FHL sources must have detected photons above 10 GeV. 
        We use the TS for the energy bands above 
        10 GeV, called TS10 below, as a measure. Out the ~11K sources with TS>10, %(numuw)d satisfy this.
        <br>Left: histogram of closested distance
        <br>Right: Venn diagram, showing number in common, selected from the 3FHL closest
        """
        from matplotlib_venn import venn2
        
        if df is None: df=self.gdf; 
        if dfuw is None: dfuw=self.df
 
        fig, axx =plt.subplots(1,2, figsize=(12,6))
        ax=axx[0]
        hist_kw=dict(bins=np.linspace(0,0.5, 26), histtype='step', log=True,lw=2)
        ax.hist(self.cl_fhl[:,1].clip(0,1800)/3600., 
            label='{} [{}]'.format('3FHL',len(self.cl_fhl)), **hist_kw);
        ax.hist(self.cl_uw[:,1].clip(0,1800)/3600., 
            label='{} [{}]'.format('uw7000',len(self.cl_uw)), color='orange',**hist_kw);
        ax.axvline(self.angle_cut, color='red', ls='--', label='angle cut')
        plt.setp(ax, ylim=(0.8,None), xlabel='distance (deg)', title='minimum distance')
        ax.legend(loc='upper left');
        ax.grid(alpha=0.5)

        # make a Venn diagram
        self.common=common = sum(df.uwok)
        self.numuw=len(dfuw)
        v=venn2(ax=axx[1], subsets=(len(dfuw)-common, len(df.uwok)-common,common)
            ,set_labels=('uw7000 with TS10>{}'.format(self.uwts10_min),'3FHL'))
        for text in v.set_labels:
            text.set_fontsize(14)
        fig.set_facecolor('white')
        return fig

    def all_plots(self):
        self.runfigures([
            self.extra_source_info,
            self.lost_source_info,
            self.comparison_plots, 
            self.quality_check_plots,
            ]
        )



class GtlikeComparison(sourcecomparison.SourceComparison):
    """Results of comparison with glike 
    Gtlike version:  %(catname)s <br>Compare the two lists of sources and spectral parameters, assuming that the skymodel names 
    correspond to the "NickName" field in the gtlike-generated FITS file.
    Assume that a special run has been made to evaluate the likelihoods for the gtlike models.

    """
    def setup(self, catpat=None, **kw):
        if catpat is None:
            catpat = self.config['gllcat']
        pat = os.path.expandvars(os.path.join('$FERMI','catalog', catpat))
        gllcats = sorted(glob.glob(pat))
        assert len(gllcats)>0, 'No gtlike catalogs found using {}'.format(pat)
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
        
    def pulsar_candidates(self, mints=16):
        """Pulsar candidate selection
        Only choose those that are NOT in the gtlike list
        %(pulsar_pivot_info)s
        """
        df = self.dfx
        cut = cut = (df.no_gtlike) & (np.abs(df.glat)>5) & (df.ts>mints)
        sedrec = self.dfx['sedrec']
        tflow = []
        tfmed = []
        for sr in sedrec:
            ts = sum(sr.ts)
            low, middle, high = [sum(sr.ts[i:i+4]/ts).round(2) for i in (0,4,8)]
            tflow.append(low)
            tfmed.append(middle)
        self.dfx['tflow']=tflow
        self.dfx['tfmed']=tfmed
        
        pcand = df[cut & (df.tfmed>0.6) & (df.locqual<5)]['ra dec glat glon ts tflow tfmed a locqual'.split()]; print len(pcand)
        pcand.index.name='name'
        pcand.to_csv('weak_pulsar_candidate.csv')
        try:
            pc=makepivot.MakeCollection('weak pulsar candidates', 'sedfig', 'weak_pulsar_candidate.csv', refresh=True)
            self.pulsar_pivot_info="""<p> These can be examined with a 
            <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
            which requires Silverlight. """% pc.cId
        except Exception, msg:
                self.pulsar_pivot_info = '<p>No pivot output; job failed %s' %msg

        fig, axx = plt.subplots(1,2, figsize=(10,4))
        
        def scat_low_med(ax):
            ax.plot(df[cut]['tflow'], df[cut]['tfmed'], '.')
            plt.setp(ax, xlabel='low TS fraction', ylabel='medium TS fraction', 
                xlim=(-0.02,1), ylim=(-0.02, 1))
         
        def hist_med(ax):
            ax.hist(df[cut]['tfmed'], np.linspace(0,1,26))
            plt.setp(ax, xlabel='medium TS fraction')
        hist_med(axx[0]); scat_low_med(axx[1])
        return fig
        
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
        <br><b>Center</b>: Photon index for sources with |b|>5 and TS>16.
        """
        df = self.dfx
        fig, ax = plt.subplots(1,3, figsize=(12,4))
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
        def plot_index(ax):
            cut = (df.no_gtlike) & (np.abs(df.glat)>5) & (df.ts>16)
            ax.hist(df[cut]['pindex'], np.linspace(1,3, 26))
            plt.setp(ax,xlabel='photon index')
            ax.grid()
            
        
        for f, ax in zip([plot1,plot_index, plot2,], ax): f(ax)
        return fig

    # def all_plots(self):
    #     self.runfigures([self.check_contents, self.missing, self.compare_fits,  
    #     self.delta_ts, self.galactic_distortion  ])


class GardianPrefactors(object):
    
    def __init__(self,
            root = '/nfs/farm/g/glast/g/diffuse/P8diffuse/results/gardian/8yr/',
            names = ['HighLatTA09', 'OuterGalaxyTA03', 'InnerGalaxyTA01']):
        
        files = [root+name+'_FinalVariables.txt' for name in names]
        regions = 'HighLat Outer Inner'.split()
        filename = files[0]
        t = open(filename).read().split('\n')
        len(t)

        pref = dict()
        all = dict()
        for region,filename in zip(regions,files):
            t = open(filename).read().split('\n')
            n=0
            for i,line in enumerate(t):
                tok = line.split()
                if len(tok)<4 or  tok[0]!='FL8Y' :continue
                sname, var = tok[1].split('_')
                if var!='pref': continue
                value = float(tok[3])
                if abs(value-1)>1e-5:
                    pref[sname] = dict(prefactor=value, region=region)
                    if sname not in all:
                        all[sname]= {region:value}
                    else:
                        all[sname].update({region:value})
                    n+=1
            print 'Found {} variable sources in region {}'.format(n,region)
        self.df =pd.DataFrame(pref).T
        self.all = all

    def add_fit_comparison(self, catpat='gll_psc*uw8011*', skymodel='$FERMI/skymodels/P305_8years/uw8503'):

        self.gtlc = FL8YComparison(skymodel, pattern=catpat)
        self.dfcmp = self.gtlc.load_pickled_info()

        self.dfcmp['prefactor'] = self.df.prefactor
        self.dfcmp['region'] = self.df.region
        #return self.dfcmp        
        #cat = fermi_catalog.GLL_PSC2(catpat)
        qq = self.gtlc.gdf
        qq['nickname']=qq.index
        qq.index= [n[5:] for n in qq.sname];
        self.df['eflux100'] = qq.eflux100*1e6 # to eV
        self.df['nickname'] = qq.nickname
    
    def prefactor_hist(self, ax=None):
        ax=plt.gca() if ax is None else ax
        x = np.array(self.df.prefactor,float)
        ax.hist(x.clip(0.8,1.2), bins=np.linspace(0.8, 1.2, 25),histtype='step', 
                lw=2, label='mean: {:.2f}'.format(x.mean()))
        ax.axvline(1.0, ls = '-', color='orange');
        ax.set(xlabel='prefactor', title='FL8Y prefactors in diffuse fits');
        ax.legend();

    def flux_hist(self, ax=None):
        df = self.df
        ax = plt.gca() if ax is None else ax
        x = np.array(df.eflux100, float)
        ax.hist(x, bins=np.logspace(1,3,20), label='all')
        ax.hist(x[df.region=='Inner'], bins=np.logspace(1,3,20), label='inner')

        ax.set(xscale='log', xlabel='Energy Flux');
        ax.legend();
        
    def prefactor_comparison(self, ax=None):
        t = self.dfcmp.query('prefactor>0');
        plt.rc('font', size=14)
        if ax is None:
            fig,ax = plt.subplots(figsize=(8,8))
        else: fig=ax.figure
        lim = (0.85,1.05)
        groups = t.groupby('region')
        for marker, (name, g) in zip('oDs', groups):
            x = np.array(g.eflux_ratio,float).clip(*lim)
            y = np.array(g.prefactor,float).clip(*lim)
            ax.plot(x,y, marker, label=name);
        ax.plot(lim, lim, '-' ,color='orange' );
        ax.set(xlabel='pointlike/FL8Y EFlux ratio', ylabel='RH01 prefactor')
        ax.legend(loc='upper left');

        import matplotlib.patches as patches
        # Add the patch to the Axes
        size=lim[1]-lim[0]
        rect = patches.Rectangle((lim[0],lim[0]),size,size,ls='--',edgecolor='grey',facecolor='none')
        ax.add_patch(rect);