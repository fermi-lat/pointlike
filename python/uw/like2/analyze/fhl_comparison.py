"""
Comparison with a FHL model

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/fhl_comparison.py,v 1.1 2016/10/28 20:48:14 burnett Exp $

"""

import os, pickle, glob
import numpy as np
import pylab as plt
import matplotlib.gridspec as gridspec
import pandas as pd

#from uw.utilities import makepivot
from uw.like import Models
from uw.like2 import sources
from uw.like2.plotting import sed
from . import  sourceinfo, fermi_catalog
from . analysis_base import html_table, FloatFormat

class FHLcomparison(sourceinfo.SourceInfo):
    """Comparison with 3FHL
    This analysis uses the 3FHL catalog version %(fhl_version)s.
    <p>This is using the %(skymodel)s model, with many more sources, and using the same 7-year data set, 
    with Source class events. There are some differences:
        <ul>
<li>The zenith cut is 100 degrees, for all energies, while 3FHL has it at 105. this loses about 3%%
<li>It restricts theta<66.4 degrees, since the IRF is not reliable above this: also about 3%% loss
<li>It uses Front/Back event types. Perhaps losing little potential localization resolution.
<li>It uses energies from 100 MeV to 1 TeV, while 3FHL is 10 GeV to 2 TeV
<li>It uses smoothed distributions for the galactic and isotropic backgrounds, with no further adjustment
 for each ROI. However, above 10 GeV, there is little effect.
</ul>
    """
    def setup(self, pattern='gll*10GeV',uwts10_min=4, angle_cut = 0.1, **kwargs):

        super(FHLcomparison, self).setup(**kwargs )
        self.plotfolder = 'fhl_comparison'        
        self.angle_cut=angle_cut
        self.uwts10_min=uwts10_min

        # Select subset of current model with ts10>4
        # add info on E>10 GeV
        self.df['r95'] = 2.50*(self.df.a * self.df.b) ** 0.5
        self.df['ts10'] = [sum(sedrec['ts'][8:]) for sedrec in self.df.sedrec]
        self.df['np10'] = [sum(sedrec['npred'][8:]) for sedrec in self.df.sedrec]
        dfuw=self.dfuw = self.df[self.df.ts10>4]
        print 'Selected {} UW sources with TS10>{}'.format(len(self.dfuw), uwts10_min)
 
        fgl = fermi_catalog.GLL_PSC(pattern)
        df=self.df = fgl.data_frame()
        self.fhl_version= fgl.version


        def differences(a,b):
            matrix = np.array([[x.difference(y) for y in b] for x in a], np.float32)
            return matrix
        def closest(t):
            n = t.argmin()
            return (n, (np.degrees(t[n])*3600).round())

        diff_array =differences(df.skydir, dfuw.skydir)
        self.cl_fhl=cl_fhl = np.array([closest(r) for r in diff_array[:,]], int)
        self.cl_uw=cl_uw = np.array([closest(c) for c in diff_array[:,].T], int) 

        # add correlation info to the GC dataframe

        close_cut = 3600*angle_cut
        df['uw_id'] = cl_fhl[:,0]
        df['dist'] = cl_fhl[:,1]
        df['uw_ts'] = [dfuw.ix[i].ts for i in cl_fhl[:,0]]
        df['uw_ts10'] = [dfuw.ix[i].ts10 for i in cl_fhl[:,0]]
        df['uw_np10'] = [dfuw.ix[i].np10 for i in cl_fhl[:,0]]
        df['uw_r95'] = [dfuw.ix[i].r95 for i in cl_fhl[:,0]]
        df['uw_pindex'] = [dfuw.ix[i].pindex for i in cl_fhl[:,0]]
        df['uw_roi'] = [int((dfuw.ix[i].roiname)[-4:]) for i in cl_fhl[:,0]]
        df['uw_name'] = [dfuw.index[i] for i in cl_fhl[:,0]]
        df['uwok'] = (df.dist<close_cut) | np.isnan(np.array(df.r95,float))
        df['nouw'] = np.logical_not(df.uwok)
        print '3FHL sources associated with uw: {}/{}'.format(sum(df.uwok), len(df))
        corr_file='UW-3FHL_correspondence.pkl'
        print 'Wrote out file with correlation info: {}'.format(corr_file)
        df.to_pickle(corr_file)

        #add correlation info to the UW dataframe

        dfuw['fhl_id'] = np.array(cl_uw[:,0])
        dfuw['dist'] = np.array(cl_uw[:,1])

        dfuw['fhl_ts'] = [df.ix[i].ts for i in cl_uw[:,0]]
        dfuw['fhl_r95'] = [df.ix[i].r95 for i in cl_uw[:,0]]
        dfuw['fhl_pindex'] = [df.ix[i].pindex for i in cl_uw[:,0]]
        dfuw['fhl_ok'] = dfuw.dist<close_cut
        dfuw['nofhl'] = np.logical_not(dfuw.fhl_ok)

        print 'UW sources associated with 3FHL: {}/{}'.format(sum(dfuw.fhl_ok), len(dfuw))

        # load pickled SED info for corresponding sources, if found
        ff = glob.glob('3FHL_correspondence/*.pkl'); 
        print 'Found {} pickle files in 3FHL_correspondence'.format(len(ff))
        if len(ff)==0:
            print 'NO SED info found: expected to find files in "3FHL_correspondendce"'
            return
        pp = [pickle.load(open(f)) for f in ff] 
        def chisq(source):
            return sum(source.sedrec['pull'][8:]**2)
        for gtsrc in pp:
            uwname = df.ix[gtsrc.name]['uw_name']
            uw_row =dfuw.ix[uwname]
            gtsrc.uwsrc = sources.PointSource(name=uwname, skydir=uw_row['skydir'],model= uw_row['model'] )
            gtsrc.uwsrc.sedrec=uw_row['sedrec']
            gtsrc.uwsrc.np10 = uw_row['np10']
            gtsrc.chisq = chisq(gtsrc)
            gtsrc.uw_chisq = chisq(gtsrc.uwsrc)
            gtsrc.uw_np10 = uw_row['np10']
            gtsrc.npred = df.ix[gtsrc.name].npred # from 3FHL 
        ppdict = dict()
        self.gtdict = dict()
        for p in pp:
            uwname=p.uwsrc.name
            ppdict[p.name]= dict(gtchisq=p.chisq, uwname=uwname, uwchisq=p.uw_chisq, 
                                uw_roi=int(dfuw.ix[uwname]['roiname'][-4:]),
                                gtnpred=p.npred, uwnpred=dfuw.ix[uwname]['np10'],

                                )
            self.gtdict[p.name]= p
        self.ppdf = pd.DataFrame(ppdict).T

    def correlation_plots(self):
        """correlation plots
        We assume that the UW sources correponding to 3FHL sources must have detected photons above 10 GeV. 
        We use the TS for the energy bands above 
        10 GeV, called TS10 below, as a measure. Out the ~11K sources with TS>10, %(numuw)d satisfy this.
        <br>Left: histogram of closested distance
        <br>Right: Venn diagram, showing number in common, selected from the 3FHL closest
        """
        from matplotlib_venn import venn2
        df=self.df; dfuw=self.dfuw
        

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

    def comparison_plots(self):
        """Comparison plots for corresponding sources

        <br>Upper Left: Test Statistic Comparison; the UW value is for the full energy range, so is nearly always greater.
        <br>Upper right: Test Statistic comparison: for the UW fits, it is "TS10", the contribution to TS from 
        energies above 10 GeV. The 
        <br>Center left: Localization radius comparison. The UW one is almost always better since it has more data
        <br>Center right: Spectral index comparison. The lack of correlation reflects the fact that the UW Value 
        is for the entire energy range
        <br>Lower left: Compares the Npred values of the UW fit above 10 GeV, with the corresponding
         Npred from 3FHL. The dashed line, at 0.93, is the relative selection efficiency
        <br>Lower right; This plot, comparing the UW quality for E>10GeV with the corresponding 3FHL fit to the UW data,
         uses the results of a special run, loading each ROI and computing
         the SED using the 3FHL model for each corresponding source. 

        """
        skymodel=self.skymodel
        df=self.df; dfuw=self.dfuw
        dfok = df[df.uwok]

        fig, axx = plt.subplots(3,2, figsize=(15,20))

        ax = axx[0,0]
        ax.loglog(dfok.ts,dfok.uw_ts, '+b');
        ax.plot([10,1e6], [10,1e6], '--g');
        plt.setp(ax, xlabel='3FHL TS', ylabel='UW TS', xlim=(10,1e5),ylim=(10,1e5))
        ax.grid(alpha=0.5)

        ax = axx[0,1]
        lim = (10,1e4)
        x,y = dfok.ts.clip(*lim),dfok.uw_ts10.clip(*lim)
        ax.semilogx(x, y/x, '+b');
        ax.plot(lim, lim, '--g');
        plt.setp(ax, xlabel='3FHL TS', ylabel='UW TS10 / 3FHL TS', xlim=lim, ylim=(0.5,1.5))
        ax.grid(alpha=0.5)

        ax = axx[1,0]
        lim = (2e-3,0.2)
        ax.loglog((dfok.r95).clip(*lim), (dfok.uw_r95).clip(*lim), '+b');
        ax.plot(lim, lim, '--g')
        plt.setp(ax, xlim=lim, ylim=lim, xlabel='3FHL R95', ylabel='UW R95');
        ax.grid(alpha=0.5);

        ax = axx[1,1]
        lim=(1.0,3.5)
        ax.plot(dfok.pindex.clip(*lim), dfok.uw_pindex.clip(*lim), '+b')
        plt.setp(ax, xlim=lim, ylim=lim, xlabel='3FHL Spectral Index', ylabel='UW pindex')
        ax.plot(lim, lim, '--g')
        ax.grid(alpha=0.5)

        ax=axx[2,0]
        gtnp = self.ppdf.gtnpred
        uwnp = self.ppdf.uwnpred
        ax.semilogx(gtnp, uwnp/gtnp, '.')
        ax.axhline(0.93, color='r', ls='--', lw=2)
        lim=(4,1000)
        plt.setp(ax, xlim=lim, ylim=(0.5,1.5), xlabel='3FHL Npred', ylabel='UW Npred / 3FHL Npred')
        ax.grid()

        ax= axx[2,1]
        lim = (0,25); limx=(0,25.2)
        gtc = self.ppdf.gtchisq
        uwc = self.ppdf.uwchisq
        poor = gtc>25
        ax.plot(gtc.clip(*lim), uwc.clip(*lim),'.', label='all')
        ax.plot(gtc.clip(*lim)[poor], uwc.clip(*lim)[poor],'or', label='Poor 3FHL fits')
        plt.setp(ax, xlabel='3FHL quality', xlim=limx,
                ylabel='UW quality', ylim=limx);
        ax.grid(); ax.legend()

        fig.suptitle('Comparison of values for common sources', fontsize=14);
        fig.set_facecolor('white')
        return fig

    def missing_from_fhl(self,min_ts=25, max_index=2.5):
        """Plots of UW sources not in FHL
       As expected, amost all UW sources with TS10 above a threshold, apparently 30, are matched. 
       The exceptions, the red histograms, need to be understood: 
       the onset represents a measure of the detection efficiency, including the 3FHL TS>25 requirement, and excluding extended sources.
       <br> %(notinfgl_html)s
        """
        df=self.df; dfuw=self.dfuw
        msg=  'UW sources not in 3FHL with TS10>{} :{}'.format(
            min_ts, sum( dfuw.nofhl & (dfuw.ts10>min_ts) ))
        msg+=  '<br>Largest UW TS10: {:.0f}'.format(max(dfuw.ts10[dfuw.nofhl]))
        self.notinfgl_html=msg
        dfuw['aname']=aname =[a['name'][0] if a is not None else None for a in dfuw.associations]
        dfuw['roi'] = [int(n[-4:]) for n in dfuw.roiname]
        dfnogt = dfuw[(dfuw.nofhl)&(dfuw.ts10>25)\
             & (np.logical_not(dfuw.isextended))\
             & (dfuw.pindex<max_index)]\
            ['ra dec glon glat ts ts10 pindex r95 locqual roi aname'.split()]
        self.notinfgl_html += html_table(dfnogt, float_format=FloatFormat(2),
                name=self.plotfolder+'/notinfhl',
                heading='<h4>List of {} {} sources with ts>{} and spec. index <{} not in 3FHL</h4>'.format(len(dfnogt), self.skymodel, min_ts, max_index),
                )

        # save a csv file as well
        outfile ='missing_from_3FHL.csv' 
        dfnogt.to_csv(os.path.join(self.plotfolder,outfile))
        print 'Wrote file {}'.format(self.plotfolder+outfile)
        self.notinfgl_html += '<br><a href="{}">CSV file of missing sources</a>'.format(outfile)

        # the plots
        fig, axx = plt.subplots(1,3, figsize=(15,6))
        ax=axx[0]
        histkw = dict(bins=np.logspace(1,4,33), histtype='step', log=True)
        ax.hist(dfuw.ts10.clip(10,1e4), label='all uw', **histkw);
        ax.hist(dfuw.ts10[dfuw.fhl_ok].clip(10,1e4), label='match', **histkw);
        ax.axvline(25, color='red', ls='--')
        ax.grid(True, alpha=0.5)
        ax.legend()
        plt.setp(ax, xscale='log', ylim=(0.8,None), xlabel='UW TS10');

        ax=axx[1]
        histkw = dict(bins=np.linspace(0,50,26), histtype='step',log=True, lw=2)
        ax.hist(dfuw.ts10.clip(0,50), label='all uw', **histkw);
        ax.hist(dfuw.ts10[dfuw.fhl_ok].clip(0,50),label='match' , **histkw);
        ax.axvline(25, color='red', ls='--')
        ax.hist(dfuw.ts10[np.logical_not(dfuw.fhl_ok)].clip(0,50),label='not 3FHL', **histkw);
        ax.grid(True, alpha=0.5)
        ax.legend()
        plt.setp(ax, xscale='linear', ylim=(0.8,None), xlabel='UW TS10');

        ax=axx[2]
        dfuw['singlat']=np.sin(np.radians(np.array(dfuw.glat,float)))
        histkw=dict(bins=np.linspace(-1,1,41),histtype='step', lw=2)
        
        tscut = dfuw.ts10>min_ts
        ax.hist(dfuw.singlat[(dfuw.fhl_ok) & tscut], label='match',  color='green', **histkw);
        histkw.update(histtype='stepfilled')
        ax.hist(dfuw.singlat[np.logical_not(dfuw.fhl_ok) & tscut], label='not 3FHL', color='red', **histkw);
        ax.text(0.05,0.9, 'TS10>{}'.format(min_ts), transform=ax.transAxes)
        plt.setp(ax, xlabel='sin(b)')
        ax.grid(alpha=0.5)
        ax.legend()

        fig.set_facecolor('white')
        return fig
    
    def missing_from_uw(self):
        """FHL Sources not in UW
        Table of FHL source not found in the UW list
        %(notinuw_html)s
        """
        dfuw=self.dfuw
        df=self.df
        dfnouw=df[df.nouw]['ra dec ts pindex r95 dist uw_r95 uw_ts uw_ts10 uw_name uw_roi'.split()]
        msg= html_table(dfnouw,float_format=FloatFormat(2), href=False,
            name=self.plotfolder+'/missing_uw',
            heading='<h4>{} FGL sources not found in {} </h4>'.format(len(dfnouw),self.skymodel)
            )
        self.notinuw_html=msg
        return None

    def poor_spectral_fits(self, qual_min=25):
        """Poor spectral fits
        This table, and corresponding SED plots, are for the corresponding 3FHL sources that have a poor
        fit using the UW model above 10 GeV. (The red dots in the quality factor comparison plot) 
        <p>Note that the fit quality should be adjusted for very strong sources to account for systematic errors 
        in flux measurements. A few of the above show changes in the spectrum, such that a log parabola fit to the
         full energy range cannot work above 10 GeV.
        <br> %(bad_fits_html)s
        """
        badfit = self.ppdf.gtchisq>qual_min
        baddf =self.ppdf[badfit]; 
        baddf.index.name='gtname'
        #baddf.to_csv('bad_spectral_check.csv')
        self.bad_fits_html = html_table(baddf,href=False, float_format=FloatFormat(1),
            name=self.plotfolder+'/poor_spectral_fit',
            heading='<h4>{} 3FHL sources with poor fits using {} model'.format(len(baddf), self.skymodel),
            )
        fig, axx = plt.subplots(4,4, figsize=(20,20), sharex=True,sharey=True)
        for name, ax in zip(baddf.index, axx.flatten()):
            sed.plot_other_source(self.gtdict[name].uwsrc, self.gtdict[name], ax=ax, emin=10000);    
        return fig

    def all_plots(self):
        self.runfigures([
            self.correlation_plots, self.comparison_plots, 
            self.missing_from_fhl, self.missing_from_uw, self.poor_spectral_fits,
            ])