"""
Comparison of a UW model with a GC analysis
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/gc_comparison.py,v 1.3 2017/11/17 22:43:17 burnett Exp $

"""

import os
from astropy.io import fits as pyfits
import numpy as np
import pylab as plt
import pandas as pd
from astropy.io import fits

from skymaps import SkyDir
# FIX LATER from matplotlib_venn import venn2
from . import sourceinfo
from . analysis_base import FloatFormat, html_table

_filenames=dict(
    GC300='GCPass8_listsources_cleaned07deg_TS10_tsmap_300MeVref_fit3TS10.fits',
    GC500='GCPass8_listsources_all_TS10_tsmap_500MeV_cleaned02deg_fit3TS10_ref.fits',
)

class GCcomparison(sourceinfo.SourceInfo):
    """Comparison with a GC analysis

    %(setup_info)s
    """

    def setup(self, gc_root = '/nfs/farm/g/glast/u/burnett/analysis/GC_check/',
        cols = 'ra dec glat glon ts pos_sigma Spectral_Index'.split(),
        gcname = 'GC300', **kw):

        super(GCcomparison, self).setup(**kw)
        self.plotfolder='comparison_%s' % gcname
        filename = _filenames[gcname]

        # load the GC analysis
        assert os.path.exists(gc_root + filename)
        hdulist =fits.open(gc_root + filename)
        data = hdulist[1].data
        df = pd.DataFrame([data.field(x) for x in cols], index=cols).T
        df.index=data.field('Source_name'); 
        df.index.name='name'
        df['sd'] = map(SkyDir, df.ra,df.dec)
        print ('Created a DataFrame from {} with {} rows'.format(filename, len(df)))
        self.dfgc = df
        self.gcname=gcname
        self.setup_info='Comparing sources in file {} with {}'\
            .format(filename, self.skymodel)
        self.select_uw()
    
    def gc_plots(self):
        """GC plots

        Using all entries
        """
        df = self.dfgc
        fig,axx =plt.subplots(2,2, figsize=(12,12))
        ax = axx[0,0]
        ax.hist(df.ts, np.logspace(1,5,41), histtype='step');
        plt.setp(ax, xscale='log', xlabel='ts');
        ax.grid(alpha=0.5)
        ax=axx[0,1]
        df.glon[df.glon>180]-=360
        ax.plot(df.glon, df.glat, '.')
        plt.setp(ax, xlabel='glon', ylabel='glat', xlim=(25,-25), ylim=(-25,25))
        ax.grid(alpha=0.5)
        ax=axx[1,0]
        ax.loglog(df.ts, df.pos_sigma, '.')
        plt.setp(ax, xlabel='ts', ylabel='pos_sigma');
        ax.grid(alpha=0.5)
        ax=axx[1,1]
        ax.hist(df.Spectral_Index, np.linspace(1,4,31), histtype='step');
        plt.setp(ax, xlabel='Spectral Index')
        ax.grid(alpha=0.5)
        fig.suptitle('{} plots'.format(self.gcname), fontsize=16)
        return fig

    def select_uw(self, angle_cut=21):
        dfuwall = self.df
        glon=dfuwall.glon
        glon[glon>180]-=360
        glat=dfuwall.glat
        dfuw = dfuwall[(abs(glon)<angle_cut) & (abs(glat)<angle_cut)]
        dfuw['sd'] = map(SkyDir, dfuw.ra, dfuw.dec)
        dfuw['pos_sigma'] = (dfuw.a * dfuw.b) ** 0.5
        text = 'Selected {} uw sources within box of {} deg about GC'.format(len(dfuw),angle_cut)
        print (text);
        self.setup_info +='<br>'+text
        self.dfuw = dfuw

    def correlate(self, angle_cut=0.1):
        """Correlate the two lists
        """
        def differences(a,b):
            matrix = np.array([[x.difference(y) for y in b] for x in a], np.float32)
            return matrix
        def closest(t):
            n = t.argmin()
            return (n, (np.degrees(t[n])*3600).round())

        diff_array =differences(self.dfgc.sd, self.dfuw.sd)
        cl_gc = np.array([closest(r) for r in diff_array[:,]], int)
        cl_uw = np.array([closest(c) for c in diff_array[:,].T], int) 

        # add correlation info to the GC dataframe
        close_cut = 3600*angle_cut
        df = self.dfgc
        dfuw = self.dfuw
        df['uw_id'] = cl_gc[:,0]
        df['dist'] = cl_gc[:,1]
        df['uw_ts'] = [dfuw.ix[i].ts for i in cl_gc[:,0]]
        df['uw_pos_sigma'] = [dfuw.ix[i].pos_sigma for i in cl_gc[:,0]]
        df['uw_pindex'] = [dfuw.ix[i].pindex for i in cl_gc[:,0]]
        df['uw_roi'] = [int((dfuw.ix[i].roiname)[-4:]) for i in cl_gc[:,0]]
        df['uwok'] = df.dist<close_cut
        df['nouw'] = np.logical_not(df.uwok)
        print ('GC sources associated with uw: {}/{}'.format(sum(df.uwok), len(df)))

        #add correlation info to the UW dataframe
        dfuw['gc_id'] = np.array(cl_uw[:,0])
        dfuw['dist'] = np.array(cl_uw[:,1])

        dfuw['gc_ts'] = [df.ix[i].ts for i in cl_uw[:,0]]
        dfuw['gc_pos_sigma'] = [df.ix[i].pos_sigma for i in cl_uw[:,0]]
        dfuw['gc_pindex'] = [df.ix[i].Spectral_Index for i in cl_uw[:,0]]
        dfuw['gcok'] = dfuw.dist<close_cut
        dfuw['nogc'] = np.logical_not(dfuw.gcok)
        print ('UW sources associated with GC: {}/{}'.format(sum(dfuw.gcok), len(dfuw)))
        
        fig, axx =plt.subplots(1,2, figsize=(12,6))
        ax=axx[0]
        hist_kw=dict(bins=np.linspace(0,0.5, 26), histtype='step', log=True)
        ax.hist(cl_gc[:,1].clip(0,1800)/3600., 
            label='{} [{}]'.format(self.gcname,len(cl_gc)), **hist_kw);
        ax.hist(cl_uw[:,1].clip(0,1800)/3600., 
            label='{} [{}]'.format(self.skymodel,len(cl_uw)), color='orange',**hist_kw);
        ax.axvline(angle_cut, color='red', ls='--', label='angle cut')
        plt.setp(ax, ylim=(0.8,None), xlabel='distance (deg)', title='minimum distance')
        ax.legend(loc='upper left');
        ax.grid(alpha=0.5)
    
        # make a Venn diagram
        common = sum(cl_uw[:,1]<close_cut)
        venn2(ax=axx[1], subsets=(len(self.dfuw)-common, len(self.dfgc)-common,common)
        ,set_labels=(self.skymodel,self.gcname))
        return fig

    def common_plots(self):
        """Common plots
        
        """
        df = self.dfgc
        fig, axx = plt.subplots(2,2, figsize=(12,12))
        ax = axx[0,0]
        dfok = df[df.uwok]
        ax.loglog(dfok.ts,dfok.uw_ts, '+b');
        ax.plot([10,1e6], [10,1e6], '--g');
        plt.setp(ax, xlabel='GC ts', ylabel='UW ts')
        ax.grid(alpha=0.5)
        ax = axx[0,1]
        lim = (1e-3,0.1)
        ax.loglog(dfok.pos_sigma.clip(*lim), dfok.uw_pos_sigma.clip(*lim), '+b');
        ax.plot(lim, lim, '--g')
        plt.setp(ax, xlim=lim, ylim=lim, xlabel='GC pos_sigma', ylabel='UW pos_sigma');
        ax.grid(alpha=0.5);
        ax = axx[1,0]
        lim=(1.5,3.5)
        ax.plot(dfok.Spectral_Index.clip(*lim), dfok.uw_pindex.clip(*lim), '+b')

        plt.setp(ax, xlim=lim, ylim=lim, xlabel='GC Spectral Index', ylabel='UW pindex')
        ax.plot(lim, lim, '--g')
        ax.grid(alpha=0.5)
        axx[1,1].set_visible(False)
        fig.suptitle('Comparison of values for common sources', fontsize=14)
        return fig

    def nouwcounterpart(self):
        """Plots for GC sources with no UW counterpart
        """
        df=self.dfgc
        fig,axx =plt.subplots(2,2, figsize=(12,12))
        ax = axx[0,0]
        ax.hist(df.ts, np.logspace(1,5,21), histtype='step', label='all');
        ax.hist(df[df.nouw].ts, np.logspace(1,5,21), histtype='step', color='red', label='no UW');
        plt.setp(ax, xscale='log', xlabel='ts');
        ax.legend()
        ax.grid(alpha=0.5)
        ax=axx[0,1]
        df.glon[df.glon>180]-=360
        ax.plot(df.glon, df.glat, '.')
        ax.plot(df[df.nouw].glon, df[df.nouw].glat, '.r')
        plt.setp(ax, xlabel='glon', ylabel='glat', xlim=(25,-25), ylim=(-25,25))
        ax.grid(alpha=0.5)
        ax=axx[1,0]
        ax.plot(df.glon, df.glat, '.')
        ax.plot(df[df.nouw].glon, df[df.nouw].glat, '.r')
        plt.setp(ax, xlabel='glon', ylabel='glat', xlim=(25,-25), ylim=(-5,5))
        ax.grid(alpha=0.5)
        ax=axx[1,1]
        ax.hist(df.Spectral_Index, np.linspace(1,4,31), histtype='step');
        ax.hist(df[df.nouw].Spectral_Index, np.linspace(1,4,31), histtype='step', color='red');
        ax.grid(alpha=0.5)
        plt.setp(ax, xlabel='Specral Index')
        fig.suptitle('{} distributions with no {} counterpart'.format(self.gcname, self.skymodel), fontsize=16); 
        return fig   

    def nogccounterpart(self):
        """Plots for UW sources with no GC counterpart
        """
        dfuw=self.dfuw
        fig,axx =plt.subplots(2,2, figsize=(12,12))
        ax = axx[0,0]
        ax.hist(dfuw.ts, np.logspace(1,5,41), histtype='step');
        ax.hist(dfuw[dfuw.nogc].ts, np.logspace(1,5,41), histtype='step', color='red', label='no GC');
        plt.setp(ax, xscale='log', xlabel='ts');
        ax.legend()
        ax.grid(alpha=0.5)
        ax=axx[0,1]
        dfuw.glon[dfuw.glon>180]-=360
        ax.plot(dfuw.glon, dfuw.glat, '.')
        ax.plot(dfuw[dfuw.nogc].glon, dfuw[dfuw.nogc].glat, '.r')
        plt.setp(ax, xlabel='glon', ylabel='glat', xlim=(25,-25), ylim=(-25,25))
        ax.grid(alpha=0.5)
        ax=axx[1,0]
        ax.plot(dfuw.glon, dfuw.glat, '.')
        ax.plot(dfuw[dfuw.nogc].glon, dfuw[dfuw.nogc].glat, '.r')
        plt.setp(ax, xlabel='glon', ylabel='glat', xlim=(25,-25), ylim=(-5,5))
        ax.grid(alpha=0.5)
        ax=axx[1,1]
        ax.hist(dfuw.pindex, np.linspace(1,4,31), histtype='step');
        ax.hist(dfuw[dfuw.nogc].pindex, np.linspace(1,4,31), histtype='step', color='red');
        ax.grid(alpha=0.5)
        plt.setp(ax, xlabel='Spectral Index')
        fig.suptitle('{} sources with no {} counterpart'.format(self.skymodel, self.gcname),
            fontsize=16);    
        return fig

    def gc_subset_table(self, tsmin=25, bmin=0,
        cols='ra dec glon glat ts pos_sigma Spectral_Index dist uw_roi'.split()
        ):
        """Table of selected strong GC sources not in the UW list
        
        %(selected_sources)s
        """
        df = self.dfgc
        selection = df.nouw & (df.ts>tsmin) & (np.abs(df.glat)>bmin)
        dfsel = df[selection][cols]
        # replace all prefixes
        split_names = [n.split()[1] for n in dfsel.index]
        dfsel.index = [self.gcname+' '+sn for sn in split_names]
        dfsel.index.name='name'
        self.selected_sources=\
            html_table(
                dfsel,
                 href=False,
                 heading = '<p>%d GC sources with TS>%.0f, |b|>%.1f' % (sum(selection),tsmin, bmin),
                    name=self.plotfolder+'/selected', maxlines=20,
                    float_format=(FloatFormat(3))
            )
        outfile = self.plotfolder+'/selected.csv'
        dfsel.to_csv(outfile)
        print ('Wrote {} seeds to file {}'.format(sum(selection), outfile))

    def all_plots(self):
        """
        %(setup_info)s
        """
        self.runfigures([ 
                self.gc_plots, self.correlate, self.common_plots, self.nouwcounterpart, self.nogccounterpart,self.gc_subset_table,
        ])
        
