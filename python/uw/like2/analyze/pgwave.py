"""
Analyze PGWAVE seeds

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/pgwave.py,v 1.1 2016/03/21 18:54:57 burnett Exp $

"""
import os, pickle
import numpy as np
import pylab as plt
import pandas as pd
from skymaps import SkyDir 
from .. import associate, seeds
from .analysis_base import FloatFormat


from . import (analysis_base, _html, sourceinfo)
from analysis_base import html_table, FloatFormat

class PGwave( sourceinfo.SourceInfo): 
    """Analysis of the PGWave seeds
    <br>
    <p>
    From file %(pgw_filename)s, month %(month)s
    """
    def setup(self, **kw):
        self.plotfolder = self.seedname='pgwave'
        self.spectral_type='power law'
        assert self.skymodel.startswith('month'), 'Current folder not a month'
        self.month=int(self.skymodel[5:]); 
        try:
            self.pgw_filename='/nfs/farm/g/glast/g/catalog/transients/TBIN_%d_all_pgw.txt'% (self.month-1)
            assert os.path.exists(self.pgw_filename), 'PGWAVE file %s not found'% pgw_filensme  

            t=self.df = seeds.read_seedfile('pgw')
            prefix = 'PGW_%02d'%self.month
        except:
            self.pgw_filename='/nfs/farm/g/glast/g/catalog/transients/P302/PGW/1m_1mp15_PGW_ALL.fits'
            t=self.df = seeds.read_seedfile('PGW')
            prefix =  'PG{:02d}'.format(self.month)
        sdirs = map(SkyDir, t.ra,t.dec)
        t['l']=[s.l() for s in sdirs]
        t.l[t.l>180]-=360
        t['b']=[s.b() for s in sdirs]
        
        # now make tables for sources that were accepted
        sdf = pickle.load(open('sources.pickle'))
        used = [n.startswith(prefix) for n in sdf.index]
        self.udf=udf = sdf[used]
        good = (udf.ts>10) & (udf.a<0.50) & (udf.locqual<8)
        print 'Used, good:',sum(used),sum(good)
        self.good=good

    def skyplot(self):
        """seed positions
        
        """
        t=self.df
        fig, axx = plt.subplots(2,2, figsize=(12,12))
        ax = axx[0,0]
        ax.plot(t.l, t.b, '.')
        plt.setp(ax, xlabel='l', ylabel='b', xlim=(180,-180), ylim=(-90,90))
        ax.axhline(0, color='k', alpha=0.3)
        ax.axvline(0, color='k', alpha=0.3)
        ax = axx[1,0]
        sinb = np.sin(np.radians(t.b))
        dom = np.linspace(-1,1,51)
        ax.hist(list(sinb),dom);
        plt.setp(ax, xlabel='sin(b)')
        ax.grid(True, alpha=0.5)
        ax = axx[0,1]
        ax.plot(t.l, t.b, '.')
        plt.setp(ax, xlabel='l', ylabel='b', xlim=(180,-180), ylim=(-25,25))
        ax.axhline(0, color='k', alpha=0.3)
        ax.axvline(0, color='k', alpha=0.3)
        ax= axx[1,1]
        ax.hist(list(t.l), np.linspace(-180, 180,31));
        plt.setp(ax, xlim=(180,-180), xlabel='l')
        ax.grid(True, alpha=0.5)
        return fig
        
    def seed_skyplot(self):
        """seed positions
        """
        t=self.df
        fig, ax = plt.subplots(1,1, figsize=(8,8))
        ax.plot(t.l, t.b, '.')
        plt.setp(ax, xlabel='l', ylabel='b', xlim=(180,-180), ylim=(-90,90))
        ax.axhline(0, color='k', alpha=0.8)
        ax.axvline(0, color='k', alpha=0.8)
        return fig
    
    def pg_cumulative_ts(self):
        """ Cumulative TS for accepted seeds
        """
        fig=self.cumulative_ts(self.udf[self.good].ts, check_localized=False);
        plt.setp(fig.get_axes()[0], xlim=(9,100), ylim=(1,100))
        return fig
    
    def skyplot_good(self):
        """ positions of accepted PGWAVE seeds
        
        """
        udf = self.udf
        good =self.good
        glon = np.array(udf.glon[good],float)
        glon[glon>180]-=360
        fig, ax = plt.subplots(1,1, figsize=(8,6))
        scat = ax.scatter(glon, udf[good].glat, c=udf[good].ts, s=50, edgecolor='none')
        plt.setp(ax, xlabel='l', ylabel='b', xlim=(180,-180), ylim=(-90,90))
        ax.grid(True, alpha=0.5)
        cb =plt.colorbar(scat)
        cb.set_label('TS')
        return fig

    def all_plots(self):
            self.runfigures([self.seed_skyplot,self.skyplot_good, self.pg_cumulative_ts])
            
