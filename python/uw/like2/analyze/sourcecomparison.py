"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/sourcecomparison.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import os, pyfits
import numpy as np
import pylab as plt
import pandas as pd

from skymaps import SkyDir
from . import sourceinfo

class SourceComparison(sourceinfo.SourceInfo):
    """Comparison with the 2FGL catalog
    """

    def setup(self, cat='gll_pscP72Y_v5r2_flags_assoc_v5r11p3.fit', #gll_psc_v06.fit', 
            catname='2FGL', **kw):
        super(SourceComparison, self).setup(**kw)
        self.catname=catname
        self.plotfolder='comparison_%s' % catname
        if cat[0]!='/':
            cat = os.path.expandvars('$FERMI/catalog/'+cat)
        assert os.path.exists(cat), 'Did not find file %s' %cat
        ft = pyfits.open(cat)[1].data
        print 'loaded FITS catalog file %s with %d entries' % (cat, len(ft))
        name = ft.Source_Name
        ra = ft.RAJ2000
        dec= ft.DEJ2000
        id_prob = [np.nan]*len(ft)
        try:
            id_prob = ft.ID_Probability[:,0]
        except: pass 
        cat_skydirs = map (lambda x,y: SkyDir(float(x),float(y)), ra,dec)
        glat = [s.b() for s in cat_skydirs]
        glon = [s.l() for s in cat_skydirs]
        # insert space to agree with my PSR name
        index = ft.NickName # note that need to squeze out blanks for comparison
        self.cat = pd.DataFrame(dict(ra=ft.RAJ2000,dec= ft.DEJ2000, ts=ft.Test_Statistic, 
                glat=glat, glon=glon, pivot=ft.Pivot_Energy, flux=ft.Flux_Density, modelname=ft.SpectrumType, id_prob=id_prob), 
            columns = 'ra dec glat glon ts pivot flux modelname id_prob'.split(), # this to order them
            index=index, ) #Source_Name )
        self.cat.index.name='name'
        
        if catname=='2FGL':
            print 'generating closest distance to catalog "%s"' % cat
            closest = np.degrees(np.array([min(map(sdir.difference, cat_skydirs))for sdir in self.df.skydir.values]))
            self.df['closest'] = closest
            closest2 = np.degrees(np.array([min(map(sdir.difference, self.df.skydir.values)) for sdir in cat_skydirs]))
            self.cat['closest']= closest2
        
            
    def distance_to_cat(self, maxdist=0.5, tscuts=[10,50,500], nbins=26):
        """Associations of sources with 2FGL
        
        """
        fig,ax = plt.subplots( figsize=(4,4))
        for tscut in tscuts:
            ax.hist(self.df.closest[self.df.ts>tscut].clip(0,maxdist), np.linspace(0,maxdist,nbins), log=True,
             label='TS>%d'%tscut)
        ax.grid()
        ax.legend(prop=dict(size=10))
        plt.setp(ax, xlabel='closest distance to %s source'%self.catname)
        return fig
    
    def lost_plots(self, close_cut=0.25, minassocprob=0.8, maxts=250):
        """2FGL sources not present in new list
        Histogram of the 2FGL catalog TS and Galactic latitude for those sources more than %(close_cut).2f deg from a skymodel source. 
        The subset of sources with associations (prob>%(minassocprob)s) is shown. <br>
        Left: Distribution vs. TS.<br>
        Right: Distribution vs sine of Galactic latitude.
        """ 
        self.minassocprob=minassocprob
        self.close_cut = close_cut
        fig,axx = plt.subplots(1,2, figsize=(8,4))
        self.lost = self.cat.closest>close_cut
        print '%d sources from %s further than %.2f deg: consider lost' % (sum(self.lost) , self.catname, close_cut )
        self.cat.ix[self.lost].to_csv('2fgl_lost.csv')
        print '\twrite to file "2fgl_lost.csv"'
        lost_assoc = self.lost * self.cat.id_prob>0.8

        def left(ax):
            space = np.linspace(0,maxts,21)
            ax.hist(self.cat.ts[self.lost].clip(0,maxts), space, label='all (%d)'%sum(self.lost))
            ax.hist(self.cat.ts[lost_assoc].clip(0,maxts), space, label='associated(%d)' %sum(lost_assoc) )
            ax.legend(prop=dict(size=10))
            ax.grid()
            plt.setp(ax, xlabel='TS of %s source' %self.catname)

        def right(ax):
            space = np.linspace(-1,1,51)
            singlat = np.sin(np.radians(self.cat.glat))
            ax.hist(singlat[self.lost], space, label='all (%d)'%sum(self.lost))
            lost_assoc = self.lost * self.cat.id_prob>0.8
            ax.hist(singlat[lost_assoc], space, label='associated(%d)' %sum(lost_assoc) )
            ax.legend(prop=dict(size=10))
            ax.grid()
            plt.setp(ax, xlabel='sin(glat) of %s source' %self.catname, xlim=(-1,1))
            return fig
        for f, ax in zip((left,right), axx.flatten()): 
            f(ax)
        return fig
        
    def all_plots(self):
        """Results of comparison with 2FGL catalog
        """
        self.runfigures([ self.distance_to_cat, self.lost_plots])