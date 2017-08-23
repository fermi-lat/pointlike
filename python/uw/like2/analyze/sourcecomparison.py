"""
Comparison with the 3FGL catalog

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/sourcecomparison.py,v 1.8 2015/12/03 17:10:03 burnett Exp $

"""

import os
import astropy.io.fits as pyfits
import numpy as np
import pylab as plt
import pandas as pd

from skymaps import SkyDir
from . import sourceinfo

class SourceComparison(sourceinfo.SourceInfo):
    """Comparison with a FITS catalog, 2FGL or beyond
    """

    def setup(self, cat='3FGL-v13r3_v6r9p1_3lacv12p1_v7.fits', #'gll_psc4yearsource_v9_assoc_v6r3p0.fit', #gll_psc_v06.fit', 
            catname='3FGL', **kw):
        super(SourceComparison, self).setup(**kw)
        self.catname=catname
        self.plotfolder='comparison_%s' % catname
        if not os.path.exists('plots/'+self.plotfolder):
            os.mkdir('plots/'+self.plotfolder)
        if cat[0]!='/':
            cat = os.path.expandvars('$FERMI/catalog/'+cat)
        assert os.path.exists(cat), 'Did not find file %s' %cat
        ft = pyfits.open(cat)[1].data
        self.ft=ft # temp
        print 'loaded FITS catalog file %s with %d entries' % (cat, len(ft))
        id_prob = [np.nan]*len(ft)
        try:
            id_prob = ft.ID_Probability_v6r9p1[:,0] ## should find that suffix
        except: 
            print 'warning: id_prob not set' 
        cat_skydirs = map (lambda x,y: SkyDir(float(x),float(y)), 
                           ft.RAJ2000, ft.DEJ2000)
        
        glat = [s.b() for s in cat_skydirs]
        glon = [s.l() for s in cat_skydirs]
        def nickfix(n):
            return n if n[:3]!='PSR' else 'PSR '+n[3:]
        index = map(nickfix, [x.strip() for x in ft.NickName_3FGL]) #Source_Name 
        self.cat = pd.DataFrame(dict(name3=ft.Source_Name_3FGL_1, 
                nickname=map(nickfix, ft.NickName_3FGL), 
                ra=ft.RAJ2000,dec= ft.DEJ2000, 
                ts=ft.Test_Statistic, 
                skydir=cat_skydirs,
                glat=glat, glon=glon, 
                #pivot=ft.Pivot_Energy, flux=ft.Flux_Density, 
                #modelname=ft.SpectrumType, 
                id_prob=id_prob,
                a95=ft.Conf_95_SemiMajor, b95=ft.Conf_95_SemiMinor, 
                                     ang95=ft.Conf_95_PosAng,
                        flags=np.asarray(ft.Flags_3FGL, int),
                ), 
            columns = """name3 nickname ra dec glat glon skydir ts 
                a95 b95 ang95 id_prob flags""".split(), # this to order them
            index=index, )
        self.cat.index.name='name'
        self.cat['pt_flags'] = self.df.flags
        self.cat['pt_ts'] = self.df.ts
        self.cat['pt_ra'] = self.df.ra
        self.cat['pt_dec'] = self.df.dec
        
        def find_close(A,B):
            """ helper function: make a DataFrame with A index containg
            columns of the
            name of the closest entry in B, and its distance
            A, B : DataFrame objects each with a skydir column
            """
            def mindist(a):
                d = map(a.difference, B.skydir.values)
                n = np.argmin(d)
                return (B.index[n], np.degrees(d[n]))
            return pd.DataFrame( map(mindist,  A.skydir.values),
                index=A.index, columns=('otherid','distance'))
                
        if catname=='2FGL' or catname=='3FGL':
            print 'generating closest distance to catalog "%s"' % cat
            closedf= find_close(self.df, self.cat)
            self.df['closest']= closedf['distance']
            self.df['close_name']=closedf.otherid
            closedf.to_csv(os.path.join('plots', self.plotfolder,
                                       'comparison_%s.csv'%catname))
            closest2 = np.degrees(np.array([min(map(sdir.difference, 
                     self.df.skydir.values)) for sdir in cat_skydirs]))
            self.cat['closest']= closest2
            
    def history_check(self):
        """use the history to get the set of 3FGL names of sources
        that were not direct ancestors. Only applies to P7R4 sources
        """
        # the guys with 3FGL ancestors
        in3fgl = [x is not None for x in self.df.cat3fgl]
        t = self.df.cat3fgl[in3fgl].values
        inset = set([x.name3 for x in t])

        #Make the set of all 3FGL names 
        fullset = set(self.cat.name3); 
        # lost: the difference
        lost = fullset.difference(inset);
        missing = [n in lost for n in self.cat.name3]
        missing_cat = self.cat[missing]
        really_missing=  [x.startswith('P7R4') for x in missing_cat.nickname]
        truly_missing = missing_cat[really_missing]
        return truly_missing
    
    def distance_to_cat(self, maxdist=0.5, tscuts=[10,50,500], nbins=26):
        """Associations of sources with 3FGL
        
        """
        fig,ax = plt.subplots( figsize=(5,5))
        hist_kw = dict(bins=np.linspace(0,maxdist,nbins), log=True,
                    histtype='step', lw=2)
        for tscut in tscuts:
            ax.hist(self.df.closest[self.df.ts>tscut].clip(0,maxdist),
                    label='TS>%d'%tscut, **hist_kw)
        ax.grid(True, alpha=0.5)
        ax.legend(prop=dict(size=10))
        plt.setp(ax, xlabel='closest distance to %s source'%self.catname,
                ylim=(0.8,None), )
        return fig
    
    def lost_plots(self, close_cut=0.25, minassocprob=0.8, maxts=250):
        """3FGL sources not present in new list
        Histogram of the 3FGL catalog TS and Galactic latitude for those sources
        more than %(close_cut).2f deg from a skymodel source. 
        The subset of sources with associations (prob>%(minassocprob)s) is shown.
        <br>
        Left: Distribution vs. TS.<br>
        Right: Distribution vs sine of Galactic latitude.
        """ 
        self.minassocprob=minassocprob
        self.close_cut = close_cut
        fig,axx = plt.subplots(1,2, figsize=(8,4))
        self.lost = self.cat.closest>close_cut
        print '%d sources from %s further than %.2f deg: consider lost' % (sum(self.lost) , self.catname, close_cut )
        self.cat.ix[self.lost].to_csv(os.path.join(self.plotfolder,'3fgl_lost.csv'))
        print '\twrite to file "%s"' % os.path.join(self.plotfolder,'3fgl_lost.csv')
        lost_assoc = self.lost & (self.cat.id_prob>0.8)

        def left(ax):
            space = np.linspace(0,maxts,21)
            ax.hist(self.cat.ts[self.lost].clip(0,maxts), space, 
                    label='all (%d)'%sum(self.lost))
            ax.hist(self.cat.ts[lost_assoc].clip(0,maxts), space, 
                    color='orange', label='associated(%d)' %sum(lost_assoc) )
            ax.legend(prop=dict(size=10))
            ax.grid()
            plt.setp(ax, xlabel='TS of %s source' %self.catname)

        def right(ax):
            space = np.linspace(-1,1,51)
            singlat = np.sin(np.radians(self.cat.glat))
            ax.hist(singlat[self.lost], space, 
                    label='all (%d)'%sum(self.lost))
            #lost_assoc = self.lost & (self.cat.id_prob>0.8)
            ax.hist(singlat[lost_assoc], space, color='orange', 
                    label='associated(%d)' %sum(lost_assoc) )
            ax.legend(prop=dict(size=10))
            ax.grid()
            plt.setp(ax, xlabel='sin(glat) of %s source' 
                     %self.catname, xlim=(-1,1))
            return fig
        for f, ax in zip((left,right), axx.flatten()): 
            f(ax)
        return fig
        
    def poorly_localized(self):
        pass
    
    def properties_of_missing(self):
        """Properties of the missing 3FGL sources
        This defines missing from the history, looking at the %(number_lost)d P7R4 sources 
        used to start the 6-year list, which were used in 3FGL, but did not 
        survive to the end
        """
        lost = self.history_check()
        self.number_lost=len(lost)
        maxdist, nbins= 2.6, 53
        fig,axx = plt.subplots(2,1, figsize=(8,12))
        ax=axx[0]
        hist_kw = dict(bins=np.linspace(0,maxdist,nbins), log=True,
                    histtype='step', lw=2)
        ax.hist(lost.closest, **hist_kw)
        plt.setp(ax, xlabel='Closest distance to 6-year source', 
                 ylim=(0.8,None),title='3FGL sources missing from 6-year list')
        ax.grid(True, alpha=0.5)

        ax = axx[1]
        maxdist, nbins= 200, 53
        hist_kw = dict(bins=np.logspace(1,3,nbins), log=False,
                    histtype='step', lw=2)
        ax.hist(lost.ts.clip(0,1000), **hist_kw)
        plt.setp(ax, xlabel='TS', ylim=(0,None),xscale='log', xlim=(25,1000))
        ax.grid(True, alpha=0.5)
        return fig
        
    def all_plots(self):
        """Results of comparison with 3FGL catalog
        """
        self.runfigures([ self.distance_to_cat, self.lost_plots,
                        self.properties_of_missing,])
        
    def lookup_3fgl(self, name3):
        if name3[-1]!=' ' and name3[-1]!='c': name3=name3+' '
        fglnames = list(self.cat.name3)
        try:
            i = fglnames.index(name3)
            nick = self.cat.ix[i].name
            j=  list(self.df.close_name).index(nick)
            return self.df.ix[j]
        except Exception, msg:
            print 'Source %s not found (%s)' % (name3, msg)
            return None