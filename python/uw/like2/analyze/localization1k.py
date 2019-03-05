"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/localization1k.py,v 1.2 2017/08/23 16:23:43 zimmer Exp $

"""

import numpy as np
import pylab as plt
import pandas as pd
import astropy.io.fits as pyfits

from . import localization

class Localization1K(localization.Localization):
    """ load and analyze a special localization-only run"""
    require='localization.zip'
    
    def setup(self, zipname='localization', ecut=1000,  **kw):
        super(Localization1K, self).setup(**kw)
        try:
            f, pk = self.load_pickles(zipname)
        except:
            print 'failed to load %s, which should have zipped pickles of sources after localization attempt'%zipname
            raise
        d = dict( (x.name, x.ellipse if hasattr(x,'ellipse') else [np.nan]*7) for x in pk)
        self.ebox1K = pd.DataFrame( d ).T  
        self.ebox1K.columns = self.ebox.columns
        self.plotfolder='sources' #needed by superclass
        
    def ellipse_ratio(self):
        """ Ratio of error ellipse for E>1GeV
        Compare the error ellipse major axis for the full fit, with a fit using energies above 1 GeV 
        <br>Left: scatter plot of the ratio vs. the value, for TS>25
        <br>Right: Histogram of the ratio, showing subset with TS>100
        """
        fig, axx = plt.subplots(1,2, figsize=(12,5))
        aratio = self.ebox1K.a/self.ebox.a
        ratio_label = 'Ratio for E>1 GeV'
        cut = self.df.ts>25
        def plot1(ax):
            ax.plot( 2.5*60*self.ebox.a[cut],aratio[cut], '.')
            plt.setp(ax, xscale='log', xlim=(0.5,10), xlabel='error ellipse major axis (arc min)', 
                         yscale='linear', ylim=(0.95,1.1), ylabel=ratio_label)
            ax.axhline(1.0, color='k')
            ax.grid()
        def plot2(ax):
            xlim = (0.90, 1.25); bins=np.linspace(xlim[0],xlim[1],36)
            tscut=100
            cut2 = self.df.ts>tscut
            ax.hist(aratio[cut].clip(*xlim),bins , label='%d sources'%sum(cut.values))
            ax.hist(aratio[cut2].clip(*xlim), bins, label='TS>%d'%tscut)
            ax.legend(loc='upper right', prop=dict(size=10))
            plt.setp(ax, xlim=xlim, xlabel =ratio_label)
            ax.grid(); ax.axvline(1.0, color='k')
            
        plot1(axx[0])
        plot2(axx[1])
        
    def all_plots(self, **kw):
        self.runfigures([self.ellipse_ratio,])