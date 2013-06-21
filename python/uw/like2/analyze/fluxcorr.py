"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/fluxcorr.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import numpy as np
import pylab as plt

from . import sourceinfo

class FluxCorr(sourceinfo.SourceInfo):

    require='fluxcorr.zip'
    def setup(self, **kwargs):
        super(FluxCorr, self).setup(**kwargs)
        self.plotfolder='fluxcorr'
        self.source_name=kwargs.pop('source_name', 'fluxcorr')
        self.title='Source-galactic diffuse flux dependence'
        self.diffuse_name='Galactic'
        self.delta = kwargs.get('delta', 0.01) # default delta is 1%

        
        self._readdata()
        
    def _readdata(self):
        # read in the flux correlation data, in DataFrame, combine to DataFrame
        fs, ps = self.load_pickles(self.source_name)
        print 'Combining with Source info...'
        #first combine the roi DataFrames
        ndf = None 
        for x in ps:
            if len(x)>0:
                ndf = ndf.append(x) if ndf is not None else x
                
        self.df = self.df.join(ndf)
        self.emins=(100, 316, 1000) # assume 
        for emin in self.emins:
            self.df['flux_dependence_ratio_%d'%emin]= \
                (10**(ndf['par_p%d'%emin] - ndf['par_m%d'%emin])-1)/(2.*self.delta)
        
    def flux_sensitivity(self, axx=None, emin=100, **kwargs):
        """ %(diffuse_name)s diffuse flux sensitivity, emin=%(emin)s
        
        Let fs be the flux sensitivity, defined as the ratio of measured flux change to change in Galactic diffuse
        left: histogram of fs<br>
        center: scatter plot, flux sensitivity vs. TS<br>
        rignt: Skymap for fs<100.
        """
        self.emin = emin
        colorbar = kwargs.pop('colorbar', True)
        if axx is None:
            fig, axx = plt.subplots(1,3, figsize=(10,5))
        flux_ratio = self.df['flux_dependence_ratio_%d'%emin]
        def plot1(ax):
            bins = np.linspace(-30,5,36)
            ax.hist(flux_ratio.clip(bins[0],bins[-1]), bins, label='%d sources'%flux_ratio.count())
            ax.hist(flux_ratio[self.df.ts<100].clip(bins[0],bins[-1]), bins, label='TS<100')
            ax.hist(flux_ratio[self.df.ts<25].clip(bins[0],bins[-1]), bins, label='TS<25')
            plt.setp(ax, xlabel='flux sensitivity')
            ax.grid()
            ax.legend(loc='upper left', prop=dict(size=10))
        def plot2(ax):
            ax.plot( self.df['ts'], flux_ratio, '.')
            plt.setp(ax, xscale='log', xlim=(10,1e5), ylim=(-30,5) , xlabel='TS')
            ax.grid()
        def plot3(ax):
            tscut = (self.df.ts<100) & (self.df.ts>10)
            self.skyplot(-flux_ratio[tscut], ax=ax ,s=15, vmax=20, vmin=0, 
                colorbar=colorbar,cbtext='abs(flux sensitivity)')
        for ax, plot in zip(axx.flatten(), (plot1,plot2, plot3)):
            plot(ax)
   
    def flux_sensitivity_all(self, colorbar=False):
        """ %(diffuse_name)s diffuse flux sensitivity
        Rows are, from the bottom, for emin=%(emins)s MeV.<br>
        Let fs be the flux sensitivity, defined as the ratio of the measured flux change 
        to the change in %(diffuse_name)s diffuse flux.
        Columns are:
        left: histogram of fs;
        center: scatter plot, flux sensitivity vs. TS; 
        right: Skymap for TS<100, showing  where the highest sensitivity is located.
        """

        #fig, axx = plt.subplots(3,3, figsize=(9,12))
        fig, axx = self.subplot_array(hsize=(1.0, 0.75, 1.2, 0.5, 2.0,0.5 ),vsize=[1.0,0.75]*3, figsize=(10,10))
        for i, emin in enumerate(self.emins):
            self.flux_sensitivity(axx[i,:], emin, colorbar=colorbar)
        return fig
    
    def ratio_vs_stat(self):
        """Compare statistical error with systematic dependence
        The red diagonal line corresponds to the relative statistical error being equal to the flux change for a 1 percent change in the 
        %(diffuse_name)s flux. Sources below this line are at risk.
        The color represent the absolute value of the Galactic latitude.
        """
        fig, axx = plt.subplots(3,1, squeeze=False, figsize=(6,8))
        plt.subplots_adjust(left=0.2,right=0.85)
        def plot1(ax, emin=100):
            fdr = -self.df['flux_dependence_ratio_%d' % emin]
            relflux_unc = self.df['unc_z%d'%emin]*100.
            scat=ax.scatter( fdr, relflux_unc, c = np.abs(np.array(self.df.glat,float)),edgecolor='none', s=10)
            plt.setp(ax, xlabel='%s diffuse dependence ratio'%self.diffuse_name, xscale='log', xlim=(0.1,40),
                 yscale='log', ylim=(1,100));
            ax.plot([1.0, 100], (1,100.0),'r-')
            ax.text(8, 1.5, 'emin=%d'%emin)
            ax.text(25,30, '1%', color='r')
            ax.grid(True)
            return scat
        scats = map(plot1, axx.flatten(), self.emins)
        cbax = fig.add_axes((0.92, 0.3, 0.025, 0.4) )
        cb = plt.colorbar(scats[0], cbax, orientation='vertical')
        fig.text(0.1, 0.5, 'relative flux uncertainty (%)', va='center', rotation='vertical')
        cb.set_label('abs(b)')
        return fig
    
    def ts_ratio(self):
        """ Ratio of TS for fits at 100 and 1000 Mev
        Compare the ratio of the TS values for fits at 100 MeV, with that at 1 GeV, vs the 100 MeV TS
        """
        fig, ax = plt.subplots(1,1, figsize=(10,5))
        ax.plot(self.df.ts_z100, self.df.ts_z1000/self.df.ts_z100, '.')
        ax.axhline(1.0, color='k')
        plt.setp(ax, xscale='log',  xlim=(10,1e5), xlabel='TS', 
                    yscale='log', ylim=(0.1, 2), ylabel='TS_1000/TS_100')
        ax.grid()
        return fig
        
    def ts_dependence(self):
        """ TS dependence on diffuse
        """
        ts100 = self.df.ts_z100
        tsdiff = ts100 - self.df.ts_p100
        fig, axx = plt.subplots(2,1, figsize=(10,10))
        def plotit(ax, tsdiff):
            scat =ax.scatter(ts100, tsdiff/np.sqrt(ts100), s=20, edgecolor='none', c = np.abs(np.asarray(self.df.glat,float)))
            plt.setp(ax, xscale='log', xlabel='delta TS/sqrt(TS)', ylim=(-0.001, 2), xlim=(10,10000))
            ax.axhline(0, color='k')
            cbar=fig.colorbar(scat, ax=ax); cbar.set_label('|b|')
            ax.grid()
        plotit(axx[0], tsdiff)
        plotit(axx[1], self.df.ts_m100-ts100)
        return fig
    
    def all_plots(self):
        self.runfigures([self.flux_sensitivity_all, self.ratio_vs_stat, self.ts_ratio,])