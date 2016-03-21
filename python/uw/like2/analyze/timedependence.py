"""
Analyze time dependence of transient sources 

$Header$

"""
import os, pickle, pyfits, glob
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from skymaps import SkyDir

class TimeDependence(object):

    def __init__(self, month_number):
        self.month_number=month_number
        
        # get the FT1 file for this month
        #ft1_file='/nfs/farm/g/glast/g/catalog/P8_P301/Source/P301_Source_%03d_zmax100.fits'% month_number
        ft1_file='/nfs/farm/g/glast/g/catalog/P8_P302/zmax105/P302_Source_%03d_zmax105.fits'% month_number
        ft1 = pyfits.open(ft1_file)
        events = ft1[1].data
        self.ra = events.RA
        self.dec = events.DEC
        self.time = events.TIME
        self.energy = events.ENERGY
        
        # load the source DataFrame from the analysis
        self.month_name='month%02d' % month_number
        self.sourcedf = pickle.load(open(os.path.expandvars('$FERMI/skymodels/P301_monthly/%s/sources.pickle') %self.month_name))
        
    
    def time_data(self, source_name, max_diff=4, min_e=400):
        adir = self.sourcedf.ix[source_name].skydir
        cut = (np.abs(self.ra-adir.ra())<10)\
            & (np.abs(self.dec-adir.dec())<10)\
            & (self.energy>min_e) 
        sdir = map(SkyDir, np.array(self.ra[cut],float), np.array(self.dec[cut],float))
        diffs = np.asarray(np.degrees( [ adir.difference(b) for b in sdir]))        
        day=(self.time[cut]-self.time.min())/(24*3600.)
        sel = diffs<max_diff
        return diffs[sel], day[sel]
        
    
    def time_plot(self, source_name):
        """
        For the given source, make plots of the angular correlation, and time structure, binned in days
        """
        diffs, day = self.time_data(source_name)
        fig, axx = plt.subplots(1,2, figsize=(12,4))

        ax = axx[0]
        dsq = diffs**2
        ax.hist(dsq, np.linspace(0,4,17), histtype='step');
        n= len(dsq[ (dsq>1) & (dsq<4)])/12. # number of bins
        ax.axhline(n, ls='--', color='gray');
   
        ax.grid(True, alpha=0.5)
        plt.setp(ax, xlabel='distance**2 [deg**2]')

        ax=axx[1]
        hist_kw = dict(bins=np.linspace(0, 30,31), histtype='step', log=True)
        ax.hist(day, **hist_kw )
        ax.hist(day[diffs**2<0.5], color='red', **hist_kw);
        ax.grid(True, alpha=0.5)
        plt.setp(ax, ylim=(0.8,None),xlabel='day of month' )
        fig.suptitle('%s for month %d' % (source_name, self.month_number) )
        return fig
        
    def ts_sources(self):
        sdf = self.sourcedf
        ts_sources = np.asarray([s.startswith('TS') for s in sdf.index], bool)
        cands = sdf.ix[ts_sources & (sdf.ts>10) & (np.abs(sdf.glat)>10)\
                        & (sdf.locqual<8)]['ts ra dec a locqual'.split()]
        cands.sort_index()
        return cands

    def create_dict(self):
        ret = dict()
        daybins=np.linspace(0, 30,31)
        for sname in self.ts_sources().index:
            diffs, day = self.time_data(sname)
            h2c=np.histogram(day[diffs**2<0.5],daybins)[0]
            ret[sname] = dict(h1= np.histogram(diffs**2, np.linspace(0,4,17))[0], 
                h2=np.histogram(day, daybins)[0],
                h2c=h2c,
                maxday=h2c.max(),
                numday=h2c.sum(),
                )
        return ret