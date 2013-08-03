"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/data.py,v 1.1 2013/06/21 20:15:30 burnett Exp $

"""

import numpy as np
import pylab as plt
import pandas as pd

from uw.like2 import dataset
from . import analysis_base

class Data(analysis_base.AnalysisBase):
    """Data summary
    """
    require='config.txt'
    """ look at binned data """
    def setup(self):
        self.plotfolder = 'data'
        config = eval(open('config.txt').read())
        datadict = config['datadict']
        self.dataset = dataset.DataSet(datadict['dataname'], interval=datadict.get('interval',None),
            irf=config['irf'],)
        t =[]
        for band in self.dataset.dmap:
            t.append(dict(emin=band.emin(), emax=band.emax(), photons=band.photons(), 
                ec=band.event_class(), nside=band.nside(), pixels=band.size()))
        self.df = pd.DataFrame(t)
        sm= np.sum(self.df.photons[self.df.emin>=100])
        self.total_counts ='{:,d}'.format(int(sm))

            
    def plot_spectrum(self):
        """Spectrum of all data
        Total above 100 MeV: %(total100)s events.
        """
        fig, ax = plt.subplots(figsize=(4,4))
        photons = self.df.photons.values
        combine = photons[0::2]+photons[1::2]
        self.total100 = np.sum(combine[3:])
        elow = self.df.emin.values[0::2]
        ehigh = self.df.emax.values[0::2]
        ax.plot(ehigh, combine, ls='steps', lw=2)
        ax.plot([elow[1],elow[1]], (100, combine[1]), '-b', lw=2) # vertical line for first bin
        plt.setp(ax, xscale='log', xlabel='Energy (MeV)', xlim=(10,1e6), yscale='log', ylabel='events/bin')
        ax.grid()
        return fig  

    def all_plots(self):
        self.runfigures( [self.plot_spectrum,])