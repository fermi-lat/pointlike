"""
Summary of the data

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/data.py,v 1.4 2013/10/11 16:34:16 burnett Exp $

"""

import numpy as np
import pylab as plt
import pandas as pd

from uw.like2 import dataset
from . import analysis_base
from .analysis_base import FloatFormat

class Data(analysis_base.AnalysisBase):
    """Data summary
    <br>look at binned data 
    """
    require='config.txt'
    def setup(self, eminmin=0, **kwargs):
        self.plotfolder = 'data'
        config = eval(open('config.txt').read())
        datadict = config['datadict']
        self.dataset = dataset.DataSet(datadict['dataname'], interval=datadict.get('interval',None),
            irf=config['irf'],)
        t =[]
        for band in self.dataset.dmap:
            if band.emin()<eminmin: continue
            if band.photons()==0: continue
            t.append(dict(emin=round(band.emin()), emax=round(band.emax()), photons=band.photons(), 
                ec=band.event_class(), nside=band.nside(), pixels=band.size()))
        self.df = pd.DataFrame(t, columns='emin emax ec nside pixels photons'.split())
        sm = np.sum(self.df.photons[self.df.emin>=100])
        sp = np.sum(self.df.pixels[self.df.emin>=100])
        self.total_counts ='{:,d}'.format(int(sm))
        self.total_pixels = '{:,d}'.format(int(sp))

            
    def plot_spectrum(self):
        """Energy Spectrum of all data
        Combined front, back. Total above 100 MeV: %(total_counts)s events.
        """
        fig, ax = plt.subplots(figsize=(6,6))
        photons = self.df.photons.values
        combine = photons[0::2]+photons[1::2]
        elow = self.df.emin.values[0::2]
        ehigh = self.df.emax.values[0::2]
        ax.plot(ehigh, combine, ls='steps', lw=2)
        ax.plot([elow[0],elow[0]], (100, combine[0]), '-b', lw=2) # vertical line for first bin
        ax.plot( (elow[0],elow[1]), (combine[0],combine[0]), '-b', lw=2) 
        plt.setp(ax, xscale='log', xlabel='Energy (MeV)', xlim=(10,1e6), yscale='log', ylabel='events/bin')
        ax.grid()
        return fig  

    def pixel_table(self):
        """Table of numbers of pixels
        Total number of pixels above 100 MeV: %(total_pixels)s
        %(pixel_table_html)s
        
        """
        pixel_table_html = self.df.to_html(float_format=FloatFormat(0))
        self.pixel_table_html = pixel_table_html.replace('<td>', '<td class="integer">')
        
    def all_plots(self):
        self.runfigures( [self.plot_spectrum,self.pixel_table])