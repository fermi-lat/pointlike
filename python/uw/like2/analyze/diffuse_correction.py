"""
Diffuse Correction

$Header$
"""

import os, pickle
import numpy as np
import pylab as plt
import pandas as pd
from skymaps import Band, SkyDir

from . import analysis_base
from . analysis_base import FloatFormat
from .. import (configuration, diffuse,)
from ..pub import healpix_map

class DiffuseCorrection(analysis_base.AnalysisBase):
    """Study the galactic diffuse correction
    This assumes that the process fitdiffuse has been run, which adjusts the input diffuse map cube
    by adding a gaussian centered on the ROI. 
    """
    
    def setup(self, **kwargs):
        """
        """
        self.plotfolder = 'diffusecorr'
        self.all_corr=np.zeros((1728,8))
        self.chisq = np.zeros(1728)
        self.chisq_init = np.zeros(1728)
        for index in range(1728):
            filename = 'diffuse/history_{:04d}.pickle'.format(index)
            try:
                t=pd.DataFrame(pickle.load(open(filename))).T;
            except Exception, msg:
                raise Exception('Failed to open file {}:'.format(filename, msg))
            c =np.array([t.ix[i]['corr'] for i in range(len(t))])
            self.all_corr[index,:] = np.prod(c,axis=0).round(3)
            self.chisq[index]=list(t.chisq)[-1].round(1)
            self.chisq_init[index]=list(t.chisq)[0].round(1)
        sdirs = map(Band(12).dir, range(1728))
        self.glon = np.array([s.l() for s in sdirs])
        self.glon[self.glon>180]-=360
        self.singlat = np.sin(np.radians([s.b() for s in sdirs]))
        # get the energy list from the current diffuse cube
        config = configuration.Configuration('.', quiet=True, postpone=True)
        t = diffuse.HealpixCube(config.diffuse['ring']['filename'])
        t.load(); self.energies = t.energies
        self.iso=diffuse.diffuse_factory(config.diffuse['isotrop'])
        print 'isotropic:\n {}'.format(self.iso)
        self.gal=diffuse.diffuse_factory(config.diffuse['ring'])
        print 'galactic:\n{}'.format(self.gal)


    def plot_chisq(self):
        """Chi squared hists
        For the count maps, before and after the new set of adjustments per ROI
        """
        fig, ax = plt.subplots(figsize=(5,5))
        ax.hist(self.chisq.clip(0,50), np.linspace(0,50, 26), lw=2,color='blue', histtype='step', label='after fit');
        ax.hist(self.chisq_init.clip(0,50), np.linspace(0,50, 26), lw=2,color='orange', histtype='step',label='before fit');
        fig.set_facecolor('white');
        ax.set_xlabel('chi squared per ROI')
        ax.grid(); ax.legend(loc='upper left')
        return fig
        
    def corr_maps(self, vmin=0.9, vmax=1.1):
        """\
        Maps of the correction factors
        The values of the corrections for each ROI
        """
        fig, axx = plt.subplots(2,4, figsize=(12,8), sharex=True, sharey=True)
        plt.subplots_adjust(right=0.9, hspace=0.15, wspace=0.1)
        for ib in range(8): #,energy in enumerate(self.energy[:8]):
            ax = axx.flatten()[ib]
            scat=self.basic_skyplot(ax, self.glon, self.singlat, self.all_corr[:,ib].clip(vmin,vmax),
                title='band %d '%ib,
                vmin=vmin,vmax=vmax, s=15, edgecolor='none', colorbar=False, labels=False)
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('normalized residual')
        fig.text(0.5,0.05, 'galactic longitude', ha='center')
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center')
        return fig
        
    def allsky_map(self, nside=128, clip=(-0.2, 0.2), saveto='smoothed_all_corr.pickle'):
        """\
        Allsky map of the combined, smoothed corrections
        Plot is for 133 MeV
        <p>Statistics: %(corr_map_stats)s
        """
        allg = np.zeros((12*nside**2,8))
        for index in range(1728):
            delta = self.all_corr[index,:]-1
            gfile = os.path.expandvars('$FERMI/misc/gausspat2/hp12_{:04d}.pickle'.format(index))
            gausspat = pickle.load(open(gfile))
            allg += np.outer(gausspat, delta) 
        self.allg=allg.clip(*clip)
        if saveto is not None:
            pickle.dump(self.allg,open(saveto,'w'))
        self.corr_map_stats ='Saved file {} with clips {} <p>'.format(saveto, clip)\
                +pd.DataFrame([self.allg.mean(axis=0),self.allg.std(axis=0)],                     
                    index='mean std'.split()).to_html(float_format=FloatFormat(4)) 
        healpix_map.HParray('band 0', allg[:,0]).plot(vmin=-0.1, vmax=0.1);
        return plt.gcf()
    
    def all_plots(self):
        self.runfigures([
            self.plot_chisq,
            self.corr_maps,
            self.allsky_map,
             ])
