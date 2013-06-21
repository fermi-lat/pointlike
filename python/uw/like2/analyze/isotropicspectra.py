"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/isotropicspectra.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import pickle
import numpy as np
import pylab as plt

from . import galacticspectra

class IsotropicSpectra(galacticspectra.GalacticSpectra):
    require = 'isofits.zip'

    def setup(self):
        self.diffuse_setup('iso')
        self.plot_functions[0] += [self.lowdiff_plots]
        self.plot_functions[1] += ['low_energy_difference']

    def lowdiff_hist(self, ax, **kwargs):
        plane=abs(self.rois.glat)<10
        v0,v1 = [self.flux['both']['values'][i] for i in (0,1)]
        delta=v1-v0
        kw = dict(histtype='stepfilled'); kw.update(kwargs)
        space = np.linspace(-0.5, 0.5, 41)
        ax.hist(delta, space, label='all: mean %.2f'%delta.mean(), **kw)
        ax.hist(delta[plane], space, color='r', label='plane: mean %.2f'%delta[plane].mean(), **kw)
        ax.legend(); ax.grid();
        ax.axvline(0, color='grey')
        ax.set_title('%s diffuse fits %s'%(self.which,self.skymodel))
        ax.set_xlabel('bin1 - bin0 difference')


    def lowdiff_scat(self, ax, vmin=-0.2, vmax=0.2, **kwargs):
        v0,v1 = [self.flux['both']['values'][i] for i in (0,1)]
        delta=v1-v0
        kw = dict(edgecolor='none');kw.update(kwargs)
        t=ax.scatter(self.rois.glon, self.rois.singlat, c=delta, s=50,vmin=vmin, vmax=vmax, **kw)
        plt.setp(ax, xlabel='glon', ylabel='sin(glat)', xlim=(180,-180), ylim=(-1,1))
        ax.axhline(0, color='k')
        ax.axvline(0, color='k')
        ax.axhline(np.sin(np.radians(10.)), lw=2, color='grey')   
        ax.axhline(np.sin(np.radians(-10.)), lw=2, color='grey')
        # draw poles outline
        try:
            poles = pickle.load(open('../../polar_circle.pickle'))
        except:
            print 'could not find the polar_circle file'
            return
        for i in range(2):
            ax.plot(poles[0,:,0]-360, np.sin(np.radians(poles[0,:,1])), '-',lw=2, color='grey')
            ax.plot(poles[1,:,0], np.sin(np.radians(poles[1,:,1])), '-', lw=2,color='grey')

    def lowdiff_plots(self):
        """ Isotropic bin0-bin1 differences 
        
        This is an indicator of the adequacy of the Limb contribution, sinc it only affects the lowest energy band.
        <br>Left: Histogram of the normalization difference between the two lowest energy bands.
        <br>Right: distribution of this over the plane, with the 45-degree Dec shown.
        """
        fig,ax=plt.subplots(1,2, figsize=(14,6))
        self.lowdiff_hist( ax[0])
        self.lowdiff_scat( ax[1])
        return fig