"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/isotropicspectra.py,v 1.1 2013/06/21 20:15:30 burnett Exp $

"""

import numpy as np
import pylab as plt

from . import galacticspectra

class IsotropicSpectra(galacticspectra.GalacticSpectra):
    """Isotropic diffuse refits
      <p>This is a set of plots to check consistency of %(title)s spectra. These result 
        from analysis of a special run that, for each ROI and each energy band, allows this diffuse component to be free.
        This is done three times: using only front, only back, and both.
        <p>There two sets of plots: using both, how consistent is it with the expected unit normalization; 
        and is the front consistent with the back?
        """

    require = 'isofits.zip'

    def setup(self, args=None):
        self.diffuse_setup('iso')
        self.plot_functions[0] += [self.lowdiff_plots]
        self.plot_functions[1] += ['low_energy_difference']

    def lowdiff_hist(self, ax, xlim=(-0.25,0.25), **kwargs):
        plane=abs(self.rois.glat)<10
        v0,v1 = [self.flux['both']['values'][i] for i in (0,1)]
        delta=v1-v0
        kw = dict(histtype='stepfilled'); kw.update(kwargs)
        space = np.linspace(xlim[0], xlim[1], 41)
        ax.hist(delta.clip(*xlim), space, label='all: mean %.2f'%delta.mean(), **kw)
        ax.hist(delta[plane].clip(*xlim), space, color='r', label='plane: mean %.2f'%delta[plane].mean(), **kw)
        ax.legend(); ax.grid();
        ax.axvline(0, color='grey')
        ax.set_title('%s diffuse fits %s'%(self.which,self.skymodel))
        ax.set_xlim(xlim)
        ax.set_xlabel('bin1 - bin0 difference')


    def lowdiff_scat(self, ax, vmin=-0.25, vmax=0.25, **kwargs):
        v0,v1 = [self.flux['both']['values'][i] for i in (0,1)]
        delta=v1-v0
        kw = dict(edgecolor='none');kw.update(kwargs)
        t=ax.scatter(self.rois.glon, self.rois.singlat, c=delta,
            marker='D', s=70,vmin=vmin, vmax=vmax, **kw)
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
        
        This is an indicator of the adequacy of the Limb contribution, since it only affects the lowest energy band.
        <br><b>Left</b>: Histogram of the normalization difference between the two lowest energy bands.
        <br><b>Right</b>: distribution of this over the sky.
        """
        fig,ax=plt.subplots(1,2, figsize=(14,6))
        self.lowdiff_hist( ax[0])
        self.lowdiff_scat( ax[1])
        return fig
        
    def all_plots(self):
        self.runfigures(*self.plot_functions)
