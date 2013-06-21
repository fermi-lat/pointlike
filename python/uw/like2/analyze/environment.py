"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/environment.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import os, pickle, types
import numpy as np
import pylab as plt
import pandas as pd
from skymaps import SkyDir #?

from . import roi_info

class Environment(roi_info.ROIinfo):
    """ Plots associated with the environment"""

    
    def setup(self, **kw):
        super(Environment, self).setup(**kw)
        self.plotfolder='environment'
        # use the fact that the isotopic diffuse compoenent is isotropic, so that
        # the ratio of the computed counts, to the fit normalization, is proportional
        # to the exposure.
        iso = self.model_counts('isotrop')
        models = self.diffuse_models('isotrop')
        norms = np.array([m.getp(0) if m is not None else np.nan for m in models])
        self.relative_exp = iso/norms/(iso/norms).mean()
        self.config = eval(open('config.txt').read())
        
    def exposure_plots(self, hsize=(1.0,1.0,2.0,1.0, 2.0, 0.7),):
        """ exposure dependence
        Examine the relative exposure, per ROI. Express in terms of the mean. Note that
        ROIs are distributed uniformly over the sky.
        <p>Use the fact that the isotopic diffuse compoenent is isotropic, so that
        the ratio of the computed counts, to the fit normalization, is proportional
        to the exposure. This involves all energies, but is weighted according to the isotropic diffuse.

        <br>Left: histogram, center: scatter plot vs. Declination; right: map on sky, in Galactic coordinates.
        
        """
        #fig, axx = plt.subplots(1,3, figsize=(15,4))
        fig, axx = self.subplot_array(hsize, figsize=(12,4))
        relative_exp = self.relative_exp
        label = 'exposure relative to mean'
        lim = (0.7, 1.6)
        def left(ax):
            ax.hist(relative_exp, np.linspace(*lim, num=25))
            plt.setp(ax, xlim=lim)# xlabel=label)
            ax.axvline(1.0, color='k')
            ax.grid()

        def center(ax):
            ax.plot(self.df.dec, relative_exp, '.')
            ax.grid()
            plt.setp(ax, xlim=(-90,90), xlabel='Dec (deg)',ylabel=label, ylim=lim)
            ax.set_xticks(range(-90,91,30))
            ax.axhline(1, color='k')
        def right(ax):
            self.skyplot(relative_exp, ax=ax, s=40)
        
        for f,ax in zip((left, center, right), axx.flatten()): f(ax)
        return fig
     
    def psf_plot(self):
        """PSF plot
        PSF files: %(psf_files)s
        <br>This is an effecive PSF size, derived from the value of the normalized function at the peak.
        """
        from uw.like import pypsf, pycaldb
        irfname=self.config['irf']
        cdm = pycaldb.CALDBManager(irf=irfname)
        psf = pypsf.CALDBPsf(cdm)
        self.psf_files=cdm.get_psf()
        #egev = np.logspace(-1,2.5,121)
        egev = np.logspace(-1.+1/8., 2.5+1/8., 3.5*4+1)
        front, back = [[np.degrees(1./np.sqrt(psf(e*1e3,ct,0)))[0] for e in egev] for ct in range(2)]
        fig,ax = plt.subplots(figsize=(5,5))
        ax.loglog(egev, front, '-g', lw=2, label='front')
        ax.plot(egev,  back, '-r', lw=2, label='back')
        plt.setp(ax, xlabel='Energy (GeV)', ylabel='PSF size (deg)',
            xlim=(0.1, 400), ylim=(0.05, 10), title='Effective PSF size')
        ax.legend(prop=dict(size=10)); ax.grid()
        self.psf_df = pd.DataFrame(dict(front=front, back=back), index=egev.round(3))
        self.psf_df.index.name='energy'
        self.psf_df.to_csv(os.path.join(self.plotfolder, 'psf.csv'))
        print 'wrote file %s' % os.path.join(self.plotfolder, 'psf.csv')
        return fig
        
    def isotropic_spectrum(self, other=None):
        """ Isotropic Spectrum from template
        
        The spectrum used to define the isotropic diffuse component.
        <br>Files for front/back: %(idfiles)s
        """
        # look up filenames used to define the isotorpic spectrum: either new or old diffuse spec; list or dict
        diffuse=self.config['diffuse']
        isokey = 'isotrop' if type(diffuse)==types.DictType else 1
        self.idfiles = [os.path.join(os.environ['FERMI'],'diffuse',diffuse[isokey][i]) for i in (0,1)]
        nf,nb = map(np.loadtxt, self.idfiles)
        energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
        fig, axs = plt.subplots(1,2, figsize=(10,5), dpi=50)
        def right(ax):
            ax.plot(energies, front/back, '-o');
            ax.axhline(1.0, color='k')
            plt.setp(ax, xscale='log', xlabel='Energy');ax.grid(True);
            ax.set_title('Isotropic flux front/back ratio', fontsize='small');
        def left(ax):
            ax.plot(energies, front*energies**2, '-g', label='front')
            ax.plot(energies, back*energies**2, '-r', label='back')
            plt.setp(ax, xlabel='Energy', ylabel='flux*e**2', xscale='log')
            ax.set_title('isotropic diffuse spectra', fontsize='small')
            ax.grid(True); ax.legend()
        for f,a in zip((left,right), axs.flatten()): f(a)
        self.iso_df = pd.DataFrame(dict(front=front*energies**2, back=back*energies**2), index=(energies/1e3).round(3))
        self.iso_df.index.name='energy'
        self.iso_df.to_csv(os.path.join(self.plotfolder, 'isotropic.csv'))
        print 'wrote file %s' % os.path.join(self.plotfolder, 'isotropic.csv')

        return fig
        
    def all_plots(self, **kw):
        self.runfigures([self.exposure_plots, self.psf_plot, self.isotropic_spectrum,])
    