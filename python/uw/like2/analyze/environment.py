"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/environment.py,v 1.3 2013/07/07 17:41:04 burnett Exp $

"""

import os, pickle, types
import numpy as np
import pylab as plt
import pandas as pd
from scipy import integrate
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
        """ Exposure dependence
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
     
    def psf_plot(self, irfname=None, outfile='psf.csv'):
        r"""PSF size
        
        <br>This is the <em>effective</em> PSF size, the measure of the shape that is relevant 
        for discrimination of a signal in the presence of a uniform background.
        Specifically, if $f(\theta)$ is the normalized PSF, then the size is 
        $1 / \sqrt{\pi \int_0^\infty f(\theta)^2 2\pi \theta \ \mathrm{d}\theta}$. For comparison, 
        the corresponding 68 percent curves are shown as dashed lines.
        <br>PSF filenames: %(psf_files)s
        """
        from uw.like import pypsf, pycaldb
        if irfname is None: irfname=self.config['irf']
        cdm = pycaldb.CALDBManager(irf=irfname)
        psf = pypsf.CALDBPsf(cdm)
        def effective_size(e, ct):
            f2 = lambda delta: psf(e,ct, delta)**2 * 2*np.pi*delta
            return np.degrees(1./np.sqrt(np.pi*integrate.quad(f2, 0, np.inf)[0]))
        self.psf_files=cdm.get_psf()
        egev = np.logspace(-1.+1/8., 2.5+1/8., 3.5*4+1)
        front, back = [[effective_size(e*1e3,ct) for e in egev] for ct in range(2)]
        f68,b68  = [[psf.inverse_integral(e*1e3, ct) for e in egev] for ct in range(2)]
        fig,ax = plt.subplots(figsize=(6,6))
        ax.loglog(egev, front, '-g', lw=2, label='front')
        ax.plot(egev,  back, '-r', lw=2, label='back')
        ax.plot(egev, f68, '--g', lw=1, label='f68')
        ax.plot(egev, b68, '--r', lw=1, label='b68')
        
        plt.setp(ax, xlabel='Energy (GeV)', ylabel='PSF size (deg)',
            xlim=(0.1, 400), ylim=(0.05, 10), title='Effective PSF size')
        ax.legend(prop=dict(size=10)); ax.grid()
        if outfile is None: return fig
        self.psf_df = pd.DataFrame(dict(front=front, back=back), index=egev.round(3))
        self.psf_df.index.name='energy'
        self.psf_df.to_csv(os.path.join(self.plotfolder, outfile))
        print 'wrote file %s' % os.path.join(self.plotfolder, outfile)
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
        self.runfigures([self.psf_plot, self.exposure_plots, self.isotropic_spectrum,])
    