"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/environment.py,v 1.5 2013/07/15 13:40:22 burnett Exp $

"""

import os, pickle, types
import numpy as np
import pylab as plt
import pandas as pd
from scipy import integrate
#from skymaps import SkyDir #?

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
     
    def get_psf(self, irfname=None, ):
        from uw.like import pypsf, pycaldb
        if irfname is None: irfname=self.config['irf']
        cdm = pycaldb.CALDBManager(irf=irfname)
        self.psf_files=cdm.get_psf()
        return pypsf.CALDBPsf(cdm)

    def psf_plot(self, irfname=None, outfile='psf.csv', title=''):
        r"""PSF size
        
        <br>Plots of two science-driven measures of the size of the PSF, compared with the standard 68 percent containment.
        If  $f(\delta)$ is the normalized PSF as a function of deviation $\delta$ for a given energy and conversion type, 
        averaged over incidence angle, then the two cases are:
        <ol>
        <li>Discrimination of a point source 
        signal in the presence of a uniform background; the <em>average</em> PSF is the inverse of the solid angle, with radius
        $$1 / \sqrt{\pi \int_0^\infty f(\delta)^2 2\pi \delta \ \mathrm{d}\delta}$$
        <li> Measurement of the position of a source. In this case the measure is the RMS of the deviation, or
        $$\sqrt{ \int_0^\infty f(\delta) 2\pi \delta^3 \ \mathrm{d}\delta}$$
        </li>
        
        <br>PSF filenames: %(psf_files)s
        """
        psf = self.get_psf(irfname)
 
        def bkg_size(e, ct):
            f2 = lambda delta: psf(e,ct, delta)**2 * 2*np.pi*delta
            return np.degrees(1./np.sqrt(np.pi*integrate.quad(f2, 0, np.inf)[0]))
        
        def loc_size(e,ct):
            f3 = lambda theta: psf(e,ct, theta) * theta**3 * 2.*np.pi
            return np.degrees(np.sqrt(integrate.quad(f3, 0, np.pi/6) [0]))

        
        egev = np.logspace(-1.+1/8., 2.5+1/8., 3.5*4+1)
        #egev = np.logspace(-1, 2, 100)
        front, back = [[bkg_size(e*1e3,ct) for e in egev] for ct in range(2)]
        floc, bloc = [[loc_size(e*1e3,ct) for e in egev] for ct in range(2)]
        
        f68,b68  = [[psf.inverse_integral(e*1e3, ct) for e in egev] for ct in range(2)]
        fig,ax = plt.subplots(figsize=(6,6))
        ax.loglog(egev, front, '-g', lw=2, label='front bkg')
        ax.plot(egev,  back, '-r', lw=2, label='back bkg')
        ax.plot(egev, floc, '--g', lw=2, label='front loc')
        ax.plot(egev,  bloc, '--r', lw=2, label='back loc')
        ax.plot(egev, f68, ':g', lw=2, label='front 68')
        ax.plot(egev, b68, ':r', lw=2, label='back 68')
        
        plt.setp(ax, xlabel='Energy (GeV)', ylabel='PSF size (deg)',
            xlim=(0.1, 100), ylim=(0.05, 10), title=title)
        ax.legend(prop=dict(size=10)); ax.grid()
        if outfile is None: return fig
        self.psf_df = pd.DataFrame(dict(front=front, floc=floc, back=back, bloc=bloc,),  index=egev.round(3))
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
    