"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/environment.py,v 1.7 2013/07/18 12:24:32 burnett Exp $

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
        # use the fact that the isotopic diffuse compoenent is isotropic, so that
        # the ratio of the computed counts, to the fit normalization, is proportional
        # to the exposure.
        iso = self.model_counts('isotrop')
        models = self.diffuse_models('isotrop')
        norms = np.array([m.getp(0) if m is not None else np.nan for m in models])
        relative_exp = iso/norms/(iso/norms).mean()
        #fig, axx = plt.subplots(1,3, figsize=(15,4))
        fig, axx = self.subplot_array(hsize, figsize=(12,4))
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
 <li> Measurement of the position of a source. This case involves the expected likelihood. The expected value for the value of the log likelihood, as a function of $\delta$
$$ w(\delta) = \int_0^\infty  \theta\ \mathrm{d} \theta\ P(\theta) \int_0^{2\pi} \mathrm{d}\phi \ln P(\sqrt{\delta^2 -2\delta\  \theta \cos(\phi)+\theta^2})
$$
The resolution is the curvature, or second derivative of this function evaluated at $\delta=0$. A simple way to evaluate this is the expression
$$ \delta \sqrt{\frac{1}{2(w(0)-w(\delta))}}
$$
for $\delta \rightarrow 0$. Note that this is represents the large $N$ limit, for which the resolution for $N$ photons is 
this value divided by $\sqrt{N}$. For small $N$, less than 10, the more appropriate measure is the RMS, typically twice 
as large.
</li>
</ol>        
<p>
        For a Gaussian PSF, $f(\delta)=\frac{1}{2\pi \sigma} \exp(-\frac{\delta^2}{2\sigma^2})$, these values are
        respectively $2\sigma$ and $\sqrt{2}\sigma$.
        <br>PSF filenames: %(psf_files)s
        """
        psf = self.get_psf(irfname)
 
        def bkg_size(e, ct):
            f2 = lambda delta: psf(e,ct, delta)**2 * 2*np.pi*delta
            return np.degrees(1./np.sqrt(np.pi*integrate.quad(f2, 0, np.inf)[0]))
        
        def loc_size(e,ct, eps=1e-4):
            def wphi(r,x):
                fun = lambda phi: np.log(psf(e,ct,np.sqrt(r**2 -2*r*x*np.cos(phi) + x**2)))
                return 2*integrate.quad(fun, 0, np.pi)[0]
 
            def we(r):
                fun = lambda rp : psf(e,ct,rp) * rp * wphi(r,rp)
                return integrate.quad(fun, 0, np.radians(3.0))[0]

            #f3 = lambda theta: psf(e,ct, theta) * theta**3 * 2.*np.pi
            #return np.degrees(np.sqrt(integrate.quad(f3, 0, np.pi/6) [0]))
            return eps * np.sqrt(1/(2*(we(0)- we(np.radians(eps)))))

        egev = np.logspace(-1.+1/8., 2.5+1/8., 3.5*4+1)
        #egev = np.logspace(-1, 2, 100)
        front, back = [[bkg_size(e*1e3,ct) for e in egev] for ct in range(2)]
        floc, bloc = [[loc_size(e*1e3,ct) for e in egev] for ct in range(2)]
        f68,b68  = [[psf.inverse_integral(e*1e3, ct) for e in egev] for ct in range(2)]
        fig,ax = plt.subplots(figsize=(6,6))
        for x, s, label in zip((front, back, floc, bloc, f68, b68),
                            ('-g', 'r', '--g', '--r', ':g', ':r'),
                            ('front bkg', 'back bkg','front loc', 'back loc', 'front 68', 'back 68')):
            ax.plot(egev, x, s, lw=2, label=label)
        
        plt.setp(ax, xlabel='Energy (GeV)', ylabel='PSF size (deg)', xscale='log', yscale='log',
            xlim=(0.1, 100), ylim=(0.05, 10), title=title)
        ax.legend(prop=dict(size=10)); ax.grid()
        if outfile is None: return fig
        self.psf_df = pd.DataFrame(dict(front=front, floc=floc, back=back, bloc=bloc,f68=f68,b68=b68), 
                index=egev.round(3))
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
    