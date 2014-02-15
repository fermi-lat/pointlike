"""
Environment plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/environment.py,v 1.15 2014/02/14 18:46:59 burnett Exp $

"""

import os, pickle, types
import numpy as np
import pylab as plt
import pandas as pd
from scipy import integrate, misc, optimize
#from skymaps import SkyDir #?

from . import roi_info
from .. import (configuration, diffuse)

class Environment(roi_info.ROIinfo):
    """ Environment plots"""

    
    def setup(self, **kw):
        super(Environment, self).setup(**kw)
        self.plotfolder='environment'
        self.config = configuration.Configuration('.',postpone=True, quiet=True)
        s=[]
        for sname in ['ring', 'isotrop']:
            z=[]
            for i,r in self.df.iterrows():
                j = r['diffuse_names'].index(sname)
                z.append(r['counts']['models'][j][1][:16])
            s.append(z)
        
        #s = [[ x[1]['counts']['models'][modelnumber][1][:16] for x in self.df.iterrows()] for modelnumber in range(2)]
        self.bdf = [pd.DataFrame(y, index=self.df.index) for y in s]

        
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
        if irfname is None: irfname=self.config.irf
        cdm = pycaldb.CALDBManager(irf=irfname)
        self.psf_files=cdm.get_psf()
        return pypsf.CALDBPsf(cdm)

    def psf_plot(self, irfname=None, outfile='psf.csv', title=''):
        r"""PSF size
        
        <br>Plots of two science-driven measures of the size of the PSF, compared with the standard 68 percent containment.
        If  $P(\delta)$ is the normalized PSF as a function of deviation $\delta$ for a given energy and conversion type, 
        in principle averaged over incidence angle, but here for normal incidence, the two cases are:
        <ol>
        <li>Discrimination of a point source 
        signal in the presence of a uniform background; the <em>average</em> PSF is the inverse of the solid angle, with radius
        $$\begin{equation}
        1 / \sqrt{\pi \int_0^\pi P(\delta)^2 2\pi \sin\delta \ \mathrm{d}\delta}
        \end{equation}$$
 <li> Measurement of the position of a source. Assuming no background, this requires the expected value for the log of the
 PSF, which is the likelihood, as a function of $\delta$
$$\begin{equation} 
w(\delta) = \int_0^\pi  \sin\theta\ \mathrm{d} \theta\ P(\theta) \int_0^{2\pi} \mathrm{d}\phi \ln P\left(\sqrt{\delta^2 -2\delta\  \theta \cos\phi+\theta^2}\right)
\end{equation}$$
The resolution is the curvature, or second derivative of this function evaluated at $\delta=0$. The curvature at $\delta=0$ is
$$\begin{equation} 
\frac{\partial^2 w(\delta)}{\partial \delta^2} = \pi \int_0^\pi \sin\theta \ \mathrm{d}\theta \frac{P'(\theta)^2}{P(\theta)}
\end{equation}$$
where $P'(\theta) = \frac{\partial P(\theta)}{\partial \theta}$. This is consistent with equation (A2) in the 2FGL paper, 
for the case with no background.
</li>
</ol>  
<p>Of course the first measure is relevant in the background-dominated case, below a GeV or so, 
the second when there is small background, above a few GeV.    
<p>
        For a Gaussian PSF, $P(\delta)=\frac{1}{2\pi \sigma^2} \exp(-\frac{\delta^2}{2\sigma^2})$, these values are
        respectively $2\sigma$ and $\sigma$. (The 68 percent is in between, at 1.5 $\sigma$.)
        <br><br>PSF filenames: %(psf_files)s
        """
        psf = self.get_psf(irfname)
 
        def bkg_size(e, ct):
            f2 = lambda delta: psf(e,ct, delta)**2 * 2*np.pi*delta
            return np.degrees(1./np.sqrt(np.pi*integrate.quad(f2, 0, np.inf)[0]))
        
#        def loc_size(e,ct, eps=1e-4):
#            def wphi(r,x):
#                fun = lambda phi: np.log(psf(e,ct,np.sqrt(r**2 -2*r*x*np.cos(phi) + x**2)))
#                return 2*integrate.quad(fun, 0, np.pi)[0]
# 
#            def we(r):
#                fun = lambda rp : psf(e,ct,rp) * rp * wphi(r,rp)
#                return integrate.quad(fun, 0, np.radians(3.0))[0]
#
#            #f3 = lambda theta: psf(e,ct, theta) * theta**3 * 2.*np.pi
#            #return np.degrees(np.sqrt(integrate.quad(f3, 0, np.pi/6) [0]))
#            return eps * np.sqrt(1/(2*(we(0)- we(np.radians(eps)))))
        def loc_size(e, ct):
            func = lambda x : psf(e,ct, x)
            fprime = lambda x : misc.derivative(func, x, dx=0.0001, order=5)
            integrand = lambda rp : rp * fprime(rp)**2/func(rp) * np.pi
            return np.degrees(1/np.sqrt(integrate.quad(integrand, 0, np.radians(5))[0]))
            
        
        egev = np.logspace(-1.+1/8., 2.5+1/8., 3.5*4+1)
        front, back = [[bkg_size(e*1e3,ct) for e in egev] for ct in range(2)]
        floc, bloc = [[loc_size(e*1e3,ct) for e in egev] for ct in range(2)]
        f68,b68  = [[psf.inverse_integral(e*1e3, ct) for e in egev] for ct in range(2)]
        fig,ax = plt.subplots(figsize=(6,6))
        for x, s, label in zip((front, back, floc, bloc, f68, b68),
                            ('-g', 'r', '--g', '--r', ':g', ':r'),
                            ('front bkg', 'back bkg','front loc', 'back loc', 'front 68', 'back 68')):
            ax.plot(egev, x, s, lw=2, label=label)
        
        plt.setp(ax, xlabel='Energy (GeV)', ylabel='PSF size (deg)', xscale='log', yscale='log',
            xlim=(0.1, 100), ylim=(0.02, 8), title=title)
        ax.legend(prop=dict(size=10)); ax.grid()
        ax.set_xticklabels('0.1 1 10 100'.split())
        ax.set_yticklabels('0.01 0.1 1'.split())
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
        # look up filenames used to define the isotropic spectrum: either new or old diffuse spec; list or dict
        #diffuse=self.config['diffuse']
        #isokey = 'isotrop' if type(diffuse)==types.DictType else 1
        #self.idfiles = [os.path.join(os.environ['FERMI'],'diffuse',diffuse[isokey][i]) for i in (0,1)]
        df=diffuse.diffuse_factory(self.config.diffuse['isotrop'])
        self.idfiles = [x.fullfilename for x in df]
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

    def limb_map(self, energy=100):
        """Limb plots
        """
        df=diffuse.diffuse_factory(config.diffuse['limb'])

        [dd.plot_map(energy, scale='linear', cbtext='flux', title='Limb %s, %d MeV'%(name,energy))\
            for dd,name in zip(df, self.config.event_type_names)]
    
    def limb_flux(self, energy=100, ra=0):
        """Limb flux 
        Note assume indpenent of RA: this is for RA=0.
        """
        from skymaps import SkyDir
        df=diffuse.diffuse_factory(self.config.diffuse['limb'])
        names = self.config.event_type_names
        sindec = np.linspace(-1,1,501)
        dec = np.degrees(np.arcsin(sindec));
        fig, ax = plt.subplots(figsize=(8,5))
        for i in range(len(df)):
            flux = np.array([df[i](SkyDir(ra,d), energy) for d in dec])
            ax.plot(sindec, flux, lw=2, label=names[i])
        ax.legend(); ax.grid();
        plt.setp(ax, xlabel='sin(Dec)', ylabel='flux', title='Limb flux @ %.0f MeV' %energy)
        return fig
    
    def get_background(self, roi):
        roiname='HP12_%04d' % roi
        return [t.ix[roiname] for t in self.bdf]
    
    def diffuse_flux(self, rois=[0,888]):
        """Diffuse flux
        Predicted counts for the low latitude and high latitude ROIs.
        """
        fig, ax = plt.subplots(1,1, figsize=(6,6), dpi=150, sharey=True)
        egev = np.array(self.energy)/1e3
        if rois is None: rois = self.rois

        for r in rois:
            gal, iso = self.get_background(r)
            ax.plot(egev, gal, '-d', label='gal %d'%r)
            ax.plot(egev, iso, '--o', label='iso %d'%r)
        plt.setp(ax, xscale='log', xlim=(0.1,300), xlabel='Energy (GeV)',
            yscale='log',  ylabel='Diffuse counts/ROI')
        ax.legend(prop=dict(size=10)); ax.grid()
        return fig
        
        
    def all_plots(self, **kw):
        self.runfigures([self.psf_plot, self.exposure_plots, 
            self.isotropic_spectrum,self.diffuse_flux, self.limb_flux, ])
    