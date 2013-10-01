"""
Limb plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/limb.py,v 1.4 2013/08/22 18:22:00 burnett Exp $

"""

import numpy as np
import pylab as plt

from skymaps import SkyDir
from . import roi_info

class Limb(roi_info.ROIinfo):
    "Limb plots"
    def setup(self, **kwargs):
        super(Limb, self).setup(**kwargs)
        self.plotfolder='limb'
        self.source_name='limb'
        self.title='Limb'
        dm = self.diffuse_models('limb')
        self.fpar,self.bpar = [np.array([m[i] if m else np.nan for m in dm] )for i in range(2)]

        self.default_plots()
        self.funcs += [self.bpar_plot, self.fpar_plot, self.flux_vs_dec, ]
        self.fnames+= ['limb_polar_back', 'limb_polar_front', 'limb_flux_vs_dec',]

     
    def polar_plots(self, values, title=None,
                vmin=0, vmax=2, vticks=5, vlabel=None, thetamax=60):
        """
        values : array of float
        Creates a Figure that must be cleared
        """
        fig, axes = plt.subplots(1,2, figsize=(8,4),subplot_kw=dict(polar=True))
        plt.subplots_adjust(bottom=0.12, top=0.90, wspace=0.3)
        plot_kw = dict(s=120, edgecolor='none', vmin=vmin, vmax=vmax)
        
        galeq = [SkyDir(float(u),0, SkyDir.GALACTIC) for u in np.arange(0,360,1)]
        galra = np.array([sd.ra() for sd in galeq])
        galdec = np.array([sd.dec() for sd in galeq])
        def radius( dec, i): return [90-dec, 90+dec][i] 
        ra = np.array(map(lambda dir: dir.ra(), self.df.skydir))
        dec = np.array(map(lambda dir: dir.dec(), self.df.skydir))
        for i, ax in enumerate(axes[:2]):
            cut = [dec>(90-thetamax), dec<(thetamax-90)][i]
            r = [90-dec, 90+dec][i]
            theta = np.radians(ra)
            c = np.array(values)[cut]
            sc =ax.scatter(theta[cut], radius(dec[cut],i), c=c, **plot_kw);
            galtheta = galra; galrad = radius(galdec,i) 
            ax.plot(np.radians(galtheta[galrad<thetamax]),galrad[galrad<thetamax], '-', color='grey', lw=2)
            ax.set_ylim(ymax=thetamax)
            ax.set_title(['North','South'][i], ha='right', fontsize='small')

        cbax = fig.add_axes((0.25,0.08,0.5, 0.04))
        cb=plt.colorbar(sc, cbax, orientation='horizontal')
        if vlabel is not None: cb.set_label(vlabel)
        if vmin is not None and vmax is not None:
            cb.set_ticks(np.linspace(vmin, vmax, vticks))
        if title is not None: plt.suptitle(title)  
        return fig
        
        
    def flux_vs_dec(self, ymax=3):
        """ front and back flux vs dec
        Plots of front and back flux normalizations, ploting ROIS with |b|>35.
        <br>The dashed lines are piecewise functions, fit to the observed pass 7 4-year mean.
        <br>The spectral function is described by the diffuse configuration entry for the limb: "%(limbspect)s"
        """
        class PieceWise(object):
            """ functiod that is a piecewise set of straight lines"""
            def __init__(self, a,b):
                self.a, self.b =a,b
                self.n = len(a)
                self.s = [(b[i+1]-b[i])/(a[i+1]-a[i]) for i in range(self.n-1)]
            def __call__(self, x):
                if x<=self.a[0]: return self.b[0]
                for i in range(self.n-1):
                    if x<self.a[i+1]:
                        return self.b[i]+(x-self.a[i])*self.s[i]
                return self.b[-1]
        limbfun = dict(front =PieceWise([-1., -0.4, 0.4, 1.0],[0.75, 0, 0, 0.75]),
                back= PieceWise([-1., -0.7, -0.5, 0.5, 0.7, 0.85, 1.0],
                                [2.0,  0.5, 0,    0,   0.5,  1.2, 0.9])    )
        ra = np.array(map(lambda dir: dir.ra(), self.df.skydir))
        dec = np.array(map(lambda dir: dir.dec(), self.df.skydir))
        dom = np.linspace(-1,1,201) 
        try:
            config = eval(open('config.txt').read())
            self.limbspect = config['diffuse']['limb']
        except Exception, msg:
            self.limbspect='(Failed to retrieve: %s)'%msg
        #dm = self.diffuse_models('limb')
        #fpar,bpar = [np.array([m[i] if m else np.nan for m in dm] )for i in range(2)]
        
        
        fig, axx = plt.subplots(2,1, figsize=(14,6), sharex=True)
        plt.subplots_adjust(right=0.9)
        c=np.abs(self.df.glat)
        cut = c>35
        for ax, par, label  in zip(axx, [self.fpar,self.bpar], 'front back'.split()):
            scat=ax.scatter(np.sin(np.radians(dec))[cut], par[cut], c=c[cut], s=25,vmin=0, vmax=90, edgecolor='none')
            plt.setp(ax, xlim=(-1,1), xlabel='sin(dec)' if label=='back' else '',  ylim=(0,ymax))
            ax.plot(dom, map(limbfun[label],dom), '--', color='k', lw=1) 
            ax.grid()
            ax.text(-0.75, ymax-0.4, label, fontsize=18)
        fig.text(0.05, 0.5, 'flux normalization factor', rotation='vertical', va='center')
        cax = fig.add_axes((0.94, 0.25, 0.02, 0.4))
        cb=plt.colorbar(scat, cax)
        cb.set_label('abs(glat)')
        plt.suptitle('Limb observed flux') 
        return fig
    
    def bpar_plot(self):
        """ Back normalization
        """
        return self.polar_plots(self.bpar)
    def fpar_plot(self):
        """ front normalization
        """
        return self.polar_plots(self.fpar)
        
    def all_plots(self):
        super(Limb, self).all_plots()