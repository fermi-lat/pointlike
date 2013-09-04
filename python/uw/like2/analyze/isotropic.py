"""
Isotropic diffuse plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/isotropic.py,v 1.2 2013/06/21 20:51:41 burnett Exp $

"""

import os, types
import numpy as np
import pylab as plt

from . import galactic

class Isotropic(galactic.Galactic):
    """Isotropic diffuse analysis
    """

    def setup(self, **kw):
        super(Isotropic, self).setup(**kw)
        self.plotfolder='iso'
        self.source_name='isotrop'
        self.title='Isotropic'
        self.default_plots()
        self.funcs += [self.isotropic_spectrum]
        self.fnames +=['isotropic_spectrum']
        # look up filenames used to define the isotorpic spectrum: either new or old diffuse spec; list or dict
        config = eval(open('config.txt').read())
        diffuse=config['diffuse']
        isokey = 'isotrop' if type(diffuse)==types.DictType else 1
        self.idfiles = [os.path.join(os.environ['FERMI'],'diffuse',diffuse[isokey][i]) for i in (0,1)]
        
    def isotropic_spectrum(self, other=None):
        """ Isotropic Spectrum from template
        
        The spectrum used to define the isotropic diffuse component.
        <br>Files for front/back: %(idfiles)s
        """
        nf,nb = map(np.loadtxt, self.idfiles)
        energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
        fig, axs = plt.subplots(1,2, figsize=(7,3), dpi=50)
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
        return fig
     
    def combined_spectra(self, other='isotrop_4years_P7_v9_repro_data_source.txt'):
        """ Special isotropic spectrum plot"""
        nf,nb = map(np.loadtxt, self.idfiles)
        fig, ax = plt.subplots(figsize=(5,5))
        energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
        ax.plot(energies, front*energies**2, '-g', label='1year_P76R_source_front')
        ax.plot(energies, back*energies**2, '-r', label='1year_P76R_source_back')
        plt.setp(ax, xlabel='Energy', ylabel='flux*e**2', xscale='log')
        ax.set_title('isotropic diffuse spectra', fontsize='small')
        f =os.path.join(os.path.expandvars('$FERMI/diffuse'),other)
        assert os.path.exists(f)
        iso = np.loadtxt(f)
        energies, fluxes =iso[:,0], iso[:,1] 
        ax.plot(energies, fluxes*energies**2, '-', label ='4years_P7-v9_source')
        plt.setp(ax, xlim=(1e2,1e4))
        ax.legend();ax.grid(True)
        return fig
    
    def all_plots(self):
        super(Isotropic, self).all_plots()