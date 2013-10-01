"""
Plots with diffuse

$Header$
"""

import os, pickle
import numpy as np
import pylab as plt
import pandas as pd

from . import roi_info
from . import analysis_base

class DiffuseInfo(roi_info.ROIinfo): #roi_info.ROIinfo):
    """Diffuse flux plots
     <br>Study of the diffuse components to the all-sky analysis. 
    """
 
    require=  'diffuse_info.zip'
    def setup(self, **kwargs):
        super(DiffuseInfo, self).setup(**kwargs)
        self.plotfolder='diffuse_info'
        f, self.pp = self.load_pickles('diffuse_info')
        self.exposure =self.get_info('exposure')
        
    def get_info(self, name, n=16):
        return [np.array([p['bdlist'][i][name] for p in self.pp]) for i in range(n)]

    def diffuse_definition(self):
        """Definition of the diffuse components
        Extracted from the config.txt file. 
             <p>Filenames are relative to %(diffuse_path)s.
        %(diffuse_def_html)s
        """
        dd= eval(open('config.txt').read())['diffuse']
        self.diffuse_def_html ='<dl>\n' \
            + '\n'.join(['<dt>%s</dt><dd>  %s</dd>' % item for item in dd.items()])\
            + '\n</dl>'
        self.diffuse_path= os.path.expandvars('$FERMI/diffuse')
        return None

    
    def flux_plot(self, diffname='limb',  nk=2, ylim=(0,4e-7), yscale='linear'):
        """plots
        """
        fig, axx = plt.subplots( 2,1, figsize=(12,8), sharex=True, sharey=True)
        diffc = self.get_info(diffname)
        for ax, k in zip(axx, range(nk)):
            ax.plot(self.df.dec, diffc[2*k]/self.exposure[2*k], '.g', label='front')
            ax.plot(self.df.dec, diffc[2*k+1]/self.exposure[2*k+1], '.r', label='back')
            plt.setp(ax, xlim=(-90,90), xlabel='Dec', 
                ylim=ylim, ylabel='band flux', yscale=yscale,
                title='%.0f MeV' % self.energy[k],  )
            ax.set_xticks(range(-90,91,30))
            ax.grid(); ax.legend()
        return fig
        
    def exposure_plot(self, k=4, ylim=(0.3, 0.8), yscale='linear'):
        """Exposure
        Points are the exposure for each ROI, normalized to the average over the sky, 
        %(exposure_average).2e cm**2 s.
        """
        fig, ax = plt.subplots( 1,1, figsize=(12,4), sharex=True, sharey=True)
        diffc = self.get_info('exposure')
        df, db = diffc[2*k:2*k+2] 
        self.exposure_average = mean = (df+db).mean()
        ax.plot(self.df.dec, df/mean, '.g', label='front')
        ax.plot(self.df.dec, db/mean, '.r', label='back')
        plt.setp(ax, xlim=(-90,90), xlabel='Dec', 
            ylim=ylim, ylabel='Exposure ratio', yscale=yscale,
            title='%.0f MeV' % self.energy[k],  )
        ax.set_xticks(range(-90,91,30))
        ax.grid(); ax.legend(loc='upper left')
        return fig

    def limb_plot(self, **kw):
        """limb flux"""
        return self.flux_plot('limb', 2, **kw)

    def isotropic_plot(self, **kw):
        """isotropic flux"""
        return self.flux_plot('isotrop', 2, **kw)

    def sunmoon_plot(self, **kw):
        """SunMoon flux"""
        return self.flux_plot('SunMoon', 2, **kw)

    def all_plots(self): #, other_html=None):
        self.runfigures([self.diffuse_definition, self.limb_plot,self.isotropic_plot, self.sunmoon_plot,self.exposure_plot,])



