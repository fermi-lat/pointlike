"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/limbrefit.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

import os, pickle
import numpy as np
import pandas as pd

from skymaps import SkyDir
from . import limb

class LimbRefit(limb.Limb):
    def setup(self, **kw):
        self.plotfolder = 'limb_refit'
        self.source_name='limb'
        self.title='Limb refit'
        files, pickles = self.load_pickles('limb')
        rdict = dict()
        for f,p in zip(files,pickles):
            name = os.path.split(f)[-1][:9]
            front, back = p['model'].parameters if 'model' in p else (np.nan, np.nan)
            if 'model' not in p: p['model']= None
            ra, dec = p['ra'], p['dec']
            skydir = SkyDir(ra,dec)
            glat,glon = skydir.b(), skydir.l()
            rdict[name] = dict(ra=ra, dec=dec, glat=glat, glon=glon,skydir=skydir, front=front, back=back)
        self.df = pd.DataFrame(rdict).transpose()
        self.fpar = self.df.front
        self.bpar = self.df.back
        
    def all_plots(self):
        """ Special run to refit Limb normalization
        %(table)s
        """
        self.table = pd.DataFrame([self.df.front, self.df.back, ], 
                index=['front', 'back']).T.describe().to_html()
        self.runfigures([self.flux_vs_dec], ['limb_fit_norm_vs_dec'] )