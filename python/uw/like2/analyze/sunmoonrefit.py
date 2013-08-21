"""
Sun/moon refit

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/sunmoonrefit.py,v 1.1 2013/06/21 20:15:31 burnett Exp $

"""

import os
import numpy as np
import pandas as pd

from skymaps import SkyDir
from . import roi_info
from .analysis_base import FloatFormat

class SunMoonRefit(roi_info.ROIinfo):
    """ Refit to SunMoon model, showing normalization, change in likelihood
    %(table)s"""

    def setup(self, **kw):
        self.plotfolder = 'sunmoon_refit'
        self.title = 'SunMoon refit'
        files, pickles = self.load_pickles('sunmoon')
        rdict = dict()
        for f,p in zip(files,pickles):
            name = os.path.split(f)[-1][:9]
            model = p['model'] if 'model' in p else None
            norm = model.parameters[0] if model is not None else np.nan
            norm_unc = np.diag(model.get_cov_matrix())[0]**0.5 if model is not None else np.nan
            if 'model' not in p: p['model']= None
            ra, dec = p['ra'], p['dec']
            skydir = SkyDir(ra,dec)
            glat,glon = skydir.b(), skydir.l()
            if glon>180: glon-=360.
            rdict[name] = dict(ra=ra, dec=dec, glat=glat, glon=glon, skydir=skydir, 
                norm=norm, norm_unc=norm_unc,
                delta_likelihood = p['delta_likelihood'])
        self.df = pd.DataFrame(rdict).transpose()
    
    def sunmoon_normalization(self):
        """ Sun/Moon normalization 
        Note that it is defined for all directions
        """
        return self.skyplot_with_hist(self.df.norm, 'norm', 0.5, 1.5, (0.5,1.5))
        
    def sunmoon_loglikelihood(self):
        """ Sun/Moon log likelihood change
        The improvement in the log likelihood when Sun/Moon normalization is freed
        """
        return self.skyplot_with_hist(self.df.delta_likelihood, 'delta loglike', 0, 5, (0,5))

    def all_plots(self):
        self.table = pd.DataFrame([self.df.norm, self.df.norm_unc,self.df.delta_likelihood ], 
                index=['norm', 'norm_unc', 'delta_likelihood']).T.describe().to_html(float_format=FloatFormat(2))
        self.runfigures([self.sunmoon_normalization, self.sunmoon_loglikelihood])