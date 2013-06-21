"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/flux_corr_iso.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

from . import fluxcorr

class FluxCorrIso(fluxcorr.FluxCorr):

    require = 'fluxcorriso.zip'
    def setup(self, **kw):
        super(FluxCorrIso,self).setup(source_name='fluxcorriso', **kw)
        self.title='Source-isotropic diffuse flux dependence'
        self.diffuse_name='Isotropic'
        
        self.plotfolder='fluxcorriso'