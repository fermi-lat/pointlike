"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/pgwseedcheck.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

from . import seedcheck

class PGWSeedCheck(seedcheck.SeedCheck):
    require='seedcheck_PGW.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname= 'seedcheck_PGW'
        self.spectral_type = 'power law'
        self.load()