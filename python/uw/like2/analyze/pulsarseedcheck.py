"""
Description here

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/pulsarseedcheck.py,v 1.1 2013/06/21 20:15:31 burnett Exp $

"""

from . import seedcheck

class PulsarSeedCheck(seedcheck.SeedCheck):
    require='pseedcheck.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname= 'pseedcheck'
        self.spectral_type = 'exponential cutoff'
        self.load()

    def all_plots(self):
        self.runfigures([self.seed_list, self.seed_cumulative_ts, self.unassoc_seed_cumulative_ts, self.spectral_parameters, self.localization],
                ('pulsar_seed_table', 'pulsar_cumulative_ts', 'pulsar_unassoc_cumulative_ts', 'pulsar_spectral_pars', 'pulsar_localization'))