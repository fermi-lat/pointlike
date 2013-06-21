"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/pulsarseedcheck.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

from . import seedcheck

class PulsarSeedCheck(seedcheck.SeedCheck):
    require='pseedcheck.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname= 'pseedcheck'
        self.spectral_type = 'exponential cutoff'
        self.load()

    def all_plots(self):
        """ Results of analysis of pulsar seeds
        %(info)s
        """
        self.info = self.df.describe().to_html()
        self.runfigures([self.seed_cumulative_ts, self.unassoc_seed_cumulative_ts, self.spectral_parameters, self.localization],
                ('pulsar_cumulative_ts', 'pulsar_unassoc_cumulative_ts', 'pulsar_spectral_pars', 'pulsar_localization'))