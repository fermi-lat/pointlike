"""
Description here

$Header: /phys/users/glast/python/uw/like2/analyze/sourcetotal.py,v 1.144 2013/06/18 12:35:36 burnett Exp $

"""

from . import roi_info

class SourceTotal(roi_info.ROIinfo):
    def setup(self, **kw):
        super(SourceTotal, self).setup(**kw)
        self.plotfolder='sourcetotal'
        self.source_name='sources'
        self.title='Sources'
        self.funcs = [self.counts_map]
        self.fnames=['source_counts']

    def all_plots(self, **kwargs):
        """ Counts for all sources, per RIO"""
    
        self.runfigures([self.counts_map], ['source_counts'], **kwargs)