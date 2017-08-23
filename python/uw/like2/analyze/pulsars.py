"""
Pulsar search analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/pulsars.py,v 1.1 2015/07/24 17:56:02 burnett Exp $

"""

import os
import astropy.io.fits as pyfits
import numpy as np
import pylab as plt
import pandas as pd

from skymaps import SkyDir
from . import sourceinfo
from analysis_base import html_table, FloatFormat

from uw.utilities import makepivot


class PulsarSearch(sourceinfo.SourceInfo):
    """Analysis related to finding pulsars
    """

    def setup(self, **kw):
        super(PulsarSearch, self).setup(**kw)
        self.plotfolder='pulsars'
        self.psr = np.asarray([s.startswith('PSR') for s in self.df.index],bool)
        self.ts_info()
    
    def ts_info(self):
        tsl=[]; tsm=[]; tsh=[]
        for i in range(len(self.df.index)):
            sedrec =self.df.ix[i].sedrec
            if sedrec is None:
                print 'Source %s has no SED' % self.df.ix[i].name
                continue
            tsa=sedrec['ts']
            tsl.append(sum(tsa[:4]))
            tsm.append(sum(tsa[4:8]))
            tsh.append(sum(tsa[8:]))
        self.df['ts_low']=tsl
        self.df['ts_med']=tsm
        self.df['ts_high']=tsh
    
    
    def efratio(self, e1=2000, e2=20000):
    
        """Examine energy flux ratio
        
        Ratio of the energy flux at 20 GeV to that at 2 GeV.
        The identified pulsar subset is shown.
        """
        df=self.df
        efr = df['eflux_ratio']=np.asarray([model(e2)/model(e1)*(e2/e1)**2 for model in self.df.model])
        fig,ax = plt.subplots(figsize=(5,5))
        xlim = (1e-2,10)
        dom = np.logspace( np.log10(xlim[0]),np.log10(xlim[1]) ,31)  
        ax.hist(efr.clip(*xlim), dom ,log=True);
        ax.hist(efr[self.psr].clip(*xlim), dom, log=True, color='orange', label='PSR');
        plt.setp(ax, xscale='log', xlabel='eflux(20 GeV)/eflux(2 GeV)')
        ax.legend()
        ax.grid();
        return fig
        
    def selection(self, curvature_cut=0.1, ts_cut=10):
        """Select candidates.
        
        %(selection_info)s
        """
        self.curvature_cut=curvature_cut
        self.ts_cut=ts_cut
        df=self.df
        probfun = lambda x: x['prob'][0] if not pd.isnull(x) else 0
        aprob = np.array([ probfun(assoc) for  assoc in self.df.associations])
        no3fgl = np.asarray([s is None for s in self.df.cat3fgl]);
        self.keep= keep = no3fgl &(~self.psr) \
            & (self.df.curvature>curvature_cut) & (self.df.ts>ts_cut) & (self.df.locqual<8) &(aprob<0.1)
        
        self.total=sum(keep)
        self.cvsname='pulsar_candidates.csv'
        t = self.df[keep]['ra dec glat ts pivot_energy pindex eflux_ratio curvature roiname'.split()]
        t.to_csv(self.cvsname)
        print 'wrote %d sources to %s' % (len(t), self.cvsname)
        self.selection_info="""\
        Cuts: non-3FGL, non-LAT PSR, association probability < 0.1, curvature>%(curvature_cut)s, TS>%(ts_cut)s<br>, 
        <br>Total:%(total)s
        <br>%(html_list)s
        <br>
         Link to csv format table:
         <a href="../../%(cvsname)s?download=true">%(cvsname)s</a></li>
        """
        self.html_list = html_table(t, name=self.plotfolder+'/candidates', 
                heading='<h4>%d Candidate pulsar sources</h4>' % len(t), 
                    float_format=FloatFormat(2))
        

    def no_curvature(self, prefix='S966', ts_high_cut=2):
        """Weak new sources with PW fits
        
        %(weak_list)s
        """
        df = self.df
        pcut = np.array([n.startswith(prefix) for n in df.index],bool); 
        cut = (df.ts>10) &  (df.locqual<8) & (df.curvature<0.01) & pcut & (df.ts_high<ts_high_cut) & (df.ts_low<5)
        t = self.df[cut]['ra dec glat ts pivot_energy pindex  fitqual locqual ts_low ts_med ts_high roiname'.split()]
        self.noc_df=t.sort_index(by='roiname')
        print 'selected %d %s sources' % (len(t), prefix)
        self.weak_list = html_table(t, name=self.plotfolder+'/weak_pl_sources', 
                heading='<h4>%d weak new power-law sources</h4>' % len(t), 
                    float_format=FloatFormat(2))

        

    def collection(self, name='pulsar candidates v2', refresh=False):
        """Pivot collection.
        
        %(pivot_link)s
        
        """
        from uw.utilities import makepivot
        cid =  263
        cid =makepivot.MakeCollection(name, 'sedfig', self.cvsname, refresh=refresh).cId
        self.pivot_link = """\
            <p>The data set, with SED images can be examined with a 
            <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
            which requires Silverlight."""  % cid

    def all_plots(self):
        self.runfigures([self.efratio,
            self.selection, 
            self.no_curvature,
            #  self.collection, 
        ]       )

    