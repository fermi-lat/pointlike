"""
Create a Pivot collection

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/collection.py,v 1.2 2013/08/01 14:24:56 burnett Exp $

"""
import os, pickle
import numpy as np
import pylab as plt
import pandas as pd

from uw.utilities import makepivot
from . import sourceinfo
from . analysis_base import FloatFormat
from . _html import HTMLindex

class Collection(sourceinfo.SourceInfo):
    """Catalog presented as a Pivot collection
    All the sources with TS>10 are shown.
    """
    def setup(self, **kw):
        super(Collection, self).setup(**kw)
        self.plotfolder='collection'
        
    def make_facets(self):    
        df = self.df[self.df.ts>10]
        pvdf = df['ra dec glat glon eflux e0'.split()]
        assoc = df.associations
        pvdf['aprob'] = [z['prob'][0] if z is not None else None for z in assoc]
        pvdf['asource'] = [':'.join([z['cat'][0], z['name'][0]]) if z is not None else None for z in assoc]
        
        try:
            catdf = pd.read_csv('plots/comparison_2FGL/comparison_2FGL.csv', index_col=0)
            pvdf['catdist']  = catdf.distance
            pvdf['catsource']= catdf.otherid
            pvdf['catok'] = catdf.distance<0.25
        except:
            print 'Did not find 2FGL comparison'
            
        self.csv_file = os.path.join(self.plotfolder, 'collection.csv')
        pvdf.to_csv(self.csv_file)
        print 'wrote special collection csv to %s' % self.csv_file

    def make_collection(self):
        """Pivot table with all sources
        Inputs are:
        <ul><li> a table of sources, a <a href="../../%(csv_file)s?download=True">csv file</a>
        <li><a href="../../%(sedfig_zip)s?download=true">zip file</a> containing all SED plots.
        </ul>
        <br>
        %(full_pivot_link)s
        """
        self.make_facets()
        
        assert os.path.exists(self.csv_file), 'Did not find the file %s'% self.csv_file
        self.sedfig_zip='sedfig.zip'
        try:
            pc =makepivot.MakeCollection('all sources %s' % self.skymodel, 'sedfig', self.csv_file)
            self.full_pivot_link = """\
                    <p>These can be examined with a 
                    <a href="http://%s/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
                    which requires Silverlight."""  % (makepivot.http_host, pc.cId)
        except Exception, msg: 
            print "**** Failed to make pivot table, perhaps need to run sedinfo first: %s" % msg

    def all_plots(self):
        self.runfigures([self.make_collection,])