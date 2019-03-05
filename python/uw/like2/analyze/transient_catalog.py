"""
Analysis of a set of monthly transients to generate a subset suitable for a catalog

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/transient_catalog.py,v 1.6 2017/11/18 22:26:38 burnett Exp $

"""

import os, pickle, glob, json
from astropy.io import fits as pyfits
from uw.like2.analyze import (transientinfo, analysis_base)
from uw.like2.tools import DateStamp
from uw.like2 import to_fits
from skymaps import Band, SkyDir

import numpy as np
import pylab as plt
import pandas as pd

from analysis_base import (html_table, FloatFormat,)

class TransientCatalog(analysis_base.AnalysisBase):
    """Transient Catalog analysis
    
    <p>This analysis checks the combined monthly transient source detections, and prepares a summary file 
    The individual monthly list is 
    <a href="http://glast-ground.slac.stanford.edu/Decorator/exp/Fermi/Decorate/groups/catalog/pointlike/skymodels/%(skymodel)s/plot_index.html?skipDecoration"/> here </a>
    <br>Output from setup.
    <pre>%(logstream)s</pre>
    """
    def setup(self, model_dir=None, saved_filename=None, **args): 
        if model_dir is None:
            model_dir=os.getcwd()
        else:
            os.chdir(model_dir)
        self.plotfolder='catalog' #needed by superclass
        assert os.path.exists('config.json'), 'File "config.json" not found in current folder {}'.format(model_dir)
        config=json.load(open('config.json'))
        self.input_model = config['input_model']
        self.nmonths = config['nmonths']
        self.startlog()
        print 'input model: {}, with {} months'.format(self.input_model, self.nmonths)
        if saved_filename is not None:
            filename='transients.csv'
            dfsaved = pd.read_csv(filename, index_col=0)
            print 'Read in file with all monthly detections {}, {} entries'.format(filename, len(dfsaved))
            self.df =df =dfsaved[dfsaved.locqual<8]; 
        
            
        print 'Running transientinfo.Analysis to collect input source detections'
        sa = transientinfo.Analysis(all_months_model=self.input_model)
        self.df=df=sa.df
        df['skydir'] = map(SkyDir, df.ra, df.dec)
        self.ts25=self.df.ts>25
        self.hilat = np.abs(self.df.glat)>10
        print 'TS>25: {}, hilat: {}, both:{}'.format(sum(self.ts25), sum(self.hilat), sum(self.ts25 & self.hilat))
 
        self.si = si= transientinfo.SourceInfo(sa)
        si.df['skydir'] = sa.df6y.skydir
        self.reference = si.df[si.df.ngood>0]
        print 'Selected reference set ({}/{}) of input sources '\
            'such that at least one month has a "good" detection.'.format(len(self.reference), len(si.df))
        self.logstream=self.stoplog()
        assert os.getcwd() == model_dir, 'Changed cwd: {}'.format(os.getcwd())
               
    def source_detections(self):
        """Source detection
        Examine how many times the input model sources are detected in individual months
        "detected" means only that the source had a TS>10.  "good" requires also a good localization. 
        This is the reference set of sources from the original model, with at least one "good" month detection.
        """
        fig, ax = plt.subplots(figsize=(6,6))
        N = self.nmonths
        histkw = dict(bins=np.linspace(0,N, N+1), histtype='step', lw=2)
        ax.hist(self.reference.nmonths, label='detected',  **histkw)
        ax.hist(self.reference.ngood, label='good',color='orange', **histkw)
        plt.setp(ax, xlabel='Number of months', xlim=(0,N))
        ax.grid(True, alpha=0.5)
        ax.legend()
        return fig
        
    def transient_summary_plots(self):
        """Transient summary
        Examination of some values for the detected transient sources.
        <br><b>Upper left</b>: TS distribution. The dashed green line shows the shape of the null distribution. 
        <br><b>Upper right</b> Check that all months are equivalent
        <br><b>Lower left</b> sine of galactic latitude. Note that a larger fraction of the TS>25 source detections
        are in the galactic plane
        <br><b>Lower Right</b> The localization quality, cut at 8.
        """
        
        fig, axx = plt.subplots(2,2, figsize=(12,12))
        histkw=dict(log=True, histtype='step', lw=2)
        axf = axx.flatten()
        df = self.df; hilat=self.hilat; ts25=self.ts25

        ax = axf[0]
        ax.hist(df.ts, np.logspace(1,3, 41), label='all', **histkw );
        ax.hist(df.ts[hilat], np.logspace(1,3, 26), label='hilat', color='orange', **histkw) ;
        plt.setp(ax, xscale='log', xlabel='TS', ylim=(0.9,None));
        ax.axvline(25, color='orange')
        x = np.logspace(1,2)
        ax.plot(x, 10*np.exp(-x/2.)/np.exp(-25/2), '--', label='null')
        ax.grid(True, alpha=0.5)
        ax.legend()

        ax=axf[1]
        ax.hist(df.month, np.linspace(0.5,48.5,49), label='all', **histkw );
        ax.hist(df.month[ts25], np.linspace(0.5,48.5,49), label='TS>25', **histkw);
        ax.hist(df.month[hilat], np.linspace(0.5,48.5,49), label='hilat', color='orange',**histkw);
        plt.setp(ax, xlabel='month', xlim=(0.5, 48.5) ,ylim=(0.9,1e4))
        ax.grid(True, alpha=0.5)
        ax.legend()

        ax=axf[2]
        df['singlat'] = np.sin(np.radians(np.array(df.glat, float))) #not sure why necessary
        ax.hist(df.singlat, np.linspace(-1,1,41),label='all', **histkw)
        ax.hist(df.singlat[ts25], np.linspace(-1,1,41), label='TS>25', **histkw)
        plt.setp(ax, xlabel='sin(glat)' ,ylim=(0.9,None));
        ax.grid(True, alpha=0.5)
        ax.legend();

        ax=axf[3]
        ax.hist(df.locqual, np.linspace(0,8,17), label='all', **histkw)
        ax.hist(df.locqual[ts25], np.linspace(0,10,21), label='TS>25', **histkw)
        ax.hist(df.locqual[hilat], np.linspace(0,10,21), label='hilat',color='orange', **histkw)
        plt.setp(ax, xlabel='localization quality',ylim=(0.9,None))
        ax.grid(True, alpha=0.5)
        ax.legend();
        return fig
   
    def close_cut(self, tsmin=[10,25], mindist=0.5, bins=np.linspace(0,4,51), log=True):
        """Closest distance
        Examine distance to the closest reference source. Since such sources are in the model used to create the TS map to 
        detect transient sources, there should be no correlation. However, when the input model source is weak, fluctuations
        on the positions of photons may be detected as separate sources. Thus we will make a %(mindist_cut)0.1f degree cut.
         
        """
        def closest(a, b):
            return np.degrees(min([a.difference(x) for x in b]))
        self.mindist_cut = mindist
        
        print 'Comparing {} sources with {} in reference set'.format(len(self.df), len(self.reference))
        self.df['closediff'] = np.array([closest(a,self.reference.skydir) for a in self.df.skydir])
  
        fig,ax=plt.subplots(figsize=(6,6))
        for tscut in tsmin:
            ax.hist(self.df.closediff[self.df.ts>tscut]**2, bins=bins, log=log, histtype='step', lw=2,
                    label='TS>{:.0f}'.format(tscut))
        plt.setp(ax, xlabel='Distance**2 [deg**2]', title='Closest distance to reference source',)
        if log: ax.set_ylim(bottom=0.9)

        ax.axvline(mindist**2, color='red')
        ax.legend()
        ax.grid(True, alpha=0.5)
        return fig
  

    def write_to_fits(self, filename='sources_4yr_transients', error_box_factor=1.1, error_box_add=5e-3, 
           ):
        """FITS output log
        <pre>%(fitslogstream)s</pre>"""
        self.startlog()
        print 'Writing all sources to file {}'.format (filename+'.csv')
        self.df.to_csv(filename+'.csv')
        cuts='(sources.ts>25) & (sources.a<0.25) &(sources.closediff>%.2f)' % self.mindist_cut
        print '\nRunning "to_fits"...'
        self.fits_file = filename+'.fits'
        to_fits.main(self.fits_file,  cuts=cuts,
                     localization_systematic = (error_box_factor, error_box_add)
                     )

        self.fitslogstream=self.stoplog()
        
      
    def all_plots(self):
        self.runfigures([self.transient_summary_plots, self.source_detections,   self.close_cut, self.write_to_fits,
            ])
        
        
class CloseCut(object):
    """
    Manage selection of a set of sources, by removing those close to another set
    """
    def __init__(self, sources, reference):
        """ sources, reference : SkyDir lists
        """
        print 'Comparing {} sources with {} in reference set'.format(len(sources), len(reference))
        def closest(a, b):
            return np.degrees(min([a.difference(x) for x in b]))
        self.closediff = np.array([closest(a,reference) for a in sources])
    
    def plot(self, marker=0.5, bins=np.linspace(0,1,101), ax=None, log=False, label=None):
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(figsize=(6,6))
        else: fig=ax.figure
        ax.hist(self.closediff**2, bins, histtype='step', lw=2, log=log, label=label);
        if marker is not None:
            ax.axvline(marker**2,color='red', label='{:.2f} deg'.format(marker))
        plt.setp(ax, xlabel='Distance**2 [deg**2]', title='Closest distance to reference source',)
        if log: ax.set_ylim(bottom=0.9)
        return fig
    
    def __call__(self, mindist=0.5):
        return self.closediff>mindist;
