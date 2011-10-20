"""
Set up an ROI (transitional?)

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roisetup.py,v 1.2 2011/10/01 13:35:06 burnett Exp $

"""
import os
import numpy as np
from . import dataset, skymodel, diffuse
from .. utilities import keyword_options, convolution
from .. like import  roi_diffuse, roi_extended


class ROIfactory(object):
    """
    combine the dataset and skymodel for an ROI
    
    """
    defaults =(
        ('analysis_kw', dict(irf='P7SOURCE_V6',minROI=7,maxROI=7, emax=316277),'roi analysis'),
        ('skymodel_kw', {}, 'skymodel keywords'),
        ('convolve_kw', dict( resolution=0.125, # applied to OTF convolution: if zero, skip convolution
                            pixelsize=0.05, # ExtendedSourceConvolution
                            num_points=25), # AnalyticConvolution
                                    'convolution parameters'),
        ('selector', skymodel.HEALPixSourceSelector,' factory of SourceSelector objects'),
        ('cache',   None, 'folder to look for cached diffuse'),
        ('quiet', False, 'set to suppress most output'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, indir, dataname, **kwargs):
        """ 
        parameters
        ----------
            dataname : string
                used to look up data specification
            indir: folder containing skymodel definition
        """
        keyword_options.process(self, kwargs)
        self.skymodel = skymodel.SkyModel(indir,  **self.skymodel_kw)
        self.dataset = dataset.DataSet(dataname, **self.analysis_kw)
        self.exposure = self.dataset.exposure #needed by roi_diffuse
        self.psf  = self.dataset.psf
        convolution.AnalyticConvolution.set_points(self.convolve_kw['num_points'])
        convolution.ExtendedSourceConvolution.set_pixelsize(self.convolve_kw['pixelsize'])

    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        show = """analysis_kw selector""".split()
        for key in show:
            s += '\t%-20s: %s\n' %(key,
                self.__dict__[key] if key in self.__dict__.keys() else 'not in self.__dict__!')
        return s

    def _diffuse_sources(self, src_sel):
        """ return  the diffuse, global and extended sources
        """
        skydir = src_sel.skydir()
        # get all diffuse models appropriate for this ROI
        globals, extended = self.skymodel.get_diffuse_sources(src_sel)
       
        global_models = [diffuse.mapper(self, src_sel.name(), skydir, source, cache=self.cache) for source in globals]

        def extended_mapper( source):
            return roi_extended.ROIExtendedModel.factory(self,source,skydir)
        extended_models = map(extended_mapper, extended)
        return global_models, extended_models

    def _local_sources(self, src_sel):
        """ return the local sources with significant overlap with the ROI
        """
        ps = self.skymodel.get_point_sources(src_sel)
        return np.asarray(ps)

    def roi(self, *pars, **kwargs):
        """ return an object based on the selector, with attributes for creating  roi analysis:
            list of ROIBand objects
            list of models: point sources and diffuse
        pars, kwargs : pass to the selector
        """
        roi_kw = kwargs.pop('roi_kw',None)
        src_sel = self.selector(*pars, **kwargs)

        class ROIdef(object):
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)
            def __str__(self):
                return 'ROIdef for %s' %self.name
        skydir = src_sel.skydir()  
        global_sources, extended_sources = self._diffuse_sources(src_sel)
        return ROIdef( name=src_sel.name() ,
                    roi_dir=skydir, 
                    bands=self.dataset(skydir), 
                    global_sources=global_sources,
                    extended_sources = extended_sources,
                    point_sources=self._local_sources(src_sel))

def main(indir='uw27', dataname='P7_V4_SOURCE_4bpd',irf='P7SOURCE_V6' ,skymodel_kw={}):
    rf = ROIfactory(indir, dataname, **skymodel_kw)
    return rf