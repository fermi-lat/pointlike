"""
Set up an ROI factory object

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roisetup.py,v 1.35 2013/10/13 13:55:58 burnett Exp $

"""
import os, sys, types
import numpy as np
import pandas as pd
import skymaps
from . import dataset, skymodel, diffusedict, configuration #, roi_bands
from .. utilities import keyword_options, convolution
from .. like import roi_extended

        


class ROIfactory(object):
    """
    combine the dataset and skymodel for an ROI
    
    """
    defaults =(
        ('analysis_kw', dict(irf=None,minROI=5,maxROI=5, emin=100, emax=316277, quiet=False),'roi analysis keywords'),
        ('skymodel_kw', {}, 'skymodel keywords'),
        ('convolve_kw', dict( resolution=0.125, # applied to OTF convolution: if zero, skip convolution
                            pixelsize=0.05, # ExtendedSourceConvolution
                            num_points=25), # AnalyticConvolution
                                    'convolution parameters'),
        ('irf', None,  'Set to override saved value with the skymodel: expect to find in custom_irf_dir'),
        ('extended', None, 'Set to override saved value with skymodel'),
        ('selector', skymodel.HEALPixSourceSelector,' factory of SourceSelector objects'),
        ('postpone', True, 'postpone file loading'),
        ('data_interval', 0, 'Data interval (e.g., month) to use'),
        ('nocreate', True, 'Do not allow creation of a binned photon file'),
        ('quiet', False, 'set to suppress most output'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, modeldir, **kwargs):
        """ 
        parameters
        ----------
        modeldir: folder containing skymodel definition
        """
        keyword_options.process(self, kwargs)
        self.analysis_kw['quiet']=self.quiet
        if not self.quiet:  'ROIfactory setup: \n\tskymodel: ', modeldir
        # extract parameters used by skymodel for defaults

        self.config = configuration.Configuration(modeldir, quiet=self.quiet, postpone=self.postpone)
        self.skymodel = skymodel.SkyModel(modeldir, **self.skymodel_kw)
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
        assert skydir is not None, 'should use the ROI skydir'
        
        # get all diffuse models appropriate for this ROI
        try:
            global_models, extended = self.skymodel.get_diffuse_sources(src_sel)

        except Exception, msg:
            print self.config.dataset, msg
            raise

        def extended_mapper( source):
            if not self.quiet:
                if not self.quiet:  'constructing extended model for "%s", spatial model: %s' \
                    %(source, source.spatial_model.__class__.__name__)
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
        roi_kw = kwargs.pop('roi_kw',dict())
        # allow parameter to be a name or a direction
        sel = pars[0]
        source_name=None
        if type(sel)==types.IntType:
            index = int(sel) # needs to be int if int type
        elif type(sel)==skymaps.SkyDir:
            index = self.skymodel.hpindex(sel)
        elif type(sel)==types.StringType:
            index = self.skymodel.hpindex(self.skymodel.find_source(sel).skydir)
            source_name=sel
        elif type(sel)==tuple and len(sel)==2: # interpret a tuple of length 2 as (ra,dec)
            index = self.skymodel.hpindex(skymaps.SkyDir(*sel))
        else:
            raise Exception( 'factory argument "%s" not recognized.' %sel)
        ## preselect the given source after setting up the ROI
        ## (not implemented here)
        #
        src_sel = self.selector(index, **kwargs)

        class ROIdef(object):
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)
            def __repr__(self):
                return 'ROIdef for %s, at %s' % (self.name, self.roi_dir)
        skydir = src_sel.skydir()  
        #global_sources, extended_sources = self._diffuse_sources(src_sel)
        global_sources, extended_sources = self.skymodel.get_diffuse_sources(src_sel)

        return ROIdef(name = src_sel.name, 
                    roi_dir=skydir,
                    bands = self.config.get_bands(skydir, **self.analysis_kw),
                    global_sources=global_sources,
                    extended_sources = extended_sources,
                    point_sources=self._local_sources(src_sel),
                    config = self.config,
                    **roi_kw)

                
                    
    def __call__(self, *pars, **kwargs):
        """ alias for roi() """
        return self.roi(*pars, **kwargs)
        
    def reload_model(self):
        """ Reload the sources in the model """
        self.skymodel._load_sources()


def main(modeldir='P202/uw29', skymodel_kw={}):
    rf = ROIfactory(modeldir, **skymodel_kw)
    return rf
