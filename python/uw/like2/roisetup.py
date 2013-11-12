"""
Set up an ROI factory object

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roisetup.py,v 1.36 2013/10/29 03:25:20 burnett Exp $

"""
import os, sys, types, pickle, zipfile
import numpy as np
import pandas as pd
import skymaps
from .. utilities import keyword_options, convolution
from .. like import roi_extended

from . import (dataset, skymodel, diffusedict as diffuse, configuration, extended , sources)
        


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


class ROImodel(object):
    """ construct the model, or list of sources, for an ROI, reading the zipped pickle format
    """
    defaults = (
        ('quiet', True, 'set False for info'),
        ('ecat',  None, 'If present, use for catalog'),
        )
    @keyword_options.decorate(defaults)
    def __init__(self, config, roi_index, **kwargs):
        """config : configuration.Configuration object
            used to find model info
        roi_index : integer
            ROI index, from 0 to 1727
        """
        keyword_options.process(self, kwargs)
        self.pickle_file = os.path.join(config.configdir, 'pickle.zip')
        assert os.path.exists(self.pickle_file)
        self.config = config
        self.index = roi_index
        self._z = zipfile.ZipFile(os.path.expandvars(self.pickle_file))
        if self.ecat is None: #speed up if already loaded
            self.ecat = extended.ExtendedCatalog(self.config.extended, quiet=self.quiet)
        self.local_sources, self.global_sources = self.load_sources(roi_index)
        self.neighbor_sources = []
        for neighbor_index in self.neighbors():
            self.neighbor_sources += self.load_sources(neighbor_index, neighbors=True)
        
    def __repr__(self):
        return '%s.%s : %d local, %d neighbor, %d global sources at ROI %d' \
            % (self.__module__, self.__class__.__name__, len(self.local_sources), len(self.neighbor_sources),
            len(self.global_sources), self.index)
            
    def prep_for_likelihood(self):
        """ return a list of sources, and a list of the ones with free parameters for BandLike 
            order as global, local, neighbors
            """
        tmp = self.global_sources
        tmp += self.local_sources
        tmp += self.neighbor_sources
        free = np.array([np.any(src.free) for src in tmp])
        return tmp, free
        
    def load_sources(self, index, neighbors=False):
        """ select and return lists of sources in the given HEALPix 
        Tags each with an index property
        if neigbors is True, return only local sources, otherwise a tuple global,local
        also set the free list to False for neighbors. (The Source constructor sets it as a 
        copy of the model's free list)
        """

        def load_local_source(name, rec):
            if not rec['isextended']:
                src = sources.PointSource(name=name, skydir=rec['skydir'], model=rec['model'])
                
            else:
                src = self.ecat.lookup(name)
                src.model = rec['model']
            if neighbors: src.free[:]=False # not sure this is necessary
            src.index = index
            return src
        
        def load_global_source(name, rec):
            if not self.quiet: print 'Loading global source %s for %d' % (name, index)
            gsrc = sources.GlobalSource(name=name, skydir=None, model=rec,
                dmodel = diffuse.diffuse_factory(self.config.diffuse[name]))
            gsrc.index =index
            return gsrc

        p = pickle.load(self._z.open('pickle/HP12_%04d.pickle' % index))
        local_sources = [load_local_source(name, rec) for name,rec in p['sources'].items()]
        if neighbors: return local_sources
        global_sources = [load_global_source(name, rec) for name, rec \
            in zip(p['diffuse_names'], p['diffuse']) if name not in self.ecat.names]
        return local_sources, global_sources
        
    def neighbors(self):
        from pointlike import IntVector
        b12 = skymaps.Band(12)
        v = IntVector()
        b12.findNeighbors(int(self.index),v) 
        return list(v)

    
def main(modeldir='P202/uw29', skymodel_kw={}):
    rf = ROIfactory(modeldir, **skymodel_kw)
    return rf
