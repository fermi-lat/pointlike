"""
Set up an ROI (transitional?)

$Header$

"""
import os
import numpy as np
from . import dataset
from .. pipeline import skymodel, dataspec
from .. utilities import keyword_options
from .. like import roi_managers, roi_diffuse, roi_extended


class ROIfactory(object):
    """
    combine the dataset and skymodel for an ROI
    """
    defaults =(
        ('analysis_kw', dict(irf='P7SOURCE_V6',minROI=7,maxROI=7, emax=316277),'roi analysis'),
        ('log',   None, ''),
        ('convolve_kw', dict( resolution=0.125, # applied to OTF convolution: if zero, skip convolution
                        pixelsize=0.05, # ExtendedSourceConvolution
                        num_points=25), 'convolution'),# AnalyticConvolution
        ('selector', skymodel.HEALPixSourceSelector,' factory of SourceSelector objects'),
        ('quiet', False, 'set to suppress most output'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, dataname, skymodel, **kwargs):
        """ dataname : string
                used to look up data specification
            skymodel: object SkyModel object
        """
        keyword_options.process(self, kwargs)
        self.skymodel = skymodel
        self.dataset = dataset.DataSet(dataname, **self.analysis_kw)
        self.exposure = self.dataset.exposure #needed by roi_diffuse
        self.psf  = self.dataset.psf
        
    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        show = """analysis_kw selector""".split()
        for key in show:
            s += '\t%-20s: %s\n' %(key,
                self.__dict__[key] if key in self.__dict__.keys() else 'not in self.__dict__!')
        return s

    def _diffuse_sources(self, src_sel):
        """ return a source manager for the diffuse,  global and extended sources
        """
        skydir = src_sel.skydir()
        # get all diffuse models appropriate for this ROI
        globals, extended = self.skymodel.get_diffuse_sources(src_sel)
       
        # perform OTF convolutions with PSFs: first diffuse, then extended if any
        # note that the wrapper could be roi_diffuse.ROIDiffuseModel_PC (for pre-convolved) it takes a tolerance
        # if resolution is zero, assume preconvolved
        def otf_diffuse_mapper( source):
            res = self.convolve_kw['resolution']
            if res>0:
                return roi_diffuse.ROIDiffuseModel_OTF(self, source, skydir, pixelsize=res)
            return roi_diffuse.ROIDiffuseModel_PC(self, source, skydir)
        def diffuse_mapper(source): #pre-convolved if isotrop
            if source.name.startswith('iso'):
                return roi_diffuse.ROIDiffuseModel_PC(self, source, skydir)
            return otf_diffuse_mapper(source)
        global_models = map( diffuse_mapper, globals)

        def extended_mapper( source):
            return roi_extended.ROIExtendedModel.factory(self,source,skydir)
        extended_models = map(extended_mapper, extended)
        
        # create and return the manager
        return roi_managers.ROIDiffuseManager(global_models+extended_models, skydir, quiet=self.quiet)
        
    def _local_sources(self, src_sel):
        """ return a manager for the local sources with significant overlap with the ROI
        """
        ps = self.skymodel.get_point_sources(src_sel)
        skydir = src_sel.skydir()
        return roi_managers.ROIPointSourceManager(ps, skydir,quiet=self.quiet)

    def roi(self, *pars, **kwargs):
        """ return an object based on the selector, with properties for creating full roi analysis
        """
        roi_kw = kwargs.pop('roi_kw',None)
        src_sel = self.selector(*pars, **kwargs)

        class ROIdef(object):
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)
            def __str__(self):
                return 'ROIdef for %s' %self.name
        skydir = src_sel.skydir()   

        bands=self.dataset(skydir)
        psm = self._local_sources(src_sel)
        dsm = self._diffuse_sources(src_sel)

        return ROIdef( name=src_sel.name() ,
                    roi_dir = skydir, 
                    bands=bands, 
                    psm=psm, 
                    dsm=dsm,)

def main(indir='uw26', dataname='P7_V4_SOURCE_4bpd', skymodel_kw={}):
    sm = skymodel.SkyModel(indir,  **skymodel_kw)
    rf = ROIfactory(dataname, sm)
    return rf.roi(749)