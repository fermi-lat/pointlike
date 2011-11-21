"""
Set up an ROI 

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roisetup.py,v 1.5 2011/11/16 14:14:20 burnett Exp $

"""
import os
import numpy as np
import skymaps
from . import dataset, skymodel, diffuse
from .. utilities import keyword_options, convolution
from .. like import roi_extended, pypsf

        
class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime."""

    def __init__(self, sa):

        skymaps.EffectiveArea.set_CALDB(sa.CALDBManager.CALDB)
        skymaps.Exposure.set_cutoff(np.cos(np.radians(sa.thetacut)))
        inst = ['front', 'back']
        aeff_files = sa.CALDBManager.get_aeff()
        ok = [os.path.exists(file) for file in aeff_files]
        if not all(ok):
            raise DataSetError('one of CALDB aeff files not found: %s' %aeff_files)
        self.ea  = [skymaps.EffectiveArea('', file) for file in aeff_files]
        if sa.verbose: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        if sa.use_weighted_livetime:
            self.exposure = [skymaps.Exposure(sa.lt,sa.weighted_lt,ea) for ea in self.ea]
        else:
            self.exposure = [skymaps.Exposure(sa.lt,ea) for ea in self.ea]

    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)


class ROIfactory(object):
    """
    combine the dataset and skymodel for an ROI
    
    """
    defaults =(
        ('analysis_kw', dict(irf='P7SOURCE_V6',minROI=7,maxROI=7, emax=316277, quiet=False),'roi analysis'),
        ('skymodel_kw', {}, 'skymodel keywords'),
        ('convolve_kw', dict( resolution=0.125, # applied to OTF convolution: if zero, skip convolution
                            pixelsize=0.05, # ExtendedSourceConvolution
                            num_points=25), # AnalyticConvolution
                                    'convolution parameters'),
        ('selector', skymodel.HEALPixSourceSelector,' factory of SourceSelector objects'),
        ('quiet', False, 'set to suppress most output'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, indir, dataname, **kwargs):
        """ 
        parameters
        ----------
        indir: folder containing skymodel definition
        dataname : string
                used to look up data specification
        """
        keyword_options.process(self, kwargs)
        self.skymodel = skymodel.SkyModel(indir,  **self.skymodel_kw)
        self.dataset = dataset.DataSet(dataname, **self.analysis_kw)
        self.exposure  = ExposureManager(self.dataset)
        self.psf = pypsf.CALDBPsf(self.dataset.CALDBManager)
 
        #self.exposure = self.dataset.exposure 
        #self.psf  = self.dataset.psf
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
        globals, extended = self.skymodel.get_diffuse_sources(src_sel)
       
        global_models = [diffuse.mapper(self, src_sel.name(), skydir, source) for source in globals]

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
        roi_kw = kwargs.pop('roi_kw',dict())
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
                    bands=self.dataset(self.psf, self.exposure, skydir), 
                    global_sources=global_sources,
                    extended_sources = extended_sources,
                    point_sources=self._local_sources(src_sel), **roi_kw)

def main(indir='uw27', dataname='P7_V4_SOURCE_4bpd',irf='P7SOURCE_V6' ,skymodel_kw={}):
    rf = ROIfactory(indir, dataname, **skymodel_kw)
    return rf