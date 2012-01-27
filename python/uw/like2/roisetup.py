"""
Set up an ROI factory object

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roisetup.py,v 1.7 2011/12/06 22:14:08 burnett Exp $

"""
import os, types
import numpy as np
import skymaps
from . import dataset, skymodel, diffuse
from .. utilities import keyword_options, convolution
from .. like import roi_extended, pypsf

        
class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime.
    
    """

    def __init__(self, dataset, **datadict): 
        """
        Parameters
        ----------
        dataset :  DataSet object
            for CALDB, aeff, some parameters
            
        datadict['exposure-correction'] : list of strings defining functions of energy
            the correction factors to apply to front, back 
        """

        skymaps.EffectiveArea.set_CALDB(dataset.CALDBManager.CALDB)
        skymaps.Exposure.set_cutoff(np.cos(np.radians(dataset.thetacut)))
        inst = ['front', 'back']
        aeff_files = dataset.CALDBManager.get_aeff()
        ok = [os.path.exists(file) for file in aeff_files]
        if not all(ok):
            raise DataSetError('one of CALDB aeff files not found: %s' %aeff_files)
        self.ea  = [skymaps.EffectiveArea('', file) for file in aeff_files]
        if dataset.verbose: print ' -->effective areas at 1 GeV: ', \
                ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        if dataset.use_weighted_livetime:
            self.exposure = [skymaps.Exposure(dataset.lt,dataset.weighted_lt,ea) for ea in self.ea]
        else:
            self.exposure = [skymaps.Exposure(dataset.lt,ea) for ea in self.ea]

        correction = datadict.pop('exposure_correction', None)
        if correction is not None:
            self.correction = map(eval, correction)
            energies = [100, 1000, 10000]
            print 'Exposure correction: for energies %s ' % energies
            for i,f in enumerate(self.correction):
                print ('\tfront:','\tback: ')[i], map( f , energies)
        else:
            self.correction = lambda x: 1.0, lambda x: 1.0
            
    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)*self.correction[event_class](energy)
        
class ExposureCorrection(object):
    """ logarithmic interpolation function
    """
    def __init__(self, a,b, ea=100, eb=300):
        self.c = (b-a)/np.log(eb/ea)
        self.d =  a -self.c*np.log(ea)
        self.a, self.b = a,b
        self.ea,self.eb = ea,eb
    def __call__(self, e):
        if e>self.eb: return self.b
        if e<self.ea: return self.a
        return self.c*np.log(e) + self.d
    def plot(self, ax=None, **kwargs):
        import pylab as plt
        if ax is None: 
            ax = plt.gca()
        dom = np.logspace(1.5, 2.5, 51) 
        ax.plot(dom, map(self, dom), **kwargs)
        ax.set_xscale('log')
            

class ROIfactory(object):
    """
    combine the dataset and skymodel for an ROI
    
    """
    defaults =(
        ('analysis_kw', dict(irf='P7SOURCE_V6',minROI=7,maxROI=7, emin=100, emax=316277, quiet=False),'roi analysis keywords'),
        ('skymodel_kw', {}, 'skymodel keywords'),
        ('convolve_kw', dict( resolution=0.125, # applied to OTF convolution: if zero, skip convolution
                            pixelsize=0.05, # ExtendedSourceConvolution
                            num_points=25), # AnalyticConvolution
                                    'convolution parameters'),
        ('selector', skymodel.HEALPixSourceSelector,' factory of SourceSelector objects'),
        ('quiet', False, 'set to suppress most output'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, modeldir, dataspec, **kwargs):
        """ 
        parameters
        ----------
        modeldir: folder containing skymodel definition
        dataspec : string or dict
                used to look up data specification
                if string, equivalent to dict(dataname=dataspec); otherwise the dict must have
                a dataname element
        """
        keyword_options.process(self, kwargs)
        print 'ROIfactory setup: \n\tskymodel: ', modeldir
        self.skymodel = skymodel.SkyModel(modeldir,  **self.skymodel_kw)
        datadict = dict(dataname=dataspec) if type(dataspec)!=types.DictType else dataspec
        print '\tdatadict:', datadict
        self.dataset = dataset.DataSet(datadict['dataname'], **self.analysis_kw)
        self.exposure  = ExposureManager(self.dataset, **datadict)
        self.psf = pypsf.CALDBPsf(self.dataset.CALDBManager)
 
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
                    point_sources=self._local_sources(src_sel), 
                    exposure=self.exposure,
                    **roi_kw)
    def __call__(self, *pars, **kwargs):
        """ alias for roi() """
        return self.roi(*pars, **kwargs)

def main(indir='uw27', dataname='P7_V4_SOURCE_4bpd',irf='P7SOURCE_V6' ,skymodel_kw={}):
    rf = ROIfactory(indir, dataname, **skymodel_kw)
    return rf