"""
Manage sources: single class SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sourcelist.py,v 1.10 2011/10/03 22:04:11 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""
import types
import numpy as np

def set_extended_property(esource):
    """ kluge until proper design for sources
        create property 'spectral_model' with get, set methods do the original smodel of an extended source
    """
    def getter(source): return source.extended_source.smodel
    def setter(source, newmodel): source.extended_source.smodel = newmodel
    esource.__class__.spectral_model = property(getter, setter,doc='spectral model')
def set_diffuse_property(dsource):
    def getter(source): return source.diffuse_source.model
    def setter(source, newmodel): source.diffuse_source.model = newmodel
    dsource.__class__.spectral_model = property(getter, setter,doc='spectral model')
def set_point_property(psource) : 
    def getter(source): return source.model
    def setter(source, newmodel): source.model = newmodel
    psource.__class__.spectral_model = property(getter, setter,doc='spectral model')

class Parameter(object):
    """ Manage a parameter """
    def __init__(self, sourcelist, source, model, parname):
        self.sourcelist=sourcelist
        self.source = source
        self.parname = parname
  
class SourceList(list):
    """ manage properties of the list of sources
        
    """

    def __init__(self, roi):
        """ initialization: get sources from the ROIAnalysis object: creates a polymorphic list of the 
            three types: point, global, and extended
            
            roi : an ROIdef object:
                expect to find attributes:
                    bgmodels : a list of global and extended sources
                    point_sources : list of point sources
                    roi_dir: the center of the ROI
        """
        def copy_defaults(a, b):
            """ if object ahas a dict 'defaults' copy corresponding stuff to b
            """
            #if not hasattr(a, 'defaults'): return
            for key,value,text in a.defaults:
                b.__dict__[key]=a.__dict__[key]
                #print key,a.__dict__[key]
                
        class DiffuseSourceFactory(object):
            """ k create new object from diffuse source
            """
            def __init__(self, roi):   self.roi = roi
            def __call__(self, **kwargs):
                cm = self.roi.global_sources[source_index]
                kwargs.update(name=self.roi.name)
                src = cm()
                src.skydir = None
                return src
        class ExtendedSourceFactory(object):
            def __init__(self, roi):   self.roi = roi
            def __call__(self, source_index, **kwargs):
                cm = self.roi.extended_sources[source_index]
                src = cm.__class__(cm.sa, cm.extended_source, self.roi.roi_dir, **kwargs)
                copy_defaults(cm, src)
                src.skydir = cm.extended_source.skydir
                return src
        class PointSourceFactory(object):
            def __init__(self, roi):    self.roi =roi
            def __call__(self, source_index):
                return roi.point_sources[source_index], self.roi.roi_dir
                
        # note that sources are added in the order diffuse, point to agree with ROIAnalysis
        # also, add two attributes to each source object and define 'spectral_model'
        for i, source in enumerate(roi.global_sources):
            self.append(source)
        for i, source in enumerate(roi.extended_sources):
            source.manager_index = i#+len(roi.global_sources)
            source.factory = ExtendedSourceFactory(roi)
            # kluge: make a getter property for the spectral model
            set_extended_property(source)
            source.skydir = source.extended_source.skydir
            source.spatial_model = source.extended_source.spatial_model
            self.append(source)
        for i,source in enumerate(roi.point_sources):
            #source.manager_index = i
            #source.factory = PointSourceFactory(roi)
            set_point_property(source)
            self.append(source)
        
    # note that properties are dynamic, in case sources or their models change
    @property
    def source_names(self): return np.array([s.name for s in self])
    @property
    def models(self): return np.array([s.spectral_model for s in self])
    @property
    def free(self): return np.array([np.any(model.free) for model in self.models])
    
    def __str__(self):
        return 'Sourcelist, with %d sources, %d free parameters' %(len(self), sum(self.free))
        
    def find_source(self, source_name):
        names = [s.name for s in self]
        try:
            return self[names.index(source_name)]
        except:
            raise Exception('source %s not found' %source_name)
    
    @property
    def parameter_names(self):
        """ array of free parameter names """
        names = []
        for source_name, model in zip(self.source_names, self.models):
            for pname in np.array(model.param_names)[model.free]:
                names.append(source_name.strip()+'_'+pname)
        return np.array(names)
    
    def get_parameters(self):
        """ array of free parameters (fitter rep)"""
        if len(self.models)==0: return []
        return np.concatenate([m.get_parameters() for m in self.models])
    
    def set_parameters(self,parameters):
        """ set the (fitter rep) parameters"""
        current_position=0
        for m in self.models:
            cp,nn = current_position, current_position+ sum(m.free)
            m.set_parameters(parameters[cp:nn])
            current_position += nn-cp

    def set_covariance_matrix(self,cov_matrix,):
        """ take over-all covariance matrix from a fit, give block-diagonal pieces to model objects"""
        current_position=0
        for m in self.models:
            cp,nn = current_position,current_position+len(m.get_parameters())
            m.set_cov_matrix(cov_matrix[cp:nn,cp:nn])
            current_position += nn-cp

    def get_model_parameters(self):
        return 10**self.get_parameters()
    
    def set_model_parameters(self, par):
        self.set_parameters(np.log10(par))
    parameters = property(get_model_parameters, set_model_parameters,
                doc='array of free model parameters')

    @property
    def uncertainties(self):
        """ return relative uncertainties 
            just scale the log10 uncertaintes
        """
        variances = np.concatenate([m.cov_matrix.diagonal()[m.free] for m in self.models])
        variances[variances<0] = 0
        return np.sqrt(variances) /np.log10(np.e)

 