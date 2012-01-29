"""
Manage sources for likelihood: single class SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sourcelist.py,v 1.17 2012/01/24 14:41:22 burnett Exp $
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

class SourceListException(Exception):pass

  
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
        def append(source):
            assert source.name not in self.source_names, 'attempt to add name %s twice' % source.name
            self.append(source)
        for source in roi.global_sources:
            append(source)
        for i, source in enumerate(roi.extended_sources):
            source.manager_index = i
            source.factory = ExtendedSourceFactory(roi)
            set_extended_property(source)
            source.skydir = source.extended_source.skydir
            source.spatial_model = source.extended_source.spatial_model
            append(source)
        for source in roi.point_sources:
            set_point_property(source)
            append(source)
        self.selected_source = None
        
    def add_source(self, source):
        """ 
        Parameters
        ---------
        source : object created by sources.PointSource, e.g.
            sources.PointSource(name='test', skydir=SkyDir(100,2)) 
        """
        set_point_property(source)
        self.append(source)
        
    # note that properties are dynamic, in case sources or their models change
    @property
    def source_names(self): return np.array([s.name for s in self])
    @property
    def models(self): return np.array([s.spectral_model for s in self])
    @property
    def free(self): 
        """ actually all global and local sources """
        return np.array([np.any(s.spectral_model.free) or s.skydir is None for s in self])
    
    def __str__(self):
        return 'Sourcelist, with %d sources, %d free parameters' %(len(self), sum(self.free))
        
    def find_source(self, source_name):
        """ Search for the source with the given name
        source_name : string or None
            if the first or last character is '*', perform a wild card search, return first match
            if None, and a source has been selected, return it
        """
        if source_name is None:
            if self.selected_source is None:
                raise SourceListException('No source is selected')
            return self.selected_source
        names = [s.name for s in self]
        def not_found():
            raise SourceListException('source %s not found' %source_name)
        if source_name[-1]=='*':
            for name in names:
                if name.startswith(source_name[:-1]): return self[names.index(name)]
            not_found()
        if source_name[0]=='*':
            for name in names:
                if name.endswith(source_name[1:]): return self[names.index(name)]
            not_found()
        try:
            selected_source = self[names.index(source_name)]
            #if self.selected_source is None or self.selected_source != selected_source:
            #    print 'selected source %s for analysis' % selected_source.name
            self.selected_source = selected_source
            return self.selected_source
        except:
            self.selected_source = None
            not_found()
   
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

    def set_covariance_matrix(self,cov_matrix, select=None):
        """ take over-all covariance matrix from a fit, give block-diagonal pieces to model objects
        if select is a list of variable indices, set diagonal elements only
        
        """
        if select is not None:
            #assert cov_matrix.shape[0] >= select[-1], 'cov_matrix.shape=%s'%cov_matrix.shape
            # upack diagonal errors to corresponding elements in 
            models = self.models
            t = [(i,np.arange(len(m.free))[m.free] ) for i,m in enumerate(models) if sum(m.free)>0]
            s = [ (i,j) for i,mf in t for j in mf]
            for n,k in enumerate(select):
                i,j = s[k]
                models[i].cov_matrix[j,j] = cov_matrix[n,n]
            return
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

    @property
    def bounds(self):
        ret = []
        for source_name, model in zip(self.source_names, self.models):
            for pname in np.array(model.param_names)[model.free]:
                pname = pname.split('_')[0]
                plim = (None,None)
                try:
                    plim = dict(
                        Index=(-3, 0.6),
                        Norm=(-15, -7),
                        Scale=np.log10(np.array([0.001, 4.0])),
                        beta=(-3, 1),
                        Cutoff=(2,6),
                        )[pname.split('_')[0]]
                except: pass
                # override for diffuse
                if source_name.startswith('ring') and pname=='Norm': 
                    plim = np.log10(np.array([0.5, 1.5]))
                ret.append(plim)
                
        return ret
        
    def add_source(self, source):
        if source.name in self.source_names:
            raise SourceListException('Attempt to add source "%s": already exists' % source.name)
        set_point_property(source)
        self.append(source)
 
    def del_source(self, source_name):
        source = self.find_source(source_name) # first get it
        self.remove(source)
        return source
        