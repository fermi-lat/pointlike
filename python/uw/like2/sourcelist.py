"""
Manage sources for likelihood: single class SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sourcelist.py,v 1.31 2013/04/09 21:26:27 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""
import types
import numpy as np
from uw.like import SpatialModels # for setting all spatial models to use template

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

## redundant?
def set_default_bounds( model, force=False):
    """
    Handy utility to set bounds for a model from like.Models
    force=True to override previously set bounds.
    """
    if not force and hasattr(model, 'bounds'):
        # model has bounds. Were they set? check to see if all are None
        notset =  np.all(np.array([np.all(b ==[None,None]) for b in model.bounds]))
        if not notset: return
    bounds=[]
    def to_internal(fun, values):
        return [fun(value) if value is not None else None for value in values]
    for pname, mp in zip(model.param_names, model.mappers):
        plim = (None,None)
        try:
            plim = dict(
                Index=(-0.5, 5), 
                Norm=(10**-16, 10**-7),
                Scale=(0.001, 4.0),
                beta=(-0.1, 5.), 
                Cutoff=(100., 1e5),
                )[pname.split('_')[0]]
        except: pass
        bounds.append( to_internal(mp.tointernal, plim) )
    model.bounds = np.array(bounds) # convert to array so can mask with free

def check_bounds(model):
    """ check that free parameters in the model are within the bounds
        could freeze the paramet
    """
    for pname, par, mp, bound, free in zip(model.param_names, model.parameters, model.mappers, model.bounds, model.free):
        if free:
            if bound[0] is not None:
                if mp.toexternal(bound[0])>=par:
                    raise SourceListException('Model %s, parameter %s (%s) fails lower bound' % (model.name, pname, par))
            if bound[1] is not None:
                if mp.toexternal(bound[1])<=par:
                    raise SourceListException('Model %s, parameter %s (%s) fails upper bound' % (model.name, pname, par))
                
    
  
class SourceList(list):
    """ manage properties of the list of sources
        
    """

    def __init__(self, roi):
        """ initialization: get sources from the ROIAnalysis object: creates a polymorphic list of the 
            three types: point, global, and extended
            
            roi : an ROIdef object:
                expect to find attributes:
                    global_sources : a list of global sources
                    extended_sources : lost of local extended sources
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
                
                
        # note that sources are added in the order diffuse, point to agree with ROIAnalysis
        # also, add two attributes to each source object and define 'spectral_model'
        def append(source):
            assert source.name not in self.source_names, 'attempt to add name %s twice' % source.name
            assert hasattr(source, 'spectral_model'), 'bug'
            self.set_default_bounds(source, True) # override until figure out more subtle
            self.append(source)
        for source in roi.global_sources:
            append(source)
        for i, source in enumerate(roi.extended_sources):
            source.manager_index = i
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
        
    # note that the following properties are dynamic, in case sources or their models change interactively
    @property
    def source_names(self): return np.array([s.name for s in self])
    @property
    def models(self): return np.array([s.spectral_model for s in self])
    
    @property
    def free(self): 
        """ mask which defines variable sources: all global and local sources with at least one variable parameter 
        """
        def free_source(s):
            return np.any(s.spectral_model.free) or s.skydir is None
        return np.array([ free_source(s) for s in self])
        
    @property 
    def bounds(self):
        """ fitter representation of applied bounds """
        return np.concatenate([m.bounds[m.free] for m in self.models])

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
            self.selected_source_index =-1
            raise SourceListException('source %s not found' %source_name)
        def found(s):
            self.selected_source=s
            self.selected_source_index = names.index(s.name)
            return s
        if source_name[-1]=='*':
            for name in names:
                if name.startswith(source_name[:-1]): 
                    return found(self[names.index(name)])
            not_found()
        if source_name[0]=='*':
            for name in names:
                if name.endswith(source_name[1:]): 
                    return found(self[names.index(name)])
            not_found()
        try:
            selected_source = self[names.index(source_name)]
            #if self.selected_source is None or self.selected_source != selected_source:
            #    print 'selected source %s for analysis' % selected_source.name
            return found(selected_source)
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

    def set_covariance_matrix(self, cov_matrix, select=None):
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
                models[i].internal_cov_matrix[j,j] = cov_matrix[n,n]
            return
        current_position=0
        for m in self.models:
            cp,nn = current_position,current_position+len(m.get_parameters())
            m.set_cov_matrix(cov_matrix[cp:nn,cp:nn])
            current_position += nn-cp

    def get_model_parameters(self):
        if len(self.models)==0: return []
        return np.concatenate([m.free_parameters for m in self.models])
        #return 10**self.get_parameters()
    
    def set_model_parameters(self, par):
        assert False, 'is this used?'
        self.set_parameters(np.log10(par))
    model_parameters = property(get_model_parameters, set_model_parameters,
                doc='array of free model parameters')

    @property
    def uncertainties(self):
        """ return relative uncertainties 
        """
        variances = np.concatenate([m.get_cov_matrix().diagonal()[m.free] for m in self.models])
        variances[variances<0]=0
        return np.sqrt(variances) / (np.abs(self.model_parameters) +1e-6) #avoid divide by zero

    def set_default_bounds(self, source, force=False):
        model = source.spectral_model
        if not force and hasattr(model, 'bounds'):
            # model has bounds. Were they set? check to see if all are None
            notset =  np.all(np.array([np.all(b ==[None,None]) for b in model.bounds]))
            if not notset: return

        #if not force and hasattr(model, 'bounds'): return
        bounds=[]
        def to_internal(fun, values):
            return [fun(value) if value is not None else None for value in values]
        for pname, mp in zip(model.param_names, model.mappers):
            plim = (None,None)
            try:
                plim = dict(
                    Index=(-0.5, 5), 
                    Norm=(10**-16, 10**-7),
                    Scale=(0.0001, 10.0),
                    beta=(-0.1, 5.), 
                    Cutoff=(100., 1e5),
                    )[pname.split('_')[0]]
            except: pass
            # override for diffuse
            if source.name.startswith('ring') and pname=='Norm': 
                plim = (0.5, 1.5)
            # finally, save internal reps
            bounds.append( to_internal(mp.tointernal, plim) )
        model.bounds = np.array(bounds) # convert to array so can mask with free
    
        
    def add_source(self, source):
        if source.name in self.source_names:
            raise SourceListException('Attempt to add source "%s": already exists' % source.name)
        set_point_property(source)
        self.set_default_bounds(source)
        self.append(source)
 
    def del_source(self, source_name):
        source = self.find_source(source_name) # first get it
        self.remove(source)
        return source
        