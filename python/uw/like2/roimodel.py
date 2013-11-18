"""
Set up and manage the model for all the sources in an ROI

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roimodel.py,v 1.1 2013/11/13 06:36:24 burnett Exp $

"""
import os ,zipfile, pickle
import numpy as np
import skymaps
from .. utilities import keyword_options
from . import (sources, 
      diffusedict as diffuse,
      )
class ROImodelException(Exception):pass

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


class ROImodel(list):
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
        self.load_sources(roi_index)
        for neighbor_index in self.neighbors():
            self.load_sources(neighbor_index, neighbors=True)
        self.selected_source = None
        self.initialize()

     
    def initialize(self):
        """For fast parameter access: must be called if any source changes
        """
        class Parameters(object):
            """ Manage the free parameters in the model, as a virtual array
            
            Note that if a parameter in a source model is changed from its current value,
            that source is marked; its 'changed' property is set True
            """
            def __init__(self, sources):
                models = sources.models
                self.free_sources = [source for source in sources if np.any(source.model.free)]
                self.clear_changed()
                self.ms = t = [(source, sum(source.model.free)) for source in self.free_sources]
                ss=[]; ii=[]
                for (s,k) in t:
                    for j in range(k):
                        ss.append(s) #free_models[s])
                        ii.append(j)
                self.index = np.array([ss, ii])
            
            def __getitem__(self, i):
                source, k = self.index[:,i]
                return source.model.get_parameters()[k]
            
            def __setitem__(self,i,x):
                source, k = self.index[:,i]
                model = source.model
                pars = model.get_parameters()
                if x==pars[k]: return
                pars[k] = x
                source.changed=True
                model.set_parameters(pars)
            
            def get_all(self):
                return np.concatenate([s.model.get_parameters() for s in self.free_sources])
                
            def set_all(self, pars):
                i =0
                for source, n in self.ms:
                    j = i+n
                    model = source.model
                    oldpars = model.get_parameters()
                    newpars = pars[i:j]
                    if np.any(oldpars != newpars):
                        source.model.set_parameters(newpars)
                        source.changed=True
                    i =j
                    
            def __repr__(self):
                return '%d parameters' % (self.index).shape[1]
            def clear_changed(self):
                for s in self.free_sources:
                    s.changed=False
            @property
            def dirty(self):
                return np.array([s.changed for s in self.free_sources])

        self.parameters = Parameters(self)
  
        
    def __repr__(self):
        return '%s.%s : %d global, %d local, %d total sources for ROI %d' \
            % (self.__module__, self.__class__.__name__, self.global_count,self.local_count,  
            len(self), self.index)
            
        
    def load_sources(self, index, neighbors=False):
        """ select and add sources in the given HEALPix to self.
        Tags each with an index property
        if neigbors is True, add only local sources.
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
            src.model.free = src.free # Don't think I still need this second copy of free
            # set parameter bounds, ignoring equivalent code in like.Models
            set_default_bounds( src.model )
            return src
        
        def load_global_source(name, rec):
            if not self.quiet: print 'Loading global source %s for %d' % (name, index)
            gsrc = sources.GlobalSource(name=name, skydir=None, model=rec,
                dmodel = diffuse.diffuse_factory(self.config.diffuse[name]))
            gsrc.index =index
            return gsrc

        p = pickle.load(self._z.open('pickle/HP12_%04d.pickle' % index))
        if not neighbors:
            global_sources = [load_global_source(name, rec) for name, rec \
                in zip(p['diffuse_names'], p['diffuse']) if name not in self.ecat.names]
            self.global_count = len(global_sources)
            for s in global_sources: self.append(s)
        local_sources = [load_local_source(name, rec) for name,rec in p['sources'].items()]
        if not neighbors: self.local_count = len(local_sources)
        for s in local_sources: self.append(s)
        
    def neighbors(self):
        from pointlike import IntVector
        b12 = skymaps.Band(12)
        v = IntVector()
        b12.findNeighbors(int(self.index),v) 
        return list(v)
    
    # note that the following properties are dynamic, in case sources or their models change interactively
    @property
    def source_names(self): return np.array([s.name for s in self])
    @property
    def models(self): return np.array([s.model for s in self])
    @property
    def free(self): 
        """ mask which defines variable sources: all global and local sources with at least one variable parameter 
        """
        return np.array([ np.any(s.model.free) for s in self])
    @property 
    def bounds(self):
        """ fitter representation of applied bounds """
        return np.concatenate([m.bounds[m.free] for m in self.models])
    
        
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

    def find_source(self, source_name):
        """ Search for the source with the given name
        
        source_name : string or None
            if the first or last character is '*', perform a wild card search, return first match
            if None, and a source has been selected, return it
        """
        if source_name is None:
            if self.selected_source is None:
                raise ROImodelException('No source is selected')
            return self.selected_source
        names = [s.name for s in self]
        def not_found():
            self.selected_source_index =-1
            raise ROImodelException('source %s not found' %source_name)
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

    def add_source(self, source):
        """Add a point source
        Must be a sources.PointSource object
        """
        
        if source.name in self.source_names:
            raise ROImodelException('Attempt to add source "%s": a source with that name already exists' % source.name)
        set_default_bounds(source.model)
        self.append(source)
 
    def del_source(self, source_name):
        source = self.find_source(source_name) # first get it
        self.remove(source)
        return source
