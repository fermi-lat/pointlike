"""
Manage a set of parameters

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/parameterset.py,v 1.6 2016/03/21 18:54:12 burnett Exp $

"""
import os, types 
import numpy as np

class ParameterSet(object):
    """ Manage the free parameters in the ROI model, as a virtual array
    
    Notes:
        if a parameter in a source model is changed from its current value,
        that source is marked; its 'changed' property is set True
        
        The values of the paramters are *internal*
    """
    def __init__(self, sources, **kw):
        """sources : set of sources.Source objects
        """
        self.free_sources = [source for source in sources if np.any(source.model.free)]
        # dangerous? self.clear_changed()
        # make two indexing arrays:
        #  ms : list of (source, npar) for each source
        #  index :  (source, index within that source) for each parameter
        self.ms = t = [(source, sum(source.model.free)) for source in self.free_sources]
        ss=[]; ii=[]
        for (s,k) in t:
            for j in range(k):
                ss.append(s) 
                ii.append(j)
        self.index = np.array([ss, ii])
        self.mask = np.ones(len(ss),bool)
    
    def __getitem__(self, i):
        """ access the ith parameter, or all parameters with [:] """
        if isinstance(i, slice):
            if i==slice(None,None,None):
                return self.get_parameters()
            else:
                raise Exception('slice format not supported')
        else:
            source, k = self.index[:,i]
        return source.model.get_parameters()[k]
    
    def __setitem__(self,i,x):
        """ set the ith parameter to x, and, if different,
            set the changed property for the source"""
        source, k = self.index[:,i]
        model = source.model
        pars = model.get_parameters()
        if x==pars[k]: return
        pars[k] = x
        source.changed=True
        model.set_parameters(pars)
    
    def setitems(self, set_dict):
        """ set a set of items by index: the dict has keys that are either the index, or the name of the variabe, and float values,
            e.g. {1:1e-14, 2:2.1, 'Source_Index': 2.0}
        """
        def par_index(self, i):
            npar = self.__len__()
            if isinstance(i,int):
                if i<0 or i>=npar:
                    raise Exception('Index, %d, out of range for %d parameters' % (i,npar ) )
                return i
            else:
                try:
                    return list(self.parameter_names).index(i)
                except:
                    raise Exception('Parameter name "%s" not found' % i)
        for key,value in set_dict.items():
            i = par_index(self,key)
            print key, i, self[i], '-->', value
            self[i]=value
            
    def __len__(self):
        return self.index.shape[1]
        
    def get_parameters(self):
        """ return array of all parameters"""
        t = [s.model.get_parameters() for s in self.free_sources]
        return np.concatenate(t) if len(t)>0 else []
        
    def set_parameters(self, pars):
        """ set parameters, checking to see if changed"""
        i =0
        for source, n in self.ms:
            model = source.model
            oldpars = model.get_parameters()
            newpars = pars[i:i+n]
            if np.any(oldpars != newpars):
                source.model.set_parameters(newpars)
                source.changed=True
            i += n
    
    def get_covariance(self, nomask=False):
        """ get the covariance matrix from the souurce models
        """
        na,nt =len(self.mask), sum(self.mask)
        cov = np.matrix( np.zeros(na*na).reshape(na,na))
        i = 0
        for source, n in self.ms:
            model = source.model
            mcov = model.internal_cov_matrix[np.outer(model.free,model.free)]
            cov[i:i+n,i:i+n] = mcov.reshape(n,n)
            i += n
        if nomask: return cov
        return np.matrix(cov[np.outer(self.mask, self.mask)].reshape(nt,nt))
    
    def set_covariance(self, cov):
        """ save the specified convariance matrix into the source models"""
        cnow = np.asarray(self.get_covariance(nomask=True)).flatten()
        
        cnow[np.outer(self.mask, self.mask).flatten()] = np.array(cov).flatten()
        na = len(self.mask)
        cnew = cnow.reshape(na,na)
        i = 0
        for source, n in self.ms:
            model = source.model
            model.set_cov_matrix(cnew[i:i+n, i:i+n])
            i += n
    
    @property
    def model_parameters(self):
        if len(self.free_sources)==0: return []
        return np.concatenate([s.model.free_parameters for s in self.free_sources])
    
    @property
    def uncertainties(self):
        """ return relative uncertainties from diagonals of individual covariance matrices 
        """
        variances = np.concatenate([s.model.get_cov_matrix().diagonal()[s.model.free] \
            for s in self.free_sources])[self.mask]
        variances[variances<0]=0
        return np.sqrt(variances) / (np.abs(self.model_parameters) +1e-20) #avoid divide by zero

    @property 
    def bounds(self):
        """ fitter representation of applied bounds """
        return np.concatenate([source.model.bounds[source.model.free] for source in self.free_sources])

    def __repr__(self):
        return '%d parameters from %d free sources' % (len(self), len(self.free_sources))
    def clear_changed(self):
        for s in self.free_sources:
            s.changed=False
    @property
    def dirty(self):
        return np.array([s.changed for s in self.free_sources])
    @property
    def parameter_names(self):
        """ array of free parameter names """
        names = []
        for source in self.free_sources:
            for pname in np.array(source.model.param_names)[source.model.free]:
                names.append(source.name.strip()+'_'+pname)                
        return np.array(names)
        
    def parameter_summary(self, out=None):
        """formatted summary of parameter names, values, gradient
        out : None or open stream
        """
        if len(self.parameter_names)==0:
            print 'No free parameters'
            return
        print >>out,'\n%-21s %8s %8s' % ('parameter', 'value', 'error(%)')
        print >>out,  '%-21s %8s %8s' % ('---------', '-----', '--------')
        for u in zip(self.parameter_names, self.get_parameters(), self.uncertainties):
            print >>out, '%-21s %8.2f %8.1f' % u


class ParSubSet(ParameterSet):
    """ adapt ParameterSet to implement a subset
    to use, set the mask property to an array of bool or call the select function 
    """
    def __init__(self, roimodel, select=None, exclude=None, mask=None):
        """
        roimodel : ROImodel object
        mask    : [array of bool | None ]
        """
        self.roimodel=roimodel
        super(ParSubSet,self).__init__(roimodel)
        self.set_mask(mask)
        self.selection_description = None
        if select is not None:
            self.select(select, exclude)
        
    def __repr__(self):
        return '%s.%s: subset of %d parameters' % (self.__module__, self.__class__.__name__, sum(self.mask))
    def set_mask(self, m=None):
        if m is None:
            self._mask = np.ones(len(self),bool)
        else: 
            assert len(m)==len(self)
            assert sum(m)>0
            self._mask = m
        self.subsetindex = np.arange(len(self))[self._mask]
    def get_mask(self): return self._mask        
    mask = property(get_mask, set_mask) 
    
    def select(self, select=None, exclude=None):
        """
        Parameters
        ----------
        select : None, item or list of items, where item is an int or a string
            if not None, it defines a subset of the parameter numbers to select
                    to define a projected function to fit
            int:  select the corresponding parameter number
            string: select parameters according to matching rules
                    The name of a source (with possible wild cards) to select for fitting
                    If initial character is '_', match the rest with parameter names
                    if initial character is not '_' and last character is '*', treat as wild card
            
        exclude : None, int, or list of int 
                if specified, will remove parameter numbers from selection
        """

        # select a list of parameter numbers, or None for all free parameters
        selected= set()
        npars = len(self)
        
        if select is not None:
            try:
                selectpar = select
                if not hasattr(select, '__iter__'): select = [select]
                for item in select:
                    if type(item)==types.IntType or type(item)==np.int64:
                        selected.add(item)
                        if item>=npars:
                            raise Exception('Selected parameter number, %d, not in range [0,%d)' %(item, npars))
                    elif type(item)==types.StringType or type(item)==np.string_:
                        if item.startswith('_'):
                            # look for parameters
                            if item[-1] != '*':
                                toadd = filter( lambda i: self.parameter_names[i].endswith(item), range(npars) )
                            else:
                                def filt(i):
                                    return self.parameter_names[i].find(item[:-1])!=-1
                                toadd = filter( filt, range(npars) )
                        elif item in self.parameter_names:
                            toadd = [list(self.parameter_names).index(item)]
                            self.selection_description = 'parameter %s' % item
                        else:
                            src = self.roimodel.find_source(item)
                            self.selection_description = 'source %s'%src.name
                            toadd = filter(lambda i: self.parameter_names[i].startswith(src.name), range(npars))
                        selected = selected.union(toadd )
                    else:
                        raise Exception('fit parameter select list item %s, type %s, must be either an integer or a string' %(item, type(item)))
            except Exception, msg:
                raise Exception('Fail parameter select: {}'.format(msg))
            select = sorted(list(selected))
            if len(select)==0:
                raise Exception('nothing selected using "%s"' % selectpar)
        
        if exclude is not None:
            if not hasattr(exclude, '__iter__'): exclude = [exclude]
            all = set(range(npars)) if select is None else set(select)
            select = list( all.difference(exclude))
        t = np.zeros(len(self), bool)
        t[select]=True
        self.set_mask( t )
        if self.selection_description is None:
            self.selection_description = 'parameters %s' % select
        
    def __getitem__(self, i):
        """ access the ith parameter"""
        return super(ParSubSet,self).__getitem__(self.subsetindex[i])
    
    def __setitem__(self,i,x):
        super(ParSubSet,self).__setitem__(self.subsetindex[i], x)
    def get_parameters(self):
        t = super(ParSubSet, self).get_parameters()
        return t[self._mask]
    def set_parameters(self, pars):
        t = super(ParSubSet, self).get_parameters()
        t[self._mask]=pars
        super(ParSubSet, self).set_parameters(t)    
    
    @property
    def spectral_model(self):
        """ access to the spectral model for the first parameter"""
        return self.index[0, self.subsetindex[0]].model
    @property
    def source(self):
        """ access to the source associated with the first parameter"""
        return self.index[0, self.subsetindex[0]]
        
    def get_model(self,i):
        """spectral model for parameter i"""
        return self.index[0, self.subsetindex[i]].model
        
    @property
    def parameter_names(self):
        t = super(ParSubSet, self).parameter_names
        return t[self._mask]
    @property
    def bounds(self):
        t = super(ParSubSet, self).bounds
        return t[self.mask]
    @property
    def model_parameters(self):
        t = super(ParSubSet, self).model_parameters
        return t[self.mask]

