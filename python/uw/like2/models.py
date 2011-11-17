"""
This is a place for modified or new spectral models needed for like2
$Header$
"""
import operator
import numpy as np
from uw.like import Models

# This is mostly copied from like.Models, but needed some changes

class CompositeModel(Models.Model):
    """ A model which joins other models. Subclasses must
        implement the __call__(e) function which says
        how to join the models. """

    def __init__(self,*models, **kwargs):
        iscopy = kwargs.pop('iscopy', False)
        if len(models) < 1:
            raise Exception("CompositeModel must be created with more than one spectral model")
        for m in models:
            if not isinstance(m,Models.Model):
                raise Exception("CompositeModel must be created with a list of models.")

        self.flux_scale = 1.
        self.models = models
        self.cov_matrix = np.zeros([self.npar,self.npar]) #default covariance matrix
        self.e0 = 1000. # not sure why, but sed_plotter needs this

    def __call__(self,e): 
        """ Must be implemented by a subclass. """
        pass

    @property
    def param_names(self):
        return reduce(operator.add,[i.param_names for i in self.models])

    @property
    def pretty_name(self):
        return self.operator.join([i.pretty_name for i in self.models])

    @property
    def background(self):
        """ Seems like a reasonable test. """
        return np.any([i.background for i in self.models])

    @property
    def n(self): 
        return np.asarray([len(i._p) for i in self.models])

    @property
    def npar(self): 
        return sum([len(i._p) for i in self.models])

    @property
    def _p(self): 
        return np.append(*[i._p for i in self.models])

    @_p.setter
    def _p(self, value):
        assert(len(self._p) == len(value))
        counter=0
        for i in xrange(len(self.n)):
            self.models[i]._p = value[counter:counter+self.n[i]]
            counter += self.n[i]

    @property
    def free(self): 
        return np.append(*[i.free for i in self.models])

    @free.setter
    def free(self, value):
        assert(len(self.free) == len(value))
        counter=0
        for i in xrange(len(self.n)):
            self.models[i].free = value[counter:counter+self.n[i]]
            counter += self.n[i]
      
    def setp(self, i, par, internal=False):
        """ set internal value, convert unless internal
        """
        i=self.mapper(i) #?
        if not internal: 
            assert par>0, 'Model external parameter cannont be negative'
        counter=0
        for k in xrange(len(self.n)):
            if counter<i:
                counter += self.n[k]
                continue
            self.models[k]._p[i-counter] = par if internal else  np.log10(par)
            return       

    def set_parameters(self,new_vals):
        """Set FREE internal parameters; new_vals should have length equal to number of free parameters."""
        assert len(new_vals)==(self.free).sum(), 'attempt to set wrong number of free parameters, model %s' %self
        pars = self._p
        pars[self.free] = new_vals.astype('f')
        self._p = pars
 

class FrontBackConstant(CompositeModel):
    """ Composite model that is either/or, for front or back
        select which constant based on value (0 or 1) of ct
    """
    name = 'FrontBackConstant'
    operator='+'
    def __init__(self, f=1, b=1, **kwargs):
        super(FrontBackConstant, self).__init__(Models.Constant(),Models.Constant(), **kwargs)
        self.models[0].param_names=['Scale_front']
        self.models[1].param_names=['Scale_back']
        self.models[0][0]=f
        self.models[1][0]=b
        self.ct = 0
        
    def __call__(self, e):
        return self.models[self.ct](e)
        
    def gradient(self, e):
        return np.hstack([(1-self.ct)*self.models[0].gradient(e).T, 
                          self.ct*self.models[1].gradient(e).T])
                          
                          