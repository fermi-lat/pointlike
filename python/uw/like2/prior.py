"""
Manage prior likelihood calculations 


$Header$
Author: T.Burnett <tburnett@uw.edu>
"""
import numpy as np

        
class NoPrior(object):
    def __call__(self,x): return 0
    def grad(self,x): return 0
    def delta(self,x): return 0
    def apply_limit(self,x): return x
    def __str__(self): return 'NoPrior'
    
class LimitPrior(NoPrior):
    def __init__(self,a,b,c=1e2):
        self.a,self.b,self.c=a,b,c
    def __str__(self):
        return 'LimitPrior(%f,%f)' % (self.a,self.b)
    def apply_limit(self,x):
        if   x<self.a: return self.a
        elif x>self.b: return self.b
        return x
    def __call__(self,x):
        if x>self.b:   return self.c*(x-self.b)**2
        elif x<self.a: return self.c*(x-self.a)**2
        return 0
    def grad(self,x):
        if   x>self.b: return 2*self.c*(x-self.b)
        elif x<self.a: return 2*self.c*(x-self.a)
        return 0
    def delta(self,x):
        if   x>self.b: return x-self.b
        elif x<self.a: return x-self.a
        return 0
           
class ModelPrior(object):
    """ manage addition of priors to the likelihood"""

    def __init__(self, sources, free=None):   
        self.sources = sources
        self.models = np.array([s.spectral_model for s in sources])
        self.initialize(free)
        self.update()
        self.enabled=True

    def initialize(self, free):
        self.free = free if free is not None else self.sources.free
        self.free_models = self.models[self.free]
        self.npar = np.sum([np.sum(m.free) for m in self.free_models])
        self.priors= []
        for model in self.free_models:
            for pname in np.array(model.param_names)[model.free]:
                if pname=='Index':
                   prior = LimitPrior(np.log10(1e-3), np.log10(3.0), 10.)
                elif pname=='Norm':
                    if model.name =='PowerLaw': #for ring
                        prior = LimitPrior(np.log10(0.5), np.log10(1.5), 10.)
                    else: prior = LimitPrior(-16, -8, 1e4)
                elif pname=='Scale': # used by iso
                   prior = LimitPrior(np.log10(1e-4), np.log10(1.5), 10.)
                else:
                    prior = NoPrior()
                self.priors.append(prior)
  
    def update(self, reset=False):
        self.pars = self.sources.get_parameters()

    def limit_pars(self, update=False):
        self.update()
        for i, prior in enumerate(self.priors):
            self.pars[i] = prior.apply_limit(self.pars[i])
        if update:
            self.sources.set_parameters(self.pars)
        
    def log_like(self):
        if not self.enabled: return 0
        t = np.sum(prior(x) for prior,x in zip(self.priors, self.pars))
        if False: # t>10.0: 
            n = np.arange(len(self.pars))[np.abs(self.gradient())>0]
            print 'likelihood penalty: %.1f %s' % (t,n)
        return -t
    
    def check(self):
        names = np.array(self.sources.parameter_names)
        deltas = np.array([prior.delta(x) for prior,x in zip(self.priors, self.pars)])
        sel= deltas!=0
        return np.vstack([names[sel], self.pars[sel], deltas[sel]])
        
    def gradient(self):
        if not self.enabled: return np.zeros(len(self.pars))
        return np.array([prior.grad(x) for prior,x in zip(self.priors, self.pars)])