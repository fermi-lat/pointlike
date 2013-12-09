"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.17 2013/12/04 05:25:51 burnett Exp $

"""
import os, sys
import numpy as np

def ufunc_decorator(f): # this adapts a bound function
    def new_ufunc(self, par):
        return np.array(map(lambda x: f(self,x),par)) if hasattr(par, '__iter__')  else f(self,par)
    return new_ufunc

# special function to replace or extend a docstring from that of another function
def decorate_with(other_func, append=False, append_init=False):
    """ append_init: If decorating with an object (which has an __init__ function),
                     then append the doc for the __init__ after the doc for the
                     overall class. """
    def decorator(func):
        if append: func.__doc__ += other_func.__doc__ 
        else:      func.__doc__  = other_func.__doc__ 

        if append_init and hasattr(other_func,'__init__'):
                func.__doc__ += other_func.__init__.__doc__
        return func
    return decorator


    
class WithMixin(object):
    """Mixin to allow simple restore of an object's state
        supports the 'with' construction, guarantees that restore is called to restore the state of the model
        example:
        -------
        with ClassName(...) as something:
            # use something ...
    """
    def __enter__(self):
        return self
        
    def __exit__(self, type, value, traceback):
        self.restore()

class OutputTee(object):
    def __init__(self, logfile):
        self.logstream = open(logfile, 'a')
        self.stdout = sys.stdout
        sys.stdout = self
    def write(self, stuff):
        self.logstream.write(stuff)
        self.stdout.write(stuff)
    def close(self):
        sys.stdout =self.stdout
        self.logstream.close()
    def flush(self):
        self.stdout.flush()
    def set_parent(self, parent):
        self.stdout.set_parent(parent) #needed??
        
class RecArray(object):
    def __init__(self, names, dtype=None):
        self.names = names
        self.fields=list([list() for i in range(len(names))])
        self.dtype=dtype

    def append(self,*arg):
        for i,x in enumerate(arg):
            self.fields[i].append(x)
    def __call__(self):
        """ return finished recarray"""
        return np.rec.fromarrays(self.fields, names=self.names, dtype=self.dtype)
    
      
