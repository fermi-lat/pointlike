"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.14 2013/11/23 16:05:33 burnett Exp $

"""
import os
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

    
      
