"""
Support for generating doc strings, and setting keyword options for class constructors
  decorate: decorator to append keyword info to the docstring
  process:  set the class dictionary from the defaults and supplied keywords
  
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/keyword_options.py,v 1.3 2010/09/14 07:49:24 lande Exp $

Author: T. Burnett <tburnett@uw.edu>
"""
import types 
import textwrap

##
#
def decorate(defaults):
    """
     Special decorator for the __init__ function of a class to append a description
     of valid keyword options.
     Assumes the class has a static variable defaults, which is a list of:
     * key, default_value, description
       where:
         key: a string that is a valid name
         default_value: any value
         description: a string describing the key: it can be broken up with newlines
    
     * a string to be used for the docstring, as a label to group the keywords perhaps
     
     usage: 
    class Myclass(object):
        defaults = ( ('key', value, 'description'), )
        @keyword_options.decorate(defaults)
        def __init__(self, *pars, **kwargs):
            keyword_options.process(self, kwargs)
    """
    indent = '\n\t\t' # do not know why two tabs seem to be needed
    hbar   = indent+60*'=' # horizontal bar
    def decorator(func):
        s= hbar+ indent+'keyword arguments'+ hbar
        for item in defaults:
            if type(item)==types.StringType:
                s+= '\n%s%s   %s'% (indent,9*'=',item.upper())
                continue

            key, value, description = item    
            if type(value)==types.StringType:
                value = "'" + value + "'"
            s += indent+'%-12s' % key
            if len(key)>=12: s += indent + 12*' '
            s += '%-10s' % value
            if len(str(value))>10: s += indent + 23*' '
            s += ' '+ (indent+23*' ').join(description.split('\n'))
        if func.__doc__ is None: func.__doc__ = ''
        func.__doc__ += s+hbar
        return func
    return decorator

# use with the above
def process(self, kwargs):
    """
    self: class instance, used to set the dictionary, and find the name
    kwargs: kwargs entry in function
    note assumes the defaults list is named 'defaults'
    
    Raises KeyError exception for any kwargs entry not in the defaults list
    """
    for item in self.defaults:
        if type(item)==types.StringType: continue
        self.__dict__[item[0].strip()] = item[1]
    for key in kwargs.keys():
        if key in self.__dict__: self.__dict__[key]=kwargs[key]
        else:
            raise KeyError, "option '%s' not recognized by %s" % (key,self.__class__.__name__)

