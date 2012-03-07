"""
Support for generating doc strings, and setting keyword options for class constructors
  decorate: decorator to append keyword info to the docstring
  process:  set the class dictionary from the defaults and supplied keywords
  
$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/utilities/keyword_options.py,v 1.9 2011/06/17 03:57:09 lande Exp $

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
            if len(item)==3:
                key, value, description = item 
            else:
                (key, value), description = item, ''
            if type(value)==types.StringType:
                value = "'" + value + "'"
            s += indent+'%-12s' % key
            if len(key)>=12: s += indent + 12*' '
            s += '%-10s' % str(value)
            if len(str(value))>10: s += indent + 23*' '
            s += ' '+ (indent+23*' ').join(description.split('\n'))
        if func.__doc__ is None: func.__doc__ = ''
        func.__doc__ += s+hbar
        return func
    return decorator

# use with the above
def process(self, kwargs, defaults=None):
    """
    self: class instance, used to set the dictionary, and find the name
    kwargs: kwargs entry in function
    defaults: None, list or tuple
        if None, use the 'default' class member
        a list item is either 
            * string, which is ignored
            * a list or 3-tuple (name, value, comment)
    Raises KeyError exception for any kwargs entry not in the defaults list
    """
    for item in self.defaults if defaults is None else defaults:
        if type(item)==types.StringType: continue
        self.__dict__[item[0].strip()] = item[1]
        
    for key in kwargs.keys():
        if key in self.__dict__: self.__dict__[key]=kwargs[key]
        else:
            raise KeyError, "option '%s' not recognized by %s" % (key,self.__class__.__name__)


def defaults_to_kwargs(obj,default_object):
    """ Take in a defaults list (used by keyword_options) and an object 
        which recognizes the keyword_options. A dictionary is
        returned with each of the defaults pointing to the value
        found in the object. This is useful for recreating 
        the object. """
    return dict([[i[0],getattr(obj,i[0])] for i in \
                default_object.defaults if isinstance(i,list) or isinstance(i,tuple)])

def change_defaults(defaults,key,value):
    """ Change a defaults dictionary's key 'key'
        to have a default value 'value'. """
    if isinstance(defaults,tuple):
        defaults=list(defaults)
        
    for i,default in enumerate(defaults):
        if default[0] == key:
            if isinstance(default,tuple): 
                defaults[i]=list(defaults[i])
                default=defaults[i]
            default[1] = value
            return defaults
    raise Exception("key %s not found in defaults" % key)

def get_default(*args, **kwargs): return get_row(*args, **kwargs)[1]

def get_row(defaults, key):
    for i,default in enumerate(defaults):
        if default[0] == key: return default
    raise Exception("key %s not found in defaults" % key)
