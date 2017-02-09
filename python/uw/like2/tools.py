"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.21 2016/10/28 21:24:32 burnett Exp $

"""
import os, sys, time
import numpy as np
import pandas as pd
from astropy.coordinates import Angle
import astropy.units as u

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

class OutputTee(WithMixin):
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
    def restore(self): self.close() #for with
        
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
    
      
class DateStamp(object):
    def _repr_html_(self):
         return '<div align=center> <h3>'+time.asctime()+'</h3> </div>'

class Profile(object):
    """ Create a "profile" of y, binned in x, according to bins
    """
    def __init__(self,x, y, bins):
        df = pd.DataFrame({'x': x, 'y': y})
        df['bin'] = np.digitize(x, bins=bins)
         # grouby bin, so we can calculate stuff
        binned = df.groupby('bin')
        # calculate mean and standard error of the mean for y in each bin
        self.pdf= binned['y'].agg(['mean', 'sem'])
        self.pdf['x'] = 0.5 * (bins[:-1] + bins[1:])
        self.pdf['xerr'] = (bins[1] - bins[0]) / 2

    def plot(self, **kwargs):
        """ Generate a plot using the profile DataFrame
        return the Axes object
        """
        return self.pdf.plot(
            x='x',       y='mean',
            xerr='xerr', yerr='sem',
            linestyle='none', capsize=0,  color='black',
            **kwargs
        )

def html_table(dataframe, float_precision=2):
    """Return an HTML table of a DataFrame with given precision
    """
    from IPython.display import display, HTML
    fmt = '{:.'+ '{:d}'.format(float_precision) + 'f}'
    return display(HTML(dataframe.to_html(float_format=lambda x:fmt.format(x))))

def parse_jname(name):
    """ return an (ra,dec) tuple from a name in J format
    J0258.9+0552 -> (44.725, 5.86667)
    """
    ra=(name[1:3]+'h'+name[3:7]+'m')
    dec = (name[7:10]+'d'+name[10:12]+'m')
    return map(lambda a: float(Angle(a, unit=u.deg).to_string(decimal=True)),(ra,dec))

def create_jname(ra,dec):
    """Create the format Jhhmm.m+ddmm
    note that last digit is truncated, not rounded
    """
    ram = np.mod(int(ra*40),1440)/10. # RA in minutes, truncate at 0.1 
    sign= '+' if dec>=0 else '-'
    dem = int(abs(dec)*60) #abs( DEC) in minutes, truncated
    return 'J' +   '{:02d}{:04.1f}'.format(int(ram/60),ram%60)\
            + sign+'{:02d}{:02.0f}'.format(int(dem/60),dem%60)
 