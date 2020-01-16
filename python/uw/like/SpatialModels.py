"""A set of classes to implement spatial models.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/SpatialModels.py,v 1.114 2018/01/27 15:36:15 burnett Exp $

   author: Joshua Lande

"""
import os
import copy
import numbers
from abc import abstractmethod

import numpy as np

from scipy import vectorize
from scipy.interpolate import interp1d, griddata, UnivariateSpline, SmoothBivariateSpline, LinearNDInterpolator
from scipy.integrate import quad
from scipy.optimize import fmin

from skymaps import PySkySpectrum,PySkyFunction,SkyDir,Hep3Vector,\
        SkyImage,SkyIntegrator,CompositeSkyFunction,PythonUtilities

from uw.utilities.quantile import Quantile
from uw.utilities.rotations import anti_rotate_equator
from uw.utilities import path


SMALL_ANALYTIC_EXTENSION=1e-10
SMALL_NUMERIC_EXTENSION=1e-3

def smart_log(p,log):
    """ This function is functionally equivalent to

            np.where(log, np.log10(p), p)

        but will not raise errors by taking the log of
        negative parameters when log=False
    """
    return np.where(log, np.log10(np.where(log, p, 1)), p)

class SpatialQuantile(object):
    """ Numerically compute the radial quantile for a spatial model, 
        assuming that space is flat. 
        
        It is easy to test this code against the values computed
        analytically for the disk and gauss.

            >>> tol = dict(rtol=3e-2, atol=3e-2)

            >>> disk = Disk(sigma=1, center=SkyDir(0, 0))
            >>> gauss = Gaussian(sigma=0.5, center=SkyDir(0,0, SkyDir.GALACTIC))
            >>> for i in [disk, gauss]:
            ...     np.allclose(i.analytic_r68(),  i.numeric_r68(), **tol)
            ...     np.allclose(i.analytic_r99(),  i.numeric_r99(), **tol)
            True
            True
            True
            True

        Good to make sure the less-analytically defined extended sources are also correct:

            >>> elliptical_gauss = EllipticalGaussian(major_axis=0.5, minor_axis=0.5, center=SkyDir(0,0, SkyDir.GALACTIC))
            >>> np.allclose(gauss.analytic_r68(),  elliptical_gauss.numeric_r68(), **tol)
            True
            >>> np.allclose(gauss.analytic_r99(),  elliptical_gauss.numeric_r99(), **tol)
            True
            >>> a,b=gauss.analytic_r68(),  elliptical_gauss.numeric_r68()
    """

    def __init__(self, spatial_model, quad_kwargs=dict()):
        self.spatial_model = spatial_model
        self.center = self.spatial_model.center

        # These tolerances get a 1% error in unittests.
        # full_output surpresses errors.
        self.quad_kwargs = dict(epsabs=np.inf, epsrel=1e-40, limit=200)
        self.quad_kwargs.update(quad_kwargs)

        rmax = self.spatial_model.effective_edge()
        self.quantile = Quantile(self.integrand, 0, rmax,
                                 quad_kwargs=self.quad_kwargs)

    def pdf(self, r,theta):
        """ Evaluate the the """
        x = r*np.cos(theta)
        y = r*np.sin(theta)

        sd=anti_rotate_equator(SkyDir(x,y),self.center)
        return self.spatial_model(sd)

    def integrand(self, r):
        return quad(lambda theta: self.pdf(r, theta), 0, 2*np.pi, **self.quad_kwargs)[0]*r

    def r68(self): return self(0.68)
    def r99(self): return self(0.99)

    def __call__(self, *args, **kwargs): return self.quantile(*args, **kwargs)



class SpatialModelException(Exception): 
    pass

class SpatialModel(object):
    """ This class represents a normalized spatial model which 
        contains a list of spatial parameters

        Any subclass must implement the __call__ function.
        __call__ takes a SkyDir 
        and returns the intensity at that
        direction. The intensity is defined such that the integral
        of intensity x solid angle = 1 when solid angle is measured
        in steradians!  Each subclass is responsible for ensuring its correct normalization.

        Any subclass must also define a list of spatial parameters by setting the following values:
    
            default_p: the values of the spatial parameters. 
                There are all assumed to be absolute.
            param_names: the names of each spatial parameters
            default_limits: the limits imposed on the parameters when fitting. 
                These values are absolute. Note that by and large, everything
                is measured in degrees. The limits on relative movement of
                longitude and latitude are measured in degrees.
            log: For each paramter, whether or not it should be mapped into log space.
            steps: Used by minuit when fitting the source. useful to make them comparable
                to distance away from true values. These are not absolute. For log parameters,
                the step is the step in the log of the parameter.

        By construction, the first two spatial parameters of all extended sources
        is the center of the source. These paramters will be set automatically
        and so do not need to be set in any subclass. For spatial models
        which define only the position in the sky, set the above paramters
        to empty lists.
        
        One slight difference between the SpatialModel and Model class
        has to do with absolute. Always, absolute=True has the same
        meanings, pass in the true value of the parameters. But
        for some functions, absolute=False has a different meanings.
        For most functions, absolute=False means the values for which log=True
        should be in log space. This is the calse for set_parameters, 
        get_parameters, and get_cov_matrix. This is different from
        how Models works where absolute=False means all parameters are
        in log space. On the other hand, the function statistical has absolute=False
        interpreted the same convention as absolute, where absolute=False
        means return the relative error (absolute error/parameter value)
        which is useful for printing percent error, etc. 
        Generally, the user should ignore these flags since it is
        easier to use the array slicing operators to get always-absolute
        paramters e.g. disk['sigma'] = 2.
       

        The center of the extended source can be set in several ways:

        By default, the center is the equatorial origin:

            >>> disk = Disk()
            >>> print (np.allclose([disk.center.ra(),disk.center.dec()],[0,0]))
            True
            >>> print (np.allclose([disk['ra'],disk['dec']],[0,0]))
            True

        But this can be overridden

            >>> disk = Disk(coordsystem=SkyDir.GALACTIC)
            >>> print (np.allclose([disk.center.l(),disk.center.b()],[0,0]))
            True
            >>> print (np.allclose([disk['l'],disk['b']],[0,0]))
            True

        The source center can be specified by defining the galactic origin

            >>> disk = Disk(l=22, b=22)
            >>> print (np.allclose([disk.center.l(),disk.center.b()],[22,22]))
            True
            >>> print (disk.coordsystem == SkyDir.GALACTIC)
            True

        Or the equatorial origin

            >>> disk = Disk(ra=22, dec=22)
            >>> print (np.allclose([disk.center.ra(),disk.center.dec()],[22,22]))
            True
            >>> print (disk.coordsystem == SkyDir.EQUATORIAL)
            True

        Or by specifying the center of the source

            >>> disk = Disk(center=SkyDir(22, 22))
            >>> print (np.allclose([disk.center.ra(),disk.center.dec()],[22,22]))
            True
            >>> print (disk.coordsystem == SkyDir.EQUATORIAL)
            True

        The names of the input paramters should be case insensitive:

            >>> disk = Disk(rA=10, dEC=-10)
            >>> print (np.allclose([disk.center.ra(),disk.center.dec()],[10,-10]))
            True

            >>> disk = Disk(L=33, B=10)
            >>> print (np.allclose([disk.center.l(),disk.center.b()],[33,10]))
            True

        Ambiguous specifications should raise an exception:

            >>> disk = Disk(center=SkyDir(), l=22, b=22)
            Traceback (most recent call last):
                ...
            Exception: Only one of center, l and b, or RA and DEC can be set

        Testing setting the free:

            >>> disk = Disk(free=[True,False,False])
            >>> print (disk.free.tolist())
            [True, False, False]

        Test setting limits:

            >>> disk = Disk(limits=[[-.1, .1], [-.1, .1], [1e-20,2]])
            >>> disk.get_limits(absolute=True).tolist()
            [[-0.1, 0.1], [-0.1, 0.1], [1e-20, 2.0]]

    """
    def __init__(self, p=None, coordsystem=SkyDir.EQUATORIAL, free=None, limits=None, **kwargs):

        assert hasattr(self,'default_p') and \
                hasattr(self,'param_names') and \
                hasattr(self,'default_limits') and \
                hasattr(self,'log') and \
                hasattr(self,'steps')

        self.name = self.pretty_name = self.__class__.__name__

        # first, set center if specified

        lower_case_kwargs = [i.lower() for i in kwargs.keys()]
        center_in_kwargs = 'center' in kwargs
        l_b_in_kwargs = 'l' in lower_case_kwargs and 'b' in lower_case_kwargs
        ra_dec_in_kwargs = 'ra' in lower_case_kwargs and 'dec' in lower_case_kwargs

        if sum([center_in_kwargs, l_b_in_kwargs, ra_dec_in_kwargs]) > 1:
            raise Exception("Only one of center, l and b, or RA and DEC can be set")

        if center_in_kwargs:
            self.center = kwargs.pop('center')
            self.coordsystem = coordsystem
        elif l_b_in_kwargs:
            # acutal center will be set later in set_parameters
            # here, we just need to get a defualt center and set 
            # the right coord system. Kludge for now - J.L.
            self.center = SkyDir() 
            self.coordsystem = SkyDir.GALACTIC
        elif ra_dec_in_kwargs:
            # same disclaimer as above - J.L.
            self.center = SkyDir() 
            self.coordsystem = SkyDir.EQUATORIAL
        else:
            self.center = SkyDir(0,0,coordsystem)
            self.coordsystem = coordsystem

        if self.coordsystem == SkyDir.EQUATORIAL:
            self.param_names = np.append(['RA','DEC'], self.param_names)
        elif self.coordsystem == SkyDir.GALACTIC:
            self.param_names = np.append(['L','B'], self.param_names)

        # The first two parameters (lon & lat) are forced to have log=False
        # cast to bool, just to be safe
        self.log=np.append([False,False],self.log).astype(bool)

        # note, self.log and self.param_names must be 
        # correctly initialized before set_parameters is called.
        self.set_parameters(p if p is not None else self.default_p, absolute=True)

        # The limits on the first two parameters are
        # defined as the allowed physical angular distance away from the source.
        if limits is not None:
            self.limits = np.asarray(limits)
        else:
            if self.default_limits == []:
                self.limits=np.asarray([[-1.,1.],[-1.,1.]])
            else:
                self.limits=np.append([[-1.,1.],[-1.,1.]],self.default_limits,axis=0)

        self.steps=np.append([0.1,0.1],self.steps)

        self.cov_matrix = np.zeros([len(self.p),len(self.p)])
        # compatibility with older numpy
        if free is not None:
            self.free = np.asarray(free)
        else:
            self.free = np.ones_like(self.p).astype(bool)

        # map the parameters/limits into log space.
        for i in range(2,len(self.log)):
            self.limits[i,:] = smart_log(self.limits[i,:], log=self.log[i])

        for k in kwargs.keys():
            if k in self:
                self[k] = kwargs.pop(k)

        if len(kwargs) > 0:
            raise Exception("Unrecognized keyword arguments %s" % kwargs.keys())

        self.cache() # overloaded by subclasses

    def cache(self):
        """ This should be inhereted by child classes to cache
            various useful quanitites after update is called. """

        self.center = SkyDir(self.p[0],self.p[1],self.coordsystem)

    def change_coordsystem(self,cs):
        """ Change the internal coordinate system. This is what is
            used when the source is displayed/what is read in
            as longitude and latitude when a parameter is set. Also
            changes what the errors are estimates of. """
        self.coordsystem = cs
        if cs  == SkyDir.EQUATORIAL:
            self.param_names[0:2] = ['RA','DEC']
            self.p[0:2] = [self.center.ra(),self.center.dec()]
        elif cs == SkyDir.GALACTIC:
            self.param_names[0:2] = ['L','B']
            self.p[0:2] = [self.center.l(),self.center.b()]

        # Errors are no longer valid, so reset cov matrix.
        self.cov_matrix = np.zeros([len(self.p),len(self.p)]) 

    def len(self):
        return len(self.p)

    def __contains__(self, item):
        try:
            self.__getitem__(item)
            return True
        except Exception as ex:
            return False
        
    def __getitem__(self, index):
        return self.getp(index)
        
    def __setitem__(self, index, value):
        self.setp(index,value)
        
    def getp(self, i, internal=False):
        """ get external value for parameter # i """
        i=self.mapper(i)
        return np.where((not internal) & self.log,10**self.p,self.p)[i]

    def numeric_r68(self, *args, **kwargs): return self.numeric_quantile(.68, *args, **kwargs)
    def numeric_r99(self, *args, **kwargs): return self.numeric_quantile(.99, *args, **kwargs)
    def numeric_quantile(self, quantile, *args, **kwargs): 
        return SpatialQuantile(self, *args, **kwargs)(quantile)


    def error(self,i, internal=False):
        """ get error for parameter # i """
        i=self.mapper(i)
        return (np.diag(self.get_cov_matrix(absolute=not internal))**0.5)[i]

    def setp(self, i, par, internal=False):
        """ set internal value, convert unless internal """
        i=self.mapper(i)
        if (not internal) and self.log[i] and par <=0:
            raise Exception("Parameter %s must be positive!" % i)

        self.p[i] = np.log10(par) if ((not internal) and self.log[i]) else par
        self.cache()

    def get_parameters(self,absolute=False):
        """Return all parameters; used for spatial fitting. 
           This is different from in Models.py """
#        return np.where(self.log,10**self.p,self.p) if absolute else self.p
#       THB: this is formally the same, but does not raise runtime overflow
        return np.array([10**x if y else x for x,y in zip(self.p,self.log)]) if absolute else self.p

    def get_param_names(self,absolute=True):
        if absolute:
            return self.param_names
        else:
            # Can't fit l or dec in log space.
            return ["log10(%s)" % n if log \
                    else n.replace('DEC','DEC (rotated)').replace('b','b (rotated)') \
                    for n,log in zip(self.param_names,self.log)]

    def set_limits(self, i, lower, upper, absolute=True):
        """ >>> sm = Disk()
            >>> print (sm.get_limits().tolist())
            [[-1.0, 1.0], [-1.0, 1.0], [1e-10, 3.0]]
            >>> sm.set_limits('dec', -2, 2)
            >>> print (sm.get_limits().tolist())
            [[-1.0, 1.0], [-2.0, 2.0], [1e-10, 3.0]]
            >>> sm.set_limits('sigma', 1e-20, 10)
            >>> print (sm.get_limits().tolist())
            [[-1.0, 1.0], [-2.0, 2.0], [1e-20, 10.0]]

        """
        i=self.mapper(i)
        #THB avoid log(0)
        if lower==0: lower=1e-6 
        lim = np.asarray([lower,upper],dtype=float)

        self.limits[i,:] = smart_log(lim, log=self.log[i])

    def get_limits(self,absolute=True):
        ret = np.asarray([10**lim if log and absolute else lim \
                         for lim,log in zip(self.limits,self.log)])
        return ret

    def get_steps(self):
        if not self.__dict__.has_key('steps'):
            raise Exception("Spatial model %s does not have fitting step sizes defined for it." % self.pretty_name)
        return self.steps

    def set_parameters(self,p,absolute=False,center=None):
        """ Set all parameters; p should have length equal to number of parameters.
            Note that this API is different from in Models.py

            If center is given as an argument, the longitude and latitude are appended to the beginning of the p
            as the first two coordinaets..
        """
        if center:
            if self.coordsystem == SkyDir.EQUATORIAL:
                p = np.append([center.ra(),center.dec()],p)
            elif self.coordsystem == SkyDir.GALACTIC:
                p = np.append([center.l(),center.b()],p)
        if isinstance(p,numbers.Real) and len(self.param_names)==3:
            return self.set_parameters([p],absolute,center=self.center)
        elif len(p)==len(self.param_names)-2:
            return self.set_parameters(p,absolute,center=self.center)

        if len(p)!=len(self.param_names):
            raise Exception("SpatialModel.set_parameters given the wrong number of parameters.")

        self.p = smart_log(p,log=self.log) if absolute else np.asarray(p,dtype=float)
        try:
            self.cache()
        except RuntimeWarning:
            print ('Overflow: parameters: before {}\n\t after {}'.format(p, self.p))

    def mapper(self,i):
        """ Maps a parameter to an index. 

                >>> d = Disk()
                >>> d.mapper('RA')
                0
                >>> d.mapper(['RA','DEC'])
                array([0, 1])

                >>> d = Disk(coordsystem=SkyDir.GALACTIC)
                >>> d.mapper('b')
                1

            """
        if isinstance(i,str):
            if i.lower() not in np.char.lower(self.param_names):
                raise Exception("Unknown parameter name %s" % i)
            return np.where(np.char.lower(self.param_names)==i.lower())[0][0]
        elif type(i) == int:
            return i
        elif type(i) == list:
            return np.asarray(map(self.mapper,i), dtype=int)
        else:
            raise Exception("Unknown parameter name %s" % i)
    
    def modify_loc(self,center):
        self.center = center
        if self.coordsystem == SkyDir.EQUATORIAL:
            self.p[0:2] = [center.ra(),center.dec()]
        elif self.coordsystem == SkyDir.GALACTIC:
            self.p[0:2] = [center.l(),center.b()]

        self.cache()

    def freeze_position(self,freeze=True):
        """ Freeze the source position. 

                >>> d = Disk()
                >>> print (d.free)
                [ True  True  True]
                >>> d.freeze_position()
                >>> print (d.free)
                [False False  True]
        """
        self.freeze([0,1],freeze)

    def freeze(self,i,freeze=True):
        """ Freeze one of the spatial parameters from fitting.
      
            i: a parameter name or index.
            freeze   : if True, freeze parameter; if False, free it 

            Previously, this function did not accept the case insensitive
            nature of input parameter names. So test both 'sigma' 
            and 'Sigma' for a disk source:

                >>> d = Disk()
                >>> print (d.free)
                [ True  True  True]
                >>> d.freeze('Sigma')
                >>> print (d.free)
                [ True  True False]
                >>> d.set_free('sigma')
                >>> print (d.free)
                [ True  True  True]
            """
        i=self.mapper(i)
        self.free[i] = not freeze
        self.cache()

    def set_free(self, i, free=True):
        self.freeze(i, freeze=not free)

    def get_free(self, i):
        """ Figure out of a parameter is free or not

                >>> d = Disk()
                >>> d.set_free('sigma', False)
                >>> print (d.get_free('sigma'))
                False
                >>> d.set_free('sigma', True)
                >>> print (d.get_free('sigma'))
                True
        """
        i=self.mapper(i)
        return self.free[i]

    def set_cov_matrix(self,new_cov_matrix):
        self.cov_matrix = new_cov_matrix

    def get_cov_matrix(self,absolute=True):
        """Return covariance matrix."""

        jac = np.log10(np.exp(1))
        p = np.where(self.log,(10**self.p)/jac,1) if absolute else np.ones_like(self.p)
        pt=p.reshape((p.shape[0],1)) #transpose
        return p*self.cov_matrix*pt

    def get_free_errors(self,absolute=False):
        """Return the diagonal elements of the covariance matrix for free parameters."""
        return np.diag(self.get_cov_matrix(absolute))**0.5

    def statistical(self,absolute=False,two_sided=False):
        """Return the parameter values and fractional statistical errors.
           If no error estimates are present, return 0 for the fractional error."""

        p = self.get_parameters(absolute=True)
        if not two_sided:
            # for one sided case, completely map covarinace matrix
            # to absolute values & then divide by p to get relative
            # errors
            errs = self.get_free_errors(absolute=True)
            return p,errs/(1. if absolute else p)
        else:
            # Perfrom conversion out of log space.
            errs = self.get_free_errors(absolute=False)
            lo_abs = np.where(self.log,p-10**(self.p-errs),errs)
            hi_abs = np.where(self.log,10**(self.p+errs)-p,errs)
            return  p, \
                    hi_abs/(1. if absolute else p), \
                    lo_abs/(1. if absolute else p)

    def copy(self): return copy.deepcopy(self)

    @abstractmethod
    def __call__(self,v,energy=None): pass

    @abstractmethod
    def has_edge(self): pass

    @abstractmethod
    def effective_edge(self,energy=None):
        """ It is useful to know an approximate edge to the image
            It is defined as a radius such that from the center to the
            enclosed circle contains (approximatly) the entire object. """
        pass

    def get_PySkyFunction(self):
        return PySkyFunction(self)

    def get_PySkySpectrum(self):
        """ Note that I am not setting the Integral function out of
            pure lazieness, since it is not needed elsewhere. """
        return PySkySpectrum(self,None)

    def save_template(self,filename,npix=150, proj='ZEA'):
        """ Saves out a template following the recommendation of
            the page LAT-Detected Extended Sources for Catalog:

                https://confluence.slac.stanford.edu/x/Qw2JBQ
        
            npix is the number of pixels in tempalate in either dimension. 
            
            
            When we save an analytic shape to a fits file and then load 
            it as a SpatialMap, we end up with a spatial model that
            is almost identical:

            First, create a uniform disk spatial model

                >>> center = SkyDir()
                >>> disk = Disk(sigma=1, center=center)

            Save it to a file:

                >>> from tempfile import NamedTemporaryFile
                >>> temp = NamedTemporaryFile()
                >>> filename = temp.name
                >>> disk.save_template(filename, npix=300)

            Then, create the spatial map:

                >>> spatial_map = SpatialMap(file = filename)
                
            In the center, the values are the same:

                >>> np.allclose(disk(center), spatial_map(center), atol=1e-3, rtol=1e-3)
                True

            Outside the disk edge, the spatial map has intensity 0:

                >>> print (spatial_map(SkyDir(1.5,0)))
                0.0
        """
        if not hasattr(self,'template_diameter'):
            raise Exception("Unable to save template because template_diameter is not defined.")

        center=self.center

        diameter=self.template_diameter()
        pixelsize=diameter/npix
        image=SkyImage(center,path.expand(filename),pixelsize,diameter,1,proj,
                       True if self.coordsystem == SkyDir.GALACTIC else False,False)
        skyfunction=self.get_PySkyFunction()
        image.fill(skyfunction)

        # Here, explicitly renormalize template to account for issues due to
        # coarse pixelization

        data = np.asarray(image.image())

        # This should be 1 if the template was correctly normalized
        normalization = np.sum(data)*np.radians(pixelsize)**2

        normed_data = data/normalization

        wsdl = image.get_wsdl()
        PythonUtilities.set_wsdl_weights(normed_data,wsdl)
        image.set_wsdl(wsdl)

        image.save()

    def __str__(self, indent=''):
        """Return a pretty print version of parameter values and errors."""
        p,hi_p,lo_p = self.statistical(absolute=False,two_sided=True)
        p,avg_p     = self.statistical(absolute=False,two_sided=False)
        p,abs_p     = self.statistical(absolute=True,two_sided=False)
        pnames      = self.param_names

        l=[]
        for name,val,hi,lo,avg,abs,log,free in zip(pnames,p,hi_p,lo_p,avg_p,abs_p,self.log,self.free):
            l += [ '%-10s: ' % name ]

            if log and avg != 0 and hi != 0 and lo !=0:
                l[-1] += '(1 + %.3f - %.3f) (avg = %.3f) ' % (hi,lo,avg)

            l[-1] += '%.4g'%(val)

            if not log and abs != 0: 
                l[-1] += ' +/- %.3g' % abs

            if not free: l[-1] += ' (FROZEN)'

        return ('\n'+indent).join(l)

    def full_name(self):
        return self.pretty_name

    def pretty_string(self):
        """ Default pretty string prints out the spatial parameters
            of the source in one terse line. This is useful to
            print during localization. """
        str = 'center = [ %.3fd, %.3fd ]' % (self.p[0],self.p[1])
        if len(self.p)>2:
            str+=', ext = [ %s ]' % (self.pretty_spatial_string())
        return str

    def pretty_spatial_string(self):
        """ Print out just the spatial part of the model, excluding
            the source location."""
        return ", ".join(["%.3f" % _ for _ in self.get_parameters(absolute=True)[2:]])

    def full_spatial_string(self):
        return '%.3f, %.3f, %s' % (self.p[0],self.p[1],self.pretty_spatial_string())

    @abstractmethod
    def _shrink(self):
        pass

    def shrink(self,**kwarg): 
        """ Update the spatial model so that the source is very small. This
            is useful for null hypothesis testing. """
        self.old_p       = self.p.copy()
        self.old_cov     = self.cov_matrix.copy()
        self.old_free    = self.free.copy()
        self._shrink(**kwarg)
        self.cache()
        pass

    def unshrink(self,**kwargs):
        if not hasattr(self,'old_p') or \
           not hasattr(self,'old_free') or \
           not hasattr(self,'old_free'): 
            raise Exception("Unable to unshrink source which was not first shrunk")
        
        self.p    = self.old_p 
        self.free = self.old_free
        self.cov_matrix = self.old_cov

        delattr(self,'old_p')
        delattr(self,'old_free')
        delattr(self,'old_cov')

        self.cache()

    @abstractmethod
    def can_shrink(self): 
        raise NotImplementedError('Cannot shrink %s!' % self.pretty_name)

    def r68(self,*args, **kwargs): 
        """ 68% containment radius, in degrees. """
        if hasattr(self, 'analytic_r68'):
            return self.analytic_r68(*args, **kwargs)
        else:
            return self.numeric_r68(*args, **kwargs)

    def r99(self, *args, **kwargs): 
        """ 99% containment radius, in degrees. """
        if hasattr(self, 'analytic_r99'):
            return self.analytic_r99(*args, **kwargs)
        else:
            return self.numeric_r99(*args, **kwargs)

    def effective_edge(self,energy=None):
        """ For analytic convolution, distance to be taken as the edge of the
            source. """
        return 5*self.r68()


class RadiallySymmetricModel(SpatialModel):

    def __call__(self,v,energy=None):
        """ v can be either a three element list which gets turned into
            a Hep3Vector and then a SkyDir, or it can be a SkyDir This is
            necessary to interface with PySkyFunction. """
        if type(v)==list and len(v)==3:
            skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        elif type(v)==SkyDir:
            skydir = v
        else:
            raise Exception("Incorrect argument to __call__ function.")

        return self.at_r(skydir.difference(self.center),energy)

    def at_r(self,r,energy=None):
        """ r is in radians. """
        return self.at_r_in_deg(np.degrees(r),energy)

    @abstractmethod
    def at_r_in_deg(self,r,energy=None):
        """ Should return the intensity at a distance r from the spatial model's center,
            r is in degrees. """
        pass

    def approximate_profile(self,numpoints=200):
        """ For outputting radial profile with sufficient accuracy. Rapidly varying spatial
            models should implement their own version of this function."""
        radius = np.linspace(0,self.effective_edge(),numpoints)
        pdf = self.at_r_in_deg(radius)
        return radius,pdf

    def template_npix(self):
        """ For outputting radial profile with sufficient accuracy. Rapidly varying spatial
            models should implement their own version of this function."""
        npix = 150
        return npix

    def save_profile(self, filename, *args, **kwarsg):
        """ Save out 1 1D text file with the profile.

            When we save an analytic shape to a text file and then load 
            it as a SpatialMap, we end up with a spatial model that
            is almost identical:

            First, create a uniform disk spatial model

                >>> center = SkyDir(0,0)
                >>> disk = Disk(sigma=1, center=center)

            Save it to a file:

                >>> from tempfile import NamedTemporaryFile
                >>> temp = NamedTemporaryFile()
                >>> filename = temp.name
                >>> disk.save_profile(filename)

            Then, create the spatial map:

                >>> profile = RadialProfile(file=filename)
                
            In the center, the values are the same:

                >>> np.allclose(disk(center), profile(center))
                True

            Outside the disk edge, the spatial map has intensity 0:

                >>> print (profile.effective_edge())
                1.0

                >>> print (profile(SkyDir(1.5,0)))
                0.0
            """
        radius,pdf = self.approximate_profile(*args, **kwarsg)
        open(path.expand(filename),'w').write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(radius,pdf)]))

class PseudoSpatialModel(SpatialModel):
    """ PseudoSpatialModel are point-like SpatialModels.
        i.e. they are SpatialModels which predict emission
        only from a very small region of the sky. """

    default_p, param_names, default_limits, log, steps = [], [], [], [], []


class Gaussian(RadiallySymmetricModel):
    """ Defined as a gaussian in x with width sigma
        times a gaussian in y with width ext.

       PDF = (1/2*pi*sigma)*exp(-|skydir-center|^2/2*sigma)

       p = [ ra, dec, sigma ]

       sigma = one dimensional size of the spatial model, measured in degrees

       The default sigma is 0.1 degrees:

           >>> gauss = Gaussian()
           >>> print (gauss['sigma'])
           0.1
       """

    # x68 and x99 are mathematical constants. 
    # They are the ratio of r68 (or r99), the %68 (or %99) containment
    # radius, to sigma, the 'size' parameter of an extended source.
    # Therefore, sigma*x68=r68
    x68,x99=1.50959219,3.03485426

    default_p = [0.1]
    param_names = ['Sigma']
    default_limits = [[SMALL_ANALYTIC_EXTENSION,3]]
    log = [True]
    steps = [0.04]

    def extension(self):
        # extension defined as a function so it is easy to overload
        # by the pseudo hypothesis.
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(Gaussian,self).cache()

        self.sigma=self.extension()
        self.sigma2=self.sigma**2 # cache this value
        self.pref=1/(2*np.pi*np.radians(self.sigma)**2)

    def at_r_in_deg(self,r,energy=None):
        return self.pref*np.exp(-r**2/(2*self.sigma2))

    def analytic_r68(self): return Gaussian.x68*self.sigma
    def analytic_r99(self): return Gaussian.x99*self.sigma
    def template_diameter(self): return 2*self.r99()

    def pretty_spatial_string(self):
        return "%.3fd" % (self.sigma)

    def has_edge(self): return False

    def _shrink(self,size=SMALL_ANALYTIC_EXTENSION): 
        self.p[2]=smart_log(size,log=self.log[2])
        self.free[2]=False

    def can_shrink(self): return True


class SunExtended(RadiallySymmetricModel):
    """Representes the 1/r emission from the sun (up to a 
       given maximum radius rmax which defaults to 20 degrees. 

       Note that there is no fittable extension parameter,
       although the position can be modified.

       PDF = (1/2*pi*rmax*r)*delta(r<rmax)

       Note that rmax & r must be in steradians in this formula.

       self.rmax in degrees
       self.rmax_in_rad in radians

       p = [ ra, dec ] """

    default_p, param_names, default_limits, log, steps = [], [], [], [], []

    def __init__(self, rmax=20, **kwargs):
        self.rmax = rmax
        super(SunExtended,self).__init__(**kwargs)

    def cache(self):
        super(SunExtended,self).cache()

        self.rmax_in_rad = np.radians(self.rmax)

        self.pref=1/(2*np.pi*self.rmax_in_rad)

    def at_r(self,r,energy=None):
        """ r in radians. """
        return self.pref*np.where(r<self.rmax_in_rad,1./r,0)

    def at_r_in_deg(self,r,energy=None):
        """ r is in radians. """
        return self.at_r(np.radians(r),energy)

    def analytic_r68(self): return 0.68*self.rmax
    def analytic_r99(self): return 0.99*self.rmax
    def template_diameter(self): return 2*self.rmax*(6.0/5.0)

    def pretty_spatial_string(self): return ""

    def effective_edge(self,energy=None):
        """ Disk has a well defined edge, so there is no reason to integrate past it. """
        return self.rmax

    def has_edge(self): return False

    def can_shrink(self): return False


class PseudoGaussian(PseudoSpatialModel,Gaussian):
    """ A PseudoGuassian is a Gaussian source with a fixed
        small radius. Useful to ensure that the null hypothesis
        of an extended source has the exact same PDF as the
        extended source."""

    def extension(self): return SMALL_ANALYTIC_EXTENSION

    def can_shrink(self): return False


class Disk(RadiallySymmetricModel):
    """ Defined as a constant value up to a distance Sigma away from the source. 

        It is easy to create a Disk Spatial Model

            >>> disk = Disk()

        The default size is 0.1 degrees and the default RA and DEC are 0 and 0:

            >>> disk['Sigma'] == 0.1 and disk['sigma'] == 0.1
            True
            >>> disk['RA'] == 0 and disk['DEC'] == 0
            True

        It is easy to create an object with custom parameters:

            >>> disk = Disk(p=[1.5],center=SkyDir(22,22,SkyDir.EQUATORIAL))
            >>> disk['Sigma'] == 1.5
            True

            >>> print (disk.center.ra())
            22.0
            >>> print (disk.center.dec())
            22.0

        You can also specify all parameters in a large list if you want:

            >>> disk = Disk(p=[22, 22, 1.5],coordsystem=SkyDir.EQUATORIAL)
            >>> print (disk['ra'])
            22.0
        
        And it is easy to set the parameters of a spatial model:

            >>> disk['Sigma'] = 0.5
            >>> print (disk['Sigma'])
            0.5

        But extensions can not be negative

            >>> disk['Sigma'] = -1
            Traceback (most recent call last):
                ...
            Exception: Parameter 2 must be positive!


        By default, the coordinate system is equatorial
            >>> np.allclose([disk['ra'],disk['dec']], [22, 22])
            True

        But it can be set to galactic

            >>> disk = Disk(center=SkyDir(22,22,SkyDir.GALACTIC),coordsystem=SkyDir.GALACTIC)

            >>> np.allclose([disk['l'],disk['b']], [22, 22])
            True

        If the galactic coodinates l and b are passed into the object, the coordiante system will be set automatically

            >>> disk = Disk(sigma=1.5, l=22, b=22)
            >>> disk.coordsystem == SkyDir.GALACTIC
            True
            >>> np.allclose([disk['sigma'], disk['l'], disk['b']], [1.5, 22, 22])
            True

        And when the coordinatsystem is galactic, you can not read out the RA or dec values

            >>> disk['ra']
            Traceback (most recent call last):
                ...
            Exception: Unknown parameter name ra

            >>> disk['ra'] = 1
            Traceback (most recent call last):
                ...
            Exception: Unknown parameter name ra

        The keywords to set a disk should be case insensitive

            >>> disk = Disk(SIGMA=1.5)
            >>> disk['SiGmA']
            1.5


        When you set instead the RA and DEC values, the coordiante system is equatorial

            >>> disk = Disk(p=1.5, ra=22, dec=22)
            >>> np.allclose([disk['sigma'],disk['ra'],disk['dec']], [1.5, 22, 22])
            True
            >>> disk['l']
            Traceback (most recent call last):
                ...
            Exception: Unknown parameter name l

            >>> disk['l'] = 1
            Traceback (most recent call last):
                ...
            Exception: Unknown parameter name l


        You can easily make deep copies of Spatial Models

            >>> disk = Disk(sigma=1)
            >>> copy = disk.copy()
            >>> disk['sigma'] = 2
            >>> copy['sigma'] == 1 
            True

        Spatial models perform correct checking to see if the input is valid
        
            >>> model=Disk(random_input=3)
            Traceback (most recent call last):
                ...
            Exception: Unrecognized keyword arguments ['random_input']

        The PDF for a disk spatial model is uniform inside the radius 'sigma' with a value 1/(pi*r^2):
        Note that the value of a spatial model is reasonable:

            >>> disk = Disk(sigma=1)
            >>> np.allclose(disk.at_r_in_deg(np.asarray([0,0.25,0.5,0.75,1])), 1/(np.pi*np.radians(1)**2))
            True

        Outside 'sigma', the disk has a value of 0

            >>> np.allclose(disk.at_r_in_deg(1.1),00)
            True
        """

    # See documentation in Disk for description
    x68,x99=0.824621125,0.994987437

    default_p = [0.1]
    param_names = ['Sigma']
    default_limits = [[SMALL_ANALYTIC_EXTENSION,3]]
    log = [True]
    steps = [0.04]

    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(Disk,self).cache()

        self.sigma=self.extension()
        self.sigma2=self.sigma**2 # cache this value
        self.pref=1/(np.pi*np.radians(self.sigma)**2)

    def at_r_in_deg(self,r,energy=None):
        return np.where(r<=self.sigma,self.pref,0)

    def analytic_r68(self): return Disk.x68*self.sigma
    def analytic_r99(self): return Disk.x99*self.sigma
    def template_diameter(self): return 2.0*(self['sigma']*6./5.)

    def effective_edge(self,energy=None):
        """ Disk has a well defined edge, so there is no reason to integrate past it. """
        return self.sigma

    def has_edge(self): return True

    def pretty_spatial_string(self):
        return "%.3fd" % (self.sigma)

    def _shrink(self,size=SMALL_ANALYTIC_EXTENSION): 
        self['sigma']=size
        self.free[2]=False
    def can_shrink(self): return True


class PseudoDisk(PseudoSpatialModel,Disk):
    """ A PseudoDisk is a Disk with a fixed
        small radius. Useful to ensure that the null hypothesis
        of an extended source has the exact same PDF as the
        extended source with small extension.
        
        A pseudodisk should be a disk with extension fixed to SMALL_ANALYTIC_EXTENSION:
        
            >>> x = PseudoDisk()
            >>> print (x.param_names)
            ['RA' 'DEC']
            >>> print (x.p)
            [ 0.  0.]
            >>> print (x.extension() == SMALL_ANALYTIC_EXTENSION)
            True

    """

    def extension(self): return SMALL_ANALYTIC_EXTENSION

    def can_shrink(self): return False


class Ring(RadiallySymmetricModel):
    """ The ring is defined as a constant value between one radius and another. """

    default_p = [0.1,0.5]
    param_names = ['Sigma', 'Fraction']
    default_limits = [[SMALL_ANALYTIC_EXTENSION,3],[0,1]]
    log = [True,False]
    steps = [0.04,0.1]

    def cache(self):
        super(Ring,self).cache()

        self.sigma,self.frac=self.get_parameters(absolute=True)[2:4]

        if self.frac < 0 or self.frac >= 1: 
            raise Exception("Ring spatial model must have 'frac' spatial parameter >=0 and < 1.")

        self.sigma2=self.sigma**2
        self.frac2=self.frac**2
        self.pref=1/(np.pi*np.radians(self.sigma)**2*(1-self.frac2))

    def at_r_in_deg(self,r,energy=None):
        return np.where((r>=self.frac*self.sigma)&(r<=self.sigma),self.pref,0)

    def analytic_r68(self): return Disk.x68*(1-self.frac2)+self.frac2
    def analytic_r99(self): return Disk.x99*(1-self.frac2)+self.frac2
    def template_diameter(self): return 2.0*(self['sigma']*6./5.)

    def effective_edge(self,energy=None):
        """ Disk has a well defined edge, so there is no reason to integrate past it. """
        return self.sigma

    def has_edge(self): return True

    def pretty_spatial_string(self):
        return "%.3fd, %.3f" % (self.sigma,self.frac)

    def _shrink(self,size=SMALL_ANALYTIC_EXTENSION): 
        self['sigma']=size
        self['fraction']=0
        self.free[2:4]=False
    def can_shrink(self): return True


class InterpProfile(RadiallySymmetricModel):

    default_p, param_names, default_limits, log, steps = [], [], [], [], []
    kind_dict = {'linear':1,'cubic':3}
    def __init__(self, r_in_degrees, pdf, kind='linear', **kwargs):
        """ Define a spatial model from a 1D interpolation of pdf 
            as a function of r_in_degrees.
            
            Note, the first column should be the offset angle from the center of the source
            in degrees while the second should be proportional to the intensity
            of the emission per steradian. Furthermore, the intensity is explicitly integrated
            and renormalized so that the integral over the whole emission radius is 1. 
            Any difference between the renormalized intensity profile and the input profile 
            is stored in the scalefactor property.

            First, define analytic shape
            
                >>> sigma = 1 # degrees
                >>> gauss = Gaussian(sigma=sigma)

            Next, define numeric shape

                >>> r = np.linspace(0,10,100) # in degrees
                >>> pdf = np.exp(-r**2/(2*sigma**2))
                >>> numeric_gauss = InterpProfile(r_in_degrees=r,pdf=pdf,kind='cubic')

            Note that the normalization of the numeric gaussian is wrong, but 
            will be renormalized anyway.

                >>> np.allclose(gauss.at_r_in_deg(r),pdf/numeric_gauss.scalefactor,
                ...             rtol=1e-5,atol=1e-5)
                True

            Note that spatial model is the same as Gaussian, even for oddly spaced points:

                >>> r_test = np.linspace(0,1,47) # sample oddly
                >>> np.allclose(gauss.at_r_in_deg(r_test), 
                ...             numeric_gauss.at_r_in_deg(r_test),rtol=1e-5,atol=1e-5)
                True
                >>> np.allclose([gauss.r68(),gauss.r99()], [numeric_gauss.r68(),numeric_gauss.r99()], atol=1e-3, rtol=1e-3)
                True

            This is just a doctest to test against a previous bug

                >>> print (numeric_gauss.log)
                [False False]
        """
        self.r_in_degrees, self.pdf, self.kind = r_in_degrees, pdf, kind
        self.r_in_radians = np.radians(self.r_in_degrees)

        super(InterpProfile,self).__init__(**kwargs)

        if self.r_in_degrees.shape != self.pdf.shape or len(self.r_in_degrees) != len(self.pdf):
            raise Exception("Size and shape of input arrays must be the same.")
        if self.r_in_degrees[0] != 0:
            pass
            #raise Exception("Profile must start at r=0")

    def cache(self):
        self.setup_spline()
        self.normalize()
        super(InterpProfile,self).cache()
        
    def setup_spline(self):
        """ Create a UnivariateSpline to interpolate in 1D. """
        # Create the 1D spline (in logspace if asked for)
        if self.kind == 'log':
            log_spline = UnivariateSpline(self.r_in_degrees,np.log10(self.pdf),k=1,s=0)
            spline = lambda r: 10**(log_spline(r))
        else:
            spline = UnivariateSpline(self.r_in_degrees,self.pdf,k=self.kind_dict[self.kind],s=0)

        self.spline = spline

    def normalize(self):
        # Keep track of the various scalings
        self.scalefactor = max(self.pdf)

        # Note that this formula assumes that space is flat, which
        # is incorrect. But the rest of the objects in this file
        # make that assumption, so this is kept for consistency.
        integrand = lambda r: self.spline(np.degrees(r))/self.scalefactor * 2*np.pi*r

        # perform integral in radians b/c the PDF must integrate
        # over solid angle (in units of steradians) to 1
        integral = quad(integrand, 0, self.r_in_radians[-1], full_output=True, epsabs=0)[0]

        self.scalefactor *= integral

    def at_r_in_deg(self,r,energy=None):
        return np.where(r<self.r_in_degrees[-1],self.spline(r)/self.scalefactor,0)

    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def quantile(self,quant, quad_kwargs=None):
        """ Calculate the quantiles of a pdf. Assumes flat space. """

        if quad_kwargs==None: quad_kwargs=dict(epsabs=0., epsrel=1e-10)

        integrand = lambda r: self.at_r_in_deg(r)*2*np.pi*np.radians(r)
        quantile=Quantile(integrand, 0, self.r_in_degrees[-1], 
                          quad_kwargs=quad_kwargs)
        return quantile(quant)

    def effective_edge(self,energy=None):
        """ Interpolation returns 0 outside of rmax, so no need to integrate past it. """
        return self.r_in_degrees[-1]

    def approximate_profile(self,numpoints=200):
        """ For outputting radial profile with sufficient accuracy. Rapidly varying spatial
            models should implement their own version of this function."""
        edge = self.effective_edge()
        radius=np.logspace(np.log10(1e-6*edge),np.log10(edge),numpoints)
        pdf = self.at_r_in_deg(radius)
        return radius,pdf

    def template_npix(self): return 1000

    def numeric_r68(self): return self.quantile(0.68)
    def numeric_r99(self): return self.quantile(0.99)

    def template_diameter(self): return 2.0*self.r_in_degrees[-1]*(6.0/5.0)

    def __getstate__(self):
        d=copy.copy(self.__dict__)
        del d['spline']
        return d

    def __setstate__(self,state):
        """ When unpickling the object, afterwords recreate the skymaps.SkyImage object. """
        self.__dict__ = state
        self.setup_spline()

class InterpProfile2D(InterpProfile):
    default_p = [1]
    param_names = ['Sigma'] # limits set in __init__
    log = [True]
    steps = [0.04]

    def __init__(self, r_in_degrees, sigmas, pdf2d, kind='log',**kwargs):
        """ Define a spatial model from a 2D interpolation of a pdf 
            as a function of the offset angle, r_in_degrees, and the extension, sigma.
            
            The interpolation in sigma is performed explicitly in this class while the
            interpolation in r_in_degrees is passed of to the InterpProfile base class.

            First, define analytic shape
                >>> sigma = 1 # degrees
                >>> gauss = Gaussian(sigma=sigma)

            Next, define numeric shape

                >>> r = np.linspace(0,10,100) # in degrees
                >>> s = np.linspace(0.5,1.5,500) # in degrees
                >>> rr,ss = np.meshgrid(r,s)
                >>> pdf2d = np.exp(-rr**2/(2*ss**2))
                >>> numeric_gauss = InterpProfile2D(r_in_degrees=r,sigmas=s,pdf2d=pdf2d,kind='cubic')

            Note that the normalization of the numeric gaussian is wrong, but 
            will be renormalized anyway.

            Note that spatial model is the same as Gaussian, even for oddly spaced points:

                >>> r_test = np.linspace(0,1,47) # sample oddly
                >>> sigma_test = 1.032
                >>> gauss.setp('sigma',sigma_test)
                >>> numeric_gauss.setp('sigma',sigma_test)
                >>> np.allclose(gauss.at_r_in_deg(r_test),
                ...             numeric_gauss.at_r_in_deg(r_test),rtol=1e-5,atol=1e-4)
                True
        """
        self.pdf2d, self.sigmas, self.kind = pdf2d, sigmas, kind
        self.default_limits = [[min(self.sigmas),max(self.sigmas)]]
        self.r_in_degrees,self.r_in_radians = r_in_degrees,np.radians(r_in_degrees)
        self.setup_interp()
        super(InterpProfile,self).__init__(**kwargs)

    def cache(self):
        self.sigma=self.extension()
        if self.sigma < self.default_limits[-1][0]:
            raise SpatialModelException("sigma=%s is lower than the lower limit=%s:" % (self.sigma,self.default_limits[-1][0]))
        if self.sigma > self.default_limits[-1][-1]:
            raise SpatialModelException("sigma=%s is larger than the upper limit=%s:" % (self.sigma,self.default_limits[-1][-1]))
        self.pdf = self.interp(self.r_in_degrees,self.sigma)
        super(InterpProfile2D,self).cache()

    def setup_interp(self):
        """ Create the 2D interpolation from a table of values. 
        This only needs to be called at initialization."""
        # Interpolate linearly in logspace
        r_mesh,sigma_mesh = np.meshgrid(self.r_in_degrees,self.sigmas)
        pts = np.array((r_mesh.ravel(),sigma_mesh.ravel())).T
        vals = np.log10(self.pdf2d.ravel())

        # Out of range, set interp to zero (log to -np.inf)
        log_interp = LinearNDInterpolator(pts,vals,fill_value=-np.inf)
        self.interp = lambda r,sigma: 10**(log_interp(r,sigma))

        return self.interp

    def __getstate__(self):
        d=super(InterpProfile2D,self).__getstate__()
        del d['interp']
        return d

    def __setstate__(self,state):
        """ When unpickling the object, afterwords recreate the skymaps.SkyImage object. """
        super(InterpProfile2D,self).__setstate__(state)
        self.setup_interp()


class RadialProfile(InterpProfile):
    r""" Define an extended source spatial model from a text file.

        Below is a simple example to demonstrate this spatial model.
        First, we will generate a text file consistent with a gaussain.

        First, define analytic shape
        
            >>> sigma = 1 # degrees
            >>> gauss = Gaussian(sigma=sigma)

        Next, define numeric shape

            >>> r = np.linspace(0,10,100) # in degrees
            >>> pdf = np.exp(-r**2/(2*sigma**2))
            >>> from StringIO import StringIO
            >>> file = StringIO('\n'.join('%s\t%s' % (i,j) for i,j in zip(r,pdf)))
            >>> numeric_gauss = RadialProfile(file=file,kind='cubic')

        Note that the normalization of the numeric gaussian is wrong, but 
        will be renormalized anyway.

        Note that spatial model is the same as Gaussian, even for oddly spaced points:

            >>> r_test = np.linspace(0,1,47) # sample odly
            >>> np.allclose(gauss.at_r_in_deg(r_test),
            ...             numeric_gauss.at_r_in_deg(r_test),rtol=1e-5,atol=1e-5)
            True
    """

    def __init__(self, file, **kwargs):

        self.file = file

        if not hasattr(file,'read') and not os.path.exists(file):
            raise Exception("RadialProfile must be passed an existing file")

        r_in_degrees,pdf=np.loadtxt(file,unpack=True)

        super(RadialProfile,self).__init__(r_in_degrees=r_in_degrees, pdf=pdf, **kwargs)


class EllipticalSpatialModel(SpatialModel):
    """  Defined as a gaussian in the major axis of width Major_Axis
         times a gaussian in the minor axis of width Minor_Axis
         where the major axis is at an angle theta from the

         The three parameters are

         p = [ Major_Axis, Minor_Axis, Theta ]

         They are all in degrees.

         sigma_x is the semi major axis and sigma_y is the semi minor axis.

         Just as with the 1FGL source catalog, Theta is defiend as the angle
         of the semimajor axis from celestial North, positive toward increasing RA (eastward).

         http://fermi.gsfc.nasa.gov/ssc/data/access/lat/1yr_catalog/1FGL_column_descriptions_v2.pdf
    """
    def extension(self):
        """ For overloading by pseudo hypothesis. """
        return self.get_parameters(absolute=True)[2:5]

    def cache(self):

        super(EllipticalSpatialModel,self).cache()

        self.sigma_x, self.sigma_y, self.theta = self.extension()

        self.call_grid = None

    def fill_grid(self,grid,energy=None,override_skydir=None):
        """ grid is a convolution.Grid object.

            It is assumed that the center of the grid is the same
            as the center of the extended source! This is ensured
            by the roi_extended object.

            Return an array of the pdf value on the rotated grid 
            (like the Grid.fill function). 
            
            Converstion from angles & lengths to a,b, and c parameters comes from
            http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

            Honestly though, there is an overall normalization of the
            angle and I canot intuit why it should have the paritcular
            value that it does. So the factor of 45 degrees and the minus
            sign in the angle formula I discovered simply by running the
            code a bunch of times with different angles and positions in
            the sky and creating skymaps of the PDF until I was convinced
            the formula was correct. 
            
            override_skydir is just used internally so the same function
            can be used by the __call__ function given a SkyDir object. """

        # NB, pixelsize is in degrees so dLons/dLats are in degres, which
        # is consistent with sigmax/sigma y in degrees.
        if override_skydir:
            # I have no idea why this works, and it really confuses me. But
            # I confirmed it just by making a lot of skyimages.
            x,y = grid.pix(override_skydir)
            x_c,y_c = grid.pix(grid.center)
            dLats=(x-x_c)*grid.pixelsize
            dLons=(y_c-y)*grid.pixelsize
        else:
            mpix = (float(grid.npix)-1)/2.
            dLons,dLats= np.meshgrid((np.arange(0,grid.npix)-mpix)*grid.pixelsize,
                                     (mpix-np.arange(0,grid.npix))*grid.pixelsize)

        # calculate the angle from the center of the grid to the celestial north
        # pole by finding the image coordinates for the center and for the sky
        # coordinate rotated up towards the north pole by one degree.
        towards_cel_north = SkyDir(self.center.ra(),self.center.dec()+1)
        x,y   = grid.pix(towards_cel_north)
        xc,yc = grid.pix(self.center)

        # Magic factor of 90 degrees still confuses me.
        angle = np.radians(90) - (np.radians(self.theta) + np.arctan2(y-yc,x-xc))

        a =  np.cos(angle)**2/(self.sigma_x**2)  + np.sin(angle)**2/(self.sigma_y**2)
        b = -np.sin(2*angle)/(2*self.sigma_x**2) + np.sin(2*angle)/(2*self.sigma_y**2)
        c =  np.sin(angle)**2/(self.sigma_x**2)  + np.cos(angle)**2/(self.sigma_y**2)

        x=a*dLons**2 + 2*b*dLons*dLats+ c*dLats**2

        return self.value_at(x)

    @abstractmethod
    def value_at(self,x):
        pass

    def __call__(self,v,energy=None):
        """ This code is very inefficient and should not be used for anything
            serious. 
            
            It is here for
            (a) completeness
            (b) for use by save_template where efficiency isn't important.
            
            """
        if type(v)==list and len(v)==3:
            skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        elif type(v)==SkyDir:
            skydir = v
        else:
            raise Exception("Incorrect argument to __call__ function.")

        # creating a grid for each function call is horribly inneficient, so
        # this function shouldn't be used for anything serious.
        if self.call_grid is None: 
            from uw.utilities.convolution import Grid
            self.call_grid=Grid(self.center)
            self.call_grid.wrap=True

        return self.fill_grid(self.call_grid,energy=None,override_skydir=skydir)

    def pretty_spatial_string(self):
        return "%.3fd, %.3fd, %.2fd" % \
                (self.sigma_x,self.sigma_y, self.theta)
    
    @abstractmethod
    def ellipse_68(self):
        """ Returns the parameters of an ellipse (sigma_x, sigma_y, theta)
            which has the same ellipticity/angle of the regular shape but
            encloses 68 percent of the intensity. """
        pass

    def _shrink(self,size=SMALL_NUMERIC_EXTENSION): 
        self['Major_Axis']=size
        self['Minor_Axis']=size
        self['Pos_Angle']=0
        self.free[2:5]=False
    def can_shrink(self): return True

    def __getstate__(self):
        """ Cannot pickle the call_grid object. """
        d=copy.copy(self.__dict__)
        if d.has_key('call_grid'): del d['call_grid']
        return d

    def __setstate__(self,state):
        """ When unpickling the object, afterwords recreate the skymaps.SkyImage object. """
        self.__dict__ = state
        self.call_grid = None

class EllipticalGaussian(EllipticalSpatialModel):

    # N,B default_limits or Non-radially symmetric sources dictated more by pixelization of grid.
    # Note, angle is in degrees
    default_p = [0.2,0.1,0]
    param_names = ['Major_Axis','Minor_Axis','Pos_Angle']
    default_limits = [[1e-6,3],
                      [1e-6,3],
                      [-45,45]]
    # Note for elliptical shapes, theta > 45 is the same as a negative angle
    log = [True,True,False]
    steps = [0.04,0.04,5]

    def effective_edge(self,energy=None):
        return 5*Gaussian.x68*max(self.sigma_x,self.sigma_y)

    def has_edge(self): return False

    def cache(self):
        super(EllipticalGaussian,self).cache()
        self.pref = 1/(2*np.pi*self.sigma_x*self.sigma_y)

    def value_at(self,x):
        return self.pref*np.exp(-x/2)

    def ellipse_68(self): return Gaussian.x68*self.sigma_x,Gaussian.x68*self.sigma_y,self.theta
    def ellipse_99(self): return Gaussian.x99*self.sigma_x,Gaussian.x99*self.sigma_y,self.theta
    def template_diameter(self): return 2.0*max(*self.ellipse_99()[0:2])


class PseudoEllipticalGaussian(PseudoSpatialModel,EllipticalGaussian):

    def extension(self):
        return SMALL_NUMERIC_EXTENSION,SMALL_NUMERIC_EXTENSION,0

    def pretty_spatial_string(self):
        return "%.3fd" % (self.sigma_x)

    def can_shrink(self): return False


class RadiallySymmetricEllipticalGaussian(EllipticalGaussian):

    def extension(self):
        sigma=self.get_parameters(absolute=True)[2]
        return sigma,sigma,0

    def pretty_spatial_string(self):
        return "%.3fd" % (self.sigma_x)


class EllipticalDisk(EllipticalSpatialModel):
    """ The elliptical disk is defined as. 
    
        Sanity check:

            >>> disk = Disk(sigma=1)
            >>> elliptical_disk = EllipticalDisk(major_axis=1, minor_axis=1, pos_angle=0)
            >>> np.allclose(disk(disk.center), elliptical_disk(disk.center))
            True
        """

    default_p = [0.2,0.1,0]
    param_names = ['Major_Axis','Minor_Axis','Pos_Angle']
    default_limits = [[SMALL_NUMERIC_EXTENSION,3],
              [SMALL_NUMERIC_EXTENSION,3],
              [-45,45]]
    log = [True,True,False]
    steps = [0.04,0.04,5]

    def effective_edge(self,energy=None):
        return max(self.sigma_x,self.sigma_y)

    def has_edge(self): return True

    def cache(self):
        super(EllipticalDisk,self).cache()
        self.pref = 1/(np.pi*np.radians(self.sigma_x)*np.radians(self.sigma_y))

    def value_at(self,x):
        return np.where(x<1,self.pref,0)

    def ellipse_68(self): return Disk.x68*self.sigma_x,Disk.x68*self.sigma_y,self.theta
    def ellipse_99(self): return Disk.x99*self.sigma_x,Disk.x99*self.sigma_y,self.theta
    def template_diameter(self): return 2.0*max(self['Major_Axis'],self['Minor_Axis'])*(6.0/5.0)


class RadiallySymmetricEllipticalDisk(EllipticalDisk):

    def extension(self):
        sigma=self.get_parameters(absolute=True)[2]
        return sigma,sigma,0


class PseudoEllipticalDisk(PseudoSpatialModel,EllipticalDisk):
    def extension(self):
        return SMALL_NUMERIC_EXTENSION,SMALL_NUMERIC_EXTENSION,0

    def can_shrink(self): return False


class EllipticalRing(EllipticalSpatialModel):

    default_p = [0.2,0.1,0,0.5]
    param_names = ['Major_Axis','Minor_Axis','Pos_Angle','Fraction']
    default_limits = [[SMALL_NUMERIC_EXTENSION,3],
             [SMALL_NUMERIC_EXTENSION,3],
             [-45,45],
             [0,1]]
    log = [True,True,False,False]
    steps = [0.04,0.04,5,0.1]

    def effective_edge(self,energy=None):
        """ sqrt(2) s for case of theta=45deg with
            respect to the imagined box's edge. """
        return max(self.sigma_x,self.sigma_y)

    def has_edge(self): return True

    def cache(self):
        super(EllipticalRing,self).cache()

        self.frac = self.get_parameters(absolute=True)[5]
        self.frac2 = self.frac**2

        self.pref = 1/(np.pi*np.radians(self.sigma_x)*np.radians(self.sigma_y)*(1-self.frac**2))

    def value_at(self,x):
        return np.where((x>self.frac2)&(x<1),self.pref,0)

    def ellipse_68(self):
        x68=Disk.x68*(1-self.frac2)+self.frac2
        return x68*self.sigma_x,x68*self.sigma_y,self.theta

    def ellipse_99(self):
        x99=DISK_X99*(1-self.frac2)+self.frac2
        return x99*self.sigma_x,x99*self.sigma_y,self.theta
    def template_diameter(self): return 2.0*max(self['Major_Axis'],self['Minor_Axis'])*(6.0/5.0)

    def pretty_spatial_string(self):
        return "%.3fd, %.3fd, %.2fd, %.2f" % \
                (self.sigma_x,self.sigma_y,
                 self.theta,self.frac)

    def _shrink(self): 
        self['fraction']=0
        self.free[5]=False

        # this calls the cache function
        super(EllipticalRing,self).shrink()
    def can_shrink(self): return True


class SpatialMap(SpatialModel):
    """ Implement an extended source not as a simple geometric shape
        but as from a 2 dimensional fits file. 
        
        This is analogous to gtlike's SpatialModel type SpatialMap. It
        is different in that this template is explicity normalized.
        
        A Template still has two spatial parameters, which represent a rotation of 
        the template away from the fits file's center.


        Testing: make sure the center of the spatial map is consistent.
        This was previously a bug:
        
            >>> gauss = Gaussian(p=[1.],center=SkyDir(41.234,42.345,SkyDir.EQUATORIAL))
            >>> from tempfile import NamedTemporaryFile
            >>> temp = NamedTemporaryFile()
            >>> filename = temp.name
            >>> f = gauss.save_template(filename,npix=150)
            >>> map = SpatialMap(file=filename)
            >>> np.allclose([map.center.ra(),map.center.dec()],[gauss.center.ra(),gauss.center.dec()])
            True
    """

    default_p, param_names, default_limits, log, steps = [], [], [], [], []

    def __init__(self, file, **kwargs):

        self.file = file

        super(SpatialMap,self).__init__(**kwargs)

        self.extension="" # use primary extension.
        self.interpolate=True # Note, interpolate=True necessary to not read outside array

        # The skyfun is not normalized. The normaliztaion happens later, after
        # the convolution step.
        
        self.skyfun=SkyImage(path.expand(self.file),self.extension,self.interpolate)

        projection = p = self.skyfun.projector()
        naxis1=self.skyfun.naxis1()
        naxis2=self.skyfun.naxis2()

        def dir(x,y): return SkyDir(x,y,projection)

        # Set the source center to the center of the image.
        self.center=SkyDir((naxis1+1.0)/2.0,(naxis2+1.0)/2.0,p)

        # the spatial parameters are just the center of the image.
        if self.coordsystem == SkyDir.EQUATORIAL:
            self.p = np.asarray([self.center.ra(),self.center.dec()],dtype=float)
        elif self.coordsystem == SkyDir.GALACTIC:
            self.p = np.asarray([self.center.l(),self.center.b()],dtype=float)

        self.init_p=self.p.copy()

        # SkyDir of image corners
        edges=[SkyDir(0,0,p),SkyDir(0,naxis2,p),
               SkyDir(naxis1,0,p),SkyDir(naxis1,naxis2,p)]

        # Find furthest corner
        self.edge=max(_.difference(self.center) for _ in edges)

    def full_name(self):
        return '%s (%s)' % (self.pretty_name,os.path.basename(self.file))

    def effective_edge(self,energy=None):
        return self.edge

    def __call__(self,v,energy=None):
        if type(v)==list and len(v)==3:
            skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        elif type(v)==SkyDir:
            skydir = v
        else:
            raise Exception("Incorrect argument to __call__ function.")

        return self.skyfun(skydir)

    def get_PySkyFunction(self):
        return self.skyfun

    def __getstate__(self):
        """ You cannot pickle a skymaps.SkyImage object. To avoid this,
            simply delete it before pickling. We are storing enough
            information to recreate it when unpickling. """
        d=copy.copy(self.__dict__)
        del d['skyfun']
        return d

    def __setstate__(self,state):
        """ When unpickling the object, afterwords recreate the skymaps.SkyImage object. """
        self.__dict__ = state
        try:
            self.skyfun=SkyImage(path.expand(self.file),self.extension,self.interpolate)
        except:
            self.skyfun=None


class EnergyDependentSpatialModel(SpatialModel):
    """ This is the base class for all spatial models where the extended source
        shape is a function of energy. """
    pass


def convert_spatial_map(spatial,filename):
    """ This function needed for xml parsing.
    
        Convert any SpatialModel object into a SpatialMap object. This is
        useful for ensuring compliance with gtlike.
        
        spatial is the name of the input spatial model
        filename is the filename for the saved template
        
        The return is a SpatialMap object with the same PDF. 
        
        When a SpatialMap object is passed in, only save out a new template if
        the source has been moved. This, of course, will degrade the quality
        of the saved template. """
    if isinstance(spatial,SpatialMap):
        return spatial

    spatial.save_template(filename,spatial.template_npix())
    new_map = SpatialMap(file=filename)
    return new_map


if __name__ == "__main__":
    import doctest
    doctest.testmod()
