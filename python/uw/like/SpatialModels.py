"""A set of classes to implement spatial models.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/SpatialModels.py,v 1.1 2010/05/24 08:10:30 lande Exp $

   author: Joshua Lande

"""

import math as M
import numpy as N
from scipy import vectorize
from skymaps import PySkySpectrum,PySkyFunction,SkyDir,Hep3Vector,SkyImage,Background,WeightedSkyDirList
from skymaps import SkyIntegrator,SkyDir
from pointlike import DoubleVector

class DefaultSpatialModelValues(object):
    models = {
       'Gaussian'           : {'p':[0,0,M.radians(.1)],                 'param_names':['lon','lat','sigma']},
       'PseudoGaussian'     : {'p':[0,0],                               'param_names':['lon','lat']},
       'EllipticalGaussian' : {'p':[0,0,M.radians(.1),M.radians(.1),0], 'param_names':['lon','lat','sigma x','sigma y','theta']},
       'Template'           : {'p':[0,0],                               'param_names':['lon','lat']},
    }

    @staticmethod
    def setup(the_model):
        classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__

        for key,val in DefaultSpatialModelValues.models[classname].items():
            exec('the_model.%s = val'%key)

        the_model.cov_matrix = N.zeros([len(the_model.p),len(the_model.p)]) #default covariance matrix
        the_model.free = N.asarray([True] * len(the_model.p))

        the_model.p=N.asarray(the_model.p)
        the_model.param_names=N.asarray(the_model.param_names)

        the_model.coordsystem = SkyDir.EQUATORIAL


#===============================================================================================#

class SpatialModel(object):
    """ This class represents a normalized spatial model which can be
        parameterized by a a simple geometric function with a list of
        free paraemters.
    
        All spatial models are assumed to be normalized such that the
        integral over solid angle of the intensity (for a given energy)
        is equal to 1.

        All SpatialModel objects must implement the __call__ function,
        which takes a skydir object and returns the intensity at that
        direction. """

    def __init__(self,**kwargs):
        DefaultSpatialModelValues.setup(self)

        self.__dict__.update(**kwargs)

        # if center is passed as a flag, add it to the paraemters.
        if 'center' in kwargs.keys(): 
            center = kwargs.pop('center')
            if not kwargs.has_key('p'):
                # remove first two elements from default parameters.
                self.p = self.p[2:]
            if self.coordsystem == SkyDir.EQUATORIAL:
                self.p = N.append([center.ra(),center.dec()],self.p)
            elif self.coordsystem == SkyDir.GALACTIC:
                self.p = N.append([center.l(),center.b()],self.p)

        # rename coordsystem to properly reflect projection.
        if self.coordsystem == SkyDir.EQUATORIAL:
            self.param_names[0:2] = ['ra','dec']
        elif self.coordsystem == SkyDir.GALACTIC:
            self.param_names[0:2] = ['l','b']

        self.cache()

    def cache(self):
        """ This should be inhereted by child classes to cache
            various useful quanitites after update is called. """

        self.center = SkyDir(self.p[0],self.p[1],self.coordsystem)

    def get_parameters(self):
        """Return FREE parameters; used for spatial fitting."""
        return self.p[self.free]

    def get_param_names(self):
        return self.param_names[self.free]

    def set_parameters(self,p,center=None):
        """ Set FREE parameters; p should have length equal to number of free parameters.

            If center is given as an argument, it is appended to the beginning of the p
            as the first two coordinaets..
        
        """
        if center:
            if self.coordsystem == SkyDir.EQUATORIAL:
                p = N.append([center.ra(),center.dec()],p)
            elif self.coordsystem == SkyDir.GALACTIC:
                p = N.append([center.l(),center.b()],p)

        if len(p)!=(self.free).sum():
            raise Exception("SpatialModel.set_parameters given the wrong number of arguments.")

        self.p[self.free] = p.astype(float) # downcast to float needed?

        self.cache()

    def freeze_position(self,freeze=True):
        """Freeze the source position. """
        self.freeze([0,1],freeze)

    def freeze(self,parameter,freeze=True):
        """Freeze one of the spatial parameters from fitting.
      
            parameter: a parameter name or index.
            freeze   : if True, freeze parameter; if False, free it """
        if type(parameter) == type(''):
            for n,name in enumerate(self.param_names):
                if parameter == name: parameter = n; break
        self.free[parameter] = not freeze

    def set_cov_matrix(self,new_cov_matrix):
        self.cov_matrix[N.outer(self.free,self.free)] = N.ravel(new_cov_matrix)

    def get_cov_matrix(self):
        """Return covariance matrix."""
        return self.cov_matrix

    def get_free_errors(self):
        """Return the diagonal elements of the covariance matrix for free parameters."""
        return N.diag(self.cov_matrix)[self.free]**0.5

    def statistical(self,absolute=False):
        """Return the parameter values and fractional statistical errors.
           If no error estimates are present, return 0 for the fractional error."""
        p = self.p 
        try: #see if error estimates are present
            err = self.get_free_errors()
            if absolute:
                return p,err
            else:
                return p,err/p
        except:
            return p,N.zeros_like(p)


    def __call__(self,v,energy=None):
        raise NotImplementedError("Subclasses should implement this!")

    def r68(self):
        """ It is useful to know the average spatial model size. """
        raise NotImplementedError("Subclasses should implement this!")

    def get_PySkyFunction(self):
        return PySkyFunction(self)

    def get_PySkySpectrum(self):
        """ Note that I am not setting the Integral function out of
            pure lazieness, since it is not needed elsewhere. """
        return PySkySpectrum(self,None)

    def save_template(self,filename,diameter=8,pixelsize=.125,galactic=True):
        center=self.center
        image=SkyImage(center,filename,pixelsize,diameter,1,"ZEA",galactic,False)
        skyfunction=self.get_PySkyFunction()
        image.fill(skyfunction)
        image.save()

    def __str__(self,absolute=False):
        """Return a pretty print version of parameter values and errors."""
        p,avg = self.statistical(absolute=absolute)
        pnames = self.param_names

        m=max([len(n) for n in pnames])
        l=[]
        if N.any(avg != 0): #if statistical errors are present   
            for i in xrange(len(pnames)):
                n=pnames[i][:m]
                t_n=n+(m-len(n))*' '
                frozen = '' if self.free[i] else '(FROZEN)'
                if not absolute:
                    l+=[t_n+': (1 +/- %.3f) %.3g %s'%(avg[i],p[i],frozen)]
                else:
                    l+=[t_n+': %.3g +/- %.3g %s'%(p[i],avg[i],frozen)]
            return '\n'.join(l)
        else: #if no errors are present
            for i in xrange(len(pnames)):
                n=pnames[i][:m]
                t_n=n+(m-len(n))*' '
                l+=[t_n+': %.3g'%(p[i])]
            return '\n'.join(l)

    def pretty_string(self):
        """ Default pretty string prints out the spatial parameters
            of the source in one terse line. This is useful to
            print during localization. """
        str = 'dir = (%.3f,%.3f)' % (self.p[0],self.p[1])
        if len(self.p)>2:
            str+=', ext = %s' % (self.pretty_spatial_string())
        return str

    def pretty_spatial_string(self):
        """ Print out just the spatial part of the model, excluding
            the source location."""
        return "[ "+" ".join(["%.3f" % _ for _ in self.p[2:]])+" ]"

#===============================================================================================#

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

        return self.at_r(skydir.difference(self.center))

    def r68(self):
        raise NotImplementedError("Subclasses should implement this!")

    def at_r(self,r):
        """ Should return the intensity at a distance r from the spatial model's center,
            where r is in radians. """
        raise NotImplementedError("Subclasses should implement this!")

#===============================================================================================#

class Gaussian(RadiallySymmetricModel):
    """ Defined as a gaussian in x with width sigma
        times a gaussian in y with width ext.

       PDF = (1/2*pi*sigma)*exp(-|skydir-center|^2/2*sigma)

       p = [ ra, dec, sigma ]

       sigma = one dimensional r68 of the spatial model, measured in radians
       """
    def extension(self):
        return self.p[2]

    def cache(self):
        super(Gaussian,self).cache()

        self.sigma=self.extension()
        self.sigma2=self.sigma**2 # cache this value
        self.pref=1/(2*M.pi*self.sigma2)

    def at_r(self,r):
        return self.pref*M.exp(-r**2/(2*self.sigma2))

    def r68(self):
        return 1.5*self.sigma

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*M.degrees(self.sigma))

#===============================================================================================#

class PseudoGaussian(Gaussian):
    """ A PseudoGuassian is a Gaussian source with a fixed
        small radius. Useful to ensure that the null hypothesis
        of an extended source has the exact same PDF as the
        extended source."""
    def extension(self): return M.radians(1e-5)

#===============================================================================================#

class EllipticalGaussian(SpatialModel):
    """  Defined as a gaussian in the major axis of width sigma_M
         times a gaussian in the minor axis of width sigma_m
         where the major axis is at an angle theta from the

         The three parameters are

         p = [ sigma_x, sigma_y, theta ]

         They are all in radians.

         sigma_x is the semi major axis and sigma_y is the semi minor axis.

         Just as with the 1FGL source catalog, Theta is defiend as the angle
         of the semimajor axis from celestial North, positive toward increasing RA (eastward).

         http://fermi.gsfc.nasa.gov/ssc/data/access/lat/1yr_catalog/1FGL_column_descriptions_v2.pdf
    """
    def cache(self):

        super(EllipticalGaussian,self).cache()

        sigma_x, sigma_y, theta = self.p

        # parameters from
        # http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
        # where I have replaced theta with -theta to get a postive angle to correspond
        # with a rotation of the semi-major axis towards positive RA.

        self.a =  M.cos(theta)**2/(2*sigma_x**2) + M.sin(theta)**2/(2*sigma_y**2)
        self.b =  M.sin(2*theta)/(4*sigma_x**2)  - M.sin(2*theta)/(4*sigma_y**2)
        self.c =  M.sin(theta)**2/(2*sigma_x**2) + M.cos(theta)**2/(2*sigma_y**2)

        self.pref = 1/(2*M.pi*sigma_x*sigma_y)


    def __call__(self,v,energy=None):
        if type(v)==list and len(v)==3:
            skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        elif type(v)==SkyDir:
            skydir = v
        else:
            raise Exception("Incorrect argument to __call__ function.")

        def sign(x,y,z):
            """ Return the sign of the difference between x and y (sign(y-x))
                with the added complication that x and y are the values mod z.
                Doesn't worry about edge case where x=y."""
            temp=(x-y) % z # calculate the difference mod z
            temp-= (temp>z/2)*z # no get the interval b/n -z/2 and z/2
            if temp>=0: return 1 # if this is > 0, then x > y
            return -1 # otherwise x < y

        # the semi major axis should point towards Celestial North, so it is ~ (dec_point - dec_center)
        delta_x = skydir.difference(SkyDir(skydir.ra(),self.center.dec()))*sign(skydir.dec(),self.center.dec(),180)

        delta_y = skydir.difference(SkyDir(self.center.ra(),skydir.dec()))*sign(skydir.ra(),self.center.ra(),360)

        return self.pref*M.exp(-(self.a*delta_x**2 +
                               2*self.b*delta_x*delta_y+
                               self.c*delta_y**2))

    def pretty_spatial_string(self):
        return "[ %.3f', %.3f', %.1d ]" % \
                (60*M.degrees(self.sigma_x),60*M.degrees(self.sigma_y), M.degrees(theta))


class Template(SpatialModel):

    def cache():
        if not self.__dict__.has_key('template'):
            raise Exception("Object Template must be initialized with template=template.fits argument.")

        self.skyfun = ConvertFitsFileToSkyFunction(self.template)

