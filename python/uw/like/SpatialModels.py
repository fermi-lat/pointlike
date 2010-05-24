"""A set of classes to implement spatial models.

   $Header$

   author: Joshua Lande

"""

import math as M
import numpy as N
from scipy import vectorize
from skymaps import PySkySpectrum,PySkyFunction,SkyDir,Hep3Vector,SkyImage,Background,WeightedSkyDirList
from skymaps import SkyIntegrator
from pointlike import DoubleVector

class DefaultSpatialModelValues(object):
    models = {
       'Gaussian'           : {'p':[M.radians(.1)],                 'param_names':['sigma']},
       'EllipticalGaussian' : {'p':[M.radians(.1),M.radians(.1),0], 'param_names':['sigma x','sigma y','theta']},
    }

    @staticmethod
    def setup(the_model):
        classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__

        for key,val in DefaultSpatialModelValues.models[classname].items():
            exec('the_model.%s = val'%key)


#===============================================================================================#

class SpatialModel(object):
    """ All spatial models are assumed to be normalized such that the
        integral over solid angle of the intensity (for a given energy)
        is equal to 1.

        All SpatialModel objects must implement the __call__ function,
        which takes a skydir object and returns the intensity at that
        direction.  

        Please do not directly modify the self.p or self.center variables
        but instead call the update() function, which will correctly
        update the internal state of the object and precompute a couple
        of useful values. """

    def init(self):
        """ This should be inhereted by child classes to setup 
            the analysis after self.p and self.center are set
            This is also called by update_params. """
        raise NotImplementedError("Subclasses should implement this!")

    def __init__(self,**kwargs):
        """ Required parameters:
              * center = A skydir of the center of the spatial model.
              * p = a list of spatial model parameters.
        """
        DefaultSpatialModelValues.setup(self)
        self.__dict__.update(**kwargs)
        if not self.p or not self.center:
            raise Exception("Spatial model %s must be initialized with parameters p and center.")

        if self.center is None:
            raise Exception('center must be pased to SpatialModel')

        if self.p is None:
            raise Exception('center must be pased to SpatialModel')

        self.init()

    def update(self,p=None,center=None):
        """ Update the spatial model parameters and center. """
        if p: self.p = p
        if center: self.center = center
        self.init()

    def __call__(self,v,energy=None):
        raise NotImplementedError("Subclasses should implement this!")

    def get_r68(self):
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

    def __str__(self):
        """Return a pretty print version of parameter values and errors."""

        p = self.p
        pnames = self.param_names

        m=max([len(n) for n in pnames])
        l=[]
        for i in xrange(len(pnames)):
            n=pnames[i][:m]
            t_n=n+(m-len(n))*' '
            l+=[t_n+': %.3g'%(p[i])]
        return '\n'.join(l)

    def pretty_string(self):
        """ Default pretty string prints out the spatial parameters
            Overload for nicer output. """
        return "[ "+" ".join(["%.3f" % _ for _ in self.p])+" ]"

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

    def get_r68(self):
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

       p = [ sigma ]

       sigma = one dimensional r68 of the spatial model, measured in radians
       """
    def init(self):
        self.sigma=self.p[0]
        self.sigma2=self.sigma**2 # cache this value
        self.pref=1/(2*M.pi*self.sigma2)

    def at_r(self,r):
        return self.pref*M.exp(-r**2/(2*self.sigma2))

    def get_r68(self):
        return 1.5*self.sigma

    def pretty_string(self):
        return "[ %.3f' ]" % (60*M.degrees(self.sigma))

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
    def init(self,*args,**kwargs):

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

    def pretty_string(self):
        return "[ %.3f', %.3f', %.1d ]" % \
                (60*M.degrees(self.sigma_x),60*M.degrees(self.sigma_y), M.degrees(theta))
