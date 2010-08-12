"""A set of classes to implement spatial models.

   $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/SpatialModels.py,v 1.11 2010/08/03 03:47:07 lande Exp $

   author: Joshua Lande

"""
import numpy as N
from scipy import vectorize
from skymaps import PySkySpectrum,PySkyFunction,SkyDir,Hep3Vector,\
        SkyImage,SkyIntegrator,CompositeSkyFunction

class DefaultSpatialModelValues(object):
    """ Spatial Parameters:
            p: the spatial parameters. There values are all assumed to be absolute.
            param_names: the names of the spatial parameters
            limits: the limits imposed on the paraemters when fitting. 
                These values are absolute. Note that by and large, everything
                is measured in radians. The limits on relative movement of
                longitude and latitude are measured in degrees.
            log: Wheter or not the parameter should be mapped into log space.
            steps: used by minuit when fitting the source. useful to make them comparable
                to distance away from true values. These are not absolute. For log parameters,
                the step is the step in the log of the parameter.

        By construction, the first two spatial parameters of all extended sources
        is the center of the source. The limits on the first two parameters are
        defined as a physical angular distance away from the source.
        The first to parametesr (lon & lat) are forced to have log=False
        The firs two parametesr are not defined in the models dict but always set
        to defaults lower in the function. """
    models = {
        'Gaussian'           : {'p':[N.radians(.1)],                 
                                'param_names':['Sigma'],                        
                                'limits':N.radians([[1e-10,3]]),
                                'log':[True],
                                'steps':[0.04]
                                # Note that the step for sigma is a step in log space!
                                # As minuit.py's doc says, a step of .04 is about 10% in log space
                               }, 
        'PseudoGaussian'     : {},
        'Disk'               : {'p':[N.radians(.1)],
                                'param_names':['Sigma'],
                                'limits':N.radians([[1e-10,3]]),
                                'log':[True],
                                'steps':[0.04]},
        'PseudoDisk'         : {},
        'NFW'                : {'p':[N.radians(.1)],
                                'param_names': ['Sigma'],
                                'limits': N.radians([[1e-10,9]]), # constrain r68 to 3 degrees.
                                'steps':[0.04],
                                'log':[True]},
        'PseudoNFW'          : {},
        'RadialProfile'      : {'p':[], 'param_names':[], 'limits':[], 'log':[], 'steps':[]},
        'EllipticalGaussian' : {'p':[N.radians(.1),N.radians(.1),0], 
                                'param_names':['Major_Axis','Minor_Axis','Position_Angle'],
                                'limits':N.radians([[1e-6,3],
                                                    [1e-6,3],
                                                    [-360,360]]),
                                'log':[True,True,False],
                                'steps':[0.04,0.04,5]},
        'SpatialMap'         : {}
    }

    @staticmethod
    def setup(the_model):
        classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__

        for key,val in DefaultSpatialModelValues.models[classname].items():
            exec('the_model.%s = val'%key)

        the_model.p=N.append([0.,0.],the_model.p) \
                if the_model.__dict__.has_key('p') else N.asarray([0.,0.])

        the_model.log=N.append([False,False,],the_model.log) \
                if the_model.__dict__.has_key('log') else N.asarray([False,False])

        the_model.param_names=N.append(['lon','lat'],the_model.param_names) \
                if the_model.__dict__.has_key('param_names') else N.asarray(['lon','lat'])

        the_model.limits=N.append([[-10.,10.],[-10.,10.]],the_model.limits,axis=0) \
                if the_model.__dict__.has_key('limits') else N.asarray([[-10.,10],[-10.,10.]])

        the_model.steps=N.append([0.1,0.1],the_model.steps) \
                if the_model.__dict__.has_key('steps') else N.asarray([0.1,0.1])

        the_model.coordsystem = SkyDir.EQUATORIAL

        the_model.cov_matrix = N.zeros([len(the_model.p),len(the_model.p)])
        the_model.free = N.asarray([True] * len(the_model.p))

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
        direction. 
        
        One slight differnece between the SpatialModel and Model class
        has to do with absolute. Always, absolute=True has the same
        meanings, pass in the true value of the parameters. But
        absolute=False has different meanings for different functions.
        
        For most functions, absolute=false means the values for which log=True
        should be in log space. This is the calse for set_parameters, 
        get_parameters, and get_cov_matrix. This is different from
        how Models works where absolute=False means all parameters are
        in log space.
        
        On the other hand, the function statistical has absolute=false
        intepreted the same convention as absolute, where absolute=false
        means return the relative error (absolute error/parameter value)
        which is useful for printing percent error, etc. """

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
            self.param_names[0:2] = ['RA','Dec']
        elif self.coordsystem == SkyDir.GALACTIC:
            self.param_names[0:2] = ['l','b']

        if self.log[0] != False or self.log[1] != False:
            raise Exception("Do not make the spatial parameters log.")

        # map the log parameters into log space.
        # careful not to take log of a negative number
        self.p = N.asarray([N.log10(p) if log else p for p,log in zip(self.p,self.log)])
        self.limits = N.asarray([N.log10(lim) if log else lim \
                                 for lim,log in zip(self.limits,self.log)])

        self.cache()

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
            self.param_names[0:2] = ['RA','Dec']
            self.p[0:2] = [center.ra(),center.dec()]
        elif cs == SkyDir.GALACTIC:
            self.param_names[0:2] = ['l','b']
            self.p[0:2] = [center.l(),center.b()]

        # Errors are no longer valid, so reset cov matrix.
        self.cov_matrix = N.zeros([len(self.p),len(self.p)]) 

    def get_parameters(self,absolute=False):
        """Return all parameters; used for spatial fitting. 
           This is different from in Models.py """
        return N.where(self.log,10**self.p,self.p) if absolute else self.p

    def get_param_names(self,absolute=True):
        if absolute:
            return self.param_names
        else:
            return ["log10(%s)" % n if log else n \
                    for n,log in zip(self.param_names,self.log)]

    def get_limits(self,absolute=False):
        ret = N.asarray([10**lim if log and absolute else lim \
                         for lim,log in zip(self.limits,self.log)])
        return ret

    def get_steps(self):
        if not self.__dict__.has_key('steps'):
            raise Exception("Spatial model %s does not have fitting step sizes defined for it." % pretty_name)
        return self.steps

    def set_parameters(self,p,absolute=False,center=None):
        """ Set all parameters; p should have length equal to number of parameters.
            Note that this API is different from in Models.py

            If center is given as an argument, the longitude and latitude are appended to the beginning of the p
            as the first two coordinaets..
        """
        if center:
            if self.coordsystem == SkyDir.EQUATORIAL:
                p = N.append([center.ra(),center.dec()],p)
            elif self.coordsystem == SkyDir.GALACTIC:
                p = N.append([center.l(),center.b()],p)

        if len(p)!=len(self.p):
            raise Exception("SpatialModel.set_parameters given the wrong number of arguments.")

        if absolute:
            # careful not to take log of a negative number
            self.p = N.asarray([N.log10(_) if log else _ for _,log in zip(p,self.log)])
        else:
            self.p = N.asarray(p)

        self.cache()
    
    def modify_loc(self,center):
        if self.coordsystem == SkyDir.EQUATORIAL:
            self.p[0:2] = [center.ra(),center.dec()]
        elif self.coordsystem == SkyDir.GALACTIC:
            self.p[0:2] = [center.l(),center.b()]

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
        self.cov_matrix = new_cov_matrix

    def get_cov_matrix(self,absolute=True):
        """Return covariance matrix."""

        jac = N.log10(N.exp(1))
        p = N.where(self.log,(10**self.p)/jac,1) if absolute else N.ones_like(self.p)
        pt=p.reshape((p.shape[0],1)) #transpose
        return p*self.cov_matrix*pt

    def get_free_errors(self,absolute=False):
        """Return the diagonal elements of the covariance matrix for free parameters."""
        return N.diag(self.get_cov_matrix(absolute))**0.5

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
            lo_abs = N.where(self.log,p-10**(self.p-errs),errs)
            hi_abs = N.where(self.log,10**(self.p+errs)-p,errs)
            return  p, \
                    hi_abs/(1. if absolute else p), \
                    lo_abs/(1. if absolute else p)


    def __call__(self,v,energy=None):
        raise NotImplementedError("Subclasses should implement this!")

    def effective_edge(self,energy=None):
        """ It is useful to know an approximate edge to the image
            It is defined as a radius such that from the center to the
            enclosed circle contains (approximatly) the entire object. """
        raise NotImplementedError("Subclasses should implement this!")

    def get_PySkyFunction(self):
        return PySkyFunction(self)

    def get_PySkySpectrum(self):
        """ Note that I am not setting the Integral function out of
            pure lazieness, since it is not needed elsewhere. """
        return PySkySpectrum(self,None)

    def save_template(self,filename,pixelsize=0.05):
        center=self.center
        image=SkyImage(center,filename,pixelsize,self.effective_edge(),1,"ZEA",
                       True if self.coordsystem == SkyDir.GALACTIC else False,False)
        skyfunction=self.get_PySkyFunction()
        image.fill(skyfunction)
        image.save()

    def __str__(self,absolute=False):
        """Return a pretty print version of parameter values and errors."""
        p,hi_p,lo_p = self.statistical(absolute=absolute,two_sided=True)
        p,avg_p     = self.statistical(absolute=absolute,two_sided=False)
        pnames      = self.param_names

        if len(pnames)==0: return 'No Spatial Parameters'# Needed for SpatialMap model.

        m=max(len(n) for n in pnames)
        l=[]
        if N.any(avg_p != 0): #if statistical errors are present   
            for i in xrange(len(pnames)):
                n=pnames[i][:m]
                t_n=n+(m-len(n))*' '
                frozen = '' if self.free[i] else '(FROZEN)'
                if not absolute:
                   l+=[t_n+': (1 + %.3f - %.3f) (avg = %.3f) %.3g %s'%(hi_p[i],lo_p[i],avg_p[i],p[i],frozen)]
                else:
                   l+=[t_n+': %.3g + %.3g - %.3g (avg = %.3g) %s'%(p[i],hi_p[i],lo_p[i],avg_p[i],frozen)]
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
        return "[ "+" ".join(["%.3f" % _ for _ in self.get_parameters(absolute=True)[2:]])+" ]"

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

        return self.at_r(skydir.difference(self.center),energy)

    def r68(self):
        raise NotImplementedError("Subclasses should implement this!")

    def at_r(self,r,energy=None):
        """ Should return the intensity at a distance r from the spatial model's center,
            where r is in radians. """
        raise NotImplementedError("Subclasses should implement this!")

    def effective_edge(self,energy=None):
        """ For analytic convolution, distance to be taken as the edge of the
            source. """
        return 5*self.r68()

#===============================================================================================#

class PseudoSpatialModel(SpatialModel):
    pass

#===============================================================================================#

class Gaussian(RadiallySymmetricModel):
    """ Defined as a gaussian in x with width sigma
        times a gaussian in y with width ext.

       PDF = (1/2*pi*sigma)*exp(-|skydir-center|^2/2*sigma)

       p = [ ra, dec, sigma ]

       sigma = one dimensional r68 of the spatial model, measured in radians
       """
    def extension(self):
        # extension defined as a function so it is easy to overload
        # by the pseudo hypothesis.
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(Gaussian,self).cache()

        self.sigma=self.extension()
        self.sigma2=self.sigma**2 # cache this value
        self.pref=1/(2*N.pi*self.sigma2)

    def at_r(self,r,energy=None):
        return self.pref*N.exp(-r**2/(2*self.sigma2))

    def r68(self):
        return 1.50959219*self.sigma

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*N.degrees(self.sigma))

#===============================================================================================#


class PseudoGaussian(PseudoSpatialModel,Gaussian):
    """ A PseudoGuassian is a Gaussian source with a fixed
        small radius. Useful to ensure that the null hypothesis
        of an extended source has the exact same PDF as the
        extended source."""
    def extension(self): return N.radians(1e-10)

#===============================================================================================#

class Disk(RadiallySymmetricModel):
    """ Defined as a constant value up to a distance Sigma away from the source. """
    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(Disk,self).cache()

        self.sigma=self.extension()
        self.sigma2=self.sigma**2 # cache this value
        self.pref=1/(N.pi*self.sigma2)

    def at_r(self,r,energy=None):
        return N.where(r<=self.sigma,self.pref,0)

    def r68(self):
        return 0.824621125*self.sigma

    def effective_edge(self,energy=None):
        """ Disk has a well defined edge, so there is no reason to integrate past it. """
        return self.sigma

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*N.degrees(self.sigma))

#===============================================================================================#

class PseudoDisk(PseudoSpatialModel,Disk):
    """ A PseudoDisk is a Disk with a fixed
        small radius. Useful to ensure that the null hypothesis
        of an extended source has the exact same PDF as the
        extended source with small extension."""
    def extension(self): return N.radians(1e-10)

#===============================================================================================#

class NFW(RadiallySymmetricModel):
    """ Ping's parameterization of the NFW Source is 
        P(x,y)=2/(pi*r*s*(1+r/s)^5) """

    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(Disk,self).cache()

        self.sigma=self.extension()
        # Ask Alex DW if you don't understand this
        self.factor=1.07
        self.scaled_sigma=self.sigma/self.factor

    def at_r(self,r,energy=None):
        return 2/(N.pi*r*self.scaled_sigma*(1+r/self.scaled_sigma)**5)

    def r68(self):
        return 0.33*self.scaled_sigma

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*N.degrees(self.sigma))

class PseudoNFW(PseudoSpatialModel,NFW):

    def extension(self): return N.radians(1e-10)


class RadialProfile(RadiallySymmetricModel):
    def __init__(self,*args,**kwargs):

        super(RadiallySymmetricModel).__init__(*args,**kwargs)

        if not self.dict.has_key('file') or not os.path.exists(self.file):
            raise Exception("RadialProfile must be passed an existing file")

        self.r,self.pdf=N.loadtxt(self.file,unpack=True)

        self.r = N.radians(self.r) # convert to radians

        # Explicitly normalize the RadialProfile.
        self.interp = N.interp1d(self.r.self.pdf,kind='cubic',bound_error=False,fill_value=0)

        r  = N.linspace(0,self.r[-1],10000)
        dr = r[1]-r[0]
        self.norm = self.interp(r)*2*N.pi*r*dr
        self.pdf /= self.norm

        # redo normalized interpolation
        self.interp = N.interp1d(self.r*self.pdf,kind='cubic',bound_error=False,fill_value=0)

    def at_r(self,r,energy=None):
        return self.interp(r)


class EllipticalGaussian(SpatialModel):
    """  Defined as a gaussian in the major axis of width Major_Axis
         times a gaussian in the minor axis of width Minor_Axis
         where the major axis is at an angle theta from the

         The three parameters are

         p = [ Major_Axis, Minor_Axis, Theta ]

         They are all in radians.

         sigma_x is the semi major axis and sigma_y is the semi minor axis.

         Just as with the 1FGL source catalog, Theta is defiend as the angle
         of the semimajor axis from celestial North, positive toward increasing RA (eastward).

         http://fermi.gsfc.nasa.gov/ssc/data/access/lat/1yr_catalog/1FGL_column_descriptions_v2.pdf
    """
    def cache(self):

        super(EllipticalGaussian,self).cache()

        self.sigma_x, self.sigma_y, self.theta = self.get_parameters(absolute=True)[2:]

        # parameters from
        # http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
        # where I have replaced theta with -theta to get a postive angle to correspond
        # with a rotation of the semi-major axis towards positive RA.

        self.a =  N.cos(theta)**2/(2*self.sigma_x**2) + N.sin(theta)**2/(2*sigma_y**2)
        self.b =  N.sin(2*theta)/(4*self.sigma_x**2)  - N.sin(2*theta)/(4*sigma_y**2)
        self.c =  N.sin(theta)**2/(2*self.sigma_x**2) + N.cos(theta)**2/(2*sigma_y**2)

        self.pref = 1/(2*N.pi*self.sigma_x*sigma_y)

    def effective_edge(self,energy=None):
        return 5*max(self.sigma_x,self.sigma_y)

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

        return self.pref*N.exp(-(self.a*delta_x**2 +
                               2*self.b*delta_x*delta_y+
                               self.c*delta_y**2))

    def pretty_spatial_string(self):
        return "[ %.3f', %.3f', %.1d ]" % \
                (60*N.degrees(self.sigma_x),60*N.degrees(self.sigma_y), N.degrees(theta))


class SpatialMap(SpatialModel):
    """ Implement an extended source not as a simple geometric shape
        but as from a 2 dimensional fits file. 
        
        This is analogous to gtlike's SpatialModel type SpatialMap. It
        is different in that this template is explicity normalized.
        
        A Template still has two spatial parameters, which represent a rotation of 
        the template away from the fits file's center."""

    def __init__(self,*args,**kwargs):

        super(SpatialMap,self).__init__(*args,**kwargs)

        if not self.__dict__.has_key('file'):
            raise Exception("Object Template must be initialized with file=template.fits keyword.")

        extension="" # use primary extension.
        interpolate=True # Note, interpolate=True necessary to not read outside array

        # The skyfun is not normalized. The normaliztaion happens later, after
        # the convolution step.
        self.skyfun=SkyImage(self.file,extension,interpolate)

        self.projection = p = self.skyfun.projector()
        naxis1=self.skyfun.naxis1()
        naxis2=self.skyfun.naxis2()

        def dir(x,y): return SkyDir(x,y,self.projection)

        self.coordsystem = SkyDir.GALACTIC if self.projection.isGalactic() else SkyDir.EQUATORIAL

        # Set the source center to the center of the image.
        self.center=SkyDir((naxis1-1)/2,(naxis2-1)/2,p)

        # Don't display any spatial parameters
        self.p=N.asarray([])
        self.param_names=N.asarray([])
        self.free=N.asarray([])
        self.log=N.asarray([])
        self.cov_matrix = N.asarray([[]])

        # SkyDir of image corners
        edges=[SkyDir(0,0,p),SkyDir(0,naxis2,p),
               SkyDir(naxis1,0,p),SkyDir(naxis1,naxis2,p)]

        # Find furthest corner
        self.edge=max(_.difference(self.center) for _ in edges)

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
