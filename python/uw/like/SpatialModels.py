"""A set of classes to implement spatial models.

   $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/SpatialModels.py,v 1.17 2010/08/27 06:02:16 lande Exp $

   author: Joshua Lande

"""
import numpy as N
from scipy import vectorize
from skymaps import PySkySpectrum,PySkyFunction,SkyDir,Hep3Vector,\
        SkyImage,SkyIntegrator,CompositeSkyFunction
from uw.utilities.convolution import Grid

# Mathematical constants. They are the ratio of r68 (the %68 containment
# radius) to sigma, the 'size' parameter of an extended source.
GAUSSIAN_X68=1.50959219
DISK_X68=0.824621125
NFW_X68=0.33

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
        'Gaussian'           : {'p':[N.radians(0.1)],                 
                                'param_names':['Sigma'],                        
                                'limits':N.radians([[1e-10,3]]),
                                'log':[True],
                                'steps':[0.04]
                                # Note that the step for sigma is a step in log space!
                                # As minuit.py's doc says, a step of .04 is about 10% in log space
                               }, 
        'PseudoGaussian'     : {},
        'Disk'               : {'p':[N.radians(0.1)],
                                'param_names':['Sigma'],
                                'limits':N.radians([[1e-10,3]]),
                                'log':[True],
                                'steps':[0.04]},
        'PseudoDisk'         : {},
        'Ring'               : {'p':[N.radians(0.1),0.5],
                                'param_names':['Sigma','Fraction'],
                                'limits':N.radians([[1e-10,3],[0,1]]),
                                'log':[True,False],
                                'steps':[0.04,0.1]},
        'NFW'                : {'p':[N.radians(.1)],
                                'param_names': ['Sigma'],
                                'limits': N.radians([[1e-10,9]]), # constrain r68 to 3 degrees.
                                'steps':[0.04],
                                'log':[True]},
        'PseudoNFW'          : {},
        'RadialProfile'      : {},
        # Note, angle is in degrees
        'EllipticalGaussian' : {'p':[N.radians(0.2),N.radians(0.1),0], 
                                'param_names':['Major_Axis','Minor_Axis','Position_Angle'],
                                'limits':N.radians([[1e-6,3],
                                                    [1e-6,3],
                                                    [N.radians(-45),N.radians(45)]]),
                                # Note for elliptical shapes, theta > radians(45) is the same as a negative angle
                                'log':[True,True,False],
                                'steps':[0.04,0.04,N.radians(5)]},
        # N,B limits or Non-radially symmetric sources dictated more by pixelization of grid.
        'EllipticalDisk'     : {'p':[N.radians(0.2),N.radians(0.1),0], 
                                'param_names':['Major_Axis','Minor_Axis','Position_Angle'],
                                'limits':N.radians([[1e-3,3],
                                                    [1e-3,3],
                                                    [N.radians(-45),N.radians(45)]]),
                                'log':[True,True,False],
                                'steps':[0.04,0.04,N.radians(5)]},
        'EllipticalRing'     : {'p':[N.radians(0.2),N.radians(0.1),0,0.5], 
                                'param_names':['Major_Axis','Minor_Axis','Position_Angle','Fraction'],
                                'limits':N.radians([[1e-3,3],
                                                    [1e-3,3],
                                                    [N.radians(-45),N.radians(45)],
                                                    [0,1]]),
                                'log':[True,True,False,False],
                                'steps':[0.04,0.04,N.radians(5),0.1]},
        'SpatialMap'         : {}
    }

    models['PseudoEllipticalGaussian'] = models['PseudoGaussian']
    models['RadiallySymmetricEllipticalGaussian'] = models['Gaussian']
    models['PseudoEllipticalDisk'] = models['PseudoDisk']
    models['RadiallySymmetricEllipticalDisk'] = models['Disk']

    @staticmethod
    def setup(the_model):
        classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__

        the_model.p,the_model.log,the_model.param_names,the_model.steps=[],[],[],[]

        for key,val in DefaultSpatialModelValues.models[classname].items():
            exec('the_model.%s = val'%key)

        # Add in point source parts.
        the_model.p=N.append([0.,0.],the_model.p)
        the_model.log=N.append([False,False],the_model.log)
        the_model.param_names=N.append(['lon','lat'],the_model.param_names)
        the_model.limits=N.append([[-10.,10.],[-10.,10.]],the_model.limits,axis=0) \
                if the_model.__dict__.has_key('limits') else N.asarray([[-10.,10],[-10.,10.]])
        the_model.steps=N.append([0.1,0.1],the_model.steps)

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

        iscopy = kwargs.pop('iscopy', False)

        self.__dict__.update(**kwargs)

        if not iscopy:
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

            if len(self.p) != len(self.param_names):
                raise Exception("SpatialModel set with wrong number of parameters.")

            # rename coordsystem to properly reflect projection.
            if self.coordsystem == SkyDir.EQUATORIAL:
                self.param_names[0:2] = ['RA','Dec']
            elif self.coordsystem == SkyDir.GALACTIC:
                self.param_names[0:2] = ['l','b']

            if self.log[0] != False or self.log[1] != False:
                raise Exception("Do not make the spatial parameters log.")

            # map the log parameters into log space.
            # careful not to take log of a negative number
            self.p = N.where(self.log,N.log10(self.p),self.p)
            log_transpose = self.log.reshape((self.log.shape[0],1))
            self.limits = N.where(log_transpose,
                                  N.log10(self.limits),self.limits)

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
            return ["log10(%s)" % n.replace('Dec','Dec (rotated)') \
                    if log else n for n,log in zip(self.param_names,self.log)]

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
            raise Exception("SpatialModel.set_parameters given the wrong number of parameters.")

        self.p = N.where(self.log,N.log10(p),p) if absolute else N.asarray(p)
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

    def copy(self):
        
        a = eval(self.name+'(iscopy=True, **self.__dict__)') #create instance of same spectral model type
        
        a.p = N.asarray(self.p).copy() #copy in parameters
        a.param_names = N.asarray(self.param_names).copy() # copy in names
        a.limits = N.asarray(self.limits).copy() #copy in limits
        a.log = N.asarray(self.log).copy() #copy in log
        a.steps = N.asarray(self.steps).copy() #copy in steps

        try: a.cov_matrix = self.cov_matrix.__copy__()
        except: pass

        return a

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

    def save_template(self,filename,pixelsize=0.025):
        if isinstance(self,EnergyDependentSpatialModel):
            raise Exception("Unable to save template for energy dependent SpatialModel.")
        center=self.center
        # NB, effective_edge is radius, SkyImage takes in diameter
        image=SkyImage(center,filename,pixelsize,2*N.degrees(self.effective_edge()),1,"ZEA",
                       True if self.coordsystem == SkyDir.GALACTIC else False,False)
        skyfunction=self.get_PySkyFunction()
        image.fill(skyfunction)
        image.save()

    def __str__(self,absolute=False, indent=''):
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
                if p[i] == 0: 
                   l+=[t_n+': (1 +       -      ) (avg =      ) %.3g %s'%(p[i],frozen)]
                elif not absolute:
                   l+=[t_n+': (1 + %.3f - %.3f) (avg = %.3f) %.3g %s'%(hi_p[i],lo_p[i],avg_p[i],p[i],frozen)]
                else:
                   l+=[t_n+': %.3g + %.3g - %.3g (avg = %.3g) %s'%(p[i],hi_p[i],lo_p[i],avg_p[i],frozen)]
            return ('\n'+indent).join(l)
        else: #if no errors are present
            for i in xrange(len(pnames)):
                n=pnames[i][:m]
                t_n=n+(m-len(n))*' '
                l+=[t_n+': %.3g'%(p[i])]
            return ('\n'+indent).join(l)

    def pretty_string(self):
        """ Default pretty string prints out the spatial parameters
            of the source in one terse line. This is useful to
            print during localization. """
        str = 'dir = [ %.3fd, %.3fd ]' % (self.p[0],self.p[1])
        if len(self.p)>2:
            str+=', ext = %s' % (self.pretty_spatial_string())
        return str

    def pretty_spatial_string(self):
        """ Print out just the spatial part of the model, excluding
            the source location."""
        return "[ "+", ".join(["%.3f" % _ for _ in self.get_parameters(absolute=True)[2:]])+" ]"

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
        return GAUSSIAN_X68*self.sigma

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
        return DISK_X68*self.sigma

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

class Ring(RadiallySymmetricModel):
    """ The ring is defined as a constant value between one radius and another. """
    def cache(self):
        super(Ring,self).cache()

        self.sigma,self.frac=self.get_parameters(absolute=True)[2:4]

        if self.frac < 0 or self.frac >= 1: 
            raise Exception("Ring spatial model must have 'frac' spatial parameter >=0 and < 1.")

        self.sigma2=self.sigma**2
        self.frac2=self.frac**2
        self.pref=1/(N.pi*self.sigma2*(1-self.frac2))

    def at_r(self,r,energy=None):
        return N.where((r>=self.frac*self.sigma)&(r<=self.sigma),self.pref,0)

    def r68(self):
        return DISK_X68*(1-self.frac2)+self.frac2

    def effective_edge(self,energy=None):
        """ Disk has a well defined edge, so there is no reason to integrate past it. """
        return self.sigma

    def pretty_spatial_string(self):
        return "[ %.3f', %.3f' ]" % (60*N.degrees(self.sigma),(60*N.degrees(self.frac*self.sigma)))

#===============================================================================================#

class NFW(RadiallySymmetricModel):
    """ Ping's parameterization of the NFW Source is 
        P(x,y)=2/(pi*r*s*(1+r/s)^5) """

    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(NFW,self).cache()

        self.sigma=self.extension()
        # Ask Alex DW if you don't understand this
        self.factor=1.07
        self.scaled_sigma=self.sigma/self.factor

    def at_r(self,r,energy=None):
        return 2/(N.pi*r*self.scaled_sigma*(1+r/self.scaled_sigma)**5)

    def r68(self):
        return NFW_X68*self.scaled_sigma

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*N.degrees(self.sigma))

#===============================================================================================#

class PseudoNFW(PseudoSpatialModel,NFW):

    def extension(self): return N.radians(1e-10)

#===============================================================================================#

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

#===============================================================================================#

class EllipticalSpatialModel(SpatialModel):
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
            value that it does. So the factor of 45degrees and the minus
            sign in the angle formula I discovered simply by running the
            code a bunch of times with different angles and positions in
            the sky and creating skymaps of the PDF until I was convinced
            the formula was correct. 
            
            override_skydir is just used internally so the same function
            can be used by the __call__ function given a SkyDir object. """

        mpix = (float(grid.npix)-1)/2.
        if override_skydir:
            # I have no idea why this works, and it really confuses me. But
            # I confirmed it just by making a lot of skyimages.
            x,y = grid.pix(override_skydir)
            dLats=(x-mpix)*N.radians(grid.pixelsize)
            dLons=(mpix-y)*N.radians(grid.pixelsize)
        else:
            dLons,dLats= N.meshgrid((N.arange(0,grid.npix)-mpix)*N.radians(grid.pixelsize),
                                     (mpix-N.arange(0,grid.npix))*N.radians(grid.pixelsize))

        # calculate the angle from the center of the grid to the celestial north
        # pole by finding the image coordinates for the center and for the sky
        # coordinate rotated up towards the north pole by one degree.
        towards_cel_north = SkyDir(self.center.ra(),self.center.dec()+1)
        x,y   = grid.pix(towards_cel_north)
        xc,yc = grid.pix(self.center)

        # Magic factor of 90 degrees still confuses me.
        angle = N.radians(90) - (self.theta + N.arctan2(y-yc,x-xc))

        a =  N.cos(angle)**2/(self.sigma_x**2)  + N.sin(angle)**2/(self.sigma_y**2)
        b = -N.sin(2*angle)/(2*self.sigma_x**2) + N.sin(2*angle)/(2*self.sigma_y**2)
        c =  N.sin(angle)**2/(self.sigma_x**2)  + N.cos(angle)**2/(self.sigma_y**2)

        x=a*dLons**2 + 2*b*dLons*dLats+ c*dLats**2

        return self.value_at(x)

    def value_at(self,x):
        raise NotImplementedError("Subclasses should implement this!")

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
            self.call_grid=Grid(self.center)
            self.call_grid.wrap=True

        return self.fill_grid(self.call_grid,energy=None,override_skydir=skydir)

    def pretty_spatial_string(self):
        return "[ %.3f', %.3f', %.2fd ]" % \
                (60*N.degrees(self.sigma_x),60*N.degrees(self.sigma_y), N.degrees(self.theta))
    
    def ellipse_68(self):
        """ Returns the parameters of an ellipse (sigma_x, sigma_y, theta)
            which has the same ellipticity/angle of the regular shape but
            encloses 68 percent of the intensity. """
        raise NotImplementedError("Subclasses should implement this!")
#===============================================================================================#

class EllipticalGaussian(EllipticalSpatialModel):

    def effective_edge(self,energy=None):
        return 5*max(self.sigma_x,self.sigma_y)

    def cache(self):
        super(EllipticalGaussian,self).cache()
        self.pref = 1/(2*N.pi*self.sigma_x*self.sigma_y)

    def value_at(self,x):
        return self.pref*N.exp(-x/2)

    def ellipse_68(self):
        return GAUSSIAN_X68*self.sigma_x,GAUSSIAN_X68*self.sigma_y,self.theta

#===============================================================================================#

class PseudoEllipticalGaussian(PseudoSpatialModel,EllipticalGaussian):

    def extension(self):
        return N.radians(1e-3),N.radians(1e-3),0

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*N.degrees(self.sigma_x))

#===============================================================================================#

class RadiallySymmetricEllipticalGaussian(EllipticalGaussian):

    def extension(self):
        sigma=self.get_parameters(absolute=True)[2]
        return sigma,sigma,0

    def pretty_spatial_string(self):
        return "[ %.3f' ]" % (60*N.degrees(self.sigma_x))

#===============================================================================================#

class EllipticalDisk(EllipticalSpatialModel):
    """ The elliptical disk is defined as. """

    def effective_edge(self,energy=None):
        return max(self.sigma_x,self.sigma_y)

    def cache(self):
        super(EllipticalDisk,self).cache()
        self.pref = 1/(N.pi*self.sigma_x*self.sigma_y)

    def value_at(self,x):
        return N.where(x<1,self.pref,0)

    def ellipse_68(self):
        return DISK_X68*self.sigma_x,DISK_X68*self.sigma_y,self.theta

#===============================================================================================#

class RadiallySymmetricEllipticalDisk(EllipticalDisk):
    def extension(self):
        sigma=self.get_parameters(absolute=True)[2]
        return sigma,sigma,0

#===============================================================================================#

class PseudoEllipticalDisk(PseudoSpatialModel,EllipticalDisk):
    def extension(self):
        return N.radians(1e-3),N.radians(1e-3),0

#===============================================================================================#

class EllipticalRing(EllipticalSpatialModel):

    def effective_edge(self,energy=None):
        """ sqrt(2) s for case of theta=45deg with
            respect to the imagined box's edge. """
        return max(self.sigma_x,self.sigma_y)

    def cache(self):
        super(EllipticalRing,self).cache()

        self.frac = self.get_parameters(absolute=True)[5]
        self.frac2 = self.frac**2

        self.pref = 1/(N.pi*self.sigma_x*self.sigma_y*(1-self.frac**2))

    def value_at(self,x):
        return N.where((x>self.frac2)&(x<1),self.pref,0)

    def ellipse_68(self):
        x68=DISK_X68*(1-self.frac2)+self.frac2
        return x68*self.sigma_x,x68*self.sigma_y,self.theta

    def pretty_spatial_string(self):
        return "[ %.3f', %.3f', %.2fd, %.2f ]" % \
                (60*N.degrees(self.sigma_x),60*N.degrees(self.sigma_y),
                 N.degrees(self.theta),self.frac)

#===============================================================================================#

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

    def get_PySkyFunction():
        return self.skyfun(skydir)

#===============================================================================================#

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
        
        The return is a SpatialMap object with the same PDF. """
    spatial.save_template(filename)
    nm = SpatialMap(file=filename)
    return nm
