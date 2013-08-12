"""A set dark matter spatial models for pointlike analyses

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/darkmatter/spatial.py,v 1.8 2012/09/27 21:36:51 lande Exp $

    author: Joshua Lande, Alex Drlica-Wagner
"""
import numpy as np
import os, pickle
from uw.utilities import path
from uw.like.SpatialModels import RadiallySymmetricModel, SMALL_ANALYTIC_EXTENSION, PseudoSpatialModel
from uw.like.SpatialModels import InterpProfile2D, smart_log

# Degrees, should come from the file...
#relpath=os.path.dirname(os.path.relpath(__file__))
#abspath=os.path.dirname(os.path.abspath(__file__))
#realpath=os.path.dirname(os.path.realpath(__file__))
#dirname=os.path.dirname(__file__)
basedir = os.path.dirname(os.path.realpath(__file__))

class NFW(InterpProfile2D):
    def __init__(self, file='%s/jvalue_table.pkl'%basedir, **kwargs):

        """ This class returns the angular profile of the line-of-sight integral
        through the NFW density profile squared. This spatial model is normalized
        such that the integral over the emission radius is 1. To retrieve a physical
        description, it is necessary to multiply the spatial distribution by the factors
        rhos**2 * rs * self.scalefactor.

        WARNING: Referring to INST_DIR is trouble if pointlike you are pointing
        to a HEAD of pointlike not in INST_DIR

        The pickle input file is created with the CVS controlled script:
        user/mdwood/python/dsph/jtable.py
        
           >>> sigma = 1
           >>> r = np.logspace(0.1,1,4)
           >>> profile = NFW(sigma=sigma)
           >>> print profile.at_r_in_deg(r)
           [  3.46679557e+01   2.73145672e+00   1.51607416e-01   6.68125546e-03]

        Note, previously this wouldn't raise an exception:
            >>> np.allclose(profile.default_limits[-1],[0.01, 15])
            True
            >>> profile['sigma'] = 0.009
            Traceback (most recent call last):
                ...
            SpatialModelException: sigma=0.009 is lower than the lower limit=0.01:
            >>> profile['sigma'] = 15.001
            Traceback (most recent call last):
                ...
            SpatialModelException: sigma=15.001 is larger than the upper limit=15.0:

        """

        infile = path.expand(file)
        data = pickle.load(open(infile,'r'))
        r_in_degrees = data['psi']
        sigmas = data['sigma']
        pdf2d = data['nfw']
        
        super(NFW,self).__init__(r_in_degrees,sigmas,pdf2d,**kwargs)

    def _shrink(self,size=None): 
        if size is None: size = min(self.sigmas)
        self.setp('Sigma',size)
        self.freeze('Sigma')

    def can_shrink(self): return True

    def effective_edge(self,energy=None):
        """ Clip at the 99% quantile """
        #return self.quantile(.99) # r99 ~ 3*self.sigma
        return min(3*self.sigma,self.r_in_degrees[-1])

class Einasto(InterpProfile2D):
    def __init__(self, file='%s/jvalue_table.pkl'%basedir, **kwargs):

        """ This class returns the angular profile of the line-of-sight integral
        through the Einasto density profile squared. This spatial model is normalized
        such that the integral over the emission radius is 1. To retrieve a physical
        description, it is necessary to multiply the spatial distribution by the factors
        rhos**2 * rs * self.scalefactor.

        The pickle input file is created with the CVS controlled script:
        user/mdwood/python/dsph/jtable.py

           >>> sigma = 1
           >>> r = np.logspace(0.1,1,4)
           >>> profile = Einasto(sigma=sigma)
           >>> print profile.at_r_in_deg(r)
           [  4.30828635e+01   2.78006838e+00   8.42642147e-02   1.01894909e-03]
        
        """

        infile = path.expand(file)
        data = pickle.load(open(infile,'r'))
        r_in_degrees = data['psi']
        sigmas = data['sigma']
        pdf2d = data['einasto']
        
        super(Einasto,self).__init__(r_in_degrees,sigmas,pdf2d,**kwargs)

    def _shrink(self,size=None): 
        if size is None: size = min(self.sigmas)
        self.setp('Sigma',smart_log(size,log=self.log[2]))
        self.freeze('Sigma')

    def can_shrink(self): return True

    def effective_edge(self,energy=None):
        """ Clip at the 99% quantile """
        #return self.quantile(.99) # r99 ~ 5.2*self.sigma
        return min(5.2*self.sigma,self.r_in_degrees[-1])

class Burkert(InterpProfile2D):
    def __init__(self, file='%s/jvalue_table.pkl'%basedir, **kwargs):

        """ This class returns the angular profile of the line-of-sight integral
        through the Burkert density profile squared. This spatial model is normalized
        such that the integral over the emission radius is 1. To retrieve a physical
        description, it is necessary to multiply the spatial distribution by the factors
        rhos**2 * rs * self.scalefactor.

        The pickle input file is created with the CVS controlled script:
        user/mdwood/python/dsph/jtable.py

           >>> sigma = 1
           >>> r = np.logspace(0.1,1,4)
           >>> profile = Burkert(sigma=sigma)
           >>> print profile.at_r_in_deg(r)
           [  1.39250630e+02   1.22157313e+01   6.14035455e-01   2.42861785e-02]
            
        """

        infile = path.expand(file)
        data = pickle.load(open(infile,'r'))
        r_in_degrees = data['psi']
        sigmas = data['sigma']
        pdf2d = data['burkert']
        
        super(Burkert,self).__init__(r_in_degrees,sigmas,pdf2d,**kwargs)

    def _shrink(self,size=None): 
        if size is None: size = min(self.sigmas)
        self.setp('Sigma',smart_log(size,log=self.log[2]))
        self.freeze('Sigma')

    def can_shrink(self): return True

    def effective_edge(self,energy=None):
        """ Clip at the 99% quantile """
        #return self.quantile(.99) # r99 ~ 2.6*self.sigma
        return min(2.6*self.sigma,self.r_in_degrees[-1])

### DEPRICATED: Spatial models using an approximate formulation.
class PingNFW(RadiallySymmetricModel):
    """ Ping's parameterization of the NFW Source is 
        P(x,y)=2/(pi*r*s*(1+r/s)^5) 
        WARNING: This only works for sources with 
        0.5 <~ sigma <~ 10. For more details, see:
        https://confluence.slac.stanford.edu/display/SCIGRPS/Pointlike+DMFit+Validation
        """

    # See documentation in Disk for description
    # I got this from Wolfram Alpha with: 'solve int 2/(pi*x/1.07*(1+x*1.07)^5)*2*pi*x for x from 0 to y=0.68'
    x68,x99=0.30801306,2.02082024

    default_p = [0.1]
    param_names = ['Sigma']
    default_limits = [[SMALL_ANALYTIC_EXTENSION,9]] # constrain r68 to 9 degrees.
    steps = [0.04]
    log = [True]

    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(PingNFW,self).cache()

        self.sigma=self.extension()
        # This factor of 1.07 is a normalization constant
        # coming from a comparison of the analytic form
        # with the l.o.s. integral.
        self.factor=1.07
        self.scaled_sigma=self.sigma/self.factor

    def at_r_in_deg(self,r,energy=None):
        return 2/(np.pi*r*self.scaled_sigma*(1+r/self.scaled_sigma)**5)

    def analytic_r68(self): return NFW.x68*self.scaled_sigma
    def analytic_r99(self): return NFW.x99*self.scaled_sigma

    def has_edge(self): return False

    def pretty_spatial_string(self):
        return "%.3fd" % (self.sigma)

    def _shrink(self,size=SMALL_ANALYTIC_EXTENSION): 
        self['sigma']=size
        self.free[2]=False
    def can_shrink(self): return True

class PseudoPingNFW(PseudoSpatialModel,PingNFW):
    """ The Pseudo variant of the NFW profile.

            >>> x = PseudoPingNFW()
            >>> print x.extension() == SMALL_ANALYTIC_EXTENSION
            True
    """

    def extension(self): return SMALL_ANALYTIC_EXTENSION

    def can_shrink(self): return False


if __name__ == "__main__":
    import doctest
    doctest.testmod()
