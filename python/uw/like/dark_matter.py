"""A set of classes to implement dark matter analysis in pointlike

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/dark_matter.py,v 1.3 2012/02/03 22:51:24 szalewsk Exp $

    author: Joshua Lande, Alex Drlica-Wagner
"""
import numpy as np
from SpatialModels import RadiallySymmetricModel, SMALL_ANALYTIC_EXTENSION, PseudoSpatialModel

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
    limits = [[SMALL_ANALYTIC_EXTENSION,9]] # constrain r68 to 9 degrees.
    steps = [0.04]
    log = [True]

    def extension(self):
        return self.get_parameters(absolute=True)[2]

    def cache(self):
        super(NFW,self).cache()

        self.sigma=self.extension()
        # This factor of 1.07 is a normalization constant
        # coming from a comparison of the analytic form
        # with the l.o.s. integral.
        self.factor=1.07
        self.scaled_sigma=self.sigma/self.factor

    def at_r_in_deg(self,r,energy=None):
        return 2/(np.pi*r*self.scaled_sigma*(1+r/self.scaled_sigma)**5)

    def r68(self): return NFW.x68*self.scaled_sigma
    def r99(self): return NFW.x99*self.scaled_sigma

    def has_edge(self): return False

    def pretty_spatial_string(self):
        return "%.3fd" % (self.sigma)

    def _shrink(self,size=SMALL_ANALYTIC_EXTENSION): 
        self['sigma']=size
        self.free[2]=False
    def can_shrink(self): return True


class PseudoPingNFW(PseudoSpatialModel,PingNFW):
    """ The Pseudo variant of the NFW profile.

            >>> x = PseudoNFW()
            >>> print x.extension() == SMALL_ANALYTIC_EXTENSION
            True
    """

    def extension(self): return SMALL_ANALYTIC_EXTENSION

    def can_shrink(self): return False


"""
from SpatialModels import InterpProfile
class NFW(InterpProfile):

    default_p = [0.1]
    param_names = ['Sigma']
    limits = [[SMALL_ANALYTIC_EXTENSION,3]]
    log = [True]
    steps = [0.04]

   file='$(INST_DIR)/Likelihood/src/dmfit/gammamc_dif.dat'),

    def __init__(self, file, **kwargs)
        super(NFW,self).__init__(**kwargs)

    def cache(self):

        # setp was just called, so self['sigma'] has now been set

        # if you want to do the calculation here
        #self.r_in_degrees, self.pdf = self.log_and_tedious_calcuation(self['sigma'])

        # if there is a table of values
        self.r_in_degrees, self.pdf = self.short_and_not_tedious_calcuation(self['sigma'])

        # This does all the fancy interpolation stuff
        super(NFW,self).cache()
"""


if __name__ == "__main__":
    import doctest
    doctest.testmod()
