""" Dark Matter spectral models

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/darkmatter/spectral.py,v 1.6 2012/06/14 23:40:12 lande Exp $

    author: Alex Drlica-Wagner, Joshua Lande
"""
import copy
import collections

import numpy as np
from scipy.optimize.minpack import check_gradient

from uw.like.Models import Model,CompositeModel,Constant
from uw.like.SpatialModels import SpatialMap

from uw.utilities.parmap import LogMapper,LimitMapper

class DMFitFunction(Model):
    """ Wrap gtlike's DMFitFunction interface. 
    
        N.B. The bug Sheridan reported that the set_flux function 
        was not working should now be fixed:
        
            >>> model = DMFitFunction()
            >>> model.set_flux(1e-7, emin=1e3, emax=1e5)
            >>> print '%g' % model.i_flux(emin=1e3, emax=1e5)
            1e-07

        Test the getters and setters

            >>> model['sigmav']=3.14
            >>> print '%g' % model['sigmav']
            3.14

        There was previously a bug in set_parameters, 
        lets see if its fixed:

            >>> model.set_parameters(np.log10([5,500]))
            >>> print '%g' % model['sigmav']
            5
            >>> print '%g' % model['mass']
            500

        Note, the parameters which are not directly fit (like bratio) get set correctly:

            >>> model = DMFitFunction(bratio=2)
            >>> print model.dmf.getParam('bratio').getTrueValue()
            2.0
            >>> model = DMFitFunction(bratio=3)
            >>> print model.dmf.getParam('bratio').getTrueValue()
            3.0

        Test a few hard coded values, to make sure the function values are correct:

            >>> model = DMFitFunction(sigmav=1e-26, mass=100,
            ... channel0=4, channel1=1, bratio=1, norm=2.5e17)

            >>> model = DMFitFunction(norm=2.5e17, sigmav=1e-26, channel0=4,channel1=1,mass=100,bratio=1.0)

        These points agree with the fortran code.

            >>> e = [1, 10, 100, 1000, 10000, 100000 , 1000000]
            >>> dnde = [ 9.55801576e-18, 2.04105211e-16,  4.43719263e-16, 1.00123992e-16, 1.44911940e-18, 0.0, 0.0 ]
            >>> print np.allclose(model(e), dnde)
            True
    """
    default_p=[1.0, 100.]
    default_extra_params=dict(norm=1, bratio=1.0, channel0=1, channel1=1, 
                              file='$(INST_DIR)/Likelihood/src/dmfit/gammamc_dif.dat')
    param_names=['sigmav','mass']
    default_mappers=[LogMapper,LogMapper]

    def full_name(self):
        return '%s, norm=%.1f, bratio=%.1f channel0=%d, channel1=%d' % (self.pretty_name,
                                                                        self.norm, self.bratio, 
                                                                        self.channel0, self.channel1)

    def __getstate__(self):
        d=copy.copy(self.__dict__)
        del d['dmf']
        return d

    def __setstate__(self,state):
        self.__dict__ = state
        self._update()

    def _update(self):
        """ Update the DMFitFunction internally.
            This function should be called
            automatically when necessary.
        """
        if not hasattr(self,'dmf'):
            import pyLikelihood
            self.dmf=pyLikelihood.DMFitFunction()

        for i,param_name in enumerate(self.param_names):
            self.dmf.setParam(param_name,self[param_name])

        # Set the parameters which are not fixed explicitly
        self.dmf.setParam('norm',self.norm)
        self.dmf.setParam('bratio',self.bratio)
        self.dmf.setParam('channel0',self.channel0)
        self.dmf.setParam('channel1', self.channel1)

    def __init__(self,  *args, **kwargs):
        import pyLikelihood

        # the DMFitFunction must exist before __init__ is called because
        # the __init__ will call setp().
        self.dmf=pyLikelihood.DMFitFunction()
        super(DMFitFunction,self).__init__(*args,**kwargs)

        # unbound all parameters in gtlike
        for n in np.append(self.param_names,['norm','bratio','channel0','channel1']):
            self.dmf.getParam(n).setBounds(-float('inf'),float('inf'))

        self.dmf.readFunction(SpatialMap.expand(self.file))
        self._update() # update all parameters in DMFitFunction

    def setp(self, *args, **kwargs):
        super(DMFitFunction,self).setp(*args, **kwargs)
        self._update()

    def set_parameters(self, *args, **kwargs):
        super(DMFitFunction,self).set_parameters(*args, **kwargs)
        self._update()

    def set_all_parameters(self, *args, **kwargs):
        super(DMFitFunction,self).set_all_parameters(*args, **kwargs)
        self._update()

    @staticmethod
    def call_pylike_spectrum(spectrum, e):
        """ Method to call a pylikelihood spectrum given
            either a python numer or a numpy array. """
        from pyLikelihood import dArg
        if isinstance(e,collections.Iterable):
            return np.asarray([spectrum(dArg(i)) for i in e])
        else:
            return spectrum(dArg(e))

    def __call__(self,e):
        """ Return energy in MeV. This could be vectorized. """
        return DMFitFunction.call_pylike_spectrum(self.dmf, e)


class ComprehensiveModel(CompositeModel):
    """ Implements a "Comprehensive Model" needed for comparing non-nested
        models using a fequentist test with a well defiend null hypthesis.

        A good reference for the test is Cox 1961, 1962:

        And a modern description is in a recent talk by Jan Conrad:
            http://www-conf.slac.stanford.edu/statisticalissues2012/talks/SLAC_Statistics_2012_Conrad.pdf


        This feature is somewhat implemented in gtlike, but compbingin the
        two models in a different way
            https://confluence.slac.stanford.edu/display/SCIGRPS/Model+Selection+Using+Likelihood+Ratios

        Using this object is easy:

            >>> from uw.like.Models import PowerLaw
            >>> dm=DMFitFunction(); dm.set_flux(1)
            >>> pl=PowerLaw(); pl.set_flux(1)
            >>> cm=ComprehensiveModel(dm,pl)

        This model has a "Scale" (the theta parameter), and the parameters
        for the dark matter & powerlaw object:

            >>> print cm.param_names
            ['sigmav', 'mass', 'Norm', 'Index', 'Scale']


        The default 'theta' parameter is 0.5
            >>> print cm.param_names[-1]
            Scale
            >>> print cm[-1]
            0.5
            >>> print cm.theta
            0.5

        And the value is defined with the strange formula:

            >>> energies=np.logspace(1,7,7)
            >>> for theta in [0, 0.25, 0.5, 0.75, 1]:
            ...    cm.theta = theta
            ...    np.all(cm(energies)==dm(energies)**theta*pl(energies)**(1-theta))
            True
            True
            True
            True
            True

        Note, make theta the last parameter because many function in pointlike (like set_flux)
        assume that the first model is the normalization.
    """
    def __init__(self,model1,model2):
        theta=Constant(Scale=0.5,mappers=[LimitMapper(0,1)])
        super(ComprehensiveModel,self).__init__(model1,model2,theta)

    @property
    def pretty_name(self):
        g,f,theta=self.models
        return '%s^theta*%s^(1-theta)' % (g.pretty_name,f.pretty_name)

    @property
    def theta(self):
        return self.models[2]['Scale']

    @theta.setter
    def theta(self,new_theta):
        self.models[2]['Scale']=new_theta

    def __call__(self,e):
        g,f=self.models[0:2]
        return g(e)**self.theta*f(e)**(1-self.theta)

    def set_prefactor(self, prefactor, energy):
        """ set_prefactor must work around non-linearities in function

                >>> from uw.like.Models import PowerLaw,ExpCutoff
                >>> model=ComprehensiveModel(PowerLaw(),ExpCutoff())
                >>> model.set_prefactor(1e-10, 100)
                >>> print model(100)
                1e-10
        """
        g,f=self.models[0:2]
        g.setp(0, 1) 
        f.setp(0, 1) 
        new_prefactor=prefactor/self(energy)
        g.setp(0, new_prefactor)
        f.setp(0, new_prefactor)

    def set_flux(self,flux,*args,**kwargs):
        """ set_flux must work around non-linearities in function

                >>> from uw.like.Models import PowerLaw,ExpCutoff
                >>> model=ComprehensiveModel(PowerLaw(),ExpCutoff())
                >>> model.set_flux(1)
                >>> print model.i_flux()
                1.0
        """
        g,f=self.models[0:2]
        g.setp(0, 1) 
        f.setp(0, 1) 
        new_prefactor = flux/self.i_flux(*args,**kwargs)
        g.setp(0,new_prefactor)
        f.setp(0,new_prefactor)

        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
        
