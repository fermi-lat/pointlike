"""Module to handle LAT exposure calculations.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/irfs/exposure.py,v 1.1 2016/06/22 17:02:51 wallacee Exp $
Author: Eric Wallace
"""
__version__='$Revision: 1.1 $'


import numpy as np

import skymaps

class Exposure(object):

    def __init__(self,livetime,aeff,
                 correction=None,
                 weighted_livetime=None,
                 exposure_cube=None,
                 cthetamin=.4):
        self.aeff = aeff 
        self.lt = livetime
        self.weighted_lt = weighted_livetime
        if correction is not None:
            self.correction = correction
        else:
            self.correction = lambda x: 1.0
        if exposure_cube is not None:
            self._cpp_exposure = skymaps.DiffuseFunction(exposure_cube,1000.,False)
        else:
            skymaps.Exposure.set_cutoff(cthetamin)
            if self.weighted_lt is None:
                self._cpp_exposure = skymaps.Exposure(self.lt,self.aeff)
            else:
                self._cpp_exposure = skymaps.Exposure(self.lt,self.weighted_lt,self.aeff)

    def value(self,skydir,energy):
        return self._cpp_exposure.value(skydir,energy)

    def __call__(self,skydir,energy):
        return self.value(skydir,energy)*self.correction(energy)

    def model_integral(self, skydir, func,  emin, emax):
        return self.integrator(skydir,  emin, emax)(func)

    def integrator(self,skydir,emin,emax):
        return ExposureIntegral(self,skydir,emin,emax)

    def band_exposure(self,energy):
        return BandExposure(self,energy)

class BandExposure(Exposure):

    def __init__(self,exp,energy):
        for attr in ('aeff','lt','weighted_lt','correction','_cpp_exposure'):
            setattr(self,attr,getattr(exp,attr))
        self.energy = energy
        self.correction = exp.correction(self.energy)

    def __call__(self,skydir,energy=None):
        if energy is None:
            energy = self.energy
        return Exposure.value(self,skydir,energy)


class ExposureCorrection(object):
    """ logarithmic interpolation function
    """
    def __init__(self, a,b, ea=100, eb=300):
        self.c = (b-a)/np.log(eb/ea)
        self.d =  a -self.c*np.log(ea)
        self.a, self.b = a,b
        self.ea,self.eb = ea,eb
    def __call__(self, e):
        if e>self.eb: return self.b
        if e<self.ea: return self.a
        return self.c*np.log(e) + self.d
    def plot(self, ax=None, **kwargs):
        import pylab as plt
        if ax is None: 
            ax = plt.gca()
        dom = np.logspace(1.5, 3.5, 501) 
        ax.plot(dom, map(self, dom), lw=2, color='r', **kwargs)
        plt.setp(ax, xscale='log', xlim=(dom[0], dom[-1]))


class ExposureIntegral(object):

    nsp_simps =16#4 # reduced from original 16

    def __init__(self, exp, skydir, emin, emax):
        """Calulate factors for  evaluating the counts under a given spectral model, 
            for integrating over the exposure within the energy limits
            note that exposure is evaluated at the skydir
        """
        self.sp_points = sp = np.logspace(np.log10(emin),np.log10(emax),self.nsp_simps+1)
        exp_points     = map(lambda e: exp(skydir, e), sp)
        # following may be marginally faster, but not implemented for exposure cube
        #exp_points      = np.asarray(self.exp.vector_value(self.sd,DoubleVector(sp)))
        simps_weights  = (np.log(sp[-1]/sp[0])/(3.*self.nsp_simps)) * \
                              np.asarray([1.] + ([4.,2.]*(self.nsp_simps/2))[:-1] + [1.])
        self.sp_vector = sp * exp_points * simps_weights
        
    def  __call__(self, model_function):
        """ return integral over exposure for function of differential flux 
        model_function : function
            if it returns an array, if it is a gradient, tnen axis should be 1
        """
        axis = 1 if hasattr(model_function(self.sp_points[0]), '__iter__') else None
        return (model_function(self.sp_points)*self.sp_vector).sum(axis=axis)
