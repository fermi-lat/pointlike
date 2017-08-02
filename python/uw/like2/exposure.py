"""
Implements exposure calcuations

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/exposure.py,v 1.5 2016/11/07 03:16:33 burnett Exp $
"""
import os
import numpy as np
import skymaps
from astropy.io import fits as  pyfits


class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime.
    
    Also handles an ad-hoc exposure correction
    """

    def __init__(self, dataset, **datadict): 
        """
        Parameters
        ----------
        dataset :  DataSet object
            for CALDB, aeff, some parameters
            
        datadict['exposure_correction'] : list of strings defining functions of energy
            the correction factors to apply to front, back 
        """

        def make_exposure():
            if dataset.exposure_cube is not None:
                ## use pregenerated gtexpcube2 cube; turn off interpolation
                return [skymaps.DiffuseFunction(f,1000.,False) for f in dataset.exposure_cube]
                
            type_names = datadict.get('type_names', ('FRONT', 'BACK'))
            if type_names=='psf':
                type_names = ['PSF{}'.format(i) for i in range(4)]
            skymaps.EffectiveArea.set_CALDB(dataset.CALDBManager.CALDB)
            skymaps.Exposure.set_cutoff(np.cos(np.radians(dataset.thetacut)))
            aeff_files = dataset.CALDBManager.get_aeff()
            ok = [os.path.exists(file) for file in aeff_files]
            if not all(ok):
                raise DataSetError('one of CALDB aeff files not found: %s' %aeff_files)
            if pyfits.open(aeff_files[0])[1].name != 'EFFECTIVE AREA':
                 # new format with combined 
                self.ea  = [skymaps.EffectiveArea('', file, 'EFFECTIVE AREA_'+name) for file,name in zip(aeff_files, type_names)]
            else:
                self.ea  = [skymaps.EffectiveArea('', file) for file in aeff_files]
            if dataset.verbose: print ' -->effective areas at 1 GeV: ', \
                    ['%s: %6.1f'% (type_names[i],self.ea[i](1000)) for i in range(len(type_names))]
            
            if dataset.use_weighted_livetime and hasattr(dataset, 'weighted_lt'):
                return [skymaps.Exposure(dataset.lt,dataset.weighted_lt,ea) for ea in self.ea]
            else:
                return  [skymaps.Exposure(dataset.lt,ea) for ea in self.ea]
                
        self.exposure = make_exposure()

        correction = datadict.pop('exposure_correction', None)
        if correction is not None:
            self.correction = map(eval, correction)
            energies = [100, 1000, 10000]
            if not self.quiet:  'Exposure correction: for energies %s ' % energies
            for i,f in enumerate(self.correction):
                if not self.quiet:  ('\tfront:','\tback: ','\tdfront:', '\tdback')[i], map( f , energies)
        else:
            self.correction = lambda x: 1.0, lambda x: 1.0
            
    def value(self, sdir, energy, event_type):
        return self.exposure[event_type].value(sdir, energy)*self.correction[event_type](energy)
        
    def __call__(self, event_type, energy=1000):
        """Return a SkySpectrum-compatible object of the exposure for the given event type (e.g., front or back)
        """
        class Exposure(object):
            def __init__(self, eman, event_type, energy, correction):
                self.eman = eman
                self.et =event_type
                self.energy=energy
                self.correction = correction
            def __repr__(self):
                return 'Exposure SkySpectrum for event_type %d, energy %.0f MeV' % (self.et, self.energy)
            def __call__(self, sdir, e=None):
                if e is None: e=self.energy
                return self.eman.value(sdir, e, self.et)
            def setEnergy(self, e):
                self.energy=e
            def model_integral(self, skydir, func,  emin, emax):
                """ return the integral of func(e)*exp(e) from emin to emax
                """
                #return ExposureIntegral(self, skydir,  emin, emax)(func)
                return self.integrator(skydir,  emin, emax)(func)
                
            def integrator(self, skydir, emin, emax):
                """ return an integrator that will  evaluate func(e)*exp(e) from emin to emax
                call with func, which may return a scalar or a 1-d array
                """
                return ExposureIntegral(self, skydir,  emin, emax)
                
        return Exposure(self, event_type, energy, self.correction[event_type](energy))
  
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

