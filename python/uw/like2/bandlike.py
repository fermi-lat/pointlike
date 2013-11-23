"""
Manage spectral and angular models for an energy band to calculate the likelihood, gradient
   Currently delegates some computation to classes in modules like.roi_diffuse, like.roi_extended
   
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/bandlike.py,v 1.38 2013/11/20 23:09:36 burnett Exp $
Author: T.Burnett <tburnett@uw.edu> (based on pioneering work by M. Kerr)
"""

import sys, types
import numpy as np
from scipy import misc, optimize
from  uw.utilities import keyword_options
from . import roimodel, tools

class FitPlotMixin(object):
    """mixin  for likelihood function to generate a plot, or set of all plots"""
    
    def estimate_solution(self):
        """ return a tuple with:
            current parameters, 
            estimated parmeters at maximum
            sigmas
        """
        parz = self.get_parameters()
        hess = self.hessian(parz)
        cov = hess.I
        sigs = np.sqrt(np.asarray(cov.diagonal()).flatten())
        parmax = parz-self.gradient(parz)*sigs**2/2.
        return parz, parmax, sigs

    def plot_fit(self, index, ax=None, nolabels=False , y2lim=(-10,10)):
        """make a plot showing the log likelihood and its derivative as a function of
        expected sigma, evaluated from the second derivative at the current point
        
        index : int
            index of the parameter
        """
        import matplotlib.pyplot as plt
        # get current parameters, gradient, and the Hessian for estimate of max liklihood position
        parz, parmax, sigs = self.estimate_solution()
        pz = parz[index]
        part = parz.copy()
        def func(x):
            part[index]=x
            return -self(part) #(restore sign for minimization)
        def gradf(x):
            part[index]=x
            return self.gradient(part)[index]
        x0, sig = (parmax)[index], sigs[index]
        ref = func(x0)
        if ax is None:
            fig, ax = plt.subplots( figsize=(3,3))
        else: fig = ax.figure
        xsig = np.linspace(-3, 3, 27)
        x =  x0 + xsig * sig 
        ax.plot(xsig, map(func,x)-ref, '-')
        ax.plot(xsig, -((x-x0)/sig)**2, '--')
        ax.plot((pz-x0)/sig, func(pz)-ref, 'db')
        plt.setp(ax, ylim=(-9,1), xlim=(-4,4))
        if not nolabels: ax.set_ylabel('log likelihood')
        j = np.arange(len(self.mask))[self.mask][index]if hasattr(self,'mask') else index
        ax.set_title('#%d: %s' %(j,self.parameter_names[index]), size=10)
        ax.text( -2,-8, 'par %.3f\nsig %.3f' % (pz,sig), size=10)
        ax.axvline(0, color='k', ls = ':')
        
        ax2 = ax.twinx()
        gradvals = -sig*np.array(map(gradf, x))
        ax2.plot(xsig, gradvals, '-r')
        ax2.axhline(0, color='r', ls=':')
        ax2.set_ylim( y2lim)
        if not nolabels: ax2.set_ylabel('derivative (sig units)')
        else: ax2.set_yticklabels([])
        
        self.set_parameters(parz) # restore when done
        return fig
    
    def plot_all_fits(self, perrow=5, figsize=(12,12)):
        """
        """
        import matplotlib.pyplot as plt
        n = len(self.parameters)
        if n==1:
            return self.plot_fit(0)
        fig, axx = plt.subplots((n+perrow-1)/perrow,perrow, 
            figsize=figsize, sharex=True, sharey=True)
        for i, ax in enumerate(axx.flatten()):
            if i>=n: break
            self.plot_fit(i, ax = ax, nolabels=True)
        return fig


class FitterMixin(object):
    def maximize(self,  **kwargs):
        """Maximize likelihood and estimate errors.
        """
        from scipy import optimize
        quiet = kwargs.pop('quiet', True)
        if not quiet: print 'using optimize.fmin_l_bfgs_b with parameter bounds %s\n, kw= %s'% (
                            self.bounds, kwargs)
        parz = self.get_parameters()
        ret = optimize.fmin_l_bfgs_b(self, parz, 
                bounds=self.bounds,  fprime=self.gradient, **kwargs)
        if ret[2]['warnflag']>0: 
            self.set_parameters(parz) #restore if error   
            raise Exception( 'Fit failure:\n%s' % ret[2])
        if not quiet:
            print ret[2]
        f = ret 
        cov = self.hessian(f[0]).I
        diag = cov.diagonal().copy()
        bad = diag<0
        if np.any(bad):
            if not quiet: print 'Minimizer warning: bad errors for values %s'\
                %np.asarray(self.parameter_names)[bad] 
            diag[bad]=np.nan
        return f[1], f[0], np.array(np.sqrt(diag)).flatten()
        
    def check_gradient(self, delta=1e-5):
        """compare the analytic gradient with a numerical derivative"""
        
        parz = self.get_parameters()
        fz = self(parz)
        grad = self.gradient(parz)
        fprime=[]
        for i in range(len(parz)):
            parz[i]+=delta
            fplus = self(parz)
            assert abs(fplus-fz)>1e-6, 'Fail consistency: variable %d not changing' % i
            parz[i]-=2*delta
            fminus = self(parz)
            parz[i]+= delta
            fzero = self(parz)
            assert abs(fzero-fz)<1e-2, 'Fail consistency: %e, %e ' % (fzero, fz)
            fprime.append((fplus-fminus)/(2*delta))
        return grad, np.array(fprime) 
   
class BandLike(object):
    """ manage the likelihood calculation for a band 
    """
    defaults = (
        ('quiet', True, 'set False for info'),
        )
    @keyword_options.decorate(defaults)
    def __init__(self, band, sources, free=None, **kwargs):
        """
           band    : ROIband object
           sources : list of sources.Source objects
           free    : [array of bool | None]
                to select models with variable parameters
                If None, select all
        """
        keyword_options.process(self, kwargs)
        # make a list of the Response objects
        self.bandsources = np.array(map(lambda s: s.response(band, quiet=self.quiet), sources))

        self.band = band 
        self.exposure_factor = band.exposure.correction
        self.data = band.pix_counts  if band.has_pixels else []# data from the band
        self.pixels=len(self.data)
        
        self.initialize(free)
        self.update()
        # special code to unweight if galactic diffuse too large
        self.unweight = self.make_unweight()
         
    def make_unweight(self):
        """ return an unweighting factor <=1.0 to use to multiply the log likelihood
        assume that first BandSource object has the galactic diffuse
        
        systematic : float
            a fraction representing the relative systematic uncertainty in the galactic diffuse
        """
        try:
            if hasattr(self['ring'].source, 'systematic'): 
                systematic = self['ring'].source.systematic
            else:  return 1.0
        except: return 1.0
        if systematic==0: return 1.0 
        n = 1/systematic**2
        # m is the number of counts from the galactic diffuse in the footprint of a point source
        m = self['ring'].counts / (self.band.psf(0)[0]*self.band.solid_angle)
        u = min(1., n/m)
        return u

    def __str__(self):
        b = self.band
        return '%s.%s: %d models (%d free) applied to band %.0f-%.0f, %s with %d pixels, %d photons'\
                % (self.__module__,self.__class__.__name__,len(self.bandsources), sum(self.free), b.emin, b.emax, 
                 ('front back'.split()[b.event_type]), self.pixels, sum(self.data), )
                 
    def __repr__(self): return self.__str__()
    
    def __getitem__(self, i): 
        """ return a BandSource object refererence, either by index, or by source name"""
        if type(i)==types.StringType:
            t = list(self.bandsources)
            for bs in t:
                if i==bs.source.name:
                    return bs
            raise Exception('Source "%s" not found in band sources' %i)
        return self.bandsources[i]
        
    def initialize(self, free):
        """ should only call if free array changes.
            Saves the combined prediction from the models with fixed parameters
        """
        self.free = free if free is not None else np.ones(len(self.bandsources), bool)
        self.free_sources = self.bandsources[self.free]
        self.counts = self.fixed_counts = sum([b.counts for b in self.bandsources[ ~ self.free]])
        if not self.band.has_pixels: return # no data, should quit?
        self.fixed_pixels = np.zeros(self.pixels)
        for m in self.bandsources[ ~ free]:
            self.fixed_pixels += m.pix_counts
        self.model_pixels = self.fixed_pixels.copy()
        for m in self.free_sources:
            self.model_pixels += m.pix_counts
        
    def update(self, reset=False, **kwargs):
        """ assume that parameters have changed. Update only contributions 
        from models with free parameters. *must* be called before evaluating likelihood.
        reset: bool, optional
            if True, need to reinitialize variable source(s), for change of position or shape
        """
        if self.band.has_pixels: self.model_pixels[:]=self.fixed_pixels
        self.counts = self.fixed_counts
        force = kwargs.get('force', False) # set to defeat
        for bandsource in self.free_sources:
            if reset: 
                bandsource.initialize()
                bandsource.source.changed=False
            elif bandsource.source.changed or force:
                bandsource.update()
            if self.band.has_pixels: self.model_pixels += bandsource.pix_counts
            self.counts+= bandsource.counts
 
        self.weights = self.data / self.model_pixels

    def log_like(self):
        """ return the Poisson extended log likelihood """
        try:
            pix = np.sum( self.data * np.log(self.model_pixels) )  if self.pixels>0 else 0
            w = pix - self.counts * self.exposure_factor
            return self.unweight * w
        except FloatingPointError, e:
            print '%s: Floating point error %s evaluating likelihood for band %s' % (e,self,)
            raise
    
    def gradient(self):
        """ gradient of the likelihood with resepect to the free parameters
        """
        if not self.band.has_pixels: return None
        return self.unweight * np.concatenate(
                [m.grad(self.weights, self.exposure_factor) for m in self.free_sources]
            )
       
    def model_counts(self, sourcemask=None):
        """ return the model predicted counts for all or a subset of the sources
        sourcemask : array of bool
            select a set of sources. In order of sources from factory
        """
        if sourcemask is not None:
            assert len(sourcemask)==len(self), 'bad input to model_counts'
        t = np.array([s.counts for s in self])
        return sum(t) if sourcemask is None else sum(t[sourcemask])

    def add_source(self, source):
        """ add a new source 
            source: sources.Source object
        """
        # ugly but compact
        t = list(self.bandsources)
        t.append(source.response(self.band))
        self.bandsources = np.array(t)
        
    def del_source(self, source):
        """ remove the source """
        t = list(self.bandsources)
        for bs in t:
            if source.name==bs.source.name:
                t.remove(bs)
                self.bandsources = np.array(t)
                return
        raise Exception('source "%s" not found to delete' % source.name)
    
    def fill_grid(self, sdirs):
        """ fill a grid with values, which are counts/sr, so must be multiplied by the pixel size
        """
        t = np.zeros(len(sdirs))
        for m in self.bandsources:
            t+= m.fill_grid(sdirs)
        return t
        
    def dataframe(self, **kw):
        """ return a pandas.DataFrame for diagnostics """
        import pandas as pd
        df = pd.DataFrame(dict([(s.source.name, 
                dict(counts=s.counts.round(), 
                    overlap=s.overlap, 
                    free=self.free[i],
                    extended=s.source.skydir is not None and hasattr(s.source, 'dmodel'),
                    diffuse=s.source.skydir is None,
                    distance=np.degrees(s.band.skydir.difference(s.source.skydir)) if s.source.skydir is not None else 0 ),
                    )
                for i,s in enumerate(self)])
            ).T
        return df

    def counts_in_pixel(self, source_index, skydir):
        """ return a tuple of predicted signal and background counts in the pixel corresponding to skydir
        Note that if the pixel has no data, it will not have been computed for the model; instead
        this will return the average background
 
        source_index : int
            the index of the source
        skydir : SkyDir object
        """
        from skymaps import WeightedSkyDir
        band = self.band
        if not hasattr(self, 'pixel_indeces'): 
            self.pixel_indeces = list([band.b.index(x) for x in band.wsdl])
        source = self[source_index]
        hp_index = band.b.index(skydir)
        try:
            pixel_index = self.pixel_indeces.index(hp_index)
            signal = source.counts * source.pixel_values[pixel_index]
            back  = self.model_pixels[pixel_index] - signal
        except:
            # missing pixel; calculate the expected source counts
            # and estimate total by mean of model in ROI
            signal = source.flux_value(skydir) * source.counts
            back = self.model_pixels.mean()
        return signal, back 
       
         
class BandLikeList(list):
    
    def __init__(self, roi_bands, roi_sources):
        """ create a list, one per band, of BandLike objects for a given ROI 
        Provide some useful methods 
        
        parameters
        ----------        
        roi_bands : list of ROIBand objects
        roi_sources :  sourcelist.SourceList object
            a list of Source objects 
        """

        self.sources = roi_sources
        self.bands = roi_bands
        for band in roi_bands:
            self.append( BandLike(band, self.sources, self.sources.free) )
            
        self.set_selected(self)# set selected for a subset?
        self.all_energies = self.energies[:]
        
    #  the band selection mechanism, used by log_like, update and gradient   
    def set_selected(self, values):
        """setter for the property selected, which must be a subset of self"""
        if values in self: # single guy
            self._selected = [values]
            return
        assert set(values).issubset(self), 'Improper selection'
        self._selected = values
    def get_selected(self):
        return self._selected
    selected = property(get_selected, set_selected)
        
    def select(self, index=None, event_type=None):
        """ Select an energy band or bands
        parameters:
        ----------
        index: None or integer
            an index into the list of energies; if None, select all bands
            and use the current spectral model, otherwise a powerlaw to 
            represent an model-independent flux over the band.
        event_type : None or integer
            if None, select both front and back, otherwise 0/1 for front/back
                
        """
        if index==None: #select all (initially selected) bands, use input model
            selected_bands = self[:]
        else:
            energy = self.all_energies[index]
            type_select = lambda x : True if event_type is None else x==event_type
            selected_bands = filter(lambda b: b.band.energy==energy and type_select(b.band.event_type), self)
            assert len(selected_bands)>0, 'did not find any bands for energy %.1f' % energy
        self.selected = selected_bands

    @property
    def free_sources(self):
        """ list of sources with currently free parameters
        """
        return self.sources.free
    
    @property
    def energies(self):
        return  np.sort(list(set([ sm.band.energy for sm in self._selected])))
    @property
    def emin(self):
        return np.array([b.band.emin for b in self.__selected]).min()
    @property
    def emax(self):
        return np.array([b.band.emax for b in self.__selected]).max()
        
    def __repr__(self):
        sel = '%d bands'%len(self.bands) if len(self.selected)==len(self) else '%d / %d selected bands'\
            %(len(self.selected),len(self))
        return '%s.%s: \n\t%s\n\t%s\n\tParameters: %d in %d/%d free sources' % (
            self.__module__, self.__class__.__name__,
            self.sources, sel, len(self.sources.parameters), 
            sum(self.free_sources),  len(self.free_sources))

    def initialize(self):
        assert False, 'needed?'
        
    # the following methods sum over the current set of bands
    def log_like(self, summed=True):
        """log likelihood for current set of bands
        summed : bool, optional
        if false, return the array of likelihods for each band
        """
        r = np.array([b.log_like() for b in self._selected])
        return  sum(r) if summed else r
        
    def update(self, **kwargs):
        for b in self._selected: 
            b.update(**kwargs)
        self.sources.parameters.clear_changed()
        
    def gradient(self):
        return np.array([blike.gradient() for blike in self._selected]).sum(axis=0) 
        
    def hessian(self, mask=None, delta=1e-6):
        """ return a hessian matrix based on the current parameter set
        This makes a numerical derivative of the analytic gradient, so not exactly
        symmetric, but the the result must be (nearly) symmetric.
        
        mask : [None, array of bool]
            If present, must have dimension of the parameters, will generate a sub matrix
        
        For sigmas and correlation coefficients, invert to covariance
                cov =  self.hessian().I
                sigs = np.sqrt(cov.diagonal())
                corr = hess / np.outer(sigs,sigs)
        """
        # get the source parameter management object
        parameters = self.sources.parameters
        parz = parameters.get_parameters()
        if mask is None: mask = np.ones(len(parz),bool)
        else:
            mask = np.asarray(mask)
            assert len(mask)==len(parz)
        # initial values for the likelihood and gradient
        fzero = self.log_like()
        glast = gzero = self.gradient()[mask]
        t = []
        for i in np.arange(len(parz))[mask]:
            # increment current variable and get new gradient
            parameters[i] = parz[i]+delta
            self.update()
            gnow = self.gradient()[mask]
            # numerical derivative of gradient with respect to this parameter
            t.append( (gnow-glast)/(2*delta))
            glast = gnow
        hess = np.matrix(t)
        parameters.set_parameters(parz) #restore all parameters, check that no difference
        self.update()
        assert abs(fzero-self.log_like())<1e-2
        return hess 
       
        
    def energy_flux_view(self, source_name, energy=None, **kw):
        """ a functor for a source, which returns log likelihood as a 
                function of the differential energy flux, in eV units, at the given energy
                
        parameters
        ----------
        source_name : string
        energy : [None | float]
            if None, use the reference energy e0
        """

        class EnergyFluxLikelihood(object):
            def __init__(self, blike, func, energy):
                self.func = func
                self.saved_pars = self.func.parameters
                source = self.func.source
                model = source.spectral_model
                if energy is None:
                    energy=model.e0
                eflux = model(energy) * energy**2 * 1e6
                self.ratio = model[0]/eflux
                assert model[0]==model['norm']
                self.tointernal = model.mappers[0].tointernal
                self.bound = model.bounds[0][0]
        
            def restore(self):
                self.func.set_parameters(self.saved_pars)

            @tools.ufunc_decorator # make this behave like a ufunc
            def __call__(self, eflux):
                if eflux<=0:
                    par = self.bound
                else:
                    par = max(self.bound, self.tointernal(eflux*self.ratio))
                return -self.func([par])
                
        try:
            func = self.fitter_view(source_name+'_Norm')
        except Exception, msg:
            raise Exception('could not create energy flux function for source %s;%s' %(source_name, msg))
        return EnergyFluxLikelihood(self, func, energy, **kw)
        
    def fitter_view(self, select=None, **kwargs):
        """ return a object to use with a fitter.
            Two versions, one with full set of parameters, other if a subset is specified
        """

        class FitFunction(FitPlotMixin, FitterMixin): 
            blike = self
            def __init__(self, blike,  **kwargs):
                self.blike = blike
                self.parameters = blike.sources.parameters
                self.sources = self.parameters.free_sources
            def get_parameters(self):
                return self.parameters.get_parameters()
            def set_parameters(self, pars):
                self.parameters.set_parameters(pars)
                self.blike.update()
            @property 
            def bounds(self):
                return np.concatenate([s.model.bounds[s.model.free] for s in self.sources]) 
            def __call__(self, pars=None):
                if pars is not None: self.set_parameters(pars)
                return -self.blike.log_like()
            def log_like(self, summed=True):
                """assume that parameters are set, possibility of individual likelihoods"""
                return self.blike.log_like(summed)

            def gradient(self,pars):
                self.set_parameters(pars)
                return self.blike.gradient()
            def hessian(self, pars):
                self.set_parameters(pars)
                return self.blike.hessian()
            @property
            def parameter_names(self):
                return self.parameters.parameter_names
             
        class LikeFunctor(roimodel.ParSubSet, FitPlotMixin, FitterMixin):
            def __init__(self, blike, select=None):
                self.blike = blike
                super(LikeFunctor, self).__init__(blike.sources, select)
            def __repr__(self):
                return '%s.%s: %s '% (self.__module__, self.__class__.__name__, self.selection_description)
            @property
            def parameters(self):
                return self.get_parameters()
            def set_parameters(self, pars):
                super(LikeFunctor,self).set_parameters(pars)
                self.blike.update()
            def __call__(self, pars=None):
                if pars is not None: self.set_parameters(pars)
                return -self.blike.log_like()
            def log_like(self, summed=True):
                """assume that parameters are set, possibility of individual likelihoods"""
                return self.blike.log_like(summed)
            def gradient(self, pars):
                self.set_parameters(pars)
                return self.blike.gradient()[self.mask]
            def hessian(self,pars):
                self.set_parameters(pars)
                return self.blike.hessian(self.mask) # needS mask

        if select is None:
            return FitFunction(self, **kwargs)
        return LikeFunctor(self, select, **kwargs)
        
    def tsmap_view(self, source_name, **kw):
        """Return TSmap function for the source
        """
        class TSmap(object):
            def __init__(self, blike, func, **kw):
                self.func = func
                self.source = self.func.source
                self.skydir = self.source.skydir
                self.wzero = func.log_like()

            def __repr__(self):
                return '%s.%s: source' % (self.__module__, self.__class__.__name__, self.source)
                
            def restore(self):
                self.source.skydir = self.skydir

            def __call__(self, skydir):
                self.source.changed=True
                self.source.skydir = skydir
                blike.update(reset=True)
                return 2*(self.func.log_like()-self.wzero)
        
        try:
            func = self.fitter_view(source_name+'_Norm')
        except Exception, msg:
            raise Exception('could not create tsmap function for source %s;%s' %(source_name, msg))
        return TSmap(self, func, **kw)
       
    
