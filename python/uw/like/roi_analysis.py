"""
Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
on the PSF.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_analysis.py,v 1.132 2013/05/07 04:17:43 lande Exp $

author: Matthew Kerr, Toby Burnett, Joshua Lande
"""

import numpy as np
import math, pickle, collections
import numbers

from . import roi_bands, roi_localize , specfitter
from . pointspec_helpers import PointSource,get_default_diffuse_mapper
from . roi_diffuse import ROIDiffuseModel,DiffuseSource
from . roi_extended import ExtendedSource,BandFitExtended,ROIExtendedModel
from . import roi_extended
from . import roi_printing
from . import roi_modify
from . import roi_save
from . import roi_image
from . import sed_plotter
from . import roi_upper_limits
from . import mapplots
from uw.utilities import keyword_options
from uw.utilities import xml_parsers
from uw.utilities import region_writer
from uw.utilities import results_writer
from scipy.optimize import fmin,fmin_powell,fmin_bfgs

from . import roi_plotting 
from . import counts_plotter

# special function to replace or extend a docstring from that of another function
def decorate_with(other_func, append=False, append_init=False):
    """ append_init: If decorating with an object (which has an __init__ function),
                     then append the doc for the __init__ after the doc for the
                     overall class. """
    def decorator(func):
        if append: func.__doc__ += other_func.__doc__ 
        else:      func.__doc__  = other_func.__doc__ 

        if append_init and hasattr(other_func,'__init__'):
                func.__doc__ += other_func.__init__.__doc__
        return func
    return decorator


class ROIAnalysis(object):

    defaults = (
        ("fit_emin",None,"""a two-element list giving front/back minimum energy. 
                            Independent energy ranges for front and back. 0th position is for event class 0.
                            If not specified, default is read from SpectralAnalysis."""),
        ("fit_emax",None,"A two-element list giving front/back maximum energy. Same qualifications as fit_emax."),
        ("conv_type",-1,"Integer specifying front(0), back(1), or all(-1) events"),
        ("quiet",False,'Set True to suppress (some) output'),
        ("catalog_aperture",-1,"Pulsar catalog analysis only -- deprecate"),
        ("phase_factor",1.,"Pulsar phase. A correction that can be made if the data has undergone a phase selection -- between 0 and 1"),
        ("bracketing_function",None,"A function by which to multiply the counts in each band, e.g. 1.2 for a 20% 'high' IRF.  It should be a function taking as arguments energy and conversion type."),
        ("skip_setup", False, "Set True to instantiate with data only" ),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi_dir,ps_manager,ds_manager,spectral_analysis,**kwargs):
        """ roi_dir    -- the center of the ROI
            ps_manager -- an instance of ROIPointSourceManager
            ds_manager -- an instance of ROIDiffuseManager
        """
        keyword_options.process(self, kwargs)

        self.__dict__.update(**kwargs)
        self.roi_dir = roi_dir
        self.psm  = ps_manager
        self.dsm  = ds_manager
        self.bgm  = self.dsm ### bgm is deprecated, set here for convenience of interactive access

        self.sa    = spectral_analysis
        self.logl = None
        self.prev_logl = None
        self.__setup_bands__()
        self.__warn_about_binning__()
        self.bin_centers = np.sort(list(set([b.e for b in self.bands])))
        self.bin_edges    = np.sort(list(set([b.emin for b in self.bands] + [b.emax for b in self.bands])))

        self.param_state, self.param_vals  = None,None
        if self.skip_setup: 
            if not self.quiet: print ('No likelihood setup done: skip_setup is set')
            return
        self.__pre_fit__()
        
        self.logl = self.prev_logl =-self.logLikelihood(self.get_parameters()) # make sure everything initialized

    def setup(self):
        """ stuff that would not be done if skip_setup was set"""
        self.__setup_counts__()
        self.__pre_fit__()
        self.logl = self.prev_logl =-self.logLikelihood(self.get_parameters()) # make sure everything initialized
        
    def __setup_bands__(self):

        if self.fit_emin is None: self.fit_emin = self.sa.emin
        if self.fit_emax is None: self.fit_emax = self.sa.emax

        if isinstance(self.fit_emin,numbers.Real):
            self.fit_emin = [self.fit_emin]*2

        if isinstance(self.fit_emax,numbers.Real):
            self.fit_emax = [self.fit_emax]*2

        self.bands = collections.deque()
        band_kwargs = {'catalog_aperture':self.catalog_aperture,
                    'bracketing_function':self.bracketing_function,
                    'phase_factor':self.phase_factor}
        for band in self.sa.pixeldata.dmap:
            evcl = band.event_class() & 1 # protect high bits

            if ((band.emin() + 1) >= self.fit_emin[evcl] and
                (band.emax() - 1) < self.fit_emax[evcl] and
                (self.conv_type==-1 or self.conv_type==evcl)) :
                self.bands.append(roi_bands.ROIBand(band,self.sa,
                                                    self.roi_dir,**band_kwargs))

        self.bands = np.asarray(self.bands)

        if not self.skip_setup: self.__setup_counts__()
        
    def __setup_counts__(self):
        self.psm.setup_initial_counts(self.bands)
        self.dsm.setup_initial_counts(self.bands)

    def __warn_about_binning__(self):
        """ Add a user friendly warning if the binning selected
            by emin & emin in DataSpecification is inconsistent with
            the energy bins selected using fit_emin and fit_emax. """
        if self.quiet: return

        for ct in [0,1]:
            if len([b for b in self.bands if b.ct==ct]) == 0:
                print ("Warning: No conversion type %s photons were selected." % ct)
                continue
            actual_emin=min(b.emin for b in self.bands if b.ct==ct)
            actual_emax=max(b.emax for b in self.bands if b.ct==ct)
            requested_emin,requested_emax=self.fit_emin[ct],self.fit_emax[ct]
            if np.abs(actual_emin-requested_emin)>1:
                print ('Warning: For ct=%d, requested emin is %d, actual emin is %d' % (ct,requested_emin,actual_emin))
            if np.abs(actual_emax-requested_emax)>1:
                print ('Warning: For ct=%d, requested emax is %d, actual emax is %d' % (ct,requested_emax,actual_emax))

    def setup_energy_bands(self,emin=[0,0]):

        groupings = dict()
        for bc in self.bin_centers:
            groupings[bc] = [band for band in self.bands if (band.e==bc and (band.emin>=emin[band.ct]))]

        self.energy_bands = [roi_bands.ROIEnergyBand(groupings[bc]) for bc in self.bin_centers]

    def mapper(self,which):
        """ Map the argument `which' passed into functions such as localize
            into the manager that which belongs to and and index in the
            manager for the particular source referred to.

            If which is an integer, it is assumed to be in the point
            source manager, so the return is the point source manager
            and the particular index.

            If which is None, it will find the closest point+extended
            to the center of the ROI and return that.

            if which is a string, the name of all the point sources and
            diffuse sources is searched and the first source with that
            name is returned.

            If which is of type PointSource, the return is the point
            source manager and the index for the particular point source.

            If which is of type DiffuseSource, the return is the diffuse
            source manager and the index is for the particular diffuse
            source. 

            If which is a list, convert each of the items in the list
            to a manager & index and if all the managers are the same,
            return a list of indices. If there are multiple different
            managers, an exception is raised. 
            
            As a practical aside, if which is a Point/Diffuse source 
            whose integral diverges, str(which) will print a warning message.
            So it is best to get rid of these cases before trying
            to cast which to a string. """
        if type(which)==int:
            if which<len(self.psm.models):
                return self.psm,which
            else:
                return self.dsm,which-len(self.psm.models)-1
        elif which is None:
            # Get closest to ROI center.
            source=self.get_sources()[0]
            if isinstance(source,PointSource):
                return self.psm,np.where(self.psm.point_sources==source)[0][0]
            else:
                return self.dsm,np.where(self.dsm.diffuse_sources==source)[0][0]
        elif isinstance(which,PointSource):
            return self.psm,int(np.where(self.psm.point_sources==which)[0])
        elif isinstance(which,DiffuseSource):
            return self.dsm,int(np.where(self.dsm.diffuse_sources==which)[0])
        elif np.any(str(which)==self.psm.names):
            return self.psm,int(np.where(str(which)==self.psm.names)[0])
        elif np.any(str(which)==self.dsm.names):
            return self.dsm,int(np.where(str(which)==self.dsm.names)[0])
        elif isinstance(which,ROIDiffuseModel):
            return self.dsm,int(np.where(self.dsm.bgmodels==which)[0])
        elif type(which) == list:
            if len(which)<1: raise Exception("Cannot pass empty list as argument for which")
            managers,indices=zip(*[list(self.mapper(_)) for _ in which])
            if np.unique(managers).shape[0]!=1:
                raise Exception("List passed as which argument must be all point or diffuse sources.")
            return managers[0],np.asarray(indices)
        else:
            raise Exception("Unknown which argument = %s" % str(which))

    def update_counts(self, parameters=None):
        """ set parameters, if specified
            in any case, update the predicted counts in all the bands from the diffuse and point source models
        """
        if parameters is not None:
            self.set_parameters(parameters)
        self.bgm.update_counts(self.bands)
        self.psm.update_counts(self.bands)

    def logLikelihood(self, parameters, *args):
        """ the total likelihood, according to model
            parameters parameters to pass to model
        """
        if np.any(np.isnan(parameters)):
            # pretty ridiculous that this check must be made, but fitter passes NaNs...
            return 1e6
            # not sure if should "set parameters" in this case

        self.update_counts(parameters)

        ll = sum(band.logLikelihood() for band in self.bands)
        return 1e6 if np.isnan(ll) else ll

    def bandFit(self,which):
        """ Perform a spectral independendent fit of the source
            specified by which and return the negative of the log
            likelihood. This is analogous to roi.fit() for a 
            spectral independent model. Note that unlike 
            roi.logLikelihood, the positive loglikelihood is returned. """
        manager,index=self.mapper(which)

        if manager == self.dsm and not \
                isinstance(self.dsm.diffuse_sources[index],ExtendedSource):
            raise Exception("Warning, bandFit only works for ExtendedSource diffuse sources.")

        self.update_counts()

        self.setup_energy_bands()
        if manager == self.psm:
            for eb in self.energy_bands:
                eb.bandFit(which=index,saveto='bandfits')

            ll = -sum(band.bandLikelihood([band.bandfits if
                                           band.bandfits > 0 else 0],index) \
                      for band in self.bands)
            return ll
        else:
            ll = 0
            for eb in self.energy_bands:
                bfe=BandFitExtended(index,eb,self)
                bfe.fit(saveto='bandfits')

                ll -= sum(bfe.bandLikelihoodExtended([band.bandfits if band.bandfits > 0 else 0], 
                                                     band, myband) \
                          for band,myband in zip(bfe.bands,bfe.mybands))
            return ll

    def gradient(self,parameters,*args):
        """ Implement the gradient of the log likelihood wrt the model parameters."""

        bands     = self.bands
        models    = self.psm.models

        # sanity check -- for efficiency, the gradient should be called with the same params as the log likelihood
        if not np.allclose(parameters,self.parameters(),rtol=0,atol=1e-6):
            self.update_counts(parameters)

        # do the point sources
        indices  = np.arange(len(models))[np.asarray([np.any(m.free) for m in models])] if len(models)>0 else []
        nparams  = np.asarray([model.free.sum() for model in models])
        gradient = np.zeros(nparams.sum())

        for b in bands:
            cp = 0
            if b.has_pixels:
                b.pix_weights = pix_weights = b.pix_counts / (b.bg_all_pix_counts + b.ps_all_pix_counts)
            else:
                pixterm = 0

            for ind,model in zip(indices,models[indices]):
                grad    = b.gradient(model)[model.free]*b.er[ind] # correct for exposure
                npar   = nparams[ind]
                apterm = b.phase_factor*b.overlaps[ind]
                if b.has_pixels:
                    pixterm = (pix_weights*b.ps_pix_counts[:,ind]).sum()
                gradient[cp:cp+npar] += grad * (apterm - pixterm)
                cp += npar

        # add in diffuse components
        gradient  = np.append(self.bgm.gradient(bands),gradient)
        
        # Note, no need to transform gradient into log space because
        # gradient now returns the gradient with respect to internal parameters.
        return gradient

    def check_gradient(self,tol=1e-3,get_approx_gradient=False):
        """ Compute numerical gradient by finite difference and compare to
            analytic gradient.  Return True if they agree numerically.

            NB -- the numerical gradient accuracy is extremely sensitive
            to step size, so the current algorithm is a crude adaptive
            scheme.  A failure to agree numerically does not necessarily
            imply an incorrect analytic gradient.
        """
        p0 = self.get_parameters().copy()
        grad = np.empty_like(p0)
        max_iter = 40
        for i in xrange(len(p0)):
            delta = 1e-2
            prev_grad = np.inf
            for j in xrange(max_iter):
                pwork = p0.copy()
                pwork[i] += delta
                lhi = self.logLikelihood(pwork)
                pwork[i] -= 2*delta
                llo = self.logLikelihood(pwork)
                g = (lhi-llo)/(2*delta)
                if abs(g-prev_grad) < tol/2:
                    grad[i] = g
                    break
                else:
                    prev_grad = g
                    delta /= 1.5
            if j == (max_iter-1):
                print ('Adaptive scheme did not converge for parameter %d.'%i)

        if get_approx_gradient: return grad
        if np.max(np.abs(grad-self.gradient(p0))) < tol:
            return True
        print (np.max(np.abs(grad-self.gradient(p0))))
        ratio = self.gradient(p0)/grad
        print ('Gradients did not agree! Ratio of gradients:')
        print (ratio)
        return False

    def _check_model_gradients(self):
        """ Determine if all spectral models support gradient."""
        for model in np.append(self.psm.models,self.dsm.models):
            if np.any(model.free) and (not hasattr(model,'external_gradient')):
                return False
        return True
         
    def parameters(self):
        """Merge parameters from background and point sources."""
        return np.asarray(self.bgm.parameters()+self.psm.parameters())

    def get_parameters(self):
        """Support for hessian calculation in specfitter module."""
        return self.parameters()

    def get_free_errors(self):
        """Return the diagonal elements of the covariance matrix -- useful for step sizes in minimization, if known."""
        return np.asarray(self.bgm.get_free_errors() + self.psm.get_free_errors())

    def set_parameters(self,parameters):
        """Support for hessian calculation in specfitter module."""
        assert len(parameters)==len(self.psm.parameters())+len(self.bgm.parameters()), 'bad parameter length, %s!=%s+%s' % (len(parameters),len(self.psm.parameters()),len(self.bgm.parameters()))
        self.bgm.set_parameters(parameters,current_position=0)
        self.psm.set_parameters(parameters,current_position=len(self.bgm.parameters()))
        self.fit_parameters = parameters

    def fit_background(self):
        old_psm_frees = []
        for m in self.psm.models:
            old_psm_frees.append(m.free.copy())
            #m.free = np.asarray([False]*len(m.free))
            m.free[:] = False
        self.fit(fit_bg_first = False,estimate_errors=False)
        for n,nm in enumerate(self.psm.models):
            nm.free[:] = old_psm_frees[n]

    def __pre_fit__(self):

        #cache frozen values
        param_state = np.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
        param_vals  = np.concatenate([m.get_all_parameters(internal=True)  for m in self.psm.models] \
                                + [m.get_all_parameters(internal=True)  for m in self.bgm.models])

        if self.param_state is None or self.param_vals is None or \
            len(param_state)  != len(self.param_state) or \
            np.any(param_state != self.param_state) or \
            np.any(param_vals  != self.param_vals):

            self.psm.cache(self.bands)

            ### NOTA BENE
            #self.bgm.cache() # remove if don't adopt with new paradigm
            #############

            self.param_state = param_state
            self.param_vals  = param_vals

    def __update_state__(self):
        """ Helper function which should, hopefully, consistently update
            all of the model predictions for an ROI, even if some of
            sources have been deleted, zeroed, or frozen. """
        self.__pre_fit__()
        self.update_counts()

    def fit(self,method='simplex', tolerance = 0.01, save_values = True, 
                     fit_bg_first = False, estimate_errors=True, error_for_steps=False,
                     use_gradient = True, gtol = 1e-1):
        """Maximize likelihood and estimate errors.

            method     -- ['simplex'] fitter; 'powell' or 'simplex' or 'minuit'
            tolerance -- (approximate) absolute tolerance of log likelihood value
        """

        if method not in ['simplex','powell','minuit']:
            raise Exception('Unknown fitting method for F.fit(): "%s"' % method)
        if use_gradient and (not self._check_model_gradients()):
            if not self.quiet:
                print ('Found a model without a gradient method.  Switching to simplex method.')
            method = 'simplex'; use_gradient = False

        if fit_bg_first:
            self.fit_background()

        self.__pre_fit__()

        if not self.quiet: print ('.....performing likelihood maximization...',)
        if method == 'minuit':
            from uw.utilities.minuit import Minuit
            temp_params = self.parameters()
            npars = self.parameters().shape[0]
            param_names = ['p%i'%i for i in xrange(npars)]
            
            if use_gradient:
                gradient         = self.gradient
                force_gradient = 1
            else:
                gradient         = None
                force_gradient = 0

            if error_for_steps:
                steps = self.get_free_errors()
                steps[steps<1e-6] = 0.04 # for models without error estimates, put in the defaults
                steps[steps > 1]  = 1     # probably don't want to step more than 100%...
                m = Minuit(self.logLikelihood,temp_params,up=.5,maxcalls=20000,tolerance=tolerance,printMode=-self.quiet,param_names=param_names,steps=steps)
            else:
                m = Minuit(self.logLikelihood,temp_params,up=.5,maxcalls=20000,tolerance=tolerance,printMode=-self.quiet,param_names=param_names)

            params,fval = m.minimize()

            if save_values:
                if estimate_errors == True:
                     self.__set_error_minuit(m,'HESSE')
                self.logLikelihood(params) # reset values to the ones found by minimization step
                self.prev_logl = self.logl if self.logl is not None else -fval
                self.logl = -fval
            #Saving this reference seems to cause a memory leak.
            #self._minuit = m
            return -fval
        else:
            ll_0 = self.logLikelihood(self.parameters())
            if use_gradient:
                f0 = fmin_bfgs(self.logLikelihood,self.parameters(),self.gradient,full_output=1,maxiter=500,gtol=gtol,disp=0)
                for i in xrange(10):
                    f = self._save_bfgs = fmin_bfgs(self.logLikelihood,self.parameters(),self.gradient,full_output=1,maxiter=500,gtol=gtol,disp=0)
                    if abs(f0[1] - f[1]) < tolerance: break # note absolute tolerance
                    if not self.quiet:
                        print ('Did not converge on this gradient iteration.  Trying again.')
                        print (f0[1],f[1],abs(f0[1]-f[1]))
                    f0 = f

            else:
                minimizer  = fmin_powell if method == 'powell' else fmin
                f = minimizer(self.logLikelihood,self.parameters(),full_output=1,
                                  maxiter=10000,maxfun=20000,ftol=0.01/abs(ll_0), disp=0 if self.quiet else 1)
            if not self.quiet: print ('Function value at minimum: %.8g'%f[1])
            if save_values:
                self.set_parameters(f[0])
                if estimate_errors: self.__set_error__(use_gradient)
                self.prev_logl = self.logl if self.logl is not None else -f[1]
                self.logl = -f[1]

            return -f[1]

        ## check for error conditions here
        #    if not self.quiet: print ('good fit!')
        #    return -f[1]

    def __set_error__(self,use_gradient=False):

        n = len(self.bgm.parameters())
        if use_gradient:
            hessian = specfitter.mycov(self.gradient,self.parameters(),full_output=True)[1]
        else:
            hessian = specfitter.SpectralModelFitter.hessian(self,self.logLikelihood)[0] #does Hessian for free parameters
        success = False
        # TODO -- check the return code

        def _has_nan(m):
            nan = np.any(np.isnan(m))
            if nan and (not self.quiet): print ('Found NaN in covariance matrix!')
            return nan

        try:
            if not self.quiet: print ('Attempting to invert full hessian...')
            self.cov_matrix = cov_matrix = np.linalg.inv(hessian)
            if _has_nan(cov_matrix): raise ValueError
            self.bgm.set_covariance_matrix(cov_matrix,current_position=0)
            self.psm.set_covariance_matrix(cov_matrix,current_position=n)
            success = True
        except Exception:
            if len(self.psm.parameters()) > 0:
                if not self.quiet: print ('Skipping full Hessian inversion, trying point source parameter subset...')
                try:
                    self.cov_matrix = cov_matrix = np.linalg.inv(hessian[n:,n:])
                    if _has_nan(cov_matrix): raise ValueError
                    self.psm.set_covariance_matrix(cov_matrix,current_position=0)
                    success = True
                except Exception:
                    if not self.quiet: print ('Unable to recover point source errors.  Any reported error is unreliable!')
            else:
                nump = len(self.get_parameters())
                self.cov_matrix = np.zeros([nump,nump])

        return success

    def __set_error_minuit(self,m,method='HESSE'):
        """Compute errors for minuit fit."""

        #Not sure yet if there will be problems with including the backgrounds.
        self.cov_matrix = m.errors(method=method)
        self.bgm.set_covariance_matrix(self.cov_matrix,current_position = 0)
        self.psm.set_covariance_matrix(self.cov_matrix,current_position = len(self.bgm.parameters()))

    def __str__(self):
        bg_header  = '======== DIFFUSE SOURCE FITS =============='
        ps_header  = '======== POINT SOURCE FITS ============'
        if (self.logl is not None) and (self.prev_logl is not None):
            ll_string  = 'Log likelihood change: %.2f'%(self.logl - self.prev_logl)
        else:
            ll_string  = ''
        return '\n\n'.join([ps_header,self.psm.__str__(),bg_header,self.bgm.__str__(),ll_string])

    def TS(self,which=0,quick=True,bandfits=False, fit_kwargs=dict(method='simplex')):
        """Calculate the significance of the central point source.

            quick -- if set True, just calculate likelihood with source flux set to 0
                        if set False, do a full refit of all other free sources

            which -- the index of source to calculate -- default to central.

            bandfits  -- if True, calcualte the likelihood using a band-by-band (model independent) spectral fit """

        if quick:
            self.zero_ps(which)
            ll_0 = self.logLikelihood(self.get_parameters())
            self.unzero_ps(which)
            if not bandfits:
                ll_1 = self.logLikelihood(self.get_parameters())
            else:
                ll_1 = -self.bandFit(which)

            if ll_0 == 1e6 or ll_1 == 1e6: 
                 print ('Warning: loglikelihood is NaN, returning TS=0')
                 return 0
            return 2*(ll_0 - ll_1)

        save_params = self.parameters().copy() # save free parameters
        self.zero_ps(which)
        self.fit(save_values=False, **fit_kwargs)
        ll_0 = -self.logLikelihood(self.parameters())

        if not self.quiet: print (self)
        self.unzero_ps(which)
        self.set_parameters(save_params) # reset free parameters
        self.__update_state__() # restore caching
        if not bandfits:
            ll = -self.logLikelihood(save_params)
        else:
            ll = self.bandFit(which)
        if ll_0 == 1e6 or ll == 1e6: 
             print ('Warning: loglikelihood is NaN, returning TS=0')
             return 0
        return 2*(ll - ll_0)

    def localize(self,which=0, tolerance=1e-3,update=False, verbose=False, bandfits=False, seedpos=None, **kwargs):
        """Localize a source using an elliptic approximation to the likelihood surface.

            which      -- index of point source; default to central
                             ***if localizing non-central, ensure ROI is large enough!***
            tolerance -- maximum difference in degrees between two successive best fit positions
            update     -- if True, update localization internally, i.e., recalculate point source contribution
            bandfits  -- if True, use a band-by-band (model independent) spectral fit; otherwise, use broadband fit
            seedpos    -- if set, use this position instead of the source position

            return fit position
        """
        rl = roi_localize.localizer(self, which=which, bandfits=bandfits,
                                   tolerance=tolerance, update=update, verbose=verbose, **kwargs)
        if seedpos is not None:
            rl.sd = seedpos  # override 
        return rl.localize()

    @decorate_with(roi_localize.DualLocalizer,append_init=True)
    def dual_localize(self,*args,**kwargs):
        return roi_localize.DualLocalizer(self,*args,**kwargs).localize()
    
    @decorate_with(ROIExtendedModel.fit_extension)
    def fit_extension(self,which,*args,**kwargs):

        manager,index=self.mapper(which)

        if manager != self.dsm or not isinstance(self.dsm.diffuse_sources[index],ExtendedSource):
            raise Exception("Can only fit the extension of extended sources.")

        return self.dsm.bgmodels[index].fit_extension(self,*args,**kwargs)

    @decorate_with(ROIExtendedModel.TS_ext)
    def TS_ext(self,which,*args,**kwargs):

        manager,index=self.mapper(which)

        if manager != self.dsm or not isinstance(self.dsm.diffuse_sources[index],ExtendedSource):
            raise Exception("Can only calculate TS_ext of extended sources.")

        return self.dsm.bgmodels[index].TS_ext(self,*args,**kwargs)

    def upper_limit(self, which, *args, **kwargs):
        f=roi_upper_limits.FluxUpperLimit(self, which, *args, **kwargs)
        return f.get_limit()

    upper_limit_quick = roi_upper_limits.upper_limit_quick
    def extension_upper_limit(self, *args, **kwargs):
        e=roi_upper_limits.ExtensionUpperLimit(self, *args, **kwargs)
        return e.results()

    def add_source(self,source,**kwargs):
         """Add a new source object to the model.

            N.B. for point sources, add a pointspec_helpers.PointSource
            object. For diffuse sources, add either a
            roi_diffuse.DiffuseSource or an ROIDiffuseModel object.
            If model is an roi_diffuse.DiffuseSource object, the
            pointspec_helpers.get_default_diffuse_mapper function is used
            to convert it to an roi_diffuse.ROIDiffuseModel object. """
         if isinstance(source,PointSource):
             manager=self.psm
         elif isinstance(source,DiffuseSource) or isinstance(source,ROIDiffuseModel):
             manager=self.dsm

             # convert DiffuseSource -> ROIDiffuseModel object
             if isinstance(source,DiffuseSource):
                 diffuse_mapper = get_default_diffuse_mapper(self.sa,self.roi_dir,self.quiet)
                 source=diffuse_mapper(source)
         else:
             raise Exception("Unable to add source %s. Only able to add PointSource, DiffuseSource, or ROIDiffuseModel objects.")
         if self.__dict__.has_key('cov_matrix'): del self.cov_matrix
         manager.add_source(source,self.bands,**kwargs)
         self.__update_state__()

    def del_source(self,which):
         """Remove the source at position given by which from the model."""
         manager,index=self.mapper(which)
         if self.__dict__.has_key('cov_matrix'): del self.cov_matrix
         source=manager.del_source(index,self.bands)
         self.__update_state__()
         return source

    def zero_source(self,which):
         """Set the flux of the source given by which to 0."""
         manager,index=self.mapper(which)
         source=manager.zero_source(index,self.bands)
         self.__update_state__()
         return source

    def unzero_source(self,which):
         """Restore a previously-zeroed flux."""
         manager,index=self.mapper(which)
         manager.unzero_source(index,self.bands)
         self.__update_state__()

    # for backwards compatability, clone functions
    add_ps = add_source
    del_ps = del_source
    zero_ps = zero_source
    unzero_ps = unzero_source

    # Get the modify functions from roi_modify
    modify_loc = roi_modify.modify_loc
    modify_spatial_model = roi_modify.modify_spatial_model
    modify_model = roi_modify.modify_model
    modify_name = roi_modify.modify_name
    modify = roi_modify.modify
        
    # get the print_summary function from roi_printing
    print_summary=roi_printing.print_summary
    print_resids=roi_printing.print_resids
    printSpectrum=roi_printing.printSpectrum

    # localization results
    print_ellipse=roi_localize.print_ellipse
    get_ellipse=roi_localize.get_ellipse
 
    # get the toXML, toRegion, and toResults function from xml_parsers
    toXML=xml_parsers.writeROI
    toRegion=region_writer.writeRegion
    toResults=results_writer.writeResults
    
    def get_model(self,which):
        """ return a reference to the model
            which : integer or string """
        manager, index = self.mapper(which) #raise exception if wrong.
        return manager.models[index]
        
    def get_source(self, which):
        """ return a reference to a source in the ROI by name, or point-source index"""
        manager, index = self.mapper(which) #raise exception if wrong.
        return manager.point_sources[index] if manager==self.psm else self.dsm.diffuse_sources[index] 
    
    def get_sources(self):
        """ Returns all localizable sources in the ROI sorted by skydir. """
        sources=self.psm.point_sources.tolist()+self.get_extended_sources()
        sources.sort(key=lambda s:s.skydir.difference(self.roi_dir))
        return sources

    def get_extended_sources(self):
        sources=[i for i in self.dsm.diffuse_sources.tolist() if hasattr(i,'skydir')]
        sources.sort(key=lambda s:s.skydir.difference(self.roi_dir))
        return sources

    def get_names(self):
        return [i.name for i in self.get_sources()]

    # get these functions from roi_save.py
    save=roi_save.save
    load=staticmethod(roi_save.load)

    @decorate_with(roi_image.ROITSMapImage,append_init=True)
    def tsmap(self,filename,**kwargs):
        i=roi_image.ROITSMapImage(self,**kwargs)
        i.get_pyfits().writeto(filename,clobber=True)
        return i

    @decorate_with(sed_plotter.plot_sed)
    def plot_sed(self,which=None,filename=None,**kwargs):
        if 'outdir' in kwargs: 
            raise Exception('Use filename, not outdir')
        return sed_plotter.plot_sed(self,which=which,outdir=filename,**kwargs)

    @decorate_with(mapplots.ROIDisplay,append_init=True)
    def plot_counts_map(self,filename=None,**kwargs):
        i=mapplots.ROIDisplay(self,**kwargs)
        i.show(filename=filename)
        return i

    @decorate_with(counts_plotter.roi_pipeline_counts_plot)
    def plot_counts_spectra(self,filename,**kwargs):
        counts_plotter.roi_pipeline_counts_plot(self,counts_dir=filename,**kwargs)

    @decorate_with(roi_plotting.ROISlice,append_init=True)
    def plot_slice(self,which=None,filename=None,datafile=None,**kwargs):
        i=roi_plotting.ROISlice(self,which=which,**kwargs)
        i.show(filename=filename,datafile=datafile)
        return i

    @decorate_with(roi_plotting.ROIRadialIntegral,append_init=True)
    def plot_radial_integral(self,which=None,filename=None,datafile=None,**kwargs):
        i=roi_plotting.ROIRadialIntegral(self,which=which,**kwargs)
        i.show(filename=filename,datafile=datafile)
        return i

    @decorate_with(mapplots.ROISmoothedSource,append_init=True)
    def plot_source(self,which=None,filename=None,**kwargs):
        i=mapplots.ROISmoothedSource(self,which=which,**kwargs)
        i.show(filename=filename)
        return i

    @decorate_with(mapplots.ROISmoothedSources,append_init=True)
    def plot_sources(self,which=None,filename=None,axes=None,**kwargs):
        i=mapplots.ROISmoothedSources(self,which=which,**kwargs)
        i.show(filename=filename,axes=axes)
        return i

    @decorate_with(mapplots.ROISmoothedResidual,append_init=True)
    def plot_residual(self,filename=None,axes=None,**kwargs):
        i=mapplots.ROISmoothedResidual(self,**kwargs)
        i.show(filename=filename,axes=axes)
        return i

    @decorate_with(mapplots.ROISignificance,append_init=True)
    def plot_significance(self,filename=None,axes=None,**kwargs):
        i=mapplots.ROISignificance(self,**kwargs)
        i.show(filename=filename,axes=axes)
        return i

    @decorate_with(mapplots.ROISmoothedBeforeAfter,append_init=True)
    def plot_before_after(self,filename=None,**kwargs):
        i=mapplots.ROISmoothedBeforeAfter(self,**kwargs)
        i.show(filename=filename)
        return i

    @decorate_with(mapplots.ROITSMapPlotter,append_init=True)
    def plot_tsmap(self,filename=None,axes=None,**kwargs):
        i=mapplots.ROITSMapPlotter(self,**kwargs)
        i.show(filename=filename, axes=axes)
        return i

    @decorate_with(mapplots.ROISmoothedDataModel,append_init=True)
    def plot_model(self,filename="model_counts.png",**kwargs):
        i=mapplots.ROISmoothedDataModel(self,**kwargs)
        i.show(filename=filename)
        return i

    def change_binning(self,fit_emin,fit_emax):
        """ This function recreates the ROI using a new energy binning.  
            Kind of inefficient, but easy. """
        kwargs=keyword_options.defaults_to_kwargs(self,ROIAnalysis)
        kwargs['fit_emin']=fit_emin
        kwargs['fit_emax']=fit_emax
        self.__init__(roi_dir=self.roi_dir,
                      ps_manager=self.psm,
                      ds_manager=self.dsm,
                      spectral_analysis=self.sa,
                      **kwargs)


load=ROIAnalysis.load
