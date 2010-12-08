"""
Module implements a binned maximum likelihood analysis with a flexible, energy-dependent ROI based
    on the PSF.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_analysis.py,v 1.55 2010/11/27 22:34:47 lande Exp $

author: Matthew Kerr
"""

import numpy as N
import math, pickle, collections

from uw.like import roi_bands, roi_localize , specfitter

from pointspec_helpers import PointSource,get_default_diffuse_mapper
from roi_diffuse import ROIDiffuseModel,DiffuseSource
from roi_extended import ExtendedSource,BandFitExtended
from uw.utilities import keyword_options

from scipy.optimize import fmin,fmin_powell,fmin_bfgs
from scipy.stats.distributions import chi2
from scipy import integrate
from numpy.linalg import inv

EULER_CONST  = N.exp(1)
LOG_JACOBIAN = 1./N.log10(EULER_CONST)

###====================================================================================================###

class ROIAnalysis(object):

    defaults = (
        ("fit_emin",[125,125],"a two-element list giving front/back minimum energy. Independent energy ranges for front and back"),
        ("fit_emax",[1e5 + 100,1e5 + 100],"A two-element list giving front/back maximum energy. 0th position for event class 0."),
        ("quiet",False,'Set True to suppress (some) output'),
        ("catalog_aperture",-1,"Pulsar catalog analysis only -- deprecate"),
        ("phase_factor",1.,"Pulsar phase. A correction that can be made if the data has undergone a phase selection -- between 0 and 1"),
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
        self.bin_centers = N.sort(list(set([b.e for b in self.bands])))
        self.bin_edges    = N.sort(list(set([b.emin for b in self.bands] + [b.emax for b in self.bands])))

        self.param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
        self.param_vals  = N.concatenate([m.p for m in self.psm.models] + [m.p for m in self.bgm.models])
        self.psm.cache(self.bands)
        #self.psm.update_counts(self.bands)

        self.logl = -self.logLikelihood(self.get_parameters()) # make sure everything initialized

    def __setup_bands__(self):

        try:
            self.fit_emin = [float(self.fit_emin)]*2
        except:
            pass

        try:
            self.fit_emax = [float(self.fit_emax)]*2
        except:
            pass

        self.bands = collections.deque()
        for band in self.sa.pixeldata.dmap:
            evcl = band.event_class() & 1 # protect high bits

            if (band.emin() + 1) >= self.fit_emin[evcl] and (band.emax() - 1) < self.fit_emax[evcl]:
                self.bands.append(roi_bands.ROIBand(band,self.sa,self.roi_dir,catalog_aperture=self.catalog_aperture))

        self.bands = N.asarray(self.bands)

        self.psm.setup_initial_counts(self.bands)
        self.dsm.setup_initial_counts(self.bands)

        for band in self.bands: band.phase_factor = self.phase_factor

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
            return self.psm,which
        elif isinstance(which,PointSource):
            return self.psm,int(N.where(self.psm.point_sources==which)[0])
        elif isinstance(which,DiffuseSource):
            return self.dsm,int(N.where(self.dsm.diffuse_sources==which)[0])
        elif N.any(str(which)==self.psm.names):
            return self.psm,int(N.where(str(which)==self.psm.names)[0])
        elif N.any(str(which)==self.dsm.names):
            return self.dsm,int(N.where(str(which)==self.dsm.names)[0])
        elif isinstance(which,ROIDiffuseModel):
            return self.dsm,int(N.where(self.dsm.bgmodels==which)[0])
        elif type(which) == list:
            if len(which)<1: raise Exception("Cannot pass empty list as argument for which")
            managers,indices=zip(*[list(self.mapper(_)) for _ in which])
            if N.unique(managers).shape[0]!=1:
                raise Exception("List passed as which argument must be all point or diffuse sources.")
            return managers[0],N.asarray(indices)
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

    def logLikelihood(self,parameters,*args):
        """ the total likelihood, according to model
            parameters parameters to pass to model
        """
        if N.any(N.isnan(parameters)):
            # pretty ridiculous that this check must be made, but fitter passes NaNs...
            return 1e6
            # not sure if should "set parameters" in this case

        self.update_counts(parameters)

        ll = sum(band.logLikelihood(phase_factor=self.phase_factor) for band in self.bands)
        return 1e6 if N.isnan(ll) else ll

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
        pf         = self.phase_factor  # can just multiply the aperture term
        models    = self.psm.models

        # sanity check -- for efficiency, the gradient should be called with the same params as the log likelihood
        if not N.allclose(parameters,self.parameters(),rtol=0,atol=1e-6):
            self.update_counts(parameters)

        # do the point sources
        indices  = N.arange(len(models))[N.asarray([N.any(m.free) for m in models])]
        nparams  = N.asarray([model.free.sum() for model in models])
        gradient = N.zeros(nparams.sum())

        for b in bands:
            cp = 0
            if b.has_pixels:
                b.pix_weights = pix_weights = b.pix_counts / (b.bg_all_pix_counts + b.ps_all_pix_counts)
            else:
                pixterm = 0

            for ind,model in zip(indices,models[indices]):
                grad    = b.gradient(model)[model.free]*b.er[ind] # correct for exposure
                np      = nparams[ind]
                apterm = pf*b.overlaps[ind]
                if b.has_pixels:
                    pixterm = (pix_weights*b.ps_pix_counts[:,ind]).sum()
                gradient[cp:cp+np] += grad * (apterm - pixterm)
                cp += np

        # add in diffuse components
        gradient  = N.append(self.bgm.gradient(bands,pf),gradient)
        
        # transform into log space and return
        return gradient * 10**parameters * LOG_JACOBIAN
         

    def parameters(self):
        """Merge parameters from background and point sources."""
        return N.asarray(self.bgm.parameters()+self.psm.parameters())

    def get_parameters(self):
        """Support for hessian calculation in specfitter module."""
        return self.parameters()

    def get_free_errors(self):
        """Return the diagonal elements of the covariance matrix -- useful for step sizes in minimization, if known."""
        return N.asarray(self.bgm.get_free_errors() + self.psm.get_free_errors())

    def set_parameters(self,parameters):
        """Support for hessian calculation in specfitter module."""
        self.bgm.set_parameters(parameters,current_position=0)
        self.psm.set_parameters(parameters,current_position=len(self.bgm.parameters()))
        self.fit_parameters = parameters

    def fit_background(self):
        old_psm_frees = []
        for m in self.psm.models:
            old_psm_frees.append(m.free.copy())
            #m.free = N.asarray([False]*len(m.free))
            m.free[:] = False
        self.fit(fit_bg_first = False,estimate_errors=False)
        for n,nm in enumerate(self.psm.models):
            nm.free[:] = old_psm_frees[n]

    def __pre_fit__(self):

        #cache frozen values
        param_state = N.concatenate([m.free for m in self.psm.models] + [m.free for m in self.bgm.models])
        param_vals  = N.concatenate([m.p for m in self.psm.models] + [m.p for m in self.bgm.models])

        if len(param_state)  != len(self.param_state) or \
            N.any(param_state != self.param_state) or \
            N.any(param_vals  != self.param_vals):

            self.psm.cache(self.bands)

            ### NOTA BENE
            #self.bgm.cache() # remove if don't adopt with new paradigm
            #############

            self.param_state = param_state
            self.param_vals  = param_vals

    def fit(self,method='simplex', tolerance = 0.01, save_values = True, do_background=True,
                     fit_bg_first = False, estimate_errors=True, error_for_steps=False,
                     use_gradient = False, gtol = 1e-1):
        """Maximize likelihood and estimate errors.

            method     -- ['simplex'] fitter; 'powell' or 'simplex' or 'minuit'
            tolerance -- (approximate) absolute tolerance of log likelihood value
        """

        if method not in ['simplex','powell','minuit']:
            raise Exception('Unknown fitting method for F.fit(): "%s"' % method)

        if fit_bg_first:
            self.fit_background()

        self.__pre_fit__()

        if not self.quiet: print '.....performing likelihood maximization...',
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
                        print 'Did not converge on first gradient iteration.  Trying again.'
                        print f0[1],f[1],abs(f0[1]-f[1])
                    f0 = f

            else:
                minimizer  = fmin_powell if method == 'powell' else fmin
                f = minimizer(self.logLikelihood,self.parameters(),full_output=1,
                                  maxiter=10000,maxfun=20000,ftol=0.01/abs(ll_0), disp=0 if self.quiet else 1)
            if not self.quiet: print 'Function value at minimum: %.8g'%f[1]
            if save_values:
                self.set_parameters(f[0])
                if estimate_errors: self.__set_error__(do_background,use_gradient)
                self.prev_logl = self.logl if self.logl is not None else -f[1]
                self.logl = -f[1]

            return -f[1]

        ## check for error conditions here
        #    if not self.quiet: print 'good fit!'
        #    return -f[1]

    def __set_error__(self,do_background=True,use_gradient=False):

        n = len(self.bgm.parameters())
        if use_gradient:
            hessian = specfitter.mycov(self.gradient,self.parameters(),full_output=True)[1]
        else:
            hessian = specfitter.SpectralModelFitter.hessian(self,self.logLikelihood)[0] #does Hessian for free parameters
        success = False
        # TODO -- check the return code

        try:
            if not do_background: raise Exception # what is this about?
            if not self.quiet: print 'Attempting to invert full hessian...'
            self.cov_matrix = cov_matrix = inv(hessian)
            if N.any(N.isnan(cov_matrix)):
                if not self.quiet: print 'Found NaN in covariance matrix!'
                raise Exception
            self.bgm.set_covariance_matrix(cov_matrix,current_position=0)
            self.psm.set_covariance_matrix(cov_matrix,current_position=n)
            success = True
        except:
            if len(self.psm.parameters()) > 0:
                if not self.quiet: print 'Skipping full Hessian inversion, trying point source parameter subset...'
                try:
                    self.cov_matrix = cov_matrix = inv(hessian[n:,n:])
                    if N.any(N.isnan(cov_matrix)):
                        if not self.quiet: print 'Found NaN in covariance matrix!'
                        raise Exception
                    self.psm.set_covariance_matrix(cov_matrix,current_position=0)
                    success = True
                except:
                    if not self.quiet: print 'Error in calculating and inverting hessian.'
            else:
                np = len(self.get_parameters())
                self.cov_matrix = N.zeros([np,np])

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

    def TS(self,quick=True,which=0,method='simplex', bandfits=False):
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
                 print 'Warning: loglikelihood is NaN, returning TS=0'
                 return 0
            return 2*(ll_0 - ll_1)

        save_params = self.parameters().copy() # save free parameters
        self.zero_ps(which)
        ll_0 = self.fit(save_values = False,method=method)
        if not self.quiet: print self
        self.unzero_ps(which)
        self.set_parameters(save_params) # reset free parameters
        self.__pre_fit__() # restore caching
        if not bandfits:
            ll = -self.logLikelihood(save_params)
        else:
            ll = self.bandFit(which)
        if ll_0 == 1e6 or ll == 1e6: 
             print 'Warning: loglikelihood is NaN, returning TS=0'
             return 0
        return -2*(ll_0 - ll)

    def localize(self,which=0, tolerance=1e-3,update=False, verbose=False, bandfits=False, seedpos=None):
        """Localize a source using an elliptic approximation to the likelihood surface.

            which      -- index of point source; default to central
                             ***if localizing non-central, ensure ROI is large enough!***
            tolerance -- maximum difference in degrees between two successive best fit positions
            update     -- if True, update localization internally, i.e., recalculate point source contribution
            bandfits  -- if True, use a band-by-band (model independent) spectral fit; otherwise, use broadband fit
            seedpos    -- if set, use this position instead of the source position

            return fit position
        """

        manager,index = self.mapper(which)

        if manager == self.dsm and not isinstance(manager.diffuse_sources[index],ExtendedSource):
            raise Exception("Can only localize Point and Extended Sources")

        f = roi_localize.ROILocalizer if manager == self.psm else roi_localize.ROILocalizerExtended
        rl = f(self,which=index,bandfits=bandfits,tolerance=tolerance,update=update,verbose=verbose)

        if seedpos is not None:
            rl.sd = seedpos  # override 
        return rl.localize()
    
    def fit_extension(self,which,*args,**kwargs):
        """ See roi_extended.ROIExtendedModelAnalytic's fit_extension function for parameters. """

        manager,index=self.mapper(which)

        if manager != self.dsm or not isinstance(self.dsm.diffuse_sources[index],ExtendedSource):
            raise Exception("Can only fit the extension of extended sources.")

        return self.dsm.bgmodels[index].fit_extension(self,*args,**kwargs)

    def TS_ext(self,which,*args,**kwargs):

        manager,index=self.mapper(which)

        if manager != self.dsm or not isinstance(self.dsm.diffuse_sources[index],ExtendedSource):
            raise Exception("Can only calculate TS_ext of extended sources.")

        return self.dsm.bgmodels[index].TS_ext(self,*args,**kwargs)

    def upper_limit(self,**kwargs):
        """Compute an upper limit on the source flux, by the "PDG Method"

        This method computes an upper limit on the flux of a specified source
        for a given time interval. The limit is computed by integrating the
        likelihood over the flux of the source, via Simpson's Rule, up to the
        desired percentile (confidence level). As such, it is essentially a
        Bayesian credible interval, using a uniform prior on the flux
        parameter.

        Note that the default integral limits are determined assuming that
        the relevant parameter is the normalization parameter of a PowerLaw
        model. For other models, especially PowerLawFlux, the limits should
        be specified appropriately.

        The limit returned is the integrated flux above 100 MeV in photons
        per square centimeter per second.

        Arguments:
            which: integer [0]
                Index of the point source for which to compute the limit.
            confidence: float [.95]
                Desired confidence level of the upper limit.
            integral_min: float [-15]
                Lower limit of the likelihood integral *in log space*.
            integral_max: float [-8]
                Upper limit of the likelihood integral *in log space*.
            simps_points: int [10]
                Number of evaluation points *per decade* for Simpson's rule.
            e_weight: float [0]
                Energy weight for the flux integral (see documentation for uw.like.Models)
            cgs: bool [False]
                If true return flux in cgs units (see documentation for uw.like.Models)
        """
        kw = dict(which=0,
                  confidence=0.95,
                  integral_min=-15,
                  integral_max =-8,
                  simps_points = 100,
                  e_weight = 0,
                  cgs = False)
        for k,v in kw.items():
            kw[k] = kwargs.pop(k,v)
        if kwargs:
            for k in kwargs.keys():
                print("Invalid keyword argument for ROIAnalysis.upper_limit: %s"%k)
        params = self.parameters().copy()
        ll_0 = self.logLikelihood(self.parameters())
        def like(norm):
            self.psm.models[kw['which']].p[0] = norm
            return N.exp(ll_0-self.logLikelihood(self.parameters()))
        npoints = kw['simps_points'] * (kw['integral_max'] - kw['integral_min'])
        points = N.log10(N.logspace(kw['integral_min'],
                                      kw['integral_max'],npoints*2+1))
        y = N.array([like(x)*10**x for x in points])
        trapz1 = integrate.cumtrapz(y[::2])
        trapz2 = integrate.cumtrapz(y)[::2]
        cumsimps = (4*trapz2 - trapz1)/3.
        cumsimps /= cumsimps[-1]
        i1 = N.where(cumsimps<.95)[0][-1]
        i2 = N.where(cumsimps>.95)[0][0]
        x1, x2 = points[::2][i1], points[::2][i2]
        y1, y2 = cumsimps[i1], cumsimps[i2]
        #Linear interpolation should be good enough at this point
        limit = x1 + ((x2-x1)/(y2-y1))*(kw['confidence']-y1)
        self.psm.models[0].p[0] = limit
        uflux = self.psm.models[0].i_flux(e_weight=kw['e_weight'],cgs=kw['cgs'])
        self.logLikelihood(params)
        return uflux

    def upper_limit_quick(self,which = 0,confidence = .95,e_weight = 0,cgs = False):
        """Compute an upper limit on the flux of a source assuming a gaussian likelihood.

        Arguments:
            which: integer [0]
                Index of the point source for which to compute the limit.
            confidence: float [.95]
                Desired confidence level of the upper limit.
            e_weight: float [0]
                Energy weight for the flux integral (see documentation for uw.like.Models)
            cgs: bool [False]
                If true return flux in cgs units (see documentation for uw.like.Models)

        The flux returned is an upper limit on the integral flux for the model
        above 100 MeV.

        The upper limit is found based on the change in the log likelihood
        from the maximum, using the Gaussian approximation. It is quicker
        than the Bayesian method performed by upper_limit, but less robust.
        In particular, fmin will sometimes get lost and return absurdly
        small values, and if the maximum is too far above zero, it will
        often find the lower of the two solutions to the appropriate
        equation.
        """

        delta_logl = chi2.ppf(2*confidence-1,1)/2.
        params = self.parameters().copy()
        #self.psm.models[which].p[0]  = -20
        zp = self.logLikelihood(self.parameters())

        def f(norm):
            self.psm.models[which].p[0] = N.log10(norm)
            ll = self.logLikelihood(self.parameters())
            return abs(ll - zp - delta_logl)

        limit = fmin(f,N.array([10**-6]),disp=0)[0]
        self.psm.models[which].p[0] = N.log10(limit)
        uflux = self.psm.models[which].i_flux(e_weight = e_weight,cgs = cgs)
        self.set_parameters(params)
        return uflux

    def printSpectrum(self,sources=None):
        """Print total counts and estimated signal in each band for a list of sources.

        Sources can be specified as PointSource objects, source names, or integers
        to be interpreted as indices for the list of point sources in the roi. If
        only one source is desired, it needn't be specified as a list. If no sources
        are specified, all sources with free fit parameters will be used."""
        if sources is None:
            sources = [s for s in self.psm.point_sources if N.any(s.model.free)]
        elif type(sources) != type([]):
            sources = [sources]
        if sources == []: return # No point sources in ROI
        bad_sources = []
        for i,s in enumerate(sources):
            if type(s) == PointSource:
                if not s in self.psm.point_sources:
                    print 'Source not found in source list:\n%s\n'%s
                    bad_sources += [s]
            elif type(s) == int:
                try:
                    sources[i] = self.psm.point_sources[s]
                except IndexError:
                    print 'No source #%i. Only %i source(s) specified.'\
                            %(s,len(self.psm.point_sources))
                    bad_sources += [s]
            elif type(s) == type(''):
                names = [ps.name for ps in self.psm.point_sources]
                try:
                    sources[i] = self.psm.point_sources[names.index(s)]
                except ValueError:
                    print 'No source named %s'%s
                    bad_sources += [s]
            else:
                print 'Unrecognized source specification:', s
                bad_sources += [s]
        sources = set([s for s in sources if not s in bad_sources])
        indices = [list(self.psm.point_sources).index(s) for s in sources]
        self.setup_energy_bands()

        fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic']\
                     +[' '*15+'Signal']*len(sources)
        outstring = 'Spectra of sources in ROI about %s at ra = %.2f, dec = %.2f\n'\
                          %(self.psm.point_sources[0].name, self.roi_dir.ra(), self.roi_dir.dec())
        outstring += ' '*54+'  '.join(['%21s'%s.name for s in sources])+'\n'
        outstring += '  '.join(fields)+'\n'
        print outstring
        for eb in self.energy_bands:
            print eb.spectralString(which=indices)


    def save_fit(self,outfile,additional_data=None):
        """Save the spectral models (and locations) for all point sources and diffuse models.

            This saves the need to refit.  A future iteration should actually save all of the
            pixel predictions to avoid lengthy recalculation, too.

            additional_data: an optional dictionary with keys to add to output; note that
                                  all entries should be serializable!"""

        d = collections.defaultdict(list)
        for ps in self.psm.point_sources:
            m = ps.model
            m.ra  = ps.skydir.ra()
            m.dec = ps.skydir.dec()
            m.source_name = ps.name
            d['point_sources'].append(m)
        for bg in self.bgm.models:
            d['backgrounds'].append(bg)
        try:
            d['localization'] = [self.ldir.ra(),self.ldir.dec(),self.lsigma,self.qform.par]
        except:
            print 'No localization to save.'
        if additional_data is not None:
            try:     d.update(additional_data)
            except: print 'Warning! Could not merge requested keys into output dictionary.'
        f = open(outfile,'w')
        pickle.dump(d,f)
        f.close()

    def __call__(self,v):

        pass #make this a TS map? negative -- spatialLikelihood does it, essentially

    def add_source(self,source):
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
                 diffuse_mapper = get_default_diffuse_mapper(self.sa,self.roi_dir)
                 source=diffuse_mapper(source)
         else:
             raise Exception("Unable to add source %s. Only able to add PointSource, DiffuseSource, or ROIDiffuseModel objects.")
         manager.add_source(source,self.bands)

    def del_source(self,which):
         """Remove the source at position given by which from the model."""
         manager,index=self.mapper(which)
         return manager.del_source(index,self.bands)

    def zero_source(self,which):
         """Set the flux of the source given by which to 0."""
         manager,index=self.mapper(which)
         return manager.zero_source(index,self.bands)

    def unzero_source(self,which):
         """Restore a previously-zeroed flux."""
         manager,index=self.mapper(which)
         manager.unzero_source(index,self.bands)

    # for backwards compatability, clone functions
    add_ps = add_source
    del_ps = del_source
    zero_ps = zero_source
    unzero_ps = unzero_source

    def modify_loc(self,skydir,which):
        """Move point source given by which to new location given by skydir."""
        manager,index=self.mapper(which)
        if manager==self.psm:
            rl = roi_localize.ROILocalizer(self,which=index,update=True,bandfits=False)
            rl.spatialLikelihood(skydir,update=True)
        elif manager==self.dsm:
            if isinstance(self.dsm.diffuse_sources[index],ExtendedSource):
                self.dsm.bgmodels[index].modify_loc(self.bands,skydir)
            else:
                raise Exception("Unable to modify_loc of diffuse source %s" % which)
        
    def print_summary(self, sdir=None, galactic=False, maxdist=5, title=''):
        """ formatted table point sources positions and parameter in the ROI
        Parameters
        ----------
             sdir : SkyDir, optional, default None for center
             galactic : bool, optional, default False
                 set true for l,b
             maxdist : float, optional, default 5
                 radius in degrees

        """
        if sdir is None: sdir = self.roi_dir
        print '\n\t Nearby sources within %.1f degrees %s' % (maxdist,title)
        colstring = 'name dist ra dec flux8 index cutoff'
        if galactic: colstring =colstring.replace('ra dec', 'l b')
        colnames = tuple(colstring.split())
        n = len(colnames)-1
        print ('%-20s'+n*'%10s')% colnames
        for ps in self.psm.point_sources:
            dist=math.degrees(sdir.difference(ps.skydir))
            if maxdist and dist>maxdist:  continue
            loc = (ps.skydir.l(),ps.skydir.b()) if galactic else (ps.skydir.ra(),ps.skydir.dec())
            par= ps.model.p
            npar = len(par)
            fmt = '%-20s'+3*'%10.3f'+' %9.2f%1s'
            freeflag = [ '*' if f else ' ' for f in ps.model.free]
            values = (ps.name, dist) +loc+( ps.model.fast_iflux()/1e-8, freeflag[0], )
            for i in range(1,npar): # parameters beyond flux
                fmt     += '%9.2f%1s'
                values += (10**par[i], freeflag[i]) 
            print fmt % values

    def print_resids(self):
        """Print out (weighted) residuals for each energy range, both in
           separate front/back columns and in a joint column.

           Useful for identifying systematic effects that distinguish between
           front and back events.
        """

        d = dict()
        for b in self.bands:
            key = (-1 if b.ct==1 else 1)*int(b.e)
            d[key] = b
        ens = N.sort(list(set([b.e for b in self.bands]))).astype(int)
        print ''
        print '        \t-------CT=0--------      -------CT=1--------      ------CT=0+1-------'
        print 'Energy\tMod      Obs      Res      Mod      Obs      Res      Mod      Obs      Res'
        print '        \t-------------------      -------------------      -------------------'
        for en in ens:
            s1 = '%-6.0f'%(en)
            tm = 0; to = 0
            for key in [en,-en]:
                if key in d.keys():
                    b  = d[key]
                    m = b.ps_all_counts + b.bg_all_counts
                    o = b.photons
                else:
                    m = o = 0
                wres = (o-m)/m**0.5 if m>0 else 0
                s1 = '\t'.join([s1,'%-6.0f\t%-6d\t%.1f'%(m,o,wres)])
            s1 = '\t'.join([s1,'%-6.0f\t%-6d\t%.1f'%(tm,to,(to-tm)/tm**0.5)])
            print s1
 
    def print_ellipse(self, label=True, line=True):
        """ print the ellipical parameters (all deg units):
                ra, dec
                a, b  : major and minor 1-sigma axes
                ang   : ellipse orientation, E of N
                qual  : measure of fit quality.
        Optional parameters:
            label [True] print a label line
            line  [True] print the line corresponding to the current fit
                  (only knows about one at a time)
        """
        if not self.qform: return
        labels = 'ra dec a b ang qual'.split()
        if label: print (len(labels)*'%10s') % tuple(labels)
        if not line: return
        p = self.qform.par[0:2]+self.qform.par[3:]
        print len(p)*'%10.4f' % tuple(p)
 
    def toXML(self,filename,**kwargs):
        """Write out a gtlike-style XML file."""
        from uw.utilities.xml_parsers import writeROI
        writeROI(self,filename,**kwargs)

    def toRegion(self,filename,**kwargs):
        """Write out a ds9 region file."""
        from uw.utilities.region_writer import writeRegion
        writeRegion(self,filename,**kwargs)

    def toResults(self,filename,**kwargs):
        """Write out a gtlike style results file."""
        from uw.utilities.results_writer import writeResults
        writeResults(self,filename,**kwargs)
    
