"""
Module to calculate flux and extension upper limits.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_upper_limits.py,v 1.29 2012/11/07 01:22:58 lande Exp $

author:  Eric Wallace <ewallace@uw.edu>, Joshua Lande <joshualande@gmail.com>
"""
import math
import numpy as np
from scipy import integrate
from scipy.stats.distributions import chi2
from scipy.optimize import fmin, fminbound, brentq

from uw.utilities import keyword_options
from uw.utilities.quantile import Quantile
from uw.utilities.parmap import LinearMapper,LimitMapper

from uw.like.roi_state import PointlikeState
from uw.like.roi_extended import ExtendedSource
from uw.like.SpatialModels import Disk
from uw.like.roi_state import PointlikeState

class FluxUpperLimit(object):
    """Compute an upper limit on the source flux, by the "PDG Method"

    This method computes an upper limit on the flux of a specified source
    for a given time interval. The limit is computed by integrating the
    likelihood over the flux of the source, via Simpson's Rule, up to the
    desired percentile (confidence level). As such, it is essentially a
    Bayesian credible interval, using a uniform prior on the flux
    parameter.

    Note that the default integral limits are taken from the spectral
    model's default parameter limits. See uw.like.Models for your
    particular model.  These should be a suitable integration range for
    any parameter. In the situation where the normalization parameter
    has a limit and those limits go outside the default limits, the
    default parameter limits are expanded by the specified parameter
    limits.  If the normalization parameter is unlimited and the current
    normalization value is outside the default parameter limits, the
    integration range will be expanded to include the current value.

    The limit returned is the integrated flux in photons
    per square centimeter per second.
    """
    defaults = (
        ('confidence',    0.95, 'Desired confidence level of the upper limit.'),
        ('simps_points',   100, 'Number of integration points (per log space).'),
        ('flux_kwargs', dict(), 'kwargs passed into i_flux function, including e_weight, cgs, emin, and emax.'),
        ('verbosity', False, 'make lots of noise')
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        """ Required arguments:
                roi - ROIAnalysis boject
                which - source for which to compute upper limit.
        """
        keyword_options.process(self, kwargs)
        self.roi = roi
        self.which = which

        self._compute()

    @staticmethod
    def get_integration_range(model):
        norm_name = model.param_names[0]
        default_norm_limits = model.default_limits[norm_name]
        norm_mapper = model.get_mapper(norm_name)
        if isinstance(norm_mapper,LimitMapper):
            lower,upper = model.get_limits(norm_name)
            norm_min = min(lower, default_norm_limits.lower)
            norm_max = max(upper, default_norm_limits.upper)
        else:
            norm_min = min(model[norm_name], default_norm_limits.lower)
            norm_max = max(model[norm_name], default_norm_limits.upper)
        return norm_min, norm_max


    def _compute(self):
        roi = self.roi
        which = self.which

        state = PointlikeState(roi)
        ll_0 = roi.logLikelihood(roi.parameters())

        source = roi.get_source(which)

        if self.verbosity:
            print ('Computing upper limit for source %s with %s spectral model' % (source.name,source.model.name))

        if not hasattr(source,'model'):
            raise Exception("upper_limit can only calculate upper limits of point and extended sources.")
        model=source.model

        integral_min, integral_max = self.get_integration_range(model)

        if self.verbosity:
            print ('For source %s, setting integration range from' % model.name )
            print (' * integration minimum = :',integral_min)
            print (' * integration maximum = :',integral_max)

        # Unbound flux temporarily to avoid parameter limits
        model.set_mapper(0,LinearMapper)

        def like(norm):
            model.setp(0,norm)
            return np.exp(ll_0-roi.logLikelihood(roi.parameters()))
        npoints = int(math.ceil(self.simps_points * (np.log10(integral_max) - np.log10(integral_min))))
        points = np.logspace(np.log10(integral_min), np.log10(integral_max),npoints*2+1)
        y = np.array([like(x)*x for x in points])
        trapz1 = integrate.cumtrapz(y[::2])
        trapz2 = integrate.cumtrapz(y)[::2]
        cumsimps = (4*trapz2 - trapz1)/3.
        cumsimps /= cumsimps[-1]
        i1 = np.where(cumsimps<.95)[0][-1]
        i2 = np.where(cumsimps>.95)[0][0]
        x1, x2 = points[::2][i1], points[::2][i2]
        y1, y2 = cumsimps[i1], cumsimps[i2]
        #Linear interpolation should be good enough at this point
        limit = x1 + ((x2-x1)/(y2-y1))*(self.confidence-y1)
        model.setp(0,limit)
        self.uflux = model.i_flux(**self.flux_kwargs)

        self.upper_limit_model = model.copy()

        state.restore(just_spectra=True)

    def get_limit(self):
        return self.uflux

def upper_limit_quick(roi,which = 0,confidence = .95,e_weight = 0,cgs = False):
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
    self=roi

    delta_logl = chi2.ppf(2*confidence-1,1)/2.
    params = self.parameters().copy()
    #self.psm.models[which]._p[0]  = -20
    zp = self.logLikelihood(self.parameters())

    def f(norm):
        self.psm.models[which]._p[0] = np.log10(norm)
        ll = self.logLikelihood(self.parameters())
        return abs(ll - zp - delta_logl)

    limit = fmin(f,np.array([10**-6]),disp=0)[0]
    self.psm.models[which]._p[0] = np.log10(limit)
    uflux = self.psm.models[which].i_flux(e_weight = e_weight,cgs = cgs)
    self.set_parameters(params)
    return uflux


class ExtensionUpperLimit(object):
    defaults = (
        ("refit_position",         False, "Refit position of source for each extension"),
        ("confidence",              0.95, "Convidence level of bayesian upper limit"),
        ("spatial_model",           None, "Spatial model to use for extnesion upper limit. Default is Disk"),
        ("delta_log_like_limits",   10,   """ delta_log_like_limits has same defintion as the parameter in
                                              pyLikelihood.IntegralUpperLimit.py function calc_int.
                                              Note, this corresponds to a change in the acutal likelihood by
                                              exp(10) ~ 20,000 which is sufficiently large that this is
                                              a pefectly fine threshold for stoping the integral."""),
        ("fit_kwargs",              dict(), "These kwargs are passed into ROIAnalysis.fit()"),
        ("spatial_model",           Disk, " Spatial model to assume during upper limit"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        """ Compute an upper limit on the source extension, by the "PDG Method". """
        keyword_options.process(self, kwargs)

        self.roi = roi
        self.which = which

        self.init_ts = roi.TS(which, quick=False)

        if self.init_ts < 4:
            # Bunt on extension upper limits for completely insignificant sources
            print ('Unable to compute extension upper limit for point-like source with too-small TS')
            self.extension_limit = None

        else:
            if not isinstance(self.spatial_model,type):
                raise Exception("The spatial model bust be a type, like Gaussian, not an instance, like Gaussian()")

            # Note, since the input is the class, not the instance, the
            # position parmaeters have not yet been added on.
            n = self.spatial_model.param_names
            assert len(n) == 1 and n[0] == 'Sigma'

            self.saved_state = PointlikeState(roi)

            self.spatial_low_lim, self.spatial_hi_lim = self.spatial_model.default_limits[0]

            results = self._compute()

            self.saved_state.restore()

    def loglike(self, extension):
        """ Perform a pointlike spectral fit for a
            given extension and return the logLikelihood.  Note,
            most robust to start with initial spectral paramaters to
            avoid situations where the previous fit totally failed to
            converge. """

        roi = self.roi
        which = self.which

        if extension < self.spatial_low_lim: extension = self.spatial_low_lim
        roi.modify(which, 
                   spatial_model=self.spatial_model(sigma=extension,
                                                    center=self.init_position,
                                                    free=np.asarray([True,True,False])),
                   keep_old_center=False)

        self.saved_state.restore(just_spectra=True)

        roi.fit(estimate_errors=False, **self.fit_kwargs)

        if self.refit_position:
            roi.fit_extension_fast(which=which, estimate_errors=False)
            roi.fit(estimate_errors=False, **self.fit_kwargs)

        ll = -roi.logLikelihood(roi.parameters()) 

        if not self.old_quiet and hasattr(self,'ll_0'):
            if self.refit_position:
                fit_position = roi.get_source(which).skydir
                position_string = ' (l,b)=(%.2f,%.2f), dist=%.2f,' % (fit_position.l(), fit_position.b(), 
                                                                 np.degrees(fit_position.difference(self.init_position)))
            else:
                position_string = ''
            print ('... sigma = %.2f,%s ll=%.2f, ll-ll_0=%.2f' % (extension, position_string, ll, ll - self.ll_0))

        return ll

    def _compute_integration_range(self):
        """ Estimate when the likelihood has fallen
            from the likelihood at Sigma=0 by an amount delta_log_like_limits. """
        roi = self.roi

        if not self.old_quiet: print ("Computing Integration range, delta_log_like_limits=%s:" % self.delta_log_like_limits)

        self.ll_0 = ll_0 = self.loglike(extension=0)

        f = lambda e: self.loglike(e) - (ll_0 - self.delta_log_like_limits)

        self.int_min = 0

        # unreasonable to have a source larger then half the ROI size.
        hi = roi.sa.maxROI/2.0 
        try:
            self.int_max = brentq(f, 0, hi, rtol=1e-4, xtol=1e-3)
        except:
            # Finding this intersect does not always work.
            print ('WARNING: Unable to find an acceptable upper limit for the integration range so defaulting to %s. Extension upper limit could be unreliable'  % hi)
            self.int_max = hi

        if not self.old_quiet: print ("Integrating range is between %s and %s" % (self.int_min, self.int_max))


    def _compute_max_loglikelihood(self):
        """ Note, it is important to evalulate the maximum loglikelihood
            so that the overall likelihood can be defined as
            exp(ll-ll_max) to avoid computing the exponential of a very
            large number in the case where ll_0 is much different from
            ll_max. Also, quad as an overall easier time since the maximum
            function value is (by defintion) equal to 1.  Since the overall
            function is normalized later, this is the most numerically
            stable normalization without causing any future trouble. """
        roi = self.roi

        if not self.old_quiet: print ("Computing maximum loglikelihood:")

        # Note, maximum loglikelihood =  minimum -1*logLikelihood
        # Note, this does not have to be very precise. The purpose of
        # this function is just to get a reasonable guess at the best extension
        # to avoid floating point issues in the likelihood l=exp(ll-ll_max).
        self.sigma_max = fminbound(lambda e: -1*self.loglike(e), self.int_min, self.int_max, disp=0, xtol=1e-3)
        self.ll_max = self.loglike(self.sigma_max)

        if not self.old_quiet: print ("Maximum logLikelihood is %.2f" % self.ll_max)

    def _compute_extension_limit(self):
        """ Compute the extnesion upper limit by

                (a) Sampling the function
                (b) Computing the CDF and normalizing
                (c) Finding when the normalized CDF is equal to the desired confidence
        
            Note, the quad accuracy parameters are roughly taken to be the same
            as in the pyLikelihood.IntegralUpperLimit.py calc_int function
            """
        roi = self.roi

        if not self.old_quiet: print ("Finding the %s quantile of the likelihood" % self.confidence)

        ll_to_l = lambda ll: np.exp(ll-self.ll_max)
        like = lambda e: ll_to_l(self.loglike(e))

        quantile = Quantile(like, self.int_min, self.int_max, 
                            quad_kwargs=dict(epsrel=1e-3, epsabs=1))
        self.extension_limit = quantile(self.confidence)

        if not self.old_quiet: print ("Extension upper limit is %.2f" % self.extension_limit)
    
    def _compute(self):

        roi = self.roi
        which = self.which

        self.old_quiet = roi.quiet
        roi.quiet = True
        self.init_position = roi.get_source(which).skydir

        self._compute_integration_range()
        self._compute_max_loglikelihood()
        self._compute_extension_limit()

    def results(self):
        return dict(extension=self.extension_limit, 
                    spatial_model = self.spatial_model.__name__,
                    confidence=self.confidence,
                    emin=self.roi.bin_edges[0],
                    emax=self.roi.bin_edges[-1],
                    delta_log_like_limits=self.delta_log_like_limits,
                    extension_units = 'degrees')
