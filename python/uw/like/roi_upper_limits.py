"""
Module to calculate flux and extension upper limits.

$Header:$

author:  Eric Wallace <ewallace@uw.edu>, Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
from scipy import integrate
from scipy.stats.distributions import chi2
from scipy.optimize import fmin

from uw.utilities import keyword_options

from uw.like.roi_state import PointlikeState
from uw.like.roi_extended import ExtendedSource
from uw.like.SpatialModels import Disk

def upper_limit(roi, which=0,
              confidence=0.95,
              integral_min=-15,
              integral_max =-8,
              simps_points = 100,
              **kwargs):
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

    All other arguments are passed into the uw.like.Models.i_flux function,
    including e_weight, cgs, emin, and emax.
    """
    self=roi
    params = self.parameters().copy()
    ll_0 = self.logLikelihood(self.parameters())

    source = self.get_source(which)
    if not source.__dict__.has_key('model'):
        raise Exception("upper_limit can only calculate upper limits of point and extended sources.")
    model=source.model

    def like(norm):
        model.setp(0,norm,internal=True)
        return np.exp(ll_0-self.logLikelihood(self.parameters()))
    npoints = simps_points * (integral_max - integral_min)
    points = np.log10(np.logspace(integral_min, integral_max,npoints*2+1))
    y = np.array([like(x)*10**x for x in points])
    trapz1 = integrate.cumtrapz(y[::2])
    trapz2 = integrate.cumtrapz(y)[::2]
    cumsimps = (4*trapz2 - trapz1)/3.
    cumsimps /= cumsimps[-1]
    i1 = np.where(cumsimps<.95)[0][-1]
    i2 = np.where(cumsimps>.95)[0][0]
    x1, x2 = points[::2][i1], points[::2][i2]
    y1, y2 = cumsimps[i1], cumsimps[i2]
    #Linear interpolation should be good enough at this point
    limit = x1 + ((x2-x1)/(y2-y1))*(confidence-y1)
    model.setp(0,limit,internal=True)
    uflux = model.i_flux(**kwargs)
    self.logLikelihood(params)
    return uflux

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
        ("confidence", 0.95, "Convidence level of bayesian upper limit"),
        ("spatial_model", None, "Spatial model to use for extnesion upper limit. Default is Disk"),
        ("integral_min", 1e-10, "Minimum extension in integration (in degrees)"),
        ("integral_max", 3, "Maximum extension in integration (in degrees)"),
        ("npoints", 100, "Number of integration points"),
        ("delta_log_like_limits", 10, """delta_log_like_limits has same defintion as the parameter in
                                      pyLikelihood.IntegralUpperLimit.py function calc_int.
                                      Note, this corresponds to a change in the acutal likelihood by
                                      exp(10) ~ 20,000 which is sufficiently large that this is
                                      a pefectly fine threshold for stoping the integral."""),
        ("fit_kwargs", dict(), 'These kwargs are passed into ROIAnalysis.fit()'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        """ Compute an upper limit on the source extension, by the "PDG Method".
            The function roi_upper_limits.upper_limit is similar, but computes
            a extension upper limit. """
        keyword_options.process(self, kwargs)

        self.roi = roi
        self.which = which


        if self.spatial_model is None: 
            self.spatial_model = Disk()

        saved_state = PointlikeState(roi)

        results = self._compute()

        saved_state.restore()


    def loglike(self, extension):
        roi = self.roi
        which = self.which
        roi.modify(which, sigma=extension)
        roi.fit(estimate_errors=False, **self.fit_kwargs)
        ll=-roi.logLikelihood(roi.parameters())
        if self.ll_0 is None: self.ll_0 = ll

        dll = ll-self.ll_0
        if not roi.old_quiet: print 'sigma = %.2f, ll=%.2f, dll=%.2f' % (extension, ll, dll)
        return dll

    def _compute(self):
        roi = self.roi
        which = self.which

        if roi.TS(which)<4:
            # Bunt on extension upper limits for completely insignificant sources
            self.extension_limit = None

        else:
            self.ll_0 = None

            roi.old_quiet = roi.quiet
            roi.quiet = True

            source = roi.get_source(which)
            if not isinstance(source,ExtendedSource):
                roi.modify(which, spatial_model=self.spatial_model, keep_old_center=True)


            extensions = np.linspace(self.integral_min, self.integral_max, self.npoints+1)
            extension_middles = 0.5*(extensions[1:] + extensions[:-1])

            ll = -np.inf*np.ones_like(extensions)
            for i in range(len(extensions)):
                ll[i] = self.loglike(extensions[i])
                if i > 0 and ll[i] < max(ll[0:i]) - self.delta_log_like_limits: break

            likelihood = np.exp(ll)

            cdf = integrate.cumtrapz(x=extensions,y=likelihood)
            cdf /= cdf[-1]
            self.extension_limit = np.interp(self.confidence, cdf, extension_middles)

            roi.quiet = roi.old_quiet

    def results(self):
        return dict(extension=self.extension_limit, 
                    spatial_model = self.spatial_model.name,
                    confidence=self.confidence,
                    emin=self.roi.bin_edges[0],
                    emax=self.roi.bin_edges[-1],
                    delta_log_like_limits=self.delta_log_like_limits,
                    extension_units = 'degree')
