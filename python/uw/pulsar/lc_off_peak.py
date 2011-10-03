"""
This module implements the class OffPeak that
can be used to calculate the off-peak region
of a pulsar.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/pulsar/lc_off_peak.py,v 1.1 2011/10/02 23:37:22 lande Exp $

author: J. Lande <joshualande@gmail.com>

"""
import numbers

import sympy
import numpy as np
from scipy.optimize import fmin, brute
from scipy.stats import poisson

        
from uw.utilities import keyword_options
from . phase_range import PhaseRange


class OffPeak(object):
    """ Object that can be used to define the off peak 
        window for a pulsar.

        This object requires an lcfitters object.
    
    
    Note, in this implementation, the pulsar is defined
        as anything above the lowest part of the light curve
        model. Even if the Pulsar model predicts emission
        at all phases, we zero the pulsar at the minimum
        value. This is, of course, a somewhat arbitrary choice,
        but so is the 1% contaminuation criteria in the first place.
    """

    defaults = (
        ("contamination", 0.01,"fraction of pulsar allowed in off peak"),
        ("quiet",False,"Don't print out"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, lcf, **kwargs):
        keyword_options.process(self, kwargs)

        self.lcf = lcf # uw.pulsar.lcfitters.LCFitter

        self.lct = lcf.template # uw.pulsar.lcfitters.LCTemplate

        self._fit()

    @staticmethod
    def get_counts(phase_range,phases):
        """ Return the number of photons in the given phase range. """
        return sum(phase in phase_range for phase in phases)

    def _find_range(self, center):
        c = center
        lct = self.lct

        # Note, by definition lct.integrate(0, 1) = 1

        def integral(dphi):
            l=self.lowest_value
            return (lct.integrate(c-dphi/2, c+dphi/2) - l*(dphi))/(1-l)
 
        # get a better starting guess
        dphi_guess=0
        while integral(dphi_guess) < self.contamination:
            dphi_guess += .01

        residual=lambda dphi: np.abs(integral(dphi) - self.contamination)

        # N.B., the integral over all phase is equal to 1
        dphi=fmin(residual, dphi_guess, disp=False)[0]

        if not self.quiet: 
            print ' center = %s, dphi = %s (int=%s)' % (center,dphi,integral(dphi))

        return PhaseRange(center-dphi/2,center+dphi/2)

    @staticmethod
    def peak_location(primitive):
        return primitive.p[np.asarray(primitive.pnames) == 'Location'][0]

    @staticmethod
    def peak_locations(lct):
        return [OffPeak.peak_location(p) for p in lct.primitives]

    @staticmethod
    def peaks_between(a,b,peaks):
        """ Calculate number of peaks between a and b
            where peaks is a list of peaks. """

        # make sure a,b are from 0-1 and b>a
        a,b=a%1,b%1
        a,b=min(a,b),max(a,b)
        if b-a < (a-0) + (1-b):
            range = PhaseRange(a,b)
        else:
            range = PhaseRange(b,a)
        return sum(peak in range for peak in peaks)

    def _fit(self):

        lct = self.lct
        lcf = self.lcf

        # (*) Calulate the lowest part of the light curve using a simple
        #     grid search.

        grid = (0,1, .001)
        lowest_phase = brute(lct, (grid,))
        self.lowest_value = lct(lowest_phase)

        # (*) Calculate the phase range for a grid of centers in phase.

        centers=np.linspace(0,1,101)
        ranges = [self._find_range(c) for c in centers]
        dphis = np.asarray([range.phase_fraction for range in ranges])
        
        i = np.argmax(dphis)
        self.first_off_peak = ranges[i]
        self.first_center = centers[i]
        self.first_dphi = centers[i]

        if not self.quiet:
            print 'Best Off Peak region is ',self.first_off_peak

        # (*) Find the largest good alternate region

        peaks = OffPeak.peak_locations(lct)
        peaks_between = np.asarray([OffPeak.peaks_between(self.first_center,c,peaks) for c in centers])
        overlaps = np.asarray([self.first_off_peak.overlaps(r) for r in ranges])

        good_region = (~overlaps) & \
                (peaks_between > 0) & \
                (len(peaks)>1)
#                (dphis > self.first_dphi/8.0) & \

        if len(peaks) < 2:
            print 'No good off pulse regions'
            self.off_peak = self.first_off_peak
            return

        i = np.argmax(np.where(~good_region,-np.inf,dphis))
        self.second_off_peak = ranges[i]

        if not self.quiet:
            print 'Second off pulse region candidate is', self.second_off_peak

        # (*) Test if second regions are statistically consistent with eachother

        counts = self.get_counts(self.first_off_peak,lcf.phases)

        # (*) Calculate total counts in the second region

        print 'Note, this may not be the optimal way to calculate number of counts (at least, when pulsar weighting...!'
        second_counts = self.get_counts(self.second_off_peak, lcf.phases)

        # counts predicted in second region if no pulsed emission!
        prediced_counts = counts*(self.second_off_peak.phase_fraction/self.first_off_peak.phase_fraction)

        poisson_likelihood = poisson.sf(second_counts,prediced_counts)

        if not self.quiet:
            print 'Poisson likelihood for second region is %.3f' % poisson_likelihood

        # keep second region if probability of obtaining more counts than observed
        # by the second region is > 5% (the second region is not unusually
        # large).
        if poisson_likelihood > 0.05:
            if not self.quiet: print 'Keeping second region'
            self.off_peak = self.first_off_peak + self.second_off_peak
        else:
            if not self.quiet: print 'Rejecting second region'
            self.off_peak = self.first_off_peak


