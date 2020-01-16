"""
This module implements the class OffPeak that
can be used to calculate the off-peak region
of a pulsar.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/pulsar/lc_off_peak.py,v 1.5 2011/10/15 23:02:06 lande Exp $

author: J. Lande <joshualande@gmail.com>

"""
import numbers

import sympy
import numpy as np
from scipy.optimize import fmin, brute, brentq
from scipy.stats import poisson

        
from uw.utilities import keyword_options
from . phase_range import PhaseRange
from . lcfitters import WeightedLCFitter


class OffPeak(object):
    """ Object that can be used to define the off peak 
        window for a pulsar.

        This object requires an lcfitters object.

        Usage:

            off_peak = OffPeak(lcf)
            region = off_peak.off_peak # phase_range.PhaseRange object
            region_list = region.tolist()
    
        N.B. in this implementation, the pulsar is defined
        as anything above the lowest part of the light curve
        model. Even if the Pulsar model predicts emission
        at all phases, we zero the pulsar at the minimum
        value. This is, of course, a somewhat arbitrary choice,
        but so is the 1% contamination criteria in the first place.
    """

    defaults = (
        ("contamination",  0.01, "fraction of pulsar allowed in off peak"),
        ("quiet",         False, "Don't print out"),
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

    def _find_range(self, center, lowest_value):
        c = center
        lct = self.lct

        def integral(dphi):
            """ Calculate the fraction of the pulsar
                light curve within the phase range
                dphi of the center.

                N.B: this formula is true because by definition 
                lct.integrate(0, 1) = 1. 
                
                N.B: the lct.integrate function does not work
                well with wraping around when c+dphi/2 > 1,
                so here we split up the integral into parts using
                the PhaseRange object. """
            l=lowest_value

            r=PhaseRange(c-dphi/2,c+dphi/2)
            integral = sum(lct.integrate(a,b) for a,b in r.tolist(dense=False))

            return (integral - l*(dphi))/(1-l)

        residual=lambda dphi: integral(dphi) - self.contamination

        dphi=brentq(residual, 0, 1, disp=False, rtol=self.contamination/10)

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
        range = PhaseRange(a,b)
        alt_range = PhaseRange(b,a)
        if alt_range.phase_fraction < range.phase_fraction:
            range = alt_range
        return sum(peak in range for peak in peaks)

    @staticmethod
    def _lowest_value(lct, exclude_phase=None):

        phases=np.linspace(0,1,1001)
        if exclude_phase is not None:
            phases = [p for p in phases if p not in exclude_phase]

        return min([lct(p) for p in phases])

    @staticmethod
    def consistent(x,xphase,y,yphase,probability = 0.05, quiet=True):
        """ Assuming x counts are observed in a phase range xphase
            and y counts are observed in phase range yphase,
            decides if the regions are consistent.
            
            The regions are consistent if the probability of obtaining 
            as many or more counts in the second region compared to the first region
            is > 5% (so there is a >5% probability that ht esecond region is
            not unusually large). """
            
        y_predicted = x*(yphase/xphase)
            
        poisson_likelihood = poisson.sf(y,y_predicted)

        if not quiet:
            print ('poisson likelihood=%.2f' % poisson_likelihood)

        if poisson_likelihood > 0.05: return True
        return False

    def _find_all_ranges(self, exclude_phase=None):

        lct = self.lct

        # (*) Calulate the lowest part of the light curve using a simple
        #     grid search.

        lowest_value = OffPeak._lowest_value(lct, exclude_phase=exclude_phase)

        # (*) Calculate the phase range for a grid of centers in phase.


        centers = np.linspace(0,1,101) 
        if exclude_phase is not None:
            centers=[c for c in centers if c not in exclude_phase]
        ranges = [self._find_range(c, lowest_value) for c in centers]
        dphis = np.asarray([range.phase_fraction for range in ranges])
        
        return centers, ranges, dphis

    def _fit(self):

        lcf = self.lcf
        lct = self.lct

        centers, ranges, dphis = self._find_all_ranges()
        
        i = np.argmax(dphis)
        self.first_off_peak = ranges[i]
        self.first_center = centers[i]
        self.first_dphi = dphis[i]

        if not self.quiet:
            print ('Best Off Peak region is ',self.first_off_peak)

        # (*) Find the largest second region

        # (*) Get a list of peaks
        peaks = OffPeak.peak_locations(lct)

        if len(peaks) < 2:
            print ('No good off pulse regions')
            self.off_peak = self.first_off_peak
            return

        # define a new set of ranges subtracting a different minimum pulsar phase
        centers, ranges, dphis = self._find_all_ranges(exclude_phase = self.first_off_peak)

        peaks_between = np.asarray([OffPeak.peaks_between(self.first_center,c,peaks) for c in centers])
        overlaps = np.asarray([self.first_off_peak.overlaps(r) for r in ranges])

        good_region = (~overlaps) & \
                (peaks_between > 0) & \
                (len(peaks)>1) & \
                (dphis > self.first_dphi/2.0)

        if np.sum(good_region) < 1:
            print ('No good off pulse regions')
            self.off_peak = self.first_off_peak
            return

        i = np.argmax(np.where(~good_region,-np.inf,dphis))
        self.second_off_peak = ranges[i]

        if not self.quiet:
            print ('Second off pulse region candidate is', self.second_off_peak)

        if isinstance(lcf,WeightedLCFitter):
            print ('Note, this may not be the optimal way to calculate number of counts (at least, when pulsar weighting...!')

        # (*) Test if second regions are statistically consistent with eachother

        first_counts  = self.get_counts(self.first_off_peak,  lcf.phases)
        second_counts = self.get_counts(self.second_off_peak, lcf.phases)

        consistent = OffPeak.consistent(first_counts,self.first_off_peak.phase_fraction,
                                        second_counts,self.second_off_peak.phase_fraction,
                                        quiet=self.quiet)

        if consistent:
            if not self.quiet: print ('Keeping second region')
            print ('phase_fraction = %.2f, %.2f' % (self.first_off_peak.phase_fraction, self.second_off_peak.phase_fraction))
            self.off_peak = self.first_off_peak + self.second_off_peak
        else:
            if not self.quiet: print ('Rejecting second region')
            self.off_peak = self.first_off_peak


