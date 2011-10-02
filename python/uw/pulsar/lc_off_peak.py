"""
This module implements the class OffPeak that
can be used to calculate the off-peak region
of a pulsar.

$Header: $

author: J. Lande <joshualande@gmail.com>

"""
import sympy
import numbers
        
from . phase_range import PhaseRange


class OffPeak(object):
    """ Object that can be used to define the off pulse 
        window for a pulsar.

        This object requires an lcfitters object.
    
    
    Note, in this implementation, the pulsar is defined
        as anything above the lowest part of the light curve
        model. Even if the Pulsar model predicts emission
        at all phases, we zero the pulsar at the minimum
        value. This is, of course, a somewhat arbitrary choice,
        but so is the 1% contaminuation criteria in the first place.
    """

    contamination = 0.01 # fraction of pulsar allowed in off peak

    def __init__(self,lct):
        self.lct = lct # uw.pulsar.lcfitters.LCTemplate

        self._fit()

    def _find_range(center):
        c = center
        lct = self.lct



        integral=lambda dphi: lct.integrate(c-dphi/2, c+dphi/2) - self.lowest_value*(dphi)
        total_integral = lct.integrate(0, 1) - self.lowest_value
 
        dphi_guess=0
        while integral(dphi_guess)/total_integral < OffPeak.contamination:
            dphi_guess += .01

        residual=lambda dphi: np.abs(integral(dphi)/total_integral - OffPeak.contamination)

        # N.B., the integral over all phase is equal to 1
        dphi=fmin(residual, dphi_guess, disp=False)[0]

        int = integral(dphi)/total_integral
        print ' center = %s, dphi = %s, int=%s' % (center,dphi,int)

        return PhaseRange(center-dphi/2,center+dphi/2)

    def _fit(self):
        lct = self.lct

        # (*) Calulate the lowest part of the light curve using a simple
        #     grid search.

        lowest_phase, lowest_value, grid, jout = scipy.optimize.brute(lct, (0,1), ns=1000)
        self.lowest_value = lowest_value

        # (*) Calculate the phase range for a grid of centers in phase.

        centers=np.linspace(0,1,101)
        ranges = np.asarray([self._find_range(c) for c in np.linspace(0,1,101)])
        
        i = np.argmax(dphis)
        self.off_pulse = ranges[i]

        # (*) Find the second largest region that does not overlap the first
        #     region.

        overlap = np.asarray([])
        j = np.argmax(dphis[ overlap == False ])
        second_off_pulse = [ min_phase[ overlap == False ][j], max_phase[ overlap == False ][j] ]

        # (*) Test if regions are statistically consistent with eachother

        self.ncounts1 = self.get_counts()
        self.ncounts2 = self.get_counts()

        if regions_are_consistent:
            self.off_pulse += second_off_pulse

    @property
    def off_pulse(self):
        return self.off_pulse



