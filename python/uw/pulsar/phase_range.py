"""
This module defines an object PhaseRange which represents a
(possibly non-contiguous) range in pulsar phase. 

See the docstring for usage information.

This object has SymPy as a dependency.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/pulsar/phase_range.py,v 1.3 2011/10/04 22:03:23 lande Exp $

author: J. Lande <joshualande@gmail.com>

"""
import sympy
import numbers
import numpy as np
from operator import isNumberType
import copy

class PhaseRange(object):
    """ This Object represents a (possibly non-contiguous) range  in pulsar 
        phase. This should automate a lot of otherwise
        tedious pulsar phase calculations.

        from phase_range import PhaseRange
            

        By default, an empty set is created, that
        is not very useful

            >>> print PhaseRange()
            EmptySet()

        The simplest thing to do is to define a single phase

            >>> print PhaseRange(0.2, 0.5)
            [0.2, 0.5]

        But you can also define several phase ranges
        
            >>> print PhaseRange([0.2, 0.5], [0.6, 0.8])
            [0.2, 0.5] U [0.6, 0.8]


        The object will automatically wrap phases numbers between -1 and 2

            >>> print PhaseRange(0.5, 1.2)
            [0, 0.2] U [0.5, 1]

            >>> print PhaseRange(0.5, -0.8)
            [0, 0.2] U [0.5, 1]

        But for simplicity, the object will crash on numbers >=2 or <=-1

            >>> PhaseRange(0, 2)
            Traceback (most recent call last):
                ...
            Exception: Error, phase range [0, 2] is outside allowed range.

            >>> PhaseRange(-1,0)
            Traceback (most recent call last):
                ...
            Exception: Error, phase range [-1, 0] is outside allowed range.


        And this code can conveniently merge multiple phases

            >>> print PhaseRange([0.2, 0.5], [0.4, 0.7])
            [0.2, 0.7]

        Or, you can pass in a sigle list if this is easir:

            >>> print PhaseRange([[0.2, 0.5], [0.4, 0.7]])
            [0.2, 0.7]

        If the first number is > the second number,
        it will interpret this as the phase range wrapinga round.

            >>> print PhaseRange(0.5, 0.2)
            [0, 0.2] U [0.5, 1]

            >>> print PhaseRange(1.5, 0.2)
            [0, 0.2] U [0.5, 1]
        
        Finally, the object fairly intellegently deals with the
        edge case of no range:

            >>> print PhaseRange(0.5, 0.5)
            {0.5}

            >>> print PhaseRange(-0.5, 0.5)
            [0, 1]

        The phase_fraction can be easily calulated:

            >>> print PhaseRange(1.5, 0.25).phase_fraction
            0.75

            >>> print PhaseRange(0, 1).phase_fraction
            1.0

        The 'in' operator has been simply defined:

            >>> print .3 in PhaseRange(0.1, 0.5)
            True
            >>> print .6 in PhaseRange(0.1, 0.5)
            False

            >>> print PhaseRange(0.2, 0.4) in PhaseRange(0.1, 0.5)
            True
            >>> print PhaseRange(0.4, 0.6) in PhaseRange(0.1, 0.5)
            False

        This function can be very usful for complicated phases

            >>> print PhaseRange([0.2, 0.3], [0.6, 0.7]) in \
                      PhaseRange([0.1, 0.4], [0.5, 0.8])
            True
            >>> print PhaseRange([0.2, 0.3], [0.6, 0.7]) in PhaseRange(0.1, 0.6)
            False


        Similarly, you can add two phases:

            >>> print PhaseRange(0.2, 0.5) + PhaseRange(0.6, 0.8)
            [0.2, 0.5] U [0.6, 0.8]


        The tolist() function can be used exported the ranges to a simple list:

            >>> print PhaseRange(0.25, 0.75).tolist()
            [0.25, 0.75]


            >>> print PhaseRange([0.25, 0.5], [0.75, 1]).tolist()
            [[0.25, 0.5], [0.75, 1.0]]

            >>> print PhaseRange(0.5, 0.25).tolist()
            [0.5, 0.25]

        This should also work in a couple of edge cases

            >>> print PhaseRange(0.5, 0.5).tolist()
            [0.5, 0.5]

            >>> print PhaseRange().tolist()
            []

        To get out a more verbose list of phases, which does not
        have assumed automatic wraping and is always a list of lists:

            >>> print PhaseRange(0.5, 0.25).tolist(dense=False)
            [[0.0, 0.25], [0.5, 1.0]]
            >>> print PhaseRange(0.25, 0.75).tolist(dense=False)
            [[0.25, 0.75]]

        Note, tolist() can be used to convert from and back to this object:

            >>> p=PhaseRange(0.5, 0.2)
            >>> PhaseRange(p.tolist()) == p
            True

        You can make a copy easily of the object:

            >>> p = PhaseRange(0.2, 0.5)
            >>> print PhaseRange(p)
            [0.2, 0.5]

        There is also the intersect and overlaps function: 

            >>> print PhaseRange(.25,.75).intersect(PhaseRange(0.5,1)) 
            [0.5, 0.75]

            >>> print PhaseRange(.25,.5).intersect(PhaseRange(0.5,.75)) 
            {0.5}

            >>> PhaseRange(.25,.75).overlaps(PhaseRange(0.5,1)) 
            True
            >>> PhaseRange(.25,.5).overlaps(PhaseRange(0.5,.75)) 
            False

    """

    # All input phases must be between 0 and 2
    allowed_phase_input = sympy.Interval(-1,2,left_open=True, right_open=True)

    def __init__(self,*args):
        if len(args) == 1 and isinstance(args[0],PhaseRange):

            self.range = copy.deepcopy(args[0].range)
            return

        if len(args) == 1 and np.alltrue(len(i)==2 for i in args):
            args = args[0]

        if len(args) == 2 and \
           isinstance(args[0],numbers.Real) and \
           isinstance(args[1],numbers.Real):
            args = [ list(args) ]

        ranges = []

        for range in args:

            if range[0] not in PhaseRange.allowed_phase_input or \
               range[1] not in PhaseRange.allowed_phase_input:
                raise Exception("Error, phase range %s is outside allowed range." % str(range))

            if np.allclose(range[0]-range[1],1) or \
               np.allclose(range[1]-range[0],1):
                range = [0,1]
            else:
                for i in [0,1]:
                    if range[i] != 1: range[i] %= 1

            if range[0] > range[1]: 
                # pulsar convention
                ranges += [
                    sympy.Interval(0, range[1]),
                    sympy.Interval(range[0], 1)
                ]
            else:
                ranges.append(sympy.Interval(*range))

        self.range = sympy.Union(*ranges)

    @property
    def phase_fraction(self):
        """ Fraction of pulsar phase that the range occupies. """
        return float(self.range.measure)

    def __str__(self):
        return self.range.__str__()

    def __in__(self): pass

    @staticmethod
    def _tolist(obj):
        """ Recursivly convert a simpy object to a simple 
            list representing the phase range. """
        if isinstance(obj,sympy.Union):
            # special case if range is [[0, x], [y, 1]] -> return [y, x]
            if len(obj.args) == 2 and \
               isinstance(obj.args[0],sympy.Interval) and \
               isinstance(obj.args[1],sympy.Interval) and \
               obj.args[0].start == 0 and \
               obj.args[1].end == 1:
                return map(float,[obj.args[1].start,obj.args[0].end])

            return map(PhaseRange._tolist, obj.args)
        elif isinstance(obj,sympy.Interval):
            return map(float,[obj.start,obj.end] )
        elif isinstance(obj, sympy.FiniteSet):
            if len(obj.args) == 1:
                return [float(obj.args[0]), float(obj.args[0])]
            else:
                return [ [float(i), float(i)] for i in obj.args ]
        elif isinstance(obj, sympy.EmptySet):
            return []
        else:
            raise Exception("Unrecognized type %s" % type(obj))

    def tolist(self,dense=True):
        r = self.range
        if dense:
            return PhaseRange._tolist(r)
        else:
            if isinstance(r,sympy.Interval):
                return [PhaseRange._tolist(r)]
            elif isinstance(r, sympy.FiniteSet):
                return []
            else:
                return [PhaseRange._tolist(i) for i in r.args]

    def __eq__(self, other):
        if isinstance(other,PhaseRange):
            return self.range == other.range
        else:
            return self.range == other

    def __contains__(self, other):
        if isNumberType(other):
            return self.range.contains(float(other))
        elif isinstance(other, sympy.Set):
            return self.range.intersect(other) == other
        elif isinstance(other, PhaseRange):
            return self.range.intersect(other.range) == other.range
        else:
            raise Exception("Unknown type %s" % type(other))

    def __add__(self, other):
        new=PhaseRange()
        new.range += self.range
        new.range += other.range
        return new

    def intersect(self, other):
        if isinstance(other, PhaseRange):
            new = PhaseRange()
            new.range = self.range.intersect(other.range)
            return new
        else:
            raise Exception("Unknown type %s" % type(other))

    def overlaps(self, other):
        """ True if the two regions have overlaping range. """
        return self.intersect(other).phase_fraction > 0

if __name__ == "__main__":
    import doctest
    doctest.testmod()
