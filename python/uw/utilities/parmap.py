from abc import abstractmethod
import copy


import numpy as np

class MapperException(Exception): pass

class ParameterMapper(object):

    @abstractmethod
    def toexternal(internal): pass

    @abstractmethod
    def tointernal(external): pass

    @abstractmethod
    def dexternaldinternal(external): 
        """ Compute the derivative of the external parametesr
            with respect to the internal parameters
        """
        pass

    def copy(self): return copy.deepcopy(self)

class LinearMapper(ParameterMapper):
    @staticmethod
    def toexternal(internal): return internal

    @staticmethod
    def tointernal(external): return external

    @staticmethod
    def dexternaldinternal(external):
        return 1

class LogMapper(ParameterMapper):
    """ Represent a logarithmic mapping:

            external = 10^internal or internal=log10(internal)
    """
    log_10 = np.log(10)
    @staticmethod
    def toexternal(internal): return 10**internal

    @staticmethod
    def tointernal(external): 
        """ Protect against bad input:

                >>> LogMapper.tointernal(-1)
                Traceback (most recent call last):
                    ...
                MapperException: Parameter value=-1 must be greater than 0
        """
        if external<0: 
            raise MapperException("Parameter value=%s must be greater than 0" % external)
        return np.log10(external)

    @staticmethod
    def dexternaldinternal(external):
        return external*LogMapper.log_10

class LimitMapper(ParameterMapper):
    """ A parameter mapper from the infinite (-inf,+inf) internal
        parameter space to a finite (a,b) external parameter space.

        The implementation uses the arctangent function.

        Note, scale doesn't acutally do anything, but is included for
        compatability with gtlike.
    """
    def __init__(self, lower, upper, scale=1):
        self.lower, self.upper = sorted([lower,upper])
        self.scale = scale
        self._c = (self.upper-self.lower)/np.pi
        self._pih = np.pi/2

    def toexternal(self,internal):
        """ Apply arctan transformation."""
        return self.lower + self._c*(np.arctan(internal)+self._pih)

    def tointernal(self,external):
        """ Equation 1.1 of "MINUIT User's Guide" 
        
            Make sure mapper protects against bad input:

                >>> mapper = LimitMapper(0, 10)
                >>> mapper.tointernal(-1)
                Traceback (most recent call last):
                    ...
                MapperException: Parameter value=-1 is smaller than lower limit=0
                >>> mapper.tointernal(11)
                Traceback (most recent call last):
                    ...
                MapperException: Parameter value=11 is larger than upper limit=10
        """
        if external < self.lower:
            raise MapperException("Parameter value=%s is smaller than lower limit=%s" % (external,self.lower))
        if external > self.upper:
            raise MapperException("Parameter value=%s is larger than upper limit=%s" % (external,self.upper))

        return np.tan( (external-self.lower)/self._c-self._pih )
        
    def dexternaldinternal(self,external):
        internal = self.tointernal(external)
        return self._c/(1.+internal**2)

    def __repr__(self):
        return 'LimitMapper(%s,%s,%s)' % (self.lower,self.upper,self.scale)

    def __eq__(self, other):
        """ Previously, equality testing was broken:

                >>> LimitMapper(1,5) == LimitMapper(1,5)
                True
                >>> LimitMapper(1,5) == LimitMapper(1,6)
                False
        """
        return self.lower == other.lower and self.upper == other.upper and self.scale == other.scale

class LimitMapperMinuit(ParameterMapper):
    """ Represents a mapper which limits parametesr between lower and upper

        See section 1.3.1 of "MINUIT User's Guide":

            http://seal.web.cern.ch/seal/documents/minuit/mnusersguide.pdf

        Note, scale doesn't acutally do anything, but is included for
        compatability with gtlike.
    """
    def __init__(self, lower, upper, scale=1):
        self.lower, self.upper = sorted([lower,upper])
        self.scale = scale

    def toexternal(self,internal):
        """ Equation 1.2 of "MINUIT User's Guide" """
        return self.lower + ((self.upper-self.lower)/2.0)*(np.sin(internal)+1.0)

    def tointernal(self,external):
        """ Equation 1.1 of "MINUIT User's Guide" 
        
            Make sure mapper protects against bad input:

                >>> mapper = LimitMapper(0, 10)
                >>> mapper.tointernal(-1)
                Traceback (most recent call last):
                    ...
                MapperException: Parameter value=-1 is smaller than lower limit=0
                >>> mapper.tointernal(11)
                Traceback (most recent call last):
                    ...
                MapperException: Parameter value=11 is larger than upper limit=10
        """
        if external < self.lower:
            raise MapperException("Parameter value=%s is smaller than lower limit=%s" % (external,self.lower))
        if external > self.upper:
            raise MapperException("Parameter value=%s is larger than upper limit=%s" % (external,self.upper))

        return np.arcsin(2.0*(external-self.lower)/(self.upper-self.lower)-1.0)
        
    def dexternaldinternal(self,external):
        return ((self.upper-self.lower)/2.0)*np.cos(self.tointernal(external))

    def __repr__(self):
        return 'LimitMapper(%s,%s,%s)' % (self.lower,self.upper,self.scale)

    def __eq__(self, other):
        """ Previously, equality testing was broken:

                >>> LimitMapper(1,5) == LimitMapper(1,5)
                True
                >>> LimitMapper(1,5) == LimitMapper(1,6)
                False
        """
        return self.lower == other.lower and self.upper == other.upper and self.scale == other.scale

if __name__ == "__main__":
    import doctest
    doctest.testmod()
