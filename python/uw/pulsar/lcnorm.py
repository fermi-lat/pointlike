"""
A module for handling normalization of light curves with an arbitrary
number of primitive components.

This is done by treating each primitives' normalization parameter as
the square of a cartesian variable lying within or on an
n-dimensional ball of unit radius.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/lcfitters.py,v 1.27 2012/03/09 02:36:12 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>
"""

import numpy as np
from math import sin,cos,asin,acos,pi

# can some of the code be reduced with inheritance here?
# TODO -- error propagation to norms

class NormAngles(object):
    """ Keep track of N angles (0 to pi/2) representing the coordinates
        inside a unit radius N-ball."""

    def __init__(self,norms):
        self.dim = len(norms)
        if not self._check_norms(norms):
            raise ValueError('Provided norms do not satisfy constraints.')
        self.p = self._get_angles(norms)
        self.free = np.asarray([True] * len(norms))

    def _check_norms(self,norms):
        ok = True
        for n in norms:
            ok = ok and ((n > 0) and (n <= 1))
        return ok and (sum(norms)<=1)

    def _get_angles(self,norms):
        """ Determine the n-sphere angles from a set of normalizations."""
        sines = sum(norms)**0.5
        angles = [asin(sines)]
        norms = np.asarray(norms)**0.5
        for i in xrange(self.dim-1):
            phi = acos(norms[i]/sines)
            sines *= sin(phi)            
            angles.append(phi)
        return np.asarray(angles)

    def set_parameters(self,p): self.p[self.free] = p

    def get_bounds(self):
        """ Angles are always [0,pi/2). """
        return [ [0,pi/2] for i in xrange(self.free.sum()) ]

    def get_parameters(self): return self.p[self.free]

    def sanity_checks(self,eps=1e-6):
        t1 = abs(self().sum() - sin(self.p[0])**2) < eps
        return t1

    def __call__(self):
        """ Return the squared value of the Cartesian coordinates.

            E.g., for a 3-sphere, return
            z^2 = cos^2(a)*cos^2(b)
            x^2 = cos^2(a)*sin^2(b)*cos^2(c)
            y^2 = cos^2(a)*sin^2(b)*sin^2(c)

            These values are guaranteed to satisfy the constraint of
            a sum <= unity and so are suitable for normalizations of
            a light curve.
        """
        p = self.p
        m = sin(p[0])
        norms = np.empty(self.dim)
        for i in xrange(1,self.dim):
            norms[i-1] = m * cos(p[i])
            m *= sin(p[i])
        norms[self.dim-1] = m
        return norms**2

    def get_grads(self):
        """ Return a matrix giving the value of the partial derivative
            of the ith normalization with respect to the jth angle.
            M_ij = dn_i/dphi_j
        """
        m = np.zeros([self.dim,self.dim],dtype=float)
        n = self()
        p = self.p
        cots = 1./np.tan(p)
        for i in xrange(self.dim):
            for j in xrange(self.dim):
                if j > (i+1): break
                if j <= i:
                    m[i,j] = n[i]*cots[j]
                else:
                    if i==(self.dim-1):
                        m[i,j] = n[i]*cots[j]
                    else:
                        m[i,j] = -n[i]/cots[j] #-cotangent
        return 2*m

    def get_total(self):
        """ Return the amplitude of all norms."""
        return sin(self.p[0])**2

            


