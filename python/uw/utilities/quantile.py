""" Module to compute quanties.
    
    Author: Joshua Lande <joshualande@gmail.com>
"""
from scipy.integrate import quad, cumtrapz
from scipy.interpolate import interp1d
import numpy as np

from . decorators import memoize

class Quantile(object):
    """ Computes the quantile for a given distribution.
        This is useful for calculating upper limits
        and containment radii (like r68).

        Example: the quantile for a gaussian distribution should be the same as the Percent point function (ppf)

            >>> import numpy as np
            >>> from scipy.stats import norm
            >>> n = norm.pdf
            >>> q = Quantile(n, -np.inf, np.inf, quad_kwargs=dict(epsabs=1e-15, epsrel=1e-15))
            >>> np.allclose(q.x68(), norm.ppf(.68), rtol=1e-4)
            True
            >>> np.allclose(q(0.5), norm.ppf(0.5), rtol=1e-4)
            True
            >>> np.allclose(q(.95), norm.ppf(0.95), rtol=1e-3)
            True
    """
    def __init__(self, integrand, xmin, xmax, quad_kwargs=dict()):

        self.integrand = integrand
        self.cache = memoize(integrand)
        self.xmin = xmin
        self.xmax = xmax
        self.quad_kwargs = quad_kwargs

        
        self._compute_cdf()

    def _compute_cdf(self):
        """ Evaluate the CDF using the quad function. """

        k = dict(full_output=True)
        k.update(self.quad_kwargs)

        quad(self.cache, self.xmin, self.xmax, **k)

        x = np.asarray(self.cache.x(), dtype=float)
        y = np.asarray(self.cache.y(), dtype=float)

        index = np.argsort(x)
        x, y = x[index], y[index]

        self.cdf = cumtrapz(y,x)
        self.cdf = np.append([0], self.cdf)
        self.cdf /= self.cdf[-1]

        self.x = x

        self.interp = interp1d(x=self.cdf, y=self.x, kind='linear')

    def __call__(self, quantile):
        return self.interp(quantile)

    def x68(self): 
        return self(0.68)
    def x99(self): 
        return self(0.99)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
