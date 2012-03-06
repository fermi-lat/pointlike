from collections import OrderedDict

def memoize(function):
    """ This decorator allows painless caching of a function.

        The function is inspired by
        
        http://programmingzen.com/2009/05/18/memoization-in-ruby-and-python/,

        But uses an object oriented return to allow eaiser access to the cached data. 

        For example:

            >>> @memoize
            ... def f(x): return x**2
            >>> f(1), f(2), f(4)
            (1, 4, 16)
            >>> print f.x()
            [1, 2, 4]
            >>> print f.y()
            [1, 4, 16]
    """
    class Mem(object):
        def __init__(self):
            self.cache = OrderedDict()
            self.function = function
        def __call__(self, *args):
            try:
                return self.cache[args]
            except KeyError:
                val = self.function(*args)
                self.cache[args] = val
                return val
        def x(self):
            """ More helpful to convert tuples of a single item into the item. """
            return [k if len(k) != 1 else k[0] for k in self.cache.keys()]
        def y(self):
            return self.cache.values()
    return Mem()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
