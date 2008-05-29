
"""

 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit.py,v 1.4 2008/03/31 09:05:05 burnett Exp $
"""
# setup to import pointlike
try: import uw.pointlike
except: pass

import os, sys, types

from pointlike import DiffuseFunction, CompositeSkySpectrum

#--------------------------------------------------------


class Background(object):
    """ manage the background
    """
    def __init__(self, diffusefilename, exposure=3e10):
        self.background=None
        if diffusefilename is None: return
        self.exposure = exposure
        
        print 'setting up background from file %s' % diffusefilename
        if 'GLAST_EXT' in os.environ and diffusefilename=='galdiffuse':
            diffusefilename = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'GP_gamma.fits')

        self._df = DiffuseFunction(diffusefilename) # need to keep this reference open!

        if type(exposure)==types.FloatType:
            print 'applying constant exposure: %.2g' % exposure
            self.background = CompositeSkySpectrum(self._df, exposure)
        else:
            print 'exposure file" %s' %exposure
            raise Exception('Use of exposure file not yet implemented')
            
    def __call__(self):
        return self.background
#--------------------------------------------------------

if __name__=='__main__':
    pass
    #test basic setup
    b = Background('galdiffuse') 
