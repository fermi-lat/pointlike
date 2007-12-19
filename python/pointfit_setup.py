#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.5 2007/11/19 01:09:36 burnett Exp $
from  pointlike_defaults import *

#  specify files with FT1 data and points to fit. pixelfile for PhotonMap, files for FT1 or merit
def test():
  " define the pixelfile for a quick test, running the pixel file created by the test program"
  import os
  path = os.environ['POINTLIKEROOT']
  return os.path.join(path, 'src', 'test', 'pointlike_test.fits')

pixelfile = test()
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'


# the troublesome triplet: the 3EG blazar is very strong, affects the HLCloud
name = ['DC2_3EGJ1635m1751', 'HLCloud_SC1_05', 'Giommi_blazar_1237']
ra   = [248.788,    248.4804, 248.34]
dec  = [-17.861,   -18.294,  -18.71]
 
