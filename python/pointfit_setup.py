#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.7 2008/01/08 22:40:15 burnett Exp $
from  pointlike_defaults import *

#  specify files with FT1 data and points to fit. pixelfile for PhotonMap, files for FT1 or merit
def test():
  " define the pixelfile for a quick test, running the pixel file created by the test program"
  import os
  path = os.environ['POINTLIKEROOT']
  return os.path.join(path, 'src', 'test', 'pointlike_test.fits')

# data selection: either "pixelfile", or "files", the latter a list of FT1 files
pixelfile = test()
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'

# if this is non-zero, lots of output
PointSourceLikelihood.verbose=0

# if this is non-zero, use the first of the list as a background for the remainder
first_is_center=0


# the troublesome triplet: the 3EG blazar is very strong, affects the HLCloud
name = ['DC2_3EGJ1635m1751', 'HLCloud_SC1_05', 'Giommi_blazar_1237', 'bogus1', 'bogus2']
ra   = [248.788,    248.4804, 248.34, 248.51, 248.27]
dec  = [-17.861,   -18.294,  -18.71  ,-17.88, -18.12]
 
PointSourceLikelihood.verbose=1
name= ['vela']
ra = [128.8359]; dec=[-45.1763]
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
#files = [r'F:\glast\data\SC2\interleave\Interleave_pruned.root']
#files = glob.glob(r'F:\glast\data\octobertest\merit\cutfile*.root')
#files = [r'F:\glast\data\octobertest\merit\cutfile8.root']
#files = glob.glob(r'F:\glast\data\octobertest\FT\*ft1.fits')

#Data.files=[r'F:\glast\data\55-day\Skimmer_pruned.root']
#Data.output_pixelfile = '55-day_pmap.fits'
#Data.pixelfile = '55-day_pmap.fits'
#PointSourceLikelihood.skip1=3 # try skipping first

# setup for selecting a range.
#tzero = 252460800
#week = 7*86400
#Data.start_time = tzero+7*week
#Data.stop_time = Data.start_time+week
#Data.output_pixelfile = 'week8_pmap.fits'
Data.pixelfile='week1_pmap.fits'

#dec=[-45.16]
ra=[128.835939]; dec=[-45.177672] # best level-10 fit?
dir =(128.836673, -45.188701)
ra=[dir[0]]; dec=[dir[1]]
PointSourceLikelihood.skip1=1 # add to minlevel for fitting position