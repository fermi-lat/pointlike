#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.7 2007/12/19 03:38:41 burnett Exp $
print 'runing setup for pointfind'

from pointlike_defaults import *

class SourceFinder:  # parameters for the SourceFinder.

    TSmin = 10           # overall minimum TS
    skipTSlevels = 2     # number of lowest energy levels to skip when localizing
    pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
    pixel_fraction = 0.4 # fraction of pixels to sample, sorted according to weighted photon count
    group_radius = 1.0   # radius for making groups
    prune_radius =0.25   # pruning radius in degrees (also needs to be tuned)

    ra, dec = 248.79, -17.86 # 3EGJ1635m1751 
    examine_radius   = 7# 179.999 # 180 for all sky

#specify data
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
PointSourceLikelihood.verbose = 0
# file to write a table of results to
outfile='../output/pointfind_3EGJ1635m1751_test.txt'
#outfile='../output/all_sky8.txt'
