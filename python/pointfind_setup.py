#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.6 2007/11/27 04:35:33 burnett Exp $
print 'runing setup for pointfind'

from pointlike_defaults import *

class SourceFinder:  # parameters for the SourceFinder.

    TSmin = 10           # overall minimum TS
    skipTSlevels = 2     # number of lowest energy levels to skip when localizing


    pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
    countThreshold = 32  # minimum count (weighted, with children) to consider examining location further
    finalPixlvl = 8      # after source is localized, pixel level to check for sufficient photons
    finalCntThresold = 2 # after source is localized, min nbr photons (with children, NOT weighted) required in surrounding pixel
    prune_radius =1.0  # pruning radius in degrees (also needs to be tuned)

    countThreshold = 3000  # override for 1-year obssim
            
    # direction and cone or radius about it to examine: either l,b, or ra,dec
    #l,b = 0, 0

    ra, dec = 248.79, -17.86 # 3EGJ1635m1751 
    radius   = 7# 179.999 # 180 for all sky

#specify data
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
PointSourceLikelihood.verbose = 0
# file to write a table of results to
outfile='../output/pointfind_3EGJ1635m1751_test.txt'
#outfile='../output/all_sky8.txt'
