#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.2 2007/08/31 04:24:41 burnett Exp $
print 'runing setup for pointfind'

# data source 

def test():
  " define the pixelfile for a quick test, running the pixel file created by the test program"
  import os
  path = os.environ['POINTLIKEROOT']
  return os.path.join(path, 'src', 'test', 'pointlike_test.fits')

pixelfile = test()

# file to write a table of results to
outfile='pointfind_test.txt'

#------------------------------------------------------------------
# Common parameters
#------------------------------------------------------------------

# direction and cone or radius about it to examine
l,b = 0, 0
radius   = 180       # 180 for all sky

eqTSmin   = 10       # minimum TS for candidates in equitorial region
midTSmin   = 10      # minimum TS for candidates in middle region
polarTSmin   = 10    # minimum TS for candidates in polar region
skipTSlevels = 2     # number of lowest energy levels to skip when localizing

pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
countThreshold = 16  # minimum count (weighted, with children) to consider examining location further
finalPixlvl = 8      # after source is localized, pixel level to check for sufficient photons
finalCntThresold = 2 # after source is localized, min nbr photons (with children, NOT weighted) required in surrounding pixel

eqBoundary = 6       # abs(b)in degrees for equatorial region < this number
polarBoundary = 40   # abs(b)in degrees for polar region > this number

plEqTSmin   = 38     # \
plMidTSmin   = 19    #  |-TS above this number is considered good, regardless of power law fit
plPolarTSmin   = 24  # /

plSlopeCutoff = -1.5 # power law slope <= this is good candidate
plFitCutoff = 0.9    # power law fit confidence factor >= this is good candidate

prune_radius = 0.25  # pruning radius in degrees (also needs to be tuned)


        
