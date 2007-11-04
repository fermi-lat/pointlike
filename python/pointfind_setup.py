#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.3 2007/09/13 22:28:43 burnett Exp $
print 'runing setup for pointfind'

# data source 

def test():
  " define the pixelfile for a quick test, running the pixel file created by the test program"
  import os
  path = os.environ['POINTLIKEROOT']
  return os.path.join(path, 'src', 'test', 'pointlike_test.fits')

pixelfile = test()
#files = [r'f:\glast\data\DC2\downlink\downlink_0169.fits']

# diffuse input file image file
diffusefile = ''
tolerance = 1.0     # average is the point in the center
import os
if 'EXTFILESSYS' in os.environ:
    diffusefile = os.path.join(os.environ['EXTFILESSYS'],'galdiffuse','GP_gamma_v0r0p1.fits')
    print 'using diffuse definition from file %s, with integral tolerance %f' % (diffusefile, tolerance)
else:
    print 'not using diffuse file'

#  specify files with FT1 data and points to fit. pixelfile for PhotonMap, files for FT1 or merit
# files = [] 
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'

#------------------------------------------------------------------
# Common parameters
#------------------------------------------------------------------

# direction and cone or radius about it to examine
#l,b = 0, 0

ra, dec = 83.6, 22.014
ra, dec = 257.5, -44.6 # PSF_B1706m44 
ra, dec = 249, -18    #HLCloud_SC2_05
radius   = 5 #  179.999 # 180 for all sky


# file to write a table of results to
#outfile='pointfind_Sc2_0930d.txt'

prune_radius = 0.25  # pruning radius in degrees (also needs to be tuned) nominal 0.25


eqTSmin   = 10       # minimum TS for candidates in equitorial region
midTSmin   = 10      # minimum TS for candidates in middle region
polarTSmin   = 10    # minimum TS for candidates in polar region
skipTSlevels = 2     # number of lowest energy levels to skip when localizing


pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
countThreshold = 32  # minimum count (weighted, with children) to consider examining location further
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

countThreshold = 346  # override for 1-year obssim
        
