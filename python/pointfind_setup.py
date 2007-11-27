#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.5 2007/11/19 01:09:36 burnett Exp $
print 'runing setup for pointfind'

# data source 

def test():
  " define the pixelfile for a quick test, running the pixel file created by the test program"
  import os
  path = os.environ['POINTLIKEROOT']
  return os.path.join(path, 'src', 'test', 'pointlike_test.fits')

pixelfile = test()
#files = [r'f:\glast\data\DC2\downlink\downlink_0169.fits']

class Diffuse:
    # diffuse input file image file
    file = ''
    exposure=3e10
    import os
    if 'EXTFILESSYS' in os.environ:
        file = os.path.join(os.environ['EXTFILESSYS'],'galdiffuse','GP_gamma.fits')

# data selection parameters
class Data:
    radius = 7.0   # radius in degrees for initial data selection
    event_type = -1 # 0, select front only; -1 no selection
    source_id =-1  # -1: all sources -- select according to Monte Carlo source id, if present



#  specify files with FT1 data and points to fit. pixelfile for PhotonMap, files for FT1 or merit
# files = [] 
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'

#------------------------------------------------------------------
# Common parameters
#------------------------------------------------------------------
class PointSourceLikelihood: #parameters for the likelihood calculation
    # HEALpix level range for energy band  fits

    minlevel=8   # minimum level to use for fits (>-6)  
    maxlevel=13  # maximum level for fits  (<=13)
    minalpha=0.15# minimum value for signal fraction to use in TS total

    # parameters governing iteration cycle of bands/position

    TSmin= 5     # minimum TS value to allow during iteration
    skip1=1      # inital number of layers to skip in localization fit
    skip2=4      # don't skip beyond this
    itermax=1    # maximum number of iterations

    verbose = 0  # set non-zero to get lots of output


class SourceFinder:  # parameters for the SourceFinder.

    TSmin = 10           # overall minimum TS
    skipTSlevels = 2     # number of lowest energy levels to skip when localizing


    pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
    countThreshold = 32  # minimum count (weighted, with children) to consider examining location further
    finalPixlvl = 8      # after source is localized, pixel level to check for sufficient photons
    finalCntThresold = 2 # after source is localized, min nbr photons (with children, NOT weighted) required in surrounding pixel
    prune_radius =1.0  # pruning radius in degrees (also needs to be tuned)

    countThreshold = 346  # override for 1-year obssim
            
    # direction and cone or radius about it to examine: either l,b, or ra,dec
    #l,b = 0, 0

    ra, dec = 248.79, -17.86 # 3EGJ1635m1751 
    radius   = 7 #  179.999 # 180 for all sky

    # file to write a table of results to
    outfile='../outputpointfind_3EGJ1635m1751.txt'
