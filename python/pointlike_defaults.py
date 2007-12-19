#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.6 2007/11/27 04:35:33 burnett Exp $


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
    maxstep = 0.2 # max step allowed during localization: abort if larger


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

