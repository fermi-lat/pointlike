# default parameters for the various parameter files
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointlike_defaults.py,v 1.2 2007/12/19 18:27:58 burnett Exp $
#
# Include this to set defaults, then override 

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

    TSmin = 8.0           # overall minimum TS
    pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
    pixel_fraction = 0.2 # fraction of pixels to sample, sorted according to weighted photon count
    group_radius = 1.0   # radius for making groups
    prune_radius =0.25   # pruning radius in degrees (also needs to be tuned)
            

