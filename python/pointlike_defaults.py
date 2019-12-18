# default parameters for the various parameter files
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointlike_defaults.py,v 1.16 2009/02/24 20:50:52 burnett Exp $
#
# Include this to set defaults, then override 
import sys
print ('running %s' % sys.argv[0])

# data selection parameters
class Data:
    pixelfile=''    # set to a photon data FITS file, or use a list of FT1 files, If set, ignore alignment, times, output
    files = []      # set to a list of FT1-like FITS files (if no pixefile)
    event_class = -1 # 0, select front only; -1 no selection
    source_id =-1   # -1: all sources -- select according to Monte Carlo source id, if present
    output_pixelfile = '' # set to create an output pixel file (if reading FT1 or ROOT files)
    start_time=0.   # select interval if non zero
    stop_time=0.    # "
    history = ''    # optional history or FT2 file, needed to correct for misalignment if reading FT1
    LATalignment=[] # alignment correction angles about x,y,z axes, in arcseconds
    energy_bins=[]  #


class Diffuse:
    file = ''      # diffuse input file image file: if not specified, try to find it below
    exposure=3e10  # this seems appropriate for 1 year.
    import os
    if 'GLAST_EXT' in os.environ and file=='':
        file = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'gll_iem_v01.fit')

class Isotropic:   # isotropic flux description
    flux = 1.5e-5
    index = 2.2
    
class Exposure:
    livetimefile=r'D:\common\pointfind\data\aug2008-jan2009_livetime.fits'
    IRF = 'P6_v1_diff'
    

class PointSourceLikelihood: #parameters for the likelihood calculation
    # HEALpix level range for energy band  fits

    emin = 200;   # only select bands including or above this energy
    minalpha=0.01 # minimum value for signal fraction to use in TS total
    minROI  = 0 # minimum size for ROI (deg)
    maxROI = 20   # prevent from bing too large

    # parameters governing iteration cycle of bands/position

    TSmin= 5     # minimum TS value to allow during iteration
    skip1=1      # inital number of layers to skip in localization fit
    skip2=9      # don't skip beyond this
    itermax=1    # maximum number of iterations

    verbose = 0  # set non-zero to get lots of output
    maxstep = 0.2 # max step allowed during localization: abort if larger
    merge=1     # set to zero to not attempt to merge bands with identical parameters
    

class SourceFinder:  # parameters for the SourceFinder.

    pass1_nside=256      # HEALpix binning for initial points.
    TSmin = 5.0          # overall minimum TS
    pixel_fraction = 1.0 # fraction of pixels to sample, sorted according to TS
    prune_radius =0.25   # pruning radius in degrees (also needs to be tuned)
    group_radius = 4.0   # maximum radius for nearby sources.
    examine_radius=180   # default is to examine full sky
            

