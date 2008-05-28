# default parameters for the various parameter files
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointlike_defaults.py,v 1.6 2008/02/19 21:00:33 burnett Exp $
#
# Include this to set defaults, then override 
import sys
print 'runing %s' % sys.argv[0]

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


class Diffuse:
    file = ''      # diffuse input file image file: if not specified, try to find it below
    exposure=3e10  # this seems appropriate for 1 year.
    import os
    if 'GLAST_EXT' in os.environ and file=='':
		file = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'GP_gamma.fits')


class PointSourceLikelihood: #parameters for the likelihood calculation
    # HEALpix level range for energy band  fits

    emin = 700;   # only select bands including or above this energy
    minalpha=0.15# minimum value for signal fraction to use in TS total

    # parameters governing iteration cycle of bands/position

    TSmin= 5     # minimum TS value to allow during iteration
    skip1=1      # inital number of layers to skip in localization fit
    skip2=4      # don't skip beyond this
    itermax=1    # maximum number of iterations

    verbose = 0  # set non-zero to get lots of output
    maxstep = 0.2 # max step allowed during localization: abort if larger
    
    # values for the gamma and sigma PSF parameters, indexed by level
    gamma_list =[0,0,0,0,0,
           2.25,  2.27,  2.22,  2.31,  2.30,  2.31,  2.16,  2.19,  2.07]
    sigma_list =[0,0,0,0,0,
           0.343, 0.335, 0.319, 0.431, 0.449, 0.499, 0.566, 0.698, 0.818]

class SourceFinder:  # parameters for the SourceFinder.

    TSmin = 8.0          # overall minimum TS
    pixel_fraction = 1.0 # fraction of pixels to sample, sorted according to weighted photon count
    prune_radius =0.25   # pruning radius in degrees (also needs to be tuned)
    group_radius = 4.0   # maximum radius for nearby sources.
    examine_radius=180   # default is to examine full sky
    pass1_nside=256      #
            

