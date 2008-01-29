# default parameters for the various parameter files
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointlike_defaults.py,v 1.4 2008/01/25 01:07:06 burnett Exp $
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
    history = ''    # optional history or FT2 file, needed to correct for misalignment if readign FT1
    Latalignment=[] # alignment correction angles about x,y,z axes, in arcseconds


class Diffuse:
    file = ''      # diffuse input file image file: if not specified, try to find it below
    exposure=3e10  # this seems appropriate for 1 year.
    import os
    if 'GLAST_EXT' in os.environ and file=='':
		file = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'GP_gamma.fits')

#    if 'EXTFILESSYS' in os.environ and file!='':
#        file = os.path.join(os.environ['EXTFILESSYS'],'galdiffuse','GP_gamma.fits')


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

    TSmin = 8.0          # overall minimum TS
    pixLevel = 8         # pixelization level for first pass.  Higher number -> more points tested for possible sources
    pixel_fraction = 0.2 # fraction of pixels to sample, sorted according to weighted photon count
    group_radius = 1.0   # radius for making groups
    prune_radius =0.25   # pruning radius in degrees (also needs to be tuned)
    examine_radius=180   # default is to examine full sky
            

