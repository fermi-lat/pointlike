#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.11 2008/06/06 17:01:59 burnett Exp $
from  pointlike_defaults import *


#Data.pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
import glob

"""
Data.files = glob.glob(r'F:\glast\data\SC2\obssim\LAT_allsky_*_V01.fits')[:3]
Data.event_class=-1
Data.bins_per_decade=4
from numpy import arange
Data.energy_bins=10.**arange(1.5,5,0.25);
"""

# specify sources as dictionary
sources = { #  ra        dec         TS (optional: used to sort)
'vela'    :  (128.8359, -45.1763,   349.1),
'Geminga' :  ( 98.476,   17.770,    287.6), 
'psr1706-44':(257.428,  -44.486,     44.0), 
'crab'  :    ( 83.633,   22.014,    140.7),
'3C454.3':   (343.491,   16.148,     99.4),
}

PointSourceLikelihood.skip1=1 # add to minlevel for fitting position
PointSourceLikelihood.verbose=0

# a chunk of the FT1 data
# Data.files=glob.glob('f:\glast\downloads\gll_ph_r02363*.fit')
#Data.files=glob.glob('f:\glast\downloads\gll_ph*.fit')
filelist = r'f:\glast\data\first_light\filelist.txt'
Data.files = [line.strip() for line in open(filelist)]
print Data.files
Data.output_pixelfile = r'F:\glast\data\first_light\binned_ph.fits'