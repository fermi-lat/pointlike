#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.12 2008/06/29 01:53:19 burnett Exp $
from  pointlike_defaults import *


# specify sources as dictionary
sources = { #  ra        dec         TS (optional: used to sort)
'vela'    :  (128.8359, -45.1763,   349.1),
'Geminga' :  ( 98.476,   17.770,    287.6), 
'psr1706-44':(257.428,  -44.486,     44.0), 
'crab'  :    ( 83.633,   22.014,    140.7),
'3C454.3':   (343.491,   16.148,     99.4),
}

PointSourceLikelihood.verbose=0

filelist = r'f:\glast\data\first_light\filelist.txt'
Data.files = [line.strip() for line in open(filelist)]

#Data.output_pixelfile = r'F:\glast\data\first_light\binned_ph.fits'
Data.pixelfile = r'F:\glast\data\first_light\binned_ph.fits'