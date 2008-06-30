#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.13 2008/06/29 21:07:17 burnett Exp $
from  pointlike_defaults import *


# specify sources as dictionary
sources = { #  ra        dec         TS (optional: used to sort)
'vela'    :  (128.8359, -45.1763,   349.1),
'Geminga' :  ( 98.476,   17.770,    287.6), 
'psr1706-44':(257.428,  -44.486,     44.0), 
'crab'    :  ( 83.633,   22.014,    140.7),
'3C454.3' :   (343.491,   16.148,     99.4),
'RX J1836.2+5925':  (279.057179, +59.425014 , 43),  # EGRET location: 278.87, 59.32 , 0.15 error circle
'3C 66A'  : (035.665048, +43.035500, 20),
'Mrk 421' : (166.113808, +38.208833, 10),
'LSI+61+303' : (40.13194,  61.22933, 6),  
}

PointSourceLikelihood.verbose=0


from runfiles import RunFiles
datapath=r'f:/glast/downloads/'
runlist = r'f:/glast/data/first_light/runlist.txt'
Data.files= RunFiles(datapath, runlist)('ph')

pixelfile = r'F:\glast\data\first_light\binned_source.fits'
#Data.output_pixelfile = pixelfile
Data.pixelfile = pixelfile