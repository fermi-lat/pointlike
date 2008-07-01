#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.14 2008/06/30 23:38:21 burnett Exp $
from  pointlike_defaults import *


# specify sources as dictionary
sources =    { #  ra        dec         TS (optional: used to sort)
'vela'       :  (128.8359, -45.1763,   349.1),
'Geminga'    :  ( 98.476,   17.770,    287.6), 
'psr1706-44' : (257.428,  -44.486,     44.0), 
'crab'       :  ( 83.633,   22.014,    140.7),
'3C454.3'    :  (343.491,   16.148,     99.4),
'RX J1836.2+5925':  (279.057179, +59.425014 , 43),  # EGRET location: 278.87, 59.32 , 0.15 error circle
'3C 66A'     : (035.665048, +43.035500, 20),
'LSI+61+303' : (40.13194,  61.22933, 6), 
'PKS 2155-304' :  (329.7169379,	-30.2255883, 37),
'GRO J0004+73' :  (002.4025, +73.1825, 42), # 42 flux, 1.85
#'3EG J1512-0849' : (228.17 ,-08.83, 10), #Gino claim
'Mrk 421'    : (166.113808, +38.208833, 10),
'x'          : (305.19,      36.95,     36),
'x2'         : (305.40,      40.50,     38),
'3EG J2021+3716    ': (305.30,      37.27,     38), 0.3 deg
'1RXS J184643.9+670636'      : (281.682921, +67.110139,     39),

}

PointSourceLikelihood.verbose=0


from runfiles import RunFiles
datapath=r'f:/glast/downloads/'
runlist = r'f:/glast/data/first_light/runlist.txt'
Data.files= RunFiles(datapath, runlist)('ph')

pixelfile = r'F:\glast\data\first_light\binned_source_01.fits'
#Data.output_pixelfile = pixelfile
Data.pixelfile = pixelfile