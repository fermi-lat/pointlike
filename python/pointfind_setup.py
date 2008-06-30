#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.14 2008/06/29 01:53:48 burnett Exp $

from pointlike_defaults import *

# modify exposure if not a year
Diffuse.exposure*=1/365. #2.5e+009

#specify data: pixels or a list of FT1 files
Data.pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
Data.pixelfile = r'F:\glast\data\obssim3_v2\uw\pixels_all_days.fits'
Data.pixelfile = r'F:\glast\data\first_light\binned_ph.fits'

pixelfile = r'F:\glast\data\first_light\binned_source.fits'
Data.pixelfile = pixelfile

# choose region to search
SourceFinder.l,SourceFinder.b= 0,0
SourceFinder.examine_radius = 180 # 180 for all sky


# files to write a table of results to
path = r'F:\glast\data\first_light/'
suffix='source'
SourceFinder.outfile=  path+'pointfind_'+suffix+'.txt'
SourceFinder.regfile = path+'pointfind_'+suffix+'.reg'
SourceFinder.logfile = path+'pointfindlog_'+suffix+'.txt' # the log file

print 'will write to file %s '% SourceFinder.outfile
SourceFinder.group_radius = 2.0
SourceFinder.TSmin = 10
SourceFinder.imagefile= path+'image_source.fits'
SourceFinder.imageresolution=0.1


print 'SourceFinder.TSmin: %s, emin: %s ' %(SourceFinder.TSmin, PointSourceLikelihood.emin)