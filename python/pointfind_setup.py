#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.13 2008/06/24 19:06:02 burnett Exp $

from pointlike_defaults import *

# modify exposure if not a year
#Diffuse.exposure*=1/12. #2.5e+009

#specify data: pixels or a list of FT1 files
Data.pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
Data.pixelfile = r'F:\glast\data\obssim3_v2\uw\pixels_all_days.fits'

# choose region to search
SourceFinder.l,SourceFinder.b= 0,0
SourceFinder.examine_radius = 180 # 180 for all sky


# files to write a table of results to
path = '../output/obssim3/'
suffix='a'
SourceFinder.outfile=  path+'pointfind_'+suffix+'.txt'
SourceFinder.regfile = path+'pointfind_'+suffix+'.reg'
SourceFinder.logfile = path+'pointfindlog_'+suffix+'.txt' # the log file

print 'will write to file %s '% SourceFinder.outfile
SourceFinder.group_radius = 2.0

print 'SourceFinder.TSmin: %s, emin: %s ' %(SourceFinder.TSmin, PointSourceLikelihood.emin)