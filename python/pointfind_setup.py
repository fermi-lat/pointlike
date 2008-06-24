#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.12 2008/05/28 21:41:52 burnett Exp $

from pointlike_defaults import *

# modify exposure if not a year
#Diffuse.exposure*=1/12. #2.5e+009

#specify data: pixels or a list of FT1 files
Data.pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'

# choose region to search
SourceFinder.l,SourceFinder.b= 0,0
SourceFinder.examine_radius = 180 # 180 for all sky


# files to write a table of results to
path = '../output/obssim2/'
SourceFinder.outfile=  path+'pointfind_f.txt'
SourceFinder.regfile = path+'pointfind_f.reg'
SourceFinder.logfile = path+'pointfindlog_f.txt' # the log file

print 'will write to file %s '% SourceFinder.outfile
SourceFinder.group_radius = 2.0

print 'SourceFinder.TSmin: %s, emin: %s ' %(SourceFinder.TSmin, PointSourceLikelihood.emin)