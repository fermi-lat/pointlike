#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.10 2008/01/25 22:35:14 burnett Exp $

from pointlike_defaults import *

# exposure appropriate for one month
#Diffuse.exposure*=1/12. #2.5e+009

#specify data
Data.pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'


# choose region to search
SourceFinder.ra, SourceFinder.dec = 248.79, -17.86 # 3EGJ1635m1751 
SourceFinder.ra, SourceFinder.dec = 137.5, 59.8 # DC2_3EGJ1048m5840 
SourceFinder.examine_radius   =  5# 180 for all sky

#Data.files=[r'F:\glast\data\SC2\obssim\month_test2a\merged_files_v1.fits']
#Data.output_pixelfile = 'test2a_photonmap.fits'

#Data.pixelfile='test2a_photonmap.fits'

# file to write a table of results to
#SourceFinder.outfile='../output/J1048m5840b.txt'
#print 'will write to file %s '% SourceFinder.outfile


SourceFinder.refit = 1  # turn on refit
SourceFinder.regfile = '../output/refit.reg'
#SourceFinder.fitsfile = '../output/refit.fits'
PointSourceLikelihood.minlevel=6