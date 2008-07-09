#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfind_setup.py,v 1.18 2008/07/07 21:57:52 burnett Exp $

from pointlike_defaults import *

suffix='v2' # digel 

Data.LATalignment=[-186,-164, -540]  # from Marshall


# modify exposure if not a year
Diffuse.exposure*=4/365. #2.5e+009

#specify data: pixels or a list of FT1 files


# specify pixelfile (BinnedPhotonData) if exists, use it: otherwise generate
pixelfile = r'D:\common\first_light\binned_source_'+suffix+'.fits'

import os
if os.path.exists(pixelfile):
  Data.pixelfile = pixelfile
elif suffix[0]=='v':
  datapath = r'F:\glast\data\first_light\digel\\'
  Data.files= [datapath+'ft1_first_src_'+suffix+'.fits']
  Data.history=datapath+'ft2_first_'+suffix+'.fits'
  print 'Using alignment: %s' % Data.LATalignment
  Data.output_pixelfile = pixelfile
else:
  from runfiles import RunFiles
  datapath=r'f:/glast/downloads/'
  Data.history=r'd:\common\first_light\ft2\merged_'+suffix+'.fits'
  print 'Using alignment: %s' % Data.LATalignment
  runlist = r'D:/common/first_light/nomsciops_runs.txt'
  Data.files= RunFiles(datapath, runlist)('ph')
  Data.output_pixelfile = pixelfile


# choose region to search
SourceFinder.l,SourceFinder.b= 0,0
SourceFinder.examine_radius = 180 # 180 for all sky


# files to write a table of results to
path = r'd:\common\first_light/'
SourceFinder.outfile=  path+'pointfind_'+suffix+'.txt'
SourceFinder.regfile = path+'pointfind_'+suffix+'.reg'
SourceFinder_regtsmin = 20
SourceFinder.logfile = path+'pointfindlog_'+suffix+'.txt' # the log file

print 'will write to file %s '% SourceFinder.outfile
SourceFinder.group_radius = 2.0
SourceFinder.TSmin = 15

imagefile= path+'image_'+suffix+'.fits'
if not os.path.exists(imagefile):
  print 'Creating FITS image file at %s' % imagefile
  SourceFinder.imagefile=imagefile
  SourceFinder.imageresolution=0.1


print 'SourceFinder.TSmin: %s, emin: %s ' %(SourceFinder.TSmin, PointSourceLikelihood.emin)

# this function, if it exists, will be called at the end of the job
def finish():
  pass  