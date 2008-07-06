#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.15 2008/07/01 23:40:58 burnett Exp $
from  pointlike_defaults import *

suffix='05'

Data.LATalignment=[-1.9*60,-2.6*60, -8.6*60]  # from Marshall
Data.history=r'd:\common\first_light\ft2\merged_'+suffix+'.fits'
print 'Using alignment: %s' % Data.LATalignment

sourcelistfile = r'D:/common/first_light/sourcelist.txt'

PointSourceLikelihood.merge=1
verbose=0

# specify pixelfile (BinnedPhotonData) if exists, use it: otherwise generate
pixelfile = r'D:\common\first_light\binned_source_'+suffix+'.fits'

import os
if os.path.exists(pixelfile):
  Data.pixelfile = pixelfile
else:
  from runfiles import RunFiles
  datapath=r'f:/glast/downloads/'
  runlist = r'D:/common/first_light/leo_runs.txt'
  Data.files= RunFiles(datapath, runlist)('ph')
  Data.output_pixelfile = pixelfile
  
outfile = r'D:/common/first_light/pointfit_'+suffix+'.txt'