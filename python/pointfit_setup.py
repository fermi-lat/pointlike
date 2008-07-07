#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.16 2008/07/06 06:41:33 burnett Exp $
from  pointlike_defaults import *

# specific for this analysis
analysis_path =r'D:/common/first_light/'
datapath=r'f:/glast/downloads/'
suffix='06'

Data.LATalignment=[-1.9*60,-2.6*60, -8.6*60]  # from Marshall
Data.history=r'd:\common\first_light\ft2\merged_'+suffix+'.fits'
print 'Using alignment: %s' % Data.LATalignment

sourcelistfile = analysis_path+'sourcelist.txt'

PointSourceLikelihood.merge=1
verbose=0

# specify pixelfile (BinnedPhotonData) if exists, use it: otherwise generate
pixelfile = analysis_path+'binned_source_'+suffix+'.fits'

import os
if os.path.exists(pixelfile):
  Data.pixelfile = pixelfile
else:
  from runfiles import RunFiles
  runlist = analysis_path+'leo_runs.txt'
  Data.files= RunFiles(datapath, runlist)('ph')
  Data.output_pixelfile = pixelfile
  
outfile = analysis_path+'pointfit_'+suffix+'.txt'

# this function, if it exists, will be called at the end of the job
def finish():
  print 'Finishing output'
  if not os.path.exists(outfile):
    print 'job failed? no output'
    return
  from confluence import Table # module in the path
  t=Table(analysis_path, suffix) 
  print 'wrote confluence-style table to %s' % t.outfilename