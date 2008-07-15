#  setup for point fit 
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.18 2008/07/09 04:01:46 burnett Exp $
from  pointlike_defaults import *
import os

# specific for this analysis
analysis_path =r'D:/common/first_light/'
datapath=r'f:/glast/downloads/'
suffix='v2d'

sourcelistfile = analysis_path+'sourcelist.txt'
sourcelistfile = analysis_path+'associated_sources.txt'

PointSourceLikelihood.merge=1
verbose=0

# specify pixelfile (BinnedPhotonData) if exists, use it: otherwise generate
Data.pixelfile = analysis_path+'binned_source_'+suffix+'.fits'

if not os.path.exists(Data.pixelfile):
  print 'Pixel file not found: run pointfind first'
  raise Exception
  
  
outfile = analysis_path+'pointfit_'+suffix+'.txt'

# this function, if it exists, will be called at the end of the job
def finish():
  print 'Finishing output'
  if not os.path.exists(outfile):
    print 'job failed? no output'
    return
  from confluence import Table # module in the path
  t=Table(analysis_path, suffix, link_sed=False) 
  print 'wrote confluence-style table to %s' % t.outfilename