#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/finder_setup.py,v 1.1 2007/07/19 13:52:13 burnett Exp $
print 'setup'
# data selection parameters

radius = 180   # radius in degrees for initial data selection
event_class = -1 # 0, select front only; -1 no selection
               # note that current algorithm only works on front-conversion events
source_id =-1  # -1: all sources -- select according to Monte Carlo source id, if present


verbose = 0  # set non-zero to get lots of output

#  specify files with FT1 data 

path = "F:/glast/data/SourceDetection"  #location on glast-ts
files = [path+"/pl_0_events_0000.fits" ,path+"/bg_low_0_events_0000.fits"]
# parameters.

# direction and cone or radius about it to examine
ra=0
dec=0
radius   = 180 # 180 for all sky

pix_level = 8  # pixel level to search
count_threshold=20 # 

TSmin   = 10   # minimum TS for candidates

# parameter for pruning: radius in degrees
prune_radius = 0.25

# file to write results to
outfile='sources.txt'
        
