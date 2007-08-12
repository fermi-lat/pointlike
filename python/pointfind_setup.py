#  setup for pointlike source finder
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/finder_setup.py,v 1.2 2007/07/19 14:07:11 burnett Exp $
print 'setup for pointfind'

# data source 

pixelfile = "F:/glast/data/SC2/obssim/sc2_obssim_map.fits"

# direction and cone or radius about it to examine
l,b = 0, -90

radius   = 60 # 180 for all sky

count_threshold=346  

TSmin   = 10   # minimum TS for candidates

# parameter for pruning: radius in degrees
prune_radius = 0.25

# file to write a table of results to
outfile='south_sources.txt'
        
