#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.1.1.1 2007/06/14 18:30:14 burnett Exp $

# data selection parameters

radius = 7.0   # radius in degrees for initial data selection
event_type = -1 # 0, select front only; -1 no selection
               # note that current algorithm only works on front-conversion events
source_id =-1  # -1: all sources -- select according to Monte Carlo source id, if present


# HEALpix level range for energy band  fits

minlevel=8   # minimum level to use for fits (>-6)  
maxlevel=13  # maximum level for fits  (<=13)

# parameters governing iteration cycle of bands/position

TSmin= 5     # minimum TS value to allow during iteration
skip1=1      # inital number of layers to skip in localization fit
skip2=4      # don't skip beyond this
itermax=2    # maximum number of iterations

verbose = 0  # set non-zero to get lots of output

#  specify files with FT1 data and points to fit.

path = "F:/glast/data/SourceDetection"  #location on glast-ts
#files = [path+"/pl_0_events_0000.fits" ,path+"/bg_low_0_events_0000.fits"]
files = ["F:/condor/alignment/scanning/jobs100-149/tuple.root"]
#,
#        "F:/condor/alignment/scanning/jobs200-299/tuple.root",
#]

        
points=[['vela',   128.73, -45.2 ]]
#        ['crab',    83.57,  22.01],
#        ['geminga', 98.49,  17.86]]

# program expects separate lists
name  = [point[0] for point in points]
ra    = [point[1] for point in points]
dec   = [point[2] for point in points]

