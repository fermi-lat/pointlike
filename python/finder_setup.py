#  setup for pointlike source finder
# $Header$

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

#  specify files with FT1 data 

path = "F:/glast/data/SourceDetection"  #location on glast-ts
files = [path+"/pl_0_events_0000.fits" ,path+"/bg_low_0_events_0000.fits"]
# parameters.

eventType=-1
center   = SkyDir(0,0, SkyDir.GALACTIC) 
radius   = 180 # for all sky

prune_radius = 0.25

TSmin   = 10
eq_TS_min =15 #25.,
mid_TS_min= 15 #19.,
polar_TS_min= 15 #18.,
pix_level = 8
count_threshold=20 # was 12
skip_TS = 2

# file to write results to
tablefile='sources.txt'
        
