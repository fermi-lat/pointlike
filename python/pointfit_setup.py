#  setup for point fit test
# $Header$

# data selection parameters

radius = 7.0   # radius in degrees for initial data selection
event_type = 0 # 0, select front only; -1 no selection
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
files = [path+"/pl_0_events_0000.fits" ,path+"/bg_low_0_events_0000.fits"]

points = [
        ["Source_54",   186.31,  -18],
        ["Source_55",   192.62,  -18],
        ["Source_56",   198.93,  -18],
        ["Source_57",   205.24,  -18],
        ["Source_58",   211.54,  -18],
        ["Source_59",   217.85,  -18],
        ["Source_60",   224.16,  -18],
        ["Source_61",   230.47,  -18]]

name  = [point[0] for point in points]
ra    = [point[1] for point in points]
dec   = [point[2] for point in points]
