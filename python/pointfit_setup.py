#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.2 2007/06/25 20:59:25 burnett Exp $

# data selection parameters

radius = 7.0   # radius in degrees for initial data selection
event_type = -1 # 0, select front only; -1 no selection
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

# test data from the pattern
path = "F:/glast/data/SourceDetection"  #location on glast-ts
files = [path+"/pl_0_events_0000.fits" ,
         path+"/bg_low_0_events_0000.fits"]

# this list is from  the sourcelist.txt, the strongest, hardest ones
sourceinfo="""\
Source_15       1.6     1.1643E-4      207.71  -30
Source_16       1.6     1.5527E-4      214.64  -30 
Source_17       1.6     2.0705E-4      221.57  -30 
Source_18       1.6     2.7611E-4       228.5  -30
Source_19       1.6     3.6820E-4      235.43  -30
Source_20       1.6     4.9100E-4      242.35  -30
Source_21       1.6     6.5476E-4      249.28  -30"""

name=[]; ra=[]; dec=[]
for line in sourceinfo.split('\n'):
   fields = line.split()
   name.append( fields[0])
   ra.append(float(fields[3]))
   dec.append(float(fields[4]))
  
