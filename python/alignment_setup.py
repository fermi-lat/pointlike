#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/alignment_setup.py,v 1.4 2008/06/17 05:26:04 mar0 Exp $

from  pointlike_defaults import *
import glob

# setup the data file, including setting alignment
Data.files = glob.glob(r'F:\glast\data\55-day\diffuse*.root') # file to test for alignment     252460800
#Data.files = [r'F:\glast\data\SC2\OneYrAllSky_startingApril\allSky_Month_01_events.fits'] # starts at 260238030
#Data.history = r'F:\glast\data\SC2\OneYrAllSky_startingApril\Gleam_survey_orbit8_long.fits'
#Data.files = [r'F:\glast\data\SC2\interleave\Interleave_pruned.root']
#Data.files = glob.glob(r'F:\glast\data\obssim3\combined_day09*.fits')
#Data.history = r'F:\glast\data\obssim3\FT2_obssim3.fits'
Data.files = glob.glob(r'f:/glast/downloads/root/d*.root')
#Data.pixelfile = r'f:/glast/downloads/root/84diffuse.fits'
#Data.files = [r'F:/glast/downloads/root/GRB.root']
#Data.pixelfile = r'D:\Users\kerrm\Python\analyses\J2021\j2021_off_75_91_bands.fits'

verbose = 0

Data.energy_bins=[]
Data.bins_per_decade=4
#Data.event_class=0
#Data.LATalignment=[72,0,-36]    # LAT alignment in arcsec
#Data.LATalignment=[0,0,0]
#PointSourceLikelihood.merge=0
#PointSourceLikelihood.emin=101
#PointSourceLikelihood.emax=56235
#PointSourceLikelihood.umax=25
dopsf=0
#verbose=1
print 'LAT alignment angles: %s' % Data.LATalignment

day=1+2*14 # weeks 5,6
day=1+3*14 # weeks 7,8

#day, ndays=57, 7 # week 9
#day, ndays=50, 7 # week 8
#day, ndays=43, 7 # week 7
#day, ndays=36, 7 # week 6
#day, ndays=29, 7 # week 5
#day, ndays=22, 7 # week 4
#day, ndays=15, 7 # week 3
#day, ndays=8, 7 # week 2
day, ndays=65, 1  # week 1


def set_days(day, nday=1):
	seconds=86400
	tzero=242053458
	tzero=236511638
	Data.start_time = tzero+(day-1)*seconds # start time in GLAST time format (252460800 for test file)
	Data.stop_time =  Data.start_time+ nday*seconds #  stop time in GLAST time format
	print 'set days %d-%d' % (day, day+nday-1)

set_days(day, ndays)
usefitdirection = False
separation = 0	# minimum separation distance of sources
maxerr = 96			# choose maximum localization error in degrees or return code, 
						# error codes: 97 - step error, 98 - lost, 99 - convergence
print 'fit direction: %s' %usefitdirection
#------------------------------------------------------------------
# Alignment parameters
#------------------------------------------------------------------
class Alignment:
	resolution = 20			# chooses the spacing of the grid points in arcseconds
	outfile = 'par2.txt'    # output file
	

#------------------------------------------------------------------
# Source position input
#------------------------------------------------------------------
sourcelist = r'D:/common/sourcelikefilter/sourcelists/bright_point_sources.txt'
sourcelist = r'D:/common/sourcelikefilter/sourcelists/associated_sources_aug.txt'
#sourcelist = r'D:/common/sourcelikefilter/sourcelists/pointfind_sources.txt'
#sourcelist = r'D:/common/sourcelikefilter/sourcelists/cgrabs_sorted.txt'
#sourcelist = r'D:/common/sourcelikefilter/sourcelists/combined_pulsar_cgrabs.txt'
#sourcelist = r'D:\common\sourcelikefilter\pointlike\v6r6p1\src\alignment\vela.txt'
print 'reading sources from file %s' % sourcelist
sourceinfo = file(sourcelist)
headings = sourceinfo.readline()
name_list=[]; ra_list=[]; dec_list=[]
for line in sourceinfo:
   fields = line.split()
   name_list.append( fields[0])
   ra_list.append(float(fields[1]))
   dec_list.append(float(fields[2]))