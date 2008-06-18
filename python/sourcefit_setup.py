#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.8 2008/01/25 22:36:11 burnett Exp $
from  pointlike_defaults import *

def source_description_lists(source):
    name, ra, dec, srctype, npar, init = [],[],[],[],[],[]
    allkeys=source.keys()
    for key in allkeys:
	name.append(key)
	ra.append(source[key][0])
	dec.append(source[key][1])
	srctype.append(source[key][2])
        if(len(source[key])<4): npar.append(0)
	else: npar.append(source[key][3])
        if(len(source[key])<5): init.append('')
	else: init.append(source[key][4])

    return name,ra,dec,srctype,npar,init    


#  specify files with FT1 data and points to fit. pixelfile for PhotonMap, files for FT1 or merit
def test():
  " define the pixelfile for a quick test, running the pixel file created by the test program"
  import os
  path = os.environ['POINTLIKEROOT']
  return os.path.join(path, 'src', 'test', 'pointlike_test.fits')

# data selection: either "pixelfile", or "files", the latter a list of FT1 files
#pixelfile = test()
#pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'

# if this is non-zero, lots of output

source={}

# if this is non-zero, use the first of the list as a background for the remainder

# the troublesome triplet: the 3EG blazar is very strong, affects the HLCloud
#name = ['DC2_3EGJ1635m1751', 'HLCloud_SC1_05', 'Giommi_blazar_1237', 'bogus1', 'bogus2']
#ra   = [248.788,    248.4804, 248.34, 248.51, 248.27]
#dec  = [-17.861,   -18.294,  -18.71  ,-17.88, -18.12]
 

#ra = [128.8359]; dec=[-45.1763]
#pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
#files = [r'F:\glast\data\SC2\interleave\Interleave_pruned.root']
#files = glob.glob(r'F:\glast\data\octobertest\merit\cutfile*.root')
#files = [r'F:\glast\data\octobertest\merit\cutfile8.root']
#files = glob.glob(r'F:\glast\data\octobertest\FT\*ft1.fits')

#Data.files=[r'F:\glast\data\55-day\Skimmer_pruned.root']
#Data.output_pixelfile = '55-day_pmap.fits'
#Data.pixelfile = '55-day_pmap.fits'
#PointSourceLikelihood.skip1=3 # try skipping first

# setup for selecting a range.
#tzero = 252460800
#week = 7*86400
#Data.start_time = tzero+7*week
#Data.stop_time = Data.start_time+week
#Data.output_pixelfile = 'week8_pmap.fits'
#Data.pixelfile='week1_pmap.fits'
#Data.source_id = 22000

first_is_center=0

SourceLikelihood.verbose=1

source['vela'] = (128.836673, -45.188701,'point',0)
#source['w28'] = (270.426, -23.4,'point',0)
#source['vela'] = (128.836673, -45.188701,'gauss',1,'0.01')
#source['vela'] = (128.836673, -45.188701,'pseudopoint')
#source['vela'] = (128.836673, -45.188701,'point')
#source['geminga'] = (98.476, 17.770,'point')
#source['crab'] = (83.633, 22.014,'point')

Data.files = [r'/u/gl/funk/ScienceData/sc2/week1_data.fits']
#Data.files = [r'/home/markusa/data/sc2/LAT_allsky_252460800.000_V01.fits']
name,ra,dec,srctype,npar,init = source_description_lists(source)

print name,ra,dec,srctype,npar,init[0]
Diffuse.exposure = 3.e10/52.
#Diffuse.file='/u/gl/markusa/ki-disk/data/pointlike/uniform_background.fits.gz'

SourceLikelihood.minlevel=8
SourceLikelihood.skip1=0 # add to minlevel for fitting position
#SourceLikelihood.NewStyle=1
#PointSourceLikelihood.maxlevel=13
SourceLikelihood.UseMinuit=1
outfile="sources.fits"

Data.output_pixelfile="test.fits"

#name= ['w28']
#dir = (270.426, -23.4)
#dec=[-45.16]
#ra=[128.835939]; dec=[-45.177672] # best level-10 fit?
