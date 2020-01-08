#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/sourcefit_setup.py,v 1.2 2008/09/11 06:40:21 markusa Exp $
from  pointlike_defaults import *
from math import sin,asin,cos,atan2,pi

def source_description_lists(source):
    name, ra, dec, srctype, npar, init, dofit = [],[],[],[],[],[],[]
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
        if(len(source[key])<6): dofit.append(1.0)
	else: dofit.append(source[key][5])

    return name,ra,dec,srctype,npar,init, dofit    

def gal2equ(ll,bb):
      bb=bb*pi/180.
      ll=ll*pi/180.
      ra_gp = 192.85948  *pi/180.
      de_gp =  27.12825  *pi/180.
      lcp   = 122.932    *pi/180.
      
      sin_d = sin(de_gp)*sin(bb)+cos(de_gp)*cos(bb)*cos(lcp-ll)
      ramragp = atan2(cos(bb)*sin(lcp-ll),cos(de_gp)*sin(bb)-sin(de_gp)*cos(bb)*cos(lcp-ll))
      
      dec=asin(sin_d)
      ra=(ramragp+ra_gp+2*pi)%(2*pi)
      return  ra*180./pi,dec*180./pi

source={}
 
first_is_center=0

SourceLikelihood.verbose=1

def ext(deg):
   return str(deg*pi/180.);

print ("SRC: ",gal2equ(0.0,0.0) )
ra,dec = 30.0, 30.0 ; source['small_00'] =(ra-0.5,dec+0.5,'gauss',1,ext(0.01),1.)
#ra,dec = 30.0, 30.0 ; source['ppoint_00'] =(ra+0.5,dec+0.5,'pseudopoint',0,ext(0.0001),1.)
#ra,dec = 30.0, 30.0 ; source['point_00'] =(ra-0.5,dec-0.5,'point',0,ext(0.0001),1.)

#ra,dec = 60.0, 30.0 ; source['small_01'] =(ra-0.5,dec+0.5,'gauss',1,ext(0.01),1.)
#ra,dec = 60.0, 30.0 ; source['ppoint_01'] =(ra+0.5,dec+0.5,'pseudopoint',0,ext(0.0001),1.)
#ra,dec = 60.0, 30.0 ; source['point_01'] =(ra-0.5,dec-0.5,'point',0,ext(0.0001),1.)

#ra,dec = 90.0, 30.0 ; source['small_02'] =(ra-0.5,dec+0.5,'gauss',1,ext(0.01),1.)
#ra,dec = 90.0, 30.0 ; source['ppoint_02'] =(ra+0.5,dec+0.5,'pseudopoint',0,ext(0.0001),1.)
#ra,dec = 90.0, 30.0 ; source['point_02'] =(ra-0.5,dec-0.5,'point',0,ext(0.0001),1.)

#ra,dec = 120.0, 30.0 ; source['small_03'] =(ra-0.5,dec+0.5,'gauss',1,ext(0.01),1.)
#ra,dec = 120.0, 30.0 ; source['ppoint_03'] =(ra+0.5,dec+0.5,'pseudopoint',0,ext(0.0001),1.)
#ra,dec = 120.0, 30.0 ; source['point_03'] =(ra-0.5,dec-0.5,'point',0,ext(0.0001),1.)

#ra,dec = 150.0, 30.0 ; source['small_04'] =(ra-0.5,dec+0.5,'gauss',1,ext(0.01),1.)
#ra,dec = 150.0, 30.0 ; source['ppoint_04'] =(ra+0.5,dec+0.5,'pseudopoint',0,ext(0.0001),1.)
#ra,dec = 150.0, 30.0 ; source['point_04'] =(ra-0.5,dec-0.5,'point',0,ext(0.0001),1.)

#ra,dec = 180.0, 30.0 ; source['small_05'] =(ra-0.5,dec+0.5,'gauss',1,ext(0.01),1.)
#ra,dec = 180.0, 30.0 ; source['ppoint_05'] =(ra+0.5,dec+0.5,'pseudopoint',0,ext(0.0001),1.)
#ra,dec = 180.0, 30.0 ; source['point_05'] =(ra-0.5,dec-0.5,'point',0,ext(0.0001),1.)

name,ra,dec,srctype,npar,init,fit = source_description_lists(source)

Data.files = [
r'/u/gl/markusa/ki-disk/data/pointlike/adam_gaussian_src_events_0000.fits']
#r'/u/gl/markusa/ki-disk/data/pointlike/extragalactic/egal_LAT_allsky_1year.000_V01.fits'
#]

Diffuse.exposure = 1.
Diffuse.file='/u/gl/markusa/ki-disk/data/pointlike/uniform_background.fits.gz'


SourceLikelihood.UseMinuit=1
SourceLikelihood.UseSimplex=0
SourceLikelihood.UseMinos=1
SourceLikelihood.UseGradient=1

SourceLikelihood.umax=100.
SourceLikelihood.minalpha=0.

SourceLikelihood.emin=500
SourceLikelihood.gammaFront=[ 1.785, 1.874, 1.981, 2.003, 2.225, 2.619, 2.527, 2.134, 1.827 ]
SourceLikelihood.sigmaFront=[ 1.77e+00, 1.06e+00, 5.69e-01, 2.85e-01, 1.54e-01, 9.16e-02, 5.15e-02, 2.52e-02, 1.55e-04]
SourceLikelihood.sigmaFront=[x*pi/180. for x in SourceLikelihood.sigmaFront]
SourceLikelihood.gammaBack =[2., 1.737, 1.742, 1.911, 2.008, 2.009, 2.207, 1.939]
SourceLikelihood.sigmaBack =[4., 2.18e+00, 1.17e+00, 6.00e-01, 3.09e-01, 1.61e-01, 8.43e-02, 3.90e-02 ]
SourceLikelihood.sigmaBack=[x*pi/180. for x in SourceLikelihood.sigmaBack]


outfile="gaussian_fit_results_bins_per_decade.fits"
bgROI  = 0.001 # rad

Data.output_pixelfile="test.fits"

#name= ['w28']
#dir = (270.426, -23.4)
#dec=[-45.16]
#ra=[128.835939]; dec=[-45.177672] # best level-10 fit?
