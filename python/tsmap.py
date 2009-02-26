import sys
import os
import re
import pointlike as pl
import skymaps as sm
import pyfits

from PeakFinder import PeakFinder
from optparse import OptionParser

gald_total = '/u/gl/markusa/ki-disk/skymaps/GP_gamma_expcounts_healpix_o8_54_59Xvarh8S.fits'
gald_front = '/u/gl/markusa/ki-disk/skymaps/GP_gamma_expcounts_front_healpix_o8_54_59Xvarh8S.fits'
gald_back  = '/u/gl/markusa/ki-disk/skymaps/GP_gamma_expcounts_back_healpix_o8_54_59Xvarh8S.fits'

iso_map    = '/u/gl/markusa/ki-disk/skymaps/extragalactic_lat_healpix_o7.pass6_irfs_5month.fits'

gald_CO      = '/u/gl/markusa/ki-disk/skymaps/CO/GP_gamma_expcounts_healpix_o8_COonly_54_59Xvarh8S.fits.gz'
COmap_galp   = '/u/gl/markusa/ki-disk/skymaps/CO/galprop_dame_co.mapcube.fits.gz'
COmap_dame   = '/u/gl/markusa/ki-disk/skymaps/CO/DAME_Wco_DHT2001.mapcube.fits.gz'
COmap_nanten = '/u/gl/markusa/ki-disk/skymaps/CO/NANTEN_all031220-to.mapcube.fits.gz'

evfiles = [
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_239557419-242000000.cut.fits',
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_242000000_244000000.cut.fits',
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_244000000_246000000.cut.fits',
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_246000000_248000000.cut.fits',
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_248000000_250000000.cut.fits',
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_250000000_252000000.cut.fits',
  '/u/gl/markusa/disk/glast_data/gammas/ft1_cut/ft1_nomSciOps_252000000-254255003.cut.fits'
]

histfile       = '/u/gl/markusa/disk/glast_data/gammas/ft2/ft2_nomSciOps_239557419-254255003.fits'

exposure_front = '/u/gl/markusa/disk/glast_data/gammas/mapcubes/fermi_diffuse_239557419-254255000_p6_v1.50bins.front.expcube.fits'
exposure_back  = '/u/gl/markusa/disk/glast_data/gammas/mapcubes/fermi_diffuse_239557419-254255000_p6_v1.50bins.back.expcube.fits'
exposure_total = '/u/gl/markusa/disk/glast_data/gammas/mapcubes/fermi_diffuse_239557419-254255000_p6_v1.100bins.expcube.fits'

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-t", "--tsfile", dest="tsfile", default='',
                  help="Name of TS output file")
parser.add_option("-o", "--outfile", dest="outfile", default='',
                  help="Name of density output file")
parser.add_option("-l", "--l", type='float', dest="l", default=0.,
                  help="L coordinate ")
parser.add_option("-b", "--b", type='float', dest="b", default=0.,
                  help="B coordinate ")
parser.add_option("-d", "--dec", type='float', dest="dec", default=-999.,
                  help="Dec coordinate ")
parser.add_option("-r", "--ra", type='float', dest="ra", default=-999.,
                  help="RA coordinate ")
parser.add_option("-f", "--fov", type='float', dest="fov", default=8.,
                  help="Field-Of-View")
parser.add_option("-s", "--pixelsize", type='float', dest="pixsize", default=0.125,
                  help="Size of pixels")
parser.add_option("-e", "--emin", type='float', dest="emin", default=200.,
                  help="Size of pixels")
parser.add_option("-B", "--bgsource", type='string', action='append', dest="bgs", 
                  help="RA coordinate ")
parser.add_option("-G", "--gasmap", dest="gasmap", default='', 
                  help="gas map to use for CO (nanten|dame|galprop)")

(options, args) = parser.parse_args()

if(options.tsfile==''):
   parser.print_help()
   sys.exit(1)

#nantenCO=sm.DiffuseFunction(COmap_dame)
#galpropCO=sm.DiffuseFunction(COmap_galp)
#img = pl.ComplexSkySpectrum("GalacticDiffuseCOvariations","([1]<=0||[0]<=0)*1.+([1]>0 &&[0]>0)*[0]/[1]",2)
#img.insert(0,nantenCO)
#img.insert(1,galpropCO)
#si=sm.SkyImage(pl.SkyDir(0,0,pl.SkyDir.GALACTIC),'dame_galprop_ratio.fits',0.125,180,1,'AIT',True)
#si.fill(img)
#del si
#sys.exit(1)
   
if(options.ra>-360.):
  sdir=pl.SkyDir(options.ra,options.dec,pl.SkyDir.EQUATORIAL)
else:  
  sdir=pl.SkyDir(options.l,options.b,pl.SkyDir.GALACTIC)

pl.SourceLikelihood.setDefaultUmax(50)
pl.SourceLikelihood.setEnergyRange(options.emin)

binner=pl.FlexibleBinner("flight/diffuse/spectrum:0/combined",0)
pl.Data.setPhotonBinner(binner)

#Set up background

galdiff = pl.HealpixDiffuseFunc(gald_total)
galCO   = pl.HealpixDiffuseFunc(gald_CO)
iso     = pl.HealpixDiffuseFunc(iso_map)


if(options.gasmap==''):
   diffuse = sm.CompositeSkySpectrum(galdiff,1.)
   diffuse.add(iso,1.);
   pl.SourceLikelihood.set_diffuse(diffuse)
else:
#   diffuse_check = sm.CompositeSkySpectrum(galdiff,1.)
#   diffuse_check.add(iso,1.);

   galpropCO=sm.DiffuseFunction(COmap_galp)
   if(options.gasmap=='nanten'): newCO=sm.DiffuseFunction(COmap_nanten)
   elif(options.gasmap=='dame'): newCO=sm.DiffuseFunction(COmap_dame)
   else:                         newCO=sm.DiffuseFunction(COmap_galp)
   diffuse = pl.ComplexSkySpectrum("GalacticDiffuseCOvariations","[0]+[1]*((([3]<=0||[2]<=0)*1.+([3]>0 &&[2]>0)*[2]/[3])-1.)+[4]",5)
   diffuse.insert(0,galdiff,1.);
   diffuse.insert(1,galCO,1.);
   diffuse.insert(2,newCO,1.);
   diffuse.insert(3,galpropCO,1.);
   diffuse.insert(4,iso,1.);

#   from numpy import *
#   for l in arange(330,334,0.1):
#     for b in arange(-2.,2.,0.1):
#        xdir=pl.SkyDir(l,b,pl.SkyDir.GALACTIC)
#	print l,b,galCO.value(xdir,200),newCO.value(xdir,200),galpropCO.value(xdir,200),diffuse.value(xdir,200.)/diffuse_check.value(xdir,200.),diffuse.value(xdir,20000.)/diffuse_check.value(xdir,20000.)
#	print l,b,diffuse.integral(xdir,200.,1000.)/diffuse_check.integral(xdir,200.,1000.),diffuse.integral(xdir,1000.,20000.)/diffuse_check.integral(xdir,1000.,20000.)
   
#   sys.exit(1)
   
pl.SourceLikelihood.set_diffuse(diffuse)

data = pl.Data(evfiles)

psl_final = pl.SourceLikelihood(data.map(),"tsmap",sdir,'point',[])
psl_final.maximize()
psl_final.printSpectrum()

drawer = pl.Draw(data.map(), diffuse, 1, options.emin, 0, 1)
if options.outfile: drawer.density(psl_final,options.outfile,options.pixsize,options.fov)
if options.tsfile: drawer.TS(psl_final,options.tsfile,options.pixsize,options.fov)

sumts_file=re.sub('\.fits$','.summed.fits',options.tsfile)
peaks_file=re.sub('\.fits$','.peaks.reg',options.tsfile)
f=pyfits.open(options.tsfile)
image=f["PRIMARY"].data

for i in range(image.shape[0]-2,-1,-1):
   image[i]+=image[i+1]
   
f["PRIMARY"].data=image

if os.path.exists(sumts_file): os.unlink(sumts_file)
f.writeto(sumts_file)   

pf=PeakFinder(25.)
pf.findPeaks(options.tsfile)
pf.writeReg(peaks_file)

