from sys import path

import pointlike as pl
import skymaps as sm
from numpy import *
import sys
import re
import os

from optparse import OptionParser

import datafiles

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-s", "--source", dest="source", default='SRC0000',
                  help="Name of the source")
parser.add_option("-o", "--outfile", dest="outfile", default='',
                  help="Name of density output file")
parser.add_option("-i", "--imagedir", dest="image_dir", default='images',
                  help="Name of density output file")
parser.add_option("-l", "--l", type='float', dest="l", default=-999.,
                  help="L coordinate ")
parser.add_option("-b", "--b", type='float', dest="b", default=-999.,
                  help="B coordinate ")
parser.add_option("-d", "--dec", type='float', dest="dec", default=-999.,
                  help="Dec coordinate ")
parser.add_option("-r", "--ra", type='float', dest="ra", default=-999.,
                  help="RA coordinate ")
parser.add_option("-B", "--bgsource", type='string', action='append', dest="bgs", default=[], 
                  help="Background source. Format <name>:<ra>,<dec> or <name>:<l>,<b>g ")
parser.add_option("-e", "--emin", type='float', dest='emin', default=500.,
                  help="Minimum energy in fit")
parser.add_option("-c", "--combine-fb",  action='store_true', dest='combined', default=False,
                  help="Combine front and back bins")
parser.add_option("-f", "--spectral-fit", type='string', dest='fit', default='plaw',
                  help="Select spectral fits to perform. Any list of plaw,bplaw,expcut, Use plaw:emin:emax to select range")
parser.add_option("-S", "--smoothed-map", action='store_true', dest='smoothmap', default=False,
                  help="Create time consuming smoothed counts map")
parser.add_option("-T", "--ts-map", action='store_true', dest='tsmap', default=False,
                  help="Create time consuming ts maps")
parser.add_option("-m", "--more-bins", action='store_true', dest='morebins', default=False,
                  help="Use more energy bins for strong sources")
parser.add_option("-p", "--point-source-only", action='store_true', dest='psonly', default=False,
                  help="Do only point source fit")
parser.add_option("-x", "--fix-position", action='store_true', dest='fixpos', default=False,
                  help="Fix the position to the initial guess")
parser.add_option("-M", "--max-size", type='float', dest='maxsize', default=180.,
                  help="Maximum source size allowed in fit in arcmin")
parser.add_option("-F", "--ts-file", dest="tsfile", default='',
                  help="Initialize fit from peak finding in TS-map ")
parser.add_option("-R", "--results-file", dest="resfile", default='',
                  help="Initialize fit from results file ")
parser.add_option("-D", "--debug", action='store_true', dest='debug', default=False,
                  help="Use only first data file")
parser.add_option("-I", "--size-scan", dest="sizescan", default='',
                  help="Scan TS values for various source sizes. Parameter: source radius range smin:smax:nsteps [arcmin]")
parser.add_option("-G", "--gasmap", dest="gasmap", default='', 
                  help="gas map to use for CO (nanten|dame|galprop)")

(options, args) = parser.parse_args()

if(options.source==''):
   parser.print_help()
   sys.exit(1)

if(options.outfile==''): options.outfile=options.source+'.fits'
if(os.path.isdir(options.outfile)): 
    options.image_dir=options.outfile+'/'+options.image_dir
    options.outfile=options.outfile+'/'+options.source+'.fits'
    
if not os.path.exists(options.image_dir): os.mkdir(options.image_dir)

#Set up environment


binner_type="flight/diffuse"
if(options.morebins): binner_type+="/spectrum:+"
else: binner_type+="/spectrum:0"

deg2rad=pi/180.

fit_models={ 'plaw': pl.PowerLaw(), 'bplaw': pl.BrokenPowerLaw(), 'expcut': pl.ExpCutoff() }
extended_source_types=['gauss','disk','radius_gauss']

bgsources={}
sources={}

scan_sourcesize=False
if(options.sizescan!=''):
  sizemin,sizemax,nsteps=options.sizescan.split(':')
  sizemin,sizemax,nsteps=float(sizemin)/60.*deg2rad,float(sizemax)/60.*deg2rad,int(nsteps)
  scanfile=re.sub('.fits','.sizescan.txt',options.outfile)
  scan_sourcesize=True

for bgs in options.bgs:
  name,coor=bgs.split(':')
  c1,c2=re.sub('g$','',coor).split(',')
  print "Adding background source: ",name,c1,c2
  if re.search('g$',coor):
     sdir=pl.SkyDir(float(c1),float(c2),pl.SkyDir.GALACTIC)
  else:   
     sdir=pl.SkyDir(float(c1),float(c2),pl.SkyDir.EQUATORIAL)
  bgsources[name] = (sdir,'point',[] )

if(options.debug): datafiles.evfiles=[datafiles.evfiles[0]]

if(options.tsfile!=''):
   from PeakFinder import PeakFinder
   pf=PeakFinder(25.)
   pf.findPeaks(options.tsfile)
   pkmap={}
   for (l,b),ts in zip(pf.peaks,pf.peakts): pkmap[ts]=(l,b)
      
   print "Found %d sources in map %s:"%(len(pf.peaks),options.tsfile)
   keys=list(sort(pkmap.keys()))
   keys.reverse()
   for i,ts in enumerate(keys):
        name="%s_PK%.2d"%(options.source,i)
	l,b =pkmap[ts]
	print "   Adding background source %s at l=%.1f,b=%.1f with TS=%.1f. "%(name,l,b,ts)
        sdir=pl.SkyDir(l,b,pl.SkyDir.GALACTIC)
        bgsources[name] = (sdir,'point',[] )
        sources[name] = (sdir,'point',[] )
        if(not options.psonly):
	    sources[name+':gauss'] = (sdir,'gauss',[0.1*deg2rad] )
            sources[name+':disk'] = (sdir,'disk',[0.1*deg2rad] )

elif(options.resfile!=''):
   from PeakFinder import ResultsFileParser
   pf=ResultsFileParser(options.resfile)
   srcs=pf.sourceList()
   pkmap={}
   srcNames=[]
   for sname,l,b,stype,ts,r in srcs: 
      pkmap[ts]=(l,b,sname,stype,r)
      srcNames.append(sname)
      
   print "Found %d sources in resultsfile %s:"%(len(srcs),options.resfile)
   keys=list(sort(pkmap.keys()))
   keys.reverse()
   for i,ts in enumerate(keys):
	l,b,sname,stype,r =pkmap[ts]
	print "   Adding background source %s (%s) at l=%.1f,b=%.1f with TS=%.1f. "%(sname,stype,l,b,ts)
        sdir=pl.SkyDir(float(l),float(b),pl.SkyDir.GALACTIC)
        if(stype in ['point','pseudopoint']): r=[]
        if(stype in extended_source_types): r=[float(r)]	
        if(not options.psonly):
	   bgsources[sname+':'+stype] = (sdir,stype,r )
           if (options.source not in srcNames) or (sname==options.source):
	      if sname==options.source and (options.l>-1 or options.ra>-1):
	         if(options.l>-1): sdir= pl.SkyDir(options.l,options.b,pl.SkyDir.GALACTIC)
		 if(options.ra>-1): sdir= pl.SkyDir(options.ra,options.dec,pl.SkyDir.EQUATORIAL)
                 if(stype in extended_source_types): r=[0.1*deg2rad]	
	      sources[sname+':'+stype] = (sdir,stype,r )
	else:   
           bgsources[sname] = (sdir,'point',[] )
           sources[sname] = (sdir,'point',[] )

else:
    if(options.dec<-900. and options.ra<-900.):
      sources[options.source] = ( pl.SkyDir(options.l,options.b,pl.SkyDir.GALACTIC),'point',[] )
      if(not options.psonly):
	  sources[options.source+':gauss'] = ( pl.SkyDir(options.l,options.b,pl.SkyDir.GALACTIC),'gauss',[0.1*deg2rad] )
	  sources[options.source+':disk'] = ( pl.SkyDir(options.l,options.b,pl.SkyDir.GALACTIC),'disk',[0.1*deg2rad] )
    else:
      sources[options.source] = ( pl.SkyDir(options.ra,options.dec,pl.SkyDir.EQUATORIAL),'point',[] )
      if(not options.psonly):
	  sources[options.source+':gauss'] = ( pl.SkyDir(options.ra,options.dec,pl.SkyDir.EQUATORIAL),'gauss',[0.1*deg2rad] )
	  sources[options.source+':disk'] = ( pl.SkyDir(options.ra,options.dec,pl.SkyDir.EQUATORIAL),'disk',[0.1*deg2rad] )

if(options.combined): binner_type+="/combined"

binner=pl.FlexibleBinner(binner_type,0)
pl.Data.setPhotonBinner(binner)

counter = 0

pl.SourceLikelihood.setMinuitMode(['simplex','hesse','grad'])
pl.SourceLikelihood.setDefaultRoI(3.5*deg2rad)
pl.SourceLikelihood.setDefaultUmax(50.)
pl.SourceLikelihood.setEnergyRange(70.)
pl.SourceLikelihood.setMaxSize(options.maxsize/60.*deg2rad)

pl.SourceLikelihood.setVerbose(1)

#pl.Data.setEnergyBins(bins)

#isotropic     = sm.IsotropicPowerLaw(2e-5,2.1)
isotropic     = pl.HealpixDiffuseFunc(datafiles.iso_map)

if(options.combined):
   frontExposure       = sm.DiffuseFunction(datafiles.exposure_total);
   backExposure        = frontExposure
   
   galDiffuse    = pl.HealpixDiffuseFunc(datafiles.gald_total);
   if (options.gasmap!=''):
       galpropCO_in  = sm.DiffuseFunction(datafiles.COmap_galp)
       galpropCO_out = pl.HealpixDiffuseFunc(datafiles.gald_CO)
       if(options.gasmap=='nanten'): newCO=sm.DiffuseFunction(datafiles.COmap_nanten)
       elif(options.gasmap=='dame'): newCO=sm.DiffuseFunction(datafiles.COmap_dame)
       else:                         newCO=sm.DiffuseFunction(datafiles.COmap_galp)
       diffuseModel_front = pl.ComplexSkySpectrum("GalacticDiffuseCOvariations",\
            "[0]+[1]*((([3]<=10.||[2]<=10.||isnan([2])||isnan([3]))*1.+([3]>10. &&[2]>10.)*[2]/[3])-1.)+[4]",5)
       diffuseModel_front.insert(0,galDiffuse,1.);
       diffuseModel_front.insert(1,galpropCO_out,1.);
       diffuseModel_front.insert(2,newCO,1.);
       diffuseModel_front.insert(3,galpropCO_in,1.);
       diffuseModel_front.insert(4,isotropic,1.);
   else:
       diffuseModel_front  = sm.CompositeSkySpectrum(galDiffuse,1.)
       diffuseModel_front.add(isotropic,1.)
   diffuseModel_back   = diffuseModel_front
#   si=sm.SkyImage(pl.SkyDir(0,0,pl.SkyDir.GALACTIC),'test.fits',0.125,180,1,'AIT',True)
#   si.fill(diffuseModel_front)
#   del si

else :
   if (options.gasmap!=''): raise RuntimeError,"Can use gasmap feature only with combined front/back bins."
   frontExposure       = sm.DiffuseFunction(datafiles.exposure_front);
   backExposure        = sm.DiffuseFunction(datafiles.exposure_back);
   
   galDiffuse_front    = pl.HealpixDiffuseFunc(datafiles.gald_front);
   galDiffuse_back     = pl.HealpixDiffuseFunc(datafiles.gald_back);
   diffuseModel_front  = sm.CompositeSkySpectrum(galDiffuse_front,1.)
   diffuseModel_back   = sm.CompositeSkySpectrum(galDiffuse_back,1.)
   diffuseModel_front.add(isotropic,1.)
   diffuseModel_back.add(isotropic,1.)

pl.SourceLikelihood.set_diffuse(diffuseModel_front,0)
pl.SourceLikelihood.set_diffuse(diffuseModel_back,1)
pl.SourceLikelihood.set_exposure(frontExposure,0)
pl.SourceLikelihood.set_exposure(backExposure,1)
   
spectral_fits={}
fit_rfiles={}
fitters=[]

if options.fit:
   fitnames=options.fit.split(',')
   for fitn in fitnames:
       fitp=fitn.split(':')
       if (fitp[0] not in fit_models.keys()):  raise RuntimeError("%s not known as spectral model"%(plaw)) 
       if(len(fitp)==1): spectral_fits[fitp[0]]=(200.,2.e5)
       elif(len(fitp)==3): spectral_fits[fitp[0]]=(float(fitp[1]),float(fitp[2]))
       else: raise RuntimeError("%s is not a valid description of a spectral fit. Format name:emin:emax"%(fitn))
   
data = pl.Data(datafiles.evfiles)

if(len(spectral_fits.keys())==0): 
    rfile= pl.ResultsFile(options.outfile,data,len(sources))
else:
    for key in spectral_fits.keys():   
         fit_rfiles[key]   = pl.ResultsFile(re.sub('.fits','.'+key+'.fits',options.outfile),data,len(sources))

for name in sort(sources.keys()):
   psl = pl.SourceLikelihood(data.map(),name,sources[name][0],sources[name][1],sources[name][2])
   psl.maximize()
   
   bgpsl=[]
   bgts=[]
   bgts_old=[]
   
   for bname in sort(bgsources.keys()):
      sbasename=re.sub(':\w+$','',name)
      bbasename=re.sub(':\w+$','',bname)
      print bbasename,'<-->',sbasename
      if(bbasename==sbasename): continue
      bgpsl.append(pl.SourceLikelihood(data.map(),bname,bgsources[bname][0],bgsources[bname][1],bgsources[bname][2]))
      bgpsl[-1].maximize()
      bgts.append(bgpsl[-1].TS())
      print "background source:",bname," TS=",bgts[-1]
    
   bgpsl.append(psl)   
   bgts.append(psl.TS())
   bgts=array(bgts)
   bgts_old=zeros_like(bgts)

   for b1 in range(len(bgpsl)):
      bgpsl[b1].clearBackgroundPointSource()
      for b2 in range(len(bgpsl)):
         if b1!=b2: bgpsl[b1].addBackgroundPointSource(bgpsl[b2],False)
      
   print "Starting background iteration with TS=",bgts
   iteration_counter=0
   while(max(abs((bgts-bgts_old)/bgts))>0.01 and iteration_counter<20):
      bgts_old=array(bgts)
      for b1 in range(len(bgpsl)):
         bgpsl[b1].setDir(bgpsl[b1].dir())
	 bgpsl[b1].maximize()
	 bgts[b1]=bgpsl[b1].TS()
      print "Iterating background sources: TS=",bgts," deltaTSmax=",max(abs(bgts-bgts_old)/bgts)
      iteration_counter+=1
#       for bpsl in bgpsl[:-1]: psl.addBackgroundPointSource(bpsl)

   psl.setDir(psl.dir())
   psl.printSpectrum()
   pl.SourceLikelihood.setEnergyRange(options.emin)
   if (options.fixpos): psl.fixPosition()
   psl.localize()
   psl.printSpectrum()
   pl.SourceLikelihood.setEnergyRange(70.)
   psl.maximize()
   psl.printSpectrum()

   if scan_sourcesize and (sources[name][1] not in ['point','pseudopoint']):
       ts_size={}
       pl.SourceLikelihood.setEnergyRange(options.emin)
       psl_scan = pl.SourceLikelihood(data.map(),name,psl.dir(),sources[name][1],psl.sourceParameters())
       for bpsl in bgpsl[:-1]: psl_scan.addBackgroundPointSource(bpsl)
       for r in arange(sizemin,sizemax+1e-10,(sizemax-sizemin)/float(nsteps-1)):
	   psl_scan.setDir(psl.dir(),[r])
	   psl_scan.maximize()
	   ts_size[r]=psl_scan.TS()
	   print "Size scan: r=",(r/deg2rad*60.),"arcmin TS=",ts_size[r]
       scan_file=open(scanfile,'w')    
       for r in sort(ts_size.keys()):
	  scan_file.write("%.4f %.4f\n"%(r/deg2rad*60.,ts_size[r]))
       scan_file.close()   
       pl.SourceLikelihood.setEnergyRange(70.)
       
       
   
#   psl_test = pl.SourceLikelihood(data.map(),name,pl.SkyDir(psl.dir().ra()+0.5,psl.dir().dec()-0.5),sources[name][1],psl.sourceParameters())
#   for bpsl in bgpsl: psl_test.addBackgroundPointSource(bpsl)
#
#   pl.SourceLikelihood.setEnergyRange(options.emin)
#   psl_test.localize()
#   pl.SourceLikelihood.setEnergyRange(70.)
#   psl_test.maximize()
#   psl_test.printSpectrum()
#   deltapsi=psl.dir().difference(psl_test.dir())*180./pi
#
#   if psl_test.TS()>25 and deltapsi>0.1: 
#       print "Warning found another source in field of view. TS=",psl_test.TS()," delta psi=",deltapsi
#       bgpsl.append(psl)
#       bgpsl[-2]=psl_test
#   else: 
#       print "No other apparent source in field of view. TS=",psl_test.TS()," delta psi=",deltapsi
#       break     

#   psl_final = pl.SourceLikelihood(data.map(),name,psl.dir(),sources[name][1],psl.sourceParameters())
#   for bpsl in bgpsl[:-1]: psl_final.addBackgroundPointSource(bpsl)
#   psl_final.maximize()
#   psl_final.printSpectrum()
  
   psl_final=psl
    
   if(len(spectral_fits.keys())==0): 
       rfile.fill(psl_final)
   else:
       for key in spectral_fits.keys():
           fitter0= pl.SpectralFitter(psl_final,fit_models[key])
           fitter0.useUpperLimit(0)
	   fitter0.setFitRange(spectral_fits[key][0],spectral_fits[key][1])
	   fitter0.specfitMinuit()
           fit_rfiles[key].fill(psl_final,fitter0)
	   fitters.append(fitter0)


   pl.SourceLikelihood.setEnergyRange(options.emin)

   print "Filling images for %s."%(name)
   image_file=options.image_dir+'/image_'+re.sub(':','_',name)
   images=["counts","model","background","count_residual","ts_residual"]
   if(options.smoothmap and sources[name][1]=='point'):images.append("smoothed_counts")
   psl_final.createImages(images,image_file,"ZEA","GAL",0.125,8.,2)
   
   if options.tsmap:
       tsdrawer = pl.Draw(data.map(), diffuseModel_front, 1, options.emin, 0, 1)

       psl_src  = pl.SourceLikelihood(data.map(),name,psl.dir(),'point',[])
       for bpsl in bgpsl[:-1]: psl_src.addBackgroundPointSource(bpsl)
       tsfile = options.image_dir+'/tsmap_src_'+re.sub(':','_',name)+'.fits'
       tsdrawer.TS(psl_src,tsfile,0.125,8.)

       psl_res  = pl.SourceLikelihood(data.map(),name,psl.dir(),'point',[])
       for bpsl in bgpsl: psl_res.addBackgroundPointSource(bpsl)
       psl_res.maximize()
       tsfile = options.image_dir+'/tsmap_res_'+re.sub(':','_',name)+'.fits'
       tsdrawer.TS(psl_res,tsfile,0.125,8.)

#       psl_nobg = pl.SourceLikelihood(data.map(),name,psl.dir(),'point',[])
#       psl_nobg.maximize()
#       tsfile = options.image_dir+'/tsmap_all_'+re.sub(':','_',name)+'.fits'
#       tsdrawer.TS(psl_nobg,tsfile,0.125,8.)

   
if(len(spectral_fits.keys())==0): 
    del rfile
else:
    for key in spectral_fits.keys(): 
#      print "writing",key
      del fit_rfiles[key]
#      print ""
            
print "Finished."
del fitters
print "deleted fitters."
del fit_models
print "deleted fit models"
del backExposure
del frontExposure
print "deleted diffuse objects"
del  diffuseModel_front
del  diffuseModel_back
print "deleted more diffuse objects"
del isotropic
print "deleted even more diffuse objects"

