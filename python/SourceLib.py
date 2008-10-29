"""Modules for the interface to the spectral fitting code.

"""
from pointlike import Data,PointSourceLikelihood,SimpleLikelihood,Draw
from skymaps import Background,DiffuseFunction,SkyDir,EffectiveArea,LivetimeCube
import pointlike as pl
from Response import *
from Fitters import *
from Models import *
from SourcePlot import *
#from fitstools import *
from psf import PSF

#Additional modules
import numpy as N
import math as M
import pylab as P
from types import *

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class GlobalData:
   """A simple wrapper to allow association of a Source with an exposure map and response functions."""
   def __init__(self,emap,response,**options):
      """Options include:

         event_class (default -1)
         mask (default None)
         method (default 'MP') <---Marginalized Poisson"""

      self.__defaults__()
      self.emap,self.response = emap,response
      self.__dict__.update(options)

   def __defaults__(self):
      self.event_class=-1
      self.fitmask=None
      self.method='MP'
      self.sources = []

   def update(self,**kwargs):
      """Update parameters.  CALL THIS RATHER THAN CHANGING MEMBER VARIABLES"""
      self.__dict__.update(kwargs)
      if 'event_class' in kwargs:
         self.response.update(self.event_class)
         for s in self.sources:
            s.__make_data__()
            for f in s.fitters.values(): f.remask()
      elif 'fitmask' in kwargs:
         for s in self.sources:
	         for f in s.fitters.values(): f.remask()
         

   def mask(self):
      """Return an array mask corresponding to the Mask object."""
      if self.fitmask is not None: return self.fitmask(self.response.bands(infinity=True),self.event_class)      
      else:
         multi = 2 if self.event_class < 0 else 1
         return N.array([True]*(len(self.response.bands(infinity=False))*multi))

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class SourceLibrary:
   """Manage a list of sources, along with global ExposureMaps, PhotonMaps, and Responses."""

   def __init__(self):
      """Initialize with a global PhotonMap, ExposureMap, and Response."""
      self.sources=[]
      self.global_data=[]
      self.index=0

   def __getitem__(self,key):
      """Return the Source object(s) corresponding to the integer/slice 'key'.
         Alternatively, search by the name associated with the psl."""
      t=type(key)
      if (t is IntType or t is SliceType):
         return self.sources[key]
      elif t is StringType:
         for x in self.sources:
            if x().name()==key:
               return x
         print key+' not found in SourceList.'
      else:
         print 'Invalid key for SourceList.  Must provide integer, slice, or valid source name.' 

   def __len__(self):
         return len(self.sources)

   def __iter__(self):
      return self   

   def next(self):
      self.index+=1
      if self.index>len(self.sources):
         self.index=0
         raise StopIteration
      else:          
         return self.sources[self.index-1]

   def ra_sort(self):
      ras = [s().dir().ra() for s in self.sources]
      self.sources = list( N.array(self.sources)[N.argsort(ras)] )

   def add_global_data(self,data_objects):
      try: list(data_objects)
      except: data_objects=[data_objects]
      self.global_data+=list(data_objects)

   def add_sourcelist(self,sourcelist,global_data_index=0):
      self.sources+=[Source(sourcelist[x].fit(),self.global_data[global_data_index]) for x in xrange(len(sourcelist))]

   def fit(self,source,model,plot=False,savepath=None,printfit=False,lsfirst=True,**kwargs):
      s=self[source]
      s.fit(model=model,printfit=printfit,lsfirst=lsfirst,**kwargs)
      if plot or savepath:
         import pylab as P
         from SourcePlot import sed
         P.clf()
         exec('sed(s,models=[s.'+s.global_data.method+'.models[-1]])')
         if savepath is not None:
            P.savefig(savepath+'/'+source+'_sed.png',dpi=75)


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class Source:
   """Provide a wrapper for PointSourceLikelihood with additional functionality."""

   def __init__(self,psl,global_data):

      self.psl,self.global_data,self.response=psl,global_data,global_data.response
      self.__make_data__()   
      self.unfolded,self.fitters=None,[]
      global_data.sources+=self
      
   def __make_data__(self):
      """Get information for both front and back from PSL, masking out an event class at the
         end of processing if necessary."""
      ec=self.global_data.event_class
      bands=N.round(self.response.bands()).astype(int)
      n=len(bands)-1
    
      event_classes=N.array([x.band().event_class() for x in self.psl])
      energies_in_psl=N.array([x.band().emin() for x in self.psl])
      if 0 in event_classes:
         iter_front=True
         front_energies=N.round(energies_in_psl[event_classes==0]).astype(int)
      else:
         front_energies = N.zeros([n])
         iter_front=False
      if 1 in event_classes:
         iter_back=True
         back_energies=N.round(energies_in_psl[event_classes==1]).astype(int)
      else:
         back_energies=N.zeros([n])
         iter_back=False
      nf,nb=len(front_energies),len(back_energies)
    
      front_mask,back_mask=N.array([False]*(n)),N.array([False]*(n))
      slikes=[0]*(2*n)
      front_index,back_index=0,0
      while iter_front and front_energies[front_index]<bands[0]: front_index+=1
      while iter_back and back_energies[back_index]<bands[0]: back_index+=1

      for i in xrange(n):
         if iter_front and front_index < nf and front_energies[front_index]==bands[i]:
            slikes[i]=self.psl[front_index+back_index]
            front_index+=1
         else: front_mask[i]=True
         if iter_back and back_index < nb and back_energies[back_index]==bands[i]:
            slikes[n+i]=self.psl[front_index+back_index]
            back_index+=1
         else: back_mask[i]=True

      for i,j in enumerate(slikes):
         if j==0: slikes[i]=self.psl[0]
      slikes=N.array(slikes)
      total_mask=N.append(front_mask,back_mask)

      photons = N.fromiter((x.photons() for x in slikes),int)
      alphas = N.array([ [x.alpha(),x.sigma_alpha()] for x in slikes ])
      umaxes = N.array( [ x.umax() for x in slikes ] )
      gammas = N.array( [ x.band().gamma() for x in slikes ] )
      ts = N.array( [x.TS() for x in slikes] )
      psl_c = 1 - (1+umaxes/gammas)**(1-gammas)

      photons[total_mask] = 0
      alphas[:,0][total_mask] = 0
      ts[total_mask] = 0

      #The 'astype(float)' is necessary; numpy doesn't protect against int overflow
      errors=(alphas[:,0]**2*photons+alphas[:,1]**2*photons.astype(float)**2)**0.5
      signals=N.column_stack([alphas[:,0]*photons,errors])

      if ec == 0 : my_slice = slice(0,n)
      elif ec == 1 : my_slice = slice(n,2*n+1)
      else : my_slice = slice(0,2*n+1)
           
      products=['photons','alphas','signals','slikes','psl_c','ts']
      for p in products: exec('self.%s = %s[my_slice]'%(p,p))
   
   def __getitem__(self,index):
      """Return the SimpleLikelihood object according to index."""
      return self.psl[index]

   def __call__(self):
      return self.psl

   def spectrum(self,energy_weight=1,model=lambda e: 1., cgs = False):
      """General purpose routine for calculating a spectrum or counts map."""
      exposure_flag=model(100)==1. #Are we doing spectrum or counts?  Crude test.
      mask = N.logical_and(self.global_data.mask(),self.photons>0) #Need to fix for upper limits
      e_weights = self.response.bands.centers()
      self.response.update(dir=self().dir(),psl_c=self.psl_c)
      exposure = self.response(model=model) #This is either exposure OR counts
      n = len(e_weights)
      #signals = N.where(mask,self.signals[:,0],0)
      exposures = N.where(mask,exposure,1e-300)
      ph = self.photons.astype(float)
      alphas,sig_alphas = self.alphas.transpose()
      signals = N.where(mask,alphas*ph,0)
      up_errs = N.where(mask,(N.minimum(1-alphas,sig_alphas)**2*ph**2 + alphas**2*ph)**0.5,0)
      down_errs = N.where(mask,(N.minimum(alphas,sig_alphas)**2*ph**2 + alphas**2*ph)**0.5,0)
   
      #errors = N.where(mask,self.signals[:,1],0)
      #errors = N.array([[up_errs[i],down_errs[i]] for i in xrange(len(up_errs))])
      if self.global_data.event_class == -1:
         signals = signals[:n]+signals[n:]
         exposures = exposures[:n]+exposures[n:]
         #errors = (errors[:n]**2+errors[n:]**2)**0.5
         up_errs = (up_errs[:n]**2+up_errs[n:]**2)**0.5
         down_errs = (down_errs[:n]**2+down_errs[n:]**2)**0.5
      if exposure_flag:
         cgs_scale = 1.60218e-6**(energy_weight-1) if cgs else 1.
         return (e_weights,self.response.bands.diffs(from_center=True),\
                  signals*cgs_scale*e_weights**energy_weight/exposures,\
                  #errors*e_weights**energy_weight/exposures
                  up_errs*cgs_scale*e_weights**energy_weight/exposures,\
                  down_errs*cgs_scale*e_weights**energy_weight/exposures)
                  #Domain, domain err, spectrum, upper error bars, lower error bars
      else:
         #exposures=N.round(exposures).astype(int)
         return (e_weights,self.response.bands.diffs(from_center=True),signals,up_errs,down_errs,exposures)

   def i_flux(self,e_weight = 0, cgs = False, scale = 1e7):
      sig,up_errs,down_errs=self.spectrum(energy_weight=e_weight,cgs=False)[2:]
      errs = N.maximum(up_errs,down_errs)
      del_e = self.response.bands.diffs(from_center=False)
      if cgs:
         #sig*=1.60218e-6**(e_weight-1)/scale
         #errs*=1.60218e-6**(e_weight-1)/scale
         sig*=1.60218e-6**(e_weight)/scale
         errs*=1.60218e-6**(e_weight)/scale
      return ( (sig*del_e)[sig>0].sum()*scale,((errs*del_e)**2)[sig>0].sum()**0.5*scale )

   def get_ts(self,accumulate=False):
      #mask = N.logical_and(self.global_data.mask(),self.photons>0)
      bands = self.response.bands.centers()
      if not accumulate:
         return bands,self.ts
      return self.ts.sum()

   def to_fits(self,filename=None):
      try:
         import pyfits as pf
      except:
         print 'Cannot find pyfits module, ergo I cannot write out a FITS file!'

      if filename is None: filename = self.name+'_spectrum.fits'
      ebins,fluxes,up_errs,down_errs = self.spectrum()
      


   def fit(self,model='PowerLaw',method=None,printfit=True,lsfirst=True,systematics=False,override=False,**kwargs):
      """Drive spectral fitting.  **kwargs are for the Model object."""

      method = method or self.global_data.method #Default here set to GlobalData default
      if printfit: print 'Fitting %s with method %s'%(self().name(),method)
      self.response.update(dir=self().dir(),psl_c=self.psl_c,random=False,override=override)

      exec('model = %s(**kwargs)'%model)

      if not method in self.fitters:
         self.fitters+=[method]
         d = {'MP':'Marginal','CP':'Conditional','EP':'Extended'} #Group Poisson fitters
         if method in d:
            exec('self.%s=PoissonFitter(source=self,method=d[method],printfit=printfit,lsfirst=lsfirst)'%(method))
         else: #Least squares fitter
            self.LS=LSFitter(source=self, printfit=printfit)
      
      exec('fitter=self.%s'%method)
      fitter.lsfirst=lsfirst
      fitter.printfit=printfit
      fitter.fit(model,systematics = systematics)
   
   def unfold(self):
      if type(self.response) is not ModelResponseDispersion:
         print '\nMust use ModelResponseDispersion to unfold!'
         return
      self.response.update(dir=self().dir())
      self.unfolded=Unfolder(self)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class SingleSourceFit(object):

   def __init__(self,datafiles,historyfile,**kwargs):

      self.init()
      self.__dict__.update(kwargs)
      self.datafiles = datafiles
      self.historyfile = historyfile

      self.psl_setup()
      #self.response_setup()
      self.open_data()
   
   def init(self):
      self.maxROI = 15.
      self.umax = 50
      self.bins = 10**N.arange(2,4.51,1./10)
      self.classlevel = 3
      self.localize = False
      self.background = None
      self.fitmask = Mask(front_ranges=[[100,1e6]],back_ranges=[[201,1e6]])
      self.tstart = 0
      self.tstop = 0
      self.pixeloutput = True
      self.ignorepix = False
      self.alignment = True
      self.theta_cut = 66.
      self.zenith_angle_cut = 105
      self.use_mc_psf = False #Use MC for PSF; if False, use estimation from data
      self.f_ecube_file = None
      self.verbose = True
      self.pixelsize = 1.
      self.irf = 'P6_v1_diff'
      self.extended = False #only for localization and signal estimation
      self.neighbours = None
      self.custom_background = None # a Background object if set
      self.use_pointlike_exposure = True

   def open_data(self):
      from pointfit import photonmap
      from skymaps import PhotonBinner
      self.pb = PhotonBinner(self.bins)
      Data.setPhotonBinner(self.pb)
      #Data.setEnergyBins(pl.DoubleVector(self.bins))
      Data.set_class_level(self.classlevel)
      Data.set_theta_cut(self.theta_cut)
      Data.set_zenith_angle_cut(self.zenith_angle_cut)
      if self.alignment:
         Data.setHistoryFile(self.historyfile)
         Data.set_alignment('default')
      self.data = photonmap(self.datafiles,pixeloutput=self.pixeloutput,\
                  tstart=self.tstart,tstop=self.tstop,ignorepix=self.ignorepix)
      self.bpd = self.data.map()

   def psl_setup(self):
      self.psf = PSF(use_mc=self.use_mc_psf)
      #self.psf.set_sigmas()
      PointSourceLikelihood.setDefaultUmax(self.umax)
      PointSourceLikelihood.set_maxROI(self.maxROI)
      PointSourceLikelihood.set_energy_range(1.01*self.bins[0])
      #if self.background is None: return
      #else:
      #   self.diffuse = DiffuseFunction(self.background)
      #   PointSourceLikelihood.set_background(self.diffuse)
   """
   def response_setup(self):
      self.bands = EnergyBands(self.bins[:-1],[self.bins[-1]])

      if self.f_ecube_file is None:
         self.emap = Exposure(self.historyfile,ft1files=self.datafiles,fovcut=self.theta_cut,zenithcut=self.zenith_angle_cut)
      else:
         fe = self.f_ecube_file
         self.emap = ExposureMap(front_emap_file = fe, back_emap_file = fe.replace('front','back'))
      self.response = ModelResponse(self.bands,self.emap)
      self.global_data = GlobalData(self.emap,self.response,fitmask = self.fitmask)
   """
   def set_source(self,name,ra,dec):

      self.src_dir = src_dir = SkyDir(ra,dec)
      self.bands = EnergyBands(self.bins[:-1],[self.bins[-1]])
      if self.f_ecube_file is None:
         if not self.use_pointlike_exposure:
            from fitstools import Exposure as my_Exposure
            self.emap = my_Exposure(self.historyfile,ft1files=self.datafiles,fovcut=self.theta_cut,zenithcut=self.zenith_angle_cut)
         else:
            self.emap = ExposureMap2(src_dir,self.historyfile,self.datafiles,self.irf,self.maxROI,self.pixelsize,self.zenith_angle_cut)
      else:
         fe = self.f_ecube_file
         self.emap = ExposureMap(front_emap_file = fe, back_emap_file = fe.replace('front','back'))
      #self.emap = ExposureMap2(src_dir,self.historyfile,self.datafiles,self.irf,self.maxROI,self.pixelsize,self.zenith_angle_cut)
      self.response = ModelResponse(self.bands,self.emap)
      self.global_data = GlobalData(self.emap,self.response,fitmask = self.fitmask)
      if self.background is not None and self.custom_background is None:
         self.diffuse = DiffuseFunction(self.background)
         self.background_obj = Background(self.diffuse, self.emap.exposure[0], self.emap.exposure[1]) # array interface does not work
         SimpleLikelihood.enable_extended_likelihood(self.extended)
         PointSourceLikelihood.set_background(self.background_obj)
      if self.custom_background is not None:
         PointSourceLikelihood.set_background(self.custom_background)
     
      
      psl = PointSourceLikelihood(self.bpd,name,src_dir)
      self.psf.set_psf_params(psl)
      """
      n_psls = []
      if self.neighbours is not None:         
         for data in self.neighbours:
            n_psl = PointSourceLikelihood(self.bpd,data[0],data[1])
            self.psf.set_psf_params(n_psl)
            n_psls += [n_psl]
            n_psl.maximize()
            psl.addBackgroundPointSource(n_psl)            
      psl.maximize()
      for n_psl in n_psls:
         n_psl.addBackgroundPointSource(psl)
         n_psl.maximize()
      psl.clearBackgroundPointSource()
      for n_psl in n_psls:
         psl.addBackgroundPointSource(n_psl)
      psl.maximize()
      self.n_psls = n_psls
      """
      if self.localize:
         self.err_radius = psl.localize()
         print 'Position difference:  %.3f +/- %.3f'%(psl.dir().difference(src_dir)*180/N.pi,self.err_radius*1.51)
      if self.verbose: psl.printSpectrum()
      self.psl = psl
      self.src = Source(psl,self.global_data)
      SimpleLikelihood.enable_extended_likelihood(False)

   def fit_spectrum(self,models=['PowerLaw']):
      for model in models:
         self.src.fit(model=model)

   def sed(self,filename=None,show=True):
      P.clf()
      sed(self.src,models=self.src.MP.models)
      if filename is not None: P.savefig(filename)
      if show: P.show()

   def i_flux(self):
      print 'Measured integrated flux (*1e7): %.2g +/- %.2g ph/cm^s/s'%self.src.i_flux()
      for model in self.src.MP.models:
         model.i_flux()

   def density_map(self,radius=10.,pix=400,smooth=True,outputname=None,equatorial=True):
      if outputname is None:
         outputname = self.src.psl.name()+'%dpix.fits'%pix
         outputname = ''.join(outputname.split())
         outputname = outputname.replace('(','_')
         outputname = outputname.replace(')','_')
      d = Draw(self.bpd)
      if equatorial: d.equatorial()
      d.region(self.psl.dir(),outputname,2.*float(radius)/pix,2*radius,smooth)
      
   def residual_image(self,scale=10,resolution=100, emin=100, emax=10000):
      try:
         self.src.MP.models[0]
      except:
         print 'No fit to calculate residuals for!'

      grid=N.linspace(-scale, scale, resolution)
      d = self.psl.dir()
      ra,dec = d.ra(),d.dec()
      cosdec = N.cos(dec*N.pi/180.)
      self.resid_image = N.array([self.__residual__(SkyDir(ra-dx/cosdec, dec-dy), emin, emax)\
                     for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
      self.image = N.array([self.__residual__(SkyDir(ra-dx/cosdec, dec-dy), emin, emax,counts=True)\
                     for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
      self.scale = scale
      #exec(';'.join( ('self.%s = %s'%x for x in ['scale','resolution'] ) )) in locals
      #self.scale, self.ra, self.dec, self.emin, self.emax, self.resolution =\
      #      scale, ra, dec, emin, emax, resolution

   def __residual__(self,sep,emin,emax,counts=False):

      delta = self.psl.dir().difference(sep)*180/N.pi #in degrees

      binc = (self.bins[:-1]*self.bins[1:])**0.5
      mask = N.array([True]*len(self.src.photons))
      n = len(mask) / 2
      for i in xrange(len(binc)):
         if binc[i] > emax or binc[i] < emin:
            mask[i] = False
            mask[i+n] = False
      binc = binc[mask[:n]]

      #psf_vals = N.array( [self.psf(e,delta,event_class=0) for e in binc] + [self.psf(e,delta,event_class=1) for e in binc] )
      #sigmas = N.array( [self.psf.sigma(e,event_class=0) for e in binc] + [self.psf.sigma(e,event_class=1) for e in binc] )
      #sigmas = sigmas*N.pi/180.
      #signals = (self.src.response(self.src.MP.models[0])*self.src.psl_c)[mask]
      #backs = self.src.photons[mask] - signals
      #densities = N.array([x.display(sep,0) for x in self.src.slikes[mask]])
      #density = self.bpd.density(sep)
      #But this is just source density!  Maybe need to use BPD.
      #back_densities = N.array([x.display(sep,2) for x in self.src.slikes[mask]])
      if not counts:
         #pre = N.array([x.display(sep,3) for x in self.src.slikes[mask]])
         resids = N.array([x.display(sep,4) for x in self.src.slikes[mask]])
         #p = pre.sum()
         #p = p if p > 0 else 1.
         return resids.sum()#/p
      else:
         return N.array([x.display(sep,1) for x in self.src.slikes[mask]]).sum()
      #back_densities = 0
      #if N.all(back_densities == 0.): back_densities = 1/self.umax #not right

      #print psf_vals
      #print signals
      #print backs
      #print densities
      #print back_densities
      #print signals*psf_vals/(N.pi*2*sigmas**2)
      #print densities

      #return (signals*psf_vals/(N.pi*2*sigmas**2) + backs*back_densities).sum() - density
      #return (signals*psf_vals/(N.pi*2*sigmas**2) + backs*back_densities - densities).sum()
      
      

   def show(self, **kwargs):
      import pylab
      scale=self.scale
      P.subplot(121)
      pylab.imshow(self.resid_image, extent=[-scale, scale, -scale, scale],interpolation='nearest', **kwargs)
      pylab.axvline(0, color='white')
      pylab.axhline(0, color='white')
      pylab.colorbar()
      P.subplot(122)
      pylab.imshow(self.image, extent=[-scale, scale, -scale, scale],interpolation='nearest', **kwargs)
      pylab.axvline(0, color='white')
      pylab.axhline(0, color='white')
      pylab.colorbar()
      #emin, emax = energy_range(self.level)
      #pylab.title('level %d (%d-%d)' %(self.level, emin, emax), size=10)


   def __call__(self):
      return self.src
      
""" A bit outdated.            


if __name__=='__main__':

   ###------------------DEMONSTRATION AND BASIC GUIDE TO SPECTRAL FITTING------------------###
   ###The following section demonstrates fitting an individual source using two binning schemes.

   ###There is method to the commend madness -- conversation comments have 3 pound signs, while
   ###optional bits of code are commented out with just one.
      
   ###Create an ExposureMap object with the output of gtexpcube (files for front, back events)
   ###Note that the minimum, maximum energies must be outside of the extremes used for fitting

   TPpath=r'f:/glast/data/SourceDetection/' #Path to TestPattern files
   emap=ExposureMap(TPpath+'HANDOFF_front_100b_30_700000_TestPattern.fits',\
      TPpath+'HANDOFF_back_100b_30_700000_TestPattern.fits')     


   ###Make two different sets of energy bands.
   ###Note that the badnds must be consistent with those used for PointSourceLikelihood.
   ###This is most easily done by specifying the bands for both simultaneously, as done here.

   b1=(10**N.arange(2,5.61,0.2)).tolist() #5 per decade
   b2=(100*2.35**N.arange(0,10)).tolist() #2.35 multiplier
   bands=[EnergyBands(b[:-1],[b[-1]]) for b in [b1,b2]]
   
   
   ###Make a ModelResponse object (calculates counts under spectral models -- needs exposure)
   ###This version does not do energy dispersion

   response1,response2=[ModelResponse(b,emap) for b in bands]
   
   ###The below version does do energy dispersion (slower, marginal difference - uncomment to use)
   ###Note that the default behaviour is to account for scattering into and out of the ROI by adding
   ###a bin to either side of the user-specified bins.  This can be overridden by setting the bins
   ###manually -- see the EnergyBins class.  Default Dispersion -- see Models.py

   #response1,response2=[ModelResponseDispersion(b,emap,Dispersion()) for b in bands]

   
   ###One can optionally set up a Mask object to specify allowed energy ranges
   
   #mask1=mask2=Mask(front_ranges=[[30,300000]],back_ranges=[[30,300000]])
   mask1=mask2=None
   

   ###The GlobalData object wraps up everything about the exposure and fitting environment
   
   global_data1=GlobalData(emap,response1,event_class=-1,fitmask=mask1)
   global_data2=GlobalData(emap,response2,event_class=-1,fitmask=mask2)
   
   
   ###Now, need photon events -- TWO WAYs TO ACCESS DATA
   
   ###1) STARTING FROM FT1 FILES
   
   #from glob import glob
   #files=glob('f:/glast/data/SourceDetection/pl*v2.fits') #events from TestPattern
   #files_bg=glob('f:/glast/data/SourceDetection/bg_low*v2.fits') #background from TestPattern
   #pl.Data.setEnergyBins(bands1(infinity=True)) #Set the energy bands in Pointlike
   #data1=pl.Data(files,0,0,0,20000) #Get all events, select on Source ID
   #data1_bg=pl.Data(files_bg,-1,0,0,-1)
   #pl.Data.setEnergyBins(bands2(infinity=True)) #Set different energy bands for the next set of data
   #data2=pl.Data(files,0,0,0,20000)
   #data2_bg=pl.Data(files_bg,-1,0,0,-1)
   #bpd1=data1.map()
   #bpd1.add(data1_bg.map())
   #bpd2=data2.map()
   #bpd2.add(data2_bg.map())
   
   ###2) STARTING FROM BINNED DATA - bins must be consistent with those specified above!!!

   bpd1=pl.BinnedPhotonData('fine.fits') #Load the finely-binned data
   #bpd1=pl.BinnedPhotonData('fine_bg.fits') #with background
   bpd2=pl.BinnedPhotonData('coarse.fits')
   #bpd2=pl.BinnedPhotonData('coarse_bg.fits')
   
   
   ###Configure PointSourceLikelihood   
   pl.PointSourceLikelihood.set_energy_range(b1[0]*1.01) #Set minimum energy in PSL

   ###Optional: set up the background -- isotropic diffuse for the TestPattern
   
   #diffuse_model=pl.IsotropicPowerLaw(2e-5,2.1)
   #diffuse=pl.Background(diffuse_model,emap.fe,emap.be)
   #pl.PointSourceLikelihood.set_diffuse(diffuse)
   
   ###Instantiate PointSourceLikelihood objects for each source   
   psl1=pl.PointSourceLikelihood(bpd1,'Source_21_Fine',sm.SkyDir(249.28,-30))
   psl2=pl.PointSourceLikelihood(bpd2,'Source_21_Coarse',sm.SkyDir(249.28,-30))
   for p in [psl1,psl2]:
      p.maximize()
      p.localize()
      p.maximize()

   ###Combine pointlike and speclike by creating Source objects
   s1=Source(psl1,global_data1)
   s2=Source(psl2,global_data2)

   ###Now, we can do some fits -- we test out and compare the three different methods of spectral fitting
   ###1) Marginal Poisson is the "best" and slowest; it uses Poisson statistics, good for low counts,
   ###     and also includes the uncertainty in the background/source separation by Pointlike.
   ###2) Conditional Poisson -- quick, uses Poisson statistics, but does *not* factor in the uncerainty
   ###     from Pointlike.
   ###3) Least Squares -- fastest, uses Gaussian approximation so only accurate for ~10 or more photons in
   ###     a bin.  Uses "propagation of errors" estimation to combine Poisson and Pointlike uncertainty.

   ###To change event class or energy masks, call "update" in the appropriate GlobalData class.
   ###We can also specify the starting position in parameter space.
   ###Note these initial positions are far from the nominal values (2e-9, 1.6) to test convergence.
   
   
   ###Test out the fit with all three methods

   printfit=True #True for verbose fitting
   methods = ['MP', 'CP', 'LS']
   for m in methods:
      s1.fit(model='PowerLaw', p = (1e-6,3.0) ,method = m , printfit=printfit)
      s2.fit(model='PowerLaw', p = (1e-6,3.0) ,method = m , printfit=printfit)


   ###This is a demonstration of accessing the fitted parameters.   
   
   method_names=['Marginal Poisson', 'Conditional Poisson', 'Least Squares']
   for i in xrange(len(methods)):
      exec('fit1=s1.%s.models[0]'%methods[i])
      exec('fit2=s2.%s.models[0]'%methods[i])
   
      print '\n\nResults for the %s fitter:'%method_names[i]
      print 'Fine binning  : Norm: %.2e +/- %.2e, Index: %.2f +/- %.2e'\
         %(fit1.p[0],fit1.cov_matrix[0,0]**0.5,fit1.p[1],fit1.cov_matrix[1,1]**0.5)
      print 'Coarse binning: Norm: %.2e +/- %.2e, Index: %.2f +/- %.2e'\
         %(fit2.p[0],fit2.cov_matrix[0,0]**0.5,fit2.p[1],fit2.cov_matrix[1,1]**0.5)



   ###-------------------------PRIMITIVE UNFOLDING TO ESTIMATE SED-----------------------###
   ###An inversion of the IRF matrix allows one to estimate the spectral energy distribution.
   ###One must use ModelResponseDispersion rather than ModelResponse.  If one unfolds a source,
   ###the unfold spectrum can be viewed by calling sed -- see section below on Plotting.
   
   s1.unfold()
   s2.unfold()

   
   
   ###---------------------PLOTTING -- REQUIRED Pylab (Matplotlib)--------------------###
   
   from SourcePlot import *

   ###Plot the spectral energy density of the Marginal Poisson fits
   
   sed(s1,models=s1.MP.models,fignum=30)
   sed(s2,models=s2.MP.models,fignum=31)

   ###We can also plot signal counts vs. model counts
   
   counts(s1,models=s1.MP.models,fignum=32)
   counts(s2,models=s2.MP.models,fignum=33)
 
 """