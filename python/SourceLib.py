###SEE MAIN METHOD BELOW FOR A DOCUMENTED EXAMPLE OF BASIC INTERACTIVE SPECTRAL ANALYSIS

try: #This try block only for glast-ts at University of Washington
   import sys
   sys.path.insert(0,'d:/users/kerrm/python/spectrum_dev4')
   import uw.pointlike
except:
   pass

import pointlike as pl

#Additional modules
import numpy as N
import math as M
from Fitters import *
from Response import *
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

   def update(self,**kwargs):
      """Update parameters.  CALL THIS RATHER THAN CHANGING MEMBER VARIABLES"""
      self.__dict__.update(kwargs)
      if 'event_class' in kwargs:
         self.response.update(self.event_class)
         for s in self.sources:
            s.__make_data__()
            for f in s.fitters.values(): f.remask()
      elif 'mask' in kwargs:
         for f in s.fitters.values(): f.remask()
      

   def mask(self):
      """Return an array mask corresponding to the Mask object."""
      if self.fitmask: return self.fitmask(self.response.bands(infinity=True),self.event_class)      
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

   def add_global_data(self,data_objects):
      try: list(data_objects)
      except: data_objects=[data_objects]
      self.global_data+=list(data_objects)

   def add_sourcelist(self,sourcelist,global_data_index=0):
      self.sources+=[Source(sourcelist[x].fit(),self.global_data[global_data_index]) for x in xrange(len(sourcelist))]

   def fit(self,source,model,plot=False,savepath=None,printfit=False):
      s=self[source]
      s.fit(model=model,printfit=printfit)
      if plot or savepath:
         import pylab as P
         from SourcePlot import sed
         P.clf()
         exec('sed(s,models=[s.'+s.global_data.method+'.models[-1]])')
         if savepath is not None:
            P.savefig(savepath+source+'_sed.png',dpi=75)


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class Source:
   """Provide a wrapper for PointSourceLikelihood with additional functionality."""

   def __init__(self,psl,global_data):

      self.psl,self.global_data,self.response=psl,global_data,global_data.response
      self.__make_data__()   
      self.unfolded,self.fitters=None,[]
      
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

      photons=N.fromiter((x.photons() for x in slikes),int)
      alphas=N.array([ [x.alpha(),x.sigma_alpha()] for x in slikes ])

      photons[total_mask]=0
      alphas[:,0][total_mask]=0

      #The 'astype(float)' is necessary; numpy doesn't protect against int overflow
      errors=(alphas[:,0]**2*photons+alphas[:,1]**2*photons.astype(float)**2)**0.5
      signals=N.column_stack([alphas[:,0]*photons,errors])

      if ec == 0 : my_slice = slice(0,n)
      elif ec == 1 : my_slice = slice(n,2*n+1)
      else : my_slice = slice(0,2*n+1)
           
      products=['photons','alphas','signals','slikes']
      for p in products: exec('self.%s = %s[my_slice]'%(p,p))
   
   def __getitem__(self,index):
      """Return the SimpleLikelihood object according to index."""
      return self.psl[index]

   def __call__(self):
      return self.psl

   def spectrum(self,energy_weight=1,model=lambda e: 1.):
      """General purpose routine for calculating a spectrum or counts map."""
      exposure_flag=model(100)==1. #Are we doing spectrum or counts?  Crude test.
      mask = N.logical_and(self.global_data.mask(),self.photons>0)
      e_weights = self.response.bands.centers()
      self.response.update(dir=self().dir())
      exposure = self.response(model=model) #This is either exposure OR counts
      n = len(e_weights)
      signals = N.where(mask,self.signals[:,0],0)
      exposures = N.where(mask,exposure,1e-300)
      errors = N.where(mask,self.signals[:,1],0)
      if self.global_data.event_class == -1:
         signals = signals[:n]+signals[n:]
         exposures = exposures[:n]+exposures[n:]
         errors = (errors[:n]**2+errors[n:]**2)**0.5
      if exposure_flag:
         return (e_weights,self.response.bands.diffs(from_center=True),\
                  signals*e_weights**energy_weight/exposures,\
                  errors*e_weights**energy_weight/exposures) #Domain, domain err, spectrum, spectrum err
      else:
         #exposures=N.round(exposures).astype(int)
         return (e_weights,self.response.bands.diffs(from_center=True),signals,errors,exposures)

   def f9(self):
      sig,errs=self.spectrum(energy_weight=0)[2:]
      del_e=self.response.bands.diffs()
      return ( (sig*del_e)[sig>0].sum()*1e9,((errs*del_e)**2)[sig>0].sum()**0.5*1e9 )


   def fit(self,model='PowerLaw',x0=None,method=None,printfit=True,lsfirst=True):

      method = method or self.global_data.method #Default here set to GlobalData default
      if printfit: print 'Fitting %s with method %s'%(self().name(),method)
      self.response.update(dir=self().dir())
      
      if x0: exec('model = %s(parameters=%s)'%(model,x0))
      else: exec('model = %s()'%model)

      if not method in self.fitters:
         self.fitters+=[method]
         if method=='MP':
            self.MP=PoissonFitter(source=self,method='Marginal',printfit=printfit,lsfirst=lsfirst)
         elif method=='CP':
            self.CP=PoissonFitter(source=self,method='Conditional',printfit=printfit,lsfirst=lsfirst)
         else:
            self.LS=LSFitter(source=self, printfit=printfit)
      
      exec('fitter=self.%s'%method)
      fitter.lsfirst=lsfirst
      fitter.printfit=printfit
      fitter.fit(model)
   
   def unfold(self):
      if type(self.response) is not ModelResponseDispersion:
         print '\nMust use ModelResponseDispersion to unfold!'
         return
      self.response.update(dir=self().dir())
      self.unfolded=Unfolder(self)


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
   psl1=pl.PointSourceLikelihood(bpd1,'Source_21_Fine',pl.SkyDir(249.28,-30))
   psl2=pl.PointSourceLikelihood(bpd2,'Source_21_Coarse',pl.SkyDir(249.28,-30))
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
      s1.fit(model='PowerLaw', x0 = (1e-6,3.0) ,method = m , printfit=printfit)
      s2.fit(model='PowerLaw', x0 = (1e-6,3.0) ,method = m , printfit=printfit)


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
 