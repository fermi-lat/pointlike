import numpy as N
import math as M
#import uw.pointlike
import pointlike as pl
from Models import *
from Fitters import *


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ExposureMap:
   """Create object wrapping an all-sky exposure map."""
   
   def __init__(self,emap_file=None,emap_file2=None,const_val=3e10):
      """Manage front/back exposure."""
      
      self.const_fe,self.const_be,self.fe,self.be,self.const_val=True,True,None,None,const_val
      if emap_file is not None:
         try:
            self.fe=pl.DiffuseFunction(emap_file)
            self.const_fe=False
         except: print 'Could not read file %s'%(emap_file)
      if emap_file2 is not None:
         try:
            self.be=pl.DiffuseFunction(emap_file2)
            self.const_be=False
         except: print 'Could not read file %s'%(emap_file2)         

   def __call__(self, dir, energies, event_class=-1):
      """Return exposure for primary (all) events or for front/back independently, at given direction."""
      f=self.const_ex if self.const_fe else self.fe.value
      b=self.const_ex if self.const_be else self.be.value
      if event_class==-1: #Return exposure for both front and back
         return N.append(\
            N.fromiter((f(dir,e) for e in energies),float),\
            N.fromiter((b(dir,e) for e in energies),float))
      elif event_class==0: #Return front exposure only
         return N.fromiter((f(dir,e) for e in energies),float)
      else: #Return back exposure
         return N.fromiter((b(dir,e) for e in energies),float)

   def const_ex(self,dir,e):
      """This could be expanded in the future to allow for constant exposure in space
         but provided energy dependence."""
      return self.const_val/2
         


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ModelResponse(object):
   """Quickly integrate a model over exposure."""
   def __init__(self,bands,emap,simps=16):
      """Bands -- an EnergyBands instance.
         Emap  -- an exposure map."""

      self.d=False #Boolean for dispersion use
      self.__simpson_setup__(bands,emap,simps)
      self.dir=pl.SkyDir(0.,0.)
      self.event_class=-1
      
   def __simpson_setup__(self,bands,emap,simps):
      """Generate the model-independent portions of the integral."""

      self.emap,self.bands=emap,bands
      self.s,s=[simps]*2
      self.s_factors=N.append(N.array([1.]+[4.,2.]*(s/2))[:-1],1.)
      
      bands = self.bands.all() if self.d else self.bands(infinity=True) 
      e_factors=bands[1:]/bands[0:-1] #Bin ratios - allows for uneven bins - check this!
      exp_factors=N.fromiter((x*(1./s) for x in xrange(s)),float)
      self.sampling_points=N.array( [bands[x]*e_factors[x]**exp_factors\
         for x in xrange(len(bands)-1) ] ).flatten()
      self.sampling_points=N.append(self.sampling_points,bands[-1])
      self.int_factors=N.log(e_factors)/(3*s)
      self.int_factors_2=N.hstack([self.int_factors,self.int_factors])
   
   def update(self,**options):
      """Update direction and/or event class."""

      self.__dict__.update(options)
      self.exposure=self.emap(self.dir,self.sampling_points,self.event_class)
   
   def __call__(self,model=lambda e: e**-2):
      """Use Simpson's Rule to integrate the model points over the exposure; logarithmic measure"""

      model_points = self.sampling_points*(model(self.sampling_points))
      n = len(model_points)
            
      if self.event_class==-1:

         model_points=N.append(model_points,model_points)*self.exposure
         return N.append(\
            self.int_factors*N.fromiter((N.dot(self.s_factors,model_points[:n][k:k+self.s+1])\
               for k in xrange(0,n-1,self.s)),float),\
            self.int_factors*N.fromiter((N.dot(self.s_factors,model_points[n:][k:k+self.s+1])\
               for k in xrange(0,n-1,self.s)),float))

      else:
         model_points*=exposure
         return self.int_factors*N.fromiter((N.dot(self.s_factors,model_points[k:k+self.s+1])\
            for k in xrange(0,n-1,self.s)),float)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ModelResponseDispersion(ModelResponse):
   """Quickly integrate a model over exposure -- use dispersion."""

   def __init__(self,bands,emap,dispersion,simps=16):
      self.d=True
      self.simpson_setup(bands,emap,simps)
      self.dispersion=N.fromiter((x*dispersion(x,y) for y in self.sampling_points for x in self.sampling_points),float)\
         .reshape([len(self.sampling_points),len(self.sampling_points)])
      #self.dispersion=N.array([ [x*dispersion(x,y) for y in self.sampling_points] for x in self.sampling_points])

   def __call__(self,model=lambda e:e**-2):

      model_points= self.sampling_points*model(self.sampling_points)
      
      if self.event_class==-1:
         model_points=N.hstack([model_points,model_points])*self.exposure
         s=len(self.sampling_points) #shorthand
         dispersed_points=N.empty_like(model_points)
         #Integrate over dispersion
         for x in xrange(s):
            vec=self.dispersion[x,:]*model_points[:s] #(row x) times vector of model_points
            dispersed_points[x]=(self.int_factors*N.fromiter( \
               (N.dot(self.s_factors,vec[k:k+self.s+1]) for k in xrange(0,s-1,self.s)) \
               ,float)).sum() #Do integral over each bin and sum
         for x in xrange(s):
            vec=self.dispersion[x,:]*model_points[s:] #(row x) times vector of model_points
            dispersed_points[s+x]=(self.int_factors*N.fromiter( \
               (N.dot(self.s_factors,vec[k:k+self.s+1]) for k in xrange(0,s-1,self.s)) \
               ,float)).sum() #Do integral over each bin and sum
         #Integrate over energy bins
         return N.append(
            (self.int_factors*N.fromiter((N.dot(self.s_factors,dispersed_points[k:k+self.s+1])\
            for k in xrange(0,len(dispersed_points)/2-1,self.s)),float))[self.bands.trim()],\
            (self.int_factors*N.fromiter((N.dot(self.s_factors,dispersed_points[k:k+self.s+1])\
            for k in xrange(len(dispersed_points)/2,len(dispersed_points)-1,self.s)),float))[self.bands.trim()])

      else:
         model_points*=self.exposure
         s=len(self.sampling_points) #shorthand
         dispersed_points=N.empty_like(model_points)
         #Integrate over dispersion
         for x in xrange(s):
            vec=self.dispersion[x,:]*model_points #(row x) times vector of model_points
            dispersed_points[x]=(self.int_factors*N.fromiter( \
               (N.dot(self.s_factors,vec[k:k+self.s+1]) for k in xrange(0,s-1,self.s)) \
               ,float)).sum() #Do integral over each bin and sum
         #Integrate over energy bins
         return (self.int_factors*N.fromiter((N.dot(self.s_factors,dispersed_points[k:k+self.s+1])\
            for k in xrange(0,len(dispersed_points)-1,self.s)),float))[self.bands.trim()] 

   def r_matrix(self,model=None):
      """Return a Simpson's rule approximation of the response matrix under given model."""
      if model is None:
         model_points=N.array([1]*len(self.sampling_points)) #D=1/E weighting
      else:
         model_points=N.fromiter((e*model(e) for e in self.sampling_points),float)
      
      dispersion=self.dispersion #shorthand
      s=len(self.sampling_points) #shorthand

      def integrate(sample):
      
         #Sum over input Simpson's points
         rm_temp=N.empty([len(self.int_factors),s]) 
         tmp=model_points*sample #Save some calculations
         for x in xrange(s):
            vec=dispersion[x,:]*tmp #Get the xth row; integrating over columns (input energy)
            rm_temp[:,x]=self.int_factors*N.fromiter((N.dot(self.s_factors,vec[k:k+self.s+1]) \
                                                      for k in xrange(0,s-1,self.s)) ,float)
            
         rm=N.empty([len(self.int_factors),len(self.int_factors)]) #The response matrix
         for x in xrange(rm.shape[0]):
            #Note self.int_factors cancels
            rm[x,:]=self.int_factors*N.fromiter( \
                     (N.dot(self.s_factors,rm_temp[x,k:k+self.s+1]) for k in xrange(0,s-1,self.s)) ,float)#/ \
                    #N.fromiter( \
                     #(N.dot(self.s_factors,model_points[k:k+self.s+1]) for k in xrange(0,len(model_points)-1,self.s)),float)
         sl=self.bands.trim()
         
         return rm[sl,sl] #Just bin for source
      n=len(self.sampling_points)
      if self.event_class==-1: return integrate(self.sample[:n])+integrate(self.sample[n:])
      else: return integrate(self.sample)


if __name__=='__main__':

   #Test routine for my integration routines

   ###WITHOUT DISPERSION###

   #Easy, analytic check with constant exposure
   e=ExposureMap() #Constant value
   b=[100*10**(0.2*x) for x in xrange(11)]
   bands=EnergyBands(b[:-1],[b[-1]])
   m=ModelResponse(bands,e)
   m.update(event_class=-1)
   b=N.array(b)
   expected=(b[1:]**-1-b[:-1]**-1)/(100**-2)*1.5e10/(-1.)*1e-9
   expected=N.append(expected,expected)
   calculated=m(PowerLaw(parameters=(1e-9,2)))

   print 'Expected:'
   print expected
   
   print 'Calculated:'
   print calculated

   print 'Percent Differences:'
   print 100*(calculated-expected)/calculated

   #Now, a numerical check with real exposure
   from scipy.integrate import quad
   TPpath=r'f:/glast/data/SourceDetection/' #Path to TestPattern files
   emap=ExposureMap(TPpath+'HANDOFF_front_100b_30_700000_TestPattern.fits',\
      TPpath+'HANDOFF_back_100b_30_700000_TestPattern.fits')
   m=ModelResponse(bands,emap)
   d=pl.SkyDir(200,45)
   m.update(event_class=-1,dir=d)

   expected = [quad(lambda e: emap(d,[e])[0]/e**2,b[i],b[i+1],full_output=1)[0] for i in xrange(len(b)-1)]
   expected = N.array(expected+
      [quad(lambda e: emap(d,[e])[1]/e**2,b[i],b[i+1],full_output=1)[0] for i in xrange(len(b)-1)])

   calculated = m()

   print 'Expected:'
   print expected
   
   print 'Calculated:'
   print calculated

   print 'Percent Differences:'
   print 100*(calculated-expected)/calculated

   ###NEED TO IMPLEMENT TESTS FOR DISPESRSION QUADRATURE###