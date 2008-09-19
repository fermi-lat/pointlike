import numpy as N
import math as M
import pointlike as pl
from Models import *
from Fitters import *
from scipy.stats import norm


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ExposureMap:
   """Create object wrapping an all-sky exposure map."""
   
   def __init__(self,front_emap_file=None,back_emap_file=None,**kwargs):
      """Manage front/back exposure."""
      
      self.__defaults__()
      self.__dict__.update(**kwargs)

      try:
         self.fe = pl.DiffuseFunction(front_emap_file)
         self.front_func = self.fe.value
      except: print 'Could not read file %s, using constant exposure!'%(front_emap_file)

      try:
         self.be = pl.DiffuseFunction(back_emap_file)
         self.back_func = self.be.value
      except: print 'Could not read file %s, using constant exposure!'%(back_emap_file)

   def __defaults__(self):
      self.front_const_val = self.back_const_val = 1.5e10
      self.ltfrac = 1.
      self.front_func,self.back_func = self.const_front_ex,self.const_back_ex

   def __call__(self, dir, energies, event_class=-1):
      """Return exposure for primary (all) events or for front/back independently, at given direction."""
      if event_class==-1: #Return exposure for both front and back
         return N.append(\
            N.fromiter((self.front_func(dir,e) for e in energies),float),\
            N.fromiter((self.back_func(dir,e) for e in energies),float))*self.ltfrac
      elif event_class==0: #Return front exposure only
         return N.fromiter((self.front_func(dir,e) for e in energies),float)*self.ltfrac
      else: #Return back exposure
         return N.fromiter((self.back_func(dir,e) for e in energies),float)*self.ltfrac

   def const_front_ex(self,dir,e): return self.front_const_val
   def const_back_ex(self,dir,e): return self.back_const_val

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ExposureCorrection(object):
   """Apply an exposure correction to the exposure, and manage systematic errors.  Note
      that this can handle *any* effect that comes in linearly to the calculation of the
      expected counts -- effective area, exposure, so on and so forth."""

   def __init__(self,edges,corrections,uncertainties = None):
      """Edges         -- the left edges of the bands used for the exposure corrections.
         Corrections   -- the factor by which to scale the previous exposure.
         Uncertainties -- relative uncertainty in the correction factor/exposure.  This
                          is equivalent to uncertainty in the exposure if the base is
                          totally correct, which makes sense to assume.
                          
         Note that we're treating bins for both front and back, since that seems to be the
         scheme for the Aeff validation.  Thank goodness everything is linear.  We'll just
         apply the same correction to front & back as necessary."""

      self.edges = N.asarray(edges).astype(float)
      self.corrections = N.asarray(corrections).astype(float)
      if uncertainties is not None:
         self.uncertainties = N.asarray(uncertainties).astype(float)
      else: self.uncertainties = uncertainties
     

   def __call__(self,energies,random=False,event_class = -1):
      """Return interpolated correction factors.  If random is set to True, return a
         random realization of the correction factors determined by self.uncertainties."""
      
      edges,corrections,uncertainties = self.edges,self.corrections,self.uncertainties
      if random and uncertainties is not None:
         multi = 1 + N.array([norm.rvs(0,x)[0] for x in uncertainties]) #why doesn't rvs work properly?
         #print multi
         multi[multi<0.2]=0.2
         multi[multi>1.8]=1.8

      else: multi = 1.
      altered_corrections = multi*corrections
      vals = N.empty_like(energies).astype(float)
      counter = 0
      n = len(edges)
      for i,e in enumerate(energies):
         if e < edges[0]:
            vals[i] = altered_corrections[0]
            continue
         if e > edges[-1]:
            vals[i] = altered_corrections[-1]
            continue
         for j in xrange(n-1):
            if edges[j]<= e and e <= edges[j+1]:
               e1,e2 = edges[j],edges[j+1]
               f1,f2 = altered_corrections[j],altered_corrections[j+1]
               vals[i] = f1*(e1/e)**(M.log(f1/f2)/M.log(e2/e1))
      if event_class >= 0: return vals
      else: return N.append(vals,vals)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ModelResponse(object):
   """Quickly integrate a model over exposure."""
   def __init__(self,bands,emap,**options):
      """Bands -- an EnergyBands instance.
         Emap  -- an exposure map."""   
      
      self.d=False #Boolean for dispersion use
      self.emap,self.bands = emap, bands
      self.init()
      self.__dict__.update(options)
      self.__simpson_setup__()
      self.update()

   def init(self):

      self.event_class = -1
      self.simps = 8 #Good for about 0.02% accuracy
      self.dir = pl.SkyDir(0,0)
      self.exposure_correction = None
      self.random = False
      self.ltfrac = 1.
      self.psl_c = 1.
            
   def __simpson_setup__(self):
      """Generate the model-independent portions of the integral."""

      self.s,s=[self.simps]*2
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
      #If there is an exposure correction, apply it.
      #If the self.random is True, perform random experiment with the effective area
      if self.exposure_correction is not None:
         self.exposure*=self.exposure_correction(self.sampling_points,random=self.random)
   
   def __call__(self,model=lambda e: e**-2):
      """Use Simpson's Rule to integrate the model points over the exposure; logarithmic measure.
         Dear Self: please remember in the future that this implies that passing the identity as 
         a model gives a logarithmic average of the exposure TIMES THE MEASURE, i.e., the result has
         dimension cm^2 s MeV."""

      model_points = self.sampling_points*(model(self.sampling_points))
      n = len(model_points)
            
      if self.event_class==-1:
         model_points=N.append(model_points,model_points)*self.exposure
         return self.ltfrac*self.psl_c*N.append(\
            self.int_factors*N.fromiter((N.dot(self.s_factors,model_points[:n][k:k+self.s+1])\
               for k in xrange(0,n-1,self.s)),float),\
            self.int_factors*N.fromiter((N.dot(self.s_factors,model_points[n:][k:k+self.s+1])\
               for k in xrange(0,n-1,self.s)),float))

      else:
         model_points*=exposure
         return self.ltfrac*self.event_class*self.int_factors*\
            N.fromiter((N.dot(self.s_factors,model_points[k:k+self.s+1])\
               for k in xrange(0,n-1,self.s)),float)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ModelResponseDispersion(ModelResponse):
   """Quickly integrate a model over exposure -- use dispersion."""

   def __init__(self,bands,emap,dispersion,simps=12):
      self.d=True
      self.__simpson_setup__(bands,emap,simps)
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
   main()

def main():

   #Test routine for my integration routines

   ###WITHOUT DISPERSION###

   #Easy, analytic check with constant exposure
   e=ExposureMap() #Constant value
   b=[100*10**(0.2*x) for x in xrange(11)]
   ec = ExposureCorrection(b,[1.0]*len(b))
   bands=EnergyBands(b[:-1],[b[-1]])
   m=ModelResponse(bands,e,exposure_correction = ec)
   #m.update(event_class=-1)
   b=N.array(b)
   expected=(b[1:]**-1-b[:-1]**-1)/(100**-2)*1.5e10/(-1.)*1e-9
   expected=N.append(expected,expected)
   calculated=m(PowerLaw(p=(0.01,2)))

   print 'Expected:'
   print expected
   
   print 'Calculated:'
   print calculated

   print 'Percent Differences (should be small):'
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

   print 'Percent Differences (should be small):'
   print 100*(calculated-expected)/calculated

   ###NEED TO IMPLEMENT TESTS FOR DISPERSION QUADRATURE###