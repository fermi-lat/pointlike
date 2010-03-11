"""
A module to manage the PSF from CALDB and handle the integration over
incidence angle and intepolation in energy required for the binned
spectral analysis.
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pypsf.py,v 1.3 2010/01/21 01:29:42 kerrm Exp $
author: M. Kerr

"""

import pyfits as pf
import numpy as N
from os.path import join
from cPickle import load
from skymaps import ExposureWeighter,SkyDir
from scipy.integrate import quad,simps
from math import cos,sin

def scale_factor(e,c0,c1,exp):
   return ( (c0*(e/100)**(exp))**2 + c1**2 )**0.5

TWOPI = (2*N.pi)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class Psf(object):

   def init(self):
      self.irf     = 'P6_v3_diff'
      self.psf_irf = None # a possible override to use an on-orbit PSF file

   def __init__(self,CALDB,**kwargs):
      self.CALDB = join(CALDB,'bcf')
      self.init()
      self.__dict__.update(kwargs)

      self.__readCALDB__()
      self.__read_params__()
      self.__calc_weights__()

   def __readCALDB__(self):

      # open CALDB files
      if self.psf_irf is not None:
         print 'Overriding default PSF; using %s'%(self.psf_irf)
      psf_files = [join(self.CALDB,'psf','psf_%s_%s.fits'%(self.psf_irf or self.irf,x)) for x in ['front','back']]
      try:
         h0,h1 = self.CALDBhandles = [pf.open(x) for x in psf_files]
      except:
         print 'Could not open CALDB files for PSF!  Aborting!'
         raise Exception

      # read in stuff that doesn't depend on conversion type
      self.scale_factors = N.asarray(h0[2].data.field('PSFSCALE')).flatten()
      self.e_los         = N.asarray(h0[1].data.field('ENERG_LO')).flatten().astype(float)
      self.e_his         = N.asarray(h0[1].data.field('ENERG_HI')).flatten().astype(float)
      self.c_los         = N.asarray(h0[1].data.field('CTHETA_LO')).flatten().astype(float)
      self.c_his         = N.asarray(h0[1].data.field('CTHETA_HI')).flatten().astype(float)

      sf = self.scale_factors
      self.scale_func = [lambda e: ( (sf[0]*(e/100)**(sf[-1]))**2 + sf[1]**2 )**0.5,
                         lambda e: ( (sf[2]*(e/100)**(sf[-1]))**2 + sf[3]**2 )**0.5 ]


   def __calc_weights__(self,livetimefile='',skydir=None):

      aeffstrs = [join(self.CALDB,'ea','aeff_%s_%s.fits'%(self.irf,ct)) for ct in ['front','back']]
      ew       = ExposureWeighter(aeffstrs[0],aeffstrs[1],livetimefile)
      dummy    = skydir or SkyDir() 

      elo,ehi,clo,chi = self.e_los,self.e_his,self.c_los[::-1],self.c_his[::-1]

      weights  = N.zeros([2,len(elo),len(clo)])
      for i in xrange(len(elo)):                     # iterator over energies
         em = (elo[i]*ehi[i])**0.5
         for j in [0,1]:                             # iterate over conversion types
            for k,(c0,c1) in enumerate(zip(clo,chi)):# iterator over cos(theta) on-axis to edge
               if c0 < 0.35: continue                # exclude bins below cos(theta)=0.4
               weights[j,i,k] = ew(c0,c1,elo[i],ehi[i],j,dummy) # exposure at energy/cos(theta)
            weights[j,i,:] /= weights[j,i,:].sum()   # normalize weights over cos(theta)

      self.weights = weights

   def set_weights(self,livetimefile,skydir):
      self.__calc_weights__(livetimefile,skydir)

   def get_p(self,e,ct,cthetabin=None):
      ind   = min(N.searchsorted(self.e_his,e),len(self.e_his) - 1)
      p     = self.tables[ct,:,ind,:]
      w     = self.weights[ct,ind,:]
      if cthetabin is None:
         return N.append(p,[w],axis=0)
      else:
         return N.append(p[:,cthetabin],1)


###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class CALDBPsf(Psf):
   """Handle reading in the parameters from CALDB for multiple formats and implement the PSF functions."""
     
   def __read_params__(self):      
      
      h        = self.CALDBhandles # a tuple with handles to front and back CALDB
      ne,nc    = len(self.e_los),len(self.c_his)
      self.newstyle = 'SCORE' in [x.name for x in h[0][1].get_coldefs()]

      if self.newstyle:

         tables = [[1./(1+N.reshape(h[i][1].data.field('NTAIL'),[nc,ne]).transpose()[:,::-1]),
                    N.reshape(h[i][1].data.field('NTAIL'),[nc,ne]).transpose()[:,::-1]/(1+N.reshape(h[i][1].data.field('NTAIL'),[nc,ne]).transpose()[:,::-1]),
                    N.reshape(h[i][1].data.field('GCORE'),[nc,ne]).transpose()[:,::-1],
                    N.reshape(h[i][1].data.field('GTAIL'),[nc,ne]).transpose()[:,::-1],
                    N.reshape(h[i][1].data.field('SCORE'),[nc,ne]).transpose()[:,::-1],
                    N.reshape(h[i][1].data.field('STAIL'),[nc,ne]).transpose()[:,::-1]] for i in [0,1]]
        
      else:
         tables = [[N.reshape(h[i][1].data.field('GCORE'),[nc,ne]).transpose()[:,::-1],
                    N.reshape(h[i][1].data.field('SIGMA'),[nc,ne]).transpose()[:,::-1]] for i in [0,1]]
      
      self.tables  = N.asarray(tables)
      
   def __call__(self,e,ct,delta, scale_sigma=True, density = False):
      sf = self.scale_func[ct](e) if scale_sigma else 1
      if self.newstyle:
         nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
         yc = self.psf_base(gc,sc*sf,delta,density)
         yt = self.psf_base(gt,st*sf,delta,density)
         return (w * (nc*yc + nt*yt)).sum(axis=1)
      else:
         # note retaining explicit old style for the nonce, for comparison testing
         gc,si,w = self.get_p(e,ct) # note each parameter is an 8-vector
         si *= sf
         us  = 0.5 * N.outer(delta,1./si)**2
         y   = (1-1./gc)*(1+us/gc)**(-gc)
         return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

   def psf_base(self,g,s,delta,density=False):
      """Implement the PSF base function; g = gamma, s = sigma, w = weight, delta = deviation in radians."""
      u = 0.5 * N.outer(delta,1./s)**2
      y = (1-1./g)*(1+u/g)**(-g)
      return y / (TWOPI*s**2 if density else 1.)
 

   """ Not sure these are in use -- if they are, update for old/new CALDB formats.
   def set_cache(self,e,ct):
      np = len(self.get_p(e[0],ct[0])[0])
      self.cache = N.empty([3,len(e),np])
      ca = self.cache; gp = self.get_p; sf = self.scale_func
      for i,(mye,myct) in enumerate(zip(e,ct)):
         ca[:,i,:]  = gp(mye,myct)
         ca[1,i,:] *= sf[myct](mye)

   def get_cache(self,delta,density = False):
      g,s,w = self.cache
      us    = 0.5 * (delta / s)**2
      y     = y  = (1-1./g)*(1+us/g)**(-g)
      return (w*(y / (TWOPI*s**2 if density else 1.))).sum(axis=1)
   """
      
   def integral(self,e,ct,dmax,dmin=0):
      """ Note -- does *not* take a vector argument right now."""
      """ Update to density-style."""
      if self.newstyle:
         nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
         #if scale_sigma: # why no test in old version? think this code may have gotten stale
         sf = self.scale_func[ct](e)
         sc *= sf
         st *= sf
         uc1 = 0.5 * (dmin / sc)**2
         uc2 = 0.5 * (dmax / sc)**2
         ut1 = 0.5 * (dmin / st)**2
         ut2 = 0.5 * (dmax / st)**2
         return (w*( nc*((1+uc1/gc)**(1-gc) - (1+uc2/gc)**(1-gc)) + nt*((1+ut1/gt)**(1-gt) - (1+ut2/gt)**(1-gt)) )).sum()
      else:
         gc,si,w = self.get_p(e,ct)
         si *= self.scale_func[ct](e)
         u1 = 0.5 * (dmin / si)**2
         u2 = 0.5 * (dmax / si)**2
         return  (w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  )).sum()

   def inverse_integral(self,e,ct,percent=68):
      """Return the radius at which the integral PSF reaches the required percentage.
         Pretty primitive, zero error checking, caveat caller.
         
         NB -- returned value is in degrees.
      """
      percent /= 100.
      def f(dmax):
         if dmax <= 0: return -percent
         return self.integral(e,ct,dmax) - percent
      from scipy.optimize import fsolve
      for seed in N.radians([1,0.5,0.1,0.05,0.01]):
         trial = fsolve(f,N.radians(1))
         if trial > 0: break
      return N.degrees(trial)
             
   def band_psf(self,band,weightfunc=None): return BandCALDBPsf(self,band,weightfunc=weightfunc,newstyle=self.newstyle)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class BandPsf(object):

   def init(self):
      self.newstyle = False

   def __init__(self,psf,band,weightfunc=None,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.par     = psf.get_p(band.e,band.ct).copy()
      if weightfunc is not None:
         dom   = N.logspace(N.log10(band.emin),N.log10(band.emax),9)
         wvals = [weightfunc(x) for x in dom]
         svals = psf.scale_func[band.ct](dom)
         num = simps(wvals*svals*dom,x=dom)
         dom = simps(wvals*dom,x=dom)
         self.scale=num/dom
         #def numfunc(loge):
         #   e = N.exp(loge)
         #   return psf.scale_func[band.ct](e) * weightfunc(e) * e
         #def domfunc(loge):
         #   e = N.exp(loge)
         #   return weightfunc(e) * e
         #e0,e1 = N.log([band.emin,band.emax])
         #num = quad(numfunc,e0,e1,epsrel=1e-3,full_output=0,limit=100)[0]
         #dom = quad(domfunc,e0,e1,epsrel=1e-3,full_output=0,limit=100)[0]
         self.scale = num/dom
         self.comp_scale = psf.scale_func[band.ct](band.e)
         print num,dom,num/dom,self.comp_scale
      else:
         self.scale   = psf.scale_func[band.ct](band.e)
      if self.newstyle:
         self.par[4] *= self.scale
         self.par[5] *= self.scale
      else:
         self.par[1] *= self.scale # scale sigma


   # ACHTUNG WITH NEW PSF -- where is this used?  I don't immediately spot a reference.
   #def sigma(self): return (self.par[1]*self.par[-1]).sum()


###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class BandCALDBPsf(BandPsf):

   def psf_base(self,g,s,delta,density=False):
      """Implement the PSF base function; g = gamma, s = sigma, delta = deviation in radians."""
      u = 0.5 * N.outer(delta,1./s)**2
      y = (1-1./g)*(1+u/g)**(-g)
      return y / (TWOPI*s**2 if density else 1.)

   def __call__(self,delta,density=False):
      if self.newstyle:
         nc,nt,gc,gt,sc,st,w = self.par
         yc = self.psf_base(gc,sc,delta,density)
         yt = self.psf_base(gt,st,delta,density)
         return (w*(nc*yc + nt*yt)).sum(axis=1)
      else:
         gc,si,w = self.par
         us  = 0.5 * N.outer(delta,1./si)**2
         y   = (1-1./gc)*(1+us/gc)**(-gc)
         return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

   def integral(self,dmax,dmin=0):
      """ Note -- does *not* take a vector argument right now."""
      if self.newstyle:
         nc,nt,gc,gt,sc,st,w = self.par
         uc1 = 0.5 * (dmin / sc)**2
         uc2 = 0.5 * (dmax / sc)**2
         ut1 = 0.5 * (dmin / st)**2
         ut2 = 0.5 * (dmax / st)**2
         return (w*( nc*((1+uc1/gc)**(1-gc) - (1+uc2/gc)**(1-gc)) + nt*((1+ut1/gt)**(1-gt) - (1+ut2/gt)**(1-gt)) )).sum()
      else:
         gc,si,w = self.par
         u1 = 0.5 * (dmin / si)**2
         u2 = 0.5 * (dmax / si)**2
         return  ( w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)) ).sum()

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class PsfOverlap(object):
   """Routines to calculate how much of the emission of a point source falls onto an ROI."""

   def init(self):
      self.quadrature_tol = 1e-4

   def __init__(self,**kwargs):
      self.init()
      self.__dict__.update(kwargs)

   def __call__(self,band,roi_dir,ps_dir,radius_in_rad=None):
      """Return an array of fractional overlap for a point source at location skydir.
         Note radius arguments are in radians."""

      roi_rad  = radius_in_rad or band.radius_in_rad
      integral = band.psf.integral
      offset   = roi_dir.difference(ps_dir)

      if offset < 1e-5:
         return band.psf.integral(roi_rad) #point source in center of ROI, symmetry

      if offset < roi_rad:

         def interior(x):

            c       = cos(x)
            s2      = (c**2 - 1)
            eff_rad = ( roi_rad**2 + offset**2*(c**2 - 1) )**0.5 - offset*c
            return integral(eff_rad)

         return quad(interior,0,N.pi,epsabs=self.quadrature_tol)[0]/N.pi

      else:

         def exterior(x):

            c    = cos(x)
            s    = (1-c**2)**0.5
            r2   = ( roi_rad**2 - (offset*s)**2 )**0.5
            de   = offset*c - r2
            return integral(dmax=de+2*r2,dmin=de)


         limit = N.arcsin(roi_rad / offset)
         return quad(exterior,0,limit,epsabs=self.quadrature_tol)[0]/N.pi















###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class NewPsf(Psf):

   def __init__(self,CALDB,my_CALDB_pickle,**kwargs):
      self.my_CALDB_pickle = my_CALDB_pickle
      super(NewPsf,self).__init__(CALDB,**kwargs)      

   def __read_params__(self):
      self.p = p = load(file(self.my_CALDB_pickle))

      # tables are in order: NCORE, SIGMA, A, GCORE, GTAIL
      ftables = N.asarray(p['ftables'])[:,:,:-1] # trim of 0-66 deg fit
      btables = N.asarray(p['btables'])[:,:,:-1] 

      # do not fit 24 MeV bin
      self.e_los    = self.e_los[1:]
      self.e_his    = self.e_his[1:]

      self.tables = N.asarray([ftables,btables])
      
   def __call__(self,e,ct,delta, scale_sigma=True, density=False,cthetabin=None):
      """Return the differential psf at given energy and conversion type.

         Parameters
         ----------
         e : float, the energy in MeV
         ct : float, the conversion type (0 = front, 1 = back)
         delta : ndarray(float), the angular separation in radians

         Returns
         -------
         output : ndarray (1d), the differential psf evaluated at the required points
      """
      nc,si,ap,gc,gt,w = self.get_p(e,ct,cthetabin) # note each parameter is an 8-vector
      if scale_sigma: si *= self.scale_func[ct](e)
      us  = 0.5 * N.outer(delta,1./si)**2
      y   = (nc*(1-1./gc))*(1+us/gc)**(-gc) + ((1-nc)*ap**(gt-1))*(1.-1./gt)*(ap+us/gt)**(-gt)
      return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

   def integral(self,e,ct,dmax,dmin=0):
      """ Note -- does *not* take a vector argument right now."""
      nc,si,ap,gc,gt,w = self.get_p(e,ct)
      si *= self.scale_func[ct](e)
      u1 = 0.5 * (dmin / si)**2
      u2 = 0.5 * (dmax / si)**2
      return  (w*(nc*           ( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  ) + \
              (1-nc)*ap**(gt-1)*( (ap+u1/gt)**(1-gt) - (ap+u2/gt)**(1-gt) ) )).sum()

   def band_psf(self,band,weightfunc=None): return BandNewPsf(self,band,weightfunc)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class OldPsf(Psf):
     
   def __read_params__(self):      
      h     = self.CALDBhandles # a tuple with handles to front and back CALDB
      ne,nc = len(self.e_los),len(self.c_his)
      tables = [[N.reshape(h[i][1].data.field('GCORE'),[nc,ne]).transpose()[:,::-1],
                 N.reshape(h[i][1].data.field('SIGMA'),[nc,ne]).transpose()[:,::-1]] for i in [0,1]]
      self.tables  = N.asarray(tables)
      
   def __call__(self,e,ct,delta, scale_sigma=True, density = False):
      gc,si,w = self.get_p(e,ct) # note each parameter is an 8-vector
      if scale_sigma: si *= self.scale_func[ct](e)
      us  = 0.5 * N.outer(delta,1./si)**2
      y   = (1-1./gc)*(1+us/gc)**(-gc)
      return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

   def set_cache(self,e,ct):
      np = len(self.get_p(e[0],ct[0])[0])
      self.cache = N.empty([3,len(e),np])
      ca = self.cache; gp = self.get_p; sf = self.scale_func
      for i,(mye,myct) in enumerate(zip(e,ct)):
         ca[:,i,:]  = gp(mye,myct)
         ca[1,i,:] *= sf[myct](mye)

   def get_cache(self,delta,density = False):
      g,s,w = self.cache
      us    = 0.5 * (delta / s)**2
      y     = y  = (1-1./g)*(1+us/g)**(-g)
      return (w*(y / (TWOPI*s**2 if density else 1.))).sum(axis=1)
      
   def integral(self,e,ct,dmax,dmin=0):
      """ Note -- does *not* take a vector argument right now."""
      gc,si,w = self.get_p(e,ct)
      si *= self.scale_func[ct](e)
      u1 = 0.5 * (dmin / si)**2
      u2 = 0.5 * (dmax / si)**2
      return  (w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  )).sum()
             
   def band_psf(self,band,weightfunc=None): return BandOldPsf(self,band,weightfunc)


###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class BandNewPsf(BandPsf):

   def __call__(self,delta,density=False):
      nc,si,ap,gc,gt,w = self.par
      us  = 0.5 * N.outer(delta,1./si)**2
      y   = (nc*(1-1./gc))*(1+us/gc)**(-gc) + ((1-nc)*ap**(gt-1))*(1.-1./gt)*(ap+us/gt)**(-gt)
      return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

   def integral(self,dmax,dmin=0):
      """ Note -- does *not* take a vector argument right now."""
      nc,si,ap,gc,gt,w = self.par
      u1 = 0.5 * (dmin / si)**2
      u2 = 0.5 * (dmax / si)**2
      return  (w*(nc*           ( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  ) + \
              (1-nc)*ap**(gt-1)*( (ap+u1/gt)**(1-gt) - (ap+u2/gt)**(1-gt) ) )).sum()


###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class BandOldPsf(BandPsf):

   def __call__(self,delta,density=False):
      gc,si,w = self.par
      us  = 0.5 * N.outer(delta,1./si)**2
      y   = (1-1./gc)*(1+us/gc)**(-gc)
      return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

   def integral(self,dmax,dmin=0):
      """ Note -- does *not* take a vector argument right now."""
      gc,si,w = self.par
      u1 = 0.5 * (dmin / si)**2
      u2 = 0.5 * (dmax / si)**2
      return  ( w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)) ).sum()