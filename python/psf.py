import numpy as N
from math import log
from skymaps import IParams

class PSF:

   def __init__(self,**kwargs):

      self.init()
      self.__dict__.update(kwargs)       

      if self.CALDB is None:
         try:
            from os import environ
            IParams.set_CALDB(environ['CALDB'])
         except:
               IParams.set_CALDB('f:\\glast\\caldb\\v0r7p1\\CALDB\\data\\glast\\lat')
      else:
         IParams.set_CALDB(self.CALDB)

      IParams.init('_'.join(self.psf_irf.split('_')[:-1]),self.psf_irf.split('_')[-1])


   def init(self):
      self.CALDB = None
      self.psf_irf = 'P6_v3_diff'
      self.umax = 50
      self.cache = False
      self.quadrature_tol = 1e-5
   
   
   def sigma(self,e,event_class):
      return (180./N.pi)*N.asarray([IParams.sigma(ie,int(iec)) for ie,iec in zip(N.asarray(e).flatten(),N.asarray(event_class).flatten())])

   def gamma(self,e,event_class):
      return N.asarray([IParams.gamma(ie,int(iec)) for ie,iec in zip(N.asarray(e).flatten(),N.asarray(event_class).flatten())])

   def conf_region(self,e,event_class, percentage=68,use_u = False):
      "Return the _percentage_ confidence region radius for the provided energy.  Unit is degrees."""
      gamma = self.gamma(e,event_class)
      u = gamma*((1-percentage/100.)**(1./(1-gamma))-1)
      if use_u: return u
      return self.u2deg(e,u,event_class)

   def u2deg(self,e,u,event_class):
      return self.sigma(e,event_class)*(2*u)**0.5

   def deg2u(self,e,deg,event_class):
      return 0.5*(deg/self.sigma(e,event_class))**2

   def rad2u(self,e,rad,event_class):
      return self.deg2u(e,rad*180/N.pi,event_class)
      
   def fmax(self,e,u,event_class,units='u'):
      if units == 'rad':
         u = self.rad2u(e,u,event_class)
      elif units == 'deg':
         u = self.deg2u(e,u,event_class)
      g = self.gamma(e,event_class)
      return 1-(1+u/g)**(1-g)

   def set_cache(self,e,event_class):
      self.cache_sigma = self.sigma(e,event_class)
      self.cache_gamma = self.gamma(e,event_class)
      self.cache = True

   def __call__(self,e,delta,event_class,radians=True,density=False):
      """Return f(u).  Input is in radians if true."""

      if self.cache:
         sigma = self.cache_sigma
         gamma = self.cache_gamma

      else:
         sigma = self.sigma(e,event_class)
         gamma = self.gamma(e,event_class)

      units = 180 / N.pi if radians else 1.

      u = 0.5 * ((delta * units)/sigma)**2
      jac = (2*N.pi/units**2)*(sigma**2) if density else 1.
      return ( (1.-1./gamma)*(1+u/gamma)**-gamma )/jac 

   def overlap(self,d1,d2,e,event_class,roi_radius):

      sigma = self.sigma(e,event_class) * N.pi / 180 #note everything in radians
      gamma = self.gamma(e,event_class)

      offset = d1.difference(d2)

      if offset < 1e-5:
         return 1-(1+0.5*(roi_radius/sigma)**2/gamma)**(1-gamma)

      from scipy.integrate import quad
      from math import cos

      if offset < roi_radius:

         def interior(x):

            c  = cos(x)
            eff_rad = ( roi_radius**2 + offset**2*(c**2 - 1) )**0.5 - offset*c
            u = 0.5 * (eff_rad / sigma)**2

            return 1.- (1. + u/gamma )**(1-gamma)

         result = quad(interior,0,N.pi,epsabs=self.quadrature_tol,epsrel=self.quadrature_tol)[0]/N.pi

      else:

         def exterior(x):

            c    = cos(x)
            s    = (1-c**2)**0.5

            r2 = ( roi_radius**2 - (offset*s)**2 )**0.5
            de = offset*c - r2

            u1 = 0.5 * (de/sigma)**2
            u2 = 0.5 * ((de + 2*r2)/sigma)**2

            return (1. + u1/gamma )**(1-gamma) - (1. + u2/gamma )**(1-gamma)

         limit = N.arcsin(roi_radius / offset)
         result = quad(exterior,0,limit,epsabs=self.quadrature_tol,epsrel=self.quadrature_tol)[0]/N.pi

      return result