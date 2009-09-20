"""Contains miscellaneous classes for background and exposure management.
   
   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointspec_helpers.py,v 1.2 2009/09/20 20:43:19 kerrm Exp $

   author: Matthew Kerr
   """

import numpy as N
from skymaps import *
from Models import Constant,PowerLaw
from roi_managers import ROIBackgroundModel
from os.path import join

###====================================================================================================###


class Singleton(object):
   """Implement a singleton class."""

   __instances = {}

   def __init__(self,constructor,key=None,*args,**kwargs):
      """Constructor -- a class name; note -- not a string!
         Key         -- an optional key; constructor name used otherwise.
                        **keys are unique, not constructors!**
         Additional args and kwargs are passed to constructor."""

      self.add(constructor,key,*args,**kwargs)

   def add(self,constructor,key=None,*args,**kwargs):
      inst = Singleton._Singleton__instances
      key  = str(constructor) if key is None else key
      if key not in inst.keys(): inst[key] = constructor(*args,**kwargs)

   def __call__(self,key):
      inst = Singleton._Singleton__instances
      key  = str(key)
      if key in inst: return inst[key]

###====================================================================================================###


class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime."""
   
    def __init__(self,sa):

        EffectiveArea.set_CALDB(sa.ae.CALDB)
        Exposure.set_cutoff(N.cos(N.radians(sa.thetacut)))
        inst = ['front', 'back']
        self.ea  = [EffectiveArea(sa.irf+'_'+x) for x in inst]
        if sa.verbose: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        self.exposure = [Exposure(sa.pixeldata.lt,ea) for ea in self.ea]

    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)

###====================================================================================================###


class ConsistentBackground(object):
   """Manage the construction of a consistent background model.

      *Notes*
      
      A background model comprises diffuse and point sources.  While
      some parameters may be refit in a spectral analysis, the default
      background model should give a good representation.
      
      The class handles matching LAT source lists with the diffuse models
      used to generate them.  The combination is a "consistent"
      background.

   """

   def __init__(self,analysis_environment,background='ems_ring',**kwargs):

      self.diffdir    = analysis_environment.diffdir
      self.catdir     = analysis_environment.catdir
      self.background = background

      eval('self.make_%s()'%background)

      if 'cat' in kwargs: cat = kwargs['cat']
      else:               cat = self.cat
      singl = Singleton(CatalogManager,'cat',cat)
      self.cm = singl('cat')

   def make_nms_galprop(self):

      singl = Singleton(DiffuseFunction,'gf',join(self.diffdir,'gas_mapcube_54_77Xvarh7S_P6_v3_diff_front.fits'))
      singl.add(        DiffuseFunction,'gb',join(self.diffdir,'gas_mapcube_54_77Xvarh7S_P6_v3_diff_back.fits'))
      singl.add(        DiffuseFunction,'if',join(self.diffdir,'ics_isotropic_mapcube_54_77Xvarh7S_P6_v3_diff_front.fits'))
      singl.add(        DiffuseFunction,'ib',join(self.diffdir,'ics_isotropic_mapcube_54_77Xvarh7S_P6_v3_diff_back.fits'))
      singl.add(        IsotropicSpectrum,'iso',join(self.diffdir,'Total_b30_EGBfree.txt'))

      self.dmodels = [ [singl('gf'),singl('gb')] , singl('iso') , [singl('if'),singl('ib')] ]
      self.smodels = [ PowerLaw(p=[1,1],free=[True,True],index_offset=1),
                       Constant(free=[True]), Constant(free[False]) ]
      self.names   = ['Gas Galactic Diffuse', 'Isotropic Diffuse', 'IC Galactic Diffuse']
      self.cat     = join(self.catdir,r'gll_psc9month_v2r2.fit')

   def make_ems_ring(self):

      singl = Singleton(DiffuseFunction,'gf',join(self.diffdir,'gll_iem_v02_P6_v3_diff_front.fits'))
      singl.add(        DiffuseFunction,'gb',join(self.diffdir,'gll_iem_v02_P6_v3_diff_back.fits'))
      singl.add(        IsotropicSpectrum,'if',join(self.diffdir,'isotropic_iem_front_v02.txt'))
      singl.add(        IsotropicSpectrum,'ib',join(self.diffdir,'isotropic_iem_back_v02.txt'))

      self.dmodels = [ [singl('gf'),singl('gb')], [singl('if'),singl('ib')] ]
      self.smodels = [ PowerLaw(p=[1,1],free=[True,True],index_offset=1),
                       Constant(free=[True]) ]
      self.names   = ['gll_iem_v02', 'Isotropic Diffuse']
      self.cat     = join(self.catdir,r'gll_psc11month_v2.fit')


   def get_bgmodels(self, models = None, lat = None):

      gal_index = True if (lat is not None and abs(lat) < 20) else False
      self.smodels[0].free[1] = gal_index # need to be careful about this

      return map(ROIBackgroundModel,self.dmodels,self.smodels,self.names)


###====================================================================================================###


class PointSource(object):
   def __init__(self,skydir,name,model=None,free_parameters=True):
      self.name   = name
      self.skydir = skydir
      self.model  = PowerLaw() if model is None else model
      if not free_parameters:
         for i in xrange(len(self.model.free)): self.model.free[i] = False
      self.duplicate = False
   def __str__(self):
      return '\n'.join(['\n',
                        ''.join(['=']*60),
                        'Name:\t\t%s'%(self.name),
                        'R.A. (J2000):\t\t%.5f'%(self.skydir.ra()),
                        'Dec. (J2000):\t\t%.5f'%(self.skydir.dec()),
                        'Model:\t\t%s'%(self.model.name)])

###====================================================================================================###


class CatalogManager(object):
   """Read a Jean Ballet catalogue and use it to set up a source list for a ROI."""

   def init(self):
      self.prune_radius  = 0.10 #deg; in a merge, consider sources closer than this duplicates
      self.min_flux      = 1e-8 #ph/cm2/s; minimum flux for sources beyond a certain proximity
      self.max_distance  = 5 #deg; distance inside which sources are returned regardless of flux
      self.min_ts        = 25

   def __init__(self,catalog_file,*args,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.__open_catalog__(catalog_file)

   def __open_catalog__(self,catalog_file):
      from pyfits import open
      f = open(catalog_file)
      ras  = f[1].data.field('RA')
      decs = f[1].data.field('DEC')
      pens = f[1].data.field('PIVOT_ENERGY')
      n0s  = f[1].data.field('FLUX_DENSITY')
      inds = f[1].data.field('SPECTRAL_INDEX')
      ts   = f[1].data.field('TEST_STATISTIC')
      inds = N.where(inds > 0, inds, -inds)

      self.dirs   = map(SkyDir,N.asarray(ras).astype(float),N.asarray(decs).astype(float))
      self.models = N.asarray([PowerLaw(p=[n0,ind],e0=pen) for n0,ind,pen in zip(n0s,inds,pens)])
      self.fluxes = N.asarray(f[1].data.field('FLUX100'))
      self.names  = N.asarray(f[1].data.field('NickName'))
      self.ts     = N.asarray(ts)

      f.close()

   def get_sources(self,skydir,radius=15.,free_radius=2.):
    
      diffs   = N.degrees(N.asarray([skydir.difference(d) for d in self.dirs]))
      mask    = ((diffs < radius)&(self.ts > self.min_ts)) & \
                ((self.fluxes > self.min_flux)|(diffs < self.max_distance))         
      diffs   = diffs[mask]
      sorting = N.argsort(diffs)

      #sort SkyDirs -- pesky vector behavior...
      dirs    = [x for i,x in enumerate(self.dirs) if mask[i]]
      dirs    = [dirs[x] for x in sorting]

      names   = self.names[mask][sorting]
      models  = self.models[mask][sorting]

      fm = diffs[sorting] < free_radius
      point_sources = map(PointSource,dirs,names,models,fm)
      return point_sources


   def merge_lists(self,skydir,radius=15,user_list=None):
      """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence."""

      cat_list = self.get_sources(skydir,radius)
      if user_list is None: return cat_list

      from collections import deque
      merged_list = deque(user_list)

      for ncps,cps in enumerate(cat_list):
         merged_list.append(cps)
         for ups in user_list:
            if N.degrees(ups.skydir.difference(cps.skydir)) < self.prune_radius:
               merged_list.pop(); break

      merged_list = list(merged_list)
      merged_list.sort(key = lambda ps: ps.skydir.difference(skydir))

      return merged_list
