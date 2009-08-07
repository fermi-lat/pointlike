
import numpy as N
from skymaps import SkyDir
from Models import PowerLaw


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

   def get_sources(self,skydir,radius=15.):
    
      diffs   = N.asarray([skydir.difference(d) for d in self.dirs])*180/N.pi
      mask    = ((diffs < radius)&(self.ts > self.min_ts)) & \
                ((self.fluxes > self.min_flux)|(diffs < self.max_distance))         
      diffs   = diffs[mask]
      sorting = N.argsort(diffs)

      #sort SkyDirs -- pesky vector behavior...
      dirs    = [x for i,x in enumerate(self.dirs) if mask[i]]
      dirs    = [dirs[x] for x in sorting]

      names   = self.names[mask][sorting]
      models  = self.models[mask][sorting]

      point_sources = map(PointSource,dirs,names,models,[False]*len(models))
      return point_sources


   def merge_lists(self,skydir,radius=15,user_list=None):
      """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence."""
      
      #implementation not very efficient, but works... can come back if it's slow.

      cat_list = self.get_sources(skydir,radius)
      if user_list is None: return cat_list

      from collections import deque
      merged_list = deque()
      
      for ups in user_list:
         merged_list.append(ups)
         for cps in cat_list:
            if ups.skydir.difference(cps.skydir)*180/N.pi < self.prune_radius:
               cps.duplicate = True
      
      for cps in cat_list:
         if not cps.duplicate: merged_list.append(cps)
      
      merged_list = list(merged_list)

      def compare(ps1,ps2):
         d1 = ps1.skydir.difference(skydir)
         d2 = ps2.skydir.difference(skydir)
         return cmp(d1,d2)

      merged_list.sort(cmp=compare)

      return merged_list
