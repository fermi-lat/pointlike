"""Contains miscellaneous classes for background and exposure management.
   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pointspec_helpers.py,v 1.18 2010/08/05 14:56:45 wallacee Exp $

   author: Matthew Kerr
   """

import numpy as N
from skymaps import *
from uw.like.Models import Model,Constant,PowerLaw,DefaultModelValues
from roi_managers import ROIBackgroundModel
from roi_diffuse import DiffuseSource
from roi_extended import ExtendedSource
from os.path import join
import os

###====================================================================================================###

class PointSource(object):
   def __init__(self,skydir,name,model=None,free_parameters=True,leave_parameters=False):
      self.name   = name
      self.skydir = skydir
      self.model  = PowerLaw() if model is None else model
      #if not free_parameters:
      if not leave_parameters:
         for i in xrange(len(self.model.free)): self.model.free[i] = free_parameters
      self.duplicate = False
   def __str__(self):
      return '\n'.join(['\n',
                        '='*60,
                        'Name:\t\t%s'%(self.name),
                        'R.A. (J2000):\t\t%.5f'%(self.skydir.ra()),
                        'Dec. (J2000):\t\t%.5f'%(self.skydir.dec()),
                        'Model:\t\t%s'%(self.model.name),
                        '\n',
                        'Model Parameters:',
                        '%s'%(str(self.model))])

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
      return inst[key]

   def __call__(self,key):
      inst = Singleton._Singleton__instances
      key  = str(key)
      if key in inst: return inst[key]

class Singleton2(Singleton):

    def __init__(self):
        pass

###====================================================================================================###


class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime."""

    def __init__(self,sa):

        caldb = sa.CALDB
        if not os.path.exists(os.path.join(caldb, 'bcf')):
            caldb = os.path.join(caldb, 'data', 'glast','lat')
            if not os.path.exists(os.path.join(caldb, 'bcf')):
                raise Exception('ExposureManager attempt to set invalid CALDB: "%s"' % sa.CALDB)
        EffectiveArea.set_CALDB(caldb)
        
        Exposure.set_cutoff(N.cos(N.radians(sa.thetacut)))
        inst = ['front', 'back']
        aeff_files = [os.path.join(caldb, 'bcf','ea', 'aeff_%s_%s.fits'%(sa.irf,fb)) for fb in inst] 
        ok = [os.path.exists(file) for file in aeff_files]
        if not all(ok):
            raise Exception('one of CALDB aeff files not found: %s' %aeff_files)
        self.ea  = [EffectiveArea('', file) for file in aeff_files]
        if sa.verbose: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        if sa.use_weighted_livetime:
            self.exposure = [Exposure(sa.pixeldata.lt,sa.pixeldata.weighted_lt,ea) for ea in self.ea]
        else:
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

   def __init__(self,analysis_environment,background='1FGL',**kwargs):

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

   def make_1FGL(self):

      singl = Singleton(DiffuseFunction,'gf',join(self.diffdir,'gll_iem_v02_P6_v3_diff_front.fits'))
      singl.add(        DiffuseFunction,'gb',join(self.diffdir,'gll_iem_v02_P6_v3_diff_back.fits'))
      singl.add(        IsotropicSpectrum,'if',join(self.diffdir,'isotropic_iem_front_v02.txt'))
      singl.add(        IsotropicSpectrum,'ib',join(self.diffdir,'isotropic_iem_back_v02.txt'))

      self.dmodels = [ [singl('gf'),singl('gb')], [singl('if'),singl('ib')] ]
      self.smodels = [ PowerLaw(p=[1,1],free=[True,True],index_offset=1),
                       Constant(free=[True]) ]
      self.names   = ['gll_iem_v02', 'Isotropic Diffuse']
      self.cat     = join(self.catdir,r'gll_psc_v03.fit')  # thb: this has one slight change from v02, involving a name


   def get_bgmodels(self, models = None, lat = None):

      gal_index = True if (lat is not None and abs(lat) < 20) else False
      self.smodels[0].free[1] = gal_index # need to be careful about this
      smodels = self.smodels if models is None else models
      return map(ROIBackgroundModel,self.dmodels,smodels,self.names)


###====================================================================================================###

class PointSourceCatalog(object):
    """ Define an interface for point source catalog that can be used to
        construct the model for an ROI analysis."""

    def get_sources(self,skydir,radius=15):
        raise NotImplementedError,'Virtual'

    def merge_lists(self,skydir,radius=15,user_point_list=None,user_diffuse_list=None):
        raise NotImplementedError,'Virtual'

###====================================================================================================###

class FermiCatalog(PointSourceCatalog):
   """Read a Jean Ballet catalogue and use it to set up a source list for a ROI.
   
      This class works (mostly) for EMS and 1FGL styles."""

   def init(self):
      self.prune_radius  = 0.10 #deg; in a merge, consider sources closer than this duplicates
      self.free_radius   = 2 #deg; sources within this distance have free spectral parameters
      self.min_flux      = 2e-9 #ph/cm2/s; minimum flux for sources beyond a certain proximity
      self.max_distance  = 5 #deg; distance inside which sources are returned regardless of flux
      self.min_ts        = 25

   def __init__(self,catalog_file,*args,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.__open_catalog__(catalog_file)

   def __open_catalog__(self,catalog_file):
      from pyfits import open
      f = open(catalog_file)
      colnames = [x.name for x in f[1].get_coldefs()]
      sname    = 'NickName' if 'NickName' in colnames else 'Source_Name'
      ras  = f[1].data.field('RA')
      decs = f[1].data.field('DEC')
      pens = f[1].data.field('PIVOT_ENERGY')
      n0s  = f[1].data.field('FLUX_DENSITY')
      inds = f[1].data.field('SPECTRAL_INDEX')
      #ts   = f[1].data.field('TEST_STATISTIC')
      inds = N.where(inds > 0, inds, -inds)

      self.dirs   = map(SkyDir,N.asarray(ras).astype(float),N.asarray(decs).astype(float))
      self.models = N.asarray([PowerLaw(p=[n0,ind],e0=pen) for n0,ind,pen in zip(n0s,inds,pens)])
      #self.fluxes = N.asarray(f[1].data.field('FLUX100'))
      self.names  = N.chararray.strip(f[1].data.field(sname))
      #self.ts     = N.asarray(ts)

      f.close()

   def get_sources(self,skydir,radius=15):

      diffs   = N.degrees(N.asarray([skydir.difference(d) for d in self.dirs]))
      #mask    = ((diffs < radius)&(self.ts > self.min_ts)) & \
      #          ((self.fluxes > self.min_flux)|(diffs < self.max_distance))
      mask    = diffs < radius
      diffs   = diffs[mask]
      sorting = N.argsort(diffs)

      #sort SkyDirs -- pesky vector behavior...
      dirs    = [x for i,x in enumerate(self.dirs) if mask[i]]
      dirs    = [dirs[x] for x in sorting]

      names   = self.names[mask][sorting]
      models  = [x.copy() for x in self.models[mask][sorting]]

      fm = diffs[sorting] < self.free_radius
      point_sources = map(PointSource,dirs,names,models,fm)
      return point_sources


   def merge_lists(self,skydir,radius=15,user_point_list=[],user_diffuse_list=[]):
      """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence."""

      user_extended_list = [i for i in user_diffuse_list if isinstance(i,ExtendedSource)]

      cat_list = self.get_sources(skydir,radius)
      if user_point_list==[] and user_extended_list==[]: return cat_list,[]

      from collections import deque
      merged_list = deque(user_point_list)


      for ncps,cps in enumerate(cat_list):
         merged_list.append(cps)
         for ups in user_point_list:
            if N.degrees(ups.skydir.difference(cps.skydir)) < self.prune_radius:
               merged_list.pop(); break

         for ues in user_extended_list:
            if N.degrees(ues.spatial_model.center.difference(cps.skydir)) < self.prune_radius:
               merged_list.pop(); break

      merged_list = list(merged_list)
      merged_list.sort(key = lambda ps: ps.skydir.difference(skydir))

      return merged_list,user_diffuse_list

class CatalogManager(FermiCatalog):
    """ For compatibility.  Temporary."""
    pass


def get_diffuse_source(spatialModel='ConstantValue',
                       spatialModelFile=None,
                       spectralModel='PowerLaw',
                       spectralModelFile=None,
                       name=None,
                       diffdir = None):

    """ Return a DiffuseSource instance suitable for
        instantiating a child of ROIDiffuseModel.

        NB -- don't support front/back distinction atm.

        The list of supported models is currently very short, but covers
        the usual cases for modeling diffuse backgrounds.  Additional
        use cases can be developed on an ad hoc basis.
        
        Arguments:
        
        spatialModel -- an XML-style keyword.  Valid options are
                        1) ConstantValue (isotropic)
                        2) MapCubeFunction (from a FITS file)
                        
        spatialModelFile -- if a mapcube is specified, its location
        
        spectralModel -- This can be either an XML-style keyword or an
                         instance of Model.
                         If an XML-style keyword, valid options are
                         1) FileFunction
                         2) PowerLaw
                         3) Constant
                        
        spectralModelFile -- if a tabular function is specified,
                             its location

        name -- a name for the ol' model

        
        diffdir -- if the XML files specify paths relative to some
                   directory, set this variable appropriately
    """

    if (diffdir is not None):
        if spatialModelFile is not None:
            spatialModelFile = os.path.join(diffdir,spatialModelFile)
        if spectralModelFile is not None:
            spectralModelFile = os.path.join(diffdir,spectralModelFile)

    # check input sanity
    if not isinstance(spectralModel,Model):
        if (spectralModelFile is not None):
            if not os.path.exists(spectralModelFile):
                raise Exception,'Could not find the ASCII file specified for FileFunction'
        elif spectralModel != 'PowerLaw':
            raise NotImplementedError,'Must provide one of the understood spectral models.'
        else:
            pass

    if spatialModel=='MapCubeFunction':
        if (spatialModelFile is None) or (not os.path.exists(spatialModelFile)):
            raise Exception,'Could not find the FITS file specified for MapCubeFunction.'
    elif spatialModel != 'ConstantValue':
        raise NotImplementedError,'Must provide one of the understood spatial models.'
    else:
        pass              

    ston = Singleton2()
    dmodel = None; smodel = None

    # deal with isotropic models
    if spatialModel=='ConstantValue':
        if spectralModelFile is not None:
            dmodel = IsotropicSpectrum(spectralModelFile)
            if isinstance(spectralModel,Model):
                smodel = spectralModel
            else:
                smodel = Constant()
        elif spectralModel == 'PowerLaw':
            # use Sreekumar-like defaults
            dmodel = IsotropicPowerLaw(1.5e-5,2.1)
            smodel = PowerLaw(p=[1,1],index_offset=1)
        else:
            # interpret the Model object as power law
            flux = spectralModel.i_flux(emin=100)
            index= 10**spectralModel.p[1]
            dmodel = IsotropicPowerLaw(flux,index)
            smodel = PowerLaw(p=[1,1],index_offset=1)

    # deal with mapcubes
    else:
        if spectralModel == 'FileFunction':
            dmodel1 = IsotropicSpectrum(spectralModelFile)
            dmodel2 = ston.add(DiffuseFunction,spatialModelFile,spatialModelFile)
            dmodel  = CompositeSkySpectrum(dmodel1,dmodel2)
            dmodel.saveme1 = dmodel1; dmodel.saveme2 = dmodel2
            smodel  = Constant()
        else:
            dmodel = ston.add(DiffuseFunction,spatialModelFile,spatialModelFile)
            if spectralModel == 'PowerLaw':
                smodel = PowerLaw(p=[1,1],index_offset=1)
            elif spectralModel == 'Constant':
                smodel = Constant()
            else:
                smodel = spectralModel

    if (dmodel is None) or (smodel is None):
        raise Exception,'Was unable to parse input.'

    return DiffuseSource(dmodel,smodel,name)


def get_default_diffuse(diffdir=None,gfile='gll_iem_v02.fit',ifile='isotropic_iem_v02.txt'):
    """ Try to get defaults for the diffuse background sources, assumed to be 
        a MapCubeFunction/Powerlaw and ConstantValue/FileFunction
        Setting gfile or ifile to None ignores that entry."""
    if diffdir is None:
        diffdir = join(os.environ['EXTFILESSYS'],'galdiffuse')
    dsources = []
    if gfile is not None:
        gfilex = join(diffdir,gfile)
        if not os.path.exists(gfilex):
            raise Exception,' Galactic diffuse file "%s" not found.' %gfilex
        else:
            dsources += [get_diffuse_source('MapCubeFunction',gfilex,'PowerLaw',None,'Galactic Diffuse (%s)'%gfile)]
    if ifile is not None:
        ifilex = join(diffdir,ifile)
        if not os.path.exists(ifilex):
            raise Exception,'isotropic diffuse file "%s" not found.'%ifilex
        else:
            dsources += [get_diffuse_source('ConstantValue',None,'FileFunction',ifilex,'Isotropic Diffuse (%s)'%ifile)]

    if len(dsources) == 0:
        raise Exception,'Unable to find any diffuse sources.'
    return dsources
