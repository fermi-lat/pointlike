"""Contains miscellaneous classes for background and exposure management.
    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/pointspec_helpers.py,v 1.49 2011/11/28 22:27:40 kerrm Exp $

    author: Matthew Kerr
    """

import numpy as N
from skymaps import SkyDir,CompositeSkySpectrum,DiffuseFunction,EffectiveArea,Exposure,IsotropicSpectrum,IsotropicPowerLaw
from uw.like.Models import Model,Constant,PowerLaw,ExpCutoff,DefaultModelValues,LogParabola
from roi_diffuse import DiffuseSource,ROIDiffuseModel_OTF
from roi_extended import ExtendedSource,ROIExtendedModel
from SpatialModels import Disk,Gaussian,EllipticalDisk,EllipticalGaussian,SpatialMap,GAUSSIAN_X68
from os.path import join
import os
from abc import abstractmethod
from uw.utilities import keyword_options

class PointSource(object):
    """ combine name, skydir, model """
    def __init__(self,skydir,name,model=None,free_parameters=True,leave_parameters=False):
        self.name    = name
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
                                'R.A. (J2000):\t%.5f'%(self.skydir.ra()),
                                'Dec. (J2000):\t%.5f'%(self.skydir.dec()),
                                'Model:\t\t%s'%(self.model.full_name()),
                                '\t'+self.model.__str__(indent='\t'), 
                                ])

    def copy(self):
        """ Create a deep copy of the point source. """
        return PointSource(SkyDir(self.skydir.ra(),self.skydir.dec()),
                           self.name,self.model.copy(),leave_parameters=True)

class Singleton(object):
    """Implement a singleton class."""

    __instances = {}

    def __init__(self,constructor,key=None,*args,**kwargs):
        """Constructor -- a class name; note -- not a string!
            Key            -- an optional key; constructor name used otherwise.
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

class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime."""

    def __init__(self,sa):

        EffectiveArea.set_CALDB(sa.CALDBManager.CALDB)
        
        Exposure.set_cutoff(N.cos(N.radians(sa.thetacut)))
        inst = ['front', 'back']
        aeff_files = sa.CALDBManager.get_aeff()
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


# For backwards compatability
from uw.like.roi_catalogs import PointSourceCatalog,FermiCatalog,ExtendedSourceCatalog,CatalogManager

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
            raise Exception,'Could not find the FITS file specified for MapCubeFunction (file = %s).' % spatialModelFile
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
            index= 10**spectralModel._p[1]
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

def get_default_diffuse_mapper(sa,roi_dir,quiet):
    """ Returns a function which maps a roi_diffuse.DiffuseSource
        object to its corresponding roi_diffuse.ROIDiffuseModel object.

        The default mapping is that an extended source is mapped into
        an ROIExtendedModel object using ROIExtendedModel's factory
        function. Otherwise, an ROIDiffuseModel_OTF object is returned. """
    return lambda x: ROIExtendedModel.factory(sa,x,roi_dir,quiet=quiet) \
                     if isinstance(x,ExtendedSource) \
                     else ROIDiffuseModel_OTF(sa,x,roi_dir,quiet=quiet)

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
