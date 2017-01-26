"""Contains miscellaneous classes for background and exposure management.
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pointspec_helpers.py,v 1.61 2013/03/21 20:09:07 lande Exp $

    author: Matthew Kerr
    """

import numpy as N
from skymaps import SkyDir,CompositeSkySpectrum,DiffuseFunction,EffectiveArea,Exposure,IsotropicSpectrum,IsotropicConstant
from uw.like.Models import Model,Constant,PowerLaw,ScalingPowerLaw,ExpCutoff,LogParabola,FileFunction
from . roi_diffuse import DiffuseSource,ROIDiffuseModel_OTF
from . roi_extended import ExtendedSource,ROIExtendedModel
from os.path import join
import os
from abc import abstractmethod
from uw.utilities import keyword_options
from uw.utilities import path

class PointlikeException(Exception): pass

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

        #self.ea  = [EffectiveArea('', file) for file in aeff_files]
        # account for new 'FB'
        self.ea  = [EffectiveArea('', file, 'EFFECTIVE AREA_'+fb.upper()) for file,fb in zip(aeff_files,inst)]

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
            if not os.path.exists(path.expand(spectralModelFile)):
                raise Exception('Could not find the ASCII file specified for FileFunction')
        elif not (spectralModel == 'PowerLaw' or spectralModel == 'Constant'):
            raise NotImplementedError,'Must provide one of the understood spectral models.'
        else:
            pass

    if spatialModel=='MapCubeFunction':
        if (spatialModelFile is None) or (not os.path.exists(path.expand(spatialModelFile))):
            raise Exception('Could not find the FITS file specified for MapCubeFunction (file = %s).' % spatialModelFile)
    elif spatialModel != 'ConstantValue':
        raise NotImplementedError,'Must provide one of the understood spatial models.'
    else:
        pass                  

    ston = Singleton2()
    dmodel = None; smodel = None


    # deal with isotropic models
    if spatialModel=='ConstantValue':
        if isinstance(spectralModel,Model):
            smodel=spectralModel
            dmodel=IsotropicConstant()
        elif spectralModelFile is not None:
            smodel = FileFunction(normalization=1, file=spectralModelFile)
            dmodel = IsotropicConstant()
        elif spectralModel == 'PowerLaw':
            # use Sreekumar-like defaults
            smodel = PowerLaw(index=2.1)
            smodel.set_flux(1.5e-5, emin=100, emax=N.inf)

            dmodel = IsotropicConstant()
        else:
            raise Exception("Unable to parse input.")

    # deal with mapcubes
    else:
        if spectralModel == 'FileFunction':
            dmodel1 = IsotropicSpectrum(spectralModelFile)
            dmodel2 = ston.add(DiffuseFunction,spatialModelFile,spatialModelFile)
            dmodel  = CompositeSkySpectrum(dmodel1,dmodel2)
            dmodel.saveme1 = dmodel1; dmodel.saveme2 = dmodel2
            smodel  = Constant()
        else:
            dmodel = ston.add(DiffuseFunction,path.expand(spatialModelFile),path.expand(spatialModelFile))
            dmodel.filename=spatialModelFile
            if spectralModel == 'PowerLaw':
                smodel = ScalingPowerLaw()
            elif spectralModel == 'Constant':
                smodel = Constant()
            else:
                smodel = spectralModel

    if (dmodel is None) or (smodel is None):
         raise Exception('Was unable to parse input.')

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

def get_default_diffuse(diffdir=None,gfile='gll_iem_v02.fit',ifile='isotropic_iem_v02.txt', limit_parameters=False):
    """ Try to get defaults for the diffuse background sources, assumed to be 
        a MapCubeFunction/Powerlaw and ConstantValue/FileFunction
        Setting gfile or ifile to None ignores that entry.
            
            >>> gal, iso = get_default_diffuse(diffdir='$(GLAST_EXT)/diffuseModels/v2r0p1/',
            ...     gfile="ring_2year_P76_v0.fits",
            ...     ifile="isotrop_2year_P76_source_v1.txt")
            >>> print gal.smodel.name
            ScalingPowerLaw
            >>> print gal.smodel['norm']
            1.0
            >>> print gal.smodel['index']
            0.0
            >>> print gal.smodel.get_mapper('norm')
            <class 'uw.utilities.parmap.LogMapper'>
            >>> print gal.smodel.get_mapper('index')
            <class 'uw.utilities.parmap.LinearMapper'>

            >>> print iso.smodel.name
            FileFunction
            >>> print iso.smodel['Normalization']
            1.0
            >>> print iso.smodel.get_mapper('Normalization')
            <class 'uw.utilities.parmap.LogMapper'>

            >>> gal, iso = get_default_diffuse(diffdir='$(GLAST_EXT)/diffuseModels/v2r0p1/',
            ...     gfile="ring_2year_P76_v0.fits",
            ...     ifile="isotrop_2year_P76_source_v1.txt",
            ...     limit_parameters=True)
            >>> print gal.smodel.get_mapper('norm')
            LimitMapper(0.1,10,1)
            >>> print gal.smodel.get_mapper('index')
            LimitMapper(-1,1,1)
            >>> print iso.smodel.get_mapper('Normalization')
            LimitMapper(0.1,10,1)



     """
    if diffdir is None:
        diffdir = join(os.environ['EXTFILESSYS'],'galdiffuse')

    if gfile is None and ifile is None:
        raise Exception('Unable to find any diffuse sources.')

    if gfile is not None:
        gfilex = join(diffdir,gfile)
        if not os.path.exists(path.expand(gfilex)):
            raise Exception(' Galactic diffuse file "%s" not found.' %gfilex)
        else:
            gfile = get_diffuse_source('MapCubeFunction',gfilex,'PowerLaw',None,'Galactic Diffuse (%s)'%gfile)
            if limit_parameters: gfile.smodel.set_default_limits()
    if ifile is not None:
        ifilex = join(diffdir,ifile)
        if not os.path.exists(path.expand(ifilex)):
            raise Exception('isotropic diffuse file "%s" not found.'%ifilex)
        else:
            ifile = get_diffuse_source('ConstantValue',None,'FileFunction',ifilex,'Isotropic Diffuse (%s)'%ifile)
            if limit_parameters: ifile.smodel.set_default_limits()

    return [gfile, ifile]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
