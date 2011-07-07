"""
Module implements New modules to read in Catalogs of sources.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_catalogs.py,v 1.3 2011/06/20 18:24:31 lande Exp $

author: Joshua Lande
"""
import os
import numpy as np
from os.path import expandvars, exists, join
from textwrap import dedent

from pyfits import open
from skymaps import SkyDir

from . pointspec_helpers import PointSourceCatalog,PointSource
from . SpatialModels import Disk,EllipticalDisk,Gaussian,EllipticalGaussian,GAUSSIAN_X68,SpatialMap
from uw.utilities import keyword_options
from . Models import PowerLaw,PowerLawFlux,LogParabola,ExpCutoff
from . roi_extended import ExtendedSource

class Catalog2FGL(PointSourceCatalog):
    """ Catalog format suitable only for 2FGL.

        This object is designed only to open the 2FGL catalog file
        gll_psc_v04.fit. It is not supporeted for use on any other
        catalog.

        Does extended sources, different spectral models.  To create
        the extended soures, this object requires either the LATEXTDIR
        environment variable to be set or the paramter latextdir to be
        passed into the object. """

    defaults = (
        ("latextdir",    None, "Directory containing the spatial model templates."),
        ("prune_radius", 0.10, "[deg] consider sources closer than this duplicates"),
        ("free_radius",     2, "[deg] sources within this distance have free spectral parameters"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,catalog,**kwargs):
        keyword_options.process(self, kwargs)

        if self.latextdir is None and \
                (not os.environ.has_key('LATEXTDIR') or not exists(expandvars(os.environ['LATEXTDIR']))):
                    raise Exception(dedent("""
                            Since environment variable $LATEXTDIR does 
                            not exist, the paramter latextdir must
                            be passed into this object."""))

        else:
            os.environ['LATEXTDIR']=expandvars(self.latextdir)

        self.catalog=catalog

        self.__make_sources__()

    def __make_sources__(self):
        """ Make a list of all the point and extended
            sources in the catalog. """

        self.__extended_models__()

        f = open(expandvars(self.catalog))
        colnames = [x.name for x in f[1].get_coldefs()]
        point = f['LAT_POINT_SOURCE_CATALOG'].data
        ras       = point.field('RAJ2000')
        decs      = point.field('DEJ2000')
        pens      = point.field('PIVOT_ENERGY')
        n0s       = point.field('FLUX_DENSITY')
        inds      = point.field('SPECTRAL_INDEX')
        names     = point.field('SOURCE_NAME')
        extended_source_names = point.field('EXTENDED_SOURCE_NAME')
        f1000s     = point.field('FLUX1000')
        cutoffs   = point.field('CUTOFF')
        betas     = point.field('BETA')
        stypes    = point.field('SPECTRUMTYPE')
        f.close()

        self.names = names = np.chararray.strip(names)

        # not sure why there is the naming inconsistency
        extended_source_names = extended_source_names.replace(' ','')

        dirs   = map(SkyDir,np.asarray(ras).astype(float),np.asarray(decs).astype(float))

        # Store both point and extended sources.
        self.sources = []

        # This is for the new 2FGL style catalogs
        for name,extended_source_name,stype,n0,ind,pen,cutoff,beta,f1000,dir in \
                zip(names,extended_source_names,stypes,n0s,inds,pens,cutoffs,betas,f1000s,dirs):

            if stype == 'PowerLaw':
                if not np.isinf(n0):
                    model=PowerLaw(p=[n0,ind],e0=pen)
                else:
                    # For some reason, the fixed extended sources don't have n0 set.
                    # So create the source from is f1000 value.
                    model=PowerLawFlux(p=[f1000,ind],emin=1e3,emax=1e5)
            elif stype == 'LogParabola':
                model=LogParabola(p=[n0,ind,beta,pen])
            elif stype == 'PLExpCutoff':
                # Following M. Kerr's suggestion, clip the index of pulsars which are too close to 0.
                if ind >= 0 and ind < 1e-10: ind=1e-10
                model=ExpCutoff(p=[n0,ind,cutoff],e0=pen)
            else:
                raise Exception("Unkown spectral model %s for source %s" % (stype,name))

            if extended_source_name:
                spatial_model=self.__get_spatial_model__(extended_source_name)
                self.sources.append(
                        ExtendedSource(name=extended_source_name,
                            model=model,
                            spatial_model=spatial_model
                        )
                )
            else:
                self.sources.append(PointSource(dir,name,model))

    def __extended_models__(self):
        """ Read in the extended source spatial models
            from the catalog. """

        f = open(expandvars(self.catalog))
        extended = f['EXTENDEDSOURCES'].data
        self.extended_names = np.chararray.strip(extended.field('Source_Name'))
        self.extended_nicknames = self.extended_names.replace(' ','')

        ras   = extended.field('RAJ2000')
        decs  = extended.field('DEJ2000')
        forms = extended.field('Model_Form')
        majors = extended.field('Model_SemiMajor')
        minors = extended.field('Model_SemiMinor')
        posangs = extended.field('Model_PosAng')

        # The xml filename for the extended sources.
        templates = extended.field('Spatial_Filename').astype(str)
        f.close()

        # SpatialModel object for each extended source.
        self.extended_models = []

        for name,form,ra,dec,major,minor,posang,template in \
                zip(self.extended_nicknames,forms,ras,decs,majors,minors,posangs,templates):
            center=SkyDir(float(ra),float(dec))
            if form == 'Disk':
                if major == minor and posang == 0:
                    self.extended_models.append(Disk(p=[major],center=center))
                else:
                    self.extended_models.append(EllipticalDisk(p=[major,minor,posang],center=center))
            elif form == '2D Gaussian':
                if major == minor and posang == 0:
                    self.extended_models.append(Gaussian(p=[major/GAUSSIAN_X68],center=center))
                else:
                    self.extended_models.append(
                        EllipticalGaussian(p=[major/GAUSSIAN_X68,minor/GAUSSIAN_X68,posang],center=center))
            else:
                self.extended_models.append(SpatialMap(file=os.path.join('$LATEXTDIR','Templates',template)))

            # remember the fits file template in case the XML needs to be saved out.
            # (for gtlike compatability)
            self.extended_models[-1].original_template = template
            self.extended_models[-1].original_parameters = self.extended_models[-1].p.copy()

        self.extended_models = np.asarray(self.extended_models)

    def __get_spatial_model__(self,name):
        """ Return the spatial model corresponding to
            the source named 'name' """
        if sum(name==self.extended_nicknames) != 1:
            raise Exception("Cannot find spatial model for extended source %s" % name)

        return self.extended_models[self.extended_nicknames==name][0]

    @staticmethod
    def sort(list,skydir):
        list.sort(key=lambda src:src.skydir.difference(skydir))

    def get_source(self,name):
        return next(source for source in self.sources if source.name == name)

    def get_sources(self,skydir,radius=15):
        """ Returns all sources (point + diffuse combined) within radius.
            Sources only allowed to vary (their spectral paramters)
            when distance from skydir is less than free_radius. """
        return_sources = []

        for source in self.sources:
            distance=np.degrees(source.skydir.difference(skydir))
            if distance > radius:
                continue

            return_sources.append(source.copy())
            return_sources[-1].model.free[:] = (distance <= self.free_radius)

        Catalog2FGL.sort(return_sources,skydir)

        return return_sources


    def merge_lists(self,skydir,radius=15,user_point_list=[],user_diffuse_list=[]):
        """ Get a list of catalog sources and merge it with an (optional)
            list of PointSource objects provided by the user.  In case
            of duplicates (closer than prune_radius), the user object
            takes precedence. """

        user_source_list = [i for i in user_point_list + user_diffuse_list if hasattr(i,'skydir')]
        user_background_list = list(set(user_point_list + user_diffuse_list).difference(user_source_list))

        catalog_sources = self.get_sources(skydir,radius)

        merged_sources = []

        for cat in catalog_sources:
            merged_sources.append(cat)

            for user in user_source_list:
                if np.degrees(user.skydir.difference(cat.skydir)) < self.prune_radius or \
                        user.name == cat.name:
                    merged_sources.pop()
                    break

        user_point_list = [i for i in merged_sources if isinstance(i,PointSource)]
        user_extended_list = list(set(merged_sources).difference(user_point_list))

        Catalog2FGL.sort(user_point_list,skydir)
        Catalog2FGL.sort(user_extended_list,skydir)

        return user_point_list,user_extended_list+user_background_list
