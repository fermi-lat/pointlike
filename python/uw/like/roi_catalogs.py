"""
Module implements New modules to read in Catalogs of sources.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_catalogs.py,v 1.30 2014/07/15 02:16:17 mdwood Exp $

author: Joshua Lande
"""
import os
from copy import copy
import glob
import numpy as np
from os.path import exists, join
from textwrap import dedent
from abc import abstractmethod

import pyfits
from skymaps import SkyDir

from . SpatialModels import Disk,EllipticalDisk,Gaussian,EllipticalGaussian,SpatialMap
from . Models import PowerLaw,PowerLawFlux,LogParabola,ExpCutoff,PLSuperExpCutoff
from . roi_extended import ExtendedSource
from uw.utilities import keyword_options
from uw.utilities.parmap import LogMapper,LimitMapper
from uw.utilities import path



class SourceCatalog(object):
    """ Define an interface for point source catalog that can be used to
         construct the model for an ROI analysis."""

    @abstractmethod
    def get_sources(self,skydir,radius=15):
        pass

    @abstractmethod
    def merge_lists(self,skydir,radius=15,user_point_list=None,user_diffuse_list=None):
        pass

PointSourceCatalog=SourceCatalog # For Backwards Compatability

class FermiCatalog(SourceCatalog):
    """Read a Jean Ballet catalogue and use it to set up a source list for a ROI.
    
        This class works (mostly) for EMS and 1FGL styles."""

    defaults = (
        ("prune_radius",0.10,"deg; in a merge, consider sources closer than this duplicates"),
        ("free_radius",2,"deg; sources within this distance have free spectral parameters"),
        ("min_ts",None,"If specified, only include sources with a larger catalog TS"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,catalog_file,**kwargs):
        keyword_options.process(self, kwargs)

        self.__open_catalog__(catalog_file)
        self.catalog_file=catalog_file

    def __open_catalog__(self,catalog_file):
        import pyfits
        f = pyfits.open(catalog_file)
        colnames = [x.name for x in f[1].get_coldefs()]
        sname     = 'NickName' if 'NickName' in colnames else 'Source_Name'
        ras  = f[1].data.field('RA' if 'RA' in colnames else 'RAJ2000')
        decs = f[1].data.field('DEC' if 'DEC' in colnames else 'DEJ2000')
        pens = f[1].data.field('PIVOT_ENERGY')
        n0s  = f[1].data.field('FLUX_DENSITY')
        inds = f[1].data.field('SPECTRAL_INDEX')
        cutoffs = f[1].data.field('CUTOFF_ENERGY') if 'Cutoff_Energy' in colnames else None
        betas = f[1].data.field('beta') if 'beta' in colnames else None       
        if "Ignore" in f[1].columns.names:
            mask = f[1].data.field("Ignore") == False
        else:
            mask = np.ones(ras.shape,dtype='bool')
        
        try:
            self.ts=np.asarray(f[1].data.field('TEST_STATISTIC'))
        except KeyError, er:
            if self.min_ts is not None:
                raise Exception("Cannot apply min_ts cut since TEST_STATISTIC column not found in fits file.")

        inds = np.where(inds > 0, inds, -inds)


        self.models = []
        for i,(n0,ind,pen) in enumerate(zip(n0s,inds,pens)):
            if ind==0:ind=1.e-5
            if cutoffs is not None and not np.isnan(cutoffs[i]) and not np.isinf(cutoffs[i]) :
                cutoff=cutoffs[i]
                self.models.append(ExpCutoff(p=[n0,ind,cutoff],e0=pen))
            elif betas is not None and not np.isnan(betas[i]) and not np.isinf(betas[i]):
                beta=betas[i]
                self.models.append(LogParabola(p=[n0,ind,beta,pen]))
            else:
                try:
                    self.models.append(PowerLaw(p=[n0,ind],e0=pen))
                except AssertionError:
                    print "Something went wrong with %s, discarding"%f[1].data.field(sname)[i]
                    mask[i] = False

        self.dirs    = map(SkyDir,np.asarray(ras).astype(float)[mask],np.asarray(decs).astype(float)[mask])
        self.models = np.asarray(self.models)
        self.names  = np.chararray.strip(f[1].data.field(sname))[mask]

        f.close()

    def get_sources(self,skydir,radius=15):

        from . pointspec_helpers import PointSource # Avoid circular imports

        diffs    = np.degrees(np.asarray([skydir.difference(d) for d in self.dirs]))
        mask     = diffs < radius
        if self.min_ts is not None: mask &= self.ts>self.min_ts

        diffs    = diffs[mask]
        sorting = np.argsort(diffs)

        #sort SkyDirs -- pesky vector behavior...
        dirs     = [x for i,x in enumerate(self.dirs) if mask[i]]
        dirs     = [dirs[x] for x in sorting]

        names    = self.names[mask][sorting]
        models  = [x.copy() for x in self.models[mask][sorting]]

        fm = diffs[sorting] < self.free_radius
        point_sources = map(PointSource,dirs,names,models,fm)
        return point_sources


    def merge_lists(self,skydir,radius=15,user_point_list=[],user_diffuse_list=[]):
        """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
            provided by the user.  In case of duplicates (closer than prune_radius), the user object
            takes precedence."""

        user_extended_list = [i for i in user_diffuse_list if isinstance(i,ExtendedSource)] \
                if user_diffuse_list is not None else []

        cat_list = self.get_sources(skydir,radius)
        if user_point_list==[] and user_extended_list==[]: return cat_list,user_diffuse_list

        from collections import deque
        merged_list = deque(user_point_list)


        for ncps,cps in enumerate(cat_list):
            merged_list.append(cps)
            for ups in user_point_list:
                if np.degrees(ups.skydir.difference(cps.skydir)) < self.prune_radius or \
                        ups.name == cps.name:
                    merged_list.pop(); break

            for ues in user_extended_list:
                if np.degrees(ues.spatial_model.center.difference(cps.skydir)) < self.prune_radius or \
                        ues.name == cps.name:
                    merged_list.pop(); break

        merged_list = list(merged_list)
        merged_list.sort(key = lambda ps: ps.skydir.difference(skydir))

        return merged_list,user_diffuse_list


class BaseCatalog2FGL(SourceCatalog):

    """ Catalog format suitable only for 2FGL.

        This object is designed only to open the 2FGL catalog file
        gll_psc_v04.fit. It is not supporeted for use on any other
        catalog.

        Does extended sources, different spectral models.  To create
        the extended soures, this object requires either the LATEXTDIR
        environment variable to be set or the paramter latextdir to be
        passed into the object. 

            >>> x=Catalog2FGL('/afs/slac/g/glast/groups/catalog/2FGL/gll_psc_v05.fit',
            ...               latextdir='/afs/slac/g/glast/groups/catalog/2FGL/gll_psc_v05_templates')
        
    """

    defaults = (
        ("latextdir",            None, "Directory containing the spatial model templates."),
        ("prune_radius",         0.10, "[deg] consider sources closer than this duplicates"),
        ("free_radius",             2, "[deg] sources within this distance have free spectral parameters"),
        ("max_free",             None, "Maximum number of sources to free (if there are more than that many sources within free_radius)."),
        ("limit_parameters",    False, "Limit fit parameters with default values"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,catalog,**kwargs):
        keyword_options.process(self, kwargs)

        print "Opening catalog from %s"%catalog

        if self.latextdir is None:
            if (not os.environ.has_key('LATEXTDIR') or
                    not exists(path.expand(os.environ['LATEXTDIR']))):
                    raise Exception(dedent("""
                            Since environment variable $LATEXTDIR does 
                            not exist, the paramter latextdir must
                            be passed into this object."""))

            else:
                self.latextdir=os.environ['LATEXTDIR']
        else:
            os.environ['LATEXTDIR']=path.expand(self.latextdir)

        self.catalog=catalog

        self.__make_sources__()

    def __make_sources__(self):
        """ Make a list of all the point and extended
            sources in the catalog. """

        from . pointspec_helpers import PointSource # Avoid circular imports

        self.__extended_models__()

        f = pyfits.open(path.expand(self.catalog))
        colnames = [x.name for x in f[1].get_coldefs()]

        if "Ignore" in colnames:
            sourceMask = f['LAT_POINT_SOURCE_CATALOG'].data.field("Ignore") == False
            point = f['LAT_POINT_SOURCE_CATALOG'].data[sourceMask]
        else:
            point = f['LAT_POINT_SOURCE_CATALOG'].data
        ras       = point.field('RAJ2000')
        decs      = point.field('DEJ2000')
        pens      = point.field('PIVOT_ENERGY')
        n0s       = point.field('FLUX_DENSITY')
        inds      = point.field('SPECTRAL_INDEX')
        f1000s     = point.field('FLUX1000')
        cutoffs   = point.field('CUTOFF')
        betas     = point.field('BETA')
        stypes    = point.field('SPECTRUMTYPE')
        f.close()

        # extended => is extended
        self.names, self.extendeds = self.kluge_names_extended(point)

        dirs = map(SkyDir,np.asarray(ras).astype(float),np.asarray(decs).astype(float))

        # Store both point and extended sources.
        self.sources = []

        # This is for the new 2FGL style catalogs
        for name,extended,stype,n0,ind,pen,cutoff,beta,f1000,dir in \
                zip(self.names,self.extendeds,stypes,n0s,inds,pens,cutoffs,betas,f1000s,dirs):

            if stype == 'PowerLaw':
                # note, np.isinf incorrectly raises annoying warning
                # http://projects.scipy.org/numpy/ticket/1500
                # This is a workaround.
                if n0 not in [-np.inf, np.inf]:
                    model=PowerLaw(norm=n0, index=ind, e0=pen)
                else:
                    # For some reason, in 2FGL the fixed extended sources don't 
                    # have n0 set. So create the source from the f1000 value.
                    model=PowerLawFlux(int_flux=f1000, index=ind, emin=1e3, emax=1e5)
            elif stype == 'PowerLaw2':
                model=PowerLawFlux(int_flux=f1000, index=ind, emin=1e3, emax=1e5)
            elif stype == 'LogParabola':
                model=LogParabola(norm=n0, index=ind, beta=beta, e_break=pen)
                model.freeze('e_break')
            elif stype == 'PLExpCutoff':
                # PLExpCutoff is defined by Eq. 1 of the 2FGL paper:
                # dN/dE = K*(E/E_pivot)**-gamma * exp(-(E - E_pivot)/E_cutoff)
                # This forces a correction factor to match pointlike ExpCutoff
                norm = n0 * np.exp(pen/cutoff)
                model=ExpCutoff(norm=norm, index=ind, cutoff=cutoff, e0=pen)
            elif stype == 'PLSuperExpCutoff':
                model=PLSuperExpCutoff(norm=n0, index=ind, cutoff=cutoff, e0=pen)
            else:
                raise Exception("Unkown spectral model %s for source %s" % (stype,name))
            
            if self.limit_parameters:
                model.set_default_limits(oomp_limits=True)

            if extended:
                spatial_model=self.__get_spatial_model__(name)
                self.sources.append(
                        ExtendedSource(name=name,
                            model=model,
                            spatial_model=spatial_model
                        )
                )
            else:
                self.sources.append(PointSource(dir,name,model,leave_parameters=True))

    def __extended_models__(self):
        """ Read in the extended source spatial models
            from the catalog. """

        f = pyfits.open(path.expand(self.catalog))

        if "Ignore" in f['EXTENDEDSOURCES'].columns.names:
            sourceMask = f['EXTENDEDSOURCES'].data.field("Ignore") == False
            extended = f['EXTENDEDSOURCES'].data[sourceMask]
        else:
            extended = f['EXTENDEDSOURCES'].data
        #self.extended_names = np.char.strip(extended.field('Source_Name'))
        #self.extended_nicknames = np.char.replace(self.extended_names,' ','')
        self.extended_names = np.asarray([x.strip() for x in extended.field('Source_Name')])
        self.extended_nicknames = np.asarray([x.replace(' ','') for x in self.extended_names])

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
                    self.extended_models.append(Gaussian(p=[major/Gaussian.x68],center=center))
                else:
                    self.extended_models.append(
                        EllipticalGaussian(p=[major/Gaussian.x68,minor/Gaussian.x68,posang],center=center))
            else:
                self.extended_models.append(SpatialMap(file=self.kluge_template_to_file_name(template)))

            # remember the fits file template in case the XML needs to be saved out.
            # (for gtlike compatability)
            self.extended_models[-1].original_template = self.kluge_template_to_file_name(template)
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

    def get_names(self):
        return self.names.tolist()

    def get_sources(self,skydir,radius=15):
        """ Returns all catalog sources (point + diffuse combined) within radius.
            Sources only allowed to vary (their spectral paramters)
            when distance from skydir is less than free_radius.
            
            A maximum of max_free sources are allowed to be specified. """
        return_sources = []

        for source in self.sources:
            distance=np.degrees(source.skydir.difference(skydir))
            if distance > radius:
                continue

            return_sources.append(source.copy())
            if distance > self.free_radius:
                return_sources[-1].model.free[:] = False

        Catalog2FGL.sort(return_sources,skydir)

        if self.max_free is not None:
            if sum(np.any(source.model.free==True) for source in return_sources) > self.max_free:
                i=0
                for source in return_sources:
                    if np.any(source.model.free==True):
                        i+=1
                        if i > self.max_free:
                            source.model.free[:] = False

        return return_sources


    def merge_lists(self,skydir,radius=15,user_point_list=[],user_diffuse_list=[]):
        """ Get a list of catalog sources and merge it with an (optional)
            list of PointSource objects provided by the user.  In case
            of duplicates (closer than prune_radius), the user object
            takes precedence. """

        from . pointspec_helpers import PointSource # Avoid circular imports

        user_source_list = [i for i in user_point_list + user_diffuse_list if hasattr(i,'skydir')]
        user_background_list = list(set(user_point_list + user_diffuse_list).difference(user_source_list))

        catalog_sources = self.get_sources(skydir,radius)

        merged_sources = copy(user_source_list)

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


class Catalog2FGL(BaseCatalog2FGL):
    """ The 2FGL and 3Y catalog are slightly different. """
    
    @staticmethod
    def kluge_names_extended(point):
        names     = point.field('SOURCE_NAME')
        names = np.asarray([x.strip() for x in names])
        extended_source_names = point.field('EXTENDED_SOURCE_NAME')
        # not sure why there is the naming inconsistency
        extended_source_names = np.asarray([x.strip().replace(' ','') for x in extended_source_names])

        extendeds = (extended_source_names != '')
        names[extendeds] = extended_source_names[extendeds]

        return names, extendeds

    @staticmethod
    def kluge_template_to_file_name(template):
        return os.path.join('$LATEXTDIR',template)
    
class Catalog3Y(BaseCatalog2FGL):

    @staticmethod
    def kluge_names_extended(point):
        names = point.field('NickName')
        names = names = np.asarray([x.strip() for x in names])
        extendeds = point.field('EXTENDED')
        return names, extendeds

    @staticmethod
    def kluge_template_to_file_name(template):

        # for some reason, Extended_Archive_v12 is inconsistent
        if 'RXJ1713_hess.fits' in template:
            return template.replace('RXJ1713_hess.fits','RXJ1713.fits')
        elif 'S147_Halpha.fits' in template:
            return template.replace('S147_Halpha.fits','S147.fits')
            pass
        else:
            return template

class ExtendedSourceCatalog(SourceCatalog):
    """ Object to read in spatially extended sources from Elizabeth Ferrara suggested
        table of extended sources.

        The input to this object is the directory contaning the extended source table
        and a folder for xml descriptions and templates of extended sources.

        For a description of the format of these sources, see:

            https://confluence.slac.stanford.edu/x/Qw2JBQ

        This code is, largely, depricated in fafor of the Catalog2FGL Object!
        
        But it *is* used by like2, which now requires map, for which there is a special kwarg to force use of a map
    '"""

    def __init__(self,archive_directory, force_map=False):
        self.force_map=force_map
        self.__open_catalog__(archive_directory)

    def __open_catalog__(self,archive_directory):
        """ Parses the LAT_extended_sources.fit table 
            to get a list of the extended sources. """
        self.archive_directory = archive_directory

        import pyfits
        filename=join(self.archive_directory,"LAT_extended_sources*.fit*")
        filename=glob.glob(filename)
        if len(filename)!=1: raise Exception("Unable to find LAT_extended_sources.fit archive file.")
        filename=filename[0]
        f = pyfits.open(filename)
        self.names = f[1].data.field('Source_Name')
        ras   = f[1].data.field('RAJ2000')
        decs  = f[1].data.field('DEJ2000')
        form = f[1].data.field('Model_Form')
        major = f[1].data.field('Model_SemiMajor')
        minor = f[1].data.field('Model_SemiMinor')
        posang = f[1].data.field('Model_PosAng')

        # The xml filename for the extended sources.
        self.xmls      = f[1].data.field('Spectral_Filename').astype(str)
        self.templates = f[1].data.field('Spatial_Filename').astype(str)

        self.dirs    = map(SkyDir,np.asarray(ras).astype(float),np.asarray(decs).astype(float))

        if self.archive_directory in [ "Extended_archive_v01", "Extended_archive_v02", 
                                       "Extended_archive_jbb02", "Extended_archive_jbb03" ]:

            # old style archives
            os.environ["LATEXTDIR"]=join(self.archive_directory,'Templates')
        else:
            os.environ["LATEXTDIR"]=self.archive_directory


        # build up a list of the analytic extended source shapes (when applicable)
        self.spatial_models = []
        for i in range(len(self.names)):
            if self.force_map:
                self.spatial_models.append(
                    SpatialMap(file=self.templates[i]))
            elif form[i] == 'Disk':
                if major[i] == minor[i] and posang[i] == 0:
                    self.spatial_models.append(Disk(p=[major[i]],center=self.dirs[i]))
                else:
                    self.spatial_models.append(EllipticalDisk(p=[major[i],minor[i],posang[i]],center=self.dirs[i]))
            elif form[i] == '2D Gaussian':
                if major[i] == minor[i] and posang[i] == 0:
                    self.spatial_models.append(Gaussian(p=[major[i]/Gaussian.x68],center=self.dirs[i]))
                else:
                    self.spatial_models.append(
                        EllipticalGaussian(p=[major[i]/Gaussian.x68,minor[i]/Gaussian.x68,posang[i]],
                                           center=self.dirs[i]))
            else:
                self.spatial_models.append(
                    SpatialMap(file=self.templates[i])
                )
            # remember the fits file template in case the XML needs to be saved out.
            self.spatial_models[-1].original_template = self.templates[i]
            self.spatial_models[-1].original_parameters = self.spatial_models[-1].p.copy()


        self.spatial_models = np.asarray(self.spatial_models)

    def get_sources(self,skydir,radius=15):
        """ Returns a list of ExtendedSource objects from the extended source
            catalog that have a center withing a distance radius of the
            position skydir. 
           
        Note that if there are none, it returns an empty list.   
        """

        from uw.utilities.xml_parsers import parse_sources

        diffs    = np.degrees(np.asarray([skydir.difference(d) for d in self.dirs]))
        mask     = diffs < radius
        if sum(mask)==0: return []
        diffs    = diffs[mask]
        sorting = np.argsort(diffs)

        names     = self.names[mask][sorting]
        xmls      = self.xmls[mask][sorting]
        spatials  = self.spatial_models[mask][sorting]

        sources = []
        for name,xml,spatial in zip(names,xmls,spatials):

            full_xml=join(self.archive_directory,'XML',os.path.basename(xml))

            # Use the built in xml parser to load the extended source.
            ps,ds=parse_sources(xmlfile=full_xml)
            if len(ps) > 0: 
                raise Exception("A point source was found in the extended source file %s" % xmlfile)
            if len(ds) > 1: 
                raise Exception("No diffuse sources were found in the extended soruce file %s" % xmlfile)
            if len(ds) < 1: 
                raise Exception("More than one diffuse source was found in the extended source file %s" % xmlfile)

            source=ds[0]
            if spatial is not None:
                # replace the SpatialMap extended source with an analytic one.
                analytic_source = ExtendedSource(name=source.name,model=source.model,
                                                 spatial_model=spatial)
                analytic_source.original_template = source.spatial_model.file # for reference
                sources.append(analytic_source)
            else:
                sources.append(source)

        return sources

    def merge_lists(self,skydir,radius=15,user_point_list=[],user_diffuse_list=[]):
        """ Get a list of extended sources from the  catalog and merge it with 
            and already existin glist of point and diffuse sources.

            Unlike the FermiCatalog, no source pruning is done. """

        cat_list = self.get_sources(skydir,radius)

        return user_point_list,user_diffuse_list+cat_list


CatalogManager=FermiCatalog # For backwards compatability.

if __name__ == "__main__":
    import doctest
    doctest.testmod()
