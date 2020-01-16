""" A module to provide simple and standard access to pointlike fitting and spectral analysis.  The
    relevant parameters are fully described in the docstring of the constructor of the SpectralAnalysis
    class.

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pointspec.py,v 1.47 2012/10/04 21:21:22 lande Exp $

    author: Matthew Kerr
"""
version='$Revision: 1.47 $'.split()[1]
import types
import os
from os.path import join
import sys
from textwrap import dedent
import collections
from glob import glob
from datetime import date,timedelta

from . pixeldata import PixelData
from . pypsf     import Psf,CALDBPsf
from . pycaldb   import CALDBManager
from . pointspec_helpers import get_default_diffuse_mapper,PointlikeException,PointSource, ExposureManager
from . roi_catalogs import PointSourceCatalog
from . roi_managers import ROIPointSourceManager,ROIBackgroundManager,ROIDiffuseManager
from . roi_analysis import ROIAnalysis
from . roi_extended import ExtendedSource,ROIExtendedModel
from . roi_diffuse import DiffuseSource

from uw.utilities.fitstools import merge_bpd,merge_lt
from uw.utilities.fermitime import MET,utc_to_met
from uw.utilities.utils import get_data
from uw.utilities import keyword_options
from uw.utilities import path
from uw.utilities.xml_parsers import parse_sources

class DataSpecification(object):
    """ Specify the data to use for an analysis."""

    defaults= (
        ('ft1files',None,"""String or list of strings
         if a single string: points to a single FT1 file, or if
         contains wild cards is expanded by glob into a list of files
         if a list of string: each string points to an FT1 file"""),
        ('ft2files',None,"""String or list of string
          same format as ft1files, but N.B. that the current
          implementation expects a one-to-one correspondence of FT1"""),
        ('ltcube',None,"""string, optional
          points to a livetime cube; if not given, the livetime will
          be generated on-the-fly.  If a file is provided but does
          not exist, the livetime cube will be generated and written
          to the specified file."""),
        ('binfile',None,"""string, optional
          points to a binned representation of the data; will be
          generated if not provided; if file specified but does not
          exist, the binned data will be written to the file
          N.B. -- this file should be re-generated if, e.g., the
          energy binning used in the later spectral analysis changes."""),
    )
    
    @keyword_options.decorate(defaults)
    def __init__(self,**kwargs):

        keyword_options.process(self, kwargs)

        if self.binfile is None: raise PointlikeException("binfile must be specified.")
        if self.ltcube is None: raise PointlikeException("ltcube must be specified.")
            
        self.binfile = path.expand(self.binfile)
        self.ltcube = path.expand(self.ltcube)

        if self.ft1files is None and not os.path.exists(self.binfile):
            raise PointlikeException('An FT1 file must be specified if the binfile does not exist.')

        ltfile_exists = self.ltcube is not None and os.path.exists(self.ltcube)
        if self.ft2files is None and not ltfile_exists:
            raise PointlikeException('No FT2 or livetime file provided! Must pass at least one of these.')

        # If string, expand it out
        if isinstance(self.ft1files,types.StringType):
            self.ft1files = [path.expand(self.ft1files)]
        elif self.ft1files is None:
            pass
        else:
            self.ft1files = map(path.expand,self.ft1files)

        if isinstance(self.ft2files,types.StringType):
            self.ft2files = [path.expand(self.ft2files)]
        elif self.ft2files is None:
            pass
        else:
            self.ft2files = map(path.expand,self.ft2files)


class SpectralAnalysis(object):
    """ Interface to the spectral analysis code."""

    defaults = (
        'keywords controlling data binning and livetime calculation',
        ('roi_dir',None,'aperture center; if None, assume all-sky analysis'),
        ('exp_radius',20,'radius (deg) to use if calculate exposure or a ROI. (180 for full sky)'),
        ('zenithcut',105,'Maximum spacecraft pointing angle with respect to zenith to allow'),
        ('thetacut',66.4,'Cut on photon incidence angle'),
        ('event_class',3,'select class level (3 - diffuse; 2 - source; 1 - transient; 0 - Monte Carlo)'),
        ('conv_type',-1,'select conversion type (0 - front; 1 - back; -1 = front + back)'),
        ('tstart',0,'Default no cut on time; otherwise, cut on MET > tstart'),
        ('tstop',0,'Default no cut on time; otherwise, cut on MET < tstop'),
        ('recalcgti',False,'if True, try to get GTI from GT1 files; otherwise, try from livetime cube or binned data file'),
        ('binsperdec',4,'energy binning granularity when binning FT1'),
        ('emin',100,'Minimum energy'),
        ('emax',1e6,'Maximum energy'),
        ('use_weighted_livetime',False,'Use the weighted livetime'),
        'keywords for monte carlo data',
        ('mc_src_id',-1,'set to select on MC_SRC_ID column in FT1'),
        ('mc_energy',False,'set True to use MC_ENERGY instead of ENERGY'),
        'keywords controlling instrument response',
        ('irf',None,'Which IRF to use'),
        ('psf_irf',None,'specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!'),
        ('CALDB',None,'override the CALDB specified by $CALDB.'),
        ('custom_irf_dir',None,'override the CUSTOM_IRF_DIR specified by the env. variable'),
        ('keywords controlling spectral analysis'),
        ('background','1FGL','a choice of global model specifying a diffuse background; see ConsistentBackground for options'),
        ('maxROI',10,'maximum ROI for analysis; note ROI aperture is energy-dependent = max(maxROI,r95(e,conversion_type))'),
        ('minROI',5,dedent("""\
                minimum ROI analysis 
                ***N.B.***
                ROI radius is energy-dependent: minROI < r95(e,conversion_type) < maxROI
                That is, the radius is given by the 95% PSF containment for the band, subject
                to the constraint that it be >= minROI and <= maxROI
                **********""")),
        'miscelleneous keywords',
        ('quiet',False,'Set True to suppress (some) output'),
        ('verbose',False,'More output'),
        ('daily_data_path','/phys/groups/tev/scratch1/users/Fermi/data/daily','Data Path'),
    )
    @keyword_options.decorate(defaults)
    def __init__(self, data_specification, **kwargs):
        """
        Create a new spectral analysis object.

        data_specification: an instance of DataSpecification with links to the FT1/FT2,
                            and/or binned data / livetime cube needed for analysis
                            (see docstring for that class) """

        self.ae = self.dataspec = data_specification

        self.__dict__.update(self.dataspec.__dict__)
        keyword_options.process(self, kwargs)

        self.CALDBManager = CALDBManager(irf=self.irf,psf_irf=self.psf_irf,
            CALDB=self.CALDB,custom_irf_dir=self.custom_irf_dir)

         #TODO -- sanity check that BinnedPhotonData agrees with analysis parameters
        self.pixeldata = PixelData(self.__dict__)
        self.exposure  = ExposureManager(self)
        self.psf = CALDBPsf(self.CALDBManager)

    def set_psf_weights(self,skydir):
        """ Set the PSF to a new position.  Weights by livetime."""

        self.psf.set_weights(self.ltcube,skydir)

    def get_sources(self, roi_dir, point_sources=[], diffuse_sources=[],
                    xmlfile=None, diffdir=None,
                    catalogs = [],
                    include_radius = None):
        """ Create a list of PointSource and ROIDiffuseModel objects
            that are needed for the ROI. """

        for ps in point_sources: 
            if not isinstance(ps, PointSource):
                raise PointlikeException("Source %s is not a point source" % ps.name)
        if diffuse_sources is not None:
            for ds in diffuse_sources:
                if not isinstance(ds, DiffuseSource):
                    raise PointlikeException("Source %s is not a diffuse source" % ds.name)

        if xmlfile is not None:
            ps,ds = parse_sources(xmlfile,diffdir=diffdir, roi_dir=roi_dir,
                                  max_roi=self.maxROI+5 if include_radius is None else include_radius)
            point_sources += ps; diffuse_sources += ds

        if not isinstance(catalogs,collections.Iterable) or \
                isinstance(catalogs,types.StringType): catalogs = [catalogs]
        for cat in catalogs:
            if not isinstance(cat,PointSourceCatalog):
                raise PointlikeException("Catalog %s must be of type PointSourceCatalog" % cat)

            point_sources,diffuse_sources = cat.merge_lists(roi_dir,
                    self.maxROI+5 if include_radius is None else include_radius,
                    point_sources,diffuse_sources)

        if point_sources == [] and not self.quiet:
            print ('WARNING!  No point sources are included in the model.')
        if diffuse_sources == [] and not self.quiet:
            print ('WARNING!  No diffuse sources are included in the model.')

        return point_sources, diffuse_sources


    def roi(self, roi_dir = None,
            point_sources = [], diffuse_sources = [],
            xmlfile=None, diffdir=None,
            catalogs = [], include_radius = None,
            **kwargs):
        """
        return an ROIAnalysis object

        Arguments:

        roi_dir          [None] A SkyDir giving the center of the ROI.  If the user
                         does not specify one, the system tries to infer it: if the
                         user has provided a list of point_sources, the position of
                         the first is chosen.  Otherwise, it defaults to the roi_dir
                         member of this object; if this element is not set, an
                         exception is raised.

        point_sources    [[]] a list of PointSource objects to merge with a Catalog list

        xmlfile          The full path to a gtlike-style XML file giving the
                         source model for the ROI.

        diffdir          [None] An optional path if the files necessary for 
                         diffuse sources are specified relative to this path.

        catalogs         [[]] a list of PointSourceCatalog objects or strings;
                         ***
                         Catalogs  are processed sequentially;
                         if duplicate objects are found (as defined within
                         the PointSourceCatalog instances), then sources found in the
                         point_sources kwarg take precendence, followed by sources in the
                         first element of the list, and so forth

                         *********** NOTA BENE******************
                         The user must set at least point_sources or catalogs
                         if he/she wants point sources in the ROI.  There isn't
                         really isn't a sensible, machine-independent default.


                         If the user provides no point sources, the user must
                         provide a center for the ROI via the roi_dir kwarg or
                         by setting the roi_dir element of this object.
                         ***************************************

        diffuse_sources  [[]] a list of DiffuseSources; if an empty list, the
                         system tries to assemble a sensible default using the GLAST_EXT
                         environment variable. If None, no diffuse sources are used
                         (this is probably only useful for MC testing).

        include_radius [None] The radius within which all catalog sources and XML sources
                         are included. The default is 5 degrees larger than maxROI

        Optional Keyword Arguments:
            ==========   =============
            keyword      description
            ==========   =============
            fit_emin     [125,125] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            conv_type    [-1] conversion type selection (0 for front, 1 for back, -1 for all)
            ==========   =============
        """

        # process kwargs
        if not kwargs.has_key('quiet'): kwargs['quiet']=self.quiet 
        if not kwargs.has_key("conv_type"): kwargs['conv_type'] = self.conv_type
        # determine ROI center
        if roi_dir is None:
            roi_dir = self.roi_dir if hasattr(self,'roi_dir') else point_sources[0].skydir
        if roi_dir is None:
            raise PointlikeException('User must provide an ROI direction!  (See docstring.)')

        # set the PSF for the ROI center (important that this happens first)
        self.set_psf_weights(roi_dir)
        # instantiate and return ROIAnalysis object

        point_sources, diffuse_sources = self.get_sources(roi_dir=roi_dir,
                                                          point_sources=point_sources, 
                                                          diffuse_sources=diffuse_sources,
                                                          xmlfile=xmlfile, diffdir=diffdir, 
                                                          catalogs=catalogs, include_radius=include_radius)

        diffuse_mapper = get_default_diffuse_mapper(self,roi_dir,self.quiet)
        diffuse_models = [diffuse_mapper(ds) for ds in diffuse_sources]  if diffuse_sources is not None else []

        psm = ROIPointSourceManager(point_sources,roi_dir,quiet=self.quiet)
        dsm = ROIDiffuseManager(diffuse_models,roi_dir,quiet=self.quiet)
        return ROIAnalysis(roi_dir,psm,dsm,self,**kwargs)

    roi_from_xml = roi


    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        ignore = ('psf', 'exposure')
        for key in sorted(self.__dict__.keys()):
            if key in self.ae.__dict__ or key in ignore: continue # avoid duplication internal functions
            s += '\t%-20s: %s\n' %(key, self.__dict__[key])
        return s



