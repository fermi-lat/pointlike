"""  A module to provide simple and standard access to pointlike fitting and spectral analysis.  The
     relevant parameters are fully described in the docstring of the constructor of the SpectralAnalysis
     class.
    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/pointspec.py,v 1.45 2012/04/03 02:27:06 lande Exp $

    author: Matthew Kerr
"""
version='$Revision: 1.45 $'.split()[1]
import types
import os
from os.path import join
import sys
from textwrap import dedent
import collections
from glob import glob
from datetime import date,timedelta

from pixeldata import PixelData
from pypsf     import Psf,CALDBPsf
from pycaldb   import CALDBManager
from pointspec_helpers import *
from roi_managers import ROIPointSourceManager,ROIBackgroundManager,ROIDiffuseManager
from roi_analysis import ROIAnalysis
from roi_extended import ExtendedSource,ROIExtendedModel
from . roi_diffuse import DiffuseSource
from uw.utilities.fitstools import merge_bpd,merge_lt
from uw.utilities.fermitime import MET,utc_to_met
from uw.utilities.utils import get_data
from uw.utilities import keyword_options
from uw.utilities import path
import numpy as N

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

        if self.binfile is None: raise Exception("binfile must be specified.")
        if self.ltcube is None: raise Exception("ltcube must be specified.")
            
        self.binfile = path.expand(self.binfile)
        self.ltcube = path.expand(self.ltcube)

        if self.ft1files is None and not os.path.exists(self.binfile):
            raise Exception,'An FT1 file must be specified if the binfile does not exist.'

        ltfile_exists = self.ltcube is not None and os.path.exists(self.ltcube)
        if self.ft2files is None and not ltfile_exists:
            raise Exception,'No FT2 or livetime file provided! Must pass at least one of these.'

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


class SavedData(DataSpecification):
    """Specify saved daily data files to use for analysis.

    Current implementation assumes that the specified directory has subfolders
    named 'daily','weekly', and 'monthly', each with subfolders 'bpd' and 'lt',
    containing BinnedPhotonData and LivetimeCube files, respectively.  It looks
    for filenames of the format 'timescale_yyyymmdd_type.fits', where timescale
    is one of 'day', 'week', or 'month', and type is one of '#bpd' or 'lt' (the
    # in the bpd filename is the number of bins per decade).  The yyyymmdd is
    the date of the first day the file covers (UTC); in the monthly case, the
    dd is dropped. So, for example, the BinnedPhotonData for Nov 6, 2009, with
    4 bins/decade, is daily/bpd/day_20091106_4bpd.fits, the BinnedPhotonData
    for the week beginning May 10, 2010, with 8 bins/decade,  is
    weekly/bpd/week_20100510_8bpd.fits, and the livetime for the month of
    December 2008 is monthly/lt/month_200812_lt.fits.

    **Parameters**

    tstart: Starting time for the analysis in MET.  Note that if this is not
            the beginning of a day, the data for the full day containing tstart
            will be used.
    tstop: Ending time for the analysis in MET.  As with tstart, if it is not
           at a boundary between days, the data for the full day containing
           tstop will be used.
    """
    try:
        default_data_dir = os.environ['DATA_DIR']
    except KeyError:
        default_data_dir = ''

    new_defaults = (('data_dir',default_data_dir,'path to the saved data products')
               ,('use_weighted_livetime',False,'''Specify whether to get
                  the weighted livetimes.''')
               ,('binsperdec',4,'''Bins per decade for the
                  BinnedPhotonData files.''')
               )
    defaults = DataSpecification.defaults + new_defaults

    @keyword_options.decorate(defaults)
    def __init__(self,tstart,tstop,**kwargs):

        keyword_options.process(self,kwargs)

        if self.data_dir =='' or not os.path.exists(self.data_dir):
            raise Exception("""No valid data directory provided. Either the DATA_DIR environment
                               variable or the data_dir keyword argument must point to a valid 
                               directory.""")

        if not (os.path.exists(str(self.binfile)) and os.path.exists(str(self.ltcube))):
            tstart = self.tstart if self.tstart else utc_to_met(2008,8,4)
            tstop = self.tstop if self.tstop else utc_to_met(*date.today().timetuple()[:6])
            bpds,lts = get_data(tstart,tstop,data_dir=self.data_dir)
            start_date = MET(tstart).time
            stop_date = MET(tstop).time
            start_date = start_date.year*10000+start_date.month*100+start_date.day
            stop_date = stop_date.year*10000+stop_date.month*100+stop_date.day
            if not self.binfile:
                self.binfile = '%i-%i_%ibpd.fits'%(start_date,stop_date,self.binsperdec)
            if not self.ltcube:
                self.ltcube = '%i-%i_lt.fits'%(start_date,stop_date)
            if not os.path.exists(self.binfile):
                merge_bpd(bpds,self.binfile)
            if not os.path.exists(self.ltcube):
                merge_lt(lts,self.ltcube,weighted = self.use_weighted_livetime)

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

    def roi(self, roi_dir = None,
                  point_sources = [], catalogs = [], catalog_mapper = None,
                  diffuse_sources = [], diffuse_mapper = None,
                  catalog_include_radius = None,
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

        catalogs         [[]] a list of PointSourceCatalog objects or strings;
                         ***IF they are strings, they will be interpreted by the
                            catalog_mapper.  If this argument is not set, the
                            strings will be interpreted as Fermi-compatible catalogs.
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

        catalog_mapper   [None] a function or instance of a class (via __call__)
                         that takes a single argmument, a filename (including
                         path), and returns an object implementing the
                         PointSourceCatalog interface.

        diffuse_sources  [[]] a list of DiffuseSources; if an empty list, the
                         system tries to assemble a sensible default using the GLAST_EXT
                         environment variable. If None, no diffuse sources are used
                         (this is probably only useful for MC testing).

        diffuse_mapper   [None] a function or instance of a class (via __call__)
                         which takes two arguments, a DiffuseSource and the ROI
                         center, and returns an object implementing the
                         ROIDiffuseModel interface.  If None, the system uses
                         the default, an on-the-fly numerical convolution.
        catalog_include_radius [None] The radius within which all catalog sources are 
                         included. The default is 5 degrees larger than maxROI

        Optional Keyword Arguments:
            ==========   =============
            keyword      description
            ==========   =============
            fit_emin     [125,125] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            conv_type    [-1] conversion type selection (0 for front, 1 for back, -1 for all)
            diffdir      [None] a directory to look for the default diffuse models (e.g. gll_iem_v02.fit)
            ==========   =============
        """

        # process kwargs
        diffdir = None
        if kwargs.has_key('diffdir'): diffdir = kwargs.pop('diffdir')
        if not kwargs.has_key('quiet'): kwargs['quiet']=self.quiet 
        if not kwargs.has_key("conv_type"): kwargs['conv_type'] = self.conv_type
        # determine ROI center
        if roi_dir is None:
            roi_dir = self.roi_dir if hasattr(self,'roi_dir') else point_sources[0].skydir
        if roi_dir is None:
            raise Exception,'User must provide an ROI direction!  (See docstring.)'

        # set the PSF for the ROI center (important that this happens first)
        self.set_psf_weights(roi_dir)

        for ps in point_sources: 
            if not isinstance(ps, PointSource):
                raise Exception("Source %s is not a point source" % ps.name)
        if diffuse_sources is not None:
            for ds in diffuse_sources:
                if not isinstance(ds, DiffuseSource):
                    raise Exception("Source %s is not a diffuse source" % ds.name)

        # process point sources
        if catalog_mapper is None:
            catalog_mapper = lambda x: FermiCatalog(x)

        if not isinstance(catalogs,collections.Iterable) or \
                isinstance(catalogs,types.StringType): catalogs = [catalogs]
        for cat in catalogs:
            if not isinstance(cat,PointSourceCatalog):
                cat = catalog_mapper(cat)
            point_sources,diffuse_sources = cat.merge_lists(roi_dir,
                    self.maxROI+5 if catalog_include_radius is None else catalog_include_radius,
                    point_sources,diffuse_sources)
        if point_sources == [] and not self.quiet:
            print 'WARNING!  No point sources are included in the model.'

        # process diffuse models
        if diffuse_sources is not None and len(diffuse_sources) == 0:
            # try to use default
            diffuse_sources = get_default_diffuse(diffdir=diffdir)
            if len(diffuse_sources) == 0 and not self.quiet:
                print 'WARNING!  No diffuse sources are included in the model.'
        if diffuse_mapper is None:
            diffuse_mapper = get_default_diffuse_mapper(self,roi_dir,self.quiet)

        diffuse_models = [diffuse_mapper(ds) for ds in diffuse_sources]  if diffuse_sources is not None else []

        # instantiate and return ROIAnalysis object
        psm = ROIPointSourceManager(point_sources,roi_dir,quiet=self.quiet)
        dsm = ROIDiffuseManager(diffuse_models,roi_dir,quiet=self.quiet)
        return ROIAnalysis(roi_dir,psm,dsm,self,**kwargs)

    def roi_from_xml(self,roi_dir,xmlfile,diffdir=None,*args,**kwargs):
        """
        return an ROIAnalysis object with a source model specified by
        a gtlike-style XML file.

        Arguments:

        roi_dir          A SkyDir giving the center of the ROI.

        xmlfile          The full path to a gtlike-style XML file giving the
                         source model for the ROI.

        diffdir          [None] An optional path if the files necessary for 
                         diffuse sources are specified relative to this path.

        Optional Arguments and Keyword Arguments:
            see docstring for SpectralAnalysis.roi
        """
        from uw.utilities.xml_parsers import parse_sources
        ps,ds = parse_sources(xmlfile,diffdir=diffdir,roi_dir=roi_dir,max_roi=self.maxROI+5)
        return self.roi(roi_dir=roi_dir,point_sources=ps,diffuse_sources=ds,
                        *args,**kwargs)


    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        ignore = ('psf', 'exposure')
        for key in sorted(self.__dict__.keys()):
            if key in self.ae.__dict__ or key in ignore: continue # avoid duplication internal functions
            s += '\t%-20s: %s\n' %(key, self.__dict__[key])
        return s



