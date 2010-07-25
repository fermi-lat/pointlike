"""  A module to provide simple and standard access to pointlike fitting and spectral analysis.  The
     relevant parameters are fully described in the docstring of the constructor of the SpectralAnalysis
     class.
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pointspec.py,v 1.17 2010/07/06 23:01:05 lande Exp $

    author: Matthew Kerr
"""
version='$Revision: 1.17 $'.split()[1]
import os
from os.path import join
import sys
from glob import glob
from datetime import date,timedelta

from pixeldata import PixelData
from pypsf     import Psf,OldPsf,NewPsf,CALDBPsf
from pointspec_helpers import *
from roi_managers import ROIPointSourceManager,ROIBackgroundManager,ROIDiffuseManager
from roi_analysis import ROIAnalysis
from roi_diffuse import ROIDiffuseModel_OTF
from roi_extended import ExtendedSource,ROIExtendedModel
from uw.utilities.fitstools import merge_bpd,sum_ltcubes
from uw.utilities.fermitime import MET,utc_to_met
from uw.utilities.utils import get_data
import numpy as N

class DataSpecification(object):
    """ Specify the data to use for an analysis.

      **Parameters**

      ft1files : string or list of strings
          if a single string: points to a single FT1 file, or if
          contains wild cards is expanded by glob into a list of files
          if a list of string: each string points to an FT1 file
      ft2files : string or list of string
          same format as ft1files, but N.B. that the current
          implementation expects a one-to-one correspondence of FT1
          and FT2 files!
      ltcube : string, optional
          points to a livetime cube; if not given, the livetime will
          be generated on-the-fly.  If a file is provided but does
          not exist, the livetime cube will be generated and written
          to the specified file.
      binfile : string, optional
          points to a binned representation of the data; will be
          generated if not provided; if file specified but does not
          exist, the binned data will be written to the file
          N.B. -- this file should be re-generated if, e.g., the
          energy binning used in the later spectral analysis changes.
      weighted_ltcube : string, optional
          points to a livetime cube with integrals weighted by livetime
          fraction; if not given, but needed, the weighted livetime will
          be generated on-the-fly.  If a file is provided but does
          not exist, the livetime cube will be generated and written
          to the specified file.
    """


    def init(self):

        self.ft1files = None
        self.ft2files = None
        self.ltcube   = None
        self.binfile  = None
        self.weighted_ltcube = None

    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

        if self.ft1files is None and self.binfile is None:
            raise Exception,'No event data (FT1 or binfile) provided!  Must pass at least one of these.'

        ltfile_exists = self.ltcube is not None and os.path.exists(self.ltcube)
        if self.ft2files is None and (self.ltcube is None or (not ltfile_exists)):
            raise Exception,'No FT2 or livetime file provided! Must pass at least one of these.'

        # make sure everything is iterable or None
        ft1 = self.ft1files; ft2 = self.ft2files
        self.ft1files = ft1 if (hasattr(ft1,'__iter__') or ft1 is None) else [ft1]
        self.ft2files = ft2 if (hasattr(ft2,'__iter__') or ft2 is None) else [ft2]

class SavedData(DataSpecification):
    """Specify saved daily data files to use for analysis.

    Current implementation assumes that the specified directory has subfolders named 'daily','weekly', and 'monthly',
    each with subfolders 'bpd','lt', and 'lt_weighted', containing BinnedPhotonData, LivetimeCube, and weighted LivetimeCube
    files, respectively.  It looks for filenames of the format 'timescale_yyyymmdd_type.fits', where timescale is one of 'day',
    'week', or 'month', and type is one of '#bpd', 'lt', or 'lt_weighted' (the # in the bpd filename is the number of bins per 
    decade).  The yyyymmdd is the date of the first day the file covers (UTC); in the monthly case, the dd is dropped. So, for 
    example, the BinnedPhotonData for Nov 6, 2009, with 4 bins/decade, is daily/bpd/day_20091106_4bpd.fits, the weighted 
    livetime cube for the week beginning May 10, 2010 is weekly/lt_weighted/week_20100510_lt_weighted.fits, and the unweighted 
    livetime for the month of December 2008 is monthly/lt/month_200812_lt.fits.

    **Parameters**

    tstart: Starting time for the analysis in MET.  Note that if this is not the beginning of a day, the data for the full day
        containing tstart will be used.
    tstop: Ending time for the analysis in MET.  As with tstart, if it is not at a boundary between days, the data for the full day
        containing tstop will be used.
    data_dir: string['/phys/groups/tev/scratch1/users/Fermi/data']
        path to the saved data products
    use_weighted_livetime: boolean[False]
        Specify whether to get the weighted livetimes.
    binsperdec: int[4] Bins/decade for the BinnedPhotonData files. Currently, we only have 4bpd files saved, but this should provide flexibility
	if we save additional binnings in the future.
    ltcube: string, optional
        File under which to save the sum of the livetime cubes for the specified time range. If not specified, the file will still
        be saved, under a name derived from the time range. 
    weighted_ltcube: string, optional 
        File under which to save the sum of the weighted livetime cubes for the specified time range. If not specified, the file will still
        be saved, under a name derived from the time range. 
    binfile: string, optional
        File under which to save the sum of the BinnedPhotonDatas for the specified time range. If not specified, the file will still
        be saved, under a name derived from the time range. 
    
    """
    
    def __init__(self,tstart,tstop,**kwargs):
        self.init()
	self.tstart = tstart
        self.tstop = tstop
        self.binsperdec = 4
        self.use_weighted_livetime = False
        self.data_dir = '/phys/groups/tev/scratch1/users/Fermi/data'
        self.__dict__.update(kwargs)


        if not (os.path.exists(str(self.binfile)) and os.path.exists(str(self.ltcube)) and 
               (os.path.exists(str(self.weighted_ltcube)) or not self.use_weighted_livetime)):
            tstart = self.tstart if self.tstart else utc_to_met(2008,8,4)
            tstop = self.tstop if self.tstop else utc_to_met(*date.today().timetuple()[:6])
            bpds,lts,weighted_lts = get_data(tstart,tstop,data_dir=self.data_dir)
            start_date = MET(tstart).time
            stop_date = MET(tstop).time
            start_date = start_date.year*10000+start_date.month*100+start_date.day
            stop_date = stop_date.year*10000+stop_date.month*100+stop_date.day
            if not self.binfile:
                self.binfile = '%i-%i_%ibpd.fits'%(start_date,stop_date,self.binsperdec)
            if not self.ltcube:
                self.ltcube = '%i-%i_lt.fits'%(start_date,stop_date)
            if self.use_weighted_livetime and not self.weighted_ltcube:
                self.weighted_ltcube = '%i-%i_lt_weighted.fits'%(start_date,stop_date)
            merge_bpd(bpds,self.binfile)
            sum_ltcubes(lts,self.ltcube)
            if self.use_weighted_livetime:
                sum_ltcubes(weighted_lts,self.weighted_ltcube)

########################################
########### DEPRECATED #################
########################################

class AnalysisEnvironment(object):
   """A class to collect locations of files needed for analysis.

      Specify the a directory in which to find the diffuse models,
      a LAT source list, the relevant version of CALDB, and the data
      files (FT1, FT2, livetime) to be used in the analysis.

      While these parameters can be specified as keyword arguments,
      it is probably best to update the default values in
      AnalysisEnvironment.init.

      **Parameters**

      ft1files : string or list of strings
          if a single string: points to a single FT1 file, or if
          contains wild cards is expanded by glob into a list of files
          if a list of string: each string points to an FT1 file
      ft2files : string or list of string
          same format as ft1files, but N.B. that the current
          implementation expects a one-to-one correspondence of FT1
          and FT2 files!
      ltcube : string, optional
          points to a livetime cube; if not given, the livetime will
          be generated on-the-fly.  If a file is provided but does
          not exist, the livetime cube will be generated and written
          to the specified file.
      weighted_ltcube : string, optional
          points to a livetime-fraction-weighted livetime_cube; if not
          given, livetime will be generated on the fly if needed.  If 
          a filename is provided and does not exist, the livetime cube
          will be generated and written to the specified file, if needed.
      binfile : string, optional
          points to a binned representation of the data; will be
          generated if not provided; if file specified but does not
          exist, the binned data will be written to the file
          N.B. -- this file should be re-generated if, e.g., the
          energy binning used in the later spectral analysis changes.
      diffdir : string, optional (but default must be correct!)
          the directory in which to find mapcubes and tables for
          the diffuse backgrounds
      catdir : string, optional (but default must be correct!)
          the directory in which to find a FITS representation of
          a LAT source list
      CALDB : string, optional (but default must be correct!)
          the directory in which to find the relevant version of CALDB.
          see the default value for the precise subdirectory required.
   """

   def init(self):

      self.diffdir = r'f:/glast/data/galprop'
      self.catdir  = r'f:/glast/data/kerr'
      #self.CALDB   = r'f:/glast/caldb/v1r1/CALDB/data/glast/lat'
      self.CALDB   = r'd:/fermi/caldb/v1r1/CALDB/data/glast/lat'

      self.ft1files = None
      self.ft2files = None
      self.ltcube   = None
      self.weighted_ltcube = None
      self.binfile  = None

   def __init__(self,**kwargs):
      self.init()
      self.__dict__.update(kwargs)

      if self.ft1files is None and self.binfile is None:
         print 'No event data (FT1 or binfile) provided!  Must pass at least one of these.'
         raise Exception

      try:
         f = open(self.ltcube)
         ltfile_exists = True
      except (IOError, TypeError):
         ltfile_exists = False

      if self.ft2files is None and (self.ltcube is None or ltfile_exists == False):
         print 'No spacecraft history (FT2) or livetime file provided! Must pass at least one of these.'
         raise Exception

      # make everything consistently a list
      ft1files = self.ft1files; ft2files = self.ft2files
      self.ft1files = ft1files if type(ft1files)==type([]) or ft1files is None else [ft1files]
      self.ft2files = ft2files if type(ft2files)==type([]) or ft2files is None else [ft2files]



class SpectralAnalysis(object):
    """
    Interface to the spectral analysis code.
Create a new spectral analysis object.

    analysis_environment: an instance of AnalysisEnvironment correctly configured with
                          the location of files needed for spectral analysis (see its
                          docstring for more information.)

Optional keyword arguments:

  ===========    =======================================================
  parameter      comments
  ===========    =======================================================
  **binning and livetime calculation**
  ----------------------------------------------------------------------
  roi_dir        [ None] aperture center; if None, assume all-sky analysis
  exp_radius     [ 20]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)
  zenithcut      [ 105]  Maximum spacecraft pointing angle with respect to zenith to allow
  thetacut       [ 66.4] Cut on photon incidence angle
  event_class    [ 3]  select class level (3 - diffuse; 2 - source; 1 - transient; 0 - Monte Carlo)
  conv_type      [ -1] select conversion type (0 - front; 1 - back; -1 = front + back)
  tstart         [0] Default no cut on time; otherwise, cut on MET > tstart
  tstop          [0] Default no cut on time; otherwise, cut on MET < tstop
  binsperdec     [4] energy binning granularity when binning FT1
  emin           [100] Minimum energy
  emax           [3e5] Maximum energy
  **Monte carlo selection**
  ----------------------------------------------------------------------
  mc_src_id      [ -1] set to select on MC_SRC_ID column in FT1
  mc_energy      [False] set True to use MC_ENERGY instead of ENERGY
  **Instrument response**
  ----------------------------------------------------------------------
  irf            ['P6_v3_diff'] Which IRF to use
  psf_irf        [None] specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!
  **spectral analysis**
  ----------------------------------------------------------------------
  background     ['ems_ring'] - a choice of global model specifying a diffuse background; see ConsistentBackground for options
  maxROI         [25] maximum ROI for analysis; note ROI aperture is energy-dependent = max(maxROI,r95(e,conversion_type))
  minROI         [0] minimum ROI analysis

                 **Note** The ROI radius is energy-dependent: minROI < r95(e,conversion_type) < maxROI
                  That is, the radius is given by the 95% PSF containment for the band, subject
                  to the constraint that it be >= minROI and <= maxROI
  **Miscellaneous**
  ----------------------------------------------------------------------
  quiet          [False] Set True to suppress (some) output
  verbose        [False] More output
  ===========    =======================================================

    """

    def __init__(self, data_specification, **kwargs):
        """

Create a new spectral analysis object.

    data_specification: an instance of DataSpecification with links to the FT1/FT2,
                        and/or binned data / livetime cube needed for analysis
                        (see docstring for that class)

Optional keyword arguments:

  =========    =======================================================
  Keyword      Description
  =========    =======================================================

  =========    KEYWORDS CONTROLLING DATA BINNING AND LIVETIME CALCULATION
  roi_dir      [ None] aperture center; if None, assume all-sky analysis
  exp_radius   [ 20]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)
  zenithcut    [ 105]  Maximum spacecraft pointing angle with respect to zenith to allow
  thetacut     [ 66.4] Cut on photon incidence angle
  event_class  [ 3]  select class level (3 - diffuse; 2 - source; 1 - transient; 0 - Monte Carlo)
  conv_type    [ -1] select conversion type (0 - front; 1 - back; -1 = front + back)
  tstart       [0] Default no cut on time; otherwise, cut on MET > tstart
  tstop        [0] Default no cut on time; otherwise, cut on MET < tstop
  recalcgti    [False] if True, try to get GTI from GT1 files;
                        otherwise, try from livetime cube or binned data file
  binsperdec   [4] energy binning granularity when binning FT1
  emin         [100] Minimum energy
  emax         [3e5] Maximum energy
  use_weighted_livetime [False] Use the weighted livetime
  use_daily_data [False] For the local cluster, use the pre-binned daily data.
  daily_data_path ['/phys/groups/tev/scratch1/users/Fermi/data/daily'] 
                  If use_daily_data, path to look for daily files.
  ***N.B.***
  If use_daily_data is True, tstart and tstop will be respected modulo one day: the
  full days containing tstart and tstop will be used. Also, it is left to the user to 
  ensure that the irf used is compatible with that used to bin the data (currently only
  P6_v3_diff on the TeV cluster).

  use_daily_data DEPRECATED in favor of using the SavedData class.
  **********

  =========    KEYWORDS FOR MONTE CARLO DATA
  mc_src_id    [ -1] set to select on MC_SRC_ID column in FT1
  mc_energy    [False] set True to use MC_ENERGY instead of ENERGY

  =========    KEYWORDS CONTROLLING INSTRUMENT RESPONSE
  irf          ['P6_v3_diff'] Which IRF to use
  psf_irf      [None] specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!
  CALDB        [environment variable] override the CALDB specified by the env. variable

  =========    KEYWORDS CONTROLLING SPECTRAL ANALYSIS
  background   ['1FGL'] - a choice of global model specifying a diffuse background; see ConsistentBackground for options
  maxROI       [25] maximum ROI for analysis; note ROI aperture is energy-dependent = max(maxROI,r95(e,conversion_type))
  minROI       [0] minimum ROI analysis

  ***N.B.***
  ROI radius is energy-dependent: minROI < r95(e,conversion_type) < maxROI
  That is, the radius is given by the 95% PSF containment for the band, subject
  to the constraint that it be >= minROI and <= maxROI
  **********

  =========    MISCELLENEOUS KEYWORDS
  quiet        [False] Set True to suppress (some) output
  verbose      [False] More output
  =========   =======================================================
  """

        self.roi_dir     = None
        self.exp_radius  = 20     # deg
        self.zenithcut   = 105    # deg
        self.thetacut    = 66.4   # deg
        self.event_class = 3      # diffuse
        self.conv_type   = -1     # front + back
        self.binsperdec  = 4
        self.tstart      = 0
        self.tstop       = 0
        self.recalcgti   = False
        self.emin        = 200    # MeV  -- note changed defaults to better deal with energy dispersion
        self.emax        = 2e5    # MeV
        self.use_weighted_livetime = False

        self.mc_src_id   = -1
        self.mc_energy   = False

        self.irf         = 'P6_v3_diff'
        self.psf_irf     = None
        self.CALDB       = os.environ['CALDB']
        self._check_CALDB()

        self.background  = '1FGL'
        self.maxROI      = 10    # deg
        self.minROI      = 5     # deg

        self.quiet       = False
        self.verbose     = False

        self.use_daily_data = False  #DEPRECATED
        self.daily_data_path = '/phys/groups/tev/scratch1/users/Fermi/data/daily'

        self.ae = self.dataspec = data_specification

        self.__dict__.update(self.dataspec.__dict__)
        self.__dict__.update(**kwargs)

        if self.use_daily_data:
            if not ((self.binfile and self.ltcube) and
                    (os.path.exists(self.binfile) and os.path.exists(self.ltcube))):
                self.setup_daily_data()

         #TODO -- sanity check that BinnedPhotonData agrees with analysis parameters
        self.pixeldata = PixelData(self.__dict__)
        self.exposure  = ExposureManager(self)
        self.psf = CALDBPsf(self.CALDB,irf=self.irf,psf_irf=self.psf_irf)

    def _check_CALDB(self):
        """ Make sure the CALDB member includes the entire path."""
        # this implementation a little crude, but there appear to be only
        # two possibilities
        toks = os.path.split(self.CALDB)
        if toks[1] != 'lat':
            self.CALDB = os.path.join(self.CALDB,'data','glast','lat')


    def setup_daily_data(self):
        """Setup paths to use saved daily data and livetime files on our local cluster.
        
        Should no longer be necossary, but keeping around for now for backward compatability."""

        data_dir = self.daily_data_path
        bpds = glob(os.path.join(data_dir,'bpd','*.fits'))
        if self.weighted_livetime:
            lts = glob(os.path.join(data_dir,'lt_weighted','*.fits'))
        else:
            lts = glob(os.path.join(data_dir,'lt','*.fits'))
        bpds.sort()
        lts.sort()
        start_date = MET(self.tstart).time if self.tstart else date(2008,8,4)
        stop_date = MET(self.tstop).time if self.tstop else date.today()-timedelta(1,0,0)
        start_date = start_date.year*10000+start_date.month*100+start_date.day
        stop_date = stop_date.year*10000+stop_date.month*100+stop_date.day
        start_ind = bpds.index(os.path.join(data_dir,'bpd','%i_%ibpd.fits'%(start_date,self.binsperdec)))
        stop_ind = bpds.index(os.path.join(data_dir,'bpd','%i_%ibpd.fits'%(stop_date,self.binsperdec)))
        bpds = bpds[start_ind:stop_ind+1]
        lts = lts[start_ind:stop_ind+1]
        if not self.binfile:
            self.binfile = '%i-%i_%ibpd.fits'%(start_date,stop_date,self.binsperdec)
        if not self.ltcube:
            if self.weighted_livetime:
                self.ltcube = '%i-%i_lt_weighted.fits'%(start_date,stop_date)
            else:
                self.ltcube = '%i-%i_lt.fits'%(start_date,stop_date)
        merge_bpd(bpds,self.binfile)
        sum_ltcubes(lts,self.ltcube)

    def set_psf_weights(self,skydir):
        """ Set the PSF to a new position.  Weights by livetime."""
        
        self.psf.set_weights(self.ltcube,skydir)

    def roi(self, roi_dir = None,
                  point_sources = [], catalogs = [], catalog_mapper = None,
                  diffuse_sources = [], diffuse_mapper = None,
                  *args,**kwargs):
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

        diffuse_sources  [[]] a list of DiffuseSources; if None, the
                         system tries to assemble a sensible default using the GLAST_EXT
                         environment variable.

        diffuse_mapper   [None] a function or instance of a class (via __call__)
                         which takes two arguments, a DiffuseSource and the ROI
                         center, and returns an object implementing the
                         ROIDiffuseModel interface.  If None, the system uses
                         the default, an on-the-fly numerical convolution.
         
        Optional Keyword Arguments:
            ==========   =============
            keyword      description
            ==========   =============
            fit_emin     [100,100] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            diffdir      [None] a directory to look for the default diffuse models (e.g. gll_iem_v02.fit)
            ==========   =============
        """

        # process kwargs
        diffdir = None
        if 'diffdir' in kwargs.keys(): diffdir = kwargs.pop('diffdir')
        
        # determine ROI center
        if roi_dir is None:
            roi_dir = self.roi_dir if len(point_sources)==0 else point_sources[0].skydir             
        if roi_dir is None:
            raise Exception,'User must provide an ROI direction!  (See docstring.)'

        # set the PSF for the ROI center (important that this happens first)
        self.set_psf_weights(roi_dir)

        # process point sources
        if catalog_mapper is None:
            catalog_mapper = lambda x: FermiCatalog(x)
        for cat in catalogs:
            if not isinstance(cat,PointSourceCatalog):
                cat = catalog_mapper(cat)
            point_sources,diffuse_sources = cat.merge_lists(roi_dir,self.maxROI+5,point_sources,diffuse_sources)
        if point_sources == []:
            print 'WARNING!  No point sources are included in the model.'

        # process diffuse models
        if len(diffuse_sources) == 0:
            # try to use default
            diffuse_sources = get_default_diffuse(diffdir=diffdir)
            if len(diffuse_sources) == 0:
                print 'WARNING!  No diffuse sources are included in the model.'
        if diffuse_mapper is None:
            diffuse_mapper = lambda x: ROIExtendedModel.factory(self,x,roi_dir) \
                    if isinstance(x,ExtendedSource) \
                    else ROIDiffuseModel_OTF(self,x,roi_dir)


        diffuse_models = [diffuse_mapper(ds) for ds in diffuse_sources]

        # instantiate and return ROIAnalysis object
        psm = ROIPointSourceManager(point_sources,roi_dir,quiet=self.quiet)
        dsm = ROIDiffuseManager(diffuse_models,roi_dir,quiet=self.quiet)           
        return ROIAnalysis(roi_dir,psm,dsm,self,**kwargs)

    def roi_from_xml(self,roi_dir,xmlfile,diffuse_mapper=None,
                          diffdir=None,*args,**kwargs):
        """
        return an ROIAnalysis object with a source model specified by
        a gtlike-style XML file.

        Arguments:
            
        roi_dir          A SkyDir giving the center of the ROI.

        xmlfile          The full path to a gtlike-style XML file giving the
                         source model for the ROI.

        diffuse_mapper   [None] a function or functor 
                         which takes two arguments, a DiffuseSource and the ROI
                         center, and returns an object implementing the
                         ROIDiffuseModel interface.  If None, the system uses
                         the default, an on-the-fly numerical convolution.

        diffdir          [None] An optional path if the files necessary for 
                         diffuse sources are specified relative to this path.

        Optional Arguments and Keyword Arguments:
            see docstring for SpectralAnalysis.roi
        """
        from uw.utilities.xml_parsers import parse_sources
        ps,ds = parse_sources(xmlfile,diffdir=diffdir)
        return self.roi(roi_dir=roi_dir,point_sources=ps,diffuse_sources=ds,
                        diffuse_mapper=diffuse_mapper,*args,**kwargs)
                

    def roi_old(self, point_sources = None, bgmodels = None, previous_fit = None, no_roi = False, **kwargs):
        """
        return an ROIAnalysis object with default settings.

        point_sources    [None] a list of PointSource objects to merge with a Catalog list
                         (if None, the nearest catalog source will be fit)

        bgmodels         a list of ROIBackgroundModels with which to override the default
                         isotropic and Galactic backgrounds (optional)

        previous_fit     [None] a file containing the results of an earlier spectral fit;
                         if not None, set spectral values to this fit
                         ***WARNING*** not tested!

        no_roi           [False] if True, return a tuple with an instance of PSManager and
                         of BGManager; can be used to instantiate an ROIAnalysis (or child)
                         object later


        Optional Keyword Arguments:
            ==========   =============
            keyword      description
            ==========   =============

            nocat        [False] if True, do not add additional sources from a catalog
            minflux      [1e-8] Minimum integral flux (ph/cm2/s) for sources more than 5 deg from ROI center to be included
            free_radius  [2] background point sources within this radius (deg) are allowed
                             to vary in the fit; others are fixed at the catalog values
            prune_radius [0.1] removes catalog sources within this distance of a user-
                               defined source; degrees
            bg_smodels   [None]  a list of spectral models to replace the default ones in ConsistentBackground
                                 i.e., a custom set of spectral scaling models
            glat         [None]  the Galactic latitude of the source; sets default free parameters in diffuse
            fit_emin     [100,100] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            no_roi       [False] If set, return a ps_manager, roi_manager instead (for separate generation of an ROIAnalysis)
            user_skydir  [None] A user-provided SkyDir that will be the center of the ROI
            ==========   =============
        """

        if point_sources is None and self.roi_dir is None:
            print 'No direction specified!  Provide a point source or set roi_dir member of this object!'
            return

        # Easier to store point_sources as an empty list then as None for later loops
        if point_sources is None: point_sources=[]

        # process kwargs
        glat,bg_smodels,nocat,minflux,free_radius,prune_radius,user_skydir = None,None,False,1e-8,2,0.1,None
        #for key in ['glat','bg_smodels','nocat','minflux','free_radius','prune_radius']
        if 'glat'        in kwargs.keys(): glat        = kwargs.pop('glat')
        if 'bg_smodels'  in kwargs.keys(): bg_smodels  = kwargs.pop('bg_smodels')
        if 'nocat'       in kwargs.keys(): nocat       = kwargs.pop('nocat')
        if 'minflux'     in kwargs.keys(): minflux     = kwargs.pop('minflux')
        if 'free_radius' in kwargs.keys(): free_radius = kwargs.pop('free_radius')
        if 'prune_radius'in kwargs.keys(): prune_radius= kwargs.pop('prune_radius')
        if 'user_skydir' in kwargs.keys(): user_skydir = kwargs.pop('user_skydir')

        # setup backgrounds and point sources
        skydir        = user_skydir or (self.roi_dir if point_sources==[] else point_sources[0].skydir)
        cb            = ConsistentBackground(self.ae,self.background)
        cb.cm.free_radius  = free_radius;
        cb.cm.prune_radius = prune_radius;
        cb.cm.min_flux     = minflux
        bgmodels      = cb.get_bgmodels(models=bg_smodels,lat=glat) if bgmodels is None else bgmodels
        point_sources = point_sources if nocat else cb.cm.merge_lists(skydir,self.maxROI+5,point_sources)

        # try to read in a previous fit
        if previous_fit is not None:
            from roi_analysis import read_fit
            source_list,bg = read_fit(previous_fit)
            for i in xrange(len(bg)):
                backgrounds[i].smodel = bg[i]

        ps_manager = ROIPointSourceManager(point_sources,skydir,quiet=self.quiet)
        bg_manager = ROIBackgroundManager(self,bgmodels,skydir,quiet=self.quiet)

        # if didn't specify a source, pick closest one and make it free -- maybe remove this?
        if point_sources==[] and (not N.any([N.any(m.free) for m in ps_manager.models])):
            if ps_manager.models.shape[0] != 0:
                ps_manager.models[0].free[:] = True

        # n.b. weighting PSF for the central point source
        self.psf.set_weights(self.ltcube,skydir)

        if no_roi: return ps_manager,bg_manager

        return ROIAnalysis(skydir,ps_manager,bg_manager,self,**kwargs)

    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        ignore = ('psf', 'exposure')
        for key in sorted(self.__dict__.keys()):
            if key in self.ae.__dict__ or key in ignore: continue # avoid duplication internal functions
            s += '\t%-20s: %s\n' %(key, self.__dict__[key])
        return s



