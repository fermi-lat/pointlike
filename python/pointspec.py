"""  A module to provide simple and standard access to pointlike fitting and spectral analysis.  The
     relevant parameters are fully described in the docstring of the constructor of the SpectralAnalysis
     class.
    
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointspec.py,v 1.33 2009/09/20 20:43:19 kerrm Exp $

    author: Matthew Kerr
"""
version='$Revision: 1.33 $'.split()[1]
import os
import sys

from pixeldata import PixelData
from pypsf     import Psf,OldPsf,NewPsf
from pointspec_helpers import ExposureManager,ConsistentBackground

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
      self.CALDB   = r'f:/glast/caldb/v1r1/CALDB/data/glast/lat'

      self.ft1files = None
      self.ft2files = None
      self.ltcube   = None
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
      except IOError:
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
  new_psf        [None] if not None, must point to a piclkle of a fit to allGamma data
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

    """
    
    def __init__(self, analysis_environment, **kwargs):
        """

Create a new spectral analysis object.  

    analysis_environment: an instance of AnalysisEnvironment correctly configured with
                          the location of files needed for spectral analysis (see its
                          docstring for more information.)

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
  binsperdec   [4] energy binning granularity when binning FT1
  emin         [100] Minimum energy                                                                                 
  emax         [3e5] Maximum energy

  =========    KEYWORDS FOR MONTE CARLO DATA
  mc_src_id    [ -1] set to select on MC_SRC_ID column in FT1
  mc_energy    [False] set True to use MC_ENERGY instead of ENERGY

  =========    KEYWORDS CONTROLLING INSTRUMENT RESPONSE
  irf          ['P6_v3_diff'] Which IRF to use
  new_psf      [None] if not None, must point to a piclkle of a fit to allGamma data
  
  =========    KEYWORDS CONTROLLING SPECTRAL ANALYSIS
  background   ['ems_ring'] - a choice of global model specifying a diffuse background; see ConsistentBackground for options
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
        self.emin        = 100    # MeV
        self.emax        = 3e5    # MeV

        self.mc_src_id   = -1
        self.mc_energy   = False

        self.irf         = 'P6_v3_diff'
        self.new_psf     = None

        self.background  = 'ems_ring'
        self.maxROI      = 10    # deg
        self.minROI      = 5     # deg

        self.quiet       = False
        self.verbose     = False

        self.ae          = analysis_environment

        self.__dict__.update(**kwargs)
        self.__dict__.update(analysis_environment.__dict__)

         #TODO -- sanity check that BinnedPhotonData agrees with analysis parameters
        self.pixeldata = PixelData(self)
        self.exposure  = ExposureManager(self)
        if self.new_psf is None:
            self.psf = OldPsf(self.ae.CALDB,irf=self.irf)
        else:
            self.psf = NewPsf(self.ae.CALDB,self.new_psf,irf=self.irf)



    def roi(self, point_sources = None, bgmodels = None, previous_fit = None, no_roi = False, **kwargs):
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
            bg_smodels   [None]  a list of spectral models to replace the default ones in ConsistentBackground
                                 i.e., a custom set of spectral scaling models
            glat         [None]  the Galactic latitude of the source; sets default free parameters in diffuse
            fit_emin     [100,100] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            ==========   =============
      """
        
        from roi_managers import ROIPointSourceManager,ROIBackgroundManager
        from roi_analysis import ROIAnalysis
        
        # process kwargs
        glat,bg_smodels,nocat = None,None,False
        if 'glat'        in kwargs.keys(): glat   = kwargs.pop('glat')
        if 'bg_smodels'  in kwargs.keys(): models = kwargs.pop('bg_smodels')
        if 'nocat'       in kwargs.keys(): nocat  = kwargs.pop('nocat')

        # setup backgrounds and point sources
        cb            = ConsistentBackground(self.ae,self.background)
        bgmodels      = cb.get_bgmodels(models=bg_smodels,lat=glat) if bgmodels is None else bgmodels
        point_sources = point_sources if nocat else cb.cm.merge_lists(self.roi_dir,self.maxROI+5,point_sources)

        # try to read in a previous fit
        if previous_fit is not None:
            from roi_analysis import read_fit
            source_list,bg = read_fit(previous_fit)
            for i in xrange(len(bg)):
                backgrounds[i].smodel = bg[i]

        ps_manager = ROIPointSourceManager(point_sources,self.roi_dir)
        bg_manager = ROIBackgroundManager(self, bgmodels)

        # if didn't specify a source, pick closest one and make it free -- maybe remove this?
        if point_sources is None and (not N.any([N.any(m.free) for m in ps_manager.models])):
            ps_manager.models[0].free[:] = True
            
        # n.b. weighting PSF for the central point source
        self.psf.set_weights(self.ltcube,point_sources[0].skydir)

        if no_roi: return ps_manager,bg_manager

        return ROIAnalysis(ps_manager,bg_manager,self,**kwargs)


