"""
Manage data and livetime information for an analysis


    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pixeldata.py,v 1.9 2009/05/09 18:39:31 burnett Exp $

"""
version='$Revision: 1.9 $'.split()[1]
import os
import psf
import math
import skymaps
import pointlike


class PixelData(object):

    def __init__(self, event_files, history_files, **kwargs):
        """
        call signature::

  data = PixelData(event_files, history_files,  **kwargs)


Create a new PixelData instance, managing data and livetime.  

    event_files: a list of event files (FT1) to process (or a single file). 
    history_files: a list of spacecraft history files (FT2) to process (or a single file).

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  roi_dir     [None] direction to use if exposure calculation will be limited; otherwise use all sky                                                                                  
  roi_radius  [ 25]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)                                                                                
  livetimefile [None] Exposure file: if specified and not found, create from FT2/FT1 info  
  zenithcut   [105]  Maximum spacecraft pointing angle with respect to zenith to allow
  align       [False] if True, perform boresigh alignment correction on pre-September flight data
  
  datafile    [None]  HEALpix data file: if specified and not found, create from FT1 info                                                                               
  class_level [ 3]  select class level (set 0 for gtobssim )
  mc_src_id   [-1] if set to != -1, select on MC_SRC_ID column in FT1 file
  use_mc_energy [False] set True to use MC_ENERGY
  binsperdecade [4] When generating Bands from the FT1 data. 
  CALDB       [None] If not specified, will use environment variable
  quiet       [False] Set True to suppress (some) output 
  verbose     [False] More output

  =========    KEYWORDS CONTROLLING INSTRUMENT RESPONSE
  CALDB       [None] If not specified, will use environment variable
  psf_irf     ['P6_v3_diff'] Which IRF to use for the point spread function

  tstart      [0] Default no cut on time; otherwise, cut on MET > tstart
  tstop       [0] Default no cut on time; otherwise, cut on MET < tstop

  emin        [100] Left edge of minimum energy bin in MeV
  emax        [5e5] Some ill-defined notion of max energy; maximum bin center will be less than this
  =========   =======================================================
"""

        self.event_files = event_files if type(event_files)==type([]) or event_files is None else [event_files]
        self.history_files=history_files if type(history_files)==type([]) or history_files is None else [history_files]

        self.roi_dir     = None
        self.roi_radius  = 25
        self.livetimefile= None
        self.datafile    = None 
        self.zenithcut   = 105
        self.binsperdecade=4
        self.tstart      = 0
        self.tstop       = 0
        self.emin        = 100.
        self.emax        = 5e5
        self.quiet       = False
        self.verbose     = False
        self.class_level = 3  # select class level
        self.mc_src_id   = -1
        self.use_mc_energy = False
        self.CALDB       = None
        self.psf_irf     = 'P6_v3_diff'
                  
        for key,value in kwargs.items():
            if key in self.__dict__:
                self.__dict__[key] = value
      	
        from numpy import arange,log10
        self.my_bins = 10**arange(log10(self.emin),log10(self.emax),1./self.binsperdecade)
        
        # check explicit files
        for filelist in  [self.event_files, self.history_files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise Exception('PixelData setup: file name or path "%s" not found'%filelist[0])

            #if filelist is None or len(filelist)==0:
	         #   raise Exception('PixelData setup: received empty list of event or history files')


        self.data= self.get_data()
        self.dmap= self.data.map()
        self.PSF_setup(self.dmap)
        self.lt =  self.get_livetime()

    def fill_empty_bands(self,bpd):
      
        from pointlike import Photon
        from skymaps import SkyDir
        dummy = SkyDir(0,0)
        bands = self.my_bins
        
        for bin_center in (bands[:-1]*bands[1:])**0.5:
             ph_f = Photon(dummy,bin_center,2.5e8,0)
             ph_b = Photon(dummy,bin_center,2.5e8,1)
             bpd.addPhoton(ph_f,0)
             bpd.addPhoton(ph_b,0)
             
    def Data_setup(self):

        from numpy import arccos,pi
        pointlike.Data.set_class_level(self.class_level)
        pointlike.Data.set_zenith_angle_cut(self.zenithcut)
        pointlike.Data.set_theta_cut(arccos(0.4)*180./pi)
        pointlike.Data.set_use_mc_energy(self.use_mc_energy)
        if not self.quiet: print 'Set Data theta cut at %.2f'%(arccos(0.4)*180./pi)

        self.binner = skymaps.PhotonBinner(self.my_bins) # per decade
        pointlike.Data.setPhotonBinner(self.binner)

    def PSF_setup(self,bpd):

        # modify the psf parameters in the band objects, which SimpleLikelihood will then use
        if self.CALDB is not None: skymaps.IParams.set_CALDB(self.CALDB)
        else: skymaps.IParams.set_CALDB(os.environ['CALDB'])
        skymaps.IParams.init('_'.join(self.psf_irf.split('_')[:-1]),self.psf_irf.split('_')[-1])
        bpd.updateIrfs()  
    
    def get_data(self):
        
        #if no binned object present, create; apply cuts
        if self.datafile is None or not os.path.exists(self.datafile):

            self.Data_setup()

            if self.verbose: print 'loading file(s) %s' % self.event_files
            data = pointlike.Data(self.event_files,-1,self.tstart,self.tstop,self.mc_src_id,'')

            self.fill_empty_bands(data.map())     # fill any empty bins           

            if self.verbose: print 'done'
            if self.datafile is not None:
                if not self.quiet: print 'saving datafile %s for subsequent use' % self.datafile
                data.write(self.datafile)
        
        else:
            data = pointlike.Data(self.datafile)
            if not self.quiet: print 'loaded datafile %s ' % self.datafile
        
        if self.verbose: data.map().info()


        if self.verbose: print '---------------------'
        return data

    def get_livetime(self,   pixelsize=1.0):
        if self.livetimefile is None or not os.path.exists(self.livetimefile):
            if self.roi_dir is None:
                # no roi specified: use full sky
                self.roi_dir = skymaps.SkyDir(0,0)
                self.roi_radius = 180

            lt = skymaps.LivetimeCube(
                cone_angle=self.roi_radius,\
                dir=self.roi_dir,\
                zcut=math.cos(math.radians(self.zenithcut)),\
                pixelsize=pixelsize,\
                quiet=self.quiet )

            #I think this code should handle cases with abritrary numbers of FT1/FT2 files
            #g = skymaps.Gti(self.event_files[0])
            #for n in xrange(1,len(self.event_files)):
            #need a union here, not intersection!g.intersection(skymaps.Gti(self.event_files[n]))
            #for hist in self.history_files:
            #   lt.load(hist,g)
            
            # if multiple history files, assume that there is one per event file, 
            # and apply the Gti (unsafe!)
            
            for n, hist in enumerate(self.history_files):
                lt.load(hist,  skymaps.Gti(self.event_files[n]))
                
            # write out ltcube if requested
            if self.livetimefile is not None: lt.write(self.livetimefile)
        else: 
            # ltcube exists: just use it! (need to bullet-proof this)
            lt = skymaps.LivetimeCube(self.livetimefile)
            print 'loaded LivetimeCube %s ' % self.livetimefile
        return lt


#--------------------------------------------------------------------
# test for this application
#--------------------------------------------------------------------

def setup(bins=8):
    data_path = r'f:\glast\data\flight'
    months = ['aug', 'sep', 'oct']
    event_files  = [os.path.join(data_path, '%s2008-ft1.fits'%m) for m in months ]
    history_files= [os.path.join(data_path, '%s2008-ft2.fits'%m) for m in months ]
    return PixelData( event_files, history_files,  
        livetimefile='../data/aug-oct_livetime.fits', 
        datafile='../data/aug-oct2008_%d.fits'%bins, 
        binsperdecade=bins)
if __name__=='__main__':
    pass
