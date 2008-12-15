"""
Manage data and livetime information for an analysis


    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pixeldata.py,v 1.2 2008/11/28 20:43:34 burnett Exp $

"""
version='$Revision: 1.2 $'.split()[1]
import os

class PixelData(object):

    def __init__(self, event_files, history_files, **kwargs):
        """
        call signature::

  data = PixelData(event_files, history_files,  **kwargs)


Create a new PixelData instance, managing data and livetime.  

    event_files: a list of event files (FT1) to process (or a single file). 
    history_files: a list of spacecraft history files (FT2) to process (or a single file).
    diffusefile:  Full path to a galactic diffuse file, like  GP_gamma_conventional.fits  (note isotropic)                                 

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
  class_level [ 3]  select class level (set 0 for gtobssim                                                             
  binsperdecade [4] When generating Bands from the FT1 data.
  use_mc_psf  [False] Use PSF determined by MC analysis of true; otherwise as defined by data
  quiet       [False] Set True to suppress (some) output 
  verbose     [False] More output
  =========   =======================================================
"""


        self.event_files = event_files if type(event_files)==type([]) else [event_files] 
        self.history_files=history_files if type(history_files)==type([]) else [history_files]
        self.roi_dir     = None
        self.roi_radius  = 25
        self.livetimefile= None
        self.datafile    = None 
        self.zenithcut   = 105
        self.binsperdecade=4
        self.quiet       = False
        self.verbose     = False
        self.class_level = 3  # select class level
        self.use_mc_psf  = False
        
        for key,value in kwargs.items():
            if key in self.__dict__:
                self.__dict__[key] = value
	    
        # check explicit files
        for filename in  [self.event_files[0], self.history_files[0]] :
            if filename is not None and not os.path.exists(filename):
                raise Exception('PixelData setup: file name or path "%s" not found'%filename)

        self.lt =  self.get_livetime()
        self.data= self.get_data()
        self.dmap= self.data.map()

    
    def get_data(self):
        from skymaps import PhotonBinner
        from pointlike import Data

        Data.set_class_level(self.class_level)
        Data.set_zenith_angle_cut(self.zenithcut)

        if self.datafile is None or not os.path.exists(self.datafile):
            self.binner = PhotonBinner(self.binsperdecade) # per decade
            Data.setPhotonBinner(self.binner)
            if self.verbose: print 'loading file(s) %s' % self.event_files
            data = Data(self.event_files)
            if self.verbose: print 'done'
            if self.datafile is not None:
                if not self.quiet: print 'saving datafile %s for subsequent use' % self.datafile
                data.write(self.datafile)
        else:
            data = Data(self.datafile)
            print 'loaded datafile %s ' % self.datafile

        # modify the psf parameters in the band objects, which SimpleLikelihood will then use
        from psf import PSF
        import math
        self.psf = PSF(use_mc=self.use_mc_psf)
        
        if not self.quiet: print 'setting PSF parameters (use_mc=%d)'%self.use_mc_psf
        if self.verbose: print '  energy class  gamma sigma(deg)'
        for band in data.map():
             e = (band.emin()*band.emax())**0.5
             cl = band.event_class()
             gamma = self.psf.gamma(e,cl)
             sigma = self.psf.sigma(e,cl)
             if self.verbose: print '%6.0f%5d%10.1f%10.2f' % (e, cl, gamma, sigma)
             band.setGamma(gamma)
             band.setSigma(math.radians(sigma))

        
        return data

    def get_livetime(self,   pixelsize=1.0):
        from skymaps import LivetimeCube, Gti, SkyDir
        from math import cos,radians
        if self.livetimefile is None or not os.path.exists(self.livetimefile):
            if self.roi_dir is None:
                # no roi specified: use full sky
                self.roi_dir = SkyDir(0,0)
                self.roi_radius = 180
            lt = LivetimeCube(
                cone_angle=self.roi_radius,\
                dir=self.roi_dir,\
                pixelsize=pixelsize,\
                zcut=cos(radians(self.zenithcut)),\
                quiet=self.quiet )
            #Not right?  Needs fixin'.
            #g = Gti(self.event_files[0])
            #for n in xrange(1,len(self.event_files)):
            #   g = g.intersect(Gti(self.event_files[n]))
            #for hist in self.history_files:
            #   lt.load(hist,g)
            for n, hist in enumerate(self.history_files):
                lt.load(hist,  Gti(self.event_files[n]))
            if self.livetimefile is not None: lt.write(self.livetimefile)
        else: 
            lt = LivetimeCube(self.livetimefile)
            print 'loaded LivetimeCube %s ' % self.livetimefile
        return lt


#--------------------------------------------------------------------
# test for this application
#--------------------------------------------------------------------

def setup(bins=8):
    import os
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
