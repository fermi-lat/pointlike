"""
Implements classes encapsulating an energy/conversion type band.  These
are the building blocks for higher level analyses.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_bands.py,v 1.28 2011/09/03 18:21:47 burnett Exp $

author: Matthew Kerr
extracted from like.roi_bands, v 1.28 for like2 by T Burnett
"""

import numpy as np
from skymaps import SkyDir, WeightedSkyDirList 
from pointlike import DoubleVector

###======================================================================###

class ROIBand(object):
    """Wrap a Band object, and provide additional functionality for like2 version of likelihood."""

    ADJUST_MEAN = True # provide a "static" interface for functionality below
    
    def init(self):
        self.umax      = 50
        self.nsp_simps = 16

    def __init__(self,band,spectral_analysis,skydir,**kwargs):
        """
        band: a skymaps.Band object
        spectral_analysis: needed for exposure, psf.band_psf, minROI, maxROI 
        skydir a SkyDir object, used to select data from the Band
        """

        self.init()
        self.__dict__.update(**kwargs)
        self.b   = band
        self.emin, self.emax = band.emin(),band.emax()
        self.e   = (self.emin*self.emax)**0.5
        self.sa  = spectral_analysis
        self.sd  = skydir
        self.ct  = band.event_class()  & 1          # note change and mask
        self.exp = self.sa.exposure.exposure[self.ct]

        self.__setup_data__()
        self.__setup_sp_simps__()

        self.psf = self.sa.psf.band_psf(self,adjust_mean=ROIBand.ADJUST_MEAN)


    def __setup_data__(self):
        """Get all pixels within the ROI in this band."""

        mi,ma                  = np.asarray([self.sa.minROI,self.sa.maxROI])*(np.pi/180.)
        # note use of band sigma for consistency in photon selection!
        self.radius_in_rad = max(min((2*self.umax)**0.5*self.b.sigma(),ma),mi)
        self.wsdl             = WeightedSkyDirList(self.b,self.sd,self.radius_in_rad,False)
        self.has_pixels     = self.photons > 0
        self.solid_angle    = self.npix*self.b.pixelArea() #ragged edge
                 

    def __setup_sp_simps__(self):
        """Cache factors for quickly evaluating the counts under a given spectral model."""

        self.sp_points = sp = np.logspace(np.log10(self.emin),np.log10(self.emax),self.nsp_simps+1)
        exp_points      = np.asarray(self.exp.vector_value(self.sd,DoubleVector(sp)))
        simps_weights  = (np.log(sp[-1]/sp[0])/(3.*self.nsp_simps)) * \
                              np.asarray([1.] + ([4.,2.]*(self.nsp_simps/2))[:-1] + [1.])
        self.sp_vector = sp * exp_points * simps_weights
     
    def pixels_from_psf(self,  skydir):
        """ return pixel array of predicted values given the PSF and direction
        """
        rvals  = np.empty(len(self.wsdl),dtype=float)
        self.psf.cpsf.wsdl_val(rvals, skydir, self.wsdl) #from C++: sets rvals
        return rvals*self.b.pixelArea()
        
    def psf_overlap(self, skydir):
        """ return the overlap of the associated PSF at the given position with this ROI
        """
        return PsfOverlap()(self, self.sd, skydir)  
        
    def exposure(self, skydir, energy=None):
        """ return the exposure at the position and the given energy, or (default) the central energy"""
        return self.exp.value(skydir, energy if energy is not None else self.e)

    def expected(self,model):
        """Integrate the passed spectral model over the exposure and return expected counts."""
        return (model(self.sp_points)*self.sp_vector).sum()


