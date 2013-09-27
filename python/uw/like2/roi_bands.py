"""
Implements classes encapsulating an energy/conversion type band.  These
are the building blocks for higher level analyses.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roi_bands.py,v 1.4 2013/02/13 03:38:54 burnett Exp $

author: Matthew Kerr
extracted from like.roi_bands, v 1.28 for like2 by T Burnett
"""

import numpy as np
from skymaps import SkyDir, WeightedSkyDirList 
from pointlike import DoubleVector
from uw.like import pypsf
from scipy import optimize # for brentq

###======================================================================###

class ROIBand(object):
    """Wrap a skymaps.Band object, and provide functionality for like2 version of likelihood."""

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
        self.ec = self.ct  = band.event_class()  & 1   #event class bit       
        self.exp = self.sa.exposure.exposure[self.ct]

        self.__setup_data__()
        self.__setup_sp_simps__()

        self.psf = self.sa.psf.band_psf(self,adjust_mean=ROIBand.ADJUST_MEAN)


    def __setup_data__(self):
        """Get all pixels within the ROI in this band."""

        mi,ma                  = np.asarray([self.sa.minROI,self.sa.maxROI])*(np.pi/180.)
        # note use of band sigma for consistency in photon selection!
        self.radius_in_rad  = max(min((2*self.umax)**0.5*self.b.sigma(),ma),mi)
        self.wsdl           = WeightedSkyDirList(self.b,self.sd,self.radius_in_rad,False)
        self.pix_counts     = np.asarray([x.weight() for x in self.wsdl]) if len(self.wsdl)>0 else []
        self.npix           = self.wsdl.total_pix()
        self.has_pixels     = self.wsdl.counts()>0
        self.solid_angle    = self.npix*self.b.pixelArea() #ragged edge
                 

    def __setup_sp_simps__(self):
        """Cache factors for quickly evaluating the counts under a given spectral model, 
            for integrating over the exposure within the energy limits
            note that exposure is evaluated at the skydir
        """

        self.sp_points = sp = np.logspace(np.log10(self.emin),np.log10(self.emax),self.nsp_simps+1)
        exp_points     = map(lambda e: self.exp.value(self.sd, e), sp)
        # following may be marginally faster, but not implemented for exposure cube
        #exp_points      = np.asarray(self.exp.vector_value(self.sd,DoubleVector(sp)))
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
        return pypsf.PsfOverlap()(self, self.sd, skydir)  
        
    def exposure(self, skydir, energy=None):
        """ return the exposure at the position and the given energy, or (default) the central energy"""
        return self.exp.value(skydir, energy if energy is not None else self.e)
        
    def exposure_integral(self, model_function, axis=None):
        """ return integral over exposure for function of differential flux """
        return (model_function(self.sp_points)*self.sp_vector).sum(axis=axis)

    def expected(self,model):
        """Integrate the passed spectral model over the exposure, at the given direction and return expected counts."""
        return (model(self.sp_points)*self.sp_vector).sum()
        
    def optimum_energy(self, f):
        """return optimum energy to determine integral over exposure as f(e)*exp(e)*(b-a)
            where exp is the exposure function
          f: function of energy
        """
        a,b = self.emin, self.emax
        alpha = np.log(f(a)/f(b))/np.log(b/a)
        fn = lambda x: x**-alpha
        t = self.expected(fn) / (b-a)
        ft = lambda e : t / (fn(e) * self.exp.value(self.sd, e) ) -1
        #print alpha, ft(a), ft(b)
        if ft(a)*ft(b)>0: return (a+b)/2. # no optimum, flat? use mean
        opte = optimize.brentq(ft, a, b)
        return opte if opte>a and opte<b else np.sqrt(a*b) # protection

