"""
Implements classes encapsulating an energy/conversion type band.  These
are the building blocks for higher level analyses.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_bands.py,v 1.28 2011/09/03 18:21:47 burnett Exp $

author: Matthew Kerr
"""

import numpy as N
from skymaps import SkyDir,WeightedSkyDirList,SkyIntegrator,PsfSkyFunction
from pointlike import DoubleVector
from Models import *
from pypsf import *
from specfitter import SpectralModelFitter
from numpy.linalg import inv
from scipy.optimize import fmin,fsolve

###======================================================================###

class ROIBand(object):
    """Wrap a Band object, and provide additional functionality for likelihood."""

    ADJUST_MEAN = True # provide a "static" interface for functionality below
    
    def init(self):
        self.umax      = 50
        self.nsp_simps = 16
        self.bracketing_function = None
        self.phase_factor = 1.

    def __init__(self,band,spectral_analysis,skydir,**kwargs):
        """
        band: a skymaps.Band object
        spectral_analysis: needed for exposure, psf.band_psf, minROI, maxROI 
        skydir a SkyDir object, used to select data from the Band
        """

        self.init()
        self.__dict__.update(**kwargs)
        self.b    = band
        self.emin, self.emax = band.emin(),band.emax()
        self.e    = (self.emin*self.emax)**0.5
        self.sa  = spectral_analysis
        self.sd  = skydir
        self.ec  = self.ct = band.event_class()  & 1          # note change and mask
        self.exp = self.sa.exposure.exposure[self.ct]

        self.__setup_data__()
        self.__setup_sp_simps__()

        self.psf = self.sa.psf.band_psf(self,adjust_mean=ROIBand.ADJUST_MEAN)

        # this is sneaky but should work just fine
        if self.bracketing_function is not None:
            self.phase_factor *= self.bracketing_function(self.e,self.ec)

    def __setup_data__(self):
        """Get all pixels within the ROI in this band."""

        mi,ma                  = N.asarray([self.sa.minROI,self.sa.maxROI])*(N.pi/180.)
        # note use of band sigma for consistency in photon selection!
        self.radius_in_rad = max(min((2*self.umax)**0.5*self.b.sigma(),ma),mi)
        self.wsdl             = WeightedSkyDirList(self.b,self.sd,self.radius_in_rad,False)
        self.pix_counts     = N.asarray([x.weight() for x in self.wsdl]) if len(self.wsdl) else 0.
        self.photons         = self.wsdl.counts()
        self.has_pixels     = self.photons > 0
        self.npix             = self.wsdl.total_pix()
        self.solid_angle    = self.npix*self.b.pixelArea() #ragged edge
        self.solid_angle_p = 2*N.pi*(1-N.cos(self.radius_in_rad)) #solid angle for a pure circle
                 

    def __setup_sp_simps__(self):
        """Cache factors for quickly evaluating the counts under a given spectral model."""

        self.sp_points = sp = N.logspace(N.log10(self.emin),N.log10(self.emax),self.nsp_simps+1)
        exp_points      = N.asarray(self.exp.vector_value(self.sd,DoubleVector(sp)))
        simps_weights  = (N.log(sp[-1]/sp[0])/(3.*self.nsp_simps)) * \
                              N.asarray([1.] + ([4.,2.]*(self.nsp_simps/2))[:-1] + [1.])
        self.sp_vector = sp * exp_points * simps_weights
     
    ### the following functions added to for functionality best performed by an ROIband
    ### they are not used by the (old) likelihood computation
    def pixels_from_psf(self,  skydir):
        """ return pixel array of predicted values given the PSF and direction
        """
        rvals  = np.empty(len(self.wsdl),dtype=float)
        self.psf.cpsf.wsdl_val(rvals, skydir, self.wsdl) #from C++: sets rvals
        return rvals*self.b.pixelArea()
        
    def exposure_integral(self, model_function, axis=None):
        """ return integral over exposure for function of differential flux """
        return (model_function(self.sp_points)*self.sp_vector).sum(axis=axis)
        
    def psf_overlap(self, skydir):
        """ return the overlap of the associated PSF at the given position with this ROI
        """
        return PsfOverlap()(self, self.sd, skydir)  
        
    def exposure(self, skydir, energy=None):
        """ return the exposure at the position and the given energy, or (default) the central energy"""
        return self.exp.value(skydir, energy if energy is not None else self.e)
    ################################################################################

    def reload_data(self,band):
        """For Monte Carlo, when everything stays same but realization of data."""
        self.b = band
        self.__setup_data__()

    def expected(self,model):
        """Integrate the passed spectral model over the exposure and return expected counts."""
        return (model(self.sp_points)*self.sp_vector).sum()

    def gradient(self,model):
        """Calculate the gradient of a spectral model (wrt its parameters) integrated over the exposure."""
        return (model.gradient(self.sp_points)*self.sp_vector).sum(axis=1)

    def bandLikelihood(self, parameters, *args):
        """Implement a model independent likelihood for the number of counts of a particular source.
            Other sources (diffuse, neighboring point sources) are treated as fixed."""

        # N.B. -- it is assumed that new_counts *does not* include exposure ratio
        # this is the easiest way to deal with the exposure ratio, but it must be explicitly set if the *counts*
        # are used elsewhere.  This comes up in, e.g., saving the values for use in localization.
        which = args[0] if len(args) > 0 else 0
        new_counts = parameters[0]*self.er[which]

        old_counts = self.ps_counts[which]

        tot_term = (self.bg_all_counts + self.ps_all_counts + self.overlaps[which]*(new_counts - old_counts))*self.phase_factor

        pix_term = (self.pix_counts * 
                            N.log(
                                self.bg_all_pix_counts + self.ps_all_pix_counts + self.ps_pix_counts[:,which]*(new_counts - old_counts)
                            )
                      ).sum() if self.has_pixels else 0.

        return tot_term - pix_term

    def logLikelihood(self):
        """ Return the (negative) log likelihood for this band.  Certain members
            of this object are set externally and are required here, specifically
            ps_all_counts, ps_all_pix_counts, bg_all_counts, and bg_all_pix_counts.
        
        """
        tot = (self.bg_all_counts + self.ps_all_counts) * self.phase_factor
               
        pix = (self.pix_counts * N.log(self.bg_all_pix_counts + self.ps_all_pix_counts)).sum()\
              if self.has_pixels else 0

        return tot - pix

###======================================================================###

class ROIEnergyBand(object):
    """Wrap 1 or 2 ROIBand objects corresponding to the same energy level 
        but different conversion classes.  Implement a likelihood as a
        function of energy.
        
        Can also accomodate multiple energies, allowing simple merge (but needs some work)
        """

    def __init__(self,bands,emin=None,emax=None):

        self.bands = bands
        self.emin = self.bands[0].emin if emin is None else emin
        self.emax = self.bands[-1].emax if emax is None else emax

    def __rois__(self):
        R2D = 180./N.pi
        froi = self.bands[0].radius_in_rad * R2D
        broi = self.bands[1].radius_in_rad * R2D
        return ( min(froi,broi), max(froi,broi) )

    def spectralString(self,which=None):
        """Return a string suitable for printSpectrum.

            which -- if not None, an array of indices for point sources to fit
        """
        which = which or [0]
        r      = []
        for w in which:
            self.bandFit(which=w)
            self.m[0] = self.uflux
            ul = sum( (b.phase_factor*b.expected(self.m)*b.er[w] for b in self.bands) )
            if self.flux is None:
                r += [0,ul,0]
            else:
                n = ul*self.flux/self.uflux
                r += [n,ul - n, n - ul*self.lflux/self.uflux]

        rois = self.__rois__()
        ph    = int(sum( (b.photons for b in self.bands) ))
        gal  = int(round(sum( (b.bg_counts[0] for b in self.bands) )))
        iso  = int(round(sum( (b.bg_counts[1] for b in self.bands) )))
        values = tuple([int(round(self.emin)),rois[0],rois[1],ph,gal,iso] + r)
        format = '  '.join(['%6i','%6.2f','%6.2f','%7i','%8i','%9i']+['%7.1f +%5.1f -%5.1f']*len(which))
        return format%values
    
    def bandLikelihood(self,parameters,*args):
        m = args[0]
        m.set_parameters(parameters)
        return sum( (b.bandLikelihood([b.expected(m)],*args[1:]) for b in self.bands) )

    def normUncertainty(self,which=0):
        # this is the uncertainty derived by taking the second derivative of the log likelihood
        # wrt the normalization parameter; it is fractional (delta_v/v)
        # this formulation is specifically for the flux density of a power law, which has such nice properties
        tot = 0
        for b in self.bands:
            if not b.has_pixels: continue
            my_pix_counts = b.ps_pix_counts[:,which]*b.expected(self.m)*b.er[which]
            all_pix_counts= b.bg_all_pix_counts + b.ps_all_pix_counts - b.ps_pix_counts[:,which]*b.ps_counts[which] + my_pix_counts
            tot += (b.pix_counts * (my_pix_counts/all_pix_counts)**2).sum()
        return tot**-0.5
            
            
    def bandFit(self,which=0,saveto=None):
        """Fit a model-independent flux to a point source.
          return value of ts for the band
        """

        bad_fit = False
        self.m = PowerLaw(free=[True,False],e0=(self.emin*self.emax)**0.5) # fix index to 2
        f = self.bandLikelihood

        self.fit = fmin(f,self.m.get_parameters(),disp=0,full_output=1,args=(self.m,which))

        def upper_limit():
            ### warning: this code depends on the log10 representation of the flux

            flux_copy = self.m[0]
            zp          = self.bandLikelihood(N.asarray([-20]),self.m,which)

            # NB -- the 95% upper limit is calculated by assuming the likelihood is peaked at
            # 0 flux and finding the flux at which it has fallen by 1.35; this is a two-sided
            # 90% limit, or a one-sided 95% limit -- that's how it works, right?
            # THB 19 June 01
            #This is true for a Gaussian distribution, but when we want limits, the likelihood
            # us usally an exponentially falling tail: for an exponential the number is -log(0.05)=3.0
            def f95(parameters):
                return abs(self.bandLikelihood(parameters,self.m,which) - zp - 3.0)
            
            # for some reason, can't get fsolve to work here.  good ol' fmin to the rescue
            self.uflux = 10**fmin(f95,N.asarray([-11.75]),disp=0)[0]
            self.lflux = None
            self.flux  = None

            self.m[0] = flux_copy

        # if flux below a certain level, set an upper limit
        if self.m[0] < 1e-20:
            bad_fit = True
            upper_limit()

        else:
            try:
                #self.m.set_cov_matrix([[self.normUncertainty,0],[0,0]])
                err = self.normUncertainty(which=which)
            except:
                bad_fit = True
                err = 0 

            self.flux  = self.m[0] 
            self.uflux = self.flux*(1 + err)
            self.lflux = max(self.flux*(1 - err),1e-30)

        if saveto is not None:
            for b in self.bands: b.__dict__[saveto] = (b.expected(self.m)*b.er[which] if not bad_fit else -1)

        if bad_fit:
            self.ts = 0
        else:
            null_ll = sum( (b.bandLikelihood([0],which) for b in self.bands) )
            alt_ll  = sum( (b.bandLikelihood([b.expected(self.m)*b.er[which]],which) for b in self.bands) )
            self.ts = 2*(null_ll - alt_ll)
            
        return self.ts
        
    def energy_flux(self):
        """ return a tuple (flux, lflux, uflux) of the energy flux in ergs cm**-2 s**-1 units
        """
        ec = self.emin*self.emax * 1.602e-6
        return (ec*self.flux, ec*self.lflux, ec*self.uflux)
