"""
Manage the psf

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/psf.py,v 1.1 2013/11/07 17:21:28 burnett Exp $

"""
import numpy as np
from uw.like import pypsf

class PSFmanager(object):
    """ manage the PSF
    This is an interface to the original pointlike code, to hide some of that complexity,
    providing a clean interface. Especially, it is the only place that the event_type is handled
    """
    def __init__(self, dataset):
        """dataset : DataSet object, used to extract CALDBmanager
        """
        cdbman = dataset.CALDBManager 
        self.cdb_psf = pypsf.CALDBPsf(cdbman)
        self.irfname = cdbman.irf
        
    def __repr__(self):
        return '%s.%s: IRF "%s"' % (self.__module__, self.__class__.__name__, self.irfname)
        
    
    def __call__(self, event_type, energy=1000):
        """Return a PSF functor of distance for the given event type (e.g., front or back)
            the energy can be changed by a setEnergy method
        """
        cdb_psf = self.cdb_psf 
        
        class PSF(object):
            def __init__(self,  event_type, energy):
                self.event_type =event_type
                self.setEnergy(energy)
                self.parent = cdb_psf # reference for clients that need it, like convolution
                
            def setEnergy(self, energy):
                self.energy=energy
                self.cpsf  = cdb_psf.get_cpp_psf(energy,event_type) # for access to C++
                class TBand(object):
                    def __init__(self, energy, event_type, **kwargs):
                        self.e, self.ct = energy, event_type
                        self.__dict__.update(kwargs)
                        
                self.bpsf = pypsf.BandCALDBPsf(cdb_psf, TBand(self.energy, self.event_type),
                        override_en=False,adjust_mean=False)

            def __repr__(self):
                return 'PSF for event_type %d, energy  %.0f MeV' % (self.event_type, self.energy)
 
            def __call__(self, delta):
                return self.bpsf(delta)
            def integral(self, dmax, dmin=0):
                """ this does not seem to be the integral over solid angle"""
                # more expressive
                #return integrate.quad(lambda r: 2*np.pi*r*self(r), dmin, dmax)[0]
                return self.bpsf.integral(dmax, dmin)
            def inverse_integral(self, percent=68,on_axis=False): 
                return self.parent.inverse_integral(self.energy, self.event_type, percent, on_axis)
            
            def overlap(self, roi_dir, radius, skydir): #not tested
                #return pypsf.PsfOverlap()(roi_dir, self.sd, skydir) 
                return self.cpsf.overlap_circle(roi_dir, np.radians(radius), skydir)
               
        return PSF(event_type, energy)

def plot_psf(roi, iband, ax=None):
    from matplotlib import pylab as plt
    if ax is None:
        fig,ax = plt.subplots(figsize=(6,6))
    else: fig = ax.figure
    q=roi[iband]
    psf = q.band.psf
    energy = q.band.energy
    evt = q.band.event_type
    pct = np.linspace(1,95,95)
    rpc = np.array([psf.inverse_integral(z) for z in pct])
    x = np.linspace(0,rpc[-1],51)
    ax.plot(x, psf(np.radians(x))/psf(0));
    rax =ax.twinx()
    rax.plot(rpc, pct)
    rax.set_ylabel('percent integral')
    ax.set_ylabel('relative density')
    ax.set_xlabel('angle [deg]')
    ax.set_title('{:.0f} MeV {}'.format(energy, ['front','back'][evt]));       
