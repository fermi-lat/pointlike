"""
Manage the psf, module in uw/irfs

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/irfs/psfman.py,v 1.4 2018/01/27 15:35:07 burnett Exp $

"""
__version__='$Revision: 1.4 $'

import os, glob
from astropy.io import fits
import numpy as np
from scipy.optimize import fmin

from skymaps import (SkyDir, ExposureWeighter, PythonPsf,)
from . import caldb


class PSFmanager(dict):
    """ manage the PSF
    Much of this code has been imported from the original uw/like/pypsf, writen by M. Kerr, 
     modified in order to expand to the PSFn event types. 
    It uses uw/irfs/caldb to interpret the CALDB information.
    Usage: 
        pman = psfman.PSFmanager() # sets up default CALDB, from $CALDB
        psf = pman(event_type, energy) # returns functor describing the PSF
        psf.plot() # look at it
    
    """
    def __init__(self, caldb, livetimefile=''):
        """caldb   : CALDB object
           livetimefile : filename of the livetime cube file
                used in the weighting, if not empty
                Could have a wild card
        """
        self.cdbman = caldb
        self.hdus = None
        if livetimefile !='':
            self.ltfiles = sorted(glob.glob(livetimefile))
            assert len(self.ltfiles)>0, 'Expected valid livetime file'
            assert len(self.ltfiles)==1 or len(self.ltfiles)==4, 'Expected 1 or 4 files to match'
        else: self.ltfiles=[]

    def _ltfile(self, event_type):
        if len(self.ltfiles)==0: return ''
        if event_type<2 or len(self.ltfiles)==1: return self.ltfiles[0]
        return self.ltfiles[event_type-2]

    def __repr__(self):
        return '{self.__class__}\nCALDB: {self.cdbman}\n\tlivetime: {self.livetimefile}'.format(self=self)

    def psf_base_integral(self,g,s,dmax,dmin=0):
        """Integral of the PSF base function; g = gamma, s = sigma (scaled),
           delta = deviation in radians."""
        u0 = 0.5 * np.outer(dmin,1./s)**2
        u1 = 0.5 * np.outer(dmax,1./s)**2
        return (1+u0/g)**(1-g)-(1+u1/g)**(1-g)

    def setup_event_type(self, et):
        """ Load info into superclass dictionary for the eventtype
        """
        if et in self.keys(): return 
        self[et]=dict()
        info =self.cdbman('psf', event_type=et)
        ext = info['extensions']
        filename = info['filename']
        if self.hdus is None or self.hdus.filename() != filename:
            self.hdus = fits.open(filename)
            #print ('opened "{}"'.format(os.path.split(filename)[-1]))
        # scale factor from PSF_SCALING hdu    
        sft = np.asarray(self.hdus[ext['PSF_SCALING']].data.field('PSFSCALE')).flatten()
        self[et]['scale_func']=\
             scale_func = lambda e: ( (sft[0]*(e/100.)**(sft[-1]))**2 + sft[1]**2 )**0.5
        data = self.hdus[ext['RPSF']].data
        
        self.e_los = np.asarray(data.field('ENERG_LO')).flatten().astype(float)
        self.e_his = np.asarray(data.field('ENERG_HI')).flatten().astype(float)
        self.c_los = np.asarray(data.field('CTHETA_LO')).flatten().astype(float)
        self.c_his = np.asarray(data.field('CTHETA_HI')).flatten().astype(float)
        ne,nc = len(self.e_los),len(self.c_his)
        
        def proc(pname):
            t = np.reshape(data.field(pname), [nc,ne])
            return t.transpose()[:,::-1]
        tables = np.asarray([proc(x) for x in 'NCORE NTAIL GCORE GTAIL SCORE STAIL'.split()])

        # normalize according to the irfs/latResponse conventions
        ens = (self.e_los*self.e_his)**0.5
        for i,en in enumerate(ens): # iterate through energy
            sf = scale_func(en)

            # vector operations in incidence angle
            nc,nt,gc,gt,sc,st = tables[:,i,:]
            normc = self.psf_base_integral(gc,sc*sf,np.pi/2)[0]
            normt = self.psf_base_integral(gt,st*sf,np.pi/2)[0]

            # NB leave scale factor out here so we can adjust norm to
            # a particular energy (effectively cancelling in integral)
            norm = (2*np.pi * (normc*sc**2 + normt*nt*st**2))**-1
            tables[0,i,:] = norm # adjust NCORE
            tables[1,i,:]*= norm # set to NTAIL*NCORE

        self[et]['tables'] =tables
        
        #### Calculate the weights -- depend in principle on direction and livetime
        self[et]['weights'] = self._calculate_weights(et)
    
    def _calculate_weights(self,event_type, skydir=None):
        """Mostly a copy of like.pypsf.Psf.__calc_weights__
        Adjusted to do a single event type
        """
        aeff_info = self.cdbman('aeff',event_type=event_type)
        filename=aeff_info['filename']
        # Have to redundantly open to get the ext name
        aeff_hdus = fits.open(aeff_info['filename'])
        extname = aeff_hdus[aeff_info['extensions']['EFF_AREA']].header['EXTNAME']

        elo,ehi,clo,chi = self.e_los,self.e_his,self.c_los[::-1],self.c_his[::-1]

        # note that it is designed to do front and back, so need to copy filenames and extname
        ew = ExposureWeighter(filename, filename, self._ltfile(event_type), extname, extname)
        
        dummy     = skydir or SkyDir() 
        weights  = np.zeros([len(elo),len(clo)])
        for i in xrange(len(elo)):                   # iterator over energies
            em = (elo[i]*ehi[i])**0.5
            for k,(c0,c1) in enumerate(zip(clo,chi)):# iterator over cos(theta) on-axis to edge
                if c0 < 0.35: continue               # exclude bins below cos(theta)=0.4
                weights[i,k] = ew(c0,c1,elo[i],ehi[i],0,dummy) # exposure at energy/cos(theta)
            weights[i,:] /= weights[i,:].sum()       # normalize weights over cos(theta)
        return weights

    def get_p(self,e,ct):
        # Identical to old except for indexing of ct
        ind    = min(np.searchsorted(self.e_his,e),len(self.e_his) - 1)
        p      = self[ct]['tables'][:,ind,:]
        w      = self[ct]['weights'][ind,:]
        return np.append(p,[w],axis=0)

    def get_cpp_psf(self,e,ct):
        sf = self[ct]['scale_func'](e)
        nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
        return PythonPsf(sc*sf,st*sf,gc,gt,nc/sf**2,nt/sf**2,w)

    def psf_base(self,g,s,delta):
        """Implement the PSF base function; 
            g = gamma, 
            s = sigma (scaled), 
            delta = deviation in radians.
        
            Operation is vectorized both in parameters and angles; return
            is an array of shape (n_params,n_angles).
        """
        u = 0.5 * np.outer(delta,1./s)**2
        return (1-1./g)*(1+u/g)**(-g)


    def __call__(self, event_type, energy=1000):
        """Return a PSF functor of distance for the given event type (e.g., front or back)
            the energy can be changed by a setEnergy method

        """
        self.setup_event_type(event_type)
        psfman = self
        
        class BandPSF(object):
            
            def __init__(self,  event_type, energy):
                self.psfman = psfman # reference for clients that need it, like convolution
                self.event_type =event_type
                self.setEnergy(energy)
                self.etname = psfman.cdbman.event_type_names[event_type]
                
            def setEnergy(self, energy):
                self.energy=energy
                self.par = psfman.get_p(energy, self.event_type)
                self.scale = psfman[self.event_type]['scale_func'](energy)
                 # scale sigma parameters
                self.par[4] *= self.scale
                self.par[5] *= self.scale
                # also correct normalization
                self.par[0] /= self.scale**2
                self.par[1] /= self.scale**2

                # Get the Python function for internal use
                nc,nt,gc,gt,sc,st,w = self.par
                self._cpsf = PythonPsf(sc,st,gc,gt,nc,nt,w)
                self.r68= self.inverse_integral()
 
            def __repr__(self):
                return '{self.__class__} PSF for event_type {self.event_type}[{self.etname}], energy  {self.energy:.0f} MeV, R68 {self.r68:.2f} deg'.format(self=self)
 
            def __call__(self, delta):
                nc,nt,gc,gt,sc,st,w = self.par
                yc = psfman.psf_base(gc,sc,delta)
                yt = psfman.psf_base(gt,st,delta)
                return (w*(nc*yc + nt*yt)).sum(axis=1)

            def integral(self, dmax, dmin=0):
                """ this does not seem to be the integral over solid angle
                Note -- does *not* take a vector argument right now.
                """
                # more expressive
                #return integrate.quad(lambda r: 2*np.pi*r*self(r), dmin, dmax)[0]
                TWOPI = 2.*np.pi
                nc,nt,gc,gt,sc,st,w = self.par
                icore = TWOPI*sc**2*nc*\
                        psfman.psf_base_integral(gc,sc,dmax,dmin=dmin)
                itail = TWOPI*st**2*nt*\
                        psfman.psf_base_integral(gt,st,dmax,dmin=dmin)
                return (w*(icore+itail)).sum()


            def inverse_integral(self, percent=68,on_axis=False): 
                """Return the radius at which the integral PSF reaches the required percentage.

                NB -- returned value is in degrees.
                """
                RAD2DEG = 180./np.pi
                e,ct = self.energy, self.event_type
                percent = float(percent)/100
                if on_axis:
                    sf = self.scale
                    nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
                    u1 = gc[0]*( (1-percent)**(1./(1-gc[0])) - 1)
                    u2 = gt[0]*( (1-percent)**(1./(1-gt[0])) - 1)
                    return RAD2DEG*sf*(nc[0]*(u1*2)**0.5*sc[0] + nt[0]*(u2*2)**0.5*st[0])           
                # off axis
                f = lambda x: abs(self.integral(x) - percent)

                seeds = np.asarray([5,4,3,2.5,2,1.5,1,0.5,0.25])*self.scale
                seedvals = np.asarray([self.integral(x) for x in seeds])
                seed = seeds[np.argmin(np.abs(seedvals-percent))]
                trial = fmin(f,seed,disp=0,ftol=0.000001,xtol=0.01)
                if trial > 0:
                    return trial[0]*RAD2DEG
                print ('Warning: could not invert integral; return best grid value.')
                return seed*RAD2DEG

            def overlap(self, roi_dir, radius, skydir): #not tested
                #return pypsf.PsfOverlap()(roi_dir, self.sd, skydir) 
                return self._cpsf.overlap_circle(roi_dir, np.radians(radius), skydir)

            def wsdl_value(self, skydir, wsdl):
                """Return a list of PSF values corresponding to the difference of skydir and the list
                    paramters:
                    skydir : SkyDir object,
                    wsdl   : WeightedSkyDirList: a list of SkyDir
                """
                rvals = np.empty(len(wsdl), dtype=float)
                self._cpsf.wsdl_val(rvals, skydir, wsdl)
                return rvals

            def plot(self, ax=None):
                from matplotlib import pylab as plt
 
                if ax is None:
                    fig,ax = plt.subplots(figsize=(5,5))
                else: fig = ax.figure
                
                pct = np.linspace(1,95,95)
                rpc = np.array([self.inverse_integral(z) for z in pct])
                x = np.linspace(0,rpc[-1],51)
                ax.plot(x, self(np.radians(x))/self(0), lw=2);
                rax =ax.twinx()
                rax.plot(rpc, pct, color='green', lw=2)
                rax.set_ylabel('percent integral', color='green')
                ax.set_ylabel('relative density', color='blue')
                ax.set_xlabel('angle [deg]')
                ax.grid(alpha=0.5)
                ax.set_xlim(0,x[-1])
                ax.set_title('{:.0f} MeV {}'.format(self.energy,self.etname));  
                fig.set_facecolor('white')
                return fig     

               
        return BandPSF(event_type, energy)


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
    ax.grid(alpha=0.5)
    ax.set_xlim(0,x[-1])
    ax.set_title('{:.0f} MeV {}'.format(energy, ['front','back'][evt]));       
