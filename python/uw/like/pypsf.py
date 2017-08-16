"""
A module to manage the PSF from CALDB and handle the integration over
incidence angle and intepolation in energy required for the binned
spectral analysis.
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pypsf.py,v 1.38 2016/06/22 17:02:52 wallacee Exp $
author: M. Kerr

"""

from astropy.io import fits as pf
import numpy as np
from os.path import join
import os
from cPickle import load
from skymaps import ExposureWeighter,SkyDir,PySkyFunction,Hep3Vector,\
    WeightedSkyDirList,PythonUtilities,PythonPsf
from scipy.integrate import quad,simps
from math import cos,sin

def scale_factor(e,c0,c1,exp):
    return ( (c0*(e/100)**(exp))**2 + c1**2 )**0.5

TWOPI = (2*np.pi)
RAD2DEG = 180./np.pi
DEG2RAD = np.pi/180
USEP7 = True

class Psf(object):

    def init(self):
        pass

    def __init__(self,CALDBManager,**kwargs):

        self.CALDBManager = CALDBManager

        self.init()
        self.__dict__.update(kwargs)

        self.__readCALDB__()
        self.__read_params__()
        self.__calc_weights__()

    def __read_params__(self):pass

    def __readCALDB__(self):

        h0,h1 = self.CALDBhandles = [pf.open(x) for x in self.CALDBManager.get_psf()]

        # read in stuff that doesn't depend on conversion type, and set references to front and back data
        
        if len(h0)>=7:
            #New format, front, back in same file (expect h0, h1 to be the same)
            sff = np.asarray(h0['PSF_SCALING_PARAMS_FRONT'].data.field('PSFSCALE')).flatten()
            sfb = np.asarray(h0['PSF_SCALING_PARAMS_BACK' ].data.field('PSFSCALE')).flatten()
            sf=self.scale_factors = [sff[0], sff[1], sfb[0], sfb[1], sff[-1]]
            self.psf_data = (h0['RPSF_FRONT'].data, h0['RPSF_BACK'].data)
        else:
            # old format, separate files for front and back
            sff = np.asarray(h0[2].data.field('PSFSCALE')).flatten()
            self.psf_data = (h0[1].data, h1[1].data)
            
            if len(sff)==5:
                # even older format, both sets in both PSFSCALE tables
                sf=self.scale_factors=sff
            elif len(sff)==3:
                # new format, individual arrays.
                sfb = np.asarray(h1[2].data.field('PSFSCALE')).flatten()
                sf=self.scale_factors = [sff[0], sff[1], sfb[0], sfb[1], sff[-1]]
            else:
                raise Exception('unexpected length of scale factors')
                
        # NB scaling functions in radians
        self.scale_func = [lambda e: ( (sf[0]*(e/100.)**(sf[-1]))**2 + sf[1]**2 )**0.5,
                           lambda e: ( (sf[2]*(e/100.)**(sf[-1]))**2 + sf[3]**2 )**0.5 ]

        # expect binning to be the same, this uses front
        self.e_los            = np.asarray(h0[1].data.field('ENERG_LO')).flatten().astype(float)
        self.e_his            = np.asarray(h0[1].data.field('ENERG_HI')).flatten().astype(float)
        self.c_los            = np.asarray(h0[1].data.field('CTHETA_LO')).flatten().astype(float)
        self.c_his            = np.asarray(h0[1].data.field('CTHETA_HI')).flatten().astype(float)


    def __calc_weights__(self,livetimefile='',skydir=None):

        aeffstrs = self.CALDBManager.get_aeff()
        if len(self.CALDBhandles[0]) >=7:
            # new format
            try:
                ew         = ExposureWeighter(aeffstrs[0],aeffstrs[1],livetimefile, 'EFFECTIVE AREA_FRONT', 'EFFECTIVE AREA_BACK')
            except Exception:
                ew         = ExposureWeighter(aeffstrs[0],aeffstrs[1],livetimefile)
        else:
            ew         = ExposureWeighter(aeffstrs[0],aeffstrs[1],livetimefile)
        
        dummy     = skydir or SkyDir() 

        elo,ehi,clo,chi = self.e_los,self.e_his,self.c_los[::-1],self.c_his[::-1]

        weights  = np.zeros([2,len(elo),len(clo)])
        for i in xrange(len(elo)):                            # iterator over energies
            em = (elo[i]*ehi[i])**0.5
            for j in [0,1]:                                      # iterate over conversion types
                for k,(c0,c1) in enumerate(zip(clo,chi)):# iterator over cos(theta) on-axis to edge
                    if c0 < 0.35: continue                     # exclude bins below cos(theta)=0.4
                    weights[j,i,k] = ew(c0,c1,elo[i],ehi[i],j,dummy) # exposure at energy/cos(theta)
                weights[j,i,:] /= weights[j,i,:].sum()    # normalize weights over cos(theta)

        self.weights = weights

    def set_weights(self,livetimefile,skydir):
        self.__calc_weights__(livetimefile,skydir)

    def get_p(self,e,ct,cthetabin=None):
        ind    = min(np.searchsorted(self.e_his,e),len(self.e_his) - 1)
        p      = self.tables[ct,:,ind,:]
        w      = self.weights[ct,ind,:]
        if cthetabin is None:
            return np.append(p,[w],axis=0)
        else:
            return np.append(p[:,cthetabin],1)

class CALDBPsf(Psf):
    """Handle reading in the parameters from CALDB for multiple formats and implement the PSF functions."""

    def __read_params__(self):        
        """ Set up the internal tables of PSF parameters."""

        h = self.CALDBhandles # tuple with handles to front and back CALDB
        ne,nc = len(self.e_los),len(self.c_his)
        self.newstyle = 'SCORE' in [x.name for x in h[0][1].columns] ### changed from get_columns()

        def proc(ct,pname):
            """ Return the parameters from CALDB for conversion type
                ct, and alter shape to match internal representation.
            """
            ### t = np.reshape(h[ct][1].data.field(pname),[nc,ne])
            ### Note change for new or old format
            t = np.reshape(self.psf_data[ct].field(pname),[nc,ne])
            return t.transpose()[:,::-1]

        if self.newstyle:

            # SHAPE = N_CT X N_PAR X N_EN X N_INC_ANGLE
            self.tables = np.asarray([
                      [proc(ct,'NCORE'),
                       proc(ct,'NTAIL'),
                       proc(ct,'GCORE'),
                       proc(ct,'GTAIL'),
                       proc(ct,'SCORE'),
                       proc(ct,'STAIL')
                      ] for ct in [0,1]])
            self._norm_pars()

        else:
            self.tables = np.asarray([
                     [proc(ct,'GCORE'),
                      proc(ct,'SIGMA')
                     ] for ct in [0,1]])

    def _norm_pars(self):
        if hasattr(self,'_normalized'):
            if self['_normalized']: return
        self._normalized = True # attempt to ensure only called once

        # recall table shape is N_CT X N_PAR X N_EN X N_INC

        # normalize according to the irfs/latResponse conventions
        ens = (self.e_los*self.e_his)**0.5
        for ct in [0,1]: # iterate through conversion type
            scale_func = self.scale_func[ct]
            for i in xrange(self.tables.shape[2]): # iterate through energy
                sf = scale_func(ens[i])
                # vector operations in incidence angle
                nc,nt,gc,gt,sc,st = self.tables[ct,:,i,:]
                normc = self.psf_base_integral(gc,sc*sf,np.pi/2)[0]
                normt = self.psf_base_integral(gt,st*sf,np.pi/2)[0]
                # NB leave scale factor out here so we can adjust norm to
                # a particular energy (effectively cancelling in integral)
                norm = (TWOPI*(normc*sc**2+normt*nt*st**2))**-1
                self.tables[ct,0,i,:] = norm # adjust NCORE
                self.tables[ct,1,i,:]*= norm # set to NTAIL*NCORE

    def __call__(self,e,ct,delta, scale_sigma=True, density=True):
        """ Return the psf at the specified energy, for events
            with conversion type ct, at a remove of delta RADIANS.

            scale_sigma -- apply the energy prescaling
            density -- divide by "jacobian" to turn into photon density
                       NB -- this can only be TRUE now
        """
        if density != True:
            raise DeprecationWarning('This function only returns PSF density now.')
            density = True
        sf = self.scale_func[ct](e) if scale_sigma else 1
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
            yc = self.psf_base(gc,sc*sf,delta)
            yt = self.psf_base(gt,st*sf,delta)
            return (w * (nc*yc + nt*yt)/sf**2).sum(axis=1)
        else:
            # note retaining explicit old style for the nonce, for comparison testing
            gc,si,w = self.get_p(e,ct) # note each parameter is an 8-vector
            si *= sf
            us  = 0.5 * np.outer(delta,1./si)**2
            y    = (1-1./gc)*(1+us/gc)**(-gc)
            return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def psf_base(self,g,s,delta):
        """Implement the PSF base function; g = gamma, s = sigma (scaled), 
           delta = deviation in radians.
           
            Operation is vectorized both in parameters and angles; return
            is an array of shape (n_params,n_angles)."""
        u = 0.5 * np.outer(delta,1./s)**2
        return (1-1./g)*(1+u/g)**(-g)

    def psf_base_integral(self,g,s,dmax,dmin=0):
        """Integral of the PSF base function; g = gamma, s = sigma (scaled),
           delta = deviation in radians."""
        u0 = 0.5 * np.outer(dmin,1./s)**2
        u1 = 0.5 * np.outer(dmax,1./s)**2
        return (1+u0/g)**(1-g)-(1+u1/g)**(1-g)

    def integral(self,e,ct,dmax,dmin=0):
        """ Return integral of PSF at given energy (e) and conversion type (ct)
            from dmin to dmax (radians).
            
            Note -- does *not* take a vector argument right now."""
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
            sf = self.scale_func[ct](e)
            icore = TWOPI*sc**2*nc*\
                    self.psf_base_integral(gc,sc*sf,dmax,dmin=dmin)[0]
            itail = TWOPI*st**2*nt*\
                    self.psf_base_integral(gt,st*sf,dmax,dmin=dmin)[0]
            return (w*(icore+itail)).sum(axis=-1)
        else:
            gc,si,w = self.get_p(e,ct)
            si *= self.scale_func[ct](e)
            u1 = 0.5 * (dmin / si)**2
            u2 = 0.5 * (dmax / si)**2
            return  (w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  )).sum()

    def inverse_integral(self,e,ct,percent=68,on_axis=False):
        """Return the radius at which the integral PSF reaches the required percentage.

           NB -- returned value is in degrees.
        """
        percent = float(percent)/100
        if on_axis:
            sf = self.scale_func[ct](e)
            if self.newstyle:
                nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
                u1 = gc[0]*( (1-percent)**(1./(1-gc[0])) - 1)
                u2 = gt[0]*( (1-percent)**(1./(1-gt[0])) - 1)
                return RAD2DEG*sf*(nc[0]*(u1*2)**0.5*sc[0] + nt[0]*(u2*2)**0.5*st[0]) #approx
            else:
                gc,sc,w = self.get_p(e,ct)
                u1 = gc[0]*( (1-percent)**(1./(1-gc[0])) - 1)
                return RAD2DEG*sf*(2*u1)**0.5*sc[0]
        f = lambda x: abs(self.integral(e,ct,x) - percent)
        from scipy.optimize import fmin
        seeds = np.asarray([5,4,3,2.5,2,1.5,1,0.5,0.25])*self.scale_func[ct](e)
        seedvals = np.asarray([self.integral(e,ct,x) for x in seeds])
        seed = seeds[np.argmin(np.abs(seedvals-percent))]
        trial = fmin(f,seed,disp=0,ftol=0.000001,xtol=0.01)
        if trial > 0:
            return trial[0]*RAD2DEG
        print 'Warning: could not invert integral; return best grid value.'
        return seed*RAD2DEG

    def band_psf(self,band,weightfunc=None,adjust_mean=False):
        return BandCALDBPsf(self,band,weightfunc=weightfunc,
            adjust_mean=adjust_mean,newstyle=self.newstyle)

    def get_cpp_psf(self,e,ct):
        sf = self.scale_func[ct](e)
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
            return PythonPsf(sc*sf,st*sf,gc,gt,nc/sf**2,nt/sf**2,w)
        else:
            gc,si,w = self.get_p(e,ct)
            return PythonPsf(si*sf,gc,w)
        

class BandPsf(object):

    def init(self):
        self.newstyle = False
        self.override_en = None
        self.adjust_mean = False
        self.weightfunc = None

    def __init__(self,psf,band,**kwargs):
        self.init()
        self.parent_psf = psf # save a reference
        self.newstyle = psf.newstyle
        self.__dict__.update(kwargs)
        self.par = psf.get_p(self.override_en or band.e,band.ct).copy()

        index = 2
        if self.adjust_mean:
            # estimate an optimal energy based on the log likelihood
            from scipy.integrate import simps
            f = lambda e: e**-index
            b = band    
            rad = min(b.radius_in_rad,np.radians(psf.inverse_integral(band.e,band.ct,percent=99,on_axis=True)))
            dom = np.linspace(0,rad**2,101)**0.5

            # calculate the actual PSF integral over sub-bands
            fopt = np.zeros_like(dom)
            for i in xrange(len(b.sp_points)):
                n = psf(b.sp_points[i],b.ct,dom,density=True)*b.sp_vector[i]*b.sp_points[i]**-index
                fopt += n
            n = (b.sp_vector*f(b.sp_points)).sum()

            # an estimate of log likelihood in infinite-N limit
            logl1 = simps(fopt) - simps(fopt*np.log(fopt))

            # now, find E_opt that gives same likelihood; linear fit
            edom = np.linspace(band.e*0.98,band.e*1.02,3)
            ecod = np.empty_like(edom)
            for ien,en in enumerate(edom):
                fe = n*psf(en,b.ct,dom,density=True)
                ecod[ien] = simps(fe) - simps(fe*np.log(fe))
            slope = (ecod[-1] - ecod[0])/(band.e*0.04)
            self.eopt = band.e + (logl1 - ecod[1])/slope
            if abs(self.eopt/band.e - 1) > 0.10: # sanity check
                self.eopt = band.e
            self.scale = psf.scale_func[band.ct](self.eopt)

        elif self.weightfunc is not None:
            weightfunc = self.weightfunc
            dom    = np.logspace(np.log10(band.emin),np.log10(band.emax),9)
            wvals = [weightfunc(x) for x in dom]
            svals = psf.scale_func[band.ct](dom)
            num = simps(wvals*svals*dom,x=dom)
            dom = simps(wvals*dom,x=dom)
            self.scale = num/dom
            self.comp_scale = psf.scale_func[band.ct](band.e)
            print num,dom,num/dom,self.comp_scale

        else:
            self.scale = psf.scale_func[band.ct](self.override_en or band.e)
        # scale sigma parameters
        if self.newstyle:
            self.par[4] *= self.scale
            self.par[5] *= self.scale
            # also correct normalization
            self.par[0] /= self.scale**2
            self.par[1] /= self.scale**2
        else:
            self.par[1] *= self.scale
        if self.newstyle:
            pass
        else:
            g,s,w = self.par
            self.int_par = w,(2*s**2*g)**-1,1-g
        self.cpsf = self.get_cpp_psf()

    def get_cpp_psf(self):
        from skymaps import PythonPsf
        if self.newstyle:
            # NB sigma parameters are prescaled by scale factor
            nc,nt,gc,gt,sc,st,w = self.par
            return PythonPsf(sc,st,gc,gt,nc,nt,w)
        else:
            gc,si,w = self.par
            return PythonPsf(si,gc,w)

    def inverse_integral_on_axis(self,frac=0.68):
        """Return an approxmation of the angle that contains frac of the
           photon.  The approximation uses the on-axis PSF and, for newstyle
           PSFs, weights the two results according to ncore and ntail.

           Return value in degrees.
        """
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.par
            u1 = gc[0]*( (1-frac)**(1./(1-gc[0])) - 1)
            u2 = gt[0]*( (1-frac)**(1./(1-gt[0])) - 1)
            return RAD2DEG*(nc[0]*(u1*2)**0.5*sc[0] + nt[0]*(u2*2)**0.5*st[0])
        else:
            gc,sc,w = self.par
            u1 = gc[0]*( (1-frac)**(1./(1-gc[0])) - 1)
            return RAD2DEG*(2*u1)**0.5*sc[0]

class BandCALDBPsf(BandPsf,CALDBPsf):

    def __call__(self,delta,density=True):
        if density != True:
            raise DeprecationWarning('This function only returns PSF density now.')
            density = True
        if self.newstyle:
            # NB, we've already corrected sigmas with scale function
            nc,nt,gc,gt,sc,st,w = self.par
            yc = self.psf_base(gc,sc,delta)
            yt = self.psf_base(gt,st,delta)
            return (w*(nc*yc + nt*yt)).sum(axis=1)
        else:
            gc,si,w = self.par
            us  = 0.5 * np.outer(delta,1./si)**2
            y    = (1-1./gc)*(1+us/gc)**(-gc)
            return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def integral(self,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.par
            icore = TWOPI*sc**2*nc*\
                    self.psf_base_integral(gc,sc,dmax,dmin=dmin)
            itail = TWOPI*st**2*nt*\
                    self.psf_base_integral(gc,sc,dmax,dmin=dmin)
            return (w*(icore+itail)).sum()

        else:
            w,a,b = self.int_par
            return ( w*( (1+a*(dmin*dmin))**b - (1+a*(dmax*dmax))**b ) ).sum()

class PsfOverlap(object):
    """Routines to calculate how much of the emission of a point source falls onto an ROI."""

    def init(self):
        self.quadrature_tol = 1e-4
        self.cache_hash = None

    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    def set_dir_cache(self,band,roi_dir,radius):
        self.cache_wsdl = WeightedSkyDirList(band.b,roi_dir,radius,True)
        self.cache_diffs = np.empty(len(self.cache_wsdl))
        self.cache_hash = hash(band)

    def num_overlap(self,band,roi_dir,ps_dir,radius_in_rad=None,override_pdf=None):
        roi_rad  = radius_in_rad or band.radius_in_rad
        if self.cache_hash != hash(band):
            # fragile due to radius dep.
            self.set_dir_cache(band,roi_dir,roi_rad) 
        if override_pdf is None:
            band.psf.cpsf.wsdl_val(self.cache_diffs,ps_dir,self.cache_wsdl)
        else:
            difference = np.empty(len(self.cache_wsdl))
            PythonUtilities.arclength(difference,self.cache_wsdl,roi_dir)
            self.cache_diffs = override_pdf(difference)
        return self.cache_diffs.sum()*band.b.pixelArea()

    def __call__(self,band,roi_dir,ps_dir,radius_in_rad=None,
                 ragged_edge=0.06,
                 override_pdf=None,override_integral=None):
        """ Return the fractional overlap for a point source at location 
            skydir.

            band           an ROIBand object
            roi_dir        the SkyDir giving the center of the ROI
            ps_dir         the SkyDir giving the point source position
            radius_in_rad  the ROI radius (radians)
            ragged_edge    tolerance for assumption that HEALPix
                               representation is circular
            override_pdf   an alternative psf density to use
            override_integral    an alternative cumulative psf to use
        """

        roi_rad  = band.radius_in_rad if radius_in_rad is None \
                    else radius_in_rad
        integral = band.psf.integral if override_integral is None \
                    else override_integral
        offset = roi_dir.difference(ps_dir)
        tol = self.quadrature_tol

        # if point source is close to center of ROI, use PSF integral
        if offset < 1e-5:
            overlap = integral(roi_rad) 

        # if point source within ROI
        elif offset < roi_rad:
            def interior(x):
                c = cos(x)
                s2 = c**2-1
                effrad = (roi_rad**2 + s2*offset**2)**0.5 - offset*c
                return integral(effrad)
            overlap = quad(interior,0,np.pi,epsabs=tol)[0]/np.pi

        # for point sources outside of ROI
        else:
            def exterior(x):
                c = cos(x)
                s2 = c**2-1
                r2 = (roi_rad**2 + s2*offset**2)**0.5
                de = offset*c - r2
                return integral(de+2*r2,dmin=de)
            limit = np.arcsin(roi_rad / offset)
            overlap = quad(exterior,0,limit,epsabs=tol)[0]/np.pi

        # check to see if ROI is circular enough for above calculation
        # or if we need instead to do a cubature with HEALPix
        if (((band.b.pixelArea()**0.5/band.radius_in_rad) > ragged_edge) 
            and (band.b.nside() < 200)):
            n_overlap = self.num_overlap(band,roi_dir,ps_dir,roi_rad,override_pdf)
            return n_overlap

        return overlap

class CPsfOverlap(PsfOverlap):
    """ Same as PsfOverlap, but using C++ implementation of the
        overlap integral."""
    def __call__(self,band,roi_dir,ps_dir,radius_in_rad=None,
                 ragged_edge=0.06):
        """ Return the fractional overlap for a point source at location 
            skydir.

            band           an ROIBand object
            roi_dir        the SkyDir giving the center of the ROI
            ps_dir         the SkyDir giving the point source position
            radius_in_rad  the ROI radius (radians)
            ragged_edge    tolerance for assumption that HEALPix
                               representation is circular
        """

        roi_rad  = band.radius_in_rad if radius_in_rad is None \
                    else radius_in_rad
        overlap = band.psf.cpsf.overlap_circle(roi_dir,roi_rad,ps_dir)

        # check to see if ROI is circular enough for above calculation
        # or if we need instead to do a cubature with HEALPix
        if (((band.b.pixelArea()**0.5/band.radius_in_rad) > ragged_edge) 
            and (band.b.nside() < 200)):
            overlap = self.num_overlap(band,roi_dir,ps_dir,roi_rad)

        return overlap

class PsfOverlapHealpix(object):
    """Routines to calculate how much of the emission of a point source falls onto an ROI."""

    def init(self):
        self.quadrature_tol = 1e-4

    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    def __call__(self,band,nside,index,ps_dir):
        """Return an array of fractional overlap for a point source at location skydir."""

        # pick an nside commensurate with PSF and that will fit evenly within base pixel
        scale = band.psf.parent_psf.scale_func[band.ct](band.e)/5
        sample_nside = (np.pi/3)**0.5/scale
        p0 = np.log(sample_nside/nside)/np.log(2)
        last_overlap = None
        overlaps = []
        for addon in [0]:
            sample_nside = min(8192,nside*2**(int(p0)+addon))
            geo = (4*np.pi/12/sample_nside**2)
            #print sample_nside,index

            # temporary code for testing!
            if True:
                from skymaps import Band
                band.b = Band(sample_nside)

            wsdl  = WeightedSkyDirList(band.b,nside,index,True)
            vals  = np.empty(len(wsdl))
            band.psf.cpsf.wsdl_val(vals,ps_dir,wsdl)
            overlap = vals.sum()*geo
            #print overlap
            overlaps.append(overlap)
            #if last_overlap is not None:
            #    new_est = (2*last_overlap+4*overlap)/6
            #    print 'Updated estimate: ',new_est
            last_overlap = overlap
        #print ('%.6f\t'*(len(overlaps)))%(tuple(overlaps))
        return overlaps[-1],sample_nside

class PretendBand(object):

     def __init__(self,energy,conversion_type,**kwargs):
          self.e = energy; self.ct = conversion_type
          self.__dict__.update(kwargs)

class ConvolutionPsf(object):
     """N.B. -- the PSF center is assumed to be at the Galactic north pole."""

     def __init__(self,CALDB,irf='P6_v3_diff',psf_irf=None):

          self.caldb_psf = CALDBPsf(CALDB,irf=irf,psf_irf=psf_irf)

     def set_band(self,energy,conversion_type):
          band = PretendBand(energy,conversion_type)
          self.psf = BandCALDBPsf(self.caldb_psf,band)

     def __call__(self,v):
          sd = SkyDir(Hep3Vector(v[0],v[1],v[2]))
          return self.psf(DEG2RAD*(90. - sd.b()))[0]          

     def get_pyskyfun(self):
          return PySkyFunction(self)

