"""
A module to manage the PSF from CALDB and handle the integration over
incidence angle and intepolation in energy required for the binned
spectral analysis.
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pypsf.py,v 1.18 2010/09/01 21:35:16 mar0 Exp $
author: M. Kerr

"""

import pyfits as pf
import numpy as N
from os.path import join
import os
from cPickle import load
from skymaps import ExposureWeighter,SkyDir,PySkyFunction,Hep3Vector,WeightedSkyDirList,PythonUtilities
from scipy.integrate import quad,simps
from math import cos,sin

def scale_factor(e,c0,c1,exp):
    return ( (c0*(e/100)**(exp))**2 + c1**2 )**0.5

TWOPI = (2*N.pi)
RAD2DEG = 180./N.pi

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class Psf(object):

    def init(self):
        self.irf      = 'P6_v3_diff'
        self.psf_irf = None # a possible override to use an on-orbit PSF file

    def __init__(self,CALDB=None,**kwargs):
        if CALDB is None:
            try: CALDB =os.environ['CALDB']
            except:
                 print 'Environment variable CALDB was not set'
                 raise
        self.CALDB = join(CALDB,'bcf')
        if not os.path.exists(self.CALDB):
            self.CALDB = join(CALDB,'data', 'glast', 'lat', 'bcf')
            if not os.path.exists(self.CALDB):
                raise Exception,  'Invalid CALDB'
        self.init()
        self.__dict__.update(kwargs)

        self.__readCALDB__()
        self.__read_params__()
        self.__calc_weights__()

    def __read_params__(self):pass

    def __readCALDB__(self):

        # open CALDB files
        if self.psf_irf is not None:
            print 'Overriding default PSF; using %s'%(self.psf_irf)
        psf_files = [join(self.CALDB,'psf','psf_%s_%s.fits'%(self.psf_irf or self.irf,x)) for x in ['front','back']]

        try:
            h0,h1 = self.CALDBhandles = [pf.open(x) for x in psf_files]
        except:
            raise Exception,'Could not open CALDB files for PSF!  Aborting!'

        # read in stuff that doesn't depend on conversion type
        self.scale_factors = N.asarray(h0[2].data.field('PSFSCALE')).flatten()
        self.e_los            = N.asarray(h0[1].data.field('ENERG_LO')).flatten().astype(float)
        self.e_his            = N.asarray(h0[1].data.field('ENERG_HI')).flatten().astype(float)
        self.c_los            = N.asarray(h0[1].data.field('CTHETA_LO')).flatten().astype(float)
        self.c_his            = N.asarray(h0[1].data.field('CTHETA_HI')).flatten().astype(float)

        sf = self.scale_factors
        self.scale_func = [lambda e: ( (sf[0]*(e/100)**(sf[-1]))**2 + sf[1]**2 )**0.5,
                                 lambda e: ( (sf[2]*(e/100)**(sf[-1]))**2 + sf[3]**2 )**0.5 ]


    def __calc_weights__(self,livetimefile='',skydir=None):

        aeffstrs = [join(self.CALDB,'ea','aeff_%s_%s.fits'%(self.irf,ct)) for ct in ['front','back']]
        ew         = ExposureWeighter(aeffstrs[0],aeffstrs[1],livetimefile)
        dummy     = skydir or SkyDir() 

        elo,ehi,clo,chi = self.e_los,self.e_his,self.c_los[::-1],self.c_his[::-1]

        weights  = N.zeros([2,len(elo),len(clo)])
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
        ind    = min(N.searchsorted(self.e_his,e),len(self.e_his) - 1)
        p      = self.tables[ct,:,ind,:]
        w      = self.weights[ct,ind,:]
        if cthetabin is None:
            return N.append(p,[w],axis=0)
        else:
            return N.append(p[:,cthetabin],1)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class CALDBPsf(Psf):
    """Handle reading in the parameters from CALDB for multiple formats and implement the PSF functions."""

    def __read_params__(self):        

        h          = self.CALDBhandles # a tuple with handles to front and back CALDB
        ne,nc     = len(self.e_los),len(self.c_his)
        self.newstyle = 'SCORE' in [x.name for x in h[0][1].get_coldefs()]

        if self.newstyle:

            tables = [[1./(1+N.reshape(h[i][1].data.field('NTAIL'),[nc,ne]).transpose()[:,::-1]),
                          N.reshape(h[i][1].data.field('NTAIL'),[nc,ne]).transpose()[:,::-1]/(1+N.reshape(h[i][1].data.field('NTAIL'),[nc,ne]).transpose()[:,::-1]),
                          N.reshape(h[i][1].data.field('GCORE'),[nc,ne]).transpose()[:,::-1],
                          N.reshape(h[i][1].data.field('GTAIL'),[nc,ne]).transpose()[:,::-1],
                          N.reshape(h[i][1].data.field('SCORE'),[nc,ne]).transpose()[:,::-1],
                          N.reshape(h[i][1].data.field('STAIL'),[nc,ne]).transpose()[:,::-1]] for i in [0,1]]

        else:
            tables = [[N.reshape(h[i][1].data.field('GCORE'),[nc,ne]).transpose()[:,::-1],
                          N.reshape(h[i][1].data.field('SIGMA'),[nc,ne]).transpose()[:,::-1]] for i in [0,1]]

        self.tables  = N.asarray(tables)

    def __call__(self,e,ct,delta, scale_sigma=True, density = False):
        sf = self.scale_func[ct](e) if scale_sigma else 1
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
            yc = self.psf_base(gc,sc*sf,delta,density)
            yt = self.psf_base(gt,st*sf,delta,density)
            return (w * (nc*yc + nt*yt)).sum(axis=1)
        else:
            # note retaining explicit old style for the nonce, for comparison testing
            gc,si,w = self.get_p(e,ct) # note each parameter is an 8-vector
            si *= sf
            us  = 0.5 * N.outer(delta,1./si)**2
            y    = (1-1./gc)*(1+us/gc)**(-gc)
            return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def psf_base(self,g,s,delta,density=False):
        """Implement the PSF base function; g = gamma, s = sigma, w = weight, delta = deviation in radians."""
        u = 0.5 * N.outer(delta,1./s)**2
        y = (1-1./g)*(1+u/g)**(-g)
        return y / (TWOPI*s**2 if density else 1.)

    """ Not sure these are in use -- if they are, update for old/new CALDB formats.
    def set_cache(self,e,ct):
        np = len(self.get_p(e[0],ct[0])[0])
        self.cache = N.empty([3,len(e),np])
        ca = self.cache; gp = self.get_p; sf = self.scale_func
        for i,(mye,myct) in enumerate(zip(e,ct)):
            ca[:,i,:]  = gp(mye,myct)
            ca[1,i,:] *= sf[myct](mye)

    def get_cache(self,delta,density = False):
        g,s,w = self.cache
        us     = 0.5 * (delta / s)**2
        y      = y  = (1-1./g)*(1+us/g)**(-g)
        return (w*(y / (TWOPI*s**2 if density else 1.))).sum(axis=1)
    """

    def integral(self,e,ct,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        """ Update to density-style."""
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.get_p(e,ct)
            #if scale_sigma: # why no test in old version? think this code may have gotten stale
            sf = self.scale_func[ct](e)
            sc *= sf
            st *= sf
            uc1 = 0.5 * (dmin / sc)**2
            uc2 = 0.5 * (dmax / sc)**2
            ut1 = 0.5 * (dmin / st)**2
            ut2 = 0.5 * (dmax / st)**2
            return (w*( nc*((1+uc1/gc)**(1-gc) - (1+uc2/gc)**(1-gc)) + nt*((1+ut1/gt)**(1-gt) - (1+ut2/gt)**(1-gt)) )).sum()
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
        seeds = N.asarray([5,4,3,2.5,2,1.5,1,0.5,0.25])*self.scale_func[ct](e)
        seedvals = N.asarray([self.integral(e,ct,x) for x in seeds])
        seed = seeds[N.argmin(N.abs(seedvals-percent))]
        trial = fmin(f,seed,disp=0,ftol=0.01,xtol=0.01)
        if trial > 0:
            return trial[0]*RAD2DEG
        print 'Warning: could not invert integral; return best grid value.'
        return seed*RAD2DEG

    def band_psf(self,band,weightfunc=None,adjust_mean=False):
        return BandCALDBPsf(self,band,weightfunc=weightfunc,adjust_mean=adjust_mean,newstyle=self.newstyle)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class BandPsf(object):

    def init(self):
        self.newstyle    = False
        self.override_en = None
        self.adjust_mean = False
        self.weightfunc  = None

    def __init__(self,psf,band,**kwargs):
        self.init()
        self.parent_psf = psf # save a reference
        self.newstyle = psf.newstyle
        self.__dict__.update(kwargs)
        self.par     = psf.get_p(self.override_en or band.e,band.ct).copy()

        index = 2
        if self.adjust_mean:
            # estimate an optimal energy based on the log likelihood
            from scipy.integrate import simps
            f = lambda e: e**-index
            b = band    
            rad = min(b.radius_in_rad,N.radians(psf.inverse_integral(band.e,band.ct,percent=99,on_axis=True)))
            dom = N.linspace(0,rad**2,101)**0.5

            # calculate the actual PSF integral over sub-bands
            fopt = N.zeros_like(dom)
            for i in xrange(len(b.sp_points)):
                n = psf(b.sp_points[i],b.ct,dom,density=True)*b.sp_vector[i]*b.sp_points[i]**-index
                fopt += n
            n = (b.sp_vector*f(b.sp_points)).sum()

            # an estimate of log likelihood in infinite-N limit
            logl1 = simps(fopt) - simps(fopt*N.log(fopt))

            # now, find E_opt that gives same likelihood; linear fit
            edom = N.linspace(band.e*0.98,band.e*1.02,3)
            ecod = N.empty_like(edom)
            for ien,en in enumerate(edom):
                fe = n*psf(en,b.ct,dom,density=True)
                ecod[ien] = simps(fe) - simps(fe*N.log(fe))
            slope = (ecod[-1] - ecod[0])/(band.e*0.04)
            self.eopt = band.e + (logl1 - ecod[1])/slope
            if abs(self.eopt/band.e - 1) > 0.10: # sanity check
                self.eopt = band.e
            self.scale = psf.scale_func[band.ct](self.eopt)

        elif self.weightfunc is not None:
            weightfunc = self.weightfunc
            dom    = N.logspace(N.log10(band.emin),N.log10(band.emax),9)
            wvals = [weightfunc(x) for x in dom]
            svals = psf.scale_func[band.ct](dom)
            num = simps(wvals*svals*dom,x=dom)
            dom = simps(wvals*dom,x=dom)
            self.scale = num/dom
            self.comp_scale = psf.scale_func[band.ct](band.e)
            print num,dom,num/dom,self.comp_scale

        else:
            self.scale    = psf.scale_func[band.ct](self.override_en or band.e)
        if self.newstyle:
            self.par[4] *= self.scale
            self.par[5] *= self.scale
        else:
            self.par[1] *= self.scale # scale sigma
        if self.newstyle:
            pass
        else:
            g,s,w = self.par
            self.int_par = w,(2*s**2*g)**-1,1-g
        self.cpsf = self.get_cpp_psf()

    def get_cpp_psf(self):
        from skymaps import PythonPsf
        if self.newstyle:
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

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
class BandCALDBPsf(BandPsf):

    def psf_base(self,g,s,delta,density=False):
        """Implement the PSF base function; g = gamma, s = sigma, delta = deviation in radians."""
        u = 0.5 * N.outer(delta,1./s)**2
        y = (1-1./g)*(1+u/g)**(-g)
        return y / (TWOPI*s**2 if density else 1.)

    def __call__(self,delta,density=False):
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.par
            yc = self.psf_base(gc,sc,delta,density)
            yt = self.psf_base(gt,st,delta,density)
            return (w*(nc*yc + nt*yt)).sum(axis=1)
        else:
            gc,si,w = self.par
            us  = 0.5 * N.outer(delta,1./si)**2
            y    = (1-1./gc)*(1+us/gc)**(-gc)
            return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def integral(self,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        if self.newstyle:
            nc,nt,gc,gt,sc,st,w = self.par
            uc1 = 0.5 * (dmin / sc)**2
            uc2 = 0.5 * (dmax / sc)**2
            ut1 = 0.5 * (dmin / st)**2
            ut2 = 0.5 * (dmax / st)**2
            return (w*( nc*((1+uc1/gc)**(1-gc) - (1+uc2/gc)**(1-gc)) + nt*((1+ut1/gt)**(1-gt) - (1+ut2/gt)**(1-gt)) )).sum()
        else:
            #gc,si,w = self.par
            #u1 = 0.5 * (dmin / si)**2
            #u2 = 0.5 * (dmax / si)**2
            #return  ( w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)) ).sum()
            w,a,b = self.int_par
            return ( w*( (1+a*(dmin*dmin))**b - (1+a*(dmax*dmax))**b ) ).sum()

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
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
        self.cache_diffs = N.empty(len(self.cache_wsdl))
        self.cache_hash = hash(band)

    def num_overlap(self,band,roi_dir,ps_dir,radius_in_rad=None,override_pdf=None):
        roi_rad  = radius_in_rad or band.radius_in_rad
        if self.cache_hash != hash(band): self.set_dir_cache(band,roi_dir,roi_rad) # fragile due to radius dep.
        if override_pdf is None:
            band.psf.cpsf.wsdl_val(self.cache_diffs,ps_dir,self.cache_wsdl)
        else:
            difference = N.empty(len(self.cache_wsdl))
            PythonUtilities.arclength(difference,self.cache_wsdl,roi_dir)
            self.cache_diffs = override_pdf(difference)
        return self.cache_diffs.sum()*band.b.pixelArea()

    def __call__(self,band,roi_dir,ps_dir,radius_in_rad=None,ragged_edge=0.06,
                 override_pdf=None,override_integral=None):
        """Return the fractional overlap for a point source at location skydir.
            Note radius arguments are in radians."""

        roi_rad  = radius_in_rad or band.radius_in_rad
        integral = band.psf.integral if override_integral is None else override_integral

        offset    = roi_dir.difference(ps_dir)

        if offset < 1e-5:
            overlap = integral(roi_rad) #point source in center of ROI, symmetry

        elif offset < roi_rad:

            def interior(x):

                c         = cos(x)
                eff_rad = ( roi_rad**2 + offset**2*(c**2 - 1) )**0.5 - offset*c
                return integral(eff_rad)

            overlap = quad(interior,0,N.pi,epsabs=self.quadrature_tol)[0]/N.pi

        else:

            def exterior(x):

                c     = cos(x)
                r2    = ( roi_rad**2 - (1-c**2)*offset**2 )**0.5
                de    = offset*c - r2
                return integral(dmax=de+2*r2,dmin=de)


            limit = N.arcsin(roi_rad / offset)
            overlap = quad(exterior,0,limit,epsabs=self.quadrature_tol)[0]/N.pi


        if ((band.b.pixelArea()**0.5/band.radius_in_rad) > ragged_edge) and (band.b.nside() < 200):
            n_overlap = self.num_overlap(band,roi_dir,ps_dir,roi_rad,override_pdf)
            #print overlap,n_overlap
            #print 'Using numerical overlap, difference is %.6f'%(n_overlap-overlap)
            return n_overlap

        return overlap

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###
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
        sample_nside = (N.pi/3)**0.5/scale
        p0 = N.log(sample_nside/nside)/N.log(2)
        last_overlap = None
        overlaps = []
        for addon in [0]:
            sample_nside = min(8192,nside*2**(int(p0)+addon))
            geo = (4*N.pi/12/sample_nside**2)
            #print sample_nside,index

            # temporary code for testing!
            if True:
                from skymaps import Band
                band.b = Band(sample_nside)

            wsdl = WeightedSkyDirList(band.b,nside,index,True)
            from pointlike import DoubleVector
            dv = DoubleVector()
            wsdl.arclength(ps_dir,dv)
            diffs = N.fromiter(dv,dtype=float)
            vals = band.psf(diffs,density=True)
            overlap = vals.sum()*geo
            #print overlap
            overlaps.append(overlap)
            #if last_overlap is not None:
            #    new_est = (2*last_overlap+4*overlap)/6
            #    print 'Updated estimate: ',new_est
            last_overlap = overlap
        #print ('%.6f\t'*(len(overlaps)))%(tuple(overlaps))
        return overlaps[-1],sample_nside


###====================================================================================================###
###deprecated==========================================================================================###
###====================================================================================================###

deg2rad = N.pi / 180.

class PretendBand(object):

     def __init__(self,energy,conversion_type):
          self.e = energy; self.ct = conversion_type

class ConvolutionPsf(object):
     """N.B. -- the PSF center is assumed to be at the Galactic north pole."""

     def __init__(self,CALDB,irf='P6_v3_diff',psf_irf=None):

          self.caldb_psf = CALDBPsf(CALDB,irf=irf,psf_irf=psf_irf)

     def set_band(self,energy,conversion_type):
          band = PretendBand(energy,conversion_type)
          self.psf = BandCALDBPsf(self.caldb_psf,band)

     def __call__(self,v):
          sd = SkyDir(Hep3Vector(v[0],v[1],v[2]))
          return self.psf(deg2rad*(90. - sd.b()))[0]          

     def get_pyskyfun(self):
          return PySkyFunction(self)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class NewPsf(Psf):

    def __init__(self,CALDB,my_CALDB_pickle,**kwargs):
        self.my_CALDB_pickle = my_CALDB_pickle
        super(NewPsf,self).__init__(CALDB,**kwargs)        

    def __read_params__(self):
        self.p = p = load(file(self.my_CALDB_pickle))

        # tables are in order: NCORE, SIGMA, A, GCORE, GTAIL
        ftables = N.asarray(p['ftables'])[:,:,:-1] # trim of 0-66 deg fit
        btables = N.asarray(p['btables'])[:,:,:-1] 

        # do not fit 24 MeV bin
        self.e_los     = self.e_los[1:]
        self.e_his     = self.e_his[1:]

        self.tables = N.asarray([ftables,btables])

    def __call__(self,e,ct,delta, scale_sigma=True, density=False,cthetabin=None):
        """Return the differential psf at given energy and conversion type.

            Parameters
            ----------
            e : float, the energy in MeV
            ct : float, the conversion type (0 = front, 1 = back)
            delta : ndarray(float), the angular separation in radians

            Returns
            -------
            output : ndarray (1d), the differential psf evaluated at the required points
        """
        nc,si,ap,gc,gt,w = self.get_p(e,ct,cthetabin) # note each parameter is an 8-vector
        if scale_sigma: si *= self.scale_func[ct](e)
        us  = 0.5 * N.outer(delta,1./si)**2
        y    = (nc*(1-1./gc))*(1+us/gc)**(-gc) + ((1-nc)*ap**(gt-1))*(1.-1./gt)*(ap+us/gt)**(-gt)
        return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def integral(self,e,ct,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        nc,si,ap,gc,gt,w = self.get_p(e,ct)
        si *= self.scale_func[ct](e)
        u1 = 0.5 * (dmin / si)**2
        u2 = 0.5 * (dmax / si)**2
        return  (w*(nc*              ( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  ) + \
                  (1-nc)*ap**(gt-1)*( (ap+u1/gt)**(1-gt) - (ap+u2/gt)**(1-gt) ) )).sum()

    def band_psf(self,band,weightfunc=None): return BandNewPsf(self,band,weightfunc)

###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class OldPsf(Psf):

    def __read_params__(self):        
        h      = self.CALDBhandles # a tuple with handles to front and back CALDB
        ne,nc = len(self.e_los),len(self.c_his)
        tables = [[N.reshape(h[i][1].data.field('GCORE'),[nc,ne]).transpose()[:,::-1],
                      N.reshape(h[i][1].data.field('SIGMA'),[nc,ne]).transpose()[:,::-1]] for i in [0,1]]
        self.tables  = N.asarray(tables)

    def __call__(self,e,ct,delta, scale_sigma=True, density = False):
        gc,si,w = self.get_p(e,ct) # note each parameter is an 8-vector
        if scale_sigma: si *= self.scale_func[ct](e)
        us  = 0.5 * N.outer(delta,1./si)**2
        y    = (1-1./gc)*(1+us/gc)**(-gc)
        return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def set_cache(self,e,ct):
        np = len(self.get_p(e[0],ct[0])[0])
        self.cache = N.empty([3,len(e),np])
        ca = self.cache; gp = self.get_p; sf = self.scale_func
        for i,(mye,myct) in enumerate(zip(e,ct)):
            ca[:,i,:]  = gp(mye,myct)
            ca[1,i,:] *= sf[myct](mye)

    def get_cache(self,delta,density = False):
        g,s,w = self.cache
        us     = 0.5 * (delta / s)**2
        y      = y  = (1-1./g)*(1+us/g)**(-g)
        return (w*(y / (TWOPI*s**2 if density else 1.))).sum(axis=1)

    def integral(self,e,ct,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        gc,si,w = self.get_p(e,ct)
        si *= self.scale_func[ct](e)
        u1 = 0.5 * (dmin / si)**2
        u2 = 0.5 * (dmax / si)**2
        return  (w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  )).sum()

    def band_psf(self,band,weightfunc=None): return BandOldPsf(self,band,weightfunc)


###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class BandNewPsf(BandPsf):

    def __call__(self,delta,density=False):
        nc,si,ap,gc,gt,w = self.par
        us  = 0.5 * N.outer(delta,1./si)**2
        y    = (nc*(1-1./gc))*(1+us/gc)**(-gc) + ((1-nc)*ap**(gt-1))*(1.-1./gt)*(ap+us/gt)**(-gt)
        return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def integral(self,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        nc,si,ap,gc,gt,w = self.par
        u1 = 0.5 * (dmin / si)**2
        u2 = 0.5 * (dmax / si)**2
        return  (w*(nc*              ( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)  ) + \
                  (1-nc)*ap**(gt-1)*( (ap+u1/gt)**(1-gt) - (ap+u2/gt)**(1-gt) ) )).sum()


###====================================================================================================###
###====================================================================================================###
###====================================================================================================###

# deprecated
class BandOldPsf(BandPsf):

    def __call__(self,delta,density=False):
        gc,si,w = self.par
        us  = 0.5 * N.outer(delta,1./si)**2
        y    = (1-1./gc)*(1+us/gc)**(-gc)
        return (w*(y / (TWOPI*si**2 if density else 1.))).sum(axis=1)

    def integral(self,dmax,dmin=0):
        """ Note -- does *not* take a vector argument right now."""
        gc,si,w = self.par
        u1 = 0.5 * (dmin / si)**2
        u2 = 0.5 * (dmax / si)**2
        return  ( w*( (1+u1/gc)**(1-gc)  - (1+u2/gc)**(1-gc)) ).sum()
