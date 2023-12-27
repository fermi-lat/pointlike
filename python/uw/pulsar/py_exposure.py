
"""
A module for calculating the exposure.  The main use case is for a calculation
at a single position on the sky, e.g., to perform a periodicity search when
the exposure varies appreciably over a cycle.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/py_exposure.py,v 1.1 2010/11/01 19:54:18 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>

"""

import numpy as N
from math import sin,cos
import astropy.io.fits as pf
import uw.utilities.fitstools as fitstools
from skymaps import Gti,Band,Healpix,SkyDir
from os.path import join

DEG2RAD = N.pi/180.

# TODO -- convert to the sqrt(1-cos(theta)) binning adopted in the ST

class Livetime(object):
    """Calculate the livetime as a function of incidence angle using the GTI
       specified in a collection of FT1 files and the livetime entries from
       an FT2 file.

       The default implementation is fully-unbinned, i.e., when the user
       requests the livetime, the exact values for the S/C z-axis and
       zenith positions are used to calculate the incidence/zenith angles.

       This executes with comparable speed (factor of ~2 slower) to the
       Science Tools application gtltcube."""

    def init(self):
        self.gti_mask  = None  # an additional mask (e.g., to remove a GRB)
        self.zenithcut = cos(DEG2RAD*105.) # 105 deg
        self.fovcut    = 0.4   # 66.4 deg
        self.nbins     = 20    # bins in cos(theta)
        self.tstart    = 0     # lower time limit in MET
        self.tstop     = 1e100 # upper time limit in MET
        self.verbose   = 1     # default verbosity
        self.mask_zero = True  # mask out bins with 0 livetime if True
        self.fast_ft2  = True  # use fast algorithm to process FT2/GTI
        self.cosbins   = N.linspace(self.fovcut,1,self.nbins+1)

    def finish(self):
        for field in self.fields:
            if ('DEC_' in field):
                self.__dict__['COS_'+field] = N.cos(self.__dict__[field])
                self.__dict__['SIN_'+field] = N.sin(self.__dict__[field])

    def __init__(self,ft2files,ft1files,**kwargs):
        self.init()
        self.__dict__.update(kwargs)
        self.prev_vals = self.prev_ra = self.prev_dec = None # initialize caching
        self.fields    = ['START','STOP','LIVETIME','RA_SCZ','DEC_SCZ','RA_ZENITH','DEC_ZENITH']

        if not hasattr(ft2files,'__len__'): ft2files = [ft2files]
        if not hasattr(ft1files,'__len__'): ft1files = [ft1files]

        self.__setup_gti(ft1files)
        self.__setup_ft2(ft2files)
        self.__process_ft2()

        for field in self.fields:
            if ('RA_' in field) or ('DEC_' in field):
                self.__dict__[field] *= DEG2RAD

        self.finish()

    def __setup_gti(self,ft1files):
        # the procedure is to take the union of all GTIs provided by FT1 files
        # then, take an intersection with the (optional) gti_mask and the time limits
        if self.verbose >= 1:
            print ('Processing GTI...')
        gti = self.gti = Gti(ft1files[0])
        if len(ft1files) > 1:
            for ft1 in ft1files[1:]: gti.combine(Gti(ft1))
        tmin = max(gti.minValue(),self.tstart)
        tmax = min(gti.maxValue(),self.tstop)
        gti = self.gti = gti.applyTimeRangeCut(tmin,tmax)
        if self.gti_mask is not None:
            before = round(gti.computeOntime())
            gti.interserction(self.gti_mask)
            if verbose >= 1:
                print ('Applied GTI mask; ontime reduced from %ds to %ds'%(before,round(gti.computerOntime())))

        ### NB -- this iteration takes way longer than it should -- put in an accessor method in SWIG
        self.gti_starts,self.gti_stops = \
                N.asarray([(x.minValue(),x.maxValue()) for x in gti]).transpose()

        N.sort(self.gti_starts);N.sort(self.gti_stops)
        if self.verbose >= 1:
            print ('Finished computing GTI from FT1 files; total ontime = %ds'%(round(gti.computeOntime())))

    def __setup_ft2(self,ft2files):
        """Load in the FT2 data.  Optionally, mask out values that will not
           contibute to the exposure."""
        if self.verbose >= 1:
            print ('Loading FT2 files...')
        handles = [pf.open(ft2,memmap=True) for ft2 in ft2files]
        ft2lens = [handle['SC_DATA'].data.shape[0] for handle in handles]
        fields  = self.fields
        arrays  = [N.empty(sum(ft2lens)) for i in range(len(fields))]
        
        counter = 0
        for ihandle,handle in enumerate(handles):
            if self.verbose > 1:
                print ('...Loading FT2 file # %d'%(ihandle))
            n = ft2lens[ihandle]
            for ifield,field in enumerate(fields):
                arrays[ifield][counter:counter+n] = handle['SC_DATA'].data.field(field)
            handle.close()
        for ifield,field in enumerate(fields):
            self.__dict__[field] = arrays[ifield]
        if self.mask_zero:
            if self.verbose >= 1: print ('Pruning values that yield 0 exposure...')
            mask = (self.LIVETIME > 0) | (self.STOP > self.gti_starts[0]) | (self.START < self.gti_stops[-1])
            self.mask_entries()

        self.mask_entries(N.argsort(self.START)) # sort the FT2 file in case it isn't
        if self.verbose > 1:
            print ('Finished loading FT2 files!')
  
    def __process_ft2(self):
        if self.verbose >= 1: print ('Processing the FT2 file (calculating overlap with GTI)...')
        if self.fast_ft2: overlaps = self.__process_ft2_fast(self.gti_starts,self.gti_stops)
        else:             overlaps = self.__process_ft2_slow(self.gti_starts,self.gti_stops)
        mask = overlaps > 0
        if self.mask_zero: self.mask_entries(mask)
        self.LTFRAC    = self.LIVETIME/(self.STOP-self.START)
        self.fields += ['LTFRAC']
        self.LIVETIME *= overlaps[mask]
        if self.verbose > 1: print ('Finished processing the FT2 file!')

    def __process_ft2_slow(self,gti_starts,gti_stops):
        """Calculate the fraction of each FT2 interval lying within the GTI.
           Uses a slow, easily-checked algorithm.
           The complexity is O(t^2) and is prohibitive
           for mission-length files."""
        t1,t2,lt = self.START,self.STOP,self.LIVETIME
        overlaps = N.zeros_like(lt)
        for i,(gti_t1,gti_t2) in enumerate(zip(gti_starts,gti_stops)):
            maxi = N.maximum(gti_t1,t1)
            mini = N.minimum(gti_t2,t2)
            overlaps += N.maximum(0,mini - maxi)
        return overlaps/(t2 - t1)
        
    def __process_ft2_fast(self,gti_starts,gti_stops):
        """Calculate the fraction of each FT2 interval lying within the GTI.
           Use binary search to quickly process the FT2 file.
           The complexity is O(t*log(t))."""
        t1,t2,lt = self.START,self.STOP,self.LIVETIME
        gti_t1,gti_t2 = gti_starts,gti_stops
        overlaps = N.zeros_like(lt)
        i1 = N.searchsorted(t2,gti_t1) # NB -- t2 in both not a typo
        i2 = N.searchsorted(t2,gti_t2)
        seps = i2 - i1
        for i,(ii1,ii2) in enumerate(zip(i1,i2)):
            overlaps[ii1:ii2+1] = 1. # fully-contained FT2 intervals
            if seps[i] > 0: # correct the endpoint FT2 intervals
                overlaps[ii1] = (t2[ii1] - gti_t1[i])/(t2[ii1] - t1[ii1])
                overlaps[ii2] = (gti_t2[i] - t1[ii2])/(t2[ii2] - t1[ii2])
            else: # edge case with exceptionally short GTI
                a = max(t1[ii1],gti_t1[i])
                b = min(t2[ii1],gti_t2[i])
                overlaps[ii1] = max(b-a,0)/(t2[ii1] - t1[ii1])
        return overlaps

    def mask_entries(self,mask=None):
        """If make is None, assume a LIVETIME > 0 cut."""
        if mask is None: mask = self.LIVETIME > 0
        for field in self.fields:
            self.__dict__[field] = self.__dict__[field][mask]

    def get_cosines(self,skydir):
        """Return the cosine of the arclength between the specified direction
           and the S/C z-axis and the zenith for each FT2 interval.  Exact."""
        ra,dec    = N.radians([skydir.ra(),skydir.dec()])
        ra_s,ra_z = self.RA_SCZ,self.RA_ZENITH
        cdec,sdec = cos(dec),sin(dec)
        scosines  = self.COS_DEC_SCZ*cdec*N.cos(ra-ra_s) + self.SIN_DEC_SCZ*sdec
        if self.zenithcut > -1:
            zcosines = self.COS_DEC_ZENITH*cdec*N.cos(ra-ra_z) + self.SIN_DEC_ZENITH*sdec
            mask = (scosines>=self.fovcut) & (zcosines>=self.zenithcut)
        else:
            mask = (scosines>=self.fovcut)
        return scosines,mask

    def __call__(self,skydir,intervals = None):
        """Return the exposure at location indicated by skydir.  Intervals is an optional list of time intervals
           that can be used to obtain a time-binned exposure for, e.g., light curves."""
      
        ra,dec = N.radians([skydir.ra(),skydir.dec()])
        if (ra == self.prev_ra) and (dec == self.prev_dec) and (intervals is None):
            return self.prev_val

        # Calculate incidence and zenith angles wrt skydir, and cuts.
        scosines,mask = self.get_cosines(skydir)

        if intervals is None:
            # Bin in cosine using livetimes as the weights
            self.prev_val = N.histogram(scosines[mask],bins=self.cosbins,weights=self.LIVETIME[mask],new=True)
            self.prev_ra,self.prev_dec = ra,dec
            return self.prev_val

        # return a time series for the livetime at the given position
        overlaps = [self.__process_ft2_fast([intervals[i][0]],[intervals[i][1]])[mask] for i in range(len(starts))]
        livetimes= [N.histogram(scosines[mask],bins=self.cosbins,weights=self.LIVETIME[mask]*ol,new=True) for ol in overlaps]
        return livetimes

#===============================================================================================#
class BinnedLivetime(Livetime):
    """See remarks for Livetime class for general information.
       
       This class provides an implementation of the livetime calculation
       in which the FT2 entries for the S/Z z-axis and zenith positions
       are binned onto a Healpix grid, allowing for a faster calculation
       with long FT2 files.
    """

    def finish(self):
        hp = Healpix(self.nside,Healpix.RING,SkyDir.EQUATORIAL)
        ras,decs = N.asarray( [hp.py_pix2ang(i) for i in range(12*self.nside**2)]).transpose()
        self.COS_HP_DEC = N.cos(decs)
        self.SIN_HP_DEC = N.sin(decs)
        self.HP_RA = ras
        ra_s,dec_s = self.RA_SCZ,self.DEC_SCZ
        ra_z,dec_z = self.RA_ZENITH,self.DEC_ZENITH
        self.S_PIX = N.fromiter((hp.py_ang2pix(ra_s[i],dec_s[i]) for i in range(len(ra_s))),dtype=int)
        self.Z_PIX = N.fromiter((hp.py_ang2pix(ra_z[i],dec_z[i]) for i in range(len(ra_z))),dtype=int)

    def __init__(self,nside=59,*args,**kwargs):
        self.nside = nside
        super(BinnedLivetime,self).__init__(*args,**kwargs)

    def get_cosines(self,skydir):
        ra,dec    = N.radians([skydir.ra(),skydir.dec()])

        # calculate the arclengths to the various Healpix
        cosines  = self.COS_HP_DEC*cos(dec)*N.cos(ra-self.HP_RA) + self.SIN_HP_DEC*sin(dec)
        scosines = cosines[self.S_PIX]
        if self.zenithcut > -1:
            zcosines = cosines[self.Z_PIX]
            mask = (scosines>=self.fovcut) & (zcosines>=self.zenithcut)
        else:
            mask = (scosines>=self.fovcut)
        return scosines,mask

#===============================================================================================#
class EfficiencyCorrection(object):
    v1 = [-1.381,  5.632, -0.830, 2.737, -0.127, 4.640]  # p0, front
    v2 = [ 1.268, -4.141,  0.752, 2.740,  0.124, 4.625]  # p1, front
    v3 = [-1.527,  6.112, -0.844, 2.877, -0.133, 4.593]  # p0, back
    v4 = [ 1.413, -4.628,  0.773, 2.864,  0.126, 4.592]  # p1, back
    
    def p(self,logE,v):
        a0,b0,a1,logEb1,a2,logEb2 = v
        b1 = (a0 - a1)*logEb1 + b0
        b2 = (a1 - a2)*logEb2 + b1
        if logE < logEb1:
            return a0*logE + b0
        if logE < logEb2:
            return a1*logE + b1
        return a2*logE + b2

    def __init__(self,e=1000):
        self.set_p(e)

    def set_p(self,e):
        loge = N.log10(e)
        for key,vec in zip(['p0f','p1f','p0b','p1b'],[self.v1,self.v2,self.v3,self.v4]):
            self.__dict__[key] = self.p(loge,vec)

    def get_efficiency(self,livetime_fraction,conversion_type=0):
        p0,p1 = (self.p0f,self.p1f) if conversion_type==0 else (self.p0b,self.p1b)
        return p0*livetime_fraction + p1

#===============================================================================================#
class EffectiveArea(object):

    def init(self):
        self.irf = 'P6_v3_diff'

    def __init__(self,CALDB,**kwargs):
        """CALDB -- path to CALDB directory"""
      
        self.init()
        self.__dict__.update(kwargs)
        self.CALDB = join(CALDB,'bcf')
        self.__read_data()

    def __read_data(self):
      
        ct0_file = join(self.CALDB,'ea','aeff_%s_front.fits'%(self.irf))
        ct1_file = join(self.CALDB,'ea','aeff_%s_back.fits'%(self.irf))
        ea = pf.open(ct0_file)
        cbins = N.append(ea['EFFECTIVE AREA'].data.field('CTHETA_LO')[0],ea['EFFECTIVE AREA'].data.field('CTHETA_HI')[0][-1])
        ebins = N.append(ea['EFFECTIVE AREA'].data.field('ENERG_LO')[0],ea['EFFECTIVE AREA'].data.field('ENERG_HI')[0][-1])
        feffarea = N.array(ea['EFFECTIVE AREA'].data.field('EFFAREA')[0])
        ea.close()
        ea = pf.open(ct1_file)
        beffarea = N.array(ea['EFFECTIVE AREA'].data.field('EFFAREA')[0])
        ea.close()
        self.cbins,self.ebins = cbins,ebins
        nc,ne = len(self.cbins),len(self.ebins)
        self.feffarea,self.beffarea = feffarea.reshape(nc-1,ne-1),beffarea.reshape(nc-1,ne-1)
        self.i_ebins,self.i_cbins = N.log((ebins[:-1]*ebins[1:])**0.5),(cbins[1:]+cbins[:-1])/2.

    def image(self,event_class=-1,logea = False):

        if event_class < 0: effarea = self.feffarea + self.beffarea
        elif event_class == 0: effarea = self.feffarea
        else: effarea = self.beffarea
        ebins,cbins = self.ebins,self.cbins

        import pylab as P

        #Generate a pseudo-color plot of the full effective area
        P.figure(2)
        P.gca().set_xscale('log')
        if logea: P.gca().set_yscale('log')
        P.pcolor((ebins[:-1]*ebins[1:])**0.5,(cbins[:-1]+cbins[1:])/2.,effarea.reshape(len(cbins)-1,len(ebins)-1))
        P.title('Effective Area')
        P.xlabel('$\mathrm{Energy\ (MeV)}$')
        P.ylabel('$\mathrm{cos( \theta)}$')
        cb = P.colorbar()
        cb.set_label('$\mathrm{Effective\ Area\ (m^2)}$')

        #Generate a plot of the on-axis effective area with and without interpolation
        energies = N.logspace(N.log10(ebins[0]),N.log10(ebins[-1]),240)
        f_vals,b_vals = N.array([self(e,.99,interpolate=True) for e in energies]).transpose()
        P.figure(4)
        P.gca().set_xscale('log')
        if logea: P.gca().set_yscale('log')
        P.plot(energies,f_vals,label='front bilinear interp.')
        P.plot(energies,b_vals,label='back bilinear interp.')
        f_vals,b_vals = N.array([self(e,.99,interpolate=False) for e in energies]).transpose()
        P.plot(energies,f_vals,label='front nearest-neighbour interp.')
        P.plot(energies,b_vals,label='back nearest-neighbour interp.')
        P.title('On-axis Effective Area')
        P.xlabel('$\mathrm{Energy\ (MeV)}$')
        P.ylabel('$\mathrm{Effective\ Area\ (cm^2)}$')
        P.legend(loc = 'lower right')
        P.grid()


    def __call__(self,e,c,event_class=-1,interpolate = True):
        """Return bilinear (or nearest-neighbour) interpolation."""
        eb,cb = self.i_ebins,self.i_cbins
        e = N.log(e)
        ne,nc = len(eb),len(cb)

        if e < eb[0]: e = eb[0]
        if e > eb[-1]: e = eb[-1]
        if c < cb[0]: c = cb[0]
        if c > cb[-1]: c = cb[-1]

        #Cute way to find nearest neighbour
        i = N.argmin(N.abs(eb-e))
        j = N.argmin(N.abs(cb-c))

        # effarea[:,-1] increasing effective area
        # effarea[-1,:] increasing effective area
        # effarea[:,-1].shape 32
        # effarea[-1,:].shape 64

        if not interpolate:
            i,j = j,i      
            if event_class < 0: return (1e4*self.feffarea[i,j],1e4*self.beffarea[i,j])
            elif event_class == 0: return 1e4*self.feffarea[i,j]
            else: return 1e4*self.beffarea[i,j]

        i = i if eb[i]<=e and i < ne-1 else i-1 #adjust nearest neighbor to lower bin
        j = j if cb[j]<=c and j < nc-1 else j-1

        def bilinear(effarea):
      
            c2,c1 = cb[j+1],cb[j]
            e2,e1 = eb[i+1],eb[i]
            f00 = effarea[j,i]
            f11 = effarea[j+1,i+1]
            f01 = effarea[j+1,i]
            f10 = effarea[j,i+1]

            return 1e4/(e2-e1)/(c2-c1)*( (e2-e)*(f00*(c2-c) + f01*(c-c1)) + (e-e1)*(f10*(c2-c) + f11*(c-c1)) )

        if event_class < 0: return (bilinear(self.feffarea),bilinear(self.beffarea))
        elif event_class == 0: return bilinear(self.feffarea)
        else: return bilinear(self.beffarea)

    def dual_avg(self,e_range,c_range,event_class=-1,steps=10):
        return N.array([self(x,y,event_class=event_class) for x in N.linspace(e_range[0],e_range[1],steps) for y in N.linspace(c_range[0],c_range[1],steps)]).sum()/steps**2

    def theta_avg(self,e,c_range,event_class=-1,steps=10):
        return N.array([self(e,ctheta,event_class=event_class) for ctheta in N.linspace(c_range[0],c_range[1],steps)]).sum()/steps

#===============================================================================================#
class Exposure(object):

   def __init__(self,ft2file,ft1files=None,fovcut=66.,zenithcut=105.,frontfile=None):
      self.ft2file = ft2file
      self.ea = EffectiveArea(frontfile = frontfile)
      self.lt = Livetime(ft2file,ft1files=ft1files,bins=self.ea.cbins,fovcut=fovcut,zenithcut=zenithcut)
      self.energies = self.event_class = None

   def __call__(self,skydir,energies,event_class=-1):

      lt = self.lt(skydir)
      if lt is None: return N.zeros(len(energies))

      #Do some caching -- effective area is the same for same energy bins!
      if N.all(energies == self.energies) and event_class == self.event_class:
         #print ('using cached values')
         vals = self.vals
      
      else:

         e_centers = N.asarray(energies)
         c_centers = (lt[1][1:]+lt[1][:-1])/2.
         
         vals = N.array([[self.ea(e_center,c_center,event_class=event_class) for c_center in c_centers] for e_center in e_centers])
         self.vals = vals    
         self.energies = energies
         self.event_class = event_class

      
      if event_class == -1:
         return N.append(N.sum(vals[:,:,0]*lt[0],axis=1),N.sum(vals[:,:,1]*lt[0],axis=1))
      else: return N.sum(vals*lt[0],axis=1)

   def change_IRF(self,frontfile = None):
      self.ea = EffectiveArea(frontfile = frontfile)
      self.energies = self.event_class = None

#===============================================================================================#
class SimpleExposureSeries(object):
    """A helper class to generate a simple exposure (evaluated at a
       single energy) time series for use in looking for periodic signals
       from sources with periods long enough to warrant exposure correction.

       The binned method is *much* faster and can be executed quickly using
       extremely fine binning, comparable to the unbinned method to much
       better than one percent, which is almost certainly smaller than the
       accuracy hit we get by ignoring interpolation of the S/C pointing.
    """

    def __init__(self,CALDB,ft1files,ft2files,
                 energy=1000,cosbins=N.linspace(0.2,1,1001),apply_correction=True,
                 **kwargs):
        self.lt = Livetime(ft2files,ft1files,**kwargs)
        self.ea = EffectiveArea(CALDB,**kwargs)
        self.en = energy
        self.ac = apply_correction
        self.apply_correction()
        self.rebin(cosbins)

    def apply_correction(self):
        if self.ac:
            # do efficiency correction
            ec = EfficiencyCorrection(self.en)
            LT0 = self.lt.LIVETIME*ec.get_efficiency(self.lt.LTFRAC,0)
            LT1 = self.lt.LIVETIME*ec.get_efficiency(self.lt.LTFRAC,1)
            self.LIVETIME = N.asarray([LT0,LT1]).transpose()
        else:
            self.LIVETIME = N.asarray([self.lt.LIVETIME,self.lt.LIVETIME]).transpose()

    def rebin(self,cosbins):
        self.cosbins = cosbins; ea = self.ea; en = self.en
        self.aeff_vals = N.asarray([ea(en,c) for c in (cosbins[1:] + cosbins[:-1])/2])

    def set_energy(self,energy):
        self.en = energy
        self.apply_correction()
        self.rebin(self.cosbins)

    def get_series(self,skydir,binned=True,energy=None):

        if energy is not None: self.set_energy(energy)

        cosines,mask = self.lt.get_cosines(skydir)
        cosines = cosines[mask]

        if binned: # use binned effective area -- faster
            indices = N.searchsorted(self.cosbins,cosines) - 1
            aeff = self.aeff_vals[indices]
        else:      # use unbinned effective area -- a little slow
            en = self.en; ea = self.ea
            aeff = N.asarray([ea(en,c) for c in cosines])
        #print (aeff.shape)
        exposure = (aeff*self.LIVETIME[mask]).sum(axis=1)
        return self.lt.START[mask],self.lt.STOP[mask],exposure

#===============================================================================================#
class SpectralExposureSeries(object):
    """Average SimpleExposureSeries over energy using the given spectral
       index and the given energy points."""

    def __init__(self,ses,index=2,emin=100,emax=1e5,nsimps=12):
        self.ses = ses
        self.e_points = sp = N.logspace(N.log10(emin),N.log10(emax),nsimps+1)
        sweights = (N.log(sp[-1]/sp[0])/(3.*nsimps)) * \
                    N.asarray([1.] + ([4.,2.]*(nsimps/2))[:-1] + [1.])
        sweights *= (self.e_points/emin)**(-index) # energy weighting
        self.s_weights = sweights / (sweights.sum())# normalization
    
    def get_series(self,skydir):
        ep,ew,ses = self.e_points,self.s_weights,self.ses
        t1,t2,results = ses.get_series(skydir,energy=ep[0])
        results *= ew[0]
        for i in range(1,len(ep)):
            results += ses.get_series(skydir,energy=ep[i])[-1]*ew[i]
        return t1,t2,results
        
        