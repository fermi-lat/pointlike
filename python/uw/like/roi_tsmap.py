"""
Module implements a TS calculation, primarily for source finding / fit verification.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_tsmap.py,v 1.3 2010/07/31 00:35:24 kerrm Exp $

author: Matthew Kerr
"""

import numpy as N
from uw.like.pypsf import PsfOverlap
from uw.like.Models import PowerLaw
from skymaps import Band,WeightedSkyDirList,PySkyFunction,Hep3Vector,SkyDir,BinnedPhotonData,SkyImage
from pointlike import IntVector,DoubleVector
from scipy.optimize import fmin
from cPickle import dump,load
from glob import glob
from collections import deque

def get_latalog(latalog_file='f:/glast/data/kerr/gll_psc_v02.fit'):
    from latalog import Latalog
    return Latalog(latalog_file)

###====================================================================================================###
class TSCalc(object):
    """Extract a TS as a function of sky position.

       The implementation is different from that in ROILocalizer in that, for each position on the sky,
       a model is fit.  (Whereas, for ROILocalizer, the same broadband model is moved about the sky
       without a refit.)  Thus, this version is appropriate for source finding.
    """
    
    def __init__(self,roi):
        self.roi = roi
        self.ro  = PsfOverlap()
        self.rd  = roi.roi_dir
        self.phase_factor = roi.phase_factor
        
        mo  = PowerLaw(p=[1e-12,2])
        self.n0 = 10**mo.p[0]
        for band in roi.bands:
            band.ts_exp = band.expected(mo) / self.n0
        
        self.seeds = N.concatenate(([-20.],N.arange(-14,-8.9,0.5)))

    def __call__(self,skydir,repeat_diffuse=False,bright_source_mask=None,no_cache=False):
        """Return the TS for the position on the sky given by the argument.
        
           bright_sources = a mask to select sources to include with the
           diffuse when generating the TS map."""

        bands = self.roi.bands
        bsm   = bright_source_mask
        dv    = DoubleVector()
        pf    = self.phase_factor

        if not repeat_diffuse or no_cache:

            for i,band in enumerate(bands):

                en,exp,pa = band.e,band.exp.value,band.b.pixelArea()

                # make a first-order correction for exposure variation
                band.ts_er = exp(skydir,en)/exp(self.rd,en)
          
                # unnormalized PSF evaluated at each data pixel
                #diffs = N.asarray([skydir.difference(w) for w in band.wsdl])
                #band.ts_pix_counts = pa * band.psf(diffs,density=True)
                band.wsdl.arclength(skydir,dv)
                band.ts_pix_counts = pa * band.psf(N.fromiter(dv,dtype=float),density=True)

                # calculate overlap
                band.ts_overlap = self.ro(band,self.rd,skydir)

        if not repeat_diffuse:

            for i,band in enumerate(bands):

                # pre-calculate the "pixel" part
                if band.has_pixels:
                    band.ts_pix_term = (band.ps_all_pix_counts + band.bg_all_pix_counts)/ \
                                       (band.ts_exp*band.ts_pix_counts)

        else:
        
            # include bright point sources and diffuse in the background model
            if bsm is not None:
                for i,band in enumerate(bands):
                    if band.has_pixels:
                        bps_term = (band.ps_counts[bsm] * band.overlaps[bsm] * band.ps_pix_counts[:,bsm]).sum(axis=1)
                        band.ts_pix_term = (bps_term + band.bg_all_pix_counts) / (band.ts_exp*band.ts_pix_counts)
            
            # include only the diffuse in the background model
            else:
                for i,band in enumerate(bands):
                    if band.has_pixels:
                        band.ts_pix_term = band.bg_all_pix_counts / (band.ts_exp*band.ts_pix_counts)


        # NB -- can save some computation by calculating f0/f1/f2 simultaneously, but it is
        # currently a minute fraction of the total time (above code dominates)
        J = N.log(10)

        def f0(n0,*args):
            n0 = 10**n0
            accum = 0
            for band in bands:
                pix_term = (band.pix_counts*N.log(1 + n0/band.ts_pix_term)).sum() if band.has_pixels else 0
                ap_term  = - n0*band.ts_overlap*band.ts_exp*pf
                accum   += pix_term + ap_term
            return accum

        def f1(n0,*args):
            n0 = 10**n0
            accum = 0
            for band in bands:
                pix_term = - (band.pix_counts*(1 + band.ts_pix_term/n0)**-1).sum() if band.has_pixels else 0
                ap_term  = n0*band.ts_exp*band.ts_overlap*pf
                accum   += pix_term + ap_term
            return J*accum

        def f2(n0,*args):
            n0 = 10**n0
            accum = 0
            for band in bands:
                if band.has_pixels:
                    quot     = band.ts_pix_term/n0
                    pix_term = - (band.pix_counts*quot/(1+quot)**2).sum()
                else:
                    pix_term = 0
                ap_term  = n0*band.ts_exp*band.ts_overlap*pf
                accum   += pix_term + ap_term
            return J*J*accum

        def TS(n0,*args): return 2*f0(n0,*args)

        # assess along a grid of seed values to make sure we have a good starting position
        vals  = [f0(x) for x in self.seeds]
        amax  = N.argmax(vals)
        if amax == 0: return 0
        seed  = self.seeds[amax] + 0.5 # for some reason, helps to start *above* the critical point

        # re-implementation of scipy version that uses half the calls!               
        def my_newton(func,x0,fprime,tol=1e-2):
            p0 = x0
            for i in xrange(30):
                fval = func(x0)
                if fval == 0: return x0,True
                gval = fprime(x0)
                delt = fval/gval
                x0  -= delt
                if (abs(delt) < tol):
                    return x0,True
            return x0,False

        n0,conv = my_newton(f1,seed,fprime=f2)
        if conv: return TS(n0)
        else:
            print 'Warning! did not converge to a value or a value consistent with 0 flux.'
            print 'Trying again...'
            n0,conv = my_newton(f1,n0,fprime=f2)
            if conv:
                print 'Converged on 2nd Try'
                return TS(n0)
            print 'DID NOT CONVERGE AFTER TWO ATTEMPTS'
            return -1

###====================================================================================================###
class TSCalcPySkyFunction(object):

    def __init__(self,tscalc,repeat_diffuse=False,bright_source_mask=None):
        self.tscalc = tscalc
        self.repeat_diffuse=repeat_diffuse
        self.bright_source_mask=bright_source_mask

    def __call__(self,v):
        sd = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        return self.tscalc(sd,repeat_diffuse=self.repeat_diffuse,bright_source_mask=self.bright_source_mask,no_cache=True)

    def get_pyskyfun(self):
        return PySkyFunction(self)

###====================================================================================================###
class HealpixTSMap(object):

    def __init__(self,hr=None,pickle=None,factor = 100,bright_sources=None):
        """hr     :   an instance of HealpixROI
           pickle :   a pickle saved from a previous HealpixTSMap
           factor :   number of intervals to divide base Healpixel side into"""

        self.bright_sources = bright_sources
        if   hr     is not None: self.setup_from_roi(hr,factor)
        elif pickle is not None: self.load(pickle)
        else:
            print 'Must provide either a HealpixROI instance or a pickle file!'
            raise Exception

    def setup_from_roi(self,hr,factor):

        band1  = hr.band
        band2  = Band(int(round(band1.nside()*factor)))
        rd     = band1.dir(hr.index) 

        # get pixels within a radius 10% greater than the base Healpix diagonal dimension
        radius = (N.pi/(3*band1.nside()**2))**0.5 * 2**0.5 * 1.1
        wsdl   = WeightedSkyDirList(band2,rd,radius,True)
        
        # then cut down to just the pixels inside the base Healpixel
        inds   = N.asarray([band1.index(x) for x in wsdl])
        mask   = inds == hr.index
        dirs   = [wsdl[i] for i in xrange(len(mask)) if mask[i]]
        inds   = N.asarray([band2.index(x) for x in dirs]).astype(int)

        # sanity check
        if abs(float(mask.sum())/factor**2 - 1) > 0.01:
            print 'Warning: number of pixels found does not agree with expectations!'

        # calculate TS for each sub-pixel
        tsc      = TSCalc(hr.roi)
        if self.bright_sources is not None:
            ts_vals = self.ts_vals = [deque(),deque(),deque()]
            bsm = N.asarray([x.source_name in self.bright_sources for x in hr.roi.psm.models])
            for d in dirs:
                ts_vals[0].append(tsc(d))
                ts_vals[1].append(tsc(d,repeat_diffuse=True))
                ts_vals[2].append(tsc(d,repeat_diffuse=True,bright_source_mask=bsm))
        else:
            self.ts_vals = [deque(),deque()]
            for d in dirs:
                ts_vals[0].append(tsc(d))
                ts_vals[1].append(tsc(d,repeat_diffuse=True))

        self.band1 = band1
        self.band2 = band2

        sorting    = N.argsort(inds)
        self.inds  = inds[sorting]
        for i in range(len(ts_vals)):
            ts_vals[i] = N.asarray(ts_vals[i])[sorting]

        self.previous_index = -1
        self.index          = hr.index
        self.mode           = 0
        self.ts             = self.ts_vals[0]

    def dump(self,filename):
        d = dict()
        d['inds'] = self.inds
        d['ts_vals']= self.ts_vals
        d['nside1'] = self.band1.nside()
        d['nside2'] = self.band2.nside()
        d['index']  = self.index
        dump(d,file(filename,'w'))
        
    def load(self,filename):
        d = load(file(filename))
        self.inds   = d['inds']
        self.band1  = Band(d['nside1'])
        self.band2  = Band(d['nside2'])
        if 'ts_vals' in d.keys():
            self.ts_vals= d['ts_vals']
        else:
            dim = (self.band2.nside()/self.band1.nside())**2
            self.ts_vals = N.empty([3,dim])
            self.ts_vals[:] = N.nan
        self.index  = d['index']
        self.ts     = self.ts_vals[0]
        self.previous_index = -1
        self.mode           = 0

    def pyskyfun(self):
        return PySkyFunction(self)

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        if self.band1.index(sd) != self.index: return N.nan
        id = self.band2.index(sd)
        if id != self.previous_index:
            self.previous_value = self.ts[N.searchsorted(self.inds,id)]
            self.previous_index = id
        return self.previous_value

    def multipix_call(self,sd):
        """A faster version of the skydir lookup where it assumed that the argument
           lies within the boundaries of this Healpixel."""
        id = self.band2.index(sd)
        if id != self.previous_index:
            self.previous_value = self.ts[N.searchsorted(self.inds,id)]
            self.previous_index = id
        return self.previous_value

    def set_mode(self,mode=1):
        if mode not in [0,1,2]:
            print 'Not a valid mode... taking no action.'
            return
        self.ts   = self.ts_vals[mode]
        self.mode = mode

class HealpixMap(object):
    """map from a set of Healpix to values, supporting a SkyFunction interface
       and interpolation.
       
       The implentation here is not really appropriate for a sparse map (much less
       than all-sky if the number of pixels are large, since the storage system is
       a vector with indices==healpix_index.
       
       An alternative (perhaps realizable with a child class) can use a sorted
       lookup and store only populated entries."""

    def __init__(self,nside,inds,vals):
        self.band = Band(nside)
        #s = N.argsort(inds)
        #self.inds = inds[s]
        #self.vals = vals[s]
        self.vals = N.empty(12*nside**2)
        self.vals[:] = N.nan
        self.vals[inds] = vals
        
        from pointlike import IntVector
        self.iv = IntVector()

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        healpix_index = self.band.index(sd)
        return self.vals[healpix_index]

    def bilinear_interp(self,v,skydir=None):
        b = self.band
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        healpix_index = b.index(sd)
        d = b.dir(healpix_index)
        b.findNeighbors(healpix_index,self.iv)
        idx   = N.asarray([x for x in iv] + [healpix_index])
        dirs  = [b.dir(idx) for idx in self.iv] + [d]
        diffs = N.asarray([sd.difference(di) for di in dirs])

        # use co-ordinate system that puts pixels closest to equator
        use_gal = abs(d.b()) < abs(d.dec())
        if use_gal:
            lons = N.asarray([d.b() for d in dirs])
            lats = N.asarray([d.l() for d in dirs])
            sdlon = sd.b();  sdlat = sd.l()
        else:
            lons = N.asarray([d.ra() for d in dirs])
            lats = N.asarray([d.dec() for d in dirs])
            sdlon = sd.ra(); sdlat = sd.dec()
        
        s = N.argsort(diffs)
        lons = lons[s][:4]
        lats = lats[s][:4]
        vals = N.asarray([self(0,skydir=dirs[x]) for x in s[:4]])
        return (N.abs(lons-sdlon)*N.abs(lats-sdlat)*vals).sum()        
    
        

class MultiHealpixTSMap(object):

    def __init__(self,glob_expansion):
        """glob_expansion: an expression that can be expanded by glob into a list of the
                           pickles of HealpixTSMaps that this class will access
                           
                           N.B. it is necessary that these objects all have been constructed
                           with the same settings of base nside and sub-grid factor!
        """

        tsmaps = [HealpixTSMap(pickle=x)for x in glob(glob_expansion)]
        nside1 = tsmaps[0].band1.nside()
        
        self.band1 = Band(nside1)
        self.tsmaps = [False] * 12 * nside1**2
        for tsmap in tsmaps:
            self.tsmaps[tsmap.index] = tsmap
        self.scale = 0

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        id = self.band1.index(sd)
        tsmap = self.tsmaps[id]
        if not tsmap: return N.nan
        return tsmap.multipix_call(sd)

    def set_mode(self,mode=1):
        """Change display mode.
           1 == TS map with all known point sources included in models
           2 == TS map using only diffuse models for background
        """
        if mode not in [0,1,2]:
            print 'Error!  Not a valid mode.  Please use 1 or 2.'
            return
        for tsmap in self.tsmaps:
            if tsmap: tsmap.set_mode(mode)

    def make_zea(center,size=10,pixelsize=0.02,galactic=False):
        from uw.utilities.image import ZEA
        z = ZEA(center,size=size,pixelsize=pixelsize,galactic=galactic)
        z.fill(PySkyFunction(self))
        self.z = z

    def get_1FGL_sources(self,z):
        """z == ZEA object"""
        fgl = get_latalog()
        sd_ul = z.skydir(0,0)        # upper left corner
        sd_lr = z.skydir(z.nx,z.ny)  # lower right corner
        if z.galactic:
            az_max,po_min = sd_ul.l(),sd_ul.b()
            az_min,po_max = sd_lr.l(),sd_lr.b()       
        else:
            az_max,po_min = sd_ul.ra(),sd_ul.dec()
            az_min,po_max = sd_lr.ra(),sd_lr.dec()
        
        # since longitude runs backward, this condition means we looped over 0
        az_zero_cross = az_min > az_max
        
        az_key,po_key = ['GLON','GLAT'] if z.galactic else ['RA','DEC']
        az = fgl.get(az_key).astype(float)
        po = fgl.get(po_key).astype(float)
        na = fgl.get('Source_Name')

        mask = (po > po_min) & (po < po_max)
        if az_zero_cross:
            mask = mask & ((az < az_max) | (az > az_min))
        else:
            mask = mask & (az > az_min) & (az < az_max)
        sds = [SkyDir(a,p,SkyDir.GALACTIC if z.galactic else SkyDir.EQUATORIAL) for a,p in zip(az[mask],po[mask])]
        return na[mask],sds

    def make_map(self,center,size=10,pixelsize=0.02,galactic=False,axes=None,
                 scale=1,thresh_low=0,thresh_high=N.inf,label_1FGL=True,mode=1,cmap=None,
                 log_transform=False,labels=None):
        """
            scale : 0 == linear, 1 == sqrt
        """
        import pylab as P
        from uw.utilities.image import ZEA
        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize

        self.set_mode(mode)
        
        if axes is None: axes = P.gca()
        cmap = cmap or P.cm.jet
        cmap.set_bad('white')
        z = ZEA(center,size=size,pixelsize=pixelsize,galactic=galactic,axes=axes)
        z.fill(PySkyFunction(self))
        im = z.image.copy()
        print 'Max Value: %.2f'%(im.max())
        im[im < thresh_low]  = N.nan
        im[im > thresh_high] = thresh_high
        z.grid()
        im = im**(0.5 if scale else 1)
        if log_transform:
            im[im < 1] = 1
            im = N.log10(im)
        P.imshow(im,cmap=cmap)
        #P.colorbar(ax=axes,cmap=cmap,orientation='horizontal')
        #P.contour(im,N.asarray([9,16,25,50,100,1000])**(0.5 if scale else 1))        
        #P.contour(im,10)
        #contours = N.asarray([9,16,25,36,49,64,81,100,225,400,625,900,1225,1600])**(0.5 if scale else 1)
        #if log_transform: contours = N.log10(contours)
        #P.contour(im,contours)
        if label_1FGL:
            names,skydirs = self.get_1FGL_sources(z)
            for na,sd in zip(names,skydirs):
                z.plot_source(na,sd,color='gray' if thresh_low > 0 else 'white')
        if labels is not None:
            try:
                names,skydirs = labels
            except:
                skydirs = labels
                names   = ['S%d'%(i) for i in xrange(len(labels))]
            for na,sd in zip(names,skydirs):
                z.plot_source(na,sd,color='gray' if thresh_low > 0 else 'white')
        return z
        #colorbar = ColorbarBase(axes,orientation='vertical',cmap=cmap)

    def make_3panel(self,center,size=10,pixelsize=0.02,galactic=False,scale=1,
                    label_1FGL=True,labels=None,separate_figures=False,fignum_base=10):
        
        import pylab as P
        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize
        from uw.utilities.image import ZEA

        cmap = P.cm.jet
        cmap.set_bad('white')
        
        if separate_figures:
            axes = []
            for i in xrange(1,4):
                P.figure(i+fignum_base)
                axes += [P.gca()]
        else:
            axes = [P.subplot(1,3,i) for i in xrange(1,4)]
        zeas = [ZEA(center,size=size,pixelsize=pixelsize,galactic=galactic,axes=ax) for ax in axes]
        mods = [1,2,0]
        for i in xrange(0,3):
            self.set_mode(mods[i])
            zeas[i].fill(PySkyFunction(self))
            zeas[i].grid()
            axes[i].imshow(zeas[i].image**(0.5 if scale else 1),cmap=cmap)
        if label_1FGL:
            names,skydirs = self.get_1FGL_sources(zeas[0])
            for na,sd in zip(names,skydirs):
                for z in zeas:
                    z.plot_source(na,sd,color='white')
        if labels is not None:
            names,skydirs = labels
            for na,sd in zip(names,skydirs):
                for z in zeas:
                    z.plot_source(na,sd,color='white')
        if separate_figures:
            # basically a "show"
            for i in xrange(1,4):
                P.figure(i+fignum_base)
                axes[i-1].set_autoscale_on(True)
                cb = ColorbarBase(axes[i-1],orientation='vertical',cmap=cmap)
                #P.colorbar()
        return axes,zeas
                                  
class DataMap(object):
    """Read from a BinnedPhotonData with a single pixel size.  Generate a counts map."""

    def __init__(self,bpdfile,**kwargs):
        self.bpd    = BinnedPhotonData(bpdfile)
        self.band   = self.bpd[0]
        self.smooth = False
        self.kernel = (4*N.pi/(12*self.band.nside()**2))**0.5

        self.iv = IntVector()
        self.previous_index = -1
        self.previous_vals  = None
        self.previous_dirs  = None

        self.__dict__.update(kwargs)

    def set_band(self,index):
        if index >= len(self.bpd):
            raise Exception
        self.band = self.bpd[index]

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        b   = self.band
        if not self.smooth: return b(sd)
        
        idx = b.index(sd)
  
        if idx == self.previous_index:
            neighbor_dirs = self.previous_dirs
            neighbor_vals = self.previous_vals
        else:
            b.findNeighbors(idx,self.iv)
            neighbor_dirs = [b.dir(idx)] + [b.dir(i) for i in self.iv]
            neighbor_vals = [b(d) for d in neighbor_dirs]
            self.previous_index = idx
            self.previous_dirs = neighbor_dirs
            self.previous_vals = neighbor_vals
        diffs   = N.asarray([sd.difference(d) for d in neighbor_dirs])
        weights = N.exp(-0.5*(diffs/self.kernel)**2)
        return (neighbor_vals*weights).sum()/weights.sum()

    def pyskyfun(self):
        return PySkyFunction(self)

    def get_AIT(self,**kwargs):
        from uw.utilities.image import AIT
        a = AIT(self.pyskyfun(),**kwargs)
        a.fill(a.skyfun)
        a.grid()
        a.axes.imshow(a.image)
        return a

    def get_ZEA(self,center,**kwargs):
        from uw.utilities.image import ZEA
        z = ZEA(center,**kwargs)
        z.fill(self.pyskyfun())
        z.grid()
        z.axes.imshow(z.image)
        return z

    def get_vector(self):
        b = self.band
        return N.asarray([self(b.dir(i)) for i in xrange(12*b.nside()**2)])

def get_seeds(self,mode=0,ts_thresh=9):
    """Do a crude test to try to find good seeds."""

    ts   = self.ts_vals[mode]
    mask = N.asarray([True]*len(ts))
    b2   = self.band2
    inds = self.inds
    iv   = IntVector()
    seeds = deque()
    seeds_ts = deque()

    n = len(inds)

    for i in xrange(n):
        my_ts  = ts[mask]
        amax   = N.argmax(my_ts)
        mts    = my_ts[amax]
        if mts < ts_thresh:
            break
        idx    = int(inds[mask][amax])
        b2.findNeighbors(idx,iv)
        neigh_ts = deque()
        for neigh in [x for x in iv]:
            arg = N.searchsorted(inds,neigh)
            if (arg < n) and (inds[arg] == neigh):
                mask[arg] = False
                neigh_ts.append(ts[arg])
        arg = N.searchsorted(inds,idx)
        mask[arg] = False
        if mts > max(neigh_ts):
            seeds.append(idx)
            seeds_ts.append(mts)

    return [b2.dir(x) for x in seeds],N.asarray(seeds_ts)

###====================================================================================================###
class HealpixKDEMap(object):

    def __init__(self,hr=None,pickle=None,factor = 100):
        """hr     :   an instance of HealpixROI
           pickle :   a pickle saved from a previous HealpixTSMap
           factor :   number of intervals to divide base Healpixel side into"""

        if   hr     is not None: self.setup_from_roi(hr,factor)
        elif pickle is not None: self.load(pickle)
        else:
            print 'Must provide either a HealpixROI instance or a pickle file!'
            raise Exception

    def setup_from_roi(self,hr,factor):

        band1  = hr.band
        band2  = Band(int(round(band1.nside()*factor)))
        rd     = band1.dir(hr.index) 

        # get pixels within a radius 10% greater than the base Healpix diagonal dimension
        radius = (N.pi/(3*band1.nside()**2))**0.5 * 2**0.5 * 1.1
        wsdl   = WeightedSkyDirList(band2,rd,radius,True)
        
        # then cut down to just the pixels inside the base Healpixel
        inds   = N.asarray([band1.index(x) for x in wsdl])
        mask   = inds == hr.index
        dirs   = [wsdl[i] for i in xrange(len(mask)) if mask[i]]
        inds   = N.asarray([band2.index(x) for x in dirs]).astype(int)

        # sanity check
        if abs(float(mask.sum())/factor**2 - 1) > 0.01:
            print 'Warning: number of pixels found does not agree with expectations!'

        # loop through the bands and image pixels to calculate the KDE
        from pointlike import DoubleVector
        dv = DoubleVector()
        img = N.zeros(len(inds))
        weights = [ N.asarray([x.weight() for x in b.wsdl]).astype(int) for b in hr.roi.bands ]
        for idir,mydir in enumerate(dirs):
            #print 'Processing pixel %d'%(idir)
            for iband,band in enumerate(hr.roi.bands):
                band.wsdl.arclength(mydir,dv)
                img[idir] += (band.psf(N.fromiter(dv,dtype=float),density=True)*weights[iband]).sum()

        self.band1 = band1
        self.band2 = band2
        self.previous_index = -1

        sorting    = N.argsort(inds)
        self.inds  = inds[sorting]
        self.img   = img[sorting]

        self.index = hr.index

    def dump(self,filename):
        d = dict()
        d['inds'] = self.inds
        d['img']  = self.img
        d['nside1'] = self.band1.nside()
        d['nside2'] = self.band2.nside()
        d['index']  = self.index
        dump(d,file(filename,'w'))
        
    def load(self,filename):
        d = load(file(filename))
        self.inds   = d['inds']
        self.band1  = Band(d['nside1'])
        self.band2  = Band(d['nside2'])
        self.img    = d['img']
        self.index  = d['index']
        self.previous_index = -1

    def pyskyfun(self):
        return PySkyFunction(self)

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        if self.band1.index(sd) != self.index: return N.nan
        id = self.band2.index(sd)
        if id != self.previous_index:
            self.previous_value = self.img[N.searchsorted(self.inds,id)]
            self.previous_index = id
        return self.previous_value

    def multipix_call(self,sd):
        """A faster version of the skydir lookup where it assumed that the argument
           lies within the boundaries of this Healpixel."""
        id = self.band2.index(sd)
        if id != self.previous_index:
            self.previous_value = self.img[N.searchsorted(self.inds,id)]
            self.previous_index = id
        return self.previous_value

class ROIWrapper(object):
    """Wrap up an ROI as if it were a HealpixROI, for testing."""

    def __init__(self,roi,nside=6):
        self.roi = roi
        self.band = Band(6)
        self.index = self.band.index(roi.roi_dir)

class MultiHealpixKDEMap(object):

    def __init__(self,glob_expansion):
        """glob_expansion: an expression that can be expanded by glob into a list of the
                           pickles of HealpixTSMaps that this class will access
                           
                           N.B. it is necessary that these objects all have been constructed
                           with the same settings of base nside and sub-grid factor!
        """

        tsmaps = [HealpixKDEMap(pickle=x)for x in glob(glob_expansion)]
        nside1 = tsmaps[0].band1.nside()
        
        self.band1 = Band(nside1)
        self.tsmaps = [False] * 12 * nside1**2
        for tsmap in tsmaps:
            self.tsmaps[tsmap.index] = tsmap
        self.scale = 0

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        id = self.band1.index(sd)
        tsmap = self.tsmaps[id]
        if not tsmap: return N.nan
        return tsmap.multipix_call(sd)

    def make_zea(self,center,size=10,pixelsize=0.02,galactic=False):
        z = ZEA(center,size=size,pixelsize=pixelsize,galactic=galactic)
        z.fill(PySkyFunction(self))
        self.z = z
