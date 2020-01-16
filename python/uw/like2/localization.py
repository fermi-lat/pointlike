"""
source localization support

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/localization.py,v 1.34 2018/01/27 15:37:17 burnett Exp $

"""
import os,sys
import numpy as np
from skymaps import SkyDir
from uw.like import quadform
from uw.utilities import keyword_options
from . import (sources, plotting )

def moment_analysis(tsmap, wcs, fudge=1.44):
    """ perform localization by a moment analysis of a TS map
        tsmap : array of float: TS values on a grid, must be square
        wcs : Projector object
            implements pix2sph function of two ints to return ra,dec
        fudge : float
            Additional factor to multiply the ellipse radii
            (Determined empirically)
    returns: 
        ra, dec, ax, bx, ang
    """
    vals = np.exp(-0.5* tsmap**2).flatten(); 
    peak_fraction = vals.max()/sum(vals)
    

    n = len(vals)
    nx = ny =int(np.sqrt(n))
    #centers of pixels have index +0.5
    ix = np.array([ i % nx for i in range(n)]) +0.5
    iy = np.array([ i //nx for i in range(n)]) +0.5
    norm = 1./sum(vals)
    t = [sum(u*vals)*norm for u in  ix,iy, ix**2, ix*iy, iy**2]
    center = (t[0],t[1])
    C = np.matrix(center)
    variance = (np.matrix(((t[2], t[3]),(t[3], t[4]))) - C.T * C)
    ra,dec = wcs.pix2sph(center[0]+0.5,center[1]+0.5)
    peak = SkyDir(ra,dec)
    # get coords of center, measure degrees/pixel
    nc = (nx+1)/2
    rac, decc = wcs.pix2sph(nc, nc)
    scale = wcs.pix2sph(nc, nc+1)[1] - decc
    size = nx*scale
    # adjust variance
    variance = scale**2 * variance
    offset = np.degrees(peak.difference(SkyDir(rac,decc)))

    # add effects of binsize
    var  = variance #NO+ np.matrix(np.diag([1,1]))*(scale/3)**2
    #Eigenvalue analysis to get ellipse coords
    u,v =np.linalg.eigh(var)
    ang =np.degrees(np.arctan2(v[1,1], -v[1,0]))
    if min(u)< 0.5* max(u): 
        print ('Too elliptical : %s, setting circular' % u)
        u[0]=u[1] = max(u)
    tt = np.sqrt(u) * fudge
    if u[1]>u[0]:
        ax,bx = tt[1], tt[0]
        ang = 90-ang
    else:
        ax,bx = tt
    return ra, dec, ax,bx, ang

class MomentAnalysis(object):
    """ localization using moment analysis
    """
    def __init__(self, tsplot, fudge=1.44):
        """tsplot : TSPlot object
        """
        self.tsp=tsplot
        zea = tsplot.zea
        wcs, tsmap = zea.projector, zea.image
        self.ellipse = moment_analysis(tsmap, wcs, fudge)
        
    def moments(self):
        tsmap = self.tsp.zea.image
        vals = np.exp(-0.5* tsmap**2).flatten(); 
        peak_fraction = vals.max()/sum(vals)
        n = len(vals)
        nx = ny =int(np.sqrt(n))
        ix = np.array([ i % nx for i in range(n)]) +0.5
        iy = np.array([ i //nx for i in range(n)]) +0.5
        norm = 1./sum(vals)
        t = [sum(u*vals)*norm for u in  ix,iy, ix**2, ix*iy, iy**2]
        return t

    def drawit(self):
        
        self.tsp.overplot(self.ellipse, color='w', lw=2, ls='-', contours=[2.45])
        self.tsp.plot(SkyDir(*self.ellipse[:2]), color='w', symbol='o' )
        return self.tsp.zea.axes.figure

        
def full_localization(roi, source_name=None, ignore_exception=False, 
            update=False, associator=None, tsmap_dir='tsmap_fail', tsfits=False, delta_ts_bad=10):
    import pylab as plt

    source = roi.sources.find_source(source_name)
    source.ellipsex = None # in case already had a moment analysis
    tsp=None
    with roi.tsmap_view(source.name) as tsm:

        loc = Localization(tsm)
        try:
            if not loc.localize():
                print ('Failed')
            if hasattr(loc, 'ellipse') and  (update or loc['qual']<1.0 and loc['a']<0.1):
                # Automatically update position if good fit.
                t = loc.ellipse
                prev = tsm.saved_skydir
                tsm.saved_skydir = SkyDir(t['ra'], t['dec'])
                print ('updated position: %s --> %s' % (prev, tsm.saved_skydir))
            else:
                print ('Failed localization')
        except Exception, msg:
            print ('Localization of %s failed: %s' % (source.name, msg))
            if not ignore_exception: raise

        if not roi.quiet and hasattr(loc, 'niter') and loc.niter>0: 
            print ('Localized %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \)
                (source.name, loc.niter, loc.delt, loc.delta_ts)
            labels = 'ra dec a b ang qual'.split()
            print ((len(labels)*'%10s') % tuple(labels))
            p = loc.qform.par[0:2]+loc.qform.par[3:7]
            print (len(p)*'%10.4f' % tuple(p))
        if associator is not None:
            try:
                make_association(source, loc.TSmap, associator, quiet=roi.quiet)
            except Exception, msg:
                print ('Exception raised associating %s: %s' %(source.name, msg))
        
        if tsmap_dir is not None : 
            if  hasattr(loc,'ellipse'): 
                a, qual, delta_ts = loc.ellipse['a'], loc.ellipse['qual'], loc.delta_ts
                tsize = min(a*15., 2.0)
                bad = a>0.25 or qual>5 or abs(delta_ts)>delta_ts_bad
                if bad:
                    print ('Flagged as possibly bad: a=%.2f>0.25 or qual=%.1f>5 or abs(delta_ts=%.1f)>%f:'% (a, qual, delta_ts,delta_ts_bad))
            else: 
                print ('no localization')
                bad = True
                tsize= 2.0
            if tsmap_dir.endswith('fail') and not bad: return

            # Make tsmap and apply moment analysis if failed fit or quality cuts
            done = False
            while not done:
                try:
                    tsp=plotting.tsmap.plot(loc, source.name, center=tsm.saved_skydir,
                                        outdir=tsmap_dir, catsig=0, size=tsize, 
                                        pixelsize= tsize/15, # was 14: desire to have central pixel
                                        assoc=source.__dict__.get('adict', None), # either None or a dictionary
                                        notitle=True, #don't do title
                                        markersize=10,
                                        primary_markersize=12,
                                        tsfits=tsfits,
                                        )
                    zea = tsp.zea
                    wcs = zea.projector
                    tsmap = zea.image
                    vals = np.exp(-0.5* tsmap**2).flatten(); 
                    peak_fraction = vals.max()/sum(vals)

                except Exception, msg:
                    print ('Plot of %s failed: %s' % (source.name, msg))
                    return None
                if peak_fraction<0.8: 
                    done = True
                else:
                    #scale is too large: reduce it
                    tsize /=2.
                    print ('peak fraction= %0.2f: setting size to %.2f' % (peak_fraction, tsize))
            ellipsex = moment_analysis(zea.image, wcs)
            source.ellipsex= list(ellipsex) + [tsize, peak_fraction] # copy to the source object
            print ('moment analysis ellipse:', np.array(ellipsex))
            rax, decx, ax,bx,phi = ellipsex
            tsp.overplot([rax,decx,ax,bx, phi], color='w', lw=2, ls='-', contours=[2.45])
            tsp.plot(SkyDir(rax,decx), color='w', symbol='o' );
            filename = source.name.replace(' ','_').replace('+','p')
            fout = os.path.join(tsmap_dir, ('%s_tsmap.jpg'%filename) )
            print ('saving updated tsplot with moment analysis ellipse to %s...' % fout ; sys.stdout.flush())
            plt.savefig(fout, bbox_inches='tight', padinches=0.2) #cuts off outherwise  
            
        return tsp
        
        
class Localization(object):
    """ manage localization of a source
    Implements a minimization interface
    see also the localize function, which uses the eliptical fitter
   
    """
    defaults = (
        ('tolerance',1e-4),
        ('verbose',False),
        ('update',False,"Update the source position after localization"),
        ('max_iteration',15,"Number of iterations"),
        #('bandfits',True,"Default use bandfits"),
        ('maxdist',1,"fail if try to move further than this"),
        ('seedpos', None, 'if set, start from this position instead of the source position'),
        ('factor', 1.0,  'factor to divide the likelihood for systmatics'),
        ('quiet', False, 'set to suppress output'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, tsm, **kwargs):
        """ 
        tsm : a TSmap object, with a source selected
            It defines a function that returns the TS, or 2x the likelihood ratio of a position with respect to the 
            source position
        """
        keyword_options.process(self, kwargs)
        
        self.tsm = tsm # roistat.tsmap_view(source_name)
        self.maxlike = self.log_like()
        self.skydir  = self.tsm.skydir
        if self.seedpos is not None: 
            if not isinstance(self.seedpos, SkyDir):
                self.seedpos = SkyDir(*self.seedpos)
            self.skydir = self.seedpos
        self.name = self.tsm.source.name
        if self.factor!=1.0: 
            print ('Applying factor {:.2f}'.format(self.factor))
    
    def log_like(self, skydir=None):
        """ return log likelihood at the given position"""
        return self.tsm(skydir)/2
   
    def TSmap(self, skydir):
        """ return the TS at given position, or 
            2x the log(likelihood ratio) from the nominal position
        """
        val= 2*(self.log_like(skydir)-self.maxlike)
        return val / self.factor

    # the following 3 functions are for a minimizer
    def get_parameters(self):
        return np.array([self.tsm.skydir.ra(), self.tsm.skydir.dec()])
    
    def set_parameters(self, par):
        self.skydir = SkyDir(par[0],par[1])
        self.tsm.skydir = self.tsm.set_dir(self.skydir)
        
    def __call__(self, par):
        # for a minimizer
        return -self.TSmap(SkyDir(par[0],par[1]))
    
    def reset(self):
        """ restore modifications to the source
        """
        self.tsm.reset()
      
    def dir(self):
        return self.skydir

    def errorCircle(self):
        return 0.05 #initial guess

    def spatialLikelihood(self, sd): #negative for legacy code below
        return -self.log_like(sd)
        
    def localize(self):
        """Localize a source using an elliptic approximation to the likelihood surface.

            return fit position, number of iterations, distance moved, delta TS
        """
        #roi    = self.roi
        #bandfits = self.bandfits
        verbose  = self.verbose
        tolerance= self.tolerance
        l   = quadform.Localize(self,verbose = verbose)
        ld  = l.dir

        ll0 = self.spatialLikelihood(self.skydir)

        if not self.quiet:
            fmt ='Localizing source %s, tolerance=%.1e...\n\t'+7*'%10s'
            tup = (self.name, tolerance,)+tuple('moved delta ra     dec    a     b  qual'.split())
            print (fmt % tup)
            print (('\t'+4*'%10.4f')% (0,0,self.skydir.ra(), self.skydir.dec()))
            diff = np.degrees(l.dir.difference(self.skydir))
            print (('\t'+7*'%10.4f')% (diff,diff, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6]))
        
        old_sigma=1.0
        for i in xrange(self.max_iteration):
            try:
                l.fit(update=True)
            except:
                #raise
                l.recenter()
                if not self.quiet: print ('trying a recenter...')
                continue
            diff = np.degrees(l.dir.difference(ld))
            delt = np.degrees(l.dir.difference(self.skydir))
            sigma = l.par[3]
            if not self.quiet: print (('\t'+7*'%10.4f')% (diff, delt, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6]))
            if delt>self.maxdist:
                l.par[6]=99 # flag very bad quality and resect position
                l.sigma =1.0
                l.par[0]=self.skydir.ra(); l.par[1]=self.skydir.dec()
                if not self.quiet: print ('\t -attempt to move beyond maxdist=%.1f' % self.maxdist)
                break 
                #self.tsm.source.ellipse = self.qform.par[0:2]+self.qform.par[3:7]
                return False # hope this does not screw things up
                #raise Exception('localize failure: -attempt to move beyond maxdist=%.1f' % self.maxdist)
            if (diff < tolerance) and (abs(sigma-old_sigma) < tolerance):
                break # converge
            ld = l.dir
            old_sigma=sigma

        self.qform    = l
        self.lsigma   = l.sigma
        q = l.par
        self.ellipse = dict(ra=float(q[0]), dec=float(q[1]),
                a=float(q[3]), b=float(q[4]),
                ang=float(q[5]), qual=float(q[6]),
                lsigma = l.sigma)

        ll1 = self.spatialLikelihood(l.dir)
        if not self.quiet: print ('TS change: %.2f'%(2*(ll0 - ll1)))

        #roi.delta_loc_logl = (ll0 - ll1)
        # this is necessary in case the fit always fails.
        delt = np.degrees(l.dir.difference(self.skydir))
        self.delta_ts = 2*(ll0-ll1)
        self.delt = delt
        self.niter = i
        # if successful, add a list representing the ellipse to the source
        self.tsm.source.ellipse = self.qform.par[0:2]+self.qform.par[3:7] +[self.delta_ts] 
        return True #success
        
    def summary(self):
        if hasattr(self, 'niter') and self.niter>0: 
            print ('Localized %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \
                (self.name, self.niter, self.delt, self.delta_ts))
            labels = 'ra dec a b ang qual'.split()
            print ((len(labels)*'%10s') % tuple(labels))
            p = self.qform.par[0:2]+self.qform.par[3:7]
            print (len(p)*'%10.4f' % tuple(p))


       
def localize_all(roi, ignore_exception=True, **kwargs):
    """ localize all variable local sources in the roi, make TSmaps and associations if requested 
        ignore if extended -- has 'spatial_model'
        kwargs can have prefix to select subset with name starting with the prefix, e.g. 'SEED'
    """
    tsmin = kwargs.pop('tsmin',10)
    prefix = kwargs.pop('prefix', None)
    source_name = kwargs.pop('source_name', None)
    update = kwargs.pop('update', False)
    def filt(s):
        ok = s.skydir is not None\
            and isinstance(s, sources.PointSource) \
            and np.any(s.spectral_model.free)
        if not ok: return False
        if not hasattr(s,'ts'): 
            s.ts = roi.TS(s.name)
        return ok and s.ts>tsmin
    if source_name is not None:
        vpsources=[roi.get_source(source_name)]
    else:
        vpsources = filter(filt, roi.sources)
    tsmap_dir = kwargs.pop('tsmap_dir', None)
    if tsmap_dir is not None:
        if tsmap_dir[0]=='$':
            tsmap_dir = os.path.expandvars(tsmap_dir)
        if not os.path.exists(tsmap_dir):
            os.makedirs(tsmap_dir)
    associator = kwargs.pop('associator', None)
    tsfits = kwargs.pop('tsfits', True) 
    if len(kwargs.keys())>0:
        print ('Warning: unrecognized args to localize_all: %s' %kwargs)
    initw = roi.log_like()
    
    for source in vpsources:
        if prefix is not None and not source.name.startswith(prefix): continue
        
        full_localization(roi, source.name, ignore_exception=ignore_exception,
            update=update, associator=associator, tsmap_dir=tsmap_dir, tsfits=tsfits)
        

    curw= roi.log_like()
    if abs(initw-curw)>1.0 and not update:
        print ('localize_all: unexpected change in roi state after localization, from %.1f to %.1f (%+.1f)'\
            %(initw, curw, curw-initw))
        return False
    else: return True

class TS_function(object):
    """ usage:
        with TS_function(roi, 'test') as tsfun:
            print (tsfun(roi.skydir))
    """
    def __init__(self, roi, source_name):
        self.loc = Localization(roi, source_name)
    def __enter__(self):
        return self.loc.TSmap
    def __exit__(self, type, value, traceback):
        self.loc.reset()
    