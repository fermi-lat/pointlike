"""
source localization support

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/localization.py,v 1.25 2014/02/12 22:38:14 burnett Exp $

"""
import os
import numpy as np
from skymaps import SkyDir
from uw.like import quadform
from uw.utilities import keyword_options
from . import (sources, plotting )

    
class Localization(object):
    """ manage localization of a source
    Implements a minimization interface
    see also the localize function, which uses the eliptical fitter
   
    """
    defaults = (
        ('tolerance',1e-2),
        ('verbose',False),
        ('update',False,"Update the source position after localization"),
        ('max_iteration',10,"Number of iterations"),
        #('bandfits',True,"Default use bandfits"),
        ('maxdist',1,"fail if try to move further than this"),
        ('seedpos', None, 'if set, start from this position instead of the source position'),
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
    
    def log_like(self, skydir=None):
        """ return log likelihood at the given position"""
        return self.tsm(skydir)/2
   
    def TSmap(self, skydir):
        """ return the TS at given position, or 
            2x the log(likelihood ratio) from the nominal position
        """
        val= 2*(self.log_like(skydir)-self.maxlike)
        return val

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
            print fmt % tup
            print ('\t'+4*'%10.4f')% (0,0,self.skydir.ra(), self.skydir.dec())
            diff = np.degrees(l.dir.difference(self.skydir))
            print ('\t'+7*'%10.4f')% (diff,diff, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
        
        old_sigma=1.0
        for i in xrange(self.max_iteration):
            try:
                l.fit(update=True)
            except:
                #raise
                l.recenter()
                if not self.quiet: print 'trying a recenter...'
                continue
            diff = np.degrees(l.dir.difference(ld))
            delt = np.degrees(l.dir.difference(self.skydir))
            sigma = l.par[3]
            if not self.quiet: print ('\t'+7*'%10.4f')% (diff, delt, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
            if delt>self.maxdist:
                if not self.quiet: print '\t -attempt to move beyond maxdist=%.1f' % self.maxdist
                break # hope this does not screw things up
                #raise Exception('localize failure: -attempt to move beyond maxdist=%.1f' % self.maxdist)
            if (diff < tolerance) and (abs(sigma-old_sigma) < tolerance):
                break
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
        if not self.quiet: print 'TS change: %.2f'%(2*(ll0 - ll1))

        #roi.delta_loc_logl = (ll0 - ll1)
        # this is necessary in case the fit always fails.
        delt = np.degrees(l.dir.difference(self.skydir))
        self.delta_ts = 2*(ll0-ll1)
        self.delt = delt
        self.niter = i
        # if successful, add a list representing the ellipse to the source
        self.tsm.source.ellipse = self.qform.par[0:2]+self.qform.par[3:7] +[self.delta_ts] 
        
    def summary(self):
        if hasattr(self, 'niter') and self.niter>0: 
            print 'Localized %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \
                (self.name, self.niter, self.delt, self.delta_ts)
            labels = 'ra dec a b ang qual'.split()
            print (len(labels)*'%10s') % tuple(labels)
            p = self.qform.par[0:2]+self.qform.par[3:7]
            print len(p)*'%10.4f' % tuple(p)


       
def localize_all(roi, **kwargs):
    """ localize all variable local sources in the roi, make TSmaps and associations if requested 
        ignore if extended -- has 'spatial_model'
        kwargs can have prefix to select subset with name starting with the prefix, e.g. 'SEED'
    """
    tsmin = kwargs.pop('tsmin',10)
    prefix = kwargs.pop('prefix', None)
    source_name = kwargs.pop('source_name', None)
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
        print 'Warning: unrecognized args to localize_all: %s' %kwargs
    initw = roi.log_like()
    
    for source in vpsources:
        if prefix is not None and not source.name.startswith(prefix): continue
        with roi.tsmap_view(source.name) as tsm:
            
            loc = Localization(tsm)
            try:
                loc.localize()
            except Exception, msg:
                print 'Localization of %s failed: %s' % (source.name, msg)
            #source.ellipse = loc.qform.par[0:2]+loc.qform.par[3:7] +[loc.delta_ts] if hasattr(loc,'qform') else None
            if not roi.quiet and hasattr(loc, 'niter') and loc.niter>0: 
                print 'Localized %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \
                    (source.name, loc.niter, loc.delt, loc.delta_ts)
                labels = 'ra dec a b ang qual'.split()
                print (len(labels)*'%10s') % tuple(labels)
                p = loc.qform.par[0:2]+loc.qform.par[3:7]
                print len(p)*'%10.4f' % tuple(p)
            if associator is not None:
                make_association(source, loc.TSmap, associator, quiet=roi.quiet)
            
            if tsmap_dir is not None : 
                if  hasattr(loc,'ellipse'): 
                    a, qual, delta_ts = loc.ellipse['a'], loc.ellipse['qual'], loc.delta_ts
                    tsize = min(a*15., 2.0)
                    bad = a>0.25 or qual>5 or abs(delta_ts)>2
                    if bad:
                        print 'Flagged as possibly bad: a>0.25 or qual>5 or abs(delta_ts)>2:', a, qual, delta_ts
                else: 
                    print 'no localization'
                    bad = True
                    tsize= 2.0
                if tsmap_dir.endswith('fail') and not bad: continue
                pixelsize= tsize/15.;
                #try:

                tsm=plotting.tsmap.plot(loc, source.name, center=tsm.saved_skydir,
                    outdir=tsmap_dir, catsig=0, size=tsize, 
                    pixelsize= pixelsize, # was 14: desire to have central pixel
                    # todo: fix this
                    assoc=source.__dict__.get('adict', None), # either None or a dictionary
                    notitle=True, #don't do title
                    markersize=10,
                    primary_markersize=12,
                    tsfits=tsfits,
                    )
                #except Exception, msg:
                #    print 'Plot of %s failed: %s' % (source.name, msg)
                #    raise
    curw= roi.log_like()
    if abs(initw-curw)>1.0:
        print 'localize_all: unexpected change in roi state after localization, from %.1f to %.1f (%+.1f)'\
            %(initw, curw, curw-initw)
        return False
    else: return True

class TS_function(object):
    """ usage:
        with TS_function(roi, 'test') as tsfun:
            print tsfun(roi.skydir)
    """
    def __init__(self, roi, source_name):
        self.loc = Localization(roi, source_name)
    def __enter__(self):
        return self.loc.TSmap
    def __exit__(self, type, value, traceback):
        self.loc.reset()
    