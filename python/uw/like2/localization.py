"""
source localization support

$Header$

"""
import os
import numpy as np
from skymaps import SkyDir
from uw.like import quadform, srcid
from uw.utilities import keyword_options
from . plotting import tsmap


def ufunc_decorator(f): # this adapts a bound function
    def new_ufunc(self, par):
        return np.array(map(lambda x: f(self,x),par)) if hasattr(par, '__iter__')  else f(self,par)
    return new_ufunc


    
class Localization(object):
    """ manage localization of a source
    Implements a minimization interface
    see also the localize function, which uses the eliptical fitter
   
    Note that it supports the context manager interface, and should be used in a with construction
    to guarantee that the reset function is called to restore the ROI state.
    
    """
    defaults = (
        ('tolerance',1e-3),
        ('verbose',False),
        ('update',False,"Update the source position after localization"),
        ('max_iteration',10,"Number of iterations"),
        #('bandfits',True,"Default use bandfits"),
        ('maxdist',1,"fail if try to move further than this"),
        ('quiet', False, 'set to suppress output'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roistat, source_name, **kwargs):
        """ roistat : an ROIstat object
            source_name : string
                the name of a source, with possible wild cards
        """
        keyword_options.process(self, kwargs)
        self.rs = roistat
        self.source = self.rs.sources.find_source(source_name)
        self.rs.select_source(self.source.name)
        self.rs.update(True)
        self.maxlike=self.rs.log_like()
        self.skydir=self.saved_skydir = self.source.skydir #saved value
        self.name = self.source.name
    
    def __enter__(self):
        """ supports the 'with' construction, guarantees that reset is called to restore the ROI
        example:
        -------
        with Localize(roi, name) as loc
            # use loc.TSmap ...
        """
        return self
        
    def __exit__(self, type, value, traceback):
        self.reset()

    def log_like(self, skydir):
        """ return log likelihood at the given position"""
        self.source.skydir =skydir
        self.rs.update(True)
        w = self.rs.log_like()
        self.source.skydir = self.saved_skydir
        return w
   
    def TSmap(self, skydir):
        """ return the TS at given position, or 
            2x the log(likelihood ratio) from the nominal position
        """
        val= 2*(self.log_like(skydir)-self.maxlike)
        return val

    
    def get_parameters(self):
        return np.array([self.source.skydir.ra(), self.source.skydir.dec()])
    
    def set_parameters(self, par):
        self.skydir = SkyDir(par[0],par[1])
        self.source.skydir = self.skydir
        
    def __call__(self, par):
        return -self.TSmap(SkyDir(par[0],par[1]))
    
    def reset(self):
        """ restore modifications to the ROIstat
        """
        self.source.skydir=self.skydir
        self.rs.select_source(None)
      
    def dir(self):
        return self.skydir

    def errorCircle(self):
        return 0.05 #initial guess

    def spatialLikelihood(self, sd, update=False):
        return -self.log_like(sd)
        
    def localize(self):
        """Localize a source using an elliptic approximation to the likelihood surface.

            return fit position, number of iterations, distance moved, delta TS
        """
        #roi    = self.roi
        #bandfits = self.bandfits
        verbose  = self.verbose
        update    = self.update
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
                raise Exception('localize failure: -attempt to move beyond maxdist=%.1f' % self.maxdist)
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
        if self.update:
            self.source.skydir = l.dir
        self.delta_ts = 2*(ll0-ll1)
        self.delt = delt
        self.niter = i

        
  
def make_association(source, tsf, associate):
    print ' %s association(s) ' % source.name,
    try:    ell = source.ellipse
    except: ell = None
    if ell is None:
        print '...no localization'
        source.adict = None
        return
    assert len(ell)>6, 'invalid ellipse for source %s' % source.name
    try:
        adict = associate(source.name, SkyDir(ell[0],ell[1]), ell[2:5]) 
    except srcid.SrcidError, msg:
        print 'Association error for %s: %s' % (source.name, msg)
        adict=None
    except Exception, msg:
        print 'Exception associating %s: %s' %( source.name, msg)
        adict=none
    source.adict = adict 
    if adict is not None:
    
        ts_local_max=tsf( SkyDir(ell[0],ell[1]) )
        adict['deltats'] = [ts_local_max-tsf(d) for d in adict['dir']]
        print '\n   cat         name                  ra        dec         ang     prob    Delta TS'
        #       15 Mrk 501               253.4897   39.7527    0.0013      0.41
        fmt = '   %-11s %-20s%10.4f%10.4f%10.4f%8.2f%8.1f' 
        for i,id_name in enumerate(adict['name']):
            tup = (adict['cat'][i], id_name, adict['ra'][i], adict['dec'][i], adict['ang'][i], 
                    adict['prob'][i],adict['deltats'][i])
            print fmt % tup
    else:
        print '...None  found'
      

def localize_all(roi, **kwargs):
    """ localize all variable local sources in the roi, make TSmaps and associations if requested 
    """
    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    tsmap_dir = kwargs.pop('tsmap_dir', None)
    associator = kwargs.pop('associator', None)
    tsfits = kwargs.pop('tsfits', False) #TODO: reimplement this to generate FITS maps
    
    initw = roi.log_like()
    
    for source in sources:
        with Localization(roi, source.name, quiet=True, **kwargs) as loc:
            try:
                loc.localize()
            except Exception, msg:
                print 'Localization of %s failed: %s' % (source.name, msg)
                continue
            source.ellipse = loc.qform.par[0:2]+loc.qform.par[3:7] +[loc.delta_ts] if hasattr(loc,'qform') else None
            if not roi.quiet and hasattr(loc, 'niter') and loc.niter>0: 
                print 'Localized %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \
                    (source.name, loc.niter, loc.delt, loc.delta_ts)
                labels = 'ra dec a b ang qual'.split()
                print (len(labels)*'%10s') % tuple(labels)
                p = loc.qform.par[0:2]+loc.qform.par[3:7]
                print len(p)*'%10.4f' % tuple(p)
            if associator is not None:
                make_association(source, loc.TSmap, associator)
            
            if tsmap_dir is not None:
                tsize = loc.ellipse['a']*15. if hasattr(loc,'ellipse') and loc.ellipse is not None else 1.1
                pixelsize= tsize/15.;
                tsm=tsmap.plot(loc, source.name, center=source.skydir, 
                    outdir=tsmap_dir, catsig=0, size=tsize, 
                    pixelsize= pixelsize, # was 14: desire to have central pixel
                    # todo: fix this
                    assoc=source.__dict__.get('adict', None), # either None or a dictionary
                    notitle=True, #don't do title
                    markersize=10,
                    primary_markersize=12,
                    )
    curw= roi.log_like()
    if abs(initw-curw)>1.0:
        print 'localize_all: unexpected change in roi state after localization, from %.1f to %.1f' %(initw, curw)
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
    