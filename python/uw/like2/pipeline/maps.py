"""
Code to generate a set of maps for each ROI
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/maps.py,v 1.10 2013/02/21 17:59:10 burnett Exp $

"""
import os, sys,  pickle, types
import numpy as np
from uw.like import Models
from uw.utilities import fitter
from skymaps import Band, SkyDir, PySkyFunction, Hep3Vector, PythonUtilities 
from uw.like2 import sourcelist, sedfuns, loglikelihood

# convenience adapters for ResidualTS model
def LogParabola(*pars):
    model = Models.LogParabola(p=pars)
    sourcelist.set_default_bounds(model)
    return model
def PowerLaw(*pars):   
    model = Models.PowerLaw(p=pars)
    sourcelist.set_default_bounds(model)
    return model
def ExpCutoff(*pars):  
    model = Models.ExpCutoff(p=pars)
    sourcelist.set_default_bounds(model)
    return model
def PLSuperExpCutoff(*pars):  
    model = Models.PLSuperExpCutoff(p=pars)
    sourcelist.set_default_bounds(model)
    return model

def make_index_table(nside, subnside, usefile=True):
    filename = 'index_table_%02d_%03d.pickle' % (nside, subnside)
    if os.path.exists(filename) and usefile:
        return pickle.load(open(filename))
    print 'generating index table for nside, subnside= %d %d' % (nside, subnside)
    band, subband = Band(nside), Band(subnside)
    npix, nsubpix = 12*nside**2, 12*subnside**2
    t=np.array([band.index(subband.dir(i)) for i in xrange(nsubpix)])
    a = np.arange(nsubpix)
    index_table = [a[t==i] for i in xrange(npix)]
    if usefile:
        pickle.dump(index_table, open(filename,'w'))
    return index_table
    

class CountsMap(dict):
    """ A map with counts per HEALPix bin """
    
    def __init__(self, roi, emin=1000., nside=256):
        """ roi : a region of interest object
            emin : float
                minimum energy for constructing the counts map
            nside : int
                the nside parameter for the HEALPix map
        Note that is is implemented as a dictionary, with the pixel index as key
        """
        self.indexfun = Band(nside).index
        self.nside = nside
        self._fill(roi, emin)
        
    def _fill(self, roi, emin):
        for bandlike in roi.selected_bands:
            band = bandlike.band
            if band.e < emin: continue
            #t = band.b.event_class() & 3 # note: mask off high bits!
            for wsd in band.wsdl:
                index, wt = self.indexfun(wsd), int(wsd.weight())
                if index in self: self[index]+=wt
                else: self[index] = wt
                
    def __call__(self,v):
        return self.get(self.indexfun(v), 0)
  
class ModelMap(CountsMap):
    """ A map with model prediction per HEALPix bin
    """
        
    def _fill(self, roi, emin):
        index_table = make_index_table(12, self.nside)
        roi_index = Band(12).index(roi.roi_dir)
        ids = index_table[roi_index]
        sdirs = [Band(self.nside).dir(int(i)) for i in index_table[roi_index]]
        grid = np.zeros(len(ids))
        
        for bandlike in roi.selected_bands:
            band = bandlike.band
            if band.e < emin: continue
            grid += bandlike.fill_grid(sdirs)
        solidangle = 4.*np.pi/(12*self.nside**2)
        grid *= solidangle
        for i,v in zip(ids, grid):
            self[i]=v
    

class KdeMap(object):
    """ Implement a SkyFunction that returns KDE data for a given ROI"""
    
    def __init__(self, roi, emin=[500,1000], **kwargs):
        """ roi: an ROIstat object
            emin: list of two minimum energies
        """
        # create local array of references corresponding to bands from ROI 
        self.bands = []
        self.r95 = []
        for  bandlike in roi.selected_bands:
            band = bandlike.band 
            if band.emin < emin[band.ct & 3]: continue
            self.r95.append( np.radians(band.psf.inverse_integral_on_axis(0.95)))
            self.bands.append(band)
    
    def __call__(self,v,skydir = None):
        """ copied from roi_tsmap.HealpixKDEMap """
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        rval = 0
        for i,band in enumerate(self.bands):
            if not band.has_pixels: continue
            rvals = np.empty(len(band.wsdl),dtype=float)
            PythonUtilities.arclength(rvals, band.wsdl, sd)
            mask = rvals < self.r95[i]
            rval +=(band.psf(rvals[mask],density=True)*band.pix_counts[mask]).sum()
        return rval


class ResidualTS(object):
    """ manage a residual TS plot 
    """

    def __init__(self, roi, **kwargs):
        """
        roi : a ROI_user object
        kwargs:
            index : float
                default 2.2 the photon index to use
            model : None, or a Models.model instance, or a string like 
                    'LogParabola(p=[6e-14, 1.44, 1.22, 4500])'
        """
        self.roi = roi
        model = kwargs.pop('model', None)
        if type(model)==types.StringType:
            print 'ResidualTS: using spectral model: %s' %model
            model = eval(model)
        self.sourcename=kwargs.pop('sourcename', 'tsmap')
        self.source =roi.add_source(name=self.sourcename, skydir = roi.roi_dir, model=model)
        self.roi.select_source(self.sourcename)
        self.index = list(roi.parameter_names).index('tsmap_Norm')
        self.model = self.source.spectral_model
        sourcelist.set_default_bounds(self.model) # in case no bounds already
        roi.summary()
        roi.fit([self.index]) #check
        print 'TS check:', roi.TS('tsmap')
        
    def __enter__(self):
        roi.summary()
        roi.fit([self.index]) #check
        print 'TS check:', roi.TS('tsmap')

        return self
    def __exit__(self,*pars):
        self.reset()
        
    def reset(self):
        #self.roi.del_source('test')
        self.roi.select_source(None)
        
    def tsfun(self, skydir):
        self.source.skydir = skydir
        self.roi.calls =0
        self.roi.update(reset=True)
        self.model[0]=1e-14 # initial value above limit
        fn = fitter.Projector(self.roi, select=[self.index])  
        mm = fitter.Minimizer(fn, quiet=True)
        mm(use_gradient=True, estimate_errors=False, use_bounds=False)
        self.model[0]=1e-15 # for TS computation
        llzero = self.roi.log_like()
        ts= 2*(-mm.fitvalue-llzero)
        return max(ts, 0)
        
    def __call__(self, v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        return self.tsfun(skydir)
        
class ResidualLikelihood(ResidualTS):
    """ save the likelihood function, as the 3-parameter representation of a shifted Poisson, plss the max dev.
    """
    def tsfun(self, skydir):
        self.source.skydir = skydir
        self.roi.update(True)
        #print 'Studying source %s at %s' % (self.sourcename, skydir) ,

        with sedfuns.SourceFlux(self.roi, self.sourcename) as sf:
            sf.select_band(None)
            try:
                pf = loglikelihood.PoissonFitter(sf, tol=1.0)
                p = pf.poiss.p+ [pf.maxdev]
                #print 'TS, maxdev: %.2f %.2f' % (pf.poiss.ts, pf.maxdev)
            except Exception, msg:
                p = [0,0,0, 1.0]
                print 'Failed at %s: %s' % (skydir,msg)
        return p 


class ROItables(object):
    """ manage one or more tables of values subdividing a HEALpix roi
    
        skyfuns : list of 3-tuples
            the 3-tuples must have: (skyfunction, tablename, dict)
            Implemented now, and default here:
                (ResidualTS,'ts',  dict(photon_index=2.0),) , 
                (KdeMap,    'kde', dict()),
            If skyfunction is a string, evaluate it 
    """
    
    def __init__(self, outdir, nside, roi_nside=12, **kwargs):
        print 'creating index table for nside %d' % nside
        self.make_index_table(roi_nside, nside)
        self.subdirfun = Band(nside).dir
        self.skyfuns = kwargs.pop('skyfuns', 
              ( (ResidualTS, 'ts', dict(photon_index=2.2),) , 
                #(KdeMap,    'kde', dict()),
              ),
            )
        self.subdirs = [os.path.join(outdir, name+'_table') for s, name, kw in self.skyfuns]
        for subdir in self.subdirs: 
            if not os.path.exists(subdir):  os.makedirs(subdir)
            
    def make_index_table(self, nside, subnside):
        band, subband = Band(nside), Band(subnside)
        npix, nsubpix = 12*nside**2, 12*subnside**2
        t=np.array([band.index(subband.dir(i)) for i in xrange(nsubpix)])
        a = np.arange(nsubpix)
        self.index_table = [a[t==i] for i in xrange(npix)]
        
    def process_table(self, skyfun,name, pos_list, outfile=None, **kwargs):
        print 'filling table %-5s with %d entries...' % (name, len(pos_list)),
        sys.stdout.flush()
        skytable = np.array([skyfun(p) for p in pos_list])
        print ' min=%6.1f, max=%6.1f, mean=%6.1f ' \
            % (skytable.min(), skytable.max(),skytable.mean()) ,
        if outfile is not None:
            print '--> %s' % outfile
            pickle.dump(skytable, open(outfile,'wb'))
        else: print  
        if hasattr(skyfun,'reset'): skyfun.reset() 
  
    
    def __call__(self, roi):
        index = int(roi.name[5:])
        pos_list = [self.subdirfun(int(i)) for i in self.index_table[index]]
          
        #if not hasattr(roi, 'bands'):
        #    roi.bands = [s.band for s in roi.selected_bands]
        #    roi.phase_factor = 1.0 
            
        for i,fun in enumerate(self.skyfuns):
            skyfun = fun[0] if type(fun[0])!=types.StringType else eval(fun[0])
            self.process_table(skyfun(roi, **fun[2]), fun[1], pos_list, 
                os.path.join(self.subdirs[i], roi.name+'.pickle'))

def countmaps(g, outdir=None, nside=256, dom=range(1728), emin=1000):
    """ generate all the roi tables for CountsMap """
    if outdir is None: outdir = g.process_kw['outdir']
    if not os.path.exists(outdir): os.mkdir(outdir)
    table = ROItables(outdir, nside, skyfuns= [[CountsMap, "counts", dict(emin=emin)]])
    quiet, g.quiet = g.quiet, True
    map(lambda i: table(g.roi(i)), dom)
    g.quiet=quiet    
    
def kdemaps(g, outdir=None, nside=256, dom=range(1728)):
    """ generate all the roi tables for KdeMap """
    if outdir is None: outdir = g.process_kw['outdir']
    if not os.path.exists(outdir): os.mkdir(outdir)
    table = ROItables(outdir, nside, skyfuns= [[KdeMap, "kde", dict()]])
    quiet, g.quiet = g.quiet, True
    map(lambda i: table(g.roi(i)), dom)
    g.quiet=quiet    

def modelmaps(g, outdir=None, nside=256, dom=range(1728)):
    """ generate all the roi tables for KdeMap """
    if outdir is None: outdir = g.process_kw['outdir']
    if not os.path.exists(outdir): os.mkdir(outdir)
    table = ROItables(outdir, nside, skyfuns= [[ModelMap, "model", dict()]])
    quiet, g.quiet = g.quiet, True
    map(lambda i: table(g.roi(i)), dom)
    g.quiet=quiet    