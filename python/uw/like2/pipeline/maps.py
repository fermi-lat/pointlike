"""
Code to generate a set of maps for each ROI
$Header$

"""
import os, sys,  pickle, types
import numpy as np
from uw.like import Models, roi_tsmap
from uw.utilities import fitter
from skymaps import Band, SkyDir, PySkyFunction, Hep3Vector 

# convenience adapters for ResidualTS model
def LogParabola(*pars):return Models.LogParabola(p=pars)
def PowerLaw(*pars):   return Models.PowerLaw(p=pars)
def ExpCutoff(*pars):  return Models.ExpCutoff(p=pars)

class CountsMap(dict):
    """ A map with counts per HEALPix bin """
    
    def __init__(self, roi, emin=1000., nside=512):
        """ roi : a region of interest object
            emin : float
                minimum energy for constructing the counts map
            nside : int
                the nside parameter for the HEALPix map
        Note that is is implemented as a dictionary, with the pixel index as key
        """
        self.indexfun = Band(nside).index
        self.emin=emin
        self._get_photons(roi, emin)
        
    def _get_photons(self, roi, emin):
        for band in roi.bands:
            e = band.e
            if e<emin: continue
            t = band.b.event_class() & 3 # note: mask off high bits!
            for wsd in band.wsdl:
                index, wt = self.indexfun(wsd), int(wsd.weight())
                if index in self: self[index]+=wt
                else: self[index] = wt
    def __call__(self,v):
        return self.get(self.indexfun(v), 0)
    
class KdeMap(dict):
    """ Implement a SkyFunction that returns KDE data for a given ROI"""
    
    def __init__(self, roi, nside=12, factor=32, **kwargs):
        """ roi: an ROIAnalysis object
            nside: int 
                HEALPix nside
            factor: int
                factor defining nside of subpixel
        """
        self.nside = nside
        # invoke code in like/roi_tsmap that will create the full map
        hr = roi_tsmap.ROIWrapper(roi,nside)
        self.kde = roi_tsmap.HealpixKDEMap(hr,factor=factor, **kwargs)
                
    def __call__(self, v, isskydir=False):
        return self.kde.call2(v)

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
        photon_index=kwargs.pop('index', 2.2) 
        self.source =roi.add_source(name='tsmap', skydir = roi.roi_dir, model=model)
        self.roi.select_source('tsmap')
        self.index = len(self.roi.get_parameters())-2
        self.model = self.source.spectral_model
        
    def reset(self):
        #self.roi.del_source('test')
        self.roi.select_source(None)
        
    def tsfun(self, skydir):
        self.source.skydir = skydir
        self.roi.calls =0
        self.roi.update(reset=True)
        fn = fitter.Projector(self.roi, select=[self.index])  
        mm = fitter.Minimizer(fn, quiet=True)
        mm(use_gradient=True, estimate_errors=False)
        self.model[0]=1e-15
        llzero = self.roi.log_like()
        ts= 2*(-mm.fitvalue-llzero)
        return max(ts, 0)
        
    def __call__(self, v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        return self.tsfun(skydir)
        

class ROItables(object):
    """ manage one or more tables of values subdividing a HEALpix roi
    
        skyfuns : list of 3-tuples
            the 3-tuples must have: (skyfunction, tablename, dict)
            Implemented now, and default here:
                (ResidualTS,'ts',  dict(photon_index=2.0),) , 
                (KdeMap,    'kde', dict()),
            If skyfunction is a string, evaluate it 
    """
    
    def __init__(self, outdir, **kwargs):
        nside = kwargs.pop('nside', 12)
        subnside = kwargs.pop('subnside', 512)
        self.make_index_table(nside,subnside)
        self.subdirfun = Band(subnside).dir
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

def countmaps(g, outdir=None):
    """ generate all the roi tables for CountsMap """
    if outdir is None: outdir = g.process_kw['outdir']
    if not os.path.exists(outdir): os.mkdir(outdir)
    table = ROItables(outdir, skyfuns= [["CountsMap", "counts", dict()]])
    quiet, g.quiet = g.quiet, True
    for index in range(0,1728):
        r = g.roi(index, roi_kw=dict(skip_setup=True)) 
        table(r)
    g.quiet=quiet    