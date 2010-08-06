"""
supplemental setup of ROI
----------------------------------
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/roi_setup.py,v 1.7 2010/06/20 13:39:30 burnett Exp $

These are near-duplicates of the classes with the same name in uw.like, but modifed for the interactive access

"""
import numpy as np
import os, pyfits, types, pickle

from uw.utilities import makerec
from uw.like import Models, roi_managers, pointspec_helpers
from skymaps import SkyDir, DiffuseFunction, IsotropicSpectrum

class ConsistentBackground(object):
    """A factory for constructing a  background model for a fit

   """

    def __init__(self, diffdir, quiet=False):
        fb = ('front','back')
        self.dmodels = ( 
           [ DiffuseFunction(os.path.join(diffdir,'gll_iem_v02_P6_v8_diff_%s.fits'%c)) for c in fb], 
           [ IsotropicSpectrum(os.path.join(diffdir,'isotropic_iem_%s_v02.txt'%c)) for c in fb],
           )
        self.names = ('gll_iem_v02', 'Isotropic Diffuse')

    def __call__(self, skydir=None, models = None, lat = None):
        """
        return a new set of background models
        """
        smodels = [ Models.PowerLaw(p=[1,1],free=[True,True],index_offset=1),
                    Models.Constant(free=[True]) ]

        # attempt to make index free depending on lat??
        smodels[0].free[1] = (lat is not None and abs(lat) < 20) 

        smodels = smodels if models is None else models
        return map(roi_managers.ROIBackgroundModel, self.dmodels, smodels, self.names)
      
class CatalogManager(object):
    """Read a  catalogue and use it to set up a source list for a ROI."""
    defaults = dict(
        prune_radius  = 0.10,   #deg; in a merge, consider sources closer than this duplicates
        free_radius   = 1,      #deg; sources within this distance have free spectral parameters
        min_flux      = 2e-9,   #ph/cm2/s; minimum flux for sources beyond a certain proximity
        max_distance  = 5,      #deg; distance inside which sources are returned regardless of flux
        min_ts        = 25,
        pulsar_dict   = None,   # optional pulsar dictionary with ExpCutoff parameters to use
        quiet         = False,
        verbose       = False,
        )

    def __init__(self,catalog_file,*args,**kwargs):
        self.__dict__.update(CatalogManager.defaults)
        for key, value in kwargs.items():
            if key in self.__dict__: self.__dict__[key]=value
            else:
                print 'warning: key %s not recognized by CatalogManager' % key
        #self.__dict__.update(kwargs)
        if not self.quiet: 
            print 'creating a CatalogManager'
        cdata = pyfits.open(catalog_file)[1].data
        try:
            ts   = np.asarray(cdata.field('Test_Statistic'),float)
        except KeyError:
            ts   = np.asarray(cdata.field('Signif_Avg'))**2
        good = ts>self.min_ts
        ras  = np.asarray(cdata.field('RA'),float)[good]
        decs = np.asarray(cdata.field('DEC'),float)[good]
        pens = cdata.field('PIVOT_ENERGY')[good]
        n0s  = cdata.field('FLUX_DENSITY')[good]
        inds = cdata.field('SPECTRAL_INDEX')[good]
        inds = np.where(inds > 0, inds, -inds)
        try:
            self.names  = np.asarray(cdata.field('Source_Name'))[good]
        except KeyError:
            self.names  = np.asarray(cdata.field('NickName'))[good]

        self.dirs   = map(SkyDir,ras,decs)
        
        pdkeys = [] 
        if self.pulsar_dict is not None:
            if type(self.pulsar_dict) == types.StringType:
                self.pulsar_dict = pickle.load(open(self.pulsar_dict))
            pdkeys = self.pulsar_dict.keys()
        def load_model(name, n0,ind, pen):
            #if name[:3]=='PSR': assert False, 'breakpoint'
            if name not in pdkeys:
                return Models.PowerLaw(p=[n0,ind],e0=pen)
            psr = self.pulsar_dict[name]
            if psr['TS']<100:return Models.PowerLaw(p=[n0,ind],e0=pen)
            stat = psr['stat'][0]
            if self.verbose: 
                print ('replacing catalog fits for source %s, par='+3*'%12.2e')\
                    % ((name,) + tuple(stat) )
            return Models.ExpCutoff(p=stat)
            
        self.models = np.asarray([load_model(name.strip(),n0,ind,pen) for name,n0,ind,pen in zip(self.names,n0s,inds,pens)])
        if not self.quiet: 
            print 'Loaded %d sources from catalog "%s" for roi backgrounds' % (len(cdata), catalog_file)
            if self.pulsar_dict is not None:
                print '\tUsing a pulsar dictionary with %d entries' % len(pdkeys)

    
    def append(self, acat):
        """ 
        append sources found in an auxilliary list of sources
        
        acat: a recarray with fields name, ra, dec, pnorm, pindex
        """
        names  = acat.name
        dirs   = map(SkyDir,acat.ra,acat.dec)
        models = [Models.PowerLaw(p=[n0,ind],e0=1000.) for n0,ind in zip(acat.pnorm,acat.pindex)]
        self.names = np.hstack((self.names, names))
        self.dirs= self.dirs+dirs
        self.models = np.hstack((self.models, models))
        if not self.quiet: print 'Added %d sources to roi background catalog, total now: %d' % (len(acat), len(self.names))
 
        
    def __call__(self,skydir, radius):

        diffs   = np.degrees(np.asarray([skydir.difference(d) for d in self.dirs]))
        #mask    = ((diffs < radius)&(self.ts > self.min_ts)) & \
        #          ((self.fluxes > self.min_flux)|(diffs < self.max_distance))         
        mask    = diffs < radius
        diffs   = diffs[mask]
        sorting = np.argsort(diffs)

        #sort SkyDirs -- pesky vector behavior...
        dirs    = [x for i,x in enumerate(self.dirs) if mask[i]]
        dirs    = [dirs[x] for x in sorting]

        names   = self.names[mask][sorting]
        models  = [x.copy() for x in self.models[mask][sorting]]

        fm = (diffs[sorting] < self.free_radius) 
        point_sources = map(pointspec_helpers.PointSource,dirs,names,models,fm)
        if not self.quiet:
            print '...selected %d sources within %.1f deg for roi' % (len(point_sources), radius)
            print '...selected %d sources within %.1f deg for refit' % (fm.sum(), self.free_radius)
        if diffs[sorting][0]< self.prune_radius: 
            if not self.quiet:
                print '...exclude catalog source %s  closer than %.1f' % (names[0], self.prune_radius)
            self.exclude = point_sources[0] 
            return point_sources[1:]
        self.exclude = None 
        return point_sources

    def merge_lists(self,skydir,radius=15,user_list=None):
        """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence.
         """

        cat_list = self(skydir,radius)
        if user_list is None: return cat_list

        from collections import deque
        merged_list = deque(user_list)

        for ncps,cps in enumerate(cat_list):
            merged_list.append(cps)
            for ups in user_list:
                if np.degrees(ups.skydir.difference(cps.skydir)) < self.prune_radius:
                   merged_list.pop(); break

        merged_list = list(merged_list)
        merged_list.sort(key = lambda ps: ps.skydir.difference(skydir))

        return merged_list
        
    def source_recarray(self):
        """ return a recarry of the full list of sources used for models
        """
        return np.rec.fromarrays(
                        [ self.names, 
                          [s.ra() for s in self.dirs], 
                          [s.dec() for s in self.dirs], 
                          [10**m.p[0] for m in self.models],
                          [10**m.p[1] for m in self.models],
                        ],
                        names = 'name ra dec pnorm pindex'.split(),
                       )

def UpdateRoiManager(object):

    """
    wrapper around the ROI class
    
    """
    def __init__(self, roi_manager):
        """
        roi_analysis
        """
        self.roi_manager = roi_manager
        self.ps_manager = roi_manager.ps
        
    def add_source(self, name, dir, model):
        rm = self.roi_manager
        rm.names.append(dir)
        rm.dirs.append(dir)
        rm.models.append(model)
        
    def add_source_array(self, recarray):
        """
        recarray: record array of sources to add
        """
        for s in recarray:
            add_source(s.name, SkyDir(s.ra,s.dec), Models.PowerLaw(p=(s.no,s.ind),e0=1000.))

if __name__=='__main__':
    pass
