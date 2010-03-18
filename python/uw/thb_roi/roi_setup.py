"""
supplemental setup of ROI
----------------------------------
$Header$

"""
import numpy as np
import os, pickle, math, glob, pyfits

from uw.utilities import makerec
from uw.like import Models, roi_managers, pointspec_helpers
from skymaps import SkyDir, DiffuseFunction, IsotropicSpectrum

class ConsistentBackground(object):
    """Manage the construction of a consistent background model.

      *Notes*
      
      A background model comprises diffuse and point sources.  While
      some parameters may be refit in a spectral analysis, the default
      background model should give a good representation.
      
      The class handles matching LAT source lists with the diffuse models
      used to generate them.  The combination is a "consistent"
      background.

    """

    def __init__(self,analysis_environment, background='1FGL', **kwargs):

        self.diffdir    = analysis_environment.diffdir
        self.catdir     = analysis_environment.catdir
        self.background = background

        eval('self.make_%s()'%background)
        ae=analysis_environment
        self.cat = os.path.join(ae.catdir, ae.catalog)
        self.cm = CatalogManager(self.cat, **kwargs)


    def make_1FGL(self):

        gf  = DiffuseFunction(os.path.join(self.diffdir,'gll_iem_v02_P6_v3_diff_front.fits'))
        gb  = DiffuseFunction(os.path.join(self.diffdir,'gll_iem_v02_P6_v3_diff_back.fits'))
        isof= IsotropicSpectrum(os.path.join(self.diffdir,'isotropic_iem_front_v02.txt'))
        isob= IsotropicSpectrum(os.path.join(self.diffdir,'isotropic_iem_back_v02.txt'))

        self.dmodels = [ [gf,gb], [isof,isob] ]
        self.smodels = [ Models.PowerLaw(p=[1,1],free=[True,True],index_offset=1),
                       Models.Constant(free=[True]) ]
        self.names   = ['gll_iem_v02', 'Isotropic Diffuse']


    def get_bgmodels(self, models = None, lat = None):

        gal_index = True if (lat is not None and abs(lat) < 20) else False
        self.smodels[0].free[1] = gal_index # need to be careful about this
        smodels = self.smodels if models is None else models
        return map(roi_managers.ROIBackgroundModel,self.dmodels,smodels,self.names)
      
class CatalogManager(object):
    """Read a  catalogue and use it to set up a source list for a ROI."""

    def init(self):
        self.prune_radius  = 0.10 #deg; in a merge, consider sources closer than this duplicates
        self.free_radius   = 1 #deg; sources within this distance have free spectral parameters
        self.min_flux      = 2e-9 #ph/cm2/s; minimum flux for sources beyond a certain proximity
        self.max_distance  = 5 #deg; distance inside which sources are returned regardless of flux
        self.min_ts        = 25
        self.quiet          =False

    def __init__(self,catalog_file,*args,**kwargs):
        self.init()
        self.__dict__.update(kwargs)
        cdata = pyfits.open(catalog_file)[1].data
        ras  = np.asarray(cdata.field('RA'),float)
        decs = np.asarray(cdata.field('DEC'),float)
        pens = cdata.field('PIVOT_ENERGY')
        n0s  = cdata.field('FLUX_DENSITY')
        inds = cdata.field('SPECTRAL_INDEX')
        inds = np.where(inds > 0, inds, -inds)

        self.dirs   = map(SkyDir,ras,decs)
        self.models = np.asarray([Models.PowerLaw(p=[n0,ind],e0=pen) for n0,ind,pen in zip(n0s,inds,pens)])
        self.names  = np.asarray(cdata.field('Source_Name'))
        if not self.quiet: print 'Loaded %d sources from catalog "%s" for roi backgrounds' % (len(cdata), catalog_file)
    
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
 
        
    def get_sources(self,skydir,radius):

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

        fm = diffs[sorting] < self.free_radius
        point_sources = map(pointspec_helpers.PointSource,dirs,names,models,fm)
        if True: #not self.quiet:
            print '...selected %d sources within %.1f deg for roi' % (len(point_sources), radius)
            print '...selected %d sources within %.1f deg for refit' % (fm.sum(), self.free_radius)
        return point_sources

    def merge_lists(self,skydir,radius=15,user_list=None):
        """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence.
         """

        cat_list = self.get_sources(skydir,radius)
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
