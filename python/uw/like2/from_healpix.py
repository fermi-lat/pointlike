"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/from_healpix.py,v 1.10 2015/07/24 17:57:06 burnett Exp $
"""
import os, pickle, zipfile
import numpy as np
import skymaps
from pointlike import IntVector

from . import (roimodel, sources, diffuse, extended)

def neighbors(index,  rings=1):
    """ return the cluster of pixel indeces around the pixel 
    Parameters
    ----------
    index : int
        specified pixel
    rings : int, optional
        number of rings around the pixel
        
    Returns a list of (index, ring_number) tuples
    """
    b12 = skymaps.Band(12)
    v = IntVector()
    outer_ring= set([index])
    cluster = set([index])
    ret = []
    for ring_number in range(rings):
        found = set([])
        for i in outer_ring:
            b12.findNeighbors(int(i),v)
            found = found.union(set(v))
        outer_ring = found.difference(cluster)
        ret = ret + [(x,ring_number+1) for x in outer_ring]
        cluster = cluster.union(found)
    return ret

class ROImodelFromHealpix(roimodel.ROImodel):
    
    def load_sources( self, roi_index, rings=1):
        """ load sources from the roi and its neighbors in the pickle file found in modeldir
        
        roi_index : integer
            HEALPix index of the ROI
        rings : integer
            number of rings of concentric pixels to search for (fixed) sources to add
            Special value: if -1, do not add any sources at all
        """

        self.pickle_file = os.path.join(self.config.modeldir, 'pickle.zip')
        if not os.path.exists(self.pickle_file):
            raise Exception('Expected file "pickle.zip" not found in %s' % config.configdir)
        if hasattr(roi_index, '__len__'):
            roi_index = skymaps.Band(12).dir(skymaps.SkyDir(*roi_index))
        self.index = roi_index
        self.roi_dir = skymaps.Band(12).dir(roi_index)
        self.name = 'HP12_%04d' % roi_index
        self._z = zipfile.ZipFile(os.path.expandvars(self.pickle_file))
        global_only = rings==-1
        self.load_sources_from_healpix((roi_index,0), global_only=global_only)
        if global_only: return
        for neighbor_index in neighbors(roi_index, rings=rings):
            self.load_sources_from_healpix(neighbor_index, neighbors=True)

        
    def load_sources_from_healpix(self, index, neighbors=False, global_only=False):
        """ select and add sources in the given HEALPix to self.
        Tags each with an index property, which is a tuple (healpix_index, ring number)
        if neigbors is True, add only local sources.
        also set the free list to False for neighbors. (The Source constructor sets it as a 
        copy of the model's free list)
        """

        def load_local_source(name, rec):
            if not rec['isextended']:
                src = sources.PointSource(name=name, skydir=rec['skydir'], 
                    model=rec['model'],
                    ellipse = rec.get('ellipse', None),
                    associations = rec.get('associations', None),
                    band_ts = rec.get('band_ts', None),
                    sedrec = rec.get('sedrec', None),
                    ts = rec.get('ts', None),
                    pivot_energy = rec.get('pivot_energy', None),
                    fixed_spectrum = rec.get('fixed_spectrum', False),
                    )
                if src.fixed_spectrum:
                    assert sum(src.model.free)==1, \
                        'Logic error? model=%s for %s' % (src.model, src.name)
            else:
                # extended source: get from extended catalog list, replace model, add info from previous fit
                src = self.ecat.lookup(name)
                src.model = rec['model']
                # don't know why this is necessary
                def tuple_check( attr):
                    t = rec.get(attr, None)
                    if isinstance(t,tuple): t = t[0]
                    src.__dict__[attr] = t
                map(tuple_check, ('ts', 'band_ts', 'sedrec', 'pivot_energy'))
                
            #if neighbors: src.free[:]=False # not sure this is necessary
            if neighbors: src.model.free[:]=False
            src.index = index
            #src.model.free[:len(src.free)] = src.free # Don't think I still need this second copy of free
            #if src.model.name=='LogParabola':
            #    src.model.free[-1] = False
            #    src.free[-1] = False    

            return src
        
        def load_global_source(name, rec):

            if not self.quiet: print 'Loading global source %s for %d' % (name, index)
            df = self.config.diffuse[name]
            if df is None: 
                return None
            gsrc = sources.GlobalSource(name=name, skydir=None, model=rec,
                dmodel = diffuse.diffuse_factory( df, self.config.event_type_names))
            gsrc.index =index
            return gsrc


        self.pickle_file = 'pickle/HP12_%04d.pickle' % index[0]
        p = pickle.load(self._z.open(self.pickle_file))
        if not neighbors:
            self.prev_logl = p.get('prev_logl', []) # make history available
            self.history = p.get('history', []) # will manage history of likelihood, stage, time, stream, cpu time
            global_sources = [load_global_source(name, rec) for name, rec \
                in zip(p['diffuse_names'], p['diffuse']) if name not in self.ecat.names]
            self.global_count = len(global_sources)
            for s in global_sources:
                if s is not None: self.append(s)
            if global_only: return
        local_sources = [load_local_source(name, rec) for name,rec in p['sources'].items()]
        if not neighbors: self.local_count = len(local_sources)
        for s in local_sources: self.append(s)
        

    def __repr__(self):
        return '%s.%s : %d global, %d local, %d total sources for ROI %d' \
            % (self.__module__, self.__class__.__name__, self.global_count,self.local_count,  
            len(self), self.index)
            
