"""
$Header$
"""
import os, pickle, zipfile
import numpy as np
import skymaps
from . import (roimodel, sources, diffuse, extended)

class ROImodelFromHealpix(roimodel.ROImodel):
    
    def load_sources( self, roi_index):
        # load sources from the roi and its neighbors in the pickle file found in modeldir
        self.pickle_file = os.path.join(self.config.modeldir, 'pickle.zip')
        if not os.path.exists(self.pickle_file):
            raise Exception('Expected file "pickle.zip" not found in %s' % config.configdir)
        self.index = roi_index
        self._z = zipfile.ZipFile(os.path.expandvars(self.pickle_file))
        self.load_sources_from_healpix(roi_index)
        for neighbor_index in self.neighbors():
            self.load_sources_from_healpix(neighbor_index, neighbors=True)

        
    def load_sources_from_healpix(self, index, neighbors=False):
        """ select and add sources in the given HEALPix to self.
        Tags each with an index property
        if neigbors is True, add only local sources.
        also set the free list to False for neighbors. (The Source constructor sets it as a 
        copy of the model's free list)
        """

        def load_local_source(name, rec):
            if not rec['isextended']:
                src = sources.PointSource(name=name, skydir=rec['skydir'], 
                    model=rec['model'],
                    ellipse = rec.get('ellipse', None),
                    band_ts = rec.get('band_ts', None),
                    sedrec = rec.get('sedrec', None),
                    ts = rec.get('ts', None),
                    associations = rec.get('associations', None),
                    pivot_energy = rec.get('pivot_energy', None),
                    )
            else:
                src = self.ecat.lookup(name)
                src.model = rec['model']
            if neighbors: src.free[:]=False # not sure this is necessary
            src.index = index
            src.model.free = src.free # Don't think I still need this second copy of free
            return src
        
        def load_global_source(name, rec):
            if not self.quiet: print 'Loading global source %s for %d' % (name, index)
            gsrc = sources.GlobalSource(name=name, skydir=None, model=rec,
                dmodel = diffuse.diffuse_factory(self.config.diffuse[name], self.config.event_type_names))
            gsrc.index =index
            return gsrc

        p = pickle.load(self._z.open('pickle/HP12_%04d.pickle' % index))
        if not neighbors:
            global_sources = [load_global_source(name, rec) for name, rec \
                in zip(p['diffuse_names'], p['diffuse']) if name not in self.ecat.names]
            self.global_count = len(global_sources)
            for s in global_sources: self.append(s)
        local_sources = [load_local_source(name, rec) for name,rec in p['sources'].items()]
        if not neighbors: self.local_count = len(local_sources)
        for s in local_sources: self.append(s)
        
    def neighbors(self):
        from pointlike import IntVector
        b12 = skymaps.Band(12)
        v = IntVector()
        b12.findNeighbors(int(self.index),v) 
        return list(v)

    def __repr__(self):
        return '%s.%s : %d global, %d local, %d total sources for ROI %d' \
            % (self.__module__, self.__class__.__name__, self.global_count,self.local_count,  
            len(self), self.index)
            
