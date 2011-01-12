"""
Manage the sky model for the UW all-sky pipeline
$Header$

"""
import os, pickle, glob, types, copy
import numpy as np
from skymaps import SkyDir, Band, IsotropicSpectrum, DiffuseFunction
from uw.utilities import keyword_options
from uw.like import pointspec_helpers

class Source(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.free = self.model.free.copy()  # save copy of free to restore
    def __str__(self):
        return self.name + ' '+ self.skydir.__str__() + self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')
 
class PointSource(Source):
   def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)
        
class GlobalSource(Source):
    pass
    def __init__(self, **kwargs):
        super(GlobalSource, self).__init__(**kwargs)
    
class ExtendedSource(Source):
    def __str__(self):
        return self.name + ' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')    
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)


class ExtendedCatalog(pointspec_helpers.ExtendedSourceCatalog):
    """ subclass to add this lookup function """

    def __init__(self, *pars, **kwargs):
        """ initialize by also filling an array with all source spectral models"""
        super(ExtendedCatalog,self).__init__(*pars, **kwargs)
        self.sources = self.get_sources(SkyDir(),180)
        assert len(self.sources)==len(self.names), 'inconsistent list lengths'
        
    def lookup(self, name):
        """ return an ExtendedSource object, or None if not found """
        for source in self.sources:
            if name==source.name:
                spatial_model = source.spatial_model
                return ExtendedSource(name=name, skydir=source.skydir,
                    model = source.model, 
                    spatial_model = source.spatial_model)
                    
        raise Exception( ' extended source %s not found' % name)
    
class SkyModel(object):
    """
    """
    
    defaults= (
        ('extended_catalog_name', '$FERMI/catalog/Extended_archive_jbb03',  'name of folder with extended info'),
        ('nside',      12,   'HEALpix nside'),
        ('quiet',   False,  'make quiet' ),
    )
    
    @keyword_options.decorate(defaults)
    def __init__(self, folder=None,  **kwargs):
        """
        folder : string or None
            name of folder to find all files defining the sky model
        """
        keyword_options.process(self, kwargs)

        if folder is None:
            folder = 'uw%02d' % int(open('version.txt').read())
        self.folder = os.path.expandvars(folder)
        if not os.path.exists(self.folder):
            raise Exception('sky model folder %s not found' % folder)
        self._setup_extended()
        self._load_sources()
        #if not self.quiet:
        #    print 'loaded %d sources from %s' % (len(self.sources), self.folder)
        
    def __str__(self):
        return 'SkyModel from  %s\n\7\textended catalog:%s' %(self.folder, self.extended_catalog)
    def _setup_extended(self):
        if self.extended_catalog_name is None: return None
        self.extended_catalog_name = os.path.expandvars(self.extended_catalog_name)
        if not os.path.exists(self.extended_catalog_name):
            raise Exception('extended source folder "%s" not found' % self.extended_catalog_name)
        self.extended_catalog= ExtendedCatalog(self.extended_catalog_name)
        print 'Loaded extended catalog %s' % self.extended_catalog_name
        
    def get_extended_sources(self,skydir, radius):
        """ add any extended sources with center within the outer radius.
            set parameters free if the center is inside the HEALpix region
        """
        if self.extended_catalog is None: return []
        ret =self.extended_catalog.get_sources(skydir, radius)
        return ret    

    def _load_sources(self):
        """
        run through the pickled roi dictionaries, create lists of point and extended sources
        """
        self.sources= []
        files = glob.glob(os.path.join(self.folder, 'pickle', '*.pickle'))
        files.sort()
        self.global_sources = []  # allocate list to index parameters for global sources
        self.extended_sources=[]  # list of unique extended sources
        self.changed=set() # to keep track of extended models that are different from catalog
        for i,file in enumerate(files):
            p = pickle.load(open(file))
            index = int(os.path.splitext(file)[0][-4:])
            assert i==index, 'logic error: file name %s inconsistent with expected index %d' % (file, i)
            sources = p['sources']
            for key,item in sources.items():
                self.sources.append( PointSource(name=key, 
                    skydir=item['skydir'], model= item['model'],
                    ts=item['ts'],band_ts=item['band_ts'], index=index)
                    )
            # make a list of extended sources used in the model   
            t = []
            for name, model in zip(p['diffuse_names'] , p['diffuse']):
                if len(t)<2: # alwasy assume first two are global
                    t.append(GlobalSource(name=name, model=model, skydir=None, index=index))
                else:
                    es = self.extended_catalog.lookup(name)
                    if es is None:
                        raise Exception( 'Extended source %s not found in extended catalog' %name)
                    if self.index(es.skydir)!=index: continue
                    
                    if es.model.name!=model.name:
                        if name not in self.changed:
                            print 'catalog model %s changed from %s for %s'% (es.model.name, model.name, name)
                        self.changed.add(name)
                    else:
                        es.model=model #update with fit values
                    self.extended_sources.append(es)
            self.global_sources.append(t)
 
    def skydir(self, index):
        return Band(self.nside).dir(index)
    def index(self, skydir):
        return Band(self.nside).index(skydir)
    
    def select(self, selection, freeze=False):
        """ return a list of point sources objects
        selection : function
            for example, lamda s : s.index=888
        freeze : bool
            determine 
        
        """
        ret =  [source for source in self.sources if selection(source)]
        for ps in ret:
            ps.model.free[:] = False if freeze else ps.free
        return ret
        
    def get_point_sources(self, index, radius=10):
        """
        return a list of PointSource objects appropriate for the ROI
        """
        ps_roi = self.select(lambda s: s.index==index, freeze=False)
        center = self.skydir(index)
        ps_rest = self.select(lambda s: s.near(center,radius) and not s.index==index, freeze=True)
        return ps_roi+ps_rest
        
    def get_diffuse_sources(self, index,radius=10):
        """return diffuse, global and extended sources needed for this HEALpix index
            always the global diffuse, and perhaps local extended sources.
            For the latter, make parameters free only if center is in the pixel
        """
        globals = self.global_sources[index]
        for s in globals:
            dfile = os.path.expandvars(os.path.join('$FERMI','diffuse', s.name))
            assert os.path.exists(dfile), 'file %s not found' % dfile
            ext = os.path.splitext(dfile)[-1]
            if ext=='.txt':
                s.dmodel = [IsotropicSpectrum(dfile)]
            elif ext=='.fits' or ext=='.fit':
                s.dmodel = [DiffuseFunction(dfile)]
            else:
                raise Exception('unrecognized diffuse file extention %s' % dfile)
            s.smodel=s.model
        center = self.skydir(index)
        extended =  [es for es in self.extended_sources if es.near(center,radius)]
        for es in extended:
            es.model.free[:] = es.free if self.index(es.skydir)==index else False
            es.dmodel = es.spatial_model # references needed 
            es.smodel = es.model
        return globals, extended
        
     