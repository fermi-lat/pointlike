"""
Manage the sky model for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/skymodel.py,v 1.2 2011/01/19 01:28:09 burnett Exp $

"""
import os, pickle, glob, types, copy
import numpy as np
from skymaps import SkyDir, Band, IsotropicSpectrum, DiffuseFunction
from uw.utilities import keyword_options, makerec
from uw.like import pointspec_helpers, Models

class Singleton(object):
    _instance={}
    def __init__(self,constructor):
        self.constructor = constructor
        self.key=str(constructor)
    def set_instance(self,  *pars):
        Singleton._instance[self.key]= self.constructor(*pars)
    def instance(self):
        return Singleton._instance[self.key]
        
class Diffuse(Singleton):
    """ manage a skymaps.DiffuseFunction object, create only when new fits file found to open
    """
    _dfile = None
    locked = False
    key = None
    def __init__(self, dfile, lock=False):
        super(Diffuse,self).__init__(DiffuseFunction)
        if Diffuse.locked or dfile==Diffuse._dfile: return
        Diffuse._dfile=dfile
        self.set_instance(dfile)
        print 'loaded new diffuse map, %s' %dfile
        Diffuse.locked = lock
        
class Isotropic(Singleton):
    """ manage a skymaps.IsotropicSpectrum object, create only when new fits file found to open
    """
    _dfile = None
    locked = False
    key = None
    def __init__(self, dfile, lock=False):
        super(Isotropic,self).__init__(IsotropicSpectrum)
        if Isotropic.locked or dfile==Isotropic._dfile: return
        Isotropic._dfile=dfile
        self.set_instance(dfile)
        print 'loaded new isotropic spectrum, %s' %dfile
        Isotropic.locked = lock

class Source(object):
    """ base class for various sources
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if 'model' not in kwargs:
            self.model=Models.PowerLaw(p=[1e-13, 2.2])
        self.free = self.model.free.copy()  # save copy of initial free array to restore
    def freeze(self, freeze):
        self.model.free[:] = False if freeze else self.free
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
                source.free = source.model.free.copy()
                if source.model.name=='LogParabola': source.free[-1]=False # E_break not free
                return source
                # since Josh checks for his ExtendedSource class, I have to modify it rather than 
                # use my own
                #spatial_model = source.spatial_model
                #return ExtendedSource(name=name, skydir=source.skydir,
                #    model = source.model, 
                #    spatial_model = source.spatial_model)
                #    
        raise Exception( ' extended source %s not found' % name)
    
class SkyModel(object):
    """
    """
    
    defaults= (
        ('extended_catalog_name', 'Extended_archive_jbb03',  'name of folder with extended info'),
        ('diffuse', ('ring_24month_P74_v1.fits', 'isotrop_21month_v2.txt'),   'pair of diffuse file names: use to locked'),
        ('auxcat', None, 'name of auxilliary catalog of point sources to append',),
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
        if self.diffuse is not None:
            assert len(self.diffuse)==2, 'expect 2 diffuse names'
            x = map(lambda f:os.path.expandvars(os.path.join('$FERMI','diffuse',f)), self.diffuse)
            Diffuse(x[0], True)
            Isotropic(x[1], True)
            
        self.load_auxcat()
        
    def __str__(self):
        return 'SkyModel %s' %self.folder\
                +'\n\t\tdiffuse: %s' %list(self.diffuse)\
                +'\n\t\textended: %s' %self.extended_catalog_name
                
    def load_auxcat(self):
        if self.auxcat is None or self.auxcat=='': return
        cat = self.auxcat 
        if not os.path.exists(cat):
            cat = os.expandvars(os.path.join('$FERMI','catalog', cat))
        if not os.path.exists(cat):
            raise Exception('auxilliary catalog %s not found locally or in $FERMI/catalog'%self.auxcat)
        ss = makerec.load(cat)
        names = set([s.name for s in self.sources])
        tagged=set()
        print 'process auxcat %s' %cat
        for s in ss:
            sname = s.name.replace('_',' ')
            if sname  not in names: 
                skydir=SkyDir(float(s.ra), float(s.dec))
                index=self.index(skydir)
                self.sources.append(PointSource(name=s.name, skydir=skydir, index=index))
                print '\tadded new source %s at ROI %d' % (s.name, index)
            else: 
                print '\t source %s is  in the model will remove if ra<0' % sname
                if s.ra<0: 
                    tagged.add(sname)
        if len(tagged)==0: return
        for s in self.sources:
            if s.name in tagged: self.sources.remove(s)
            
    def _setup_extended(self):
        if self.extended_catalog_name is None: return None
        extended_catalog_name = \
            os.path.expandvars(os.path.join('$FERMI','catalog',self.extended_catalog_name))
        if not os.path.exists(extended_catalog_name):
            raise Exception('extended source folder "%s" not found' % extended_catalog_name)
        self.extended_catalog= ExtendedCatalog(extended_catalog_name)
        #print 'Loaded extended catalog %s' % self.extended_catalog_name
        
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
                if len(t)<2: # always assume first two are global
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
            ps.freeze(freeze) 
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
                s.dmodel = [Isotropic(dfile).instance()]
            elif ext=='.fits' or ext=='.fit':
                s.dmodel = [Diffuse(dfile).instance()]
            else:
                raise Exception('unrecognized diffuse file extention %s' % dfile)
            s.smodel=s.model
        center = self.skydir(index)
#        extended =  [es for es in self.extended_sources if es.near(center,radius)]
        extended =  [es for es in self.extended_sources if center.difference(es.skydir)<np.radians(radius)]
        for es in extended:
            es.model.free[:] = es.free if self.index(es.skydir)==index else False
            es.dmodel = es.spatial_model # references needed 
            es.smodel = es.model
        return globals, extended
        
     