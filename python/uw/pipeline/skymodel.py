"""
Manage the sky model for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/skymodel.py,v 1.3 2011/01/20 16:08:04 burnett Exp $

"""
import os, pickle, glob, types, copy
import numpy as np
from skymaps import SkyDir, Band, IsotropicSpectrum, DiffuseFunction
from uw.utilities import keyword_options, makerec
from uw.like import  Models
from uw.like import pointspec_helpers

class Singleton(object):
    _instance={}
    def __init__(self,constructor):
        self.constructor = constructor
        self.key=str(constructor)
    def set_instance(self,  *pars):
        Singleton._instance[self.key]= self.constructor(*pars)
    def instance(self):
        try:
            return Singleton._instance[self.key]
        except KeyError:
            print 'SkyModel: Global source %s not initialized' % self.key
            raise
            
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
        print 'loaded new diffuse map, %s with lock=%s' %(dfile, lock)
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
        print 'loaded new isotropic spectrum, %s, with lock=%s' %(dfile, lock)
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
    def __init__(self, **kwargs):
        super(GlobalSource, self).__init__(**kwargs)
    
class ExtendedSource(Source):
    def __str__(self):
        return self.name + ' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')    
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)


class ExtendedCatalog( pointspec_helpers.ExtendedSourceCatalog):
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
                s.name = os.path.split(Isotropic._dfile)[-1]
            elif ext=='.fits' or ext=='.fit':
                s.dmodel = [Diffuse(dfile).instance()]
                s.name = os.path.split(Diffuse._dfile)[-1]
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
        
def get_association_class(adict):
    """Given association dictionary, decide what class to ascribe the source to.  Partly guesswork!
        original version by Eric Wallace
        added tev
    """
    if adict is None: return '   '
    cat_list=adict['cat']
    priority = '''bllac bzcat cgrabs crates crates_fom seyfert seyfert_rl qso agn 
                vcs galaxies pulsar_lat snr snr_ext pulsar_high pulsar_low pulsar_fom
                msp pwn hmxb lmxb globular tev ibis lbv dwarfs
               '''.split()
    others = ['ostar', 'starbursts', 'ocl']
    ass_class = ['bzb','bzcat']+['bzq']*3+['agn']*6+['LAT psr']+\
                ['snr']*2 + ['psr']*4 + ['pwn'] + ['hmxb'] + ['lmxb']+ ['glc'] +['tev'] + 3*['None']
    cls = None
    for c,a in zip(priority,ass_class):
        if c in cat_list:
            cls = a
            break
    if cls is None:
        if cat_list[0] not in others: print 'warning: ''%s'' not recognized' % cat_list
        return '   '
    if cls == 'bzcat': #special, use first 3 chars of source name
        cls = adict['name'][cat_list.index(c)][:3].lower()
    return cls
        
     
def create_catalog(outdir, **kwargs):
    """ make source and roi rec files in the skymodel directory, containing updated values for all the 
        sources found in the pickle folder
            note that order of first 6 parameters is assumed above
        also make a diffuse parameter dictionary
    """
    assert os.path.exists(os.path.join(outdir,'pickle')), 'pickle folder not found under %s' %outdir
    filelist = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
    assert len(filelist)>0, 'no .pickle files found in %s/pickle' % outdir 
    failed,maxfail = 0,kwargs.pop('maxfail',10)
    ignore_exception = kwargs.pop('ignore_exception',False)
    save_local = kwargs.pop('save_local',True) 
    assert len(kwargs.keys())==0, 'unrecognized kw %s' %kwargs 
    filelist.sort()
    
    class CatalogRecArray(object):
        def __init__(self, minflux=1e-16, update_position=False, ts_min=25):
            self.minflux=minflux
            self.update_position = update_position
            self.count=self.reject=self.moved=0
            self.rejected=[]
            self.ts_min=ts_min
            self.moved = 0
            self.colnames ="""name ra dec 
                pnorm pindex cutoff 
                pnorm_unc pindex_unc cutoff_unc
                e0 pivot_energy 
                flux flux_unc
                beta beta_unc
                modelname 
                ts band_ts bts10
                fit_ra fit_dec a b ang qual delta_ts
                id_prob aclass hp12
                """.split() 
            self.rec =makerec.RecArray(self.colnames) 

        def process(self,pk, cnan=np.nan):
            sources = pk['sources']
            for name,entry in sources.items():
                self.count+=1
                skydir = entry['skydir']
                data = [name, skydir.ra(), skydir.dec()]
                model  = entry['model']
                p,p_relunc = model.statistical()
                if p[0]< self.minflux or np.any(np.isnan(p[:2])):
                    self.reject+=1
                    self.rejected.append(name)
                    continue
                p_unc = p*p_relunc
                psr_fit =  model.name=='ExpCutoff'
                data += [p[0],     p[1],     p[2] if psr_fit else cnan, ]
                data += [p_unc[0], p_unc[1] ,p_unc[2] if psr_fit else cnan,]
                pivot_energy = entry.get('pivot_energy',model.e0)
                if pivot_energy=='None': pivot_energy=model.e0
                e0 = model.e0 if model.name!='LogParabola' else p[3]
                flux = model(e0)
                flux_unc = flux*p_relunc[0]
                data += [e0, pivot_energy]
                data += [flux, flux_unc]
                if psr_fit:
                    data += [cnan,cnan, 'ExpCutoff']
                else:
                    data += [cnan,cnan, 'PowerLaw'] if len(p)<3 or p[2]<=0.01 else [p[2], p_unc[2], 'LogParabola']
                ts = entry['ts']
                data += [ts, entry['band_ts']]
                data += [sum(entry['band_info']['ts'][7:])] ### note assumption that 10 GeV starts at 7
                ellipse = entry.get('ellipse', None)
                if ellipse is None or np.any(np.isnan(ellipse)):
                    data += [np.nan]*7
                else:
                    data += ellipse 
                    if self.update_position and ts>self.ts_min:
                        fit_ra, fit_dec, a, b, ang, qual, delta_ts = ellipse
                        #assert False, 'break'
                        if qual<5 and delta_ts<9 and a < 0.2 :
                            data[1:3] = [fit_ra, fit_dec]
                            self.moved +=1
                adict =  entry.get('associations', None)
                data.append( cnan if adict is None else adict['prob'][0])
                data.append('%-7s'%get_association_class(adict))
                data.append(Band(12).index(skydir))

                assert len(data)==len(self.colnames), 'mismatch between names, data'
                #assert not np.any(np.isinf(data[1:])), 'found inf for source %s' % name 
                self.rec.append(*data)

        def __call__(self): 
            print 'processed %d sources, rejected %d' %(self.count, self.reject)
            if self.reject>0:
                print 'rejected: flux<%.1e ' %self.minflux, self.rejected
            if self.update_position:
                print '\tmoved %d sources according to localization' % self.moved
            return self.rec()
        
    class DiffuseRecArray(object):
    
    
        def __init__(self, ndiff=3, nside=12):
            self.ndiff=ndiff
            if ndiff==3:
                self.colnames = """name galnorm galindex isonorm 
                                galnorm_unc galindex_unc isonorm_unc 
                                 loglike chisq""".split()
            else:
                self.colnames = """name galnorm galindex isonorm  isonorm2
                                galnorm_unc galindex_unc isonorm_unc isonorm2_unc
                                loglike chisq""".split()
            self.rec = makerec.RecArray(self.colnames)
            self.nside=nside
            
        def process(self, pk):
            name = pk['name']
            p, p_relunc = [np.hstack([m.statistical()[i] for m in pk['diffuse']] ) for i in range(2)]
            if len(p)!=self.ndiff:
                msg = 'unexpected number of diffuse parameters, %d vs %d, processing %s' % (len(p),self.ndiff,name)
                #print msg
                p = p[:self.ndiff]; p_relunc = p_relunc[:self.ndiff]
            data = [name] + list(np.hstack((p, p*p_relunc)))
            data += [pk['logl']]
            counts = pk['counts']
            obs,mod = counts['observed'], counts['total']
            data += [((obs-mod)**2/mod).sum()]
            assert len(data)==len(self.colnames), 'mismatch between names, data'
            self.rec.append(*data)
            
        def __call__(self):
            t = self.rec()
            n = 12*self.nside**2
            if len(t)!=n: 
                q = np.array([ 'HP12_%04d' % i not in t.name for i in range(n)])
                msg  = 'pixel file missing entries %s' % np.arange(n)[q]
                print msg
                raise Exception, msg
            return t

        
    crec = CatalogRecArray(**kwargs)
    drec = DiffuseRecArray()
    for fname in filelist:
        try:
            p = pickle.load(open(fname))
            crec.process(p)
            drec.process(p)
        except Exception, arg:
            print 'Failed to load file  %s: %s' % (fname, arg)
            if not ignore_exception or failed>maxfail: raise
            failed +=1
    print 'read %d entries from %s (%d failed)' % (len(filelist),outdir,failed)
    for r,name in ((crec, 'sources'), (drec, 'rois')):
        fname = '%s_%s.rec'%(name, outdir) if not save_local else '%s/%s.rec' % (outdir, name)
        rec = r()
        pickle.dump(rec, open(fname,'wb'))
        print 'saved %d entries to pickle file %s' % (len(rec), fname)
        
 
     