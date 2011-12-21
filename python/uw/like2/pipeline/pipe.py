"""
Main entry for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipe.py,v 1.6 2011/12/09 16:07:44 burnett Exp $
"""
import os, types, glob, time, copy
import cPickle as pickle
import numpy as np
from . import processor, parallel
from .. import main, roisetup, skymodel
from uw.utilities import makerec
from .parallel import  AssignTasks, get_mec, kill_mec, free

class Pipe(roisetup.ROIfactory):
    """ This is a subclass of ROIfactory,
    It implements the capabilities needed for 
    IPython paralellization: the constructor defines the environment, and if called with
    an index, it will run a full analysis for that ROI.
    
    """
    def __init__(self, indir, dataset,  **kwargs):
        """
        parameters
        ----------
        indir : string
            name of a folder containing the sky model description, passed to 
               skymodel.SkyModel
        dataset : instance of a DataSpec object or a string used to lookup
              
        keyword arguments:
            processor : if set, it is a function
        """
        self.nside      = kwargs.pop('nside', 12)
        self.process_kw = kwargs.pop('process_kw', dict())
        self.fit_kw     = kwargs.pop('fit_kw', dict())
        self.skymodel_kw= kwargs.pop('skymodel_kw', dict())
        self.analysis_kw =kwargs.pop('analysis_kw', dict())
        self.roi_kw      =kwargs.pop('roi_kw', dict())
 
        # TODO: get association back in
        #associator = kwargs.pop('associate', 'all_but_gammas')
        #if associator is not None and associator[0]=='$':
        #    print 'associations will come from the catalog %s' %associator
        #    associator = associate.Association(associator)
        #elif associator is not None and associator!='None':
        #    associator = associate.SrcId('$FERMI/catalog', associator)
            #print 'will associate with catalogs %s' % associator.classes
        #self.process_kw['associate'] = associator
        self.selector = skymodel.HEALPixSourceSelector
        self.selector.nside = self.nside

        # the function that will be called for each ROI: default is processor.process
        self.processor = kwargs.pop('processor', processor.process)
        if type(self.processor)==types.StringType:
            self.processor = eval(self.processor)

        super(Pipe, self).__init__(indir, dataset, 
            analysis_kw=self.analysis_kw, skymodel_kw=self.skymodel_kw)
       
    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        show = """analysis_kw selector skymodel dataset fit_kw convolve_kw process_kw""".split()
        for key in show:
            s += '\t%-20s: %s\n' %(key,
                self.__dict__.get(key,'not in self.__dict__!'))
        return s
        
    def __call__(self, index, **kwargs):
        """ perform analysis for the ROI
        Note that it is done by a function in the processor module
        requires a setup call
        """
        roi_kw = self.roi_kw.copy()
        roi_kw.update(**kwargs)
        roi = self.roi(index , roi_kw=roi_kw )
        
        return self.processor(roi, **self.process_kw)

    def roi(self, *pars, **kwargs):
        try:
            t = super(Pipe, self).roi(*pars, **kwargs)
        except Exception, e:
            print 'trying again after exception %s' %e
            t = super(Pipe, self).roi(*pars, **kwargs)

        return main.ROI_user(t)
        
    def names(self):
        return [self.selector(i).name() for i in range(12*self.nside**2)]
        
    def setup_mec(self, block=True):
        mec = parallel.get_mec()
        assert len(mec)>0, 'no engines to setup'
        print 'setting up %d engines...' % len(mec)
        ar = mec[:].execute(self(), block=block)
        return ar

class Setup(dict):
    """ Setup is a dictionary with variable run parameters
    also will setup and run parallel
    """

    def __init__(self, version=None,  **kwargs):
        """ generate setup string"""
        if version is None: version = int(open('version.txt').read()) if os.path.exists('version.txt') else -1
        indir=kwargs.pop('indir', 'uw%02d'%(version)) 
        self.outdir=outdir=kwargs.pop('outdir', 'uw%02d'%(version+1))
        if not os.path.exists(outdir): os.mkdir(outdir)
        self.version=version
        if os.name=='nt':
            os.system('title %s %s'% (os.getcwd(), indir))
        self.update(dict( cwd='$HOME/analysis/'+os.path.split(os.getcwd())[-1], 
                indir = indir,
                auxcat='',
                outdir=outdir,
                pass_number=0,
                datadict = dict(dataname='P7_V4_SOURCE'),
                diffuse = ('ring_24month_P76_v1.fits', 'isotrop_21month_P76_v2.txt'),
                extended= None,
                alias= {}, 
                skymodel_extra = '',
                irf = 'P7SOURCE_V6',
                associator ='all_but_gammas',
                sedfig_dir = None,
                tsmap_dir  = None, tsfits=False,
                localize=True,
                emin=100, emax=316227,
                processor='processor.process',
                fix_beta=False, dofit=True,
                source_kw=dict(),
                fit_kw=dict(ignore_exception=False, use_gradient=True, call_limit=1000),
                repivot = False,
                update_positions=None,
                free_index=None,
                tables = None,  #roi_maps.ROItables("%(outdir)s", skyfuns=(
                                # (roi_tsmap.TSCalc, 'ts', dict(photon_index=2.0),) 
                                #  (ts_map.KdeMap, "kde", dict()),))
                dampen = 1.0,
                setup_cmds='',
                ))
        self.update(kwargs)
        # first-order replace
        for key in self:
            if type(self[key])== types.StringType and self[key].find('%(')>=0: 
                self[key]=self[key]%self
                self['tables'] = self['tables']%self
                #print 'fix key %s: %s' % (key, self[key])
        self.setup_string =  """\
import os, pickle; os.chdir(os.path.expandvars(r"%(cwd)s"));%(setup_cmds)s
from uw.like2.pipeline import pipe; from uw.like2 import skymodel
g=pipe.Pipe("%(indir)s", %(datadict)s, 
        skymodel_kw=dict(auxcat="%(auxcat)s",diffuse=%(diffuse)s,
            extended_catalog_name="%(extended)s", update_positions=%(update_positions)s,
            free_index=%(free_index)s,
            alias =%(alias)s, %(skymodel_extra)s), 
        analysis_kw=dict(irf="%(irf)s", minROI=%(minROI)s, maxROI=%(maxROI)s, emin=%(emin)s,emax=%(emax)s),
        cache="",
        processor="%(processor)s",
        process_kw=dict(outdir="%(outdir)s", pass_number=%(pass_number)s,
            fit_kw = %(fit_kw)s,
            tsmap_dir=%(tsmap_dir)s,  tsfits=%(tsfits)s,
            sedfig_dir=%(sedfig_dir)s,
            localize=%(localize)s,
            associate="%(associator)s",
            fix_beta= %(fix_beta)s, dofit=%(dofit)s,
            tables= %(tables)s,
            repivot=%(repivot)s, dampen=%(dampen)s,
            ),
    ) 
n,chisq = len(g.names()), -1
""" % self
        self.mecsetup=False
        print 'Pipeline input, output: %s -> %s ' % (indir, outdir)
        
    def __call__(self): return self.setup_string
    
    def g(self):
        exec(self(), globals()) # should set g
        return g

    def setup_mec(self, block=True):
        """ send our setup string to the current engines
        """
        mec = parallel.get_mec()
        assert len(mec)>0, 'no engines to setup'
        print 'setting up %d engines...' % len(mec)
        mec[:].clear()
        ar = mec[:].execute(self(), block=block)
        self.mecsetup=True
        self.dump()
        return ar
        
    def dump(self):
        """ save the setup string, and dictionary """
        open(self.outdir+'/setup_string.txt', 'w').write(self.setup_string)
        open(self.outdir+'/config.txt', 'w').write(self.__str__())
    
    def run(self, fn=None, tasklist=None, **kwargs):
        """ note: fn should be  lambda x:(x,g(x))
        """
        if not self.mecsetup:
            self.setup_mec()
        outdir = self.outdir
        if outdir is not None:
            if not os.path.exists(outdir): 
                os.mkdir(outdir)
        if fn is None:
            fn = lambda x:(x,g(x))
        print 'writing results to %s' %outdir
        local = kwargs.pop('local', False)
        ignore_exception = kwargs.pop('ignore_exception', False)
        sleep_interval = kwargs.pop('sleep_interval', 60)
        if tasklist is None:
            mytasks =[ (864+i/2) if i%2==0 else (863-i/2) for i in range(1728)]
        else: 
            mytasks = [t for t in tasklist]
        lc= AssignTasks(fn, mytasks, 
            timelimit=5000, local=local, ignore_exception=ignore_exception,  )
        if sleep_interval>0:
            lc(sleep_interval)
        return lc


def getg(setup_string):
    """ for testing"""
    exec(setup_string, globals()) # should set g
    return g

def roirec(version=None, nside=12):
    if version is None:
        version = int(open('version.txt').read())
    roi_files = glob.glob('uw%02d/pickle/*.pickle'%version)
    roi_files.sort()
    n = 12*nside**2
    if len(roi_files)<n:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(n) if x not in t]
        print 'misssing roi files: %s' % missing
        raise Exception('misssing roi files: %s' % missing)
        
    recarray = makerec.RecArray('name chisq loglike'.split())
    for fname in roi_files:
        p = pickle.load(open(fname))
        counts = p['counts']
        obs,mod = counts['observed'], counts['total']
        chisq =  ((obs-mod)**2/mod).sum()

        recarray.append(p['name'], chisq, p['logl'])
    return recarray()

def check_missing_files(folder):
    roi_files = sorted(glob.glob(os.path.join(folder, '*.pickle')))
    n = 1728
    missing = []
    if len(roi_files)<n:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(n) if x not in t]
        if len(missing)<10:
            print '\tmisssing roi files: %s' % missing
        else:
            print '\tMissing %d roi files: %s...' %( len(missing), missing[:10])
    return missing    



def roirec(outdir, check=False):
    roi_files = sorted(glob.glob(os.path.join(outdir, 'pickle/*.pickle')))
    n = 1728
    if len(roi_files)<n:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(n) if x not in t]
        if len(missing)<10:
            print '\tmisssing roi files: %s' % missing
        else:
            print '\tMissing %d roi files: %s...' %( len(missing), missing[:10])
        if check: return missing    
    if check: return[]#    return None # raise Exception('misssing roi files: %s' % missing)
        
    recarray = makerec.RecArray('name chisq loglike prevlike niter'.split())
    bad = []
    for fname in roi_files:
        p = pickle.load(open(fname))
        if 'counts' not in p.keys():
            bad.append(fname)
            continue
        counts = p['counts']
        obs,mod = counts['observed'], counts['total']
        chisq =  ((obs-mod)**2/mod).sum()
        logl_list = p.get('prev_logl', [0])
        recarray.append(p['name'], chisq, p['logl'], logl_list[-1], len(logl_list))
    if len(bad)>0:
        print 'no fit info in file(s) %s' % bad if len(bad)<10 else (str(bad[:10])+'...')
    return recarray()

def check_converge(month, tol=10, add_neighbors=True, log=None):
    """ check for convergence, ROI that have been updated
    month: int or string
        if int, intrepret as a month, else a folder
        
    """
    from pointlike import IntVector
    from skymaps import Band
    outdir = 'month%02d'%month if type(month)==types.IntType else month
    print '%s:' %outdir
    r = roirec(outdir)
    if r is None: return
    diff = r.loglike-r.prevlike
    dmin,dmax = diff.min(), diff.max()
    rmin,rmax = list(diff).index(dmin), list(diff).index(dmax)
    changed = set(np.arange(1728)[np.abs(diff)>tol])
    print >>log, '\tpass %d:  %d changed > %d, min, max: %d(#%d) %d(#%d)' % (max(r.niter),len(changed), tol, dmin,rmin,dmax,rmax),
    if not add_neighbors: return list(changed)
    nbrs = set()
    b12 = Band(12)
    for x in changed: 
        v = IntVector()
        b12.findNeighbors(int(x),v) # int is tricky
        for n in v:
            nbrs.add(n)
    q =list(changed.union( nbrs))
    print >>log, ' (total %d)' %len(q)
    if log is not None: log.flush()
    return q
    
       
def pmain( setup, fn, taskids=None, local=False,
        ignore_exception=False, 
        logpath='log',
        sleep_interval=60,
        ):
        
    """
    Parameters
        setup : Setup object, implements
            setup_string = setup()
                Python code, must define an object "g". It must implement:
                a function g(n), n an integer
                g.names() must return a list of task names
            outdir =setup.outidr : directory to save files
        taskids : None or a list of integers
            the integers should be in range(0,len(g.names())
    """
    setup_string = setup()
    outdir = setup.outdir
    if outdir is not None:
        if not os.path.exists(outdir): 
            os.mkdir(outdir)
        print 'writing results to %s' %outdir
            
    lc= AssignTasks(fn, taskids, 
        timelimit=5000, local=local, ignore_exception=ignore_exception,
         )
    if sleep_interval>0:
        lc(sleep_interval)
       
    return lc
   

