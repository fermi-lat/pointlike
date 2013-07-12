"""
Main entry for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipe.py,v 1.37 2013/05/04 14:44:22 burnett Exp $
"""
import os, types, glob, time, copy
import cPickle as pickle
import numpy as np
from . import processor, associate
from .. import main, roisetup, skymodel
from uw.like import Models 
from uw.utilities import makerec

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
            if None, get from the model configuration
              
        keyword arguments:
            processor : if set, it is a function
        """
        self.nside      = kwargs.pop('nside', 12)
        self.process_kw = kwargs.pop('process_kw', dict())
        self.fit_kw     = kwargs.pop('fit_kw', dict())
        self.skymodel_kw= kwargs.pop('skymodel_kw', dict())
        self.analysis_kw =kwargs.pop('analysis_kw', dict())
        self.roi_kw      =kwargs.pop('roi_kw', dict())
 
        self.selector = skymodel.HEALPixSourceSelector
        self.selector.nside = self.nside

        # the function that will be called for each ROI: default is processor.process
        # note that if the name has arguments, extract them as process_kw
        self.processor = kwargs.pop('processor', processor.process)
        if type(self.processor)==types.StringType:
            self.processor = eval(self.processor)

        super(Pipe, self).__init__(indir, dataset, 
            analysis_kw=self.analysis_kw, skymodel_kw=self.skymodel_kw, **kwargs)
       
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
        requires a setup call to define the function on a remote machine
        returns a tuple: (index, processor call return)
        """
        roi_kw = self.roi_kw.copy()
        roi_kw.update(**kwargs)
        roi = self.roi(index , roi_kw=roi_kw )
        
        return index, self.processor(roi, **self.process_kw)

    def roi(self, *pars, **kwargs):
        try:
            t = super(Pipe, self).roi(*pars, **kwargs)
        except Exception, e:
            #if not e.startswith('Root'): 
            #    raise
            print 'trying again after exception "%s"' %e
            t = super(Pipe, self).roi(*pars, **kwargs)

        return main.ROI_user(t)
        
    def names(self):
        return [self.selector(i).name() for i in range(12*self.nside**2)]
        

class Setup(dict):
    """ Setup is a dictionary with variable run parameters.
    It manages creation of Pipe object, either locally, or in parallel
    also will setup and run parallel
    """

    def __init__(self, indir,  **kwargs):
        """ generate setup string"""
        self.outdir=outdir=kwargs.pop('outdir', indir)
        self.indir=indir
        self.quiet = kwargs.pop('quiet', False)
        if not os.path.exists(outdir): os.mkdir(outdir)
        if os.name=='nt':
            os.system('title %s %s'% (os.getcwd(), indir))
        self.update(dict( cwd=os.getcwd(), #'$HOME/analysis/'+os.path.split(os.getcwd())[-1], 
                indir = indir,
                auxcat='',
                outdir=outdir,
                datadict = None,
                diffuse = None,
                emin=100, emax=316227, minROI=5, maxROI=5,
                extended= None, #flag to get from model
                skymodel_extra = '',
                process_extra='',
                irf = None, # flag to get from model 'P7SOURCE_V6',
                associator = None,
                sedfig_dir = None,
                tsmap_dir  = None, tsfits=False,
                localize=False,
                processor='processor.process',
                fix_beta=False, 
                source_kw=dict(),
                fit_kw=dict(ignore_exception=False, use_gradient=True, call_limit=1000),
                repivot = False,
                tables = None,  #roi_maps.ROItables("%(outdir)s", skyfuns=(
                                # (roi_tsmap.TSCalc, 'ts', dict(photon_index=2.0),) 
                                #  (ts_map.KdeMap, "kde", dict()),))
                dampen = 1.0,
                setup_cmds='',
                ))
        self.update(kwargs)
        # first-order replace
        #for key in self:
        #    if type(self[key])== types.StringType and self[key].find('%(')>=0: 
        #        self[key]=self[key]%self
        #        self['tables'] = self['tables']%self
        #        print 'fix key %s: %s' % (key, self[key])
        # Turn apparent args in processor kw to process_exta list.
        if '(' in self['processor']:
            p,q = self['processor'].split('(')
            self['processor'] = p
            self['process_extra']=q[:-1]

        self.create_string()
        self.mecsetup=False
        self.mypipe = None
        if not self.quiet:
            print 'Pipeline input, output: %s -> %s ' % (indir, outdir)
        #if os.path.exists(indir+'/config.txt'):
        #    input_config = eval(open(indir+'/config.txt').read())
        #    for key in 'extended diffuse irf'.split():
        #        if self[key] is None: 
        #            self[key]=input_config[key]
        #            print 'updating %s from skymodel' %key

        
    def create_string(self):
        self.setup_string= """\
import os, pickle; os.chdir(os.path.expandvars(r"%(cwd)s"));%(setup_cmds)s
from uw.like2.pipeline import pipe,associate; from uw.like2 import skymodel;
g=pipe.Pipe("%(indir)s", %(datadict)s, 
        skymodel_kw=dict(auxcat="%(auxcat)s", 
            extended_catalog_name=%(extended)s,
            %(skymodel_extra)s), 
        analysis_kw=dict(irf="%(irf)s", minROI=%(minROI)s, maxROI=%(maxROI)s, emin=%(emin)s,emax=%(emax)s),
        irf="%(irf)s",
        diffuse=%(diffuse)s,
        processor="%(processor)s",
        process_kw=dict(outdir="%(outdir)s", dampen=%(dampen)s,
            fit_kw = %(fit_kw)s,
            sedfig_dir=%(sedfig_dir)s,
            tsmap_dir=%(tsmap_dir)s,  tsfits=%(tsfits)s,
            localize=%(localize)s, associate=%(associator)s,
            fix_beta= %(fix_beta)s, repivot=%(repivot)s,
            tables= %(tables)s,
            %(process_extra)s
            ),
    ) 
n,chisq = len(g.names()), -1
""" % self
        
    def _repr_html_(self):
        """ used by the notebook for display """
        return '<h3>'+self.__class__.__name__+'</h3>'+'<ul><li>Setup: %s -> %s' % (self.indir, self.outdir)

    def __call__(self): return self.setup_string
    
    def g(self):
        """
        return an exec of self, the Pipe object
        also save reference as mypipe
        """
        exec(self(), globals()) # should set g
        self.mypipe = g
        return g
        
    def process(self, hp12):
        """ process locally
        """
        if self.mypipe is None: self.g()
        roi =  self.mypipe.roi(hp12)
        processor.process(roi, **self.mypipe.process_kw)

    def setup_mec(self, profile=None, sleep_interval=60, **kwargs):
        """ send our setup string to the engines associated with profile
        """
        self.mec = engines.Engines(profile=profile, sleep_interval=sleep_interval, **kwargs)
        self.mec.clear()
        self.mec.execute('import gc; gc.collect();'+self())
        self.mecsetup=True
        self.dump()
        
    def dump(self, override=False):
        """ save the setup string, and dictionary 
            TODO: make sure not to overwrite important info
        """
        open(self.outdir+'/setup_string.txt', 'w').write(self.setup_string)
        config_file=self.outdir+'/config.txt'
        if not os.path.exists(config_file) or override:
            open(config_file, 'w').write(self.__str__())
    
    def run(self, tasklist=None,  **kwargs):
        """ 
        kwargs
        ------

        """
        profile = kwargs.pop('profile','default')
        sleep_interval = kwargs.pop('sleep_interval', 60)

        if not self.mecsetup:
            self.setup_mec(profile=profile, sleep_interval=sleep_interval)
        outdir = self.outdir
        if outdir is not None:
            if not os.path.exists(outdir): 
                os.mkdir(outdir)

        print 'writing results to %s' %outdir
        if tasklist is None:
            mytasks =[ (864+i/2) if i%2==0 else (863-i/2) for i in range(1728)]
        else: 
            mytasks = [int(t) for t in tasklist]
        self.mec.submit( 'g', mytasks, **kwargs)




def getg(setup_string):
    """ for testing"""
    exec(setup_string, globals()) # should set g
    return g

def roirec(version, nside=12):
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
    #print '%s:' %outdir
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
    
       

class NotebookPipe(object):

    def __init__(self, analysisdir, indir, **kwargs):
        self.analysisdir = os.path.expandvars(analysisdir)
        self.indir=indir
        outdir = kwargs.pop('outdir', None)
        self.outdir = indir if outdir is None else outdir
        if not os.path.exists(self.analysisdir):
            os.makedirs(self.analysisdir)
        os.chdir(self.analysisdir)
        self.setup = Setup(self.indir, outdir=self.outdir, **kwargs)
       
    def __str__(self):
        return 'Pipeline processing in %s: %s->%s' % (self.analysisdir,self.indir, self.outdir)

    def _repr_html_(self):
        return '<h3>'+time.asctime()+'</h3>'+'<ul><li>Pipeline %s %s -> %s' % (self.__class__.__name__,self.indir, self.outdir)
   
    def runcheck(self, older=3600):
        """ check for jobs that did not finish within range
        """
        t = sorted(glob.glob(os.path.join(self.analysisdir, self.outdir, 'log', '*.txt')))
        assert len(t)==1728, 'did not find all 1728 files'
        ftimes = map(os.path.getmtime, t)
        return [n for n,f in enumerate(t) if os.path.getmtime(f)<time.time() -older]

  
class Update(NotebookPipe):
    """ Update a model 
        Note that dampen=0.5; may need adjustment
    """
    def __init__(self, analysisdir='.', indir='.', outdir=None, **kwargs):
        """
        
        """
        self.analysisdir = os.path.expandvars(analysisdir)
        self.indir=indir
        self.outdir = indir if outdir is None else outdir
        if not os.path.exists(self.analysisdir):
            os.makedirs(self.analysisdir)
        os.chdir(self.analysisdir)
        kw = self.defaults()
        kw.update(skymodel_extra = self.extra_check())
        kw.update(kwargs)
        self.setup = Setup(self.indir, outdir=self.outdir, **kw)
     
    def defaults(self):
        return dict( dampen=0.5, quiet=True)
    
    def extra_check(self):
        # if config has skymodel_kw, will add to the skymodel keywords by interpreting it and setting skymodel_extra
        # especially, for a filter entry
        config = eval(open(os.path.join(self.analysisdir, self.indir, 'config.txt')).read())
        skykw = config.pop('skymodel_kw', None)
        return  ','.join(['%s="%s"' % item for item in skykw.items() ]) if skykw is not None else ''

     
    def __call__(self): return self.setup
    def check(self):
        self.t = check_converge(self.outdir)
        
    def g(self): return self.setup.g()
    def run(self, **kwargs): 
        rc=self.setup.run(**kwargs)
        self.check()
        return rc
    def iterate(self, full=False):
        self.setup.mecsetup=False
        rc=self.setup.run(tasklist=self.t if not full else None)
        self.check()
        return rc

class Finish(Update):
    """ finish processing with localization and association
       Make tsmaps for problem fits
    """
    def defaults(self):
        return dict(dampen=0,
            localize=True, 
            tsmap_dir='"tsmap_fail"',  
            setup_cmds = 'from uw.like2.pipeline import associate ',
            associator="associate.SrcId('$FERMI/catalog','all_but_gammas')",quiet=True)
            
class Tables(Update):
    """ create standard Tables """
    skyfuns=[("CountsMap", "counts", {}), ("KdeMap", "kde", {}), ("ResidualTS", "ts", dict(photon_index=2.2))]
    nside=512
    def defaults(self):
        return dict(dampen=0,
            processor="processor.table_processor",
            tables="""maps.ROItables("%(outdir)s", nside=%(nside)s, skyfuns=%(skyfuns)s)""" %\
                        dict(skyfuns=self.skyfuns, outdir=self.outdir, nside=self.nside),
            setup_cmds= 'from uw.like2.pipeline import maps', quiet=True
            )
            
class PulsarLimitTables(Tables):
    """ create tables with pulsar fits (used for 2PC limits) """
    skyfuns=[ ("ResidualLikelihood", "pulsar_like", dict(model="Models.ExpCutoff(p=[1e-15,1.7, 2000.])"),)]

class PulsarDetection(Tables):
    """ create TS table with pulsar-like function"""
    skyfuns = [("ResidualTS", "pts", dict(model="Models.ExpCutoff(p=[6e-14, 1.2, 2000])"),)]
            
class ModelTables(Tables):
    """ create Model prediction tables """
    nside=64
    skyfuns = [("ModelMap", "model", dict(nside=64)),]
    
class Create(Update):
    """ create a new model, assuming appropriate config.txt
        model_dir points to the new model, which must have an entry "input_model" in its config.txt
        also look for keys auxcat and skymodel_kw
    """
    def __init__(self, analysisdir='.', model_dir='.',  **kwargs):
        """
        
        """
        self.analysisdir = os.path.expandvars(analysisdir)
        os.chdir(self.analysisdir)
        config = eval(open(os.path.join(analysisdir,model_dir, 'config.txt')).read())
        try:
            input_model = config['input_model']
            if type(input_model)==types.StringType:
                self.indir=os.path.join(analysisdir,model_dir,input_model)
                skykw = config.pop('skymodel_kw', None)
                auxcat= config.pop('auxcat', '')
            else:
                self.indir=os.path.join(analysisdir,model_dir,input_model['path'])
                skykw = input_model.get('skymodel_kw', None)
                auxcat= input_model.get('auxcat', '')

        except KeyError:
            print '"input_model" not found in config: %s' %config
            raise
        self.outdir=model_dir
        kw = self.defaults()
        
        # if config has skymodel_kw, will add to the skymodel keywords by interpreting it and setting skymodel_extra
        # especially, for a filter entry. Also possible: update_positions=10 (or other TS)
        skymodel_extra = ','.join(['%s=%s' % item for item in skykw.items() ]) if skykw is not None else ''
        
        kw.update(outdir=model_dir, datadict=config['datadict'],
            irf=config['irf'],
            diffuse=config['diffuse'],
            auxcat=auxcat,
            skymodel_extra=skymodel_extra,
            )
        kw.update(kwargs)
        self.setup = Setup(self.indir,  **kw)

    def defaults(self):
        return dict(dampen=1.0,
            )
            
