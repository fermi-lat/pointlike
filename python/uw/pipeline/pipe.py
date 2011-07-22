"""
Main entry for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pipe.py,v 1.14 2011/06/24 04:54:49 burnett Exp $
"""
import os, types, glob, time
import cPickle as pickle
import numpy as np
from . import skymodel, skyanalysis, processor, associate, bandlike
from uw.utilities import makerec
from uw.utilities.assigntasks import setup_mec, AssignTasks, get_mec, kill_mec, free

class Pipe(skyanalysis.SkyAnalysis):
    """ This is a subclass of SkyAnalysis, itself a subclass of the basic pointlike
    setup class, SpectralAnalysis. It implements the capabilities needed for 
    IPython paralellization: the constructur defines the environment, and if called with
    an index, it will run a full analysis for that ROI.
    
    indir : string
        name of a folder containing the sky model description, passed to 
           skymodel.SkyModel
           
    """
    def __init__(self, indir, dataset,  **kwargs):
        """
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
 
        associator = kwargs.pop('associate', 'all_but_gammas')
        if associator is not None and associator[0]=='$':
            print 'associations will come from the catalog %s' %associator
            associator = associate.Association(associator)
        elif associator is not None and associator!='None':
            associator = associate.SrcId('$FERMI/catalog', associator)
            #print 'will associate with catalogs %s' % associator.classes
        self.process_kw['associate'] = associator
        self.selector = skymodel.HEALPixSourceSelector
        self.selector.nside = self.nside

        # the function that will be called for each ROI: default is processor.process
        self.processor = kwargs.pop('processor', processor.process)
        if type(self.processor)==types.StringType:
            self.processor = eval(self.processor)

        # create the SkyModel object passing any additional keywords to it, and start the superclass
        sm = skymodel.SkyModel(indir,  **self.skymodel_kw)
        super(Pipe, self).__init__(sm, dataset, **kwargs)
        
    def __call__(self, index):
        """ perform analysis for the ROI
        Note that it is done by a function in the processor module
        """
        roi = self.roi(index )
        self.processor(roi, **self.process_kw)

    def names(self):
        return [self.selector(i).name() for i in range(12*self.nside**2)]

class Setup(dict):
    """ Setup is a dictionary with variable run parameters
    """

    def __init__(self, version=None,  **kwargs):
        """ generate setup string"""
        if version is None: version = int(open('version.txt').read())
        indir=kwargs.pop('indir', 'uw%02d'%(version)) 
        self.outdir=outdir=kwargs.pop('outdir', 'uw%02d'%(version+1))
        self.version=version
        if os.name=='nt':
            os.system('title %s %s'% (os.getcwd(), indir))
        self.update(dict( cwd =os.getcwd(), 
                indir = indir,
                auxcat='',
                nside = 12,
                outdir=outdir,
                pass_number=0,
                dataset = 'P7_V4_SOURCE',
                diffuse = ('ring_24month_P76_v1.fits', 'isotrop_21month_P76_v2.txt'),
                extended= None,
                alias= {}, 
                skymodel_extra = '',
                irf = 'P7SOURCE_V4PSF',
                associator ='all_but_gammas',
                sedfig = None,
                tsmap  = None, tsfits=False,
                localize=True,
                fit_emin=100,
                fit_emax=800000,
                processor='processor.process',
                fix_beta=False, dofit=True,
                source_kw=dict(),
                fit_kw=dict(),
                repivot = True,
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
import os, pickle; os.chdir(r"%(cwd)s");%(setup_cmds)s
from uw.pipeline import pipe, maps, skymodel;
g=pipe.Pipe("%(indir)s", "%(dataset)s",  event_class=0, 
        skymodel_kw=dict(auxcat="%(auxcat)s",diffuse=%(diffuse)s,
            extended_catalog_name="%(extended)s", update_positions=%(update_positions)s,
            free_index=%(free_index)s,
            alias =%(alias)s, %(skymodel_extra)s), 
        irf="%(irf)s", nside=%(nside)s,
        fit_emin=%(fit_emin)s, fit_emax=%(fit_emax)s, minROI=%(minROI)s, maxROI=%(maxROI)s,
        processor="%(processor)s",
        associate="%(associator)s",
        process_kw=dict(outdir="%(outdir)s", pass_number=%(pass_number)s,
            tsmap_dir=%(tsmap)s,  tsfits=%(tsfits)s,
            sedfig_dir=%(sedfig)s,
            localize=%(localize)s,
            fix_beta= %(fix_beta)s, dofit=%(dofit)s,
            tables= %(tables)s,
            repivot=%(repivot)s, dampen=%(dampen)s,
            ),
        fit_kw = %(fit_kw)s,
    ) 
""" %self
        print 'Pipeline input, output: %s -> %s ' % (indir, outdir)
    def __call__(self): return self.setup_string
    
    def g(self):
        exec(self(), globals()) # should set g
        return g

    def increment_version(self, increment=1):
        out = open('version.txt','w')
        out.write('%2d'%(self.version+increment))
        out.close()
        print 'Increment version.txt to %d' % (self.version+increment) 
 
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
def converge_test(version=None, thresh=10,nside=12):
    if version is None:
        version = int(open('version.txt').read())
    old_drec = roirec(version-1, nside)
    old_loglike = sum(old_drec.loglike)
    new_drec = roirec (version,nside)
    new_loglike= sum(new_drec.loglike)
    deltas = new_drec.loglike - old_drec.loglike
    moved = sum(np.abs(deltas)>thresh)
    big = new_drec.name[(deltas==deltas.max())+(deltas==deltas.min())]
    s='iteration %2d: total log likelihood change: %6.0f\n number changed>%2.0f:'\
        +'%4d, min, max changes %6.1f,%6.1f (%s,%s)'
    print s % (version, (new_loglike - old_loglike),thresh, moved,  deltas.min(),deltas.max(),big[0],big[1])
    chisq=new_drec.chisq
    print 'number with chisq>50: %d, max value %.0f' % (sum(chisq>50), chisq.max())
 
def iterate(setup, **kwargs):
    """ 
    """
    version = setup.version
    main(setup, **kwargs)
    converge_test(version+1, setup['nside'])
    setup.increment_version()

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

def check_converge(month, tol=10, log=None):
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
    
def iterate_months(month, mec, make_setup, minpass=0, maxpass=7, log=None, tol=10, dampen=0.5, **kwargs):
    """
    run main for a given month or dataset. minpass is starting pass number:
        0 : inital run from default skymodel
        1 : update all ROIs in place
        >1: update only changed ROIs and neighbors
    loops until maxpass, or no ROIs to update
    """
   
    outdir = 'month%02d'%month if type(month)==types.IntType else month
    if log is None:
        log = open(os.path.join(outdir, 'log.txt'), 'a')
    converged=False
    iterfilename = os.path.join(outdir, 'iteration.txt')
    if  not os.path.exists(iterfilename):
        iterfile = open(iterfilename, 'w').write('%d'%minpass)
    minpass = int(open(iterfilename).read())
    for pass_number in range(minpass, maxpass):
        print >>log, time.asctime(), 'start iteration %d, damping factor=%.2f' % (pass_number, dampen)
        log.flush()
        if pass_number==0:
            setup = make_setup(outdir, False, pass_number)
            main(setup, mec=mec, **kwargs)
        elif pass_number==1:
            setup = make_setup(outdir, True, pass_number)
            main(setup, mec=mec, **kwargs)
        else:
            taskids = check_converge(outdir, log=log, tol=tol)
            if len(taskids)==0: 
                converged=True; break
            setup = make_setup(outdir, True, pass_number, dampen=dampen)
            main(setup, mec=mec, taskids=taskids, **kwargs
            )
        open(iterfilename, 'w').write('%d'%pass_number);
    
    print >>log, time.asctime(), 'finished'
    log.flush()
    print 'done!' if converged else 'Warning: did not converge' 
 
#--------------------        
def main( setup, mec=None, taskids=None, local=False,
        ignore_exception=False, 
        logpath='log',
        progress_bar=False):
    """
    Parameters
        setup : Setup object, implements
            setup_string = setup()
                Python code, must define an object "g". It must implement:
                a function g(n), n an integer
                g.names() must return a list of task names
            outdir =setup.outidr : directory to save files
        mec : None or a MultiEngineClient instance
        taskids : None or a list of integers
            the integers should be in range(0,len(g.names())
    """
    setup_string = setup()
    outdir = setup.outdir
    if outdir is not None:
        if not os.path.exists(outdir): 
            os.mkdir(outdir)
        print 'writing results to %s' %outdir
    # create an instance as a test, and to get the number of tasks that it defines.
    # executing the setup string must create an object g, where g(n) for 0<=n<len(g.names())
    # will execute the task n
    g = getg(setup_string)
    if logpath is not None and outdir is not None:
        print >> open(os.path.join(outdir, 'config.txt'), 'w'), str(g)
        print >> open(os.path.join(outdir, 'setup_string.txt'), 'w'), setup_string
    names = g.names(); 

    if taskids is None: taskids = range(len(names))
    if len(names)==0:
        raise InvalidArgument('no tasks defined by Pipeline object')
    else:
        print 'found %d tasks, will process %d' % (len(names), len(taskids))

    
    # list of function calls to exec for all the tasks     
    tasks = ['g(%d)'% i for i in taskids]
    post = 'del g' #for cleanup
    del g #do not need this instance
    
    def callback(id, result):
        if logpath is None: return
        try:
            name = names[taskids[id]]
            logdir = os.path.join(outdir, logpath) if outdir is not None else logpath
            if not os.path.exists(logdir): os.mkdir(logdir)
            out = open(os.path.join(logdir, '%s.txt' % name.strip().replace(' ','_').replace('+','p')), 'a')
            print >>out, '='*80
            print >>out, '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (name,))
            print >>out, result
            out.close()
        except:
            print 'got exception writing %s' % name
            raise
            
    if not local and mec is None: 
        setup_mec()
        time.sleep(20)
        mec=get_mec()
        assert len(mec.get_ids())>10, 'Failed to setup mec?'
    lc= AssignTasks(setup_string, tasks, post=post, mec=mec, 
        timelimit=5000, local=local, callback=callback, 
        ignore_exception=ignore_exception,
         progress_bar=progress_bar)
    lc(15)
    if not local:
        if mec is None: 
            get_mec().kill(True)
        else:
            get_mec().reset()
        
    return lc
   

