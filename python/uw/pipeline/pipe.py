"""
Main entry for the UW all-sky pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pipe.py,v 1.3 2011/01/24 20:25:58 kerrm Exp $
"""
import os, types, glob, time, pickle
import numpy as np
from . import skymodel, skyanalysis, processor, associate
from uw.utilities import makerec
from uw.utilities.assigntasks import setup_mec, AssignTasks, get_mec, kill_mec

class Pipe(skyanalysis.SkyAnalysis):
    """ This is a subclass of SkyAnalysis, itself a subclass of the basic pointlike
    setup class, SpectralAnalysis. It implements the capabilities needed for 
    IPython paralellization: the constructur defines the environment, and if called with
    an index, it will run a full analysis for that ROI.
    
    indir : string
        name of a folder containing the sky model description, passed to 
           skymodel.SkyModel
           
    """
    def __init__(self, indir, dataset, **kwargs):
        """
        indir : string
            name of a folder containing the sky model description, passed to 
               skymodel.SkyModel
        dataset : instance of a DataSpec object
               
        """
        self.nside      = kwargs.pop('nside', 12)
        self.process_kw = kwargs.pop('process_kw', dict())
        self.fit_kw     = kwargs.pop('fit_kw', dict())
        self.skymodel_kw= kwargs.pop('skymodel_kw', dict())
 
        associator = kwargs.pop('associate', 'all_but_gammas')
        if associator is not None and associator!='None':
            associator = associate.SrcId('$FERMI/catalog', associator)
            #print 'will associate with catalogs %s' % associator.classes
        self.process_kw['associate'] = associator

        # create the SkyModel object passing any additional keywords to it, and start the superclass
        sm = skymodel.SkyModel(indir,  **self.skymodel_kw)
        super(Pipe, self).__init__(sm, dataset, **kwargs)
        
    def __call__(self, index):
        """ perform analysis for the ROI
        Note that it is done by a function in the processor module
        """
        roi = self.roi(index)
        processor.process(roi, **self.process_kw)

    def names(self):
        return ['HP%d_%04d'%(self.nside,i) for i in range(12*self.nside**2)]

class Setup(dict):
    """ Setup is a dictionary with variable run parameters
    """

    def __init__(self, version=None,  **kwargs):
        """ generate setup string"""
        if version is None: version = int(open('version.txt').read())
        indir='uw%02d'%(version) 
        outdir=self.outdir='uw%02d'%(version+1)
        self.version=version
        if os.name=='nt':
            os.system('title %s %s'% (os.getcwd(), indir))
        self.update(dict( cwd =os.getcwd(), 
                indir = indir,
                auxcat='',
                outdir=outdir,
                dataset = 'P7_V4_SOURCE',
                diffuse = ('ring_24month_P74_v1.fits', 'isotrop_21month_v2.txt'),
                irf = 'P7SOURCE_V4PSF',
                associator ='all_but_gammas',
                sedfig = None,
                tsmap  = None,
                localize=True,
                fit_emin=100,
                fit_emax=800000,
                fix_beta=False, dofit=True,
                source_kw=dict(),
                fit_kw=dict(use_gradient=False,),
                repivot = True,
                tables = None,  #roi_maps.ROItables("%(outdir)s", skyfuns=(
                                # (roi_tsmap.TSCalc, 'ts', dict(photon_index=2.0),) 
                                #  (ts_map.KdeMap, "kde", dict()),))
                ))
        self.update(kwargs)
        # first-order replace
        for key in self:
            if type(self[key])== types.StringType and self[key].find('%(')>=0: 
                self[key]=self[key]%self
                self['tables'] = self['tables']%self
                #print 'fix key %s: %s' % (key, self[key])
        self.setup_string =  """\
import os; os.chdir(r"%(cwd)s");
from uw.pipeline import pipe, maps;
g=pipe.Pipe("%(indir)s", "%(dataset)s",
        skymodel_kw=dict(auxcat="%(auxcat)s",diffuse=%(diffuse)s,),
        irf="%(irf)s",
        fit_emin=%(fit_emin)s, fit_emax=%(fit_emax)s, minROI=5, maxROI=5,
        associate="%(associator)s",
        process_kw=dict(outdir="%(outdir)s",
            tsmap_dir=%(tsmap)s,  sedfig_dir=%(sedfig)s,
            localize=%(localize)s,
            fix_beta= %(fix_beta)s, dofit=%(dofit)s,
            tables= %(tables)s,
            repivot=%(repivot)s,
            ),
        fit_kw = %(fit_kw)s,
    ) 
""" %self
        print 'new version: %d, output to %s ' % (version+1, self.outdir)
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

def roirec(version=None):
    if version is None:
        version = int(open('version.txt').read())
    roi_files = glob.glob('uw%02d/pickle/*.pickle'%version)
    roi_files.sort()
    if len(roi_files)<1728:
        t = map(lambda x : int(x[-11:-7]), roi_files)
        missing = [x for x in xrange(1728) if x not in t]
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
def converge_test(version=None, thresh=10,):
    if version is None:
        version = int(open('version.txt').read())
    old_drec = roirec(version-1)
    old_loglike = sum(old_drec.loglike)
    new_drec = roirec (version)
    new_loglike= sum(new_drec.loglike)
    deltas = new_drec.loglike - old_drec.loglike
    moved = sum(np.abs(deltas)>thresh)
    print 'iteration %2d: total log likelihood change: %6.0f, number changed>%2.0f: %4d, min, max changes %8.1f,%8.1f'\
        % (version, (new_loglike - old_loglike),thresh, moved,  deltas.min(),deltas.max())
    chisq=new_drec.chisq
    print 'number with chisq>50: %d, max value %.0f' % (sum(chisq>50), chisq.max())
 
def iterate(setup, **kwargs):
    """ 
    """
    version = setup.version
    main(setup(), setup.outdir, **kwargs)
    converge_test(version+1)
    setup.increment_version()

 
#--------------------        
def main( setup_string, outdir, mec=None, startat=0, n=0, local=False,
        machines=[], engines=None,
        ignore_exception=False, 
        logpath='log',
        progress_bar=False):
        
    if not os.path.exists(outdir): 
        os.mkdir(outdir)
    print 'writing results to %s' %outdir
    # create an instance as a test, and to get the number of tasks that it defines.
    # executing the setup string must create an object g, where g(n) for 0<=n<len(g.names())
    # will execute the task n
    g = getg(setup_string)
    print >> open(os.path.join(outdir, 'config.txt'), 'w'), str(g)
    print >> open(os.path.join(outdir, 'setup_string.txt'), 'w'), setup_string
    names = g.names(); 
    print 'Start at source %d' % startat
    if len(names)==0:
        raise InvalidArgument('no tasks defined by Pipeline object')
    else:
        print 'found %d entries to process' % len(names)

    if n==0: endat = len(names)
    else: endat = min(startat+n, len(names))
    # list of function calls to exec for all the tasks     
    tasks = ['g(%d)'% i for i in range(len(names))]
    del g #do not need this instance
    
    def callback(id, result):
        try:
            name = names[id+startat]
            logdir = os.path.join(outdir, logpath)
            if not os.path.exists(logdir): os.mkdir(logdir)
            out = open(os.path.join(logdir, '%s.txt' % name.strip().replace(' ','_').replace('+','p')), 'a')
            print >>out, '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (name,))
            print >>out, result
            out.close()
        except:
            print 'got exception writing %s' % name
            raise
            
    if not local and mec is None: 
        setup_mec(machines=machines, engines=engines)
    time.sleep(10)
    lc= AssignTasks(setup_string, tasks[startat:endat], mec=mec, 
        timelimit=2000, local=local, callback=callback, 
        ignore_exception=ignore_exception,
         progress_bar=progress_bar)
    lc(5)
    if not local and mec is None: 
        get_mec().kill(True)
    return lc
   

