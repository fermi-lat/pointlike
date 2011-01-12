"""
Main entry for the UW all-sky pipeline
$Header$
"""
import os, types, glob, time
import skymodel, skyanalysis, processor, associate
from uw.utilities.assigntasks import setup_mec, AssignTasks, get_mec, kill_mec

class Pipe(skyanalysis.SkyAnalysis):
    """ This is a subclass of SkyAnalysis, itself a subclass of the basic pointlike
    setup class, SpectralAnalysis. It implements the capabilities needed for 
    IPython paralization: the constructur defines the environment, and if called with
    an index, it will run a full analysis for that ROI.
    
    indir : string
        name of a folder containing the sky model description, passed to 
           skymodel.SkyModel
           
    """

    def __init__(self, indir, **kwargs):
        self.nside      = kwargs.pop('nside', 12)
        self.process_kw = kwargs.pop('process_kw', dict())
        self.fit_kw     = kwargs.pop('fit_kw', dict())
 
        associator = kwargs.pop('associate', 'all_but_gammas')
        if associator is not None:
            associator = associate.SrcId('$FERMI/catalog', associator)
            print 'will associate with catalogs %s' % associator.classes
        self.process_kw['associate'] = associator

        sm = skymodel.SkyModel(indir)
        super(Pipe, self).__init__(sm, **kwargs)
        
    def __call__(self, index):
        """ perform analysis for the ROI
        Note that it is done by a function in the processor module
        """
        roi = self.roi(index)
        processor.process(roi, **self.process_kw)

    def names(self):
        return ['HP12_%04d'%i for i in range(12*12**2)]

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
                outdir=outdir,
                dataset = '24MP7source',
                #diffuse = (config.diffuse_path, 'ring_24month_P74_v1.fits', 'isotrop_21month_v2.txt'),
                irf = 'P7SOURCE_V4PSF',
                associator ='all_but_gammas',
                sedfig = None,
                tsmap  = None,
                localize=True,
                fit_emin=100,
                fit_emax=800000,
                fix_beta=False, dofit=True,
                source_kw=dict(),
                fit_kw=dict(use_gradient=True,),
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
from uw.pipeline.pipe import Pipe;
g=Pipe("%(indir)s", 
        irf="%(irf)s",
        fit_emin=%(fit_emin)s, fit_emax=%(fit_emax)s, minROI=5, maxROI=5,
        associate="%(associator)s",
        process_kw=dict(outdir="%(outdir)s",
            tsmap_dir=%(tsmap)s,  sedfig_dir=%(sedfig)s,
            localize=%(localize)s,
            fix_beta= %(fix_beta)s, dofit=%(dofit)s,
            tables= %(tables)s,
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
    
def iterate(setup, **kwargs):
    """ 
    """
    version = setup.version
    main(setup(), setup.outdir, **kwargs)
    setup.increment_version()
 
#--------------------        
def main( setup_string, outdir, mec=None, startat=0, n=0, local=False,
        machines=[], engines=None,
        ignore_exception=False, 
        logpath='log',
        progress_bar=False):
        
    if not os.path.exists(outdir): 
        os.mkdir(outdir)
        #raise Exception('outdir folder, "%s", not found' % outdir)
    print 'writing results to %s' %outdir
    if setup_string[0]=='@':
        setup_string = open(setup_string).read()
    g = getg(setup_string)
    print >> open(os.path.join(outdir, 'config.txt'), 'w'), g
    print >> open(os.path.join(outdir, 'setup_string.txt'), 'w'), setup_string
    names = g.names(); 
    print 'Start at source %d' % startat
    if len(names)==0:
        raise InvalidArgument('no tasks defined by Pipeline object')
    else:
        print 'found %d sources to process' % len(names)

    if n==0: endat = len(names)
    else: endat = min(startat+n, len(names))
         
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
            
    if not local: setup_mec(machines=machines, engines=engines)
    time.sleep(10)
    lc= AssignTasks(setup_string, tasks[startat:endat], mec=mec, 
        timelimit=2000, local=local, callback=callback, 
        ignore_exception=ignore_exception,
         progress_bar=progress_bar)
    lc(5)
    if not local: get_mec().kill(True)
    return lc
   

