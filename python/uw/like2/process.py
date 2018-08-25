"""
Classes for pipeline processing
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/process.py,v 1.31 2018/01/27 15:37:17 burnett Exp $

"""
import os, sys, time, glob
import cPickle as pickle
import numpy as np
import pandas as pd
from scipy import optimize
from skymaps import SkyDir, Band
from uw.utilities import keyword_options
from uw.like2 import (main, tools, sedfuns, maps, sources, localization, roimodel, seeds,)



class Process(main.MultiROI):

    defaults=(
        ('outdir',        '.',   'output folder'),
        ('load_kw', {'rings':2}, 'a dict specific for the loading'),
        ('config_kw',    {},     'additional configuration keywords'),
        ('localize_flag', False, 'perform localiation'),
        ('localize_kw',   {},    'keywords for localization'),
        ('repivot_flag',  False, 'repivot the sources'),
        ('betafix_flag',  False, 'check betas, refit if needed'),
        ('ts_beta_zero',  None,   'set to threshold for converting LP<->PL'),
        ('ts_min',        5,      'Minimum TS for saving after an iteration'),
        ('dampen',        1.0,   'damping factor: set <1 to dampen, 0 to not fit'),
        ('selected_pars', None,   'Apply to fit, to select subset of parameters for global fit'),
        ('counts_dir',    None,  'folder for the counts plots'),
        ('norms_only',    False, 'fit to norms only'),
        ('fix_spectra_flag',False,  'set variable sources to fit norms only, for this and subsequent iterations'),
        ('countsplot_tsmin', 100, 'minimum for souces in counts plot'),
        ('source_name',   None,   'for localization?'),
        ('fit_kw',        dict(ignore_exception=True),     'extra parameters for fit'),
        ('associate_flag',False,  'run association'),
        ('tsmap_dir',     None,   'folder for TS maps'),
        ('sedfig_dir',    None,   'folder for sed figs'),
        ('quiet',         False,  'Set false for summary output'),
        ('finish',        False,  'set True to turn on all "finish" output flags'),
        ('residual_flag', False,  'set True for special residual run; all else ignored'),
        ('diffuse_key',   None,   'set to "gal" or "iso" to evaluate diffuse spectral corrections'),
        ('profile_flag',  False,    'create profile entries for all free sources'),
        ('tables_flag',   False,  'set True for tables run; all else ignored'),
        #('xtables_flag',  False,  'set True for special tables run; all else ignored'),
        ('tables_nside',  512,    'nside to use for table generation'),
        ('table_keys',    None,   'list of keys for table generation: if None, all else ignored'),
        ('seed_key',      None,   'set to name of key for seed check run'),
        ('update_positions_flag',False,  'set True to update positions before fitting'),
        ('add_seeds_flag', False,    'Add seeds found within the ROI, from the table, "plots/seedcheck/good_seeds.csv"'), 
        ('special_flag',  False,  'set for special processing: invoke member func "special"'), 
        ('model_counts',  None,   'set to run model counts'),
        
    )
    
    @keyword_options.decorate(defaults)
    def __init__(self, config_dir, roi_list=None, **kwargs):
        """ process the roi object after being set up
        """
        keyword_options.process(self, kwargs)
        if self.finish:
            self.__dict__.update(dampen=0,
                                 localize_flag=True, associate_flag=True,
                                 sedfig_dir='sedfig',
                                 counts_dir='countfig', 
                                 tsmap_dir='tsmap_fail',
                                 profile_flag=True,
                                 )
        super(Process,self).__init__(config_dir,quiet=self.quiet,)
        self.stream = os.environ.get('PIPELINE_STREAMPATH', 'interactive')
        #if self.xtables_flag:
        #    # suppress loading all point sources for from-scratch source finding
        #    #self.load_kw={'rings':-1} # not now
        #    pass
        if roi_list is not None:
            for index in roi_list:
                self.process_roi(index)
        
    def process_roi(self, index):
        """ special for batch: manage the log file, make sure an exception closes it
        """
        print 'Setting up ROI #%04d ...' % index,
        sys.stdout.flush()
        self.setup_roi(index)
        if self.outdir is not None: 
            logpath = os.path.join(self.outdir, 'log')
            if not os.path.exists(logpath): os.mkdir(logpath)
            outtee = tools.OutputTee(os.path.join(logpath, self.name+'.txt'))
        else: outtee=None
        print 'Processing...'
        sys.stdout.flush()
        try:
            self.process()
        finally:
            if outtee is not None: outtee.close()
        
        
    def process(self):
        roi=self
        dampen=self.dampen 
        outdir = self.outdir
        print  '='*80
        print '%4d-%02d-%02d %02d:%02d:%02d - %s - %s' %(time.localtime()[:6]+ (roi.name,)+(self.stream,))

        # special processing flags
        if self.diffuse_key is not None and self.diffuse_key!='post_gal':
            if self.diffuse_key=='iso':
                fit_isotropic(self)
                return
            elif self.diffuse_key=='gal':
                if not fit_galactic(self):return
            elif self.diffuse_key=='gal_only':
                fit_galactic(self)
                return
            else:
                raise Exception('Unexpected key: {}'.format(self.diffuse_key))
             #fit_diffuse()
            
        if self.residual_flag:
            self.residuals()
            return        
        if self.tables_flag:
            self.tables()
            return
        if self.table_keys is not None:
            self.tables(mapkeys=self.table_keys)
            return
        #if self.xtables_flag:
        #    self.tables(special=True)
        #    return
        if self.seed_key is not None:
            key = self.seed_key
            if not seeds.add_seeds(self, key, config=self.config) :
                # nothing added, so nothing to do with the model for this ROI
                write_pickle(self) # make sure to update anyway
                return
        if self.model_counts is not None:
            maps.ModelCountMaps(self, bandlist=self.model_counts, subdir='model_counts' )
            return

        if self.special_flag:
            """ special processing, converting something."""
            # if not self.model_count_maps(): #fit_second_order() : 
            #     print '-----special processing had nothing to to'
            psc_check(self)
            return
        if self.add_seeds_flag:
            self.add_sources()
        
        if self.fix_spectra_flag:
            # special for monthly or smaller processing
            fix_spectra(self)
            
        if self.counts_dir is not None and not os.path.exists(self.counts_dir) :
            try: os.makedirs(self.counts_dir) # in case some other process makes it
            except: pass
        sys.stdout.flush()
        init_log_like = roi.log_like()
        if self.update_positions_flag:
            self.update_positions()
            
        roi.print_summary(title='before fit, logL=%0.f'% init_log_like)
        fit_sources = [s for s in roi.free_sources if not s.isglobal]
        if len(roi.sources.parameters[:])==0 or dampen==0:
            print '===================== not fitting ========================'
        else:
            fit_kw = self.fit_kw
            try:
                if self.norms_only:
                    print 'Fitting parameter names ending in "Norm"'
                    roi.fit('_Norm',  **fit_kw)
                roi.fit(select=self.selected_pars, update_by=dampen, **fit_kw)
                if self.fix_spectra_flag:
                    # Check for bad errors, 
                    diag = np.asarray(self.hessian().diagonal())[0]
                    if np.any(diag<0):
                        print 'Retrying bad fits, reset '
                        for i,v in enumerate(diag):
                            if v>0: continue
                            self.fit([i], setpars={i: -13}, **fit_kw)
                        self.fit(**fit_kw)
                
                change =roi.log_like() - init_log_like 
                if  abs(change)>1.0 :
                    roi.print_summary(title='after global fit, logL=%0.f, change=%.1f'% (roi.log_like(), change))
                    
                if self.repivot_flag:
                    # repivot, iterating a few times
                    n = 3
                    while n>0:
                        if not self.repivot( fit_sources, select=self.selected_pars): break
                        n-=1
                if self.betafix_flag:
                    if not self.betafix(ts_beta_zero=self.ts_beta_zero):
                        print 'betafix requested, but no refit needed, quitting'

                if self.diffuse_key=='post_gal':
                    fit_galactic(self)
            except Exception, msg:
                print '============== fit failed, no update!! %s'%msg
                raise
        
        def getdir(x ):
            if x is None or outdir is None: return None
            t = os.path.join(outdir, x)
            if not os.path.exists(t): os.mkdir(t)
            return t
        sedfig_dir = getdir(self.sedfig_dir)
        if sedfig_dir is not None:
            print '------------ creating seds, figures ---------------'; sys.stdout.flush()
            skymodel_name = os.path.split(os.getcwd())[-1]
            roi.plot_sed('all', sedfig_dir=sedfig_dir, suffix='_sed_%s'%skymodel_name, )
        
        if self.profile_flag:
            print '------creating profile entries for all free sources'; sys.stdout.flush()
            self.profile('all')

        if self.localize_flag:
            print '------localizing all local sources------'; sys.stdout.flush()
            tsmap_dir = getdir(self.tsmap_dir)
            skymodel = os.getcwd().split('/')[-1]
            if skymodel.startswith('month') or skymodel.startswith('year'): 
                print 'Not running tsmap analysis since data subset'
                tsmap_dir=None
            roi.localize('all', tsmap_dir=tsmap_dir)
        
        if self.associate_flag:
            print '-------- running associations --------'; sys.stdout.flush()
            self.find_associations('all')

        print '-------- analyzing counts histogram, ',; sys.stdout.flush()
        cts=roi.get_count_dict() # always do counts
        print 'chisquared: %.1f ----'% cts['chisq']

        counts_dir = getdir(self.counts_dir)
        if counts_dir is not None:
            print '------- saving counts plot ------'; sys.stdout.flush()
            try:
                fig = roi.plot_counts( tsmin=self.countsplot_tsmin)
                fout = os.path.join(counts_dir, ('%s_counts.jpg'%roi.name) )
                print '----> %s' % fout ; sys.stdout.flush()
                fig.savefig(fout, dpi=60)
            except Exception,e:
                print '***Failed to analyze counts for roi %s: %s' %(roi.name,e)
                chisq = -1
        
        if outdir is not None:  
            write_pickle(self)
            #pickle_dir = os.path.join(outdir, 'pickle')
            #if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
            #roi.to_healpix( pickle_dir, dampen, 
            #    counts=cts,
            #    stream=self.stream,
            #    )
    
   
    def repivot(self, fit_sources=None, min_ts = 10, max_beta=3.0, emin=200., emax=100000.,
             dampen=1.0, tolerance=0.10, test=False, select=None):
        """ invoked  if repivot flag set;
        returns True if had to refit, allowing iteration
        """
        roi = self
        if fit_sources is None:
            fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
        if len(fit_sources)>1:
            print '\ncheck need to repivot % sources with TS>%.0f, beta<%.1f: \n'\
            'source                     TS        e0      pivot' % (len(fit_sources), min_ts, max_beta)
        
        need_refit =False


        for source in fit_sources:
            model = source.spectral_model
            try:
                ts, e0, pivot = roi.TS(source.name),model.e0, model.pivot_energy()
            except Exception, e:
                print 'source %s:exception %s' %(source.name, e)
                continue
                
            if pivot is None: 
                print 'pivot is none'
                continue
            if model.name=='LogParabola': e0 = model[3]
            elif model.name=='ExpCutoff' or model.name=='PLSuperExpCutoff':
                e0 = model.e0
            else:
                print 'Model %s is not repivoted' % model.name
                continue
            print '%-20s %8.0f %9.0f %9.0f '  % (source.name, ts, e0, pivot),
            if ts < min_ts: 
                print 'TS too small'
                continue
            if model.name=='LogParabola':
                if model[2]>max_beta: 
                    print 'beta= %.2f too large' %(model[2])
                    continue #very 
            if pivot < emin or pivot > emax:
                print 'pivot energy, not in range (%.0f, %.0f): setting to limit' % (emin, emax)
                pivot = min(emax, max(pivot,emin))
                model.set_e0(pivot)
                limited = True
            else: limited=False
            if abs(pivot/e0-1.)<tolerance: #0.05:
                print 'converged'; continue
            print 'will refit'
            need_refit=True
            if not test and not limited: model.set_e0(pivot*dampen+ e0*(1-dampen))
        if need_refit and not test:
            roi.fit(select=select, tolerance=0, ignore_exception=True)
        return need_refit

    def betafix(self, ignore_exception=True, ts_beta_zero=9):
        """ evalute ts_beta for all sources, add to source info
            ts_beta_zero: float or None
                if a float, convert source so that log parabola are greater
        """
        for source in self.free_sources:
            if source.isglobal: continue #skip global
            print '----------------- %s (%.1f)-------------' % (source.name, source.ts)
            t=source.ts_beta = self.ts_beta(source.name, ignore_exception=ignore_exception)
            if t is None: continue
            print 'ts_beta ', t,
            if ts_beta_zero is None: 
                print ' -- not checking'
                continue
            changed=True
            powerlaw = not source.model.free[2]
            if t<ts_beta_zero and not powerlaw:
                self.freeze('beta', source.name, 0.)
                print '--> PowerLaw'
            elif t>ts_beta_zero and powerlaw:
                self.thaw('beta', source.name)
                print '--> LogParabola'
            else:
                print ': OK'
                changed=False
            if changed:
                self.fit(source.name)
        self.fit(ignore_exception=ignore_exception) # seems necessary
        return False # nofollowup

        
    def residuals(self, tol=0.3):
        print 'Creating tables of residuals'
        if not os.path.exists('residuals'):
            os.mkdir('residuals')
        resids = sedfuns.residual_tables(self, tol)
        filename = 'residuals/%s_resids.pickle' %self.name
        with open(filename, 'w') as out:
            pickle.dump(resids, out)
            print 'wrote file %s' %filename
            
    def tables(self, special=False, mapkeys=['ts', 'kde']):
        """create a set of tables"""
        if mapkeys==['rmap']:
            # residual maps
            maps.residual_maps(self)
            return
        tinfo = [maps.table_info[key] for key in mapkeys]
        skyfuns = [(entry[0], key, entry[1]) for key,entry in zip(mapkeys, tinfo)]  
        rt = maps.ROItables(self.outdir, nside=self.tables_nside, skyfuns=skyfuns )
        rt(self)
        
    def update_positions(self, tsmin=10, qualmax=8):
        """ use the localization information associated with each source to update position
            require ts>tsmin, qual<qualmax
            
        """
        print '---Updating positions---'
        sources = [s for s in self.sources if s.skydir is not None and np.any(s.spectral_model.free)]
        #print 'sources:', [s.name for s in sources]
        print '%-15s%6s%8s %8s' % ('name','TS','qual', 'delta_ts')
        for source in sources:
            has_ts= hasattr(source, 'ts')
            print '%-15s %6.0f' % (source.name, source.ts if has_ts else -1.0) , 
            if not hasattr(source, 'ellipse') or source.ellipse is None:
                print ' no localization info'
                continue
            if not has_ts:   
                print '    no TS'
                continue

            if source.ts<tsmin:
                print '    TS<%.0f' % (tsmin)
                continue
            newdir = SkyDir(*source.ellipse[:2]); qual, delta_ts = source.ellipse[-2:]
            print '%6.1f%6.1f' % (qual, delta_ts) ,
            if qual>qualmax:
                print ' qual>%.1f' % qualmax
                continue
            print ' %s -> %s, moved %.2f' % (source.skydir,newdir, np.degrees(newdir.difference(source.skydir)))
            source.skydir = newdir
    
    
    def fit_second_order(self, summarize=False):
        """
        Fit the second order parameter (beta or Cutoff) for all variable sources
        Leave them frozen.
        """
    
        def fit2(source_name, parname='beta', fmt='{:6.2f}'):
            s=self.get_source(source_name)
            sources.set_default_bounds(s.model, True)# force change of bounds
            self.thaw(parname)
            parval = s.model[parname]
            try:
                self.fit(s.name, summarize=summarize, estimate_errors=True, ignore_exception=False)
            except Exception, msg:
                print 'Failed to fit {} for {}: {}'.format(parname,source_name, msg)
                s.model[parname]=parval
            self.freeze(parname)
            print ('{:15}{:8.0f}'+fmt+fmt).format(source_name, s.ts, s.model[parname], s.model.error(parname))
        
        LP_sources = [s.name for s in self.free_sources 
            if not s.isglobal and s.model.name=='LogParabola']
        PLEX_sources = [s.name for s in self.free_sources 
            if not s.isglobal and s.model.name=='PLSuperExpCutoff']

        print '{:15}{:>8} {:6} {}'.format('LP source', 'TS', 'beta', 'error')
        map( lambda s: fit2(s,'beta'), LP_sources) 
        print '{:15}{:>8} {:6} {:10}'.format('PLEX source', 'TS','Cutoff', 'error')
        map( lambda s: fit2(s,'Cutoff', '{:6.0f}'), PLEX_sources)
        return True

    def model_count_maps(self):
        maps.ModelCountMaps(self, nbands=12, subdir='model_counts' )
        return False

def fix_spectra(roi):
    for src in roi.free_sources:
        m=src.model
        if src.name=='isotrop':
            print 'Freezing isotrop'
            roi.freeze('Scale', src.name, 1.0)
            continue
        
        for i,parname in enumerate(m.param_names[1:]):
            if m.free[i+1]:
                roi.freeze(parname, src.name)
        src.fixed_spectrum=True


class BatchJob(Process):
    """special interface to be called from uwpipeline
    Expect current dir to be output dir.
    """
    def __init__(self, **kwargs):
        config_dir= kwargs.pop('config_dir', '.')
        roi_list = kwargs.pop('roi_list', range(1,3)) 
        kwargs['outdir']= os.getcwd()
        super(BatchJob, self).__init__(config_dir, **kwargs)
        
    def __call__(self, roi_index):
        self.process_roi(roi_index)
    

def fit_isotropic(roi, nbands=8, folder='isotropic_fit'):
    """ fit only the front and back"""
    from uw.like import Models
    iso_source = roi.get_source('isotrop')
    old_model=iso_source.model.copy()
    roi.sources.set_model(Models.FrontBackConstant(), 'isotrop')
    iso_model=iso_source.model
    roi.reinitialize()

    cx = []
    for eband in range(nbands):
        print '*** Energy Band {}: iso counts {}'.format( eband,
                                [t[1].counts.round() for t in roi[2*eband:2*eband+2]])
        roi.select(eband)
        iso_model[0]=1; iso_model[1]=1
        roi.fit([0,1]); 
        u = iso_model.get_all_parameters();
        cx.append(u)
    roi.select()
    roi.sources.set_model(old_model)
    if folder is not None:
        if not os.path.exists(folder): os.mkdir(folder)
        filename= '{}/{}.pickle'.format(folder, roi.name)
        pickle.dump(np.array(cx), open(filename, 'w'))
        print 'wrote file {}'.format(filename)
    return np.array(cx)

class FitGalactic(object):
    """Manage the galactic correction fits
    """
    def __init__(self, roi, nbands=8, folder=None, upper_limit=5.0):
        """ fit only the galactic normalization"""
        if folder is not None and not os.path.exists(folder):
            os.mkdir(folder)
        self.roi=roi
        gal_model = roi.get_source('ring').model
        roi.thaw('Norm')
        gal_norm = gal_model[0]
        #gal_model.set_limits('Norm', 0.2, 5.0) # override default
        if upper_limit is not None:
            gal_model.bounds[0][1]=np.log10(upper_limit)# set upper limit by hand
        roi.reinitialize()
        cx = []
        for eband in range(nbands):
            print '*** Energy Band {}: gal counts {}'.format( eband,
                                    [t[0].counts.round() for t in roi[2*eband:2*eband+2]])
            roi.select(eband)
            gal_model[0]=gal_norm
            roi.fit([0], ignore_exception=True); 
            cx.append((gal_model[0], gal_model.error(0)))
        self.fitpars= cx = np.array(cx) # convert to array, shape (nbands, 2)
        x,s = cx.T
        self.chisq= sum(((x-1)/s)**2)
        # re-select all bands, freeze galactic again    
        roi.select()
        roi.freeze('Norm')

        if folder is not None:
            filename= '{}/{}.pickle'.format(folder, roi.name)
            pickle.dump(cx[:,0], open(filename, 'w'))
            print 'wrote file {}'.format(filename)
    
    def update(self):
        from uw.like2 import (response,diffuse)
        r = self.roi
        if r.sources.diffuse_normalization is None:
            print 'FitGalactic: Setting up diffuse normalization'
            roi_index= int(r.name[-4:])
            dn = self.create_corr_dict(r.config['diffuse'], roi_index)
            r.sources.diffuse_normalization = diffuse.normalization = dn

        a = self.roi.sources.diffuse_normalization
        b = self.fitpars[:,0]
        before = a['gal']
        a['gal'] = before * self.fitpars[:,0]
        print before, '\n',a['gal']
        # update the Galactic Response objects
        for gr in self.roi[:16]:
            gr[0].initialize(force=True)

    def create_corr_dict(self, diffuse_dict,  roi_index, event_type_names=('front','back')):
        import response
        corr_dict = {}
        galf = diffuse_dict['ring']['correction']
        corr_dict['gal'] = response.DiffuseCorrection(galf).roi_norm(roi_index)

        isof =  diffuse_dict['isotrop']['correction']
        corr_dict['iso']= dict()
        for x in event_type_names:
            isoc = response.DiffuseCorrection(isof.replace('*',x));
            corr_dict['iso'][x]= isoc.roi_norm(roi_index)
        return corr_dict

def fit_galactic(roi, nbands=8, folder=None, upper_limit=5.0):
    t = FitGalactic(roi, nbands, folder, upper_limit)
    print 'Chisq: {:.1f}'.format(t.chisq)
    t.update()
    return True
    
def fit_diffuse(roi, nbands=8, select=[0,1,2], restore=False, folder='diffuse_fit'):
    """
    Perform indpendent fits to the gal, iso_front, and iso_back for each of the first nbands bands
    select: None or list of variables
    """
    from uw.like import Models
    # freeze all free sources, thaw gal and iso
    roi.thaw('Norm', 'ring')
    iso_model =roi.sources.find_source('isotrop').model.copy()
    roi.sources.set_model(Models.FrontBackConstant(), 'isotrop')
    roi.reinitialize()
    
    # do the fitting
    dpars=[]
    energies = []
    for ie in range(nbands):
        roi.select(ie); 
        roi.fit(select, ignore_exception=True)
        energies.append(int(roi.energies[0]))
        dpars.append( roi.sources.parameters.get_parameters()[:3])
    t =np.power(10, dpars)
    df = pd.DataFrame(t, columns=['gal iso_front iso_back'.split()])
    df.index=energies
    
    if restore:
        # does not seem to work, comment this out for now
        # restore sources
        roi.sources.diffuse_normalization *= df
        for s,f in zip(free_sources, saved_free):
            s.model.free=f
        roi.sources.find_source('ring').model.free[0]=False
        roi.sources.set_model(iso_model, 'isotrop')
        roi.reinitialize()
        roi.select()
        roi.fit() # needed to restore gradient, at least.
        
        # update the pickle file
        write_pickle(roi)
    elif folder is not None:
        # simply save results
        if not os.path.exists(folder):
            os.mkdir(folder)
        filename= '{}/{}.pickle'.format(folder, roi.name)
        pickle.dump(df, open(filename, 'w'))
        print 'wrote file {}'.format(filename)
    return df
    
def write_pickle(roi):
    pickle_dir = os.path.join(roi.outdir, 'pickle')
    if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
    roi.to_healpix( pickle_dir, dampen=1.0, 
        counts=roi.get_count_dict(),
        stream=os.environ.get('PIPELINE_STREAMPATH', 'interactive'),
        ts_min = roi.ts_min,
        )


def psc_check(roi, psc_name='gll_psc*uw8011*' , outdir='psc_check'):
    """Compare the spectra of sources from a "gll" file with the corresponding
    pointlike original fits.
    """

    from uw.like2.analyze import fermi_catalog
    from uw.like2.plotting import sed

    #load the catalog
    fgl = roi.config.get('fgl', None)
    if fgl is None:
        fgl = fermi_catalog.GLL_PSC2(psc_name)
        roi.config['fgl']= fgl
    
    # Replace sources, keep list of (old,new)
    changed = []
    source_pairs=[]
    for s in roi.free_sources:
        cname=s.name.replace(' ','')
        if cname not in fgl.df.index: 
            print '{:14s} {:6.0f} '.format(s.name, s.ts)
            continue
        fl8y = fgl.df.loc[cname]
        trunc = fl8y.sname[5:]
        print '{:14s} {:6.0f} --> {:14s} {:6.0f}'.format(s.name, s.ts, trunc, fl8y.ts),
        if fl8y.extended:
            sx = sources.ExtendedSource(name=trunc, skydir=(fl8y.ra,fl8y.dec), 
                model=fl8y.model, dmodel=s.dmodel)
        else:
            sx = sources.PointSource(name=trunc, skydir=(fl8y.ra,fl8y.dec), model=fl8y.model)
        roi.del_source(s.name)
        roi.add_source(sx)
        newts=roi.TS(trunc)
        print ' -->{:6.0f}'.format(newts)
        changed.append(( (s.name,s.model,s.ts),
            (sx.name,sx.model,sx.ts))) #avoid saving dmodel?
        source_pairs.append((s, sx))
        
    if outdir is None: return
    # save info for comparison
    def path_check(x):
        if not os.path.exists(x): os.mkdir(x)

    map( path_check, [outdir, outdir+'/info', outdir+'/sed'])        
    pickle.dump(changed, open('psc_check/info/HP12_{}.pickle'.format(roi.name[5:]),'w'))

    # save plots
    for old,new in source_pairs:
        sed.plot_pair(old,new).savefig('psc_check/sed/{}_sed.jpg'.format(new.name.replace('+','p')))
    

#### TODO
#def run(rois, **kw):
#    Process('.', rois, **kw )
#    
#if __name__=='__main__':
#    run()
