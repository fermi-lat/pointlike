"""
roi and source processing used by the roi pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/processor.py,v 1.65 2013/09/25 12:51:12 burnett Exp $
"""
import os, time, sys, types, glob
import cPickle as pickle
import numpy as np
import pylab as plt
import pandas as pd
from skymaps import SkyDir, Hep3Vector
from uw.like import srcid  
from uw.utilities import image
from . import associate
from ..plotting import sed, counts 
from .. import localization, sedfuns
#np.seterr(invalid='raise', divide='raise', over='raise')
np.seterr(invalid='raise', divide='raise', over='ignore')
def isextended(source):
    return source.__dict__.get(  'spatial_model', None) is not None
  
class OutputTee(object):
    def __init__(self, logfile):
        self.logstream = open(logfile, 'a')
        self.stdout = sys.stdout
        sys.stdout = self
    def write(self, stuff):
        self.logstream.write(stuff)
        self.stdout.write(stuff)
    def close(self):
        sys.stdout =self.stdout
        self.logstream.close()
    def flush(self):
        self.stdout.flush()
    def set_parent(self, parent):
        self.stdout.set_parent(parent) #needed??
        
def fix_beta(roi, bts_min=20, qual_min=15,):
    """ invoked by process() if fix_beta flag set, but can be run standalone for testing
    if beta=0.001, it has not been tested.
    if beta<0.01 or error exists and  >0.1
    """
    refit=candidate=False
    print 'checking for beta fit: minTS %s, min qual %s...'% (bts_min, qual_min)
    models_to_fit=[]
    for source in roi.sources:
        model = source.spectral_model
        which = source.name
        if not np.any(model.free) or model.name!='LogParabola': continue
        if not candidate:
            print 'name                      beta         band_ts  fitqual'
        candidate=True
        beta = model[2]
        try:
            beta_unc = np.sqrt(model.get_cov_matrix()[2,2])
        except:
            beta_unc=0
        band_ts, ts = roi.band_ts(which), roi.TS(which)
        sbeta = '%5.3f+/-%5.3f' % (beta, beta_unc) if beta_unc>0 else '%5.3f        '%beta
        print '%-20s %s %10.1f %10.1f ' %(which, sbeta, band_ts, band_ts-ts),
        # beta is free: is it a good fit? check value, error if any
        if model.free[3]: # shouldn't be free
           model.free[3]=False
           model.internal_cov_matrix[3,:]=0
           model.internal_cov_matrix[:,3]=0
           print 'freezing E_break' ,
           refit=True
        if model.free[2] or beta>0.001:
            # free: is the fit ok?
            if beta>0.001 and beta_unc>0.001 and beta > 2*beta_unc:
                print ' fit is ok'
                if not model.free[2]:
                    model.free[2]=True # make sure free, since was once.
                    refit=True
                continue
            else:
                print '<--- reseting to PowerLaw' 
                model.free[2]=False
                # this should be done by the freeze? dangerous direct access, oh well.
                model.internal_cov_matrix[2,:]=0
                model.internal_cov_matrix[:,2]=0
                model[2]=0.
                models_to_fit.append(model)
                refit=True
                continue
        if beta>0 and beta<=0.001:
            print '<--- freezing at zero'
            model[2]=0
            model.internal_cov_matrix[2,:]=0
            model.internal_cov_matrix[:,2]=0
            continue

        if beta>=3.0: print 'beta>1 too large'; continue
        if beta==0: print 'frozen previously'; continue
        if band_ts<bts_min: print 'band_ts< %.1f'%bts_min; continue # not significant
        if band_ts-ts < qual_min and beta<=0.01:print 'qual<%.1f' %qual_min; continue # already a good fit
        print '<-- select to free beta' # ok, modify
        model.free[2]=True
        models_to_fit.append(model) # save for later check
        refit = True
    if refit:    
        print 'start refit with beta(s) freed or refixed...'
        roi.fit()
        roi.print_summary(title='after freeing one or more beta parameters')
        # now check for overflow
        refit = False
        for model in models_to_fit:
            beta = model[2]
            if beta < 3.0: continue
            print 'reseting model: beta =%.1f too large' % model[2]
            model.free[2]=False
            model[0:3]=(1e-15, 2.0, 3.0) #reset everything
            model.cov_matrix[:]=0 
            refit=True
        if refit:
            print 'need to re-refit: beta too large'
            roi.fit()
            roi.print_summary(title='re-refit after freeing/fixing one or more beta parameters')
    else:
        print 'none found'
    return refit    
        

def pickle_dump(roi, fit_sources, pickle_dir, dampen, failed=False, **kwargs):
    """ dump the source information from an ROI constructed from the sources here
    """
    
    name = roi.name.strip()
    fname = kwargs.pop('fname', name)
    filename=os.path.join(pickle_dir,fname+'.pickle')
    if os.path.exists(filename):
        # if the file exists
        output=pickle.load(open(filename))
        #curp = output.get('passes',None)
        #if curp is None: curp = []
        #curp.append(pass_number)
        #output['passes']=curp
        # update history of previous runs
        prev_logl = output.get('prev_logl', [])
        if prev_logl is None: prev_logl = []
        if dampen>0 : # update list only if a new fit
            last_logl = output.get('logl')
            prev_logl.append(last_logl)
            output['prev_logl'] = prev_logl
        oldsrc = output.get('sources', dict())
        print 'updating pickle file: log likelihood history:', \
             ''.join(map(lambda x: '%.1f, '%x, prev_logl)) 
    else:
        output = dict()
        oldsrc = dict()
        #output['passes']=[pass_number]
    diffuse_sources =  [s for s in roi.sources if s.skydir is None]\
                        +[s for s in roi.sources if isextended(s)]
    output['name'] = name
    output['skydir']  = roi.roi_dir
    # add extended sources to diffuse for backwards compatibilty: extended alos in the sources dict.
    output['diffuse'] = [s.spectral_model for s in diffuse_sources] #roi.dsm.models
    output['diffuse_names'] = [s.name for s in diffuse_sources]  #roi.dsm.names
    output['logl'] = -roi(roi.get_parameters()) if not failed else 0
    output['parameters'] = roi.get_parameters() 
    output['gradient'] = np.asarray(roi.gradient(), np.float32)
    output['time'] = time.asctime()

    if failed:
        
        f = open(filename,'wb') #perhaps overwrite
        pickle.dump(output,f)
        f.close()
        print 'saved (failed) pickle file to %s' % filename
        return
        
    sources=dict()
    output['sources']= sources
    def getit(s, key, savekey=None):
        """ key: current key
            savekey: key to save, expect to find in saved pickle
        """
        if savekey is None: savekey=key
        t = s.__dict__.get(key, None)
        if t is None and s.name in oldsrc:
            t = oldsrc[s.name].get(savekey, None)
            if t is not None:
                print 'ROI pickle: keeping previous calculation of %s for %s' % (savekey, s.name)
        return t
    def getit(s, key, savekey=None):
        """ key: current key
            savekey: key to save, expect to find in saved pickle
        """
        if savekey is None: savekey=key
        t = s.__dict__.get(key, None)
        if t is not None: return t #use a new version
        if s.name in oldsrc: #was it measured before?
            t = oldsrc[s.name].get(savekey, None)
            if t is not None: #yes, make a note and return it
                print 'ROI pickle: keeping previous calculation of %s for %s' % (key, s.name)
        return t
    
    for s in fit_sources:
        try:
            pivot_energy = s.spectral_model.pivot_energy()
        except: # if not fit
            pivot_energy = None 
        if pivot_energy is None: pivot_energy=s.spectral_model.e0 # could be a pulsar?
        
        if 'sedrec' not in s.__dict__ or s.sedrec is None:
            print 'warning: no sedrec in %s' %s.name
            s.sedrec = None
        sedrec = s.sedrec
        loc = s.__dict__.get('loc', None)
        sources[s.name]=dict(
            skydir=s.skydir, 
            model=s.spectral_model,
            isextended=isextended(s),
            #extent = s.__dict__.get('spatial_model', None),
            ts = s.ts if hasattr(s, 'ts') else roi.TS(s.name),
            sedrec = sedrec,
            band_ts=0 if sedrec is None else sedrec.ts.sum(),
            pivot_energy = pivot_energy,
            # if ellipse or adict not done, but already in pickle, keep them
            ellipse= getit(s, 'ellipse'), #s.__dict__.get('ellipse', None), 
            associations = getit(s, 'adict', 'associations'), #s.__dict__.get('adict',None),
            )
    output.update(kwargs) # add additional entries from kwargs
    f = open(filename,'wb') #perhaps overwrite
    pickle.dump(output,f)
    f.close()
    print 'saved pickle file to %s' % filename
        

def repivot(roi, fit_sources=None, min_ts = 10, max_beta=3.0, emin=200, emax=20000.):
    """ invoked by process() if repivot flag set; can be run separately to test
    
    returns True if had to refit, allowing iteration
    """
    
    print '\ncheck need to repivot sources with TS>%.0f, beta<%.1f: \n'\
    'source                     TS        e0      pivot' % (min_ts, max_beta)
    need_refit =False
    if fit_sources is None:
        fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    print 'processing %d sources' % len(fit_sources)
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
        if abs(pivot/e0-1.)<0.05:
            print 'converged'; continue
        print 'will refit'
        need_refit=True
        model.set_e0(pivot)
    if need_refit:
        roi.fit()
    return need_refit
        #roi.print_summary(sdir=roi.roi_dir,maxdist=5, title='after pivot refit(s): logL=%0.f' % roi.logl)

class Damper(object):
    """ manage adjustment of parameters for damping """
    def __init__(self, roi, dampen, tol=2.0):
        self.dampen=dampen
        self.roi = roi
        self.ipar = roi.get_parameters()[:] # get copy of initial parameters
        self.initial_logl = -roi(self.ipar)
        self.tol = tol
        self.change = 0
    def __call__(self):
        self.roi.logl = -self.roi(self.roi.get_parameters())
        self.change = abs(self.initial_logl-self.roi.logl)
        if self.change< self.tol or self.dampen==1.0 : return False
        # change more than tol and dampen specified: dampen
        fpar = self.roi.get_parameters()
        if len(fpar)==len(self.ipar): 
            dpar = self.ipar+self.dampen*(fpar-self.ipar)
            self.roi.set_parameters(dpar)
            self.roi.logl=-self.roi(dpar)
            # check for change, at least 0.5
        return True
    

def process(roi, **kwargs):
    """ process the roi object after being set up
    """
    outdir   = kwargs.get('outdir', None)
    localize = kwargs.pop('localize', False)
    repivot_flag  = kwargs.pop('repivot', False)
    fixbeta  = kwargs.pop('fix_beta', False)
    dampen   = kwargs.pop('dampen', 1.0)  # factor to adjust parameters before writing out: if zero, no fit
    counts_dir = kwargs.pop('counts_dir', 'counts_dir')
    pass_number= 0
    tables = kwargs.pop('tables', None)
    localize_kw = kwargs.pop('localize_kw', {}) # could have bandfits=False
    diffuse_only = kwargs.pop('diffuse_only', False)
    norms_first = kwargs.pop('norms_first', True)
    freeze_iem = kwargs.pop('freeze_iem', 1.0)
    freeze_iso = kwargs.pop('freeze_iso', None)
    countsplot_tsmin = kwargs.pop('countsplot_tsmin', 100) # minimum for counts plot
    source_name = kwargs.pop('source_name', None) # for localize perhaps
    damp = Damper(roi, dampen)
    
    if outdir is not None: 
        counts_dir = os.path.join(outdir,counts_dir)
        logpath = os.path.join(outdir, 'log')
        if not os.path.exists(logpath): os.mkdir(logpath)
        outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
        print  '='*80
        print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    else: outtee=None

    if counts_dir is not None and not os.path.exists(counts_dir) :
        try: os.makedirs(counts_dir) # in case some other process makes it
        except: pass
    if freeze_iem is not None:
        print 'Freezeing IEM to %f' % freeze_iem
        roi.freeze('Norm', 'ring', freeze_iem)
    if freeze_iso is not None:
        print 'Freezeing isotropic to %f' % freeze_iso
        roi.freeze('Scale', 'iso*', freeze_iso)

    init_log_like = roi.log_like()
    roi.print_summary(title='before fit, logL=%0.f'% init_log_like)
    fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    if len(roi.get_parameters())==0:
        print '===================== nothing to fit========================'
    else:
        if dampen>0:
            fit_kw = kwargs.get('fit_kw', {})
            try:
                if norms_first:
                    t = np.array([n.endswith('Norm') for n in roi.parameter_names])
                    if sum(t)>0:
                        print 'Fitting parameter names ending in "Norm"'
                        roi.fit(np.arange(len(t))[t], **fit_kw)

                if diffuse_only:
                    ndiff = len([n for n in roi.parameter_names if n.split('_')[0] in ('ring','isotrop', 'limb')])
                    roi.summary(range(ndiff), title='Before fit to diffuse components')
                    fit_kw.update(select=range(ndiff))
                    roi.fit(**fit_kw)
                else:
                    roi.fit(**fit_kw)
                    change =roi.log_like() - init_log_like 
                    if  abs(change)>1.0 :
                        roi.print_summary(title='after global fit, logL=%0.f, change=%.1f'% (roi.log_like(), change))
                if repivot_flag:
                    # repivot, iterating a few times
                    n = 3
                    while n>0:
                        if not repivot(roi, fit_sources): break
                        n-=1
                if fixbeta:
                    if not fix_beta(roi):
                        print 'fixbeta requested, but no refit needed, quitting'
                if damp():
                    change =roi.log_like() - init_log_like
                    roi.print_summary(title='after damping with factor=%.2f, logL=%0.f, change=%.1f'\
                        %(dampen,roi.log_like(), change))
                else:
                    print 'No damping requested (dampen=%.1f), or difference in log Likelihood (%.1f) < tol (%.1f)'\
                           % (dampen, damp.change,damp.tol )
            except Exception, msg:
                print '============== fit failed, aborting!! %s'%msg
                #pickle_dump(roi, fit_sources, os.path.join(outdir, 'pickle'), pass_number, failed=True)
                return False
        else:
            print 'No fit requested'
    
    outdir     = kwargs.pop('outdir', '.')
    associator=  kwargs.pop('associate', None) #None #### disable for now ####
    if type(associator)==types.StringType and associator=='None': associator=None
    def getdir(x ):
        if not kwargs.get(x,False): return None
        t = os.path.join(outdir, kwargs.get(x))
        if not os.path.exists(t): os.mkdir(t)
        return t
    skymodel_name = os.path.split(os.getcwd())[-1]
    sedfuns.makesed_all(roi, sedfig_dir=getdir('sedfig_dir'), suffix='_sed_%s'%skymodel_name, )
    if localize:
        print 'localizing and associating all sources with variable...'
        q, roi.quiet = roi.quiet,False
        tsmap_dir = None ####### getdir('tsmap_dir') 
        localization.localize_all(roi, tsmap_dir=tsmap_dir, associator = associator, source_name=source_name)
        roi.quiet=q

    try:
        ax = counts.stacked_plots(roi,None, tsmin=countsplot_tsmin)
        cts=counts.get_counts(roi)
        chisq = cts['chisq']
        print 'chisquared for counts plot: %.1f'% chisq
 
        if counts_dir is not None and dampen>0:
            fout = os.path.join(counts_dir, ('%s_counts.png'%roi.name) )
            ax[1].figure.savefig(fout, dpi=60)
            print 'saved counts plot to %s' % fout
    except Exception,e:
        print 'Failed to analyze counts for roi %s: %s' %(roi.name,e)
        chisq = -1
    
    if tables is not None:
        tables(roi)
    
    if outdir is not None:  
        pickle_dir = os.path.join(outdir, 'pickle')
        if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
        pickle_dump(roi, fit_sources, pickle_dir, dampen, 
            initial_logl=damp.initial_logl, 
            counts=cts,
            )
    if outtee is not None: outtee.close() 
    return chisq

def localize(roi, **kwargs):
    """ special processor to perform localization for tests 
    also do tsmaps for failed fits if kw 'poor_loc' is name of csv file
    (sort of wired in now)
    """
    outdir = kwargs.get('outdir', '.')
    tsmin = kwargs.get('tsmin', 10)
    logpath = os.path.join(outdir, 'log')
    locdir = kwargs.pop('locdir', 'localization')
    tsmap_dir = 'tsmap_fail' #kwargs.get('tsmap_dir', None)
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    poor_loc = kwargs.pop('poor_loc', 'poorly_localized.csv')
    if poor_loc is not None:
       poor = pd.read_csv(poor_loc, index_col=0)
       names = list(poor[poor.roiname==roi.name].index.values)
       if len(names)==0:
           print '***no source names to fit'
           outtee.close()
           return
       print 'Localizing and perhaps ts_maps for: ', names
       for source_name in names:
           localization.localize_all(roi, source_name=source_name, tsmap_dir=tsmap_dir)
       outtee.close()
       return
       
    emin = kwargs.pop('emin', None)
    if emin is not None:
        roi.select_bands(emin=emin)
    # extract only relevant keywords
    loc_kw =dict() 
    for x in localization.Localization.defaults:
        kw = x[0]
        if kw in kwargs:
            loc_kw[kw]=kwargs[kw]
    loc_kw.update(tsmap_dir=kwargs.get('tsmap_dir',None), source_name=kwargs.get('source_name',None))
    localization.localize_all(roi,  **loc_kw)
    if locdir is None:
       outtee.close()
       return
    if not os.path.exists(locdir): os.mkdir(locdir)
    sources = [s for s in roi.sources if s.skydir is not None\
        and s.__dict__.get(  'spatial_model', None) is None \
        and np.any(s.spectral_model.free) and s.ts>tsmin]

    for source in sources:
        outfile = os.path.join(locdir, source.name.replace(' ', '_').replace('+','p')+'.pickle')
        pickle.dump(source, open(outfile, 'w'))
        print 'wrote file %s' % outfile
    outtee.close()

    
def table_processor(roi, **kwargs):
    """ do only tables """
    outdir = kwargs.get('outdir')
    tables = kwargs.get('tables')
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    tables(roi)
    outtee.close()
    
def residual_processor(roi, **kwargs):
    """ a roi processor that creates tables of residuals from the bands """
    from uw.like2 import roistat, bandplot
    from . import maps
    outdir = kwargs['outdir']
    iband_range = kwargs.get('iband_range', range(4))
    print 'generating residual plots for %s, bands %s' % (roi.name, iband_range)
    class ResidualCounts(object):
        def __init__(self, roi,  **kwargs):
            iband = kwargs.get('iband')
            band = roi.bands[iband]
            energy, event_class = band.e, band.b.event_class()
            s = roistat.ROIstat(roi, lambda b: b.e==energy and b.b.event_class()==event_class)
            bm = s.all_bands[0]
            resids = (bm.data-bm.spectral_model)/np.sqrt(bm.spectral_model)
            self.skyfun = bandplot.BandSkyFunction(band, resids)
        def __call__(self,v):
            skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
            return self.skyfun(skydir)
    table = maps.ROItables(outdir, skyfuns= [
        [ResidualCounts, "residual_%02d"%iband, dict(iband=iband)] for iband in iband_range])
    table(roi)    

def background_density(roi, source_name):
    bgsources = [(j,s) for j,s in enumerate(roi.sources) if s.skydir is None]
    bgdensity = []
    b0 = roi.selected_bands[0]
    pti = list(roi.sources.source_names).index(source_name)
    try:
        for j, bgsource in bgsources:
            bgdensity.append(np.array([ 
                ( band[j].pix_counts * band[pti].pix_counts ).sum() / band[pti].counts\
                            for band in roi.selected_bands],np.float32) )
    except Exception, e:
        print 'failed bgdensity for source %s: %s' % (source.name, e)
    return bgdensity

def full_sed_processor(roi, **kwargs):
    """ roi processor to include front and back sed info 
    """
    print 'processing ROI %s: creating full sedinfo ' %roi.name
    outdir   = kwargs.get('outdir')
    ts_min   = kwargs.get('ts_min', 10)
    sedinfo = os.path.join(outdir, 'sedinfo')
    make_seds = kwargs.get('make_seds', True)
    if not os.path.exists(sedinfo): os.mkdir(sedinfo)
    ptsources = [(i,s) for i,s in enumerate(roi.sources) if s.skydir is not None and np.any(s.spectral_model.free)]
    bgsources = [(j,s) for j,s in enumerate(roi.sources) if s.skydir is None]
    print 'evaluating %d sources' %len(ptsources)
    def getdir(x ):
        if not kwargs.get(x,False): return None
        t = os.path.join(outdir, kwargs.get(x))
        if not os.path.exists(t): os.mkdir(t)
        return t
    if make_seds:
       sedfuns.makesed_all(roi, sedfig_dir=getdir('sedfig_dir'))
    for pti,source in ptsources:
        print 'processing %s, ts=%.1f' % (source.name, source.ts)
        try:
            ts = source.ts #roi.TS(source.name)
            if ts<ts_min: continue
            sf = sedfuns.SourceFlux(roi, source.name, )
            sedrecs = [sedfuns.SED(sf, event_class).rec for event_class in (0,1,None)]
                     
        except Exception,e:
            print '*** source %s failed flux measurement: %s' % (source.name, e)
            raise
        fname = source.name.replace('+','p').replace(' ', '_')+'_sedinfo.pickle'
        fullfname = os.path.join(sedinfo, fname)
        bgdensity = []
        roi.update() # make sure likelihood initialize
        try:
            for j, bgsource in bgsources:
                bgdensity.append(np.array([ 
                    ( band[j].pix_counts * band[pti].pix_counts ).sum() / band[pti].counts\
                                for band in roi.selected_bands],np.float32) )
        except Exception, e:
            print 'failed bgdensity for source %s: %s' % (source.name, e)
        
        bgcounts = []
        for j, bgsource in bgsources:
            bgcounts.append(np.array(
                [band[j].counts for band in roi.selected_bands], np.float32))
        try:
            counts = np.array([band[source.name].counts for band in roi.selected_bands], np.float32)
        except:
            counts = None 
        sdir = source.skydir
        pickle.dump( dict(
                ra=sdir.ra(), dec=sdir.dec(), 
                glat=sdir.b(), glon=sdir.l(),
                ts=ts, 
                elow =sedrecs[0].elow,
                ehigh=sedrecs[0].ehigh,
                flux = np.vstack([sr.flux   for sr in sedrecs]),
                lflux= np.vstack([sr.lflux  for sr in sedrecs]),
                uflux= np.vstack([sr.uflux  for sr in sedrecs]),
                mflux= np.vstack([sr.mflux  for sr in sedrecs]),
                bts   =np.vstack([sr.ts     for sr in sedrecs]),
                counts = counts,
                bgdensity = bgdensity,
                bgcounts = bgcounts,
                bgnames = [s.name for ii,s in bgsources],
                ),
            open(fullfname, 'w'))
        print 'wrote sedinfo pickle file to %s' % fullfname
    sys.stdout.flush()
            

def covariance(roi, **kwargs):
    """Compute covariance for first two backgrounds, low energy bands variable point sources
    """
    outdir   = kwargs.get('outdir', '.')
    ts_min   = kwargs.get('ts_min', 25)
    covinfo = os.path.join(outdir, 'covariance')
    if not os.path.exists(covinfo): os.mkdir(covinfo)
    ptsources = [(i,s) for i,s in enumerate(roi.sources) if s.skydir is not None and np.any(s.spectral_model.free)]
    print 'evaluating %d sources' %len(ptsources)
    
    for pti,source in ptsources:
        try: 
            ts = source.ts 
        except: 
            ts = roi.TS(source.name)
        if ts<ts_min: continue
        fname = source.name.replace('+','p').replace(' ', '_')+'_covinfo.pickle'
        fullfname = os.path.join(covinfo, fname)

        def cov3(bd):
            """ for source at index j, in bandlike bd, return errors, correlation coefficients with
            respect to first two diffuse sources
            """
            ix = (pti,0,1)
            cix = [bd[i].pix_counts for i in ix] #counts
            q = bd.data/(bd.model_pixels)**2
            t = [[sum(x*y*q) for x in cix] for y in cix] 
            v = np.matrix(t).I
            sigs = np.array(np.sqrt(v.diagonal()),np.float32)[0]
            return dict(sigs=sigs,
                        cc = np.array([v[i,j]/(sigs[i]*sigs[j]) for i,j in ((0,1),(0,2),(1,2))], np.float32) ,
                        unweight = bl.unweight,)

        covariance = []
        try:
            covariance.append(np.array([cov3(bl) for bl in roi.selected_bands[:16]]))
                
        except Exception, msg:
            print '*** failed covariance for source %s: %s' % (source.name,msg)
        pickle.dump( dict(name=source.name,
                        covariance=covariance, 
                        ts=ts,
                        skydir=source.skydir ,), 
                open(fullfname, 'w'))
        print 'wrote covariance pickle file to %s' % fullfname
        
    sys.stdout.flush()



def roi_refit_processor(roi, **kwargs):
    """ roi processor to refit diffuse for each energy band
    
    kwargs
    ------
    outdir : string
    emax : float 
        maximum energy to fit (default 10000)
    diffuse_names : list of diffuse names to refit
        default ('ring',)
    plot_dir : string
        sub folder in outdir to save plots ['galfits_plots']
    fit_dir : string
        sub folder in outdir to save pickled fits ['galfits_all']
    """
    outdir= kwargs.get('outdir', '.')
    emax =  kwargs.pop('emax', 10000)
    names = kwargs.get('diffuse_names', ('ring',))
    plot_dir =os.path.join(outdir, kwargs.pop('plot_dir', 'galfit_plots'))
    fit_dir = os.path.join(outdir, kwargs.pop('fit_dir', 'galfits_all'))
    ylim = kwargs.get('ylim', (0.75,1.25))
    if not os.path.exists(fit_dir): os.mkdir(fit_dir)
    if not os.path.exists(plot_dir): os.mkdir(plot_dir)
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))


    print 'processing ROI %s: refitting diffuse %s' % (roi.name, names)

    def make_plot( d, labels=('front','back','both'), outdir=None, ylim=(0.75, 1.25)):
        plt.figure(99, figsize=(5,5));plt.clf()
        data = [np.array(d[label]['values'])[:,0] for label in labels]
        errors= [np.array(d[label]['errors'])[:,0] for label in labels]
        energy = np.array(d['both']['energies'])
        name = d['roiname']
        diffuse_name = d['diffuse_names'][0] #assume only one?
        for y,yerr,label in zip(data,errors,labels):
            plt.errorbar(energy[:8]/1e3, y[:8], yerr=yerr[:8], 
                fmt='o-' , ms=10, lw=2, label=label)
        plt.xscale('log'); plt.ylim(ylim)
        plt.title('Diffuse fits to %s for ROI %s'%(diffuse_name, name), fontsize='medium')
        plt.axhline(d['norms'][0], color='k', linestyle='--', label='full fit')
        plt.xlabel('Energy (GeV)');plt.ylabel('normalization factor')
        plt.grid(True);plt.legend(loc='upper right');
        if outdir is not None:
            filename = os.path.join(outdir,  '%s_%s.png'%(name, diffuse_name))
            plt.savefig(filename, dpi=60)
            print 'saved figure to %s'%filename
        
    dlike = sedfuns.DiffuseLikelihood(roi, names=names)
    r = dict(roiname=roi.name, glon=roi.roi_dir.l(), glat=roi.roi_dir.b(), 
        diffuse_names=names, norms=dlike.saved_pars)
    for event_class,name in ((0,'front'),(1,'back'),(None,'both')):
        r[name] = dlike.multifit(event_class=event_class, emax=emax)
    make_plot(r, outdir=plot_dir, ylim = ylim)    
    fullfname = os.path.join(fit_dir, roi.name+'_diffuse.pickle')
    pickle.dump( r, open(fullfname, 'w'))
    print 'wrote diffuse refit pickle file to %s' % fullfname
    outtee.close()
    
def iso_refit_processor(roi, **kwargs):
    kwargs.update(diffuse_names=('isotrop',),
            plot_dir='isofit_plots', fit_dir='isofits', ylim=(0,2))
    return roi_refit_processor(roi, **kwargs)

def limb_processor(roi, **kwargs):
    """ report on limb fit, perhaps refit"""
    outdir= kwargs.get('outdir')
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s limb processor' %(time.localtime()[:6]+ (roi.name,))

    limbdir = os.path.join(outdir, kwargs.get('limbdir', 'limb'))
    if not os.path.exists(limbdir): os.mkdir(limbdir)
    refit = kwargs.get('refit', True)
    try:
        limb = roi.get_model('limb')
    except:
        print 'No limb source: no pickle file saved'
        outtee.close()
        return
    if refit:
        for m in limb.models:
            m.free[:]=True
        roi.initialize()
        names = roi.parameter_names
        u = np.arange(len(names))
        i = np.array(['limb' in x for x in names])
        
        roi.fit(u[i])
        
    t = roi.all_bands[1][2]
    fname = os.path.join(limbdir,'%s_limb.pickle' %roi.name)
    pickle.dump( dict(ra=roi.roi_dir.ra(), dec=roi.roi_dir.dec(), model=limb,
        counts=t.counts, pixel_values=t.pixel_values), open(fname,'w'))
    print 'wrote file to %s' %fname
    outtee.close()


def sunmoon_processor(roi, **kwargs):
    """ report on sunmon refit """
    outdir= kwargs.get('outdir')
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s sunmoon processor' %(time.localtime()[:6]+ (roi.name,))

    sunmoondir = os.path.join(outdir, kwargs.get('sunmoondir', 'sunmoon'))
    if not os.path.exists(sunmoondir): os.mkdir(sunmoondir)
    refit = kwargs.get('refit', True)
    try:
        sunmoon = roi.get_model('SunMoon')
    except:
        print 'No sunmoon source: no pickle file saved'
        outtee.close()
        return
    sunmoon.free[:] = True
    roi.initialize()
    names = roi.parameter_names
    u = np.arange(len(names))
    i = np.array(['SunMoon' in x for x in names])
    before = roi.log_like()
    roi.fit(u[i])
        
    fname = os.path.join(sunmoondir,'%s_sunmoon.pickle' %roi.name)
    pickle.dump( dict(ra=roi.roi_dir.ra(), dec=roi.roi_dir.dec(), model=sunmoon, skydir=roi.roi_dir,
      delta_likelihood=roi.log_like()-before ), open(fname,'w'))
    print 'wrote file to %s' %fname
    outtee.close()

    
def check_seeds(roi, **kwargs):
    """ Evaluate a set of seeds: fit, localize with position update, fit again
    """
    outdir = kwargs.get('outdir')
    prefix = kwargs.get('prefix')
    tsmap_dir=kwargs.get('tsmap', None)
    tsmin = kwargs.pop('tsmin', 10)
    repivot_flag = True ##########kwargs.pop('repivot', True)
    seedcheck_dir = kwargs.get('seedcheck_dir', 'seedcheck')
    if not os.path.exists(seedcheck_dir): os.mkdir(seedcheck_dir)
    associator= kwargs.pop('associate', None)
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    print 'Checking seeds, if any'
    fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    seed_sources = [ s for s in fit_sources if s.name.startswith(prefix)]
    if len(seed_sources)==0: 
        print 'no sources to refit, locate'
        outtee.close()
        return
    # initial fit to norm only
    seednorms = np.arange(len(roi.parameter_names))[np.array([s.startswith(prefix) and s.endswith('_Norm') for s in roi.parameter_names])]
    roi.fit(seednorms)
    for s in seed_sources:
        s.ts=ts = roi.TS(s.name)
    localization.localize_all(roi, prefix=prefix, tsmap_dir=tsmap_dir, associator = associator, update=True, tsmin=tsmin)
    roi.fit() 
    if repivot_flag: 
        repivot(roi) # one iteration
    for s in seed_sources: 
        s.ts=ts = roi.TS(s.name)
        sfile = os.path.join(seedcheck_dir, s.name.replace(' ', '_').replace('+','p')+'.pickle')
        pickle.dump(s, open(sfile, 'w'))
        print 'wrote file %s' %sfile
    outtee.close()

class GtlikeCatalog(object):
    """ read in a gll catalog FITS file, then make spectral model available
    """
    def __init__(self, name=None): #'gll_psc4yearclean_v4.fit'):
        import pyfits
        catfile = sorted(glob.glob(os.path.expandvars('$FERMI/catalog/gll*.fit')))[-1]
        print 'opening catalog file %s' % catfile
        self.cat=pyfits.open(catfile)[1].data
        
    def __call__(self, nickname):
        """ Return a Model corresponding to the nickname field
        """
        from uw.like import Models
        cselect = np.array([np.any(s.field('NickName')==(nickname.replace(' ',''))) for s in self.cat])
        if sum(cselect)!=1:
            print 'did not find a source with NickName %s' %nickname
            return None
        catsource = self.cat[cselect][0]
        try:
            st = catsource.field('SpectrumType')
            flux,index,cutoff,b,pivot,beta=[catsource.field(f) for f in 'Flux_Density Spectral_Index Cutoff Index2 Pivot_Energy beta'.split()]
            if st=='PowerLaw':
                return Models.PowerLaw(p=[flux, index], e0=pivot)
            elif st=='PLSuperExpCutoff':
                prefactor = flux*np.exp((pivot/cutoff)**b)
                return Models.PLSuperExpCutoff(p=[prefactor, index, cutoff,b], e0=pivot)
            elif st=='LogParabola':
                return Models.LogParabola(p=[flux, index, beta, pivot])
            elif st=='PowerLaw2':   ### same as PowerLaw in table
                return Models.PowerLaw(p=[flux, index], e0=pivot)
            else:
                raise Exception('unexpected spectrum type %s'%st)
                catmodel = cat_model(catsource); 
        except:
            print 'Source %s failed' % nickname
            return None


class CompareOtherModel(object):

    def __init__(self, roi, other='uw22c'):
        self.others = pickle.load(open('../%s/sources.pickle' % other))
        self.other_models=[]
        self.roi=roi
    
    def __call__(self, source):
        roi = self.roi
        source_name = source.name
        ptts = roi.TS(source_name)
        #    sed = roi.plot_sed(source_name, butterfly=False, annotate=None, fit_kwargs=dict(label='this: %.0f'%ptts, color='orange', lw=2))
        plot_kw=dict(energy_flux_unit=kwargs.pop('energy_flux_unit','eV'),
                 gev_scale=kwargs.pop('gev_scale',True))
        roi.get_sed(update=True)
        ps = sed.Plot(source, **plot_kw)
        #annotation =(0.05,0.9, 'TS=%.0f'% self.TS(source.name))
        plot_kw = dict(label= 'current: %.0f'%ptts )#annotate=annotation)
        ps(fit_kwargs=plot_kw)
        saved_model = source.spectral_model
        axes = plt.gca()
        if source_name in self.others.index:
            othermodel = self.others.ix[source_name]['model']
            source.spectral_model = othermodel
            gtts = roi.TS(source_name)
            roi.get_sed(update=True)
            ps.plot_model( othermodel, butterfly=False, label='other: %.0f'%gtts, color='g', lw=2)
            source.spectral_model = saved_model
        else:
            gtts=-1
            othermodel=None
        axes.legend(prop=dict(size=10))
        plt.setp(axes, xlim=(100, 31.6e3), ylim=(0.1,1000))
        outfile = os.path.join(sed_dir, '%s_sed.png' % (source_name.replace(' ','_').replace('+','p')))
        plt.savefig(outfile)
        print 'wrote file %s' % outfile
        self.other_models.append(dict(name=source_name, m_other=othermodel, m_pt=saved_model, ts_pt=ptts, ts_other=gtts))
 
class GtlikeModels(object):
    def __init__(self, catpath=os.path.expandvars('$FERMI/catalog/gll_psc4year*.fit'), otherversion='v7'):
        import pyfits, glob
        catfile = sorted(glob.glob(catpath))[-1]
        print 'Loading gtlike file %s' % catfile
        self.cat = pyfits.open(catfile)[1].data
        self.otherversion=otherversion
        self.uwversion = os.path.split(os.getcwd())[-1]
                 
    def lookup(self, source):
        """ return the corresponding model, or None"""
        from uw.like import Models
        cselect = np.array([np.any(s.field('NickName')==(source.name.replace(' ',''))) for s in self.cat])
        if sum(cselect)!=1: return None
        cs = self.cat[cselect][0]
        st = cs.field('SpectrumType')
        flux,index,cutoff,b,pivot,beta=[cs.field(f) for f in 'Flux_Density Spectral_Index Cutoff Exp_Index Pivot_Energy beta'.split()]
        try:
            if st=='PowerLaw':
               return Models.PowerLaw(p=[flux, index], e0=pivot)
            elif st=='PLSuperExpCutoff':
                prefactor = flux*np.exp((pivot/cutoff)**b)
                return Models.PLSuperExpCutoff(p=[prefactor, index, cutoff,b], e0=pivot)
            elif st=='LogParabola':
                return Models.LogParabola(p=[flux, index, beta, pivot])
            elif st=='PowerLaw2':   ### same as PowerLaw in table
                return Models.PowerLaw(p=[flux, index], e0=pivot)
            else:
                raise Exception('unexpected spectrum type %s'%st)
        except Exception, msg:
            print 'Failed to create model for %s: "%s"' % (source.name, msg)
            return None
                
    def sed_plot(self, roi, source, **kwargs):
        source_name = source.name
        othermodel = self.lookup(source)
        if othermodel is None: 
            print 'Source %s not found' %source_name
            return None
        ptts = roi.TS(source_name)
        plot_kw=dict(energy_flux_unit=kwargs.pop('energy_flux_unit','eV'),
                 gev_scale=kwargs.pop('gev_scale',True))
        roi.get_sed(update=True)
        ps = sed.Plot(source, **plot_kw)
        #annotation =(0.05,0.9, 'TS=%.0f'% self.TS(source.name))
        plot_kw = dict(label= '%s: %.0f'%(self.uwversion,ptts) )#annotate=annotation)
        ps(fit_kwargs=plot_kw)
        saved_model = source.spectral_model
        axes = plt.gca()
        source.spectral_model = othermodel
        gtts = roi.TS(source_name)
        roi.get_sed(update=True)
        ps.plot_model( othermodel, butterfly=False, label='%s: %.0f'%(self.otherversion,gtts), color='g', lw=2)
        source.spectral_model = saved_model
        axes.legend(prop=dict(size=10))
        plt.setp(axes, xlim=(100, 31.6e3), ylim=(0.1,1000))
        outfile = os.path.join(self.sed_dir, '%s_sed.png' % (source_name.replace(' ','_').replace('+','p')))
        plt.savefig(outfile)
        print 'wrote file %s' % outfile
        return dict(name=source.name, uwts=ptts, other_ts=gtts, othermodel=othermodel)
   
    def __call__(self, roi, **kwargs):
        outdir = kwargs.pop('outdir', '.')
        self.sed_dir = os.path.join(outdir, kwargs.pop('sed_dir', 'gtlike/sed'))
        model_dir=os.path.join(outdir, kwargs.pop('model_dir', 'gtlike/models'))
        for t in (self.sed_dir, model_dir):
            if not os.path.exists(t): os.makedirs(t)
        sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
        catmodels=[]
        for source in sources:
            catmodel = self.sed_plot(roi, source)
            if catmodel is not None:
                catmodels.append(  catmodel )
        outfile = os.path.join(model_dir, '%s.pickle'%roi.name)
        pickle.dump(catmodels,open(outfile, 'w'))
        print 'wrote file %s' %outfile
        

    
gtm = None
def gtlike_compare(roi, **kwargs):
    """
    compare with gtlike
    """
    global gtm
    if gtm is None:
        gtm = GtlikeModels()
    gtm(roi, **kwargs)
    
    
others=None
def UW_compare(roi, **kwargs):
    global others
    other = kwargs.pop('other', 'uw25') #wire in default for now
    if others is None:
        others = pickle.load(open('../%s/sources.pickle' % other))
    outdir = kwargs.pop('outdir', '.')
    sed_dir = os.path.join(outdir, kwargs.pop('sed_dir', 'compare_%s/sed'%other))
    model_dir=os.path.join(outdir, kwargs.pop('model_dir', 'compare_%s/models' %other))
    for t in (sed_dir,model_dir):
        if not os.path.exists(t): os.makedirs(t)
    other_models=[]
    def check_othersource(source):
        source_name = source.name
        ptts = roi.TS(source_name)
        #    sed = roi.plot_sed(source_name, butterfly=False, annotate=None, fit_kwargs=dict(label='this: %.0f'%ptts, color='orange', lw=2))
        plot_kw=dict(energy_flux_unit=kwargs.pop('energy_flux_unit','eV'),
                 gev_scale=kwargs.pop('gev_scale',True))
        roi.get_sed(update=True)
        ps = sed.Plot(source, **plot_kw)
        #annotation =(0.05,0.9, 'TS=%.0f'% self.TS(source.name))
        plot_kw = dict(label= 'current: %.0f'%ptts )#annotate=annotation)
        ps(fit_kwargs=plot_kw)
        saved_model = source.spectral_model
        axes = plt.gca()
        if source_name in others.index:
            othermodel = others.ix[source_name]['model']
            source.spectral_model = othermodel
            gtts = roi.TS(source_name)
            roi.get_sed(update=True)
            ps.plot_model( othermodel, butterfly=False, label='other: %.0f'%gtts, color='g', lw=2)
            source.spectral_model = saved_model
        else:
            gtts=-1
            othermodel=None
        axes.legend(prop=dict(size=10))
        plt.setp(axes, xlim=(100, 31.6e3), ylim=(0.1,1000))
        outfile = os.path.join(sed_dir, '%s_sed.png' % (source_name.replace(' ','_').replace('+','p')))
        plt.savefig(outfile)
        print 'wrote file %s' % outfile
        other_models.append(dict(name=source_name, m_other=othermodel, m_pt=saved_model, ts_pt=ptts, ts_other=gtts))

    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    map( check_othersource, sources)
    outfile = os.path.join(model_dir, '%s.pickle'%roi.name)
    pickle.dump(other_models,open(outfile, 'w'))
    print 'wrote file %s' %outfile
 
    
    
def flux_correlations(roi, **kwargs):

    class DiffuseDependence(object):
        
        def __init__(self, roi, diffuse='ring'):
            self.roi = roi        
            normpar = np.array([name.endswith('_Norm') for name in roi.parameter_names])
            normpar[0]=False
            self.names = np.array([name.split('_')[0] for name in roi.parameter_names])[normpar]
            self.ipar = map(int, np.arange(len(normpar))[normpar])
            self.ring = roi.get_model(diffuse)
            self.rnorm = self.ring.getp(0)
            self.df = pd.DataFrame(index=self.names)
            
        def pars(self):
            return self.roi.model_parameters[self.ipar]
        def errs(self):
            return self.roi.sources.uncertainties[self.ipar]
        
        def ts(self):
            return map(self.roi.TS, self.names)
        
        def fit(self, suffix, delta):
            self.ring.setp(0, self.rnorm*(1.+delta))
            t = self.roi.fit(self.ipar)
            self.df['par_'+suffix]=t.get_parameters()
            self.df['unc_'+suffix]=self.errs()
            self.df['ts_'+suffix]=self.ts()
            roi.set_parameters(roi.saved_pars) # always restore
                           
        def run(self, emin=100, delta=0.01,):
            """ 
            """
            if len(self.names)==0:
                return 
            if emin >100:
                roi.select_bands(emin=emin)
            for suffix, delt in zip('z p m'.split(), (0,delta,-delta)):
                self.fit(suffix+'%d'%emin, delt)
    outdir= kwargs.get('outdir')
    eminlist = kwargs.get('eminlist', (100, 316, 1000) )
    flux_corr_dir = os.path.join(outdir, kwargs.get('fluxcorr', 'fluxcorr'))
    if not os.path.exists(flux_corr_dir): os.mkdir(flux_corr_dir)
    
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    diffuse=kwargs.get('diffuse', 'ring')
    print 'Running diffuse dependence for %s'%diffuse
    t = DiffuseDependence(roi, diffuse=diffuse)
    for emin in eminlist:
        t.run(emin)
    
    fname = os.path.join(flux_corr_dir,'%s_fluxcorr.pickle' %roi.name)
    pickle.dump( t.df, open(fname,'w'))
    print 'wrote file to %s' %fname
    outtee.close()
  

