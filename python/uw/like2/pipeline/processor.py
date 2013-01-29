"""
roi and source processing used by the roi pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/processor.py,v 1.40 2013/01/29 18:40:36 burnett Exp $
"""
import os, time, sys, types
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
        sbeta = '%5.2f+/-%5.2f' % (beta, beta_unc) if beta_unc>0 else '%5.2f        '%beta
        print '%-20s %s %10.1f %10.1f ' %(which, sbeta, band_ts, band_ts-ts),
        # beta is free: is it a good fit? check value, error if any
        if model.free[2]:
            # free: is the fit ok?
            if beta>0.01 and beta_unc>0.001 and beta_unc< 2*beta:
                print ' fit is ok'
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
    output['logl'] = roi.logl if not failed else 0
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
    def getit(s, key):
        t = s.__dict__.get(key, None)
        if t is None and s.name in oldsrc:
            t = oldsrc[s.name].get(key, None)
            if t is not None:
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
            associations = getit(s, 'adict'), #s.__dict__.get('adict',None),
            )
    output.update(kwargs) # add additional entries from kwargs
    f = open(filename,'wb') #perhaps overwrite
    pickle.dump(output,f)
    f.close()
    print 'saved pickle file to %s' % filename
        

def repivot(roi, fit_sources=None, min_ts = 10, max_beta=3.0, emin=200, emax=10000.):
    """ invoked by process() if repivot flag set; can be run separately to test
    
    returns True if had to refit, allowing iteration
    """
    
    print '\ncheck need to repivot sources with TS>%.0f, beta<%.1f: \n'\
    'source                     TS        e0      pivot' % (min_ts, max_beta)
    need_refit =False
    if fit_sources is None:
        fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]

    for source in fit_sources:
        model = source.spectral_model
        try:
            ts, e0, pivot = roi.TS(source.name),model.e0, model.pivot_energy()
        except Exception, e:
            print 'source %s:exception %s' %(source.name, e)
            continue
            
        if pivot is None: continue
        if model.name=='LogParabola': e0 = model[3]
        elif model.name=='ExpCutoff':
            e0 = model.e0
        else:
            continue
        print '%-20s %8.0f %9.0f %9.0f '  % (source.name, ts, e0, pivot),
        if ts < min_ts: 
            print 'TS too small'
            continue
        if model.name=='LogParabola':
            if model[2]>max_beta: 
                print 'beta= %.2f too large' %(model[2])
                continue #very 
            else:
                model.free[2]=False # make sure beta fixed?

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
    def __init__(self, roi, dampen, tol=0.5):
        self.dampen=dampen
        self.roi = roi
        self.ipar = roi.get_parameters()[:] # get copy of initial parameters
        self.initial_logl = roi.logl =  -roi(self.ipar)
        self.tol = tol
    def __call__(self):
        fpar = self.roi.get_parameters()
        if self.dampen!=1.0 and  len(fpar)==len(self.ipar): 
            dpar = self.ipar+self.dampen*(fpar-self.ipar)
            self.roi.set_parameters(dpar)
            self.roi.logl=-self.roi(dpar)
            # check for change, at least 0.5
            return abs(self.initial_logl-self.roi.logl)>self.tol
        return False
    

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
    countsplot_tsmin = kwargs.pop('countsplot_tsmin', 100) # minimum for counts plot
    damp = Damper(roi, dampen)
    
    if outdir is not None: 
        counts_dir = os.path.join(outdir,counts_dir)
        logpath = os.path.join(outdir, 'log')
        if not os.path.exists(logpath): os.mkdir(logpath)
        outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
        print  '='*80
        print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    else: outtee=None

    if not os.path.exists(counts_dir):
        try: os.makedirs(counts_dir) # in case some other process makes it
        except: pass
    init_log_like = roi.log_like()
    roi.print_summary(title='before fit, logL=%0.f'% init_log_like)
    fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    if len(roi.get_parameters())==0:
        print '===================== nothing to fit========================'
    else:
        if dampen>0:
            fit_kw = kwargs.get('fit_kw', {})
            try:
                if diffuse_only:
                    ndiff = len([n for n in roi.parameter_names if n.split('_')[0] in ('ring','isotrop', 'limb')])
                    roi.summary(range(ndiff), title='Before fit to diffuse components')
                    fit_kw.update(select=range(ndiff))
                    roi.fit(**fit_kw)
                else:
                    roi.fit(**fit_kw)
                    if  abs(roi.log_like() - init_log_like)>1.0 :
                        roi.print_summary(title='after global fit, logL=%0.f'% roi.log_like())
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
                    roi.print_summary(title='after damping with factor=%.2f, logL=%0.f'%(dampen,roi.log_like()))
                else:
                    print 'No damping requested, or difference in log Likelihood < %f' % damp.tol
            except Exception, msg:
                print '============== fit failed, aborting!! %s'%msg
                #pickle_dump(roi, fit_sources, os.path.join(outdir, 'pickle'), pass_number, failed=True)
                return False
        else:
            print 'No fit requested'
    
    outdir     = kwargs.pop('outdir', '.')
    associator= kwargs.pop('associate', None)
    if type(associator)==types.StringType and associator=='None': associator=None
    def getdir(x ):
        if not kwargs.get(x,False): return None
        t = os.path.join(outdir, kwargs.get(x))
        if not os.path.exists(t): os.mkdir(t)
        return t
    sedfuns.makesed_all(roi, sedfig_dir=getdir('sedfig_dir'))
    if localize:
        print 'localizing and associating all sources with variable...'
        q, roi.quiet = roi.quiet,False
        localization.localize_all(roi, tsmap_dir=getdir('tsmap_dir'), associator = associator)
        roi.quiet=q

    try:
        ax = counts.stacked_plots(roi,None, tsmin=countsplot_tsmin)
        cts=counts.get_counts(roi)
        chisq = cts['chisq']
        print 'chisquared for counts plot: %.1f'% chisq
 
        if counts_dir is not None:
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
    """ special processor to perform localization for tests """
    outdir = kwargs.get('outdir')
    tsmin = kwargs.get('tsmin', 10)
    logpath = os.path.join(outdir, 'log')
    locdir = os.path.join(outdir, 'localization')
    if not os.path.exists(locdir): os.mkdir(locdir)
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
    emin = kwargs.pop('emin', None)
    if emin is not None:
        roi.select_bands(emin=emin)
    # extract only relevant keywords
    loc_kw =dict() 
    for x in localization.Localization.defaults:
        kw = x[0]
        if kw in kwargs:
            loc_kw[kw]=kwargs[kw]
    localization.localize_all(roi, **loc_kw)
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
    ts_min   = kwargs.get('ts_min', 25)
    sedinfo = os.path.join(outdir, 'sedinfo')
    if not os.path.exists(sedinfo): os.mkdir(sedinfo)
    ptsources = [(i,s) for i,s in enumerate(roi.sources) if s.skydir is not None and np.any(s.spectral_model.free)]
    bgsources = [(j,s) for j,s in enumerate(roi.sources) if s.skydir is None]
    print 'evaluating %d sources' %len(ptsources)
    def getdir(x ):
        if not kwargs.get(x,False): return None
        t = os.path.join(outdir, kwargs.get(x))
        if not os.path.exists(t): os.mkdir(t)
        return t
    sedfuns.makesed_all(roi, sedfig_dir=getdir('sedfig_dir'))

    for pti,source in ptsources:
        try:
            ts = roi.TS(source.name)
            if ts<ts_min: continue
            sf = sedfuns.SourceFlux(roi, source.name, )
            sedrecs = [sedfuns.SED(sf, event_class).rec for event_class in (0,1,None)]
                     
        except Exception,e:
            print 'source %s failed flux measurement: %s' % (source.name, e)
            raise
        fname = source.name.replace('+','p').replace(' ', '_')+'_sedinfo.pickle'
        fullfname = os.path.join(sedinfo, fname)
        bgdensity = []
        try:
            for j, bgsource in bgsources:
                bgdensity.append(np.array([ 
                    ( band[j].pix_counts * band[pti].pix_counts ).sum() / band[pti].counts\
                                for band in roi.selected_bands],np.float32) )
        except Exception, e:
            print 'failed bgdensity for source %s: %s' % (source.name, e)
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
                bts   =np.vstack([sr.ts     for sr in sedrecs]),
                bgdensity = bgdensity,
                bgnames = [s.name for ii,s in bgsources],
                ),
            open(fullfname, 'w'))
        print 'wrote sedinfo pickle file to %s' % fullfname
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
    outdir= kwargs.get('outdir')
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
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))

    limbdir = os.path.join(outdir, kwargs.get('limbdir', 'limb'))
    if not os.path.exists(limbdir): os.mkdir(limbdir)
    refit = kwargs.get('refit', True)
    limb = roi.get_model('limb')
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

    
def check_seeds(roi, **kwargs):
    """ Evaluate a set of seeds: fit, localize with position update, fit again
    """
    outdir = kwargs.get('outdir')
    prefix = kwargs.get('prefix')
    tsmap_dir=kwargs.get('tsmap', None)
    tsmin = kwargs.pop('tsmin', 10)
    seedcheck_dir = kwargs.get('seedcheck_dir', 'seedcheck')
    if not os.path.exists(seedcheck_dir): os.mkdir(seedcheck_dir)
    associator= kwargs.pop('associate', None)
    logpath = os.path.join(outdir, 'log')
    outtee = OutputTee(os.path.join(logpath, roi.name+'.txt'))
    print  '='*80
    print '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (roi.name,))
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
    for s in seed_sources: 
        s.ts=ts = roi.TS(s.name)
        sfile = os.path.join(seedcheck_dir, s.name.replace(' ', '_').replace('+','p')+'.pickle')
        pickle.dump(s, open(sfile, 'w'))
        print 'wrote file %s' %sfile
    outtee.close()

    
cat = None
def gtlike_compare(roi, **kwargs):
    """
    compare with gtlike
    """
    def cat_model(cs):
        from uw.like import Models
        st = cs.field('SpectrumType')
        flux,index,cutoff,b,pivot,beta=[cs.field(f) for f in 'Flux_Density Spectral_Index Cutoff Index2 Pivot_Energy beta'.split()]
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
    global cat
    if cat is None:
        import pyfits
        cat = pyfits.open(os.path.expandvars('$FERMI/catalog/gll_psc3yearclean_v1.fit'))[1].data
    outdir = kwargs.pop('outdir')
    sed_dir = os.path.join(outdir, kwargs.pop('sed_dir', 'gtlike/sed'))
    model_dir=os.path.join(outdir, kwargs.pop('model_dir', 'gtlike/models'))
    for t in (sed_dir,model_dir):
        if not os.path.exists(t): os.makedirs(t)
    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    catmodels=[]
    for source in sources:
        cselect = np.array([np.any(s.field('NickName').startswith(source.name.replace(' ',''))) for s in cat])
        if sum(cselect)!=1:
            print 'did not find source %s' %source.name
            continue
        catsource = cat[cselect][0]
        try:
            catmodel = cat_model(catsource); 
        except:
            print 'Source %s failed' %source.name
            catmodel = None
        
        ptts = roi.TS(source.name)
        sed = roi.plot_sed(source.name, butterfly=False, annotate=None, fit_kwargs=dict(label='pointlike: %.0f'%ptts, color='orange', lw=2))
        axes = plt.gca()
        if catmodel is not None:
            s = roi.get_source(source.name)
            t=s.spectral_model
            s.spectral_model = catmodel
            gtts = roi.TS(source.name)
            s.spectral_model = t
            sed.plot_model(axes, catmodel, np.logspace(2, np.log10(3.16e4)), False, label='gtlike: %.0f'%gtts, color='g', lw=2)
        else:
            gtts=-1
        axes.legend()
        plt.setp(axes, xlim=(100, 31.6e3), ylim=(0.1,1000))
        outfile = os.path.join(sed_dir, '%s_sed.png' % (source.name.replace(' ','_').replace('+','p')))
        plt.savefig(outfile)
        print 'wrote file %s' % outfile
        catmodels.append(dict(name=source.name, m_gt=catmodel, m_pt=t, info=np.array(catsource), 
            catinfo=catsource, ts_pt=ptts, ts_gt=gtts))
    outfile = os.path.join(model_dir, '%s.pickle'%roi.name)
    pickle.dump(catmodels,open(outfile, 'w'))
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
  

