"""
roi and source processing used by the roi pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/processor.py,v 1.13 2012/01/11 14:03:25 burnett Exp $
"""
import os, time, sys, types
import cPickle as pickle
import numpy as np
import pylab as plt
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
        
def fix_beta(roi, bts_min=30, qual_min=15,):
    refit=candidate=False
    print 'checking for beta fit...'
    models_to_fit=[]
    for source in roi.sources:
        model = source.spectral_model
        which = source.name
        if not np.any(model.free): break
        if model.name!='LogParabola': continue
        if model.free[2]: continue  # already free 
        if not candidate:
            print 'name            beta   band_ts  fitqual'
        candidate=True
        beta = model[2]
        band_ts, ts = roi.band_ts(which=which), roi.TS(which=which)
        print '%-20s %10.2f %10.1f %10.1f ' %(roi.psm.point_sources[which].name,beta, band_ts, band_ts-ts),
        if beta>=3.0: print 'beta>1 too large'; continue
        if band_ts<bts_min: print 'band_ts< %.1f'%bts_min; continue # not significant
        if band_ts-ts < qual_min and beta<=0.01:print 'qual<%.1f' %qual_min; continue # already a good fit
        print '<-- select to free beta' # ok, modify
        model.free[2]=True
        models_to_fit.append(model) # save for later check
        refit = True
    if refit:    
        print 'start refit with beta(s) freed...'
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
            roi.print_summary(title='re-refit after freeing one or more beta parameters')
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
            ts = roi.TS(s.name),
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
        

def repivot(roi, fit_sources, min_ts = 16, max_beta=3.0):
        print '\ncheck need to repivot sources with TS>%.0f, beta<%.1f: \n'\
        'source                     TS        e0      pivot' % (min_ts, max_beta)
        need_refit =False

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
 
            if pivot < roi.fit_emin[0] or pivot > roi.fit_emax[0]:
                print 'bad pivot energy, not in range (%.0f, %.0f)' % (roi.fit_emin[0], roi.fit_emax[0])
                continue
            if abs(pivot/e0-1.)<0.05:
                print 'converged'; continue
            print 'will refit'
            need_refit=True
            model.set_e0(pivot)
        if need_refit:
            roi.fit()
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
    localize = kwargs.pop('localize', True)
    repivot_flag  = kwargs.pop('repivot', True)
    fixbeta  = kwargs.pop('fix_beta', False)
    dampen   = kwargs.pop('dampen', 1.0)  # factor to adjust parameters before writing out
    counts_dir = kwargs.pop('counts_dir', 'counts_dir')
    dofit   =  kwargs.pop('dofit', True)
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

    if not os.path.exists(counts_dir):os.makedirs(counts_dir)
    init_log_like = roi.log_like()
    roi.print_summary(title='before fit, logL=%0.f'% init_log_like)
    fit_sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    if len(roi.get_parameters())==0:
        print '===================== nothing to fit========================'
    else:
        if dofit:
            fit_kw = kwargs.get('fit_kw')
            try:
                if diffuse_only:
                    ndiff = len([n for n in roi.parameter_names if n.split('_')[0] in ('ring','isotrop')])
                    roi.summary(range(ndiff), title='Before fit to diffuse components')
                    fit_kw.update(select=range(ndiff))
                    roi.fit(**fit_kw)
                else:
                    roi.fit(**fit_kw)
                    if  abs(roi.log_like() - init_log_like)>1.0 :
                        roi.print_summary(title='after global fit, logL=%0.f'% roi.log_like())
                if repivot_flag:
                    repivot(roi, fit_sources)
                
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
    print 'processing ROI %s ' %roi.name
    outdir   = kwargs.get('outdir')
    ts_min   = kwargs.get('ts_min', 25)
    sedinfo = os.path.join(outdir, 'sedinfo')
    if not os.path.exists(sedinfo): os.mkdir(sedinfo)
    ptsources = [(i,s) for i,s in enumerate(roi.sources) if s.skydir is not None and np.any(s.spectral_model.free)]
    bgsources = [(j,s) for j,s in enumerate(roi.sources) if s.skydir is None]
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
    outdir= kwargs.get('outdir')
    limbdir = os.path.join(outdir, kwargs.get('limbdir', 'limb'))
    if not os.path.exists(limbdir): os.mkdir(limbdir)
    t = roi.all_bands[1][2]
    fname = os.path.join(limbdir,'%s_limb.pickle' %roi.name)
    pickle.dump( dict(ra=roi.roi_dir.ra(), dec=roi.roi_dir.dec(),
        counts=t.counts, pixel_values=t.pixel_values), open(fname,'w'))
    print 'wrote file to %s' %fname
    
    
    
def sed_roi(roi, **kwargs):
    roi_refit_processor(roi, **kwargs)
    full_sed_processor(roi, **kwargs)
    
    
def cache_diffuse(roi, **kwargs):
    """ roi processor saves the convolved diffuse maps for an roi"""
    print 'processing %s' %roi.name
    outdir = kwargs.get('outdir')
    diffuse_maps = os.path.join(outdir, 'diffuse_maps')
    if not os.path.exists(diffuse_maps): os.mkdir(diffusemaps)
    
    bands = roi.all_bands
