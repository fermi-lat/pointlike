"""
roi and source processing used by the roi pipeline
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/processor.py,v 1.11 2011/04/05 21:42:08 burnett Exp $
"""
import os, pickle
import numpy as np
import pylab as plt
from skymaps import SkyDir
from uw.like import counts_plotter, srcid
from uw.utilities import image
from . import plot_sed  
  
def isextended(source):
    return source.__dict__.get(  'spatial_model', None) is not None
    

        
def fix_beta(roi, bts_min=30, qual_min=15,):
    refit=candidate=False
    print 'checking for beta fit...',
    models_to_fit=[]
    for which, model in enumerate(roi.psm.models):
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
        print 'start fit'
        roi.fit()
        roi.print_summary(title='after freeing one or more beta parameters')
        # now check for overflow
        refit = False
        for which,model in enumerate(roi.psm.models):
            if not np.any(model.free): break
            if model.name!='LogParabola': continue
            beta = model[2]
            if beta < 2.5: continue
            model.free[2]=False
            model[2]=2.5
            model.cov_matrix[:]=0 # can screw up fit if bad.
        if refit:
            print 'need to refit: beta too large'
            roi.fit()
            refit=False
            for model in models_to_fit():
                if model[2]>3.0 :
                    print 'reseting model: beta =%.1f too large' % model[2]
                    model[2]=2.; model.free[2]=False
                    refit = True
            if refit:
                print 'refitting again'
                roi.fit()
            
    else:
        print 'none found'
    return refit    
        
def source_band_info(roi, which):
    """ return dictionary of the band ts values and photons, diffuse  
        which: index of source in the ROI, or source name
    """
    # assume done: roi.setup_energy_bands()
    man,i = roi.mapper(which)
    bands = roi.energy_bands
    return dict(
        ts =      np.asarray([band.bandFit(i) for band in bands],np.float32),
        fit_ts =  roi.fit_ts_list(which),
        signal_counts= roi.signal_counts(which),
        )

def pickle_dump(roi, fit_sources, pickle_dir, pass_number, **kwargs):
    """ dump the source information from an ROI constructed from the sources here
    """
    
    name = roi.name.strip()
    fname = kwargs.pop('fname', name)
    filename=os.path.join(pickle_dir,fname+'.pickle')
    if os.path.exists(filename):
        # if the file exists
        output=pickle.load(open(filename))
        curp = output.get('passes',[])
        if curp is None: curp = []
        output['passes']=curp.append(pass_number)
        output['prev_logl'] = output['logl']
        print 'updating pickle file with pass %d' %pass_number
    else:
        output = dict()
        output['passes']=[pass_number]
    output['name'] = name
    output['skydir']  = roi.roi_dir
    output['diffuse'] = roi.dsm.models
    output['diffuse_names'] = roi.dsm.names
    output['logl'] = roi.logl
    output['parameters'] = roi.get_parameters() 
    
    sources=dict()
    output['sources']= sources
    for s in fit_sources:
        try:
            pivot_energy = s.model.pivot_energy()
        except: # if not fit
            pivot_energy = None 
        if pivot_energy is None: pivot_energy=s.model.e0 # could be a pulsar?
        
        roi.setup_energy_bands() 
        sedrec = s.sedrec
        band_info = source_band_info(roi, s.name)
        sources[s.name]=dict(skydir=s.skydir, 
            model=s.model,
            extent = s.__dict__.get('spatial_model', None),
            ts = roi.TS(which=s.name),
            sedrec = sedrec,
            band_info = band_info, 
            band_ts=np.sum(band_info['ts']),
            pivot_energy = pivot_energy,
            ellipse= np.array(s.ellipse,np.float32) if 'ellipse' in s.__dict__ else None,
            associations = s.adict if 'adict' in s.__dict__ else None,
            )
    output.update(kwargs) # add additional entries from kwargs
    f = open(filename,'wb') #perhaps overwrite
    pickle.dump(output,f)
    f.close()
    

def fname(name):
    return name.replace(' ','_').replace('+','p')

def make_sed(roi, source, sedfig_dir, **kwargs):
    plot_sed.PlotSED(source.sedrec)(source.model, source.name)
    # add a galactic map if requested
    image.galactic_map(roi.roi_dir, color='lightblue', marker='s', markercolor='r')
    fout = os.path.join(sedfig_dir, ('%s_sed.png'%fname(source.name)) )
    plt.title(source.name)
    plt.savefig(fout)
    print 'saved SED plot to %s' % fout 

def make_tsmap(roi, source, tsmap_dir, **kwargs):
    """ generate a tsmap for source name in the roi
    """
    adict = source.adict if 'adict' in source.__dict__ else None #kwargs.pop('adict', None)
    tsize = source.ellipse[2]*15. if source.ellipse is not None else 1.1
    name = source.name
    fout = os.path.join(tsmap_dir, ('%s_tsmap.png'%fname(name)) )
    which = source.name 
    try:
        if not isextended(source): 
            roi.qform= None #kluge to ignore the last fit
            tsm=roi.plot_tsmap( which=which, center=source.skydir, name=roi.name,
                outdir=None, catsig=0, size=tsize, 
                pixelsize= tsize/14,
                # todo: fix this
                assoc=adict if adict is not None else None, # either None or a dictionary
                notitle=True, #don't do title
                markersize=10,
                primary_markersize=12,
                )
            if source.ellipse is not None:
                tsm.overplot(source.ellipse)
            tsm.zea.axes.set_title('%s'% name, fontsize=16)  # big title
        else: 
            # yucky branch to avoid trying to generate TS maps for extended sources
            print 'dummy TS map for extended source %s ' % name
            fig = plt.figure(figsize=(5,5))
            plt.title('%s'% name, fontsize=16)
            plt.text(0.2,0.5, '(extended source)', fontsize=16)
            plt.gca().set_axis_off()
        plt.savefig(fout)
        print 'saved tsplot to %s' % fout 
    except Exception, e:
        print 'Failed to calculate tsmap: %s' % e
        fig = plt.figure(figsize=(5,5))
        plt.title('%s'% name, fontsize=16) 
        plt.savefig(fout)
 
def make_association(source, tsf, associate):
    print ' %s association(s) ' % source.name,
    try:    ell = source.ellipse
    except: ell = None
    if ell is None:
        print '...no localization'
        source.adict = None
        return
    assert len(ell)>6, 'invalid ellipse for source %s' % source.name
    try:
        adict = associate(source.name, SkyDir(ell[0],ell[1]), ell[2:5]) 
    except srcid.SrcidError, msg:
        print 'Association error %s' % msg
        raise
        adict=None
    source.adict = adict 
    if adict is not None:
    
        ts_local_max=tsf(source.tsmaxpos)
        adict['deltats'] = [ts_local_max-tsf(d) for d in adict['dir']]
        print '\n   cat         name                  ra        dec         ang     prob    Delta TS'
        #       15 Mrk 501               253.4897   39.7527    0.0013      0.41
        fmt = '   %-10s %-20s%10.4f%10.4f%10.4f%8.2f%8.1f' 
        for i,id_name in enumerate(adict['name']):
            tup = (adict['cat'][i], id_name, adict['ra'][i], adict['dec'][i], adict['ang'][i], 
                    adict['prob'][i],adict['deltats'][i])
            print fmt % tup
    else:
        print '...None  found'

def process_sources(roi, sources, **kwargs):
    outdir     = kwargs.pop('outdir', '.')
    associate= kwargs.pop('associate', None)
    
    if associate is not None and associate!='None':
        for source in sources:
            make_association(source, roi.tsmap(which=source.name), associate)
 
    # add band info to the source objects
    for source in sources:
        source.sedrec = plot_sed.SEDflux(roi, which=source.name).rec
        
    dirs = [kwargs.pop(dirname,None) for dirname in ('sedfig_dir', 'tsmap_dir')]
    procs = (make_sed, make_tsmap)
    for dir,proc in zip(dirs, procs):
        if dir is None: continue
        full_path = os.path.join(outdir,dir)
        if not os.path.exists(full_path): os.makedirs(full_path)
        for source in sources:
            proc(roi, source, full_path, **kwargs)
     

def localize_all(roi,sources):

    roi.qform=None
    for source in sources:
        if isextended(source):
            print 'source %s is extended: not localizing' % source.name
            source.ellipse=None
        else:
            source.tsmaxpos, delta_ts =roi.localize(which=source.name)
            source.ellipse = roi.qform.par[0:2]+roi.qform.par[3:7] +[delta_ts] if roi.qform is not None else None
        
def repivot(roi, fit_sources, min_ts = 16, max_beta=3.0):
        print '\ncheck need to repivot sources with TS>%.0f, beta<%.1f: \n'\
        'source                     TS        e0      pivot' % (min_ts, max_beta)
        need_refit =False

        for source in fit_sources:
            model = source.model
            try:
                ts, e0, pivot = roi.TS(which=source.name),model.e0, model.pivot_energy()
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
    def __init__(self, roi, dampen):
        self.dampen=dampen
        self.roi = roi
        self.ipar = roi.get_parameters()[:] # get copy of initial parameters
        self.initial_logl = roi.logl =  -roi.logLikelihood(self.ipar)
    def __call__(self):
        fpar = self.roi.get_parameters()
        if self.dampen!=1.0 and  len(fpar)==len(self.ipar): 
            dpar = self.ipar+self.dampen*(fpar-self.ipar)
            self.roi.set_parameters(dpar)
            self.roi.logl=-self.roi.logLikelihood(dpar)
            # check for change, at least 0.5
            return abs(self.initial_logl-self.roi.logl)>0.5
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
    pass_number=kwargs.pop('pass_number', 0)
    tables = kwargs.pop('tables', None)

    damp = Damper(roi, dampen)
    
    if outdir is not None: counts_dir = os.path.join(outdir,counts_dir)
    if not os.path.exists(counts_dir):os.makedirs(counts_dir)
    
    roi.print_summary(sdir=roi.roi_dir, title='before fit, logL=%0.f'%roi.logl)
    fit_sources = [s for s in roi.psm.point_sources if np.any(s.model.free)]\
                + [s for s in roi.dsm.diffuse_sources if isextended(s) and np.any(s.model.free)]
    if len(roi.get_parameters())==0:
        print '===================== nothing to fit========================'
    else:
        if dofit:
            try:
                roi.fit(ignore_exception=False)
                roi.print_summary(sdir=roi.roi_dir,title='after global fit, logL=%0.f'% roi.logl)

                if repivot_flag:
                    repivot(roi, fit_sources)
                
                if fixbeta:
                    if not fix_beta(roi):
                        print 'fixbeta requested, but no refit needed, quitting'
                if damp():
                    roi.print_summary(sdir=roi.roi_dir, title='after damping with factor=%.2f, logL=%0.f'%(dampen,roi.logl))
                else:
                    print 'No damping requested, or no difference in log Likelihood'
            except Exception, msg:
                print '============== fit failed, aborting!! %s'%msg
                return False
        else:
            print 'Refit not requested'
    
    if localize: 
        localize_all(roi, fit_sources)
    if tables is not None:
        tables(roi)
    
    process_sources(roi, fit_sources, **kwargs)
    
    ax = counts_plotter.roi_pipeline_counts_plot(roi,None, merge_all=True)
    counts=counts_plotter.get_counts(roi, merge_non_free=True, merge_all=True)
    obs,mod = counts['observed'], counts['total']
    chisq = ((obs-mod)**2/mod).sum()
    print 'chisquared for counts plot: %.1f'% chisq

    ax[1].text(roi.sa.fit_emin*1.1, 0.2,'chisq=%.1f'% chisq)
    if counts_dir is not None:
        fout = os.path.join(counts_dir, ('%s_counts.png'%roi.name) )
        ax[1].figure.savefig(fout)
        print 'saved counts plot to %s' % fout

    
    if outdir is None: return True
    
    pickle_dir = os.path.join(outdir, 'pickle')
    if not os.path.exists(pickle_dir): os.makedirs(pickle_dir)
    pickle_dump(roi, fit_sources, pickle_dir, pass_number,
        initial_logl=damp.initial_logl, 
        counts=counts,
        )
    return True

