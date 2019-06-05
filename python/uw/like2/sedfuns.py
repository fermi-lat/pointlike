# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for ROI analysis - Spectral Energy Distribution functions

"""
import os, pickle
from collections import OrderedDict
import numpy as np
import pandas as pd

from uw.utilities import ( keyword_options)
from . import ( plotting, tools, loglikelihood, sources, bands)
# 2/decade above 31.6 GeV
energybins=np.concatenate( [np.logspace(2,4.25,10), np.logspace(4.5,6,4)])
       
class SED(tools.WithMixin):
    """ measure the energy flux vs. energy for a given source
    It can be set for all bands, just the band(s) at given energy, or a single band:
        see select_band
    Note that it supports the 'with' 
    """
    def __init__(self, rstat, source_name, quiet=False):
        """ rstat : ROIstat object
            source_name : name of one of the sources in the SourceList, or None
        """
        self.rs = rstat
        self.rs.quiet=quiet
        self.func = self.rs.energy_flux_view(source_name, bound=-20) # note very low bound
        self.source_name = source_name
        # make a list of energies with data; only have info if there is data in the ROI
        self.energies = list(set(np.array([b.band.energy for b in self.rs if b.pixels>0])))  

        # combines the bands above 100 GeV 
        emax = rstat[-1].band.emax
        global energybins
        self.energybins=filter(lambda e: e<=emax, energybins)
    
    def full(self):
        try:
            fp = self.select(None)
            self.full_poiss = fp.poiss
        except Exception, msg:
            print 'Failed poisson fit to source %s: "%s"' % (self.source_name, msg)
            raise
        return fp.poiss, fp.maxdev
        
    def __repr__(self):
        return '{}.{} : SED analysis for source "{}", using {:d} selected bands, {:.0f}-{:.0f} MeV'.format(
                    self.__module__, self.__class__.__name__,self.source_name,len(self.rs.selected), 
                    self.rs.emin, self.rs.emax)
    
    def select(self, index, event_type=None, poisson_tolerance=0.20, 
        elow=None, ehigh=None,**kwargs):
        """ Select an energy band or bands
        parameters:
            index: None or integer
                an index into the list of energies; if None, select all bands
                and use the current spectral model, otherwise a powerlaw to 
                represent an model-independent flux over the band.
            event_type : None or integer
                if None, select all, otherwise 0/1 for front/back, 2-5 for PSF0-3
            elow, ehigh : None or float
                If set, use to select bands, allowing combined
                For this case only, check to see if there is any data, return None if not 

        returns an equivalent Poisson object.
        """
        if index is None and elow is None and event_type is None:
            self.rs.select()
            self.func.set_energy(None)# =func = self.rs.energy_flux_view(self.source_name)
        elif index is not None:
            self.rs.select(index, event_type)
            energies = self.rs.energies
            assert len(energies)==1
            energy = self.rs.energies[0]
            self.func.set_energy(energy)
            assert self.func(0) != self.func(1), 'Function not variable? energy %.0f' % energy
        else:
            # case to perhaps combine bands
            self.rs.select(event_type=event_type, elow=elow, ehigh=ehigh)
            has_data = np.any([b.band.has_pixels for b in self.rs.selected])
            if not has_data: return None
            energy = np.sqrt(elow*ehigh)
            self.func.set_energy(energy)
        pf = loglikelihood.PoissonFitter(self.func, tol=poisson_tolerance, **kwargs)
        return pf

    def all_poiss(self, event_type=None, tol=0.1, debug=False):
        """ return array of Poisson objects for each energy band """
        pp = []
        for i,e in enumerate(self.energies):
            if debug: print '%3i %8.0f' % (i,e),
            try:
                pf = self.select(i, event_type=event_type,poisson_tolerance=tol)
                pp.append(pf.poiss)
                if debug: print pf
            except Exception, msg:
                print 'Fail poiss fit for %.0f MeV: %s ' % (e,msg)
                pp.append(None)
                
        self.restore()
        return np.array(pp)
        
    def sed_rec(self, event_type=None, tol=0.1):
        """ return a numpy.recarray with values for each band (or set of bands)
           elow ehigh       -- energy limits
           flux lflux uflux -- flux at max, upper and lower 1-sigma
           ts               -- Test Statistic for the band
           mflux            -- Flux predicted by the model
           npred            -- number predicted photons for max L 
           pindex           -- photon index at central energy
           delta_ts         -- TS difference, fit-model
           pull             -- signed square root of delta_ts
           zero_fract       -- predicted fraction of the time expect zero flux
        """
        names = 'elow ehigh flux lflux uflux npred pindex ts mflux  delta_ts pull maxdev zero_fract'.split()
        rec = tools.RecArray(names, dtype=dict(names=names, formats=['>f4']*len(names)) )
        
        ebins = self.energybins
        if event_type is not None:
            ebins = filter(lambda e: e>=bands.event_type_min_energy[event_type], ebins)
        for i,(elow,ehigh) in enumerate(zip(ebins[:-1], ebins[1:])):
        #for i,energy in enumerate(self.energies):
                
            try:
                pf = self.select(None, elow=elow,ehigh=ehigh, 
                    event_type=event_type, poisson_tolerance=tol)
                xlo,xhi = self.rs.emin, self.rs.emax
            except Exception, msg:
                print 'Fail poiss fit for %.0f-%.0f MeV: %s ' % (elow,ehigh,msg)
                rec.append(elow,ehigh, 0, 0, np.nan, 0,0,0, np.nan, np.nan, np.nan, np.nan,     np.nan )
                continue
            if pf is None: # no data
                rec.append(elow,ehigh, 0, 0, np.nan, 0,0,0, np.nan, np.nan, np.nan, np.nan,     np.nan )
                continue
            elif np.isnan(pf.wprime):
                print 'Fail poiss fit for %.0f-%.0f MeV: %s ' % (elow,ehigh,'bad poiss')
                rec.append(elow,ehigh, 0, 0, np.nan, 0,0,0, np.nan, np.nan, np.nan, np.nan, np.nan )
                continue
            
            w = pf.poiss
            err = pf.maxdev
            lf,uf = w.errors
            maxl  = w.flux
            mf    = self.func.eflux
            self.func(maxl) # set to maxl for npred
            npred = sum([bs[self.source_name].counts for bs in self.rs.selected])
            
            # get spectral function evaluate exponential slope by finite difference
            m = self.rs.get_model(self.source_name)
            x = np.sqrt(xlo*xhi)
            delta=0.01 # 1%
            pindex= (1-m((1+delta)*x)/m(x))/delta
            
            delta_ts = 2.*(self(maxl) - self(mf) )
            zf =   w.zero_fraction()
            if lf>0 :
                pull = np.sign(maxl-mf) * np.sqrt(max(0, delta_ts))
                assert not np.isnan(pull), 'row {}: pull = {}'.format(i,pull)
                rec.append(xlo, xhi, maxl, lf, uf, npred, pindex, w.ts, mf, delta_ts, pull, err, zf)
            else:
                pull = -np.sqrt(max(0, delta_ts))
                rec.append(xlo, xhi, 0, 0, w.cdfcinv(0.05), 0,pindex, 0, mf, delta_ts, pull, err, zf )
            
        self.restore()
        return rec()

    def data_frame(self, event_type=None, tol=0.1):
        """DataFrame summary of the sed_rec"""
        si = self.sed_rec(event_type,tol)
        r =pd.DataFrame(
            dict(elow=si.elow, ehigh=si.ehigh, 
                flux=si.flux.round(2), TS=si.ts.round(1), lflux=si.lflux.round(2),
                npred= si.npred.round(1),
                pindex=si.pindex.round(2),
                uflux=si.uflux.round(2), mflux=si.mflux.round(2), pull=si.pull.round(2), zf=si.zero_fract.round(3)), 
            index=np.array(np.sqrt(si.elow*si.ehigh),int), 
            columns='elow ehigh flux lflux uflux mflux npred pindex TS pull zf'.split())
        r.index.name='energy'
        return r

        
    def restore(self):
        self.rs.select()
        self.func.restore()
        
    def __call__(self, eflux):
        """eflux : float or array of float
            energy flux in eV units
        """
        return self.func(eflux)
        
    def plots(self):
        import matplotlib.pylab as plt
           
        fig, axx = plt.subplots(4,4, figsize=(12,12), sharey=True)
        for i, ax in enumerate(axx.flatten()):
            if i >= len(self.energies):
                ax.set_visible(False)
                continue
            pf = self.select(i)
            pf.plot(ax)
            ax.set_title('%d MeV' %( int(self.func.energy),), size=10)
            ax.set_ylim(0,1)
        self.restore()
        fig.suptitle('Binned likelihood plots for '+self.source_name, size=14)
        return fig
 
def sed_table(roi, source_name=None, event_type=None, tol=0.1):
    """
    Return a DataFrame
    """
    source = roi.sources.find_source(source_name)
    
    if isinstance(event_type,str):
        etname = event_type.lower()
        if etname=='all': event_type=None
        elif etname in roi.config.event_type_names:
            event_type = roi.config.event_type_names.index(etname)
        else:
            raise Exception('event type name %s not recognized' %event_type_name)
            
    with SED(roi, source.name) as sf:
        return sf.data_frame(event_type=event_type, tol=tol)


def norm_table(roi, source_name=None, event_type=None, tol=0.25, ignore_exception=True):
    """
    Return a DataFrame table for the given source with Poisson results
    """
    source = roi.sources.find_source(source_name)
    #print 'table for {}'.format(source.name)
    roi.select()
    energies = roi.energies
    poiss_list = dict()
    with roi.normalization_view(source.name) as nv:
        for i,energy  in enumerate(energies):
            roi.select(i, event_type)
            try:
                p = loglikelihood.PoissonFitter(nv, tol=tol)
                poiss_list[int(energy)] = p.normalization_summary()
            except Exception, msg:
                print 'Fail for %.f: %s' % (energy, msg)
                if not ignore_exception: raise
                poiss_list[int(energy)]= {}
    roi.select()
    ret = pd.DataFrame(poiss_list, index='maxl lower upper ts err'.split() ).T
    ret['pull'] = (ret.maxl-1)/ret.err
    ret.index.name='energy'
    return ret
                
def residual_tables(roi, tol=0.3, types=None, globals=None, locals=None):
    """ make residual tables for global and local sources
    """
    if types is None: types = ['all']+ list(roi.config.event_type_names)
    if globals is None: globals = filter(lambda s: s.isglobal,roi.sources)
    residuals = dict()
    for source in globals:
        yy = residuals[source.name] = dict()
        yy['model'] = source.model    
        for et in types:
            print source.name, et
            yy[et] = norm_table(roi, source.name,et, tol)
            
    if locals is None: locals = filter(lambda s: np.any(s.model.free) and not s.isglobal, roi.sources)
    for source in locals:
        yy = residuals[source.name] = dict()
        yy['model'] = source.model    
        for et in ('all', 'front', 'back'):
            print source.name, et
            yy[et] = sed_table(roi, source.name, et, tol)
    return residuals



def print_sed(roi, source_name=None):
    source = roi.get_source(source_name)
    t = pd.get_option('display.float_format')
    pd.set_option('display.float_format', lambda x: '%.1f'%x)
    print sed_table(roi, source_name)
    pd.set_option('display.float_format', t)
               

def makesed_all(roi, source_name='all', **kwargs):
    """ add sed information to each free local source.
    Except: can give an individual source name
    
    kwargs:
        sedfig_dir : string or None
            if string, a folder name in which to put the figures
        showts : bool
        ndf : int
            default 10, for fit quality
    other kwargs passed to sed.Plot().__call__
    """
    from scipy import stats # for chi2 
    sedfig_dir = kwargs.pop('sedfig_dir', None)
    if sedfig_dir is not None and sedfig_dir[0]=='$':
        sedfig_dir = os.path.expandvars(sedfig_dir)
    ndf = kwargs.pop('ndf', 10) 
    if sedfig_dir is not None and not os.path.exists(sedfig_dir): os.mkdir(sedfig_dir)
    showts = kwargs.pop('showts', True)
    poisson_tolerance = kwargs.pop('poisson_tolerance', 0.50)
    initw = roi.log_like()

    if source_name=='all':
        sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    else:
        sources = [roi.get_source(source_name)]
    print 'sources:', [s.name for s in sources]
    for source in sources:
        with SED(roi, source.name, ) as sf:
            print source.name,':',
            try:
                source.sedrec = sf.sed_rec( tol=poisson_tolerance)
                source.ts = roi.TS(source.name)
                qual = sum(source.sedrec.pull**2)
                pval = 1.- stats.chi2.cdf(qual, ndf)
                if sedfig_dir is not None:
                    annotation =(0.04,0.88, 'TS=%.0f\npvalue %.1f%%'% (source.ts,pval*100.)) if showts else None 
                    plotting.sed.stacked_plots(sf,  #gev_scale=True, energy_flux_unit='eV',
                         galmap=source.skydir, outdir=sedfig_dir, 
                            annotate=annotation, **kwargs)
                        
            except Exception,e:
                print '***Warning: source %s failed flux measurement: %s' % (source.name, e)
                #raise
                source.sedrec=None
    curw= roi.log_like()
    assert abs(initw-curw)<0.1, \
        'makesed_all: unexpected change in roi state after spectral analysis, from %.1f to %.1f' %(initw, curw)

def add_flat_sed(roi, source_name=None, cols='flux lflux uflux ts'.split()):
    """For the given source (or 'ALL'), look for additional unmodelled photons by defining a flat source
        at the same position, measuring its SED
        Add fit flux with limits, and TS to the sedrec.
    """

    def do_one(s):
        roi.add_source(name='temp', skydir=s.skydir, model='LogParabola(1e-12, 2, 0, 1e4)')
        ss = roi.get_source('temp')
        t = roi.get_sed()
        roi.del_source('temp')
        df = pd.DataFrame(s.sedrec)
        df_flat=pd.DataFrame(OrderedDict( [(n, t[n].astype(float)) for n in cols]))
        for col in cols:
            df[col+'_flat'] = df_flat[col]
        s.sedrec = rec = df.to_records(index=False)
        rec.dtype.names = map(str, rec.dtype.names) # numpy dtype names must be str 
        return sum(df_flat.ts)

    if source_name!='ALL':
        s = roi.get_source(source_name)
        return do_one(s)
    else:
        sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
        return [(s.name, do_one(s)) for s in sources]


def normalization_poiss(roi, source_name, event_type=None):
    """ return a list of Poisson objects for each energy band
    
    parameters
    ----------
    source_name : str
        name of a source
    event_type : [None | int]
        None for all, int to select, say, front or back
        
    """
    roi.select()
    energies = roi.energies
    poiss_list = []
    with roi.normalization_view(source_name) as nv:
        for i,energy  in enumerate(energies):
            roi.select(i, event_type)
            p = loglikelihood.PoissonFitter(nv, tol=0.25).poiss
            poiss_list.append(p)
    roi.select()
    return poiss_list

def alternate_source(roi, source, name, skydir, model):
    """
    Create a PointSource object, coincident with a current one, but with a different 
    
    roi    : a ROI object
    source : Current Source object, with sedrec
    name   : Name of alternative source
    skydir : its position (presumably very close to source.skydir)
    model  : alternative Model object
    
    Will create, and return a Source object, containing a different sedrec.
    
    """
    roi.get_source(source.name) # make sure selected
    altsrc = sources.PointSource(name=name, skydir=skydir, model=model)
    saved_model = source.spectral_model
    saved_sed = source.sedrec.copy()
    source.spectral_model = model
    altsrc.sedrec = roi.get_sed(update=True)
    source.spectral_model = saved_model
    source.sedrec = saved_sed
    return altsrc
    