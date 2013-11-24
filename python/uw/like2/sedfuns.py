"""
Tools for ROI analysis - Spectral Energy Distribution functions

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sedfuns.py,v 1.20 2013/11/23 16:05:33 burnett Exp $

"""
import os, pickle
import numpy as np
from uw.utilities import fitter,  makerec, keyword_options
from . import ( plotting, tools, loglikelihood)

#class EnergyFluxConverter(object):
#    """ a functor for any Model, which returns the flux internal parameter as a 
#        function of the differential energy flux, in eV units, at the given energy
#    """
#    def __init__(self, model, energy):
#        eflux = model(energy) * energy**2 * 1e6
#        self.ratio = model[0]/eflux
#        assert model[0]==model['norm']
#        self.tointernal = model.mappers[0].tointernal
#        self.bound = model.bounds[0][0]
#        
#    def __call__(self, eflux):
#        if eflux<=0:
#            return self.bound
#        t = self.tointernal(eflux*self.ratio)
#        return max(self.bound, t)
        
class SourceFlux(tools.WithMixin):
    """ measure the energy-dependent flux for a given source
    An object can be passed to LogLikelihood to measure errors, limit, TS
    It can be set for all bands, just the band(s) at given energy, or a single band:
        see select_band
    But beware: it alters the source's spectral_model if so: be sure to call reset when done!
    Note that it supports the 'with' 
    """
    def __init__(self, rstat, source_name, quiet=False):
        """ rstat : ROIstat object
            source_name : name of one of the sources in the SourceList, or None
        """
        self.rs = rstat
        self.rs.quiet=quiet
        self.func = self.rs.energy_flux_view(source_name)
        self.source_name = source_name
        self.full_poiss = self.select(None)
        self.energies = self.rs.energies # the original list
    
    def __repr__(self):
        return '%s.%s : %d bands selected, energy range %.0f-%.0f'% (
                    self.__module__, self.__class__.__name__,len(self.rs.selected), self.emin, self.emax)
    
    def select(self, index, event_type=None, poisson_tolerance=0.05):
        """ Select an energy band or bands
        parameters:
            index: None or integer
                an index into the list of energies; if None, select all bands
                and use the current spectral model, otherwise a powerlaw to 
                represent an model-independent flux over the band.
            event_type : None or integer
                if None, select both front and back, otherwise 0/1 for front/back
                
        returns an equivalent Poisson object
        """
        if index is None:
            self.rs.select()
            self.func =func = self.rs.energy_flux_view(self.source_name)
        else:
            self.rs.select(index, event_type)
            energies = self.rs.energies
            assert len(energies)==1
            energy = self.rs.energies[0]
            func = self.rs.energy_flux_view(self.source_name, energy)
        return loglikelihood.PoissonFitter(func).poiss

    def all_poiss(self, event_type=None):
        """ return array of Poisson objects for each energy band """
        pp = []
        for i,e in enumerate(self.energies):
            pp.append(self.select(i, event_type=event_type))
        self.restore()
        return np.array(pp)
        
    def restore(self):
        self.rs.select()
        self.func.restore()
        
    def __call__(self, eflux):
        """eflux : float or array of float
            energy flux in eV units
        """
        return self.func(eflux)
        
    def plots(self, full=True, x = np.linspace(0, 80, 25)):
        import matplotlib.pylab as plt
           
        fig, axx = plt.subplots(4,4, figsize=(12,12), sharex=True, sharey=True)
        
        for i,ax in enumerate(axx.flatten()):
            if i>=len(self.energies): 
                ax.set_visible(False)
                continue
            ll = self.select(i)
            if full:
                y = self(x)
                ax.plot(x, y-y.max(), '-')
            ax.set_title('%d' % self.selected_energy, size=10)
            y2 = np.array(ll(x))
            ax.plot(x, y2-y2.max() ,'+', label='Poisson fit')
            plt.setp(ax, ylim=(-9,1))
        self.restore()
        return fig

    
        
class SED(object):
    """     
    generates a recarray (member rec) with fields:
        elow ehigh : energy range
        flux lflux uflux : energy flux (eV units) for peak, upper and lower limits
        ts : Test Statistic value for the band
        mflux delta_ts : model flux, and 2*(logl(flux)-logl(mflux), where logl is the log likelihood
        pull : sign(mflux) * sqrt(delta_ts)
    """

    def __init__(self, source_flux, event_class=None, scale_factor=1.0, merge=False,):
        """ 
            source_flux:  SourceFlux object
            event_class : None or int
                None for all, 0/1 for front/back
            merge: bool
                flag to merge adjacent upper and lower bands with upper limits only
                (not implemented yet)
                
        """
        sf = source_flux
        self.scale_factor=scale_factor
        rec = makerec.RecArray('elow ehigh flux lflux uflux ts mflux delta_ts pull'.split())
        self.loglikes = []
        for i,energy in enumerate(sf.energies):
            sf.select_band(i, event_class)
            xlo,xhi = sf.emin,sf.emax
            try:
                w = tools.LogLikelihood(sf)
                lf,uf = w.errors()
                mf = sf.model_flux()
                delta_ts = 2.*(sf(w.maxl) - sf(mf) )
            except Exception, e:
                print 'Failed likelihood analysis for source %s energy %.0f: %s' %(sf.source.name, energy, e)
                rec.append(xlo, xhi, np.nan, np.nan, np.nan, 0, np.nan, np.nan)
                continue
            if lf>0 :
                pull = np.sign(w.maxl-mf) * np.sqrt(max(0, delta_ts))
                rec.append(xlo, xhi, w.maxl, lf, uf, w.TS(), mf, delta_ts, pull)
            else:
                pull = -np.sqrt(max(0, delta_ts))
                rec.append(xlo, xhi, 0, 0, w.upper_limit(), 0, mf, delta_ts, pull )
            
        self.rec = rec()
        sf.restore() # restore model normalization
    
    def __str__(self):
        n  = len(self.rec.dtype)-2
        return ((2*'%10s'+n*'%8s'+'\n') % self.rec.dtype.names)\
             +'\n'.join( [(2*'%10.0f'+n*'%8.1f') % tuple(row) for row in self.rec])


class DiffuseLikelihood(fitter.Fitted):
    """ implement likelihood function of diffuse normalization
    """
    def __init__(self, rstat, names=('ring','isotrop'), quiet=True):
        """ rstat : ROIstat object
            source_name : name of one of the sources in the SourceList
        """
        self.rs = rstat
        self.rs.quiet=quiet
        self.models = [rstat.get_model(name) for name in names]
        self.energies = np.sort(list(set([ sm.band.e for sm in rstat.all_bands])))
        self.saved_pars = self.get_parameters()
        self.selected_bands = self.rs.selected_bands.copy()
        self.eopt = [bl[0].energy for bl in self.selected_bands]
        
    def restore(self):
        self.set_parameters(self.saved_pars)

    def select_band(self, index, event_class=None):
        """ Select an energy band or bands
        
        Parameters
        ----------
        index: None or integer
                an index into the list of energies; if None, select all bands
                and use the current spectral model, otherwise a powerlaw to 
                represent an model-independent flux over the band.
            event_class : None or integer
                if None, select both front and back, otherwise 0/1 for front/back
        """
        if index==None: #select all bands, use input model
            self.rs.selected_bands = self.selected_bands
        else:
            # selected band(s) at a given energy: use a powerlaw 
            self.selected_energy = energy =self.energies[index]
            class_select = lambda x : True if event_class is None else x==event_class
            self.rs.select_bands(lambda b: b.e==energy and class_select(b.ec))
            assert len(self.rs.selected_bands)>0, 'did not find any bands for energy %.1f' % energy
        self.emin = np.min([bandlike.band.emin for bandlike in self.rs.selected_bands])
        self.emax = np.max([bandlike.band.emax for bandlike in self.rs.selected_bands])
    
    def log_like(self):
        return self.rs.log_like()
    def __call__(self, pars):
        self.set_parameters(pars)
        return -self.log_like()
        
    def get_parameters(self):
        return np.array([m[0] for m in self.models])
    def set_parameters(self, pars):
        def setpar(model,par): model[0]=par
        map(setpar, self.models, pars)
        
    def fit(self, index, event_class=None, **kwargs):
        """ 
        perform a fit
        
        parameters
            index : None or integer
                use to select energy band, or all
        """
        self.select_band(index, event_class)
        self.restore()
        initial_loglike= self.log_like()
        mm = fitter.Minimizer(self)
        t = mm(**kwargs)
        self.restore() 
        return (-t[0] - initial_loglike, t[1], t[2])
    
    def multifit(self, event_class=None, emax=None, **kwargs):
        """ run fit for all energy bands, make dict to return with results
        """
        loglike=[]
        values = []
        errors = []
        energies = filter(lambda e:e<emax, self.energies) if emax is not None else self.energies
        #print ' fitting energies', energies
        for index,energy in enumerate(energies):
            try:
                t = self.fit(index, event_class, **kwargs)
            except Exception, e:
                # TODO: refit with only one parameter?
                print 'Failed fit for e=%.0f: %s' % (energy , e)
                nm = len(self.models)
                t = [0, np.zeros(nm), np.zeros(nm)]
            loglike.append(t[0])
            values.append(t[1])
            errors.append(t[2])
        # get optimized energies for reference NOTE: replace that 0 with source name
        return dict(energies=energies, loglike=loglike, event_class=event_class, values=values, errors=errors,
            eopt=self.eopt)

def makesed_all(roi, **kwargs):
    """ add sed information to each free local source
    
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
    ndf = kwargs.pop('ndf', 10) 
    if sedfig_dir is not None and not os.path.exists(sedfig_dir): os.mkdir(sedfig_dir)
    showts = kwargs.pop('showts', True)
    initw = roi.log_like()

    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    for source in sources:
        try:
            sf = SourceFlux(roi, source.name, )
            source.sedrec = SED(sf, merge=False).rec
            source.ts = roi.TS(source.name)
            qual = sum(source.sedrec.pull**2)
            pval = 1.- stats.chi2.cdf(qual, ndf)
            if sedfig_dir is not None:
                annotation =(0.04,0.88, 'TS=%.0f\npvalue %.1f%%'% (source.ts,pval*100.)) if showts else None 
                plotting.sed.stacked_plots(roi, source.name, #gev_scale=True, energy_flux_unit='eV',
                     galmap=source.skydir, outdir=sedfig_dir, 
                        annotate=annotation, **kwargs)
                    
        except Exception,e:
            print 'source %s failed flux measurement: %s' % (source.name, e)
            raise
            source.sedrec=None
    curw= roi.log_like()
    assert abs(initw-curw)<0.1, \
        'makesed_all: unexpected change in roi state after spectral analysis, from %.1f to %.1f' %(initw, curw)


def test_model(roi, model, source_name=None, **kwargs): 
    """
    
    """
    figdir,seddir = [kwargs.pop(x, None) for x in ('figdir', 'seddir')]
    for x in figdir,seddir:
        if x is not None and not os.path.exists(x):
            os.mkdir(x)
            
    alternate = kwargs.pop('alternate', 'alternate')
   
    sedrec = roi.get_sed(source_name)
    fitqual=sum(sedrec.delta_ts)
    source = roi.get_source()    
    plx = plotting.sed.Plot(source, gev_scale=True, energy_flux_unit='eV')
    plx( galmap=source.skydir, fit_kwargs=dict(color='r', lw=2, label='best fit (%.1f)' % fitqual),
                    **kwargs)
    old_model = roi.set_model(model)
    with SourceFlux(roi, source_name) as sf:
        tsedrec = SED(sf).rec
    alt_fitqual = sum(tsedrec.delta_ts)
    plx.plot_model(model, color='green', lw=2,  ls='--', label='%s (+%.1f)' %(alternate,alt_fitqual-fitqual) )
    plx.axes.legend(loc='upper left', prop=dict(size=8) )
    if figdir is not None:
        plx.savefig(figdir)
    if seddir is not None:
        pickle.dump(tsedrec, open(os.path.join(seddir, 
            source.name.replace(' ','_').replace('+','p')+'.pickle'), 'w'))
    roi.set_model(old_model)
    return alt_fitqual
    
def sed_table(roi, source_name=None, emax=1e6, **kwargs):
    """
    Return a pandas DataFrame with spectral information for the source
    Columns, with energy fluxes in eV units, are:
        flux, lflux, uflux : measured flux, lower and upper uncertainty, or 0,0, 95% limit
        mflux : predicted flux for this bin
        TS : Test Statistic for the signal
        pull : signed square root of the TS difference for the model
    Index is mean energy in MeV
    """
    import pandas as pd

    sedrec = roi.get_sed(source_name, **kwargs)
    si = sedrec[sedrec.elow<emax]
    pull = np.sign(si.flux-si.mflux) * np.sqrt( si.delta_ts.clip(0,100) )
    return pd.DataFrame(dict(flux=si.flux.round(1), TS=si.ts.round(1), lflux=si.lflux.round(1), 
        uflux=si.uflux.round(1), model=si.mflux.round(1), pull=pull.round(2) ),
            index=np.array(np.sqrt(si.elow*si.ehigh),int), columns='flux lflux uflux model TS pull'.split())
    