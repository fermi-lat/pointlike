"""
Tools for ROI analysis - Spectral Energy Distribution functions

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sedfuns.py,v 1.9 2013/03/28 22:14:14 burnett Exp $

"""
import os
import numpy as np
from uw.utilities import fitter,  makerec, keyword_options
from . import plotting 
from . import tools
from uw.like import Models

class SourceFlux(object):
    """ measure the energy-dependent flux for a given source
    An object can be passed to LogLikelihood to measure errors, limit, TS
    It can be set for all bands, just the band(s) at given energy, or a single band:
        see select_band
    But beware: it alters the source's spectral_model if so: be sure to call reset when done!
    Note that it supports the 'with' 
    """
    def __init__(self, rstat, source_name=None, quiet=True):
        """ rstat : ROIstat object
            source_name : name of one of the sources in the SourceList, or None
        """
        self.rs = rstat
        self.rs.quiet=quiet
        self.source = rstat.sources.find_source(source_name)
        parname = self.source.name+'_Norm'
        try:
            self.pindex = list(rstat.parameter_names).index(parname)
        except:
            raise Exception('did not find parameter name, %s, for source flux'%parname)
        self.saved_flux = self.source.spectral_model[0]
        self.saved_model = self.source.spectral_model
        self.energies = np.sort(list(set([ sm.band.e for sm in rstat.all_bands])))
        self.all_bands = self.rs.selected_bands.copy() # save list of initial bands selected
        #self.select_band(None)
      
    def __str__(self):
        return 'SourceFlux of %s in ROI %s' %(self.source.name, self.rs.name)
        
    def restore(self):
        """ restore the selected model """
        self.source.spectral_model = self.saved_model
        self.source.spectral_model[0] = self.saved_flux
        self.rs.selected_bands = self.all_bands

    @tools.ufunc_decorator # make this behave like a ufunc
    def __call__(self, eflux):
        """ eflux : double
                energy flux in eV units
        """
        self.source.spectral_model[0] = max(eflux,1e-3)*self.factor
        self.rs.update()
        return self.rs.log_like()
        
    def model_flux(self):
        """ return the energy flux predicted by the saved model, for the current energy
        (select_band must have been called)
        """
        return self.saved_model(self.selected_energy)/self.factor

    def select_band(self, index, event_class=None):
        """ Select an energy band or bands
        parameters:
            index: None or integer
                an index into the list of energies; if None, select all bands
                and use the current spectral model, otherwise a powerlaw to 
                represent an model-independent flux over the band.
            event_class : None or integer
                if None, select both front and back, otherwise 0/1 for front/back
        Sets self.factor as conversion factor from flux to eflux in eV    
        """
        if index==None: #select all (initially selected) bands, use input model
            self.selected_energy = self.source.spectral_model.e0
            self.rs.selected_bands = self.all_bands
            self.factor = self.saved_model[0]/(self.saved_model.eflux*1e6)
        else:
            # selected band(s) at a given energy: use a powerlaw 
            # (but beware: this replaces the spectral model, which must be restored)
            self.selected_energy = energy =self.energies[index]
            self.source.spectral_model = Models.PowerLaw(free=[True,False],p=[1e-11,2.1], 
                        e0=self.selected_energy) 
            class_select = lambda x : True if event_class is None else x==event_class
            self.rs.select_bands(lambda b: b.e==energy and class_select(b.ec))
            assert len(self.rs.selected_bands)>0, 'did not find any bands for energy %.1f' % energy
            self.factor = 1.0/(energy**2*1e6) 
        self.emin = np.min([bandlike.band.emin for bandlike in self.rs.selected_bands])
        self.emax = np.max([bandlike.band.emax for bandlike in self.rs.selected_bands])
        w = tools.LogLikelihood(self)
        return w
        
    def __enter__(self):
        """ supports the 'with' construction, guarantees that restore is called to restore the ROI
        example:
        -------
        with SourceFlux(roi, name) as sf:
            # use sf ...
        """
        return self
        
    def __exit__(self, type, value, traceback):
        self.restore()


        
class SED(object):
    """     
    generates a recarray (member rec) with fields:
        elow ehigh : energy range
        flux lflux uflux : energy flux (eV units) for peak, upper and lower limits
        ts : Test Statistic value for the band
        mflux delta_ts : model flux, and 2*(logl(flux)-logl(mflux), where logl is the log likelihood
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
        rec = makerec.RecArray('elow ehigh flux lflux uflux ts mflux delta_ts'.split())
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
                rec.append(xlo, xhi, w.maxl, lf, uf, w.TS(), mf, delta_ts)
            else:
                rec.append(xlo, xhi, 0, 0, w.upper_limit(), 0, mf, delta_ts )
            
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
    other kwargs passed to sed.Plot().__call__
    """
    sedfig_dir = kwargs.pop('sedfig_dir', None)
    if sedfig_dir is not None and not os.path.exists(sedfig_dir): os.mkdir(sedfig_dir)
    showts = kwargs.pop('showts', True)
    initw = roi.log_like()

    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    for source in sources:
        try:
            sf = SourceFlux(roi, source.name, )
            source.sedrec = SED(sf, merge=False).rec
            source.ts = roi.TS(source.name)
            if sedfig_dir is not None:
                annotation =(0.05,0.9, 'TS=%.0f'% source.ts) if showts else None 
                plotting.sed.Plot(source, gev_scale=True, energy_flux_unit='eV')\
                    ( galmap=source.skydir, outdir=sedfig_dir, 
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
    outdir = kwargs.pop('outdir', None)
    if outdir is not None:
        if not os.path.exists(outdir): os.mkdir(outdir)
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
    plx.axes.legend(loc='upper left', prop=dict(size=10) )
    if outdir is not None:
        plx.savefig(outdir)
    roi.set_model(old_model)
    return alt_fitqual
    
    
    