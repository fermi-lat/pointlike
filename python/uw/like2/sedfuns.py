"""
Tools for ROI analysis - Spectral Energy Distribution functions

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sedfuns.py,v 1.30 2013/12/22 20:05:21 burnett Exp $

"""
import os, pickle
import numpy as np
from uw.utilities import ( keyword_options)
from . import ( plotting, tools, loglikelihood)

       
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
        self.func = self.rs.energy_flux_view(source_name)
        self.source_name = source_name
        # make a list of energies with data; only have info if there is data in the ROI
        hasdata = np.array([b.pixels>0 for b in self.rs])
        self.energies = self.rs.energies[hasdata[0::2] | hasdata[1::2]] ## note assumption about energies 
    
    def full(self):
        try:
            self.full_poiss = self.select(None).poiss
        except Exception, msg:
            print 'Failed poisson fit to source %s: "%s"' % (source_name, msg)
            raise

    def __repr__(self):
        return '%s.%s : %d bands selected for source, %s energy range %.0f-%.0f'% (
                    self.__module__, self.__class__.__name__,len(self.rs.selected), self.source_name, 
                    self.rs.emin, self.rs.emax)
    
    def select(self, index, event_type=None, poisson_tolerance=0.10):
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
            self.func.set_energy(None)# =func = self.rs.energy_flux_view(self.source_name)
        else:
            self.rs.select(index, event_type)
            energies = self.rs.energies
            assert len(energies)==1
            energy = self.rs.energies[0]
            self.func.set_energy(energy)
            assert self.func(0) != self.func(1), 'Function not variable? energy %.0f' % energy
        pf = loglikelihood.PoissonFitter(self.func, tol=poisson_tolerance)
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
                print msg
                pp.append(None)
                
        self.restore()
        return np.array(pp)
        
    def sed_rec(self, event_type=None, tol=0.1):
        """ return a numpy.recarray with values for each band
           elow ehigh       -- energy limits
           flux lflux uflux -- flux at max, upper and lower 1-sigma
           ts               -- Test Statistic for the band
           mflux            -- Flux predicted by the model
           delta_ts         -- TS difference, fit-model
           pull             -- signed square root of delta_ts
        """
        names = 'elow ehigh flux lflux uflux ts mflux delta_ts pull maxdev'.split()
        rec = tools.RecArray(names, dtype=dict(names=names, formats=['>f4']*len(names)) )
        for i,energy in enumerate(self.energies):
            pf = self.select(i, event_type=event_type, poisson_tolerance=tol)
            w = pf.poiss
            err = pf.maxdev
            xlo,xhi = self.rs.emin,self.rs.emax
            lf,uf = w.errors
            maxl  = w.flux
            mf    = self.func.eflux
            delta_ts = 2.*(self(maxl) - self(mf) )
            if lf>0 :
                pull = np.sign(maxl-mf) * np.sqrt(max(0, delta_ts))
                rec.append(xlo, xhi, maxl, lf, uf, w.ts, mf, delta_ts, pull, err)
            else:
                pull = -np.sqrt(max(0, delta_ts))
                rec.append(xlo, xhi, 0, 0, w.cdfcinv(0.05), 0, mf, delta_ts, pull, err )
            
        self.restore()
        return rec()

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
           
        fig, axx = plt.subplots(4,4, figsize=(12,12), sharex=True, sharey=True)
        for i, ax in enumerate(axx.flatten()):
            if i >= len(self.energies):
                ax.set_visible(False)
                continue
            pf = self.select(i)
            pf.plot(ax)
            ax.set_title('%s@ %d MeV' %( self.source_name, int(self.func.energy),), size=10)
        self.restore()
        return fig
        

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
    if sedfig_dir is not None and sedfig_dir[0]=='$':
        sedfig_dir = os.path.expandvars(sedfig_dir)
    ndf = kwargs.pop('ndf', 10) 
    if sedfig_dir is not None and not os.path.exists(sedfig_dir): os.mkdir(sedfig_dir)
    showts = kwargs.pop('showts', True)
    poisson_tolerance = kwargs.pop('poisson_tolerance', 0.25)
    initw = roi.log_like()

    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    for source in sources:
        with SED(roi, source.name, ) as sf:
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


