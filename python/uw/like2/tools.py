"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.1 2011/08/31 23:13:18 burnett Exp $

"""
import numpy as np
from uw.utilities import fitter
from scipy import optimize, interpolate
from skymaps import SkyDir
from uw.utilities import keyword_options


class ProjectToLinear(fitter.Projector):
    """ accept a function object and convert from log to linear
    """

    def __init__(self, fun, select, par=None, logpar=None, log=None):
        """ fun : minimizable function obejct to wrap
            select: list of variable numbers to use
            par :  initial parameter values or None
                if None, assume [1,1..]
            logpar : parameters to pass to function, or None
                if None, get from the fun.get_parameters()
            """
        if par is None: par=np.ones(len(select))
        super(ProjectToLinear, self).__init__( fun, select=select, par=logpar)
        self.logpar = self.par.copy()
        self.linpar = par
        self.log = log#logopen('log.txt', 'w')
        
    def __call__(self, pf):
        """ convert to normalization 
        """
        pfc = np.array(pf) #makes a copy
        #assert np.all(pfc>0), 'invalid flux value'
        if np.any(pfc<1e-3):
            pfa = pfc.copy()
            pfa[pfc<1e-3]=1e-3
            pfc[pfc<1e-6]=1e-6
            penalty = 1e2*np.sum(((pfa-pfc)/pfc)**2) 
        else: penalty=0
        pfc[pfc<1e-8] = 1e-8
        if self.log is not None: print >>self.log , pfc, penalty;  self.log.flush()
        return super(ProjectToLinear,self).__call__(self.logpar+np.log10(pfc)) + penalty

    def get_parameters(self):
        return 10**(super(ProjectToLinear,self).get_parameters()-self.logpar)
    def set_parameters(self,par):
        super(ProjectToLinear, self).set_parameters(self.logpar+np.log10(par))
        
    def gradient(self, pf):
        raise Exception( 'gradient not implemented' )

class FluxFitter(object):
    """ adapt a roistat.ROIstat object to fit only fluxes
    """
    def __init__(self, roi_stat, **kwargs):
        self.rs = roi_stat
        self.likelihood = kwargs.get('likeliood',True)
        self.sourcenames=[source.name for source in self.rs.sources if source.model.free[0]]
        print 'found %d free models in %s' % (len(self.sourcenames), roi_stat.name)
        select =[]
        for i,parname in enumerate(self.rs.parameter_names): 
            sname, p = parname.split('_')
            if p in ('Norm','Scale') and sname in self.sourcenames:
                select.append(i)
        assert len(select)==len(self.sourcenames), 'Missing model!'
        print 'selected parameters:', select
        
        self.setup_bands()
        self.setup_parameters(select)
        
    def setup_parameters(self,select):    
        # setup parameters. 
        self.projector = fitter.Projector(self.rs, select=select)
        self.logpar = self.projector.get_parameters().copy() # log parameters for reference
        self.par = np.ones_like(self.logpar) # inital are ones.
        self.linear_projector = ProjectToLinear(self.rs, select, self.par, self.logpar)
        self.dpar = np.zeros_like(self.par)
        self.fit_type = 'none yet'
        self.val = 0.

 
    def select_bands(self, selector=lambda b : True):
        self.rs.select_bands( selector)
        self.setup_bands()
    
    def setup_bands(self):
        bands = self.rs.selected_bands
        self.pixels = np.sum([band.pixels for band in bands])
        self.data = np.array([(band.data).sum() for band in bands])
        self.emin = np.min([band.band.emin for band in bands])
        self.emax = np.max([band.band.emax for band in bands])
        self.label = 'flux-only fits to ROI %s, %d bands, energies %0.f to %0.f' \
                % (self.rs.name, len(bands), self.emin, self.emax)
        
    def reset(self):
        self.projector.set_parameters(self.logpar)

    def __call__(self,par):
        return self.linear_projector(par)
        
    def fit(self, par0=None, quiet=True, **fit_kw):
        if par0 is None: par0=self.par

        try:
            if self.likelihood:
                self.fitter = fitter.Minimizer(self.linear_projector, par0, quiet=quiet)
                self.fit_type ='-logL'
            else:
                assert False, 'chisq not implemented yet'
                #self.fitter = fitter.Minimizer(lambda p: 0.5*self.bandstat.chisq(p), par0, quiet=True, **fit_kw)
                self.fit_type ='(chi squared)/2'
            self.val, self.par, self.dpar = self.fitter(**fit_kw)
        except Exception, e:
            print 'fit failed: %s, resetting' %e
            self.par = np.ones_like(self.logpar)
            self.dpar = np.zeros_like(self.par)
            raise
        self.reset() # restore original parameters, I hope
     
    def __str__(self):
        allcounts = np.array([[model.counts for model in band.freemodels]\
                            for band in self.rs.selected_bands])
        self.counts = allcounts.sum(axis=0)
        def fitval(i):
            if self.dpar[i]>0:
                return (' '+self.sourcenames[i], '%6.0f (%5.2f +/-%5.2f)'\
                    % (self.counts[i]*self.par[i], self.par[i], self.dpar[i]))
            return (' '+self.sourcenames[i], '%6.0f (%5.2f   ----  )'\
                    % (self.counts[i]*self.par[i], self.par[i], ))
        outtups = [
            (self.fit_type, '%6.1f'%self.val), 
            ('pixels', '%6d'%self.pixels), 
            ('data', '%6d'%self.data.sum()),
            ('model components', ''),
            ('  name            fit  measured/predicted', '')]\
            + map(fitval, range(len(self.par))) \
            + [('sum', '%6.0f' % np.sum(self.counts*self.par))]
            
        return self.label+'\n' + '\n'.join('%-15s %s' % t for t in  outtups  )

    def get_parameters(self):
        return self.par
    def analyze(self, index=0):
        t = fitter.Projector(self, select=[index], par=[self.par[index]])
        return lambda x : -t([x])
        

class SourceFlux(object):
    """ measure the energy-dependent flux for a given source
    """
    def __init__(self, rstat, source_name):
        """ rstat : ROIstat object
            source_name : name of one of the sources in the SourceList
        """
        self.rs = rstat
        names = [s.name for s in rstat.sources]
        try:
            self.source = rstat.sources[names.index(source_name)]
        except:
            raise Exception('source %s not found in ROI' %source_name)
        parname = source_name+'_Norm'
        try:
            self.pindex = list(rstat.parameter_names).index(parname)
        except:
            raise Exception('did not find parameter name, %s, for source flux'%parname)
        self.model=self.source.model
        self.energies = np.sort(list(set([ sm.band.e for sm in rstat.all_bands])))
        self.projector = ProjectToLinear(self.rs, [self.pindex])
        
    def fitall(self):
        """ fits for only this source, adjusting only the normalization, for each energy"""
        self.norm= [self.fitone(energy) for energy in self.energies]
    def fitone(self, energy):
        self.rs.select_bands(lambda b: b.e==energy)
        assert len(self.rs.selected_bands)==2, 'did not find two bands for energy %.1f' % energy
        mm = fitter.Minimizer(self.projector)
        #mm =self.rs.fit(select=[self.pindex], estimate_errors=False, quiet=True)
        return mm(estimate_errors=False)
    
def make_fits(self, **kwargs):
    self.rs.select_bands()
    projfn = fitter.Projector(self.rs, select=[self.pindex])
    par = projfn.get_parameters()
    print 'parameter: %.3f'% par
    print 'likelihood, gradient = %.0f, %.0f' % (projfn(par), projfn.gradient(par))
    fitfn  = fitter.Minimizer(projfn, **kwargs)
    a,b = fitfn(estimate_errors=False)
    print 'full fit:', b
    for energy in self.energies:

        print 'processing energy %.0f' % energy ,
        self.rs.select_bands(lambda b: b.e==energy)
        print 'likelihood, gradient = %.0f, %.0f' % (projfn(par), projfn.gradient(par))
        a,b,c = fitfn(estimate_errors=True)
        print '\t fit:', b,c
        
      
class NormFitter(object):
    """ code from ROIEnergyBand for reference"""
    def __init__(self, f):
        bad_fit = False
        self.m = PowerLaw(free=[True,False],e0=(self.emin*self.emax)**0.5) # fix index to 2

        self.fit = optimize.fmin(f,self.m.get_parameters(),disp=0,full_output=1,args=(self.m,which))

        def upper_limit():
            ### warning: this code depends on the log10 representation of the flux

            flux_copy = self.m[0]
            zp          = self.f(N.asarray([-20]),self.m,which)

            def f95(parameters):
                return abs(f(parameters,self.m,which) - zp - 3.0)
            
            # for some reason, can't get fsolve to work here.  good ol' fmin to the rescue
            self.uflux = 10**optimize.fmin(f95,N.asarray([-11.75]),disp=0)[0]
            self.lflux = None
            self.flux  = None

            self.m[0] = flux_copy

        # if flux below a certain level, set an upper limit
        if self.m[0] < 1e-20:
            bad_fit = True
            upper_limit()

        else:
            try:
                #self.m.set_cov_matrix([[self.normUncertainty,0],[0,0]])
                err = self.normUncertainty(which=which)
            except:
                bad_fit = True
                err = 0 

            self.flux  = self.m[0] 
            self.uflux = self.flux*(1 + err)
            self.lflux = max(self.flux*(1 - err),1e-30)

        if saveto is not None:
            for b in self.bands: b.__dict__[saveto] = (b.expected(self.m)*b.er[which] if not bad_fit else -1)

        if bad_fit:
            self.ts = 0
        else:
            null_ll = sum( (b.bandLikelihood([0],which) for b in self.bands) )
            alt_ll  = sum( (b.bandLikelihood([b.expected(self.m)*b.er[which]],which) for b in self.bands) )
            self.ts = 2*(null_ll - alt_ll)
            
        return self.ts
        
    def normUncertainty(self,which=0):
        # this is the uncertainty derived by taking the second derivative of the log likelihood
        # wrt the normalization parameter; it is fractional (delta_v/v)
        # this formulation is specifically for the flux density of a power law, which has such nice properties
        tot = 0
        for b in self.bands:
            if not b.has_pixels: continue
            my_pix_counts = b.ps_pix_counts[:,which]*b.expected(self.m)*b.er[which]
            all_pix_counts= b.bg_all_pix_counts + b.ps_all_pix_counts - b.ps_pix_counts[:,which]*b.ps_counts[which] + my_pix_counts
            tot += (b.pix_counts * (my_pix_counts/all_pix_counts)**2).sum()
        return tot**-0.5
 

        
class Localization(object):
    """ manage localization of a source
    Implements a minimization interface
    
    """

    def __init__(self, roistat, source_name):
        """ roistat : an ROIstat object
            source_name : string
                the name of a source
        """
        self.rs = roistat
        self.source_mask = np.array([source.name==source_name for source in roistat.sources])
        if sum(self.source_mask)!=1:
            raise Exception('Localization: source %s not found'%source_name)
        self.rs.initialize(self.source_mask)
        self.rs.update(True)
        self.maxlike=self.rs.log_like()
        self.source = np.array(roistat.sources)[self.source_mask][0]
        self.skydir = self.source.skydir
   
    def log_like(self, skydir):
        """ return log likelihood at the given position"""
        self.source.skydir =skydir
        self.rs.update(True)
        return self.rs.log_like()
   
    def TS(self, skydir):
        """ return the TS at given position, or 
            2x the log(likelihood ratio) from the nominal position
        """
        return 2*(self.log_like(skydir)-self.maxlike)

    def get_parameters(self):
        return np.array([self.source.skydir.ra(), self.source.skydir.dec()])
    
    def set_parameters(self, par):
        self.skydir = SkyDir(par[0],par[1])
        self.source.skydir = self.skydir
        
    def __call__(self, par):
        return -self.TS(SkyDir(par[0],par[1]))
    
    def reset(self):
        """ restore modifications to the ROIstat
        """
        self.source.skydir=self.skydir
        self.rs.update(True)
        self.rs.initialize()
        
        
class PoissLogLikelihood(object):
    defaults = (
        ('atol',0.5,'absolute maximum difference allowed between poly and interp'),
        #('remove_offenders',False,'try to remove bad likelihood points'),
        ('tstart',0,'optional start time'),
        ('tstop',0,'optional stop time'),
    )
    
    @keyword_options.decorate(defaults)
    def _poiss(self,dom):
        sp,e,b = self.poissp
        b = abs(b) # this is a bit of a kluge to keep log argument positive
        r = e*(dom+b)
        r_peak = e*(sp+b)
        if sp > 0:
            const = r_peak*np.log(r_peak) - r_peak
        else:
           # const = r_peak*np.log(r_peak) - r_peak
            t = e*b
            const = r_peak*np.log(t) - t
        f = r_peak*np.log(r) - r
        return f - const

    def _pre_process(self):
        diffs = self.cod0[1:]-self.cod0[:-1]
        signs = diffs > 0
        divider = len(self.dom0)/2 - 1
        mask = np.asarray([True]*len(self.cod0))
        for i in xrange(len(diffs)-1):
            if (signs[i+1]!=signs[i]) and (signs[i-1]!=signs[i]):
                if i < divider:
                    mask[:i] = False
                else: mask[i:] = False
        if mask.sum()!=len(mask):
            print 'Removed %d entries'%(len(mask)-mask.sum())
        self.premodcod0 = self.cod0.copy()
        self.premoddom0 = self.dom0.copy()
        self.cod0 = self.cod0[mask]
        self.dom0 = self.dom0[mask]
            
    def _do_fit(self,exp_seed=3e11,bg_fac=1):
        dom0 = self.dom0; cod0 = self.cod0
        amax = np.argmax(cod0)
        smax = dom0[amax]
        # use just three (best?) points
        #dom0 = np.asarray([dom0[1],dom0[amax],dom0[-1]])
        #cod0 = np.asarray([cod0[1],cod0[amax],cod0[-1]])
        def fitfunc(p):
            self.poissp = p
            return self._poiss(dom0)-cod0
        if amax==0:
            smax = -dom0[1]
            #smax = 0.
        #    e = -(cod0[-1]-cod0[0])/(dom0[-1]-dom0[0])
        #    return leastsq(fitfunc,[1e-13,e,1e-13])[0]
        return optimize.leastsq(fitfunc,[smax,exp_seed,dom0[-1]/bg_fac],maxfev=3600)[0]
        
    def __init__(self,dom0,cod0,**kwargs):
        """
        dom0 -- the sampled domain of the log likelihood surface
        cod0 -- the log likelihood as a function of dom0
        """
        keyword_options.process(self,kwargs)
        self.dom0 = dom0; self.cod0 = cod0-cod0.max();
        self._pre_process()
        self.interp = interpolate.interp1d(self.dom0,self.cod0,bounds_error=False,fill_value=self.cod0[-1])
        for seed in [1e14,3e10,1e9]:
            for bg_fac in [1.,10.,0.1]:
                self.poissp = self._do_fit(exp_seed=seed,bg_fac=bg_fac)
                self.ok = self._check_agreement()
                if self.ok: break
            if self.ok: break # use return instead?
        if not self.ok:
            print 'WARNING fit and interpolation are significantly different (%.2f)'%(self.agreement)
    def _check_agreement(self):
        v1 = self(self.dom0,use_fit=True); v2 = self(self.dom0,use_fit=False)
        self.agreement = max_diff = (np.abs(v2-v1)).max()
        return max_diff < self.atol

    def __call__(self,dom,use_fit=True):
        """ Return the likelihood as a function of the domain."""
        if use_fit: return self._poiss(dom)
        return self.interp(dom)

    def show(self,use_fit=True):
        import pylab as pl
        pl.clf()
        #pl.plot(self.dom0,self(self.dom0,use_fit=use_fit))
        dom = np.linspace(0,self.dom0[-1],100)
        pl.plot(dom,self(dom,use_fit=use_fit))
        pl.plot(self.dom0,self.cod0,marker='o',ls=' ',color='red')
        pl.axvline(self.find_max(use_fit=use_fit))

    def find_max(self,use_fit=True):
        if use_fit: return max(0,self.poissp[0])
        return self.dom0[np.argmax(self.cod0)]

    def ts(self,use_fit=True):
        ll1 = self(self.find_max(use_fit=use_fit))
        ll0 = self(0,use_fit=use_fit)#.cod0[0] # is this best?
        return 2*(ll1-ll0)

    def error_bars(self,delta_logl=0.5):
        xmax = self.find_max()
        ymax = self(xmax,use_fit=True)
        def fitfunc(x):
            return self(x)+delta_logl-ymax
        r_lo = optimize.bisect(fitfunc,0,xmax)
        r_hi = optimize.bisect(fitfunc,xmax,10*xmax)
        roots = np.asarray([r_lo,r_hi])
        mask = np.isreal(roots)
        if (mask.sum() != 2) or (r_hi==r_lo):
            print 'Could not find 2 roots!'
            return None
        return xmax-r_lo,r_hi-xmax
        