
"""
A module implementing an unbinned maximum likelihood fit of phases from some aperture
to a profile template.  The template is formed from an arbitrary number of components,
each represented by its own class.

LCPrimitives are combined to form a light curve (LCTemplate).  LCFitter then performs
a maximum likielihood fit to determine the light curve parameters.

LCFitter also allows fits to subsets of the phases for TOA calculation.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/lcfitters.py,v 1.6 2011/04/22 01:01:27 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>

"""

#TODO -- remove TOA stuff from LCFitter

import numpy as np

from lcprimitives import *
from operator import mul
from scipy.special import gamma
from scipy.optimize import fmin,leastsq
from uw.pulsar.stats import z2mw

SECSPERDAY = 86400.

def prim_io(template):
    """ Read files and build LCPrimitives. """

    def read_gaussian(toks):
        primitives = []
        for i,tok in enumerate(toks):
            if tok[0].startswith('phas'):
                g = LCGaussian()
                g.p[2] = float(tok[2])
                g.errors[2] = float(tok[4])
                primitives += [g]
            elif tok[0].startswith('fwhm'):
                g = primitives[-1]
                g.p[1] = float(tok[2])/2.3548200450309493      # kluge for now
                g.errors[1] = float(tok[4])/2.3548200450309493
            elif tok[0].startswith('ampl'):
                g = primitives[-1]
                g.p[0] = float(tok[2])
                g.errors[0] = float(tok[4])
        return primitives

    toks = [line.strip().split() for line in file(template) if len(line.strip()) > 0]
    if 'gauss' in toks[0]:     return read_gaussian(toks[1:])
    elif 'kernel' in toks[0]:  return [LCKernelDensity(input_file=toks[0])]
    elif 'fourier' in toks[0]: return [LCEmpiricalFourier(input_file=toks[0])]
    raise ValueError,'Template format not recognized!'


#==============================================================================#
class LCTemplate(object):
    """Manage a lightcurve template (collection of LCPrimitive objects).
   
       IMPORTANT: a constant background is assumed in the overall model, so there is
                 no need to furnish this separately.
    """

    def __init__(self,primitives=None,template=None):
        """primitives -- a list of LCPrimitive instances."""
        self.primitives = primitives
        if template is not None:
            self.primitives = prim_io(template)
        if self.primitives is None:
            raise ValueError,'No light curve components or template provided!'
        self.shift_mode = np.any([p.shift_mode for p in self.primitives])

    def set_parameters(self,p):
        start = 0
        for prim in self.primitives:
            n = int(prim.free.sum())
            prim.set_parameters(p[start:start+n])
            start += n

    def set_errors(self,errs):
        start = 0
        for prim in self.primitives:
            n = int(prim.free.sum())
            prim.errors = np.zeros_like(prim.p)
            prim.errors[prim.free] = errs[start:start+n]
            start += n

    def get_parameters(self):
        return np.concatenate( [prim.get_parameters() for prim in self.primitives] )

    def set_overall_phase(self,ph):
        """Put the peak of the first component at phase ph."""
        if self.shift_mode:
            self.primitives[0].p[0] = ph
            return
        shift = ph - self.primitives[0].get_location()
        for prim in self.primitives:
            new_location = (prim.get_location() + shift)%1
            prim.set_location(new_location)

    def get_location(self):
        return self.primitives[0].get_location()

    def norm(self):
        self.last_norm = sum( (prim.integrate() for prim in self.primitives) )
        return self.last_norm

    def max(self,resolution=0.01):
        return self(np.arange(0,1,resolution)).max()

    def __call__(self,phases,ignore_cache=False,suppress_bg=False):
        rval = np.zeros_like(phases)
        for prim in self.primitives:
            #rval += prim.cache_vals if (prim.cache and not ignore_cache) else prim(phases)
            rval += prim(phases)
        if suppress_bg: return rval/self.norm()
        else          : return (1.-self.norm()) + rval

    def gradient(self,phases):
        r = np.empty([len(self.get_parameters()),len(phases)])
        c = 0
        for prim in self.primitives:
            n = prim.free.sum()
            r[c:c+n,:] = prim.get_gradient(phases)
            c += n
        return r

    def __sort__(self):
        def cmp(p1,p2):
            if   p1.p[-1] <  p2.p[-1]: return -1
            elif p1.p[-1] == p2.p[-1]: return 0
            else: return 1
        self.primitives.sort(cmp=cmp)

    def __str__(self):
        self.__sort__()
        s1 = '\n'+'\n\n'.join( ['P%d -- '%(i+1)+str(prim) for i,prim in enumerate(self.primitives)] ) + '\n'
        if len(self.primitives)==2:
            p1,e1 = self.primitives[0].get_location(error=True) 
            p2,e2 = self.primitives[1].get_location(error=True)
            s1 +=  '\nDelta   : %.4f +\- %.4f\n'%(p2-p1,(e1**2+e2**2)**0.5)
        return s1

    def prof_string(self,outputfile=None):
        """ Return a string compatible with the format used by pygaussfit.
            Assume all primitives are gaussians."""
        rstrings = []
        dashes = '-'*25
        norm,errnorm = 0,0
      
        for nprim,prim in enumerate(self.primitives):
            phas = prim.get_location(error=True)
            fwhm = prim.get_width(error=True,fwhm=True)
            ampl = prim.get_norm(error=True)
            norm += ampl[0]
            errnorm += (ampl[1]**2)
            for st,va in zip(['phas','fwhm','ampl'],[phas,fwhm,ampl]):
                rstrings += ['%s%d = %.5f +/- %.5f'%(st,nprim+1,va[0],va[1])]
        const = 'const = %.5f +/- %.5f'%(1-norm,errnorm**0.5)
        rstring = [dashes] + [const] + rstrings + [dashes]
        if outputfile is not None:
            f = open(outputfile,'w')
            f.write('# gauss\n')
            for s in rstring: f.write(s+'\n')
        return '\n'.join(rstring)
       
    def random(self,n):
        """ Return n random variables drawn from the distribution given
            by this light curve template.

            Uses a mulitinomial to divvy the n phases up amongs the various
            components, which then each generate MC phases from their own
            distributions.

        """
        # multinomial implementation
        if len(self.primitives)==0: return np.random.rand(n)
        norms = [prim.get_norm() for prim in self.primitives]
        norms = np.append(norms,[1-sum(norms)])
        a = np.argsort(norms)[::-1]
        boundaries = np.cumsum(norms[a])
        components = np.searchsorted(boundaries,np.random.rand(n))
        rvals = np.empty(n)
        counter = 0
        for mapped_comp,comp in zip(a,np.arange(len(norms))):
            n = (components==comp).sum()
            if mapped_comp == len(norms)-1:
                rvals[counter:counter+n] = np.random.rand(n) 
            else:
                rvals[counter:counter+n] = self.primitives[mapped_comp].random(n)
            counter += n
        return rvals

    def swap_primitive(self,index,ptype='LCLorentzian'):
       """ Swap the specified primitive for a new one with the parameters
           that match the old one as closely as possible."""
       self.primitives[index] = convert_primitive(self.primitives[index],ptype)

#==============================================================================#

def LCFitter(template,phases,weights=None,**kwargs):
    if weights is None:
        return UnweightedLCFitter(template,phases,**kwargs)
    return WeightedLCFitter(template,phases,weights,**kwargs)

#==============================================================================#

class UnweightedLCFitter(object):
    """Perform the maximum likelihood fit template to the unbinned phases."""

    def __init__(self,template,phases,times=1):
        """
        template -- an instance of LCTemplate or a file with a pref-fit Gaussian template
        phases   -- a list or array of phase values
        """
        if type(template) == type(""):
            self.template = LCTemplate(gaussian_template=template)
        else: self.template = template
        self.phases   = np.asarray(phases)
        self.times    = np.asarray(times) #times (MET) for phases

        self.__hist_setup__()
        self.loglikelihood = self.unbinned_loglikelihood

    def __hist_setup__(self):
        """ Setup binning for a quick chi-squared fit."""
        nbins = min(max( int( len(self.phases)/50. ) , 32),101)
        hist = np.histogram(self.phases,bins=np.linspace(0,1,nbins))
        if len(hist[0])==nbins: raise ValueError,'Histogram too old!'
        x = ((hist[1][1:] + hist[1][:-1])/2.)[hist[0]>0]
        counts = hist[0][hist[0]>0]
        y    = counts / (x[1] - x[0]) / counts.sum()
        yerr = counts**0.5 / (x[1] - x[0]) / counts.sum()
        self.chistuff = x,y,yerr
        # now set up binning for binned likelihood
        nbins = 101
        hist = np.histogram(self.phases,bins=np.linspace(0,1,nbins))
        self.counts_centers = ((hist[1][1:] + hist[1][:-1])/2.)[hist[0]>0]
        self.counts = hist[0][hist[0]>0]

    def unbinned_loglikelihood(self,p,*args):
        if (not self.template.shift_mode and np.any(p < 0)):
         #guard against negative parameters
            return 2e20
        args[0].set_parameters(p)
        rvals = -np.log(self.template(self.phases)).sum()
        if np.isnan(rvals): return 2e20 # NB need to do better accounting of norm
        return rvals

    def binned_loglikelihood(self,p,*args):
        if not self.template.shift_mode and np.any(p<0):
            return 2e20
        args[0].set_parameters(p)
        return -(self.counts*np.log(self.template(self.counts_centers))).sum()
      
    def gradient(self,p,*args):
        args[0].set_parameters(p); t = self.template
        return -(t.gradient(self.phases)/t(self.phases)).sum(axis=1)

    def chi(self,p,*args):
        x,y,yerr = self.chistuff
        if not self.template.shift_mode and np.any(p < 0):
            return 2e100*np.ones_like(x)/len(x)
        args[0].set_parameters(p)
        chi = (self.template(x,ignore_cache=True) - y)/y**0.5
        if self.template.last_norm > 1:
            return 2e100*np.ones_like(x)/len(x)
        else: return chi

    def quick_fit(self):
        f = leastsq(self.chi,self.template.get_parameters(),args=(self.template))
         
    def fit(self,quick_fit_first=False, estimate_errors=True, unbinned=True):
      # an initial chi squared fit to find better seed values
        if quick_fit_first: self.quick_fit()

        self.loglikelihood = self.unbinned_loglikelihood if unbinned else self.binned_loglikelihood
        f = fmin(self.loglikelihood,self.template.get_parameters(),args=(self.template,),disp=0,ftol=1e-6,full_output=True)
        self.fitval = f[0]
        self.ll  = -f[1]
        if estimate_errors:
            self.__errors__()
            self.template.set_errors(np.diag(self.cov_matrix)**0.5)
            return self.fitval,np.diag(self.cov_matrix)**0.5
        else: return self.fitval

    def fit_cg(self):
        from scipy.optimize import fmin_cg
        fit = fmin_cg(self.loglikelihood,self.template.get_parameters(),fprime=self.gradient,args=(self.template,),full_output=1,disp=1)
        return fit

    def __errors__(self):
        from numpy.linalg import inv
        h1 = hessian(self.template,self.loglikelihood)
        try: 
            c1 = inv(h1)
            h2 = hessian(self.template,self.loglikelihood,delt=np.diag(c1)**0.5)
            c2 = inv(h2)
            if np.all(np.diag(c2)>0): self.cov_matrix = c2
            elif np.all(np.diag(c1)>0): self.cov_matrix = c1
            else: raise ValueError
        except:
            print 'Unable to invert hessian!'
            self.cov_matrix = np.zeros_like(h1)

    def __str__(self):
        if 'll' in self.__dict__.keys():
            return '\nLog Likelihood for fit: %.2f\n'%(self.ll) + str(self.template)
        return str(self.template)

    def write_template(self,outputfile='template.gauss'):
        s = self.template.prof_string(outputfile=outputfile)

    def plot(self,nbins=50,fignum=2):
        import pylab as pl
        dom = np.linspace(0,1,200)
        cod = self.template(dom,ignore_cache=True)
        pl.figure(fignum)
        pl.hist(self.phases,bins=np.linspace(0,1,nbins+1),histtype='step',ec='red',normed=True,lw=1)
        pl.plot(dom,cod,color='blue',lw=1)
        pl.xlabel('Phase')
        pl.grid(True)

#==============================================================================#

class WeightedLCFitter(UnweightedLCFitter):
    """Perform the maximum likelihood fit template to the unbinned phases."""

    def __init__(self,template,phases,weights,times=1):
        """
        template -- an instance of LCTemplate or a file with a pref-fit Gaussian template
        phases   -- a list or array of phase values
        """
        self.weights = weights
        super(WeightedLCFitter,self).__init__(template,phases,times=times)

    def unbinned_loglikelihood(self,p,*args):
        if not self.template.shift_mode and np.any(p < 0):
         #guard against negative parameters
            return 2e20
        args[0].set_parameters(p)
        return -np.log(1+self.weights*(self.template(self.phases)-1)).sum()

#==============================================================================#

class WeightedLCFitter_Approx(WeightedLCFitter):
    """Perform the maximum likelihood fit template to the unbinned phases."""

    def __init__(self,template,phases,weights,times=1,stride=1000,min_chi=15):
        """
        template -- an instance of LCTemplate or a file with a pref-fit Gaussian template
        phases   -- a list or array of phase values
        """
        super(WeightedLCFitter_Approx,self).__init__(template,phases,weights,times=times)
        self.stride = stride
        self.min_chi = min_chi
        #self.logl_accuracy = logl_accuracy
        #self.bes = np.linspace(0,1,nbins+1)
        #self.bcs = (self.bes[:-1]+self.bes[1:])/2
        a = np.argsort(weights)
        self.weights = weights[a]
        self.phases = phases[a]
        self.prep_approx()

    def fit(self,*args,**kwargs):
        #self.prep_approx()
        super(WeightedLCFitter_Approx,self).fit(*args,**kwargs)

    def prep_approx(self):
        min_n = self.weights.sum()
        for i in xrange(int(float(len(self.phases))/self.stride)+1):
            stat = z2mw(self.phases[:(i+1)*self.stride],self.weights[:(i+1)*self.stride])
            if stat[0] > self.min_chi:
                break
        nphotons = max(min_n,len(self.phases)-(i+1)*self.stride)
        self.sig_weights = self.weights[-nphotons+1:]
        self.sig_phases  = self.phases[-nphotons+1:]
        print 'using %d photons'%(nphotons)

    def unbinned_loglikelihood(self,p,*args):
        if not self.template.shift_mode and np.any(p < 0):
         #guard against negative parameters
            return 2e20
        #self.approx_ok()
        args[0].set_parameters(p)
        sig_logl = np.log(1+self.sig_weights*(self.template(self.sig_phases)-1)).sum()
        return -sig_logl
 
#==============================================================================#


def hessian(m,mf,*args,**kwargs):
   """Calculate the Hessian; mf is the minimizing function, m is the model,args additional arguments for mf."""
   p = m.get_parameters().copy()
   if 'delt' in kwargs.keys():
      delta = kwargs['delt']
   else: delta = [0.01]*len(p)
   hessian=np.zeros([len(p),len(p)])
   for i in xrange(len(p)):
      delt = delta[i]
      for j in xrange(i,len(p)): #Second partials by finite difference; could be done analytically in a future revision
         
         xhyh,xhyl,xlyh,xlyl=p.copy(),p.copy(),p.copy(),p.copy()
         xdelt = delt if p[i] >= 0 else -delt
         ydelt = delt if p[j] >= 0 else -delt

         xhyh[i]*=(1+xdelt)
         xhyh[j]*=(1+ydelt)

         xhyl[i]*=(1+xdelt)
         xhyl[j]*=(1-ydelt)

         xlyh[i]*=(1-xdelt)
         xlyh[j]*=(1+ydelt)

         xlyl[i]*=(1-xdelt)
         xlyl[j]*=(1-ydelt)

         hessian[i][j]=hessian[j][i]=(mf(xhyh,m,*args)-mf(xhyl,m,*args)-mf(xlyh,m,*args)+mf(xlyl,m,*args))/\
                                       (p[i]*p[j]*4*delt**2)

   mf(p,m,*args) #call likelihood with original values; this resets model and any other values that might be used later
   return hessian

#==============================================================================#

def get_gauss2(pulse_frac=1,x1=0.1,x2=0.55,ratio=1.5,width1=0.01,width2=0.02):
    """Return a two-gaussian template.  Convenience function."""
    n1,n2 = np.asarray([ratio,1.])*(pulse_frac/(1.+ratio))
    return LCTemplate(primitives=[LCGaussian(p=[n1,width1,x1]),
                                  LCGaussian(p=[n2,width2,x2])])

def get_gauss1(pulse_frac=1,x1=0.5,width1=0.01):
    """Return a one-gaussian template.  Convenience function."""
    return LCTemplate(primitives=[LCGaussian(p=[pulse_frac,width1,x1])])

def make_twoside_gaussian(one_side_gaussian):
    """ Make a two-sided gaussian with the same initial shape as the
        input one-sided gaussian."""
    g2 = LCGaussian2() 
    g1 = one_side_gaussian
    g2.p[0] = g1.p[0]
    g2.p[-1] = g1.p[-1]
    g2.p[1:3] = g1.p[1]
    return g2

def get_errors(template,total,n=100):
    from scipy.optimize import fmin
    ph0 = template.get_location()
    def logl(phi,*args):
        phases = args[0]
        template.set_overall_phase(phi%1)
        return -np.log(template(phases)).sum()
    errors = np.empty(n)
    fitvals = np.empty(n)
    errors_r = np.empty(n)
    delta = 0.01
    mean = 0
    for i in xrange(n):
        template.set_overall_phase(ph0)
        ph = template.random(total)
        results = fmin(logl,ph0,args=(ph,),full_output=1,disp=0)
        phi0,fopt = results[0],results[1]
        fitvals[i] = phi0
        mean += logl(phi0+delta,ph)-logl(phi0,ph)
        errors[i] = (logl(phi0+delta,ph)-fopt*2+logl(phi0-delta,ph))/delta**2
        my_delta = errors[i]**-0.5
        errors_r[i] = (logl(phi0+my_delta,ph)-fopt*2+logl(phi0-my_delta,ph))/my_delta**2
    print 'Mean: %.2f'%(mean/n)
    return fitvals-ph0,errors**-0.5,errors_r**-0.5

def make_err_plot(template,totals=[10,20,50,100,500],n=1000):
    import pylab as pl
    fvals = []; errs = []
    bins = np.arange(-5,5.1,0.25)
    for tot in totals:
        f,e = get_errors(template,tot,n=n)
        fvals += [f]; errs += [e]
        pl.hist(f/e,bins=np.arange(-5,5.1,0.5),histtype='step',normed=True,label='N = %d'%tot);
    g = lambda x: (np.pi*2)**-0.5*np.exp(-x**2/2)
    dom = np.linspace(-5,5,101)
    pl.plot(dom,g(dom),color='k')
    pl.legend()
    pl.axis([-5,5,0,0.5])


