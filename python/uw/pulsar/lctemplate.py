"""
A module implementing a mixture model of LCPrimitives to form a
normalized template representing directional data.
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/lctemplate.py,v 1.2 2012/03/11 23:00:52 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>

"""

import numpy as np
from lcnorm import NormAngles
from lcprimitives import *

class LCTemplate(object):
    """Manage a lightcurve template (collection of LCPrimitive objects).
   
       IMPORTANT: a constant background is assumed in the overall model, 
       so there is no need to furnish this separately.

       The template is such that 
    """

    #def __init__(self,primitives=None,template=None,norms=None):
    def __init__(self,primitives,norms):
        """primitives -- a list of LCPrimitive instances."""
        self.primitives = primitives
        #if template is not None:
        #    self.primitives = prim_io(template)
        #if self.primitives is None:
        #    raise ValueError,'No light curve components or template provided!'
        self.shift_mode = np.any([p.shift_mode for p in self.primitives])
        # new stuff for normalization
        #self.free = np.asarray([True]*len(self.primitives)])
        if norms is None: norms = np.ones(len(primitives))/len(primitives)
        self.norms = NormAngles(norms)

    def __getitem__(self,index): return self.primitives[index]
    def __setitem__(self,index,value): self.primitives[index]=value

    def set_parameters(self,p):
        start = 0
        for prim in self.primitives:
            n = int(prim.free.sum())
            prim.set_parameters(p[start:start+n])
            start += n
        self.norms.set_parameters(p[start:start+self.norms.free.sum()])

    def set_errors(self,errs):
        start = 0
        for prim in self.primitives:
            n = int(prim.free.sum())
            prim.errors = np.zeros_like(prim.p)
            prim.errors[prim.free] = errs[start:start+n]
            start += n

    def get_parameters(self):
        return np.append(np.concatenate( [prim.get_parameters() for prim in self.primitives]) , self.norms.get_parameters())

    def get_bounds(self):
        b1 = np.concatenate([prim.get_bounds() for prim in self.primitives])
        b2 = self.norms.get_bounds()
        return np.concatenate((b1,b2)).tolist()

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
        #self.last_norm = sum( (prim.integrate() for prim in self.primitives) )
        #return self.last_norm
        return self.norms.get_total()
        

    def integrate(self,phi1,phi2, suppress_bg=False):
        #norm = self.norm()
        norms = self.norms()
        t = norms.sum()
        dphi = (phi2-phi1)
        if suppress_bg: return sum( (n*prim.integrate(phi1,phi2) for n,prim in zip(norms,self.primitives)) )/t

        return (1-t)*dphi + sum( (n*prim.integrate(phi1,phi2) for n,prim in zip(norms,self.primitives)) )

    def max(self,resolution=0.01):
        return self(np.arange(0,1,resolution)).max()

    def __call__(self,phases,suppress_bg=False):
        #n = self.norm()
        norms = self.norms()
        rval = np.zeros_like(phases)
        for n,prim in zip(norms,self.primitives):
            rval += n*prim(phases)
        if suppress_bg: return rval/norms.sum()
        return (1-norms.sum()) + rval

    def gradient(self,phases):
        r = np.empty([len(self.get_parameters()),len(phases)])
        c = 0
        norms = self.norms()
        prim_terms = np.empty([len(phases),len(self.primitives)])
        for i,(norm,prim) in enumerate(zip(norms,self.primitives)):
            n = prim.free.sum()
            r[c:c+n,:] = norm*prim.get_gradient(phases)
            c += n
            prim_terms[:,i] = prim(phases)-1
        m = self.norms.get_grads()
        for j in xrange(len(norms)):
            if not self.norms.free[j]: continue
            r[c,:] = (prim_terms*m[:,j]).sum(axis=1)
            c += 1
        return r

    def approx_gradient(self,phases,eps=1e-5):
        orig_p = self.get_parameters().copy()
        g = np.zeros([len(orig_p),len(phases)])
        weights = np.asarray([-1,8,-8,1])/(12*eps)

        def do_step(which,eps):
            p0 = orig_p.copy()
            p0[which] += eps
            self.set_parameters(p0)
            return self(phases)

        for i in xrange(len(orig_p)):
            # use a 4th-order central difference scheme
            for j,w in zip([2,1,-1,-2],weights):
                g[i,:] += w*do_step(i,j*eps)

        self.set_parameters(orig_p)
        return g

    def check_gradient(self,tol=1e-6,quiet=False):
        """ Test gradient function with a set of MC photons."""
        ph = self.random(1000)
        g1 = self.gradient(ph)
        g2 = self.approx_gradient(ph)
        anyfail = False
        for i in xrange(g1.shape[0]):
            d1 = np.abs(g1[i]-g2[i])
            d2 = d1/g1
            fail = np.any((d1>tol) | (d2>tol))
            if not quiet:
                pass_string = 'FAILED' if fail else 'passed'
                print '%d (%s) %.3g (abs) %.3g (frac)'%(i,pass_string,d1.max(),d2.max())
            anyfail = anyfail or fail
        return not anyfail

    def delta(self,index=None):
        """ Return radio lag -- reckoned by default as the posittion of the            first peak following phase 0."""
        if (index is not None) and (index<=(len(self.primitives))):
            return self[index].get_location(error=True)
        return self.Delta(delta=True)

    def Delta(self,delta=False):
        """ Report peak separation -- reckoned by default as the distance
            between the first and final component locations.
            
            delta [False] -- if True, return the first peak position"""
        if len(self.primitives)==1: return -1
        prim0,prim1 = self[0],self[-1]
        for p in self.primitives:
            if p.get_location() < prim0.get_location(): prim0 = p
            if p.get_location() > prim1.get_location(): prim1 = p
        p1,e1 = prim0.get_location(error=True) 
        p2,e2 = prim1.get_location(error=True)
        if delta: return p1,e1
        return (p2-p1,(e1**2+e2**2)**0.5)

    def _sorted_prims(self):
        def cmp(p1,p2):
            if   p1.p[-1] <  p2.p[-1]: return -1
            elif p1.p[-1] == p2.p[-1]: return 0
            else: return 1
        return sorted(self.primitives,cmp=cmp)

    def __str__(self):
        #prims = self._sorted_prims()
        prims = self.primitives
        def norm_string(i):
            fstring = '' if self.norms.free[i] else ' [FIXED]'
            return 'P%d : %.4f +\- %.4f%s'%(i+1,self.norms()[i],0,fstring)
        s0 = '\nMixture Amplitudes\n------------------\n'+\
             '\n'.join([norm_string(i) for i in xrange(len(prims))])+\
             '\nDC : %.4f +\- %.4f'%(1-self.norms.get_total(),0)
        s1 = '\n\n'+'\n\n'.join( ['P%d -- '%(i+1)+str(prim) for i,prim in enumerate(prims)] ) + '\n'
        s1 +=  '\ndelta   : %.4f +\- %.4f'%self.delta()
        s1 +=  '\nDelta   : %.4f +\- %.4f'%self.Delta()
        return s0+s1

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
       
    def random(self,n,weights=None):
        """ Return n pseudo-random variables drawn from the distribution 
            given by this light curve template.

            Uses a mulitinomial to divvy the n phases up amongs the various
            components, which then each generate MC phases from their own
            distributions.

            If a vector of weights is provided, the weight in interpreted
            as the probability that a photon comes from the template.  A
            random check is done according to this probability to
            determine whether to draw the photon from the template or from
            a uniform distribution.
        """
        rvals = np.empty(n)
        if weights is not None:
            if len(weights) != n:
                raise Exception('Provided weight vector does not provide a weight for each photon.')
            t = np.random.rand(n)
            m = t <= weights
            t_indices = np.arange(n)[m]
            n = m.sum()
            rvals[~m] = np.random.rand(len(m)-n)
        else:
            t_indices = np.arange(n)

        # multinomial implementation -- draw from the template components
        if len(self.primitives)==0: return np.random.rand(n)
        #norms = [prim.get_norm() for prim in self.primitives]
        norms = self.norms()
        norms = np.append(norms,[1-sum(norms)])
        a = np.argsort(norms)[::-1]
        boundaries = np.cumsum(norms[a])
        components = np.searchsorted(boundaries,np.random.rand(n))
        counter = 0
        for mapped_comp,comp in zip(a,np.arange(len(norms))):
            n = (components==comp).sum()
            if mapped_comp == len(norms)-1:
                rvals[t_indices[counter:counter+n]] = np.random.rand(n) 
            else:
                rvals[t_indices[counter:counter+n]] = self.primitives[mapped_comp].random(n)
            counter += n
        return rvals

    def swap_primitive(self,index,ptype=LCLorentzian):
       """ Swap the specified primitive for a new one with the parameters
           that match the old one as closely as possible."""
       self.primitives[index] = convert_primitive(self.primitives[index],ptype)

def get_gauss2(pulse_frac=1,x1=0.1,x2=0.55,ratio=1.5,width1=0.01,width2=0.02,lorentzian=False,bridge_frac=0,skew=False):
    """Return a two-gaussian template.  Convenience function."""
    n1,n2 = np.asarray([ratio,1.])*(1-bridge_frac)*(pulse_frac/(1.+ratio))
    if skew:
        prim = LCLorentzian2 if lorentzian else LCGaussian2
        p1,p2 = [width1,width1*(1+skew),x1],[width2*(1+skew),width2,x2]
    else:
        if lorentzian:
            prim = LCLorentzian; width1 *= (2*np.pi); width2 *= (2*np.pi)
        else:
            prim = LCGaussian
        p1,p2 = [width1,x1],[width2,x2]
    if bridge_frac > 0:
        nb = bridge_frac*pulse_frac
        b = LCGaussian(p=[0.1,(x2+x1)/2])
        return LCTemplate([prim(p=p1),b,prim(p=p2)],[n1,nb,n2])
    return LCTemplate([prim(p=p1),prim(p=p2)],[n1,n2])

def get_gauss1(pulse_frac=1,x1=0.5,width1=0.01):
    """Return a one-gaussian template.  Convenience function."""
    return LCTemplate(primitives=[LCGaussian(p=[pulse_frac,width1,x1])])

def get_2pb(pulse_frac=0.5,lorentzian=False):
    """ Convenience function to get a 2 Lorentzian + Gaussian bridge template."""
    prim = LCLorentzian if lorentzian else LCGaussian
    p1 = prim(p=[0.3*pulse_frac,0.03,0.1])
    b = LCGaussian(p=[0.4*pulse_frac,0.15,0.3])
    p2 = prim(p=[0.3*pulse_frac,0.03,0.55])
    return LCTemplate(primitives=[p1,b,p2])

def make_twoside_gaussian(one_side_gaussian):
    """ Make a two-sided gaussian with the same initial shape as the
        input one-sided gaussian."""
    g2 = LCGaussian2() 
    g1 = one_side_gaussian
    g2.p[0] = g1.p[0]
    g2.p[-1] = g1.p[-1]
    g2.p[1:3] = g1.p[1]
    return g2
