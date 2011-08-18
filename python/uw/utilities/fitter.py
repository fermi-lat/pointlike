"""
Basic fitter utilities

Authors: Matthew Kerr, Toby Burnett
$Header$
"""

import numpy as np
from scipy import optimize #for fmin,fmin_powell,fmin_bfgs
from numpy import linalg  #for inv


class Minimizer(object):
    """ this is mostly extracted as is from uw.like.specfitter and turned into a utility
    """

    def __init__(self, fn, parameters=None, args=(), quiet=False):
        """ fn : function object
                note that it will be minimized, so should be negative of log likelihood
        """
        self.quiet = quiet
        self.par = parameters
        self.args = args
        self.fn = fn 
        npar = len(self.get_parameters())
        self.cov_matrix=np.zeros([npar,npar])
    
    def gradient(self,parameters,*args):
        """ access gradient if defined by the function
        """
        assert hasattr(self.fn, 'gradient'), 'Minimize: use_gradient set, but function did not define a gradient'
        return self.fn.gradient(parameters)

    def get_parameters(self):
        return self.fn.get_parameters() if self.par is None else self.par

    def set_parameters(self, par):
        if self.par is None:
            self.fn.set_parameters(par)
        else:
            self.par = par

    def get_free_errors(self):
        """Return the diagonal elements of the covariance matrix -- useful for step sizes in minimization, if known.
        """
        assert False, 'get_free_errors not implemented yet'
    
    def __call__(self,method='simplex', tolerance = 0.01, save_values = True, 
                      estimate_errors=True, error_for_steps=False,
                     use_gradient = False, gtol = 1e-1):
        """Maximize likelihood and estimate errors.
            method     -- ['simplex'] fitter; 'powell' or 'simplex' or 'minuit'
            tolerance -- (approximate) absolute tolerance 
        """
        if method not in ['simplex','powell','minuit']:
            raise Exception('Unknown fitting method for F.fit(): "%s"' % method)
        if method == 'minuit':
            from uw.utilities.minuit import Minuit
            temp_params = self.parameters()
            npars = self.parameters().shape[0]
            param_names = ['p%i'%i for i in xrange(npars)]
            
            if use_gradient:
                gradient         = self.gradient
                force_gradient = 1
            else:
                gradient         = None
                force_gradient = 0

            if error_for_steps:
                steps = self.get_free_errors()
                steps[steps<1e-6] = 0.04 # for models without error estimates, put in the defaults
                steps[steps > 1]  = 1     # probably don't want to step more than 100%...
                m = Minuit(self.fn,temp_params,up=.5,maxcalls=20000,tolerance=tolerance,printMode=-self.quiet,param_names=param_names,steps=steps)
            else:
                m = Minuit(self.fn,temp_params,up=.5,maxcalls=20000,tolerance=tolerance,printMode=-self.quiet,param_names=param_names)

            params,fval = m.minimize()

            if save_values:
                if estimate_errors == True:
                    self.__set_error_minuit(m,'HESSE')
                self.fn(params) # reset values to the ones found by minimization step
            self.fitvalue= fval
            return fval
        # scipy
        ll_0 = self.fn(self.get_parameters(), *self.args) 
        if ll_0==0: ll_0=1.0
        if use_gradient:
            f0 = optimize.fmin_bfgs(self.fn,self.get_parameters(),self.gradient,full_output=1,maxiter=500,gtol=gtol,disp=0)
            for i in xrange(10):
                f = self._save_bfgs = optimize.fmin_bfgs(self.fn,self.get_parameters(),self.gradient,
                        full_output=1,maxiter=500,gtol=gtol,disp=0)
                if abs(f0[1] - f[1]) < tolerance: break # note absolute tolerance
                if not self.quiet:
                    print 'Did not converge on first gradient iteration.  Trying again.'
                    print f0[1],f[1],abs(f0[1]-f[1])
                f0 = f
        else:
            minimizer  = optimize.fmin_powell if method == 'powell' else optimize.fmin
            f = minimizer(self.fn, self.get_parameters(),full_output=1,
                              maxiter=10000,maxfun=20000, ftol=0.01/abs(ll_0), disp=0 if self.quiet else 1)
        if not self.quiet: print 'Function value at minimum: %.8g'%f[1]
        self.set_parameters(f[0])
        self.fitvalue=f[1]
        if estimate_errors: 
            self.__set_error__(use_gradient)
        return f[1], f[0], np.sqrt(self.cov_matrix.diagonal()) if estimate_errors else None
    
    def correlations(self):
        """Return the linear correlation coefficients for the estimated covariance matrix."""
        sigmas = np.diag(self.cov_matrix)**0.5
        return self.cov_matrix / np.outer(sigmas,sigmas)

    def __set_error__(self,use_gradient=False):

        if use_gradient:
            hessian = Minimizer.mycov(self.gradient,self.get_parameters(),full_output=True)[1]
        else:
            hessian = Minimizer.hessian(self.fn, self.get_parameters())[0] 
        success = False
        # TODO -- check the return code
        full = True
        npar = len(self.get_parameters())
        try:
            if linalg.det(hessian)<=0:
                full=False
                #print 'singular? %s' %hessian
                hessian = hessian[1:,1:] #for now, first one
                
            if not self.quiet: print 'Attempting to invert full hessian...'
            self.cov_matrix =t = linalg.inv(hessian)
            if not full:
                self.cov_matrix = np.zeros((npar,npar))
                self.cov_matrix[1:,1:] = t
            if np.any(np.isnan(self.cov_matrix)):
                if not self.quiet: print 'Found NaN in covariance matrix!'
                raise Exception('Found NaN in covariance matrix!')
            if np.any(self.cov_matrix.diagonal()<0):
                print 'bad covariance matrix: diagonals:', self.cov_matrix.diagonal()
            success = True
        except linalg.LinAlgError:
            raise
            self.cov_matrix = np.zeros([npar,npar])
            success = False
        return success

    @staticmethod
    def hessian(mf, pars, *args):
        """Calculate the Hessian matrix using finite differences (adapted from specfitter.SpectralModelFitter.hessian)
        
         mf:   minimizing function
         pars: parameters at the minimum,
         args: additional arguments for mf.
        
        returns matrix, error code array
        """
        p  = pars.copy()
        npar = len(pars)
        deltas    = np.abs(0.01 * p) #initial guess
        hessian  = np.zeros([npar,npar])
        bad_mask = np.asarray([False] * npar)
        return_code = np.zeros(npar)

        l0 = mf(p, *args)

        #find good values with which to estimate the covariance matrix -- look at diagonal deviations
        #iterate until change in function consistent with ~1 sigma conditional error
        for i in xrange(npar):
            #print 'Working on parameter %d'%(i)
            h,l = p.copy(),p.copy()
            for j in xrange(10):
                h[:] = p[:]; l[:] = p[:];
                h[i] += deltas[i]
                l[i] -= deltas[i]

                delta_f_1 = mf(h, *args) - l0
                delta_f_2 = mf(l, *args) - l0
                delta_f = max(delta_f_1 + delta_f_2,0) #twice difference, really
                deltas[i] /= max(delta_f**0.5,0.33) # can change by half decade
                if delta_f < 5 and delta_f > 0.5: break

            if delta_f < 5e-3:
                # no constraint on parameter -- ignore it in further fittingor :
                bad_mask[i] = True
                return_code[i] = 1
            if (delta_f_1/delta_f_2 > 10 or delta_f_1/delta_f_2 < 1./10):
                # significant asymmetry in likelihood             
                bad_mask[i] = True
                return_code[i] = 2
            if (delta_f_2 < 5e-3 and delta_f_1 > 0.5):
                # not actually at maximum of likelihood -- upper limit condition
                bad_mask[i] = True
                return_code[i] = 3

        for i in xrange(npar):
            if bad_mask[i]:
                hessian[i,:] = 0 #no correlation?
                hessian[:,i] = 0
                continue
            for j in xrange(i,npar): #Second partials by finite difference
                
                xhyh,xhyl,xlyh,xlyl=p.copy(),p.copy(),p.copy(),p.copy()
                xdelt = deltas[i]
                ydelt = deltas[j]
                xhyh[i] += xdelt;  xhyh[j] += ydelt
                xhyl[i] += xdelt;  xhyl[j] -= ydelt
                xlyh[i] -= xdelt;  xlyh[j] += ydelt
                xlyl[i] -= xdelt;  xlyl[j] -= ydelt
                hessian[i][j]=hessian[j][i]=(mf(xhyh, *args)-mf(xhyl, *args)-mf(xlyh, *args)+mf(xlyl,*args))/\
                                                        (4*xdelt*ydelt)

        mf(p, *args) #call likelihood with original values; this resets model and any other values that might be used later
        return hessian,return_code

    @staticmethod
    def mycov(grad,par,full_output=False,init_step=0.04,min_step=1e-6,max_step=1,max_iters=5,target=0.5,min_func=1e-4,max_func=4):
        """Perform finite differences on the _analytic_ gradient provided by user to calculate hessian/covariance matrix.

        Positional args:
            grad                : a function to return a gradient
            par                 : vector of parameters (should be function minimum for covariance matrix calculation)

        Keyword args:

            full_output [False] : if True, return information about convergence, else just the covariance matrix
            init_step   [1e-3]  : initial step size (0.04 ~ 10% in log10 space); can be a scalar or vector
            min_step    [1e-6]  : the minimum step size to take in parameter space
            max_step    [1]     : the maximum step size to take in parameter sapce
            max_iters   [5]     : maximum number of iterations to attempt to converge on a good step size
            target      [0.5]   : the target change in the function value for step size
            min_func    [1e-4]  : the minimum allowable change in (abs) function value to accept for convergence
            max_func    [4]     : the maximum allowable change in (abs) function value to accept for convergence
        """

        nparams   = len(par)
        step_size = np.ones(nparams)*init_step
        step_size = np.maximum(step_size,min_step*1.1)
        step_size = np.minimum(step_size,max_step*0.9)
        hess      = np.zeros([nparams,nparams])
        min_flags = np.asarray([False]*nparams)
        max_flags = np.asarray([False]*nparams)

        def revised_step(delta_f,current_step,index):
            if (current_step == max_step):
                max_flags[i] = True
                return True,0
            elif (current_step == min_step):
                min_flags[i] = True
                return True,0
            else:
                adf = abs(delta_f)
                if adf < 1e-8:
                    # need to address a step size that results in a likelihood change that's too
                    # small compared to precision
                    pass
                    
                if (adf < min_func) or (adf > max_func):
                    new_step = current_step/(adf/target)
                    new_step = min(new_step,max_step)
                    new_step = max(new_step,min_step)
                    return False,new_step
                else:
                    return True,0
        
        iters = np.zeros(nparams)
        for i in xrange(nparams):
            converged = False
            for j in xrange(max_iters):
                iters[i] += 1
                di = step_size[i]
                par[i] += di
                g_up    = grad(par)
                par[i] -= 2*di
                g_dn    = grad(par)
                par[i] += di
                delta_f = (g_up - g_dn)[i]
                converged,new_step = revised_step(delta_f,di,i)
                #print 'Parameter %d -- Iteration %d -- Step size: %.2e -- delta: %.2e'%(i,j,di,delta_f)
                if converged: break
                else: step_size[i] = new_step
            hess[i,:] = (g_up - g_dn) / (2*di)  # central difference
            if not converged:
                print 'Warning: step size for parameter %d (%.2g) did not result in convergence.'%(i,di)
        
        try:
            cov = np.linalg.inv(hess)
        except:
            print 'Error inverting hessian.'
            #cov = np.zeros([nparams,nparams])
            raise Exception('Error inverting hessian')
        if full_output:
            return cov,hess,step_size,iters,min_flags,max_flags
        else:
            return cov

class AdaptFunc(object):
    """ adapt a function to make it a function of a subset of its parameters
    """
    def __init__(self, fn, par, select=[0]):
        """ par: array type
                default parameters to use
            fn: function of par: should be minimizable if use minimize
            select: list of free parameters
            TODO: use mask instead or optionally
        """
        self.fn=fn
        self.par = par[:]
        self.select = select
    def expand(self, x):
        p = self.par
        for i,j in enumerate(self.select): p[j]=x[i]
        return p
    def __call__(self, x):
        """ len of x must be number of selected parameters"""
        for i,j in enumerate(self.select): self.par[j]=x[i]
        return self.fn(self.par)
        
    def minimize(self, par0, **fit_kw):
        self.fit = Minimizer(self, par0, **fit_kw)
        c2, par, dpar = self.fit()
        return c2, par, dpar

def test(x0=1.1, pars=[1.0, 1.5], **kwargs):
    """ test with a parabola corresponding to a Gaussian with mean, sigma in pars
    """
    testf = lambda p: 1.+ 0.5*((p[0]-pars[0])/pars[1])**2
    class Func(object):
        def __init__(self, fn, pars):
            self.fn = fn
            self.pars = pars
        def __call__(self, pars):
            return self.fn(self.pars)
        def parameters(self): return self.pars
        def set_parameters(self,pars):
            self.pars = pars
    print 'input parameters:', pars
    m =  Minimizer(testf, [x0], )
    f = m(use_gradient=False)
    print 'solution at %.2f, +/- %.2f ' % (m.par[0], np.sqrt(m.cov_matrix.diagonal()[0]))
    return m, f
    
    