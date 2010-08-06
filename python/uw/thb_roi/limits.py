"""
Study limits

"""
from uw.thb_roi import variation 
from setup import *
import pylab as plt
import numpy as np
from uw.utilities import tsmap_grid, makerec
from scipy import stats
import types, glob
def getlikes(months=5, source=1):
    md = MonthlyData()
    lc = [pipeline.LightCurve(outdir,md(i)) for i in range(months)]
    roi = [lc[i].roi(source) for i in range(len(lc))]
    ww = [LogLikelihood(r, 1e-11) for r in roi]
    return ww

def plotlikes(ww):
    plt.figure(10); plt.clf(); ax = plt.gca()
    for i,w in enumerate(ww): 
        w.plot(axes=ax, label='month %d'% (i+1) )
        sig = np.exp(-0.5)
        a, b = w.errors()
        ax.plot( [a,b], [sig,sig], '-k')
        ax.plot( [w.maxl, w.maxl], [0.975,1.025], '-k')
        f95 = w.upper_limit()
        ax.plot( [f95,f95], [0, w(f95)], '-k')
    ax.legend()
    ax.grid(True)
    ax.set_xlabel('Flux density at pivot')
    plt.title('flux likelihoods')

def point(axis, x, yvals, xerr=0, TS=None, arrow_size=0.4, **kwargs):
    """
        x, x-value
        yvals - a tuple (y, yl,yu, y95)
    """
    y, yl,yu,y95 = yvals[:4]
    a,b = axis.get_ylim(); 
    fmt =kwargs.pop('fmt','o')
    if yl>2*a or TS>9: 
        # lower limit is clearly above lower limit: plot the point
        axis.plot([x], [y], fmt, **kwargs)
        axis.plot([x-xerr,x+xerr], [y,y], '-', **kwargs)
        axis.plot([x,x], [yl, yu], '-', **kwargs)
        #axis.errorbar( [x], [y], xerr =xerr, yerr=[[yu-y],[y-yl]] if yu>y else 0, fmt=fmt, capsize=0,**kwargs)
    elif y95>0:
        axis.plot( [x-xerr,x+xerr], [y95,y95], **kwargs)
        dx = arrow_size
        dy = arrow_size*(b-a) 
        if axis.get_yscale()=='linear':
            ya,yb,yc = y95, y95-dy, y95-2*dy
        else:
            ya, yb, yc = y95, y95/1.2, y95/1.2**2
        kwargs.pop('ms', None)
        axis.fill( [x,  x, x+dx, x,  x-dx, x], 
                   [ya, yb,yb,   yc, yb,  yb], **kwargs)
    else:
        axis.plot( [x-xerr,x+xerr], [y95,y95], '*', **kwargs)
  
class LightCurve(object):
    """ manage the light curve info contained in the source pickle file
    """
    def __init__(self, pk, N=None, **kwargs):
        if type(pk)==types.StringType:
            raise Exception, 'String not supported'
        self.pk =pk
        self.norm =pk['src_par'][0]
        self.pivot = pk['pivot_energy']
        self.months =pk['months']
        self.name = pk['name']
        # allow for negative keys, used for sum
        self.keys = np.sort(np.array([k for k in self.months.keys() if k>=0]))
        self.N = len(self.keys) if N is None else N
        self.like = np.asarray([self.months[i]['likelihood'] for i in self.keys[:self.N]] )

        # compute check on assumption that monthly sums correspond to full fit
        # sigma should be 1.0, delta 0
        b,a,c  = [self.delta_loglike(i) for i in range(3)]
        self.sigma = np.sqrt(a+c-2*b) # corrected sigma (in original units)
        self.delta = (c-a)/self.sigma # corrected peak, in sigma units
        
    def delta_loglike(self, i=0):
        """ determine the log likelihood difference for the constant assumption
            i [0 default] has meaning: 0 max, 1, -sigma, 2 +sigma
            Expect that the 1 sigma points will be 0.5 less
        """
        return -sum([np.log(v[5+i]) for v in self.like]) 
     
    def survival_prob(self):
        """ calcualte the "survival probability, assuming distributed as chi**2 with N-1 deg of f"""
        return stats.chi2.sf(2*self.delta_loglike(),self.N-1 )
        
    def plot(self, *pars, **kwargs):
        plotmon(self.pk, *pars, **kwargs)
  
    def print_line(self,i, q):
        z = np.array(q[:4])/self.norm
        t = (i+1,)+tuple(z)+(q[4], -np.log(q[5]),)
        return ('\n%8d'+ (len(t)-1)*'  %8.2f') % t

    def __str__(self, **kwargs):
        s =  '\nLight curve for %s, flux=%.2e at %.2f GeV\n' % (self.name, self.norm, self.pivot/1e3)
        s += ('   '+7*'%-10s') % tuple('month peak low high 95% TS deltaL'.split()) 
        for i,q in enumerate(self.like):
            s += self.print_line(i,q)
        s += '\n   Total'+ '%50.2f%10.2f' % (sum(self.like[:,4]), -sum(np.log(self.like[:,5])))
        s += '\nlog likelihood difference =%.1f for %d months: prob=%.1e\n'\
            %( self.delta_loglike(), len(self.like), self.survival_prob())
            
        return s
            
def plotmon(p, fignum=11, axes=None, ymax=10, saveto=None, **kwargs):
    if type(p)==types.StringType:
        p = pickle.load(open(os.path.join(outdir, 'pickle', '%s.pickle' % p.replace(' ', '_').replace('+','p'))))
    norm =p['src_par'][0]
    pivot = p['pivot_energy']
    months =p['months']
    name = p['name']
    q = [months[i]['likelihood'] for i in months.keys() if i>=0] 
    oldlw = plt.rcParams['axes.linewidth']
    plt.rcParams['axes.linewidth'] = 2
    if axes is None:
        fig=plt.figure(fignum, figsize=(4,4)); plt.clf()
        fig.add_axes((0.22,0.15,0.75,0.72))
        axes = plt.gca()
    if 'color' not in kwargs: kwargs['color']='k'
    xmax = len(q)+1
    axes.set_xlim((0,xmax+1))
    axes.set_yscale('log')
    axes.set_ylim((0.1,ymax))
    axes.set_autoscale_on(False)
    for i,r in enumerate(q):
        point(axes, i+1, np.array(r)/norm, 0.5, TS=r[4],  **kwargs)
    axes.set_ylim((0.1,ymax))
    axes.set_xlim((0,xmax+1))

    axes.set_xlabel('month')   
    axes.axhline(1.0, color='k', lw=2)
    axes.set_ylabel('flux / 18-month average')
    axes.set_title('%s' % name, fontsize=12)
    axes.text(0.2*xmax,0.8*ymax, '<flux>=%.2e at %d' % (norm, pivot), fontsize=10)
    axes.grid(True)
    plt.rcParams['axes.linewidth'] = oldlw

    if saveto is not None:
        plt.savefig(os.path.join(saveto, name.replace(' ', '_').replace('+','p')+'__lc.png'))
    return p
    
def make_plots(outdir=outdir, imagedir='light_curves', **kwargs):
    pklist = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
    saveto = os.path.join(outdir, imagedir)
    if not os.path.exists(saveto): os.mkdir(saveto)
    failed= dict()
    for p in pklist:
        pk = pickle.load(open(p))
        print pk['name'], '...',
        try:
            plotmon(pk, color='r', ms=8, saveto=saveto)
            print 'ok!'
        except Exception, arg:
            print 'fail: "%s %s"' % (p,arg)
        #    raise
            failed[p]=arg
    if len(failed.keys())>0: print 'Some failed!'
    return failed    

def test_variability(outdir=outdir, **kwargs):
    """ process the monthly stuff, add info to the pickle files if update==True
    
    """
    pklist = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
    rec = makerec.RecArray('name ts prob sigma delta w new_prob'.split())
    update = kwargs.pop('update', False)
    assert len(kwargs.keys()) ==0, 'unrecognized: %s' % kwargs.keys()
    for p in pklist:
        try:
            pk = pickle.load(open(p))
        except:
            print 'Could load file %s, skipping it.' % p
            continue
        ts = pk['ts2']
        #if ts<25 : continue
        name = pk['name']
        N = kwargs.pop('N', 18)
        try:
            lc = LightCurve(pk, N)
            varindex = 2*lc.delta_loglike()
            prob = lc.survival_prob()
            b,a,c  = [lc.delta_loglike(i) for i in range(3)]
            sigma = np.sqrt(a+c-2*b) # corrected sigma
            delta = (c-a)/sigma
            new_prob = stats.chi2.sf(2*(b-delta**2), lc.N-1 )
            rec.append(name, ts, prob, sigma,  delta, b, new_prob)
            if update:
                pk['variability'] = dict(sigma=sigma, delta=delta, prob=prob, cprob=new_prob, varindex=varindex, N=lc.N)
                print 'updating file %s' % p
                pickle.dump(pk,open(p,'wb'))
        except KeyError:
            print '-------------- fail for %s------------------------' % name
    ret = rec()
    print 'read %d, analyzed %d sources' % (len(pklist), len(ret))
    return ret
        
    
def make_grid(outdir=outdir, title='18M light curves',imagedir='light_curves',  
        html_file=None, **kwargs):
    """create a grid of thumbnails from the files in image_dir, in alphabetical order
    """
    if html_file is None: html_file='%s_lc_grid.htm' %outdir
    cur = os.getcwd()
    os.chdir(outdir)
    images = glob.glob(os.path.join(imagedir,'*.png'))
    assert len(images)>0, 'check that images exist'
    images.sort()
    names = ['%d: %s' %(i+1,(os.path.split(img)[1][:-3]).replace('p','+')) for (i,img) in enumerate(images)]
    tsmap_grid.make_map(images, names,title=title, html_file=html_file, **kwargs)
    os.chdir(cur)
     