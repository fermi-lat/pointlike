""" Manage fit quality  for sourceinfo.SourceInfo

"""
import numpy as np
import pylab as plt
import pandas as pd
from .analysis_base import (html_table, FloatFormat)
from scipy import stats

class FitQualityPlots(object):
    """ Spectral fit quality
    
    This is the difference between the TS from the fits in the individual energy bands, and that for the spectral fit.
    It should be distributed approximately as chi squared of the number of degrees of freedom. 
    However, high energy bins usually do not contribute, there is a range of effective bins that peaks at 8.
    So we compare with ndf= 9 or 8.
    All sources with TS>%(tscut)d are shown.<br>
    <b>Left</b>: Power-law fits. Tails in this distribution perhaps could be improved by changing the curvature 
    parameter. 
    <br><b>Center</b>: Log parabola fits.
    <br><b>Right</b>: Fits for the pulsars, showing high latitude subset.
    <br>
    %(badfit_check)s

    %(poorfit_table)s

    """
    
    def __init__(self, sourceinfo, xlim=(0,30), ndf=(9,8,8), tscut=25, grid_flag=True, make_table=True, legend_flag=True):
        
        from scipy import stats
        self.sourceinfo = sourceinfo
        s = self.df = sourceinfo.df
        psr = np.asarray(s.psr, bool)
        fq = np.array(s.fitqual, float)
        fq[pd.isnull(fq)]=0
        beta = s.beta
        logparabola = (~psr) & (beta>0.01)
        powerlaw = (~psr) & (beta.isnull() | (beta<0.01) )

        cut = np.array((s.ts>tscut) & (fq>0) , bool)
        
        dom = np.linspace(xlim[0],xlim[1],26)
        d = np.linspace(xlim[0],xlim[1],51); delta=dom[1]-dom[0]
        chi2 = lambda x,ndf: stats.chi2.pdf(x,ndf) * delta

        hilat = np.abs(s.glat)>5
        def tobool(a): return np.array(a, bool)
        
        fig, axx = plt.subplots(1,3, figsize=(15,6))
        plt.subplots_adjust(left=0.1)
        self.average=[]
        def mystats(cut):
            return '\nmean: {:.1f} tail: {:d}'.format(fq[cut].mean(), sum(fq[cut]>xlim[1]))
        
        def make_hist(ax, n, label, cut_expr):
            mycut=tobool(cut & (cut_expr))
            count = sum(mycut)
            if count==0:
                print ('Not generating plot for %s' % label)
                return

            ax.hist(fq[mycut].clip(*xlim), dom, histtype='stepfilled', 
                    color='lightblue',edgecolor='blue', 
                    label=label+ mystats(mycut))
            self.average.append( fq[mycut].mean() )
            if label=='exp. cutoff':
                ax.hist(fq[mycut & hilat].clip(*xlim), dom,
                        histtype='stepfilled', color='palegoldenrod', edgecolor='orange',
                        label='|b|>5 subset'+mystats(mycut&hilat))
            self.average.append(fq[mycut & hilat].mean())
            ax.plot(d, chi2(d,n)*count, 'r--', lw=2, label=r'$\mathsf{\chi^2\ ndf=%d}$'%n)
            ax.grid(alpha=0.25); ax.set_xlabel('fit quality')
            if legend_flag: ax.legend(loc='upper right', prop=dict(size=10))
            else: ax.set_title(label)
                
        for ax, n, label, cut_expr in zip(axx, ndf, 
                        ('power law', 'log-normal','exp. cutoff'), 
                        (powerlaw,logparabola, psr)):
            make_hist(ax, n, label, cut_expr)
        self.sourceinfo.ndf = ndf
        self.sourceinfo.tscut=tscut


    def badfit(self, make_table=True):
        
        self.df['badfit2'] =np.array(self.df.badfit.values, bool)
        t = self.df.loc[(self.df.badfit2) & (self.df.ts>10)].sort_values(by='roiname')
        print ('%d sources with bad fits' %len(t))
        if len(t)>0:
            print ('%d sources with missing errors' % len(t))
            self.badfit = t[['ts', 'freebits', 'badbits', 'pindex', 'beta', 'e0','roiname']]
            self.badfit_check = html_table(self.badfit, name=self.sourceinfo.plotfolder+'/badfits', 
                heading='<h4>%d Sources with missing errors</h4>' % len(t), float_format=FloatFormat(1))
            ids = np.array(sorted([int(name[-4:]) for name in set(self.badfit.roiname)]))
            self.badfit_check +='<p>List of roi numbers: {}'.format(np.array(ids)) + """
                <br>Missing errors means that the Hessian matrix in second derivates of the log likelihood could not 
                be completely inverted, which can happen if it is singular due to too high correlations.
                """
        else: self.badfit_check = '<p>All sources fit ok.'
        self.fit_quality_average =  ', '.join( map(lambda x,n :'%s: %.1f' %(n,x) ,
                            self.average, 'powerlaw logparabola expcutoff(hilat) expcutoff(lolat)'.split()) )
        
        print ('fit quality averages:', self.fit_quality_average)
        if make_table:
            s = self.df
            # Make table of the poor fits
            s['pull0'] = np.array([x.pull[0] if x is not None else np.nan for x in s.sedrec ])
            t =self.df.loc[((s.fitqual>30) | (np.abs(s.pull0)>3)) & (s.ts>10) ]\
                ['ra dec glat fitqual pull0 ts modelname freebits index2 roiname'.split()].sort_values(by='roiname')
            if len(t)==0:
                self.poorfit_table= '<p>No poor fits found'
            else:
            
                self.poorfit_table  = html_table(t, name=self.sourceinfo.plotfolder+'/poorfit', 
                        heading='<h4>Table of %d poor spectral fits</h4>'%len(t),
                        float_format=FloatFormat(2),
                        formatters=dict(ra=FloatFormat(3), dec=FloatFormat(3), ts=FloatFormat(0),index2=FloatFormat(3)))
        self.sourceinfo.badfit_check = self.badfit_check
        self.sourceinfo.poorfit_table = self.poorfit_table
