"""
Count plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/counts.py,v 1.14 2015/12/03 17:10:03 burnett Exp $

"""

import os, pickle
import numpy as np
import pylab as plt
import pandas as pd

from . import analysis_base
from . analysis_base import FloatFormat
from ..pipeline import stream

from uw.utilities import makepivot

class CountPlots(analysis_base.AnalysisBase): 
    """Count plots
    <br> Plots generated after each iteration, checking quality of counts histogram
    %(iteration_info)s
    %(missing_info)s
    """ 
    require='pickle.zip'
    def setup(self, **kwargs):
        self.plotfolder = 'counts'

         # get the basic pickles with the model
        files, pkls = self.load_pickles()
        self.pkls = pkls # for development
        assert len(pkls)==1728, 'expect to find 1728 pickled roi files'
        def chisq10(counts):
            total, observed = counts['total'], counts['observed']
            return ((observed-total)**2/total)[:8].sum()
        def lat180(l): return l if l<180 else l-360
        rdict = dict()
        skipped = 0
        for r in pkls:
            if 'counts' not in r or r['counts'] is None:
                print '***No counts in %s: skipping' % r['name']
                skipped +=1
                continue
            cnts = r['counts']
            rdict[r['name']]=dict(
                glon = lat180(r['skydir'].l()),
                glat = r['skydir'].b(),
                chisq = cnts['chisq'],
                chisq10= chisq10(r['counts']),
                uchisq = cnts['uchisq'] if 'uchisq' in cnts else 0,
                )
        if skipped>0:
            self.missing_info = '<p>%d missing ROIS' % skipped
            print '***Skipped %d ROIs' % skipped
        else: self.missing_info = ''
        self.rois = pd.DataFrame(rdict).transpose()
        self.rois['singlat'] = np.sin(np.radians(np.asarray(self.rois.glat,float)))
        self.rois['glon'] = np.asarray(self.rois.glon, float)
        
        # dict of dataframes with count info. columns are energies
        self.energy = self.pkls[0]['counts']['energies'] # extract list from first pickle
        counts = [p['counts'] for p in self.pkls if p['counts'] is not None]
        self.counts=dict()
        for key in ['observed', 'total']:
            self.counts[key]= pd.DataFrame([x[key] for x in counts], index=self.rois.index)
        #try:
        self.add_model_info()
        #except Exception, msg:
        #    print msg
            
        if 'history' in pkls[0].keys():
            print 'Extracting history info from the ROI analyses'
            self.sinfo =t= self.history_info()
            # check for previous creation, ignore them: look for last "monthly*" or "create"
            y= [(x.startswith('monthly') or (x.startswith('create'))) for x in t.stage]; 
            lc = t.index[y][-1]; print lc
            self.toshow = t[t.index>= lc]; print self.toshow
            skipped = t.index[y][:-1]; print 'Skipped starts:\n{}'.format(t.ix[skipped])
            
            cfile = 'config.txt' if os.path.exists('config.txt') else '../config.txt'
            try:
                config = eval(open(cfile).read()) 
            except Exception, msg:
                raise Exception('Could not read config file, %s' %msg)
            input_model=config['input_model']['path']
            self.iteration_info = """<p>Input model: <a href="../../%s/plots/index.html?skipDecoration">%s</a>
                <p>Iteration history: log likelihood change for each step: \n%s
                """ % (input_model,input_model, 
                    self.toshow.to_html(float_format=FloatFormat(1)) )
        else:
            self.iteration_info=''

    def history_info(self):
        from uw.like2.pipeline import stream
        model = '/'.join(os.getcwd().split('/')[-2:])
        sinfo = stream.StreamInfo(model)
        
        #look up info from stream creation
        model_streams =sorted(sinfo.keys()); # streams in this model
        hd= dict()
        for s in model_streams:
            hd[s]=np.zeros(1728)
        interactive = 0
        for i,r in enumerate(self.pkls):
            hh = r['history']
            for h in hh:
                s = h['stream']
                if s=='interactive': 
                    interactive +=1
                    continue
                stream = int(h['stream'].split('.')[0])
                if stream not in model_streams: continue
                hd[stream][i]=h['logl']
        p = model_streams[0]

        sinfo[p].update(delta_min=0, delta_max=0, delta_sum=0,gt10=0, nroi=1728)
        for s  in model_streams[1:]:
            nroi = sum(hd[s]!=0)
            for k in range(1728):
                if hd[s][k]==0:
                    hd[s][k]= hd[p][k]
            delta = hd[s]-hd[p]
            sinfo[s].update(delta_sum=sum(delta), delta_min=delta.min(), 
                            gt10=sum(abs(delta)>10), delta_max=delta.max(), nroi=nroi)
            p = s
        return pd.DataFrame(sinfo).T['stage date nroi gt10 delta_sum delta_min delta_max'.split()]
    
    def loglikelihood(self):
        """log likelihood

        Progression of the log likelihood (left scale, blue cirles) and the number of ROI's changing by more than 10 (right scale, red diamonds) 
        for the processing stages.
        """
        fig, ax = plt.subplots(1,1, figsize=(8,4))
        sinfo = self.toshow #subset
        ds = sinfo.delta_sum
        x = range(len(sinfo))
        ax.plot( x, np.cumsum(ds), 'o--');
        ax.set_xticks(x)
        ax.set_xticklabels(list(sinfo.stage), rotation=45, va='top', ha='right')
        plt.setp(ax, xlim=(x[0]-0.5, x[-1]+0.5), 
                 title='change in total log likelihood', xlabel='processing stage', ylabel='log likelihood')
        ax.grid(True, alpha=0.5)
        axr = ax.twinx()
        plt.setp(axr, ylabel='changed ROIs', ylim=(0,sinfo.gt10.max()+5) )
        axr.plot( x, sinfo.gt10, 'Dr--')
        return fig    
    
    def add_model_info(self):
        for i,key in enumerate(['ring','isotrop', 'SunMoon', 'limb',]): # the expected order
            t = []
            for j,p in enumerate(self.pkls):
                if p['counts'] is None: continue
                if key in p['diffuse_names']:
                    k =  p['diffuse_names'].index(key) 
                    y=p['counts']['models'][k]
                    assert y[0]==key, 'bad key, roi %d: %s!=%s; list is %s'% (j,key, y[0], p['diffuse_names'])
                    t.append(y[1])
                else:
                    t.append(np.zeros(len(self.energy)))
            self.counts[key]= pd.DataFrame(t, index=self.rois.index)

    def counts_map(self):
        """ Sum, for E>100 Mev
        """
        obs = self.counts['observed']
        total = np.array([sum(x[1]) for x in obs.iterrows()])

    def residual(self, ib):
        """ residual DF array for energy band ib 
        """
        obs   = self.counts['observed'].transpose().ix[ib]
        model = self.counts['total'].transpose().ix[ib]
        resid = (obs-model)/np.sqrt(model)
        return resid
     
    def residual_hists(self):
        """ histograms of normalized residuals 
        subset for ridge (|b|<10, |l|<60) shown
        """
        fig,axx = plt.subplots(3,4, figsize=(12,9), sharex=True,)
        ridge = ( np.abs(self.rois.glat)<10) & ( np.abs(self.rois.glon)<60 )

        for ib,ax in enumerate(axx.flatten()):
            resid = self.residual(ib)
            hkw = dict(bins=np.linspace(-5,5,21), histtype='stepfilled')
            ax.hist(resid.clip(-5,5), **hkw)
            ax.hist(resid[ridge].clip(-5,5), color='orange', **hkw)
            ax.set_title('%.0f MeV'% self.energy[ib], fontsize=10)
            ax.axvline(0, color='k')
            plt.setp(ax, xlim=(-5,5))
            ax.grid(True)
        return fig
    
    def residual_plot(self):
        """ plot of the average normalized residual
        """
        res = [self.residual(ib) for ib in range(len(self.energy))]
        means = [x.mean() for x in res]
        stds  = [x.std() for x in res]
        fig,ax= plt.subplots(figsize=(4,4))
        ax.errorbar(self.energy, y=means, yerr=stds, fmt='o')
        plt.setp(ax, xscale='log', xlabel='energy', ylabel='average residual')
        ax.set_title('count residuals')
        ax.grid()
        ax.axhline(0, color='k')
        return fig
        
    def chisq_plots(self, use10=False, unweight=True, hsize=(1.0, 0.7, 1.5, 0.7), 
            vmin=0, vmax=50, bcut=10, grid_flag=True):
        """ chi squared plots
        chi squared distribution
        <p>Note that this chi squared is modified by the unweighting factors.
        %(bad_roi_link)s
        """
        
        
        fig, axs = self.subplot_array( hsize, figsize=(11,5))
        if unweight:
            chisq = self.rois.uchisq
        else:
            chisq = self.rois.chisq if not use10 else self.rois.chisq10
        chisqtxt= r'$\chi^2$'
        
        # make a table of the bad ones
        
        bad_rois = self.rois[chisq>vmax]['glat glon chisq'.split()]
        bad_rois['chisq'] = chisq
        bad_rois.to_csv('bad_rois.csv')
        if not self.skymodel.startswith('month'):
            try:
                pc =makepivot.MakeCollection('bad rois %s' % os.path.split(os.getcwd())[-1], 'countfig', 'bad_rois.csv')
            
                self.bad_roi_link = """\
                    <p>A list of %d bad ROIs, with chisq>%.0f, can examined with a 
                    <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
                    which requires Silverlight."""  % (len(bad_rois), vmax, pc.cId)
            except Exception, msg:
                self.bad_roi_link = 'Failed to create Pivot: {}'.format(msg)
                print self.bad_roi_link
            
        else:
            self.bad_roi_link=''
        
        def chisky(ax):
            self.basic_skyplot(ax, self.rois.glon, self.rois.singlat, chisq, 
                s=55, marker='D', vmin=vmin, vmax=vmax,  edgecolor='none', 
                colorbar=True, cbtext=chisqtxt);
                
        def chihist(ax, htype='stepfilled'):
            bins = np.linspace(0, vmax, 26)
            lolat = np.abs(self.rois.glat)<bcut
            ax.hist(chisq.clip(0,vmax), bins, label='all: mean=%.0f'%chisq.mean(), histtype=htype)
            ax.hist(chisq.clip(0,vmax)[lolat], bins, color='orange', 
                label='|b|<%d (%.0f)'%(bcut, chisq[lolat].mean()), histtype=htype)
            ax.legend(loc='upper right', prop=dict(size=10)) 
            plt.setp(ax, xlabel=chisqtxt, xlim=(0,vmax))
            ax.grid(grid_flag)
            
        for f, ax in zip( (chihist, chisky), axs.flatten()): f(ax)
        return fig
        
    def residual_maps(self, vmin=-5, vmax=5):
        """ Maps of the residuals 
        """
        fig, axx = plt.subplots(3,4, figsize=(12,10), sharex=True, sharey=True)
        plt.subplots_adjust(right=0.9, hspace=0.15, wspace=0.1)
        for ib,energy in enumerate(self.energy[:12]):
            ax = axx.flatten()[ib]
            scat=self.basic_skyplot(ax, self.rois.glon, self.rois.singlat, self.residual(ib).clip(vmin,vmax),
                 title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, s=15, edgecolor='none', colorbar=False, labels=False)
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('normalized residual')
        fig.text(0.5,0.05, 'galactic longitude', ha='center')
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center')
        return fig
        
    def resid_vs_dec(self, ib=0, ax=None, ylim=(-8,8), labels=True):
        """ residual vs. declination angle
        """
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        r = self.residual(ib).clip(*ylim)
        ax.plot(self.rois.dec, r, '.', color='gray')
        galplane = np.abs(self.rois.glat)<5
        ax.plot(self.rois.dec[galplane], r[galplane], 'or', label='|b|<5')
        ax.axhline(0, color='k')
        plt.setp(ax, ylim=ylim, xlim=(-90,90))
        if labels: plt.setp(ax, xlabel='Dec', ylabel='normalized residual')
        ax.set_xticks((-90,-45,0,45,90))
        ax.set_title('%d MeV' % self.energy[ib], size=10)
        ax.legend(prop=dict(size=10))
        ax.grid()
        return fig
        
    def ridge_spectral_residuals(self, ax=None, glat=(-10,10), glon=(-60,60), nolabels=False):
        """ Spectral residuals along the Galactic ridge
       
        Designed to match, except for different ROI definitions, the Saclay standard plot.
        """
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        def cut( x, range):
            return (x>range[0])  & (x<range[1])
        #ridge = ( np.abs(self.rois.glat)<10) * ( np.abs(self.rois.glon)<60 )
        ridge = cut(self.rois.glat, glat) & cut(self.rois.glon, glon)
        data =self.counts['observed'][ridge].sum()
        model = self.counts['total'][ridge].sum()
        x = self.energy
        y = data/model-1
        yerr = 1/np.sqrt(model) # needs correlation factor
        self.ridge_correction = pd.DataFrame(dict(energy=x, corr=y, error=yerr))
        ax.errorbar(x, y, yerr=yerr, fmt='o')
        plt.setp(ax, xscale='log')
        if not nolabels:
            plt.setp(ax, xlabel='energy (MeV)', ylabel='(counts-model)/model')
        ax.axhline(0, color='gray')
        ax.grid()
        return fig

    def all_plots(self):
        self.runfigures([
            self.loglikelihood,
            self.chisq_plots,
            self.residual_maps, 
            self.residual_plot, 
            self.residual_hists, 
            self.ridge_spectral_residuals,
            ])
