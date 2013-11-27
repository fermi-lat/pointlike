"""
Code to generate an ROI counts plot 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/plotting/counts.py,v 1.7 2013/11/26 04:39:59 burnett Exp $

Authors M. Kerr, T. Burnett

"""
import os
import numpy as np
import pylab as plt
from matplotlib import font_manager


def get_counts(roi, event_class=None, tsmin=10):
    """
    return a dictionary with counts information for plotting
    
    roi : an ROIstat object
    event_class : None or integer
        if integer, 0/1 for front/back
    tsmin : float
        used to select sources if bandts info in sedrec exists
        the bandts values are in the dict, None if no info
    """
    bands = roi.selected
    if event_class is not None:
        bands = [b for b in bands if b.band.ec==event_class]
    assert len(bands)>0, 'get_counts: no bands found'
    all_sources = np.array(roi.sources)
    global_mask = np.array([s.skydir is None for s in roi.sources])
    free_mask = np.array([s.skydir is not None and np.any(s.model.free) for s in roi.sources])
    
    global_sources = np.array([ s for s in roi.sources if s.skydir is None])
    free_sources = all_sources[free_mask] 
    strong_filter= lambda s: s.sedrec.ts.sum()>tsmin if hasattr(s,'sedrec')\
            else True if s in free_sources else False
    strong_mask = np.array([strong_filter(s) for s in all_sources])
    weak_mask = free_mask * (~strong_mask)
    
    names = np.array([s.name for s in np.hstack([global_sources,free_sources])])
    bandts = np.array([s.sedrec.ts.sum() if hasattr(s,'sedrec') else None for s in free_sources])
    energies = roi.energies #sorted(list(set([b.energy for b in bands])))
    nume = len(energies)
    numsrc = sum(global_mask | strong_mask)
    observed = np.zeros(nume)
    fixed = np.zeros(nume)
    weak = np.zeros(nume)
    model_counts = np.zeros((nume, numsrc))
    total = np.zeros(nume)
    for i, energy in enumerate(energies):
        for b in bands:
            if b.band.energy!=energy: continue
            counts = np.array([s.counts for s in b]) # array of all model counts
            observed[i] += sum(b.data)
            fixed[i] += b.fixed_counts            
            weak[i]  += sum(counts[weak_mask]) 
            model_counts[i,:]  += counts[global_mask | strong_mask]
            total[i] += b.counts
    models = [(names[j] , model_counts[:,j]) for j in range(numsrc)]
    if sum(weak)>0: models.append( ('TS<%.0f'%tsmin, weak))
    models.append( ('(fixed)', fixed))
    chisq = ((observed-total)**2/total).sum()
    return dict(energies=energies, observed=observed, models=models, 
        names=names, total=total, bandts=bandts, chisq=chisq)
    

def plot_counts(roi,fignum=1, event_class=None, outfile=None,
        max_label=10, #merge_non_free=True,
        #merge_all=False,
        axes = None,
        **kwargs
        ):
    """Make counts and residual plots for the ROI roi
    keyword parameters
        fignum    : integer 
        outfile   : filename, or None,
        max_label : integer
        merge_non_free : bool
            if True (default) combine all sources that are not free
        merge_all: bool
            if True (False default) combine *all* sources
        axes : array of one or two Axes instances, or None 
            If None (default), create standard horizontal layout with figure number fignum
            First is filled with the log plot, second, if there, with the residual 
    """

    def plot_counts_and_models(ax, count_data,
                model_kw = dict(linestyle='-', marker=''),
                total_kw = dict(linestyle='steps-mid', color='black', linewidth=2),
                obs_kw= dict(linestyle=' ', marker='o', color='black',)
                ):
        ax.set_xscale('log')
        ax.set_yscale('log')
        en,obs,tot = count_data['energies'],count_data['observed'],count_data['total']
        for name, data in count_data['models']:
            if np.any(data<=0): continue # ignore models with no predicted counts
            assert len(en)==len(data), 'energy, data mismatch'
            if len(name)>20: name=name[:17]+'...'
            tmodel_kw = model_kw.copy()
            if name.startswith('iso'): tmodel_kw.update(linestyle=':', lw=2)
            elif name.startswith('ring'):tmodel_kw.update(linestyle='--', lw=2)
            ax.loglog(en, data, label=name, **tmodel_kw)
        ax.loglog( en, tot, label='Total Model', **total_kw)
        err = obs**0.5
        low_err = np.where(obs-err <= 0, 0.99*obs, err)
        ax.errorbar(en,obs,yerr=[low_err,err], label='Counts', **obs_kw )
        ax.set_ylabel('Counts per Bin')
        def gevticklabel(x):
            if x<100 or x>1e5: return ''
            elif x==100: return '0.1'
            return '%d'% (x/1e3)
        """ make it look nicer """
        ax.set_xticklabels(map(gevticklabel, ax.get_xticks()))
        ax.set_xlabel(r'$\mathsf{Energy\ (GeV)}$')

        ax.legend(loc=0,prop=dict(size=8))
        ax.grid(b=True)
        ax.set_ylim(ymin=0.3)
        
    def plot_residuals(ax, count_data, 
                ylabel='fract. dev',
                plot_kw=dict( linestyle=' ', marker='o', color='black',),
                show_chisq=True, plot_pulls=True,
                ):
        energy, obs,tot = count_data['energies'],count_data['observed'],count_data['total']
        ax.set_xscale('log')
        if not plot_pulls:
            ax.errorbar(energy, (obs-tot)/(tot), yerr=tot**-0.5, **plot_kw)
            ybound = min( 0.5, np.abs(np.array(ax.get_ylim())).max() )

        else:
            ylabel = 'pull'
            ybound = 3.5
            pull = (obs-tot)/np.sqrt(tot)
            axes[1].plot(energy, pull, 'ko')

            nhigh = sum(pull>3)
            if nhigh>0:  ax.plot(energy[pull>3], [3]*nhigh, '^r', markersize=10) 
            nlow = sum(pull<-3)
            if nlow >0:  ax.plot(energy[pull<-3], [-3]*nlow, 'vr', markersize=10) 
            plt.setp(ax, xscale='log', ylabel='pull', ylim=(-3.5,3.5) )

        ax.axhline(0, color = 'black')
        ax.set_ylim((-ybound,ybound))
        ax.set_ylabel(ylabel)
        def gevticklabel(x):
            if x<100 or x>1e5: return ''
            elif x==100: return '0.1'
            return '%d'% (x/1e3)
        """ make it look nicer """
        ax.set_xticklabels(map(gevticklabel, ax.get_xticks()))
        ax.set_xlabel(r'$\mathsf{Energy\ (GeV)}$')
        ax.grid(b=True)
        if show_chisq :
            ax.text(0.8, 0.8,'chisq=%.0f'% count_data['chisq'], transform = ax.transAxes, fontsize=8)


    if axes is None:
        plt.close(fignum) # close it if exists
        fig, axes = plt.subplots(1,2, sharex=True, num=fignum, figsize=(12,6))

    tsmin = kwargs.pop('tsmin', 10)
    count_data = get_counts(roi, event_class, tsmin=tsmin
    ) #, integral=integral, merge_non_free=merge_non_free, merge_all=merge_all)

    plot_counts_and_models(axes[0], count_data)
    if len(axes>1): plot_residuals(axes[1], count_data, kwargs.pop('show_chisq', True))
    if outfile is not None: plt.savefig(outfile)


def stacked_plots(roi, counts_dir=None, fignum=6, title=None, **kwargs):
    """ 
    Make stacked plots
        
        roi : A ROIstat object
            Uses the name as a title unless title specified
        counts_dir : None or the name of a folder
            In the folder case, makes a file name from the ROI name
            
        Creates the two Axes objects, and returns them
    """
    plt.close(fignum)
    oldlw = plt.rcParams['axes.linewidth']
    plt.rcParams['axes.linewidth'] = 2
    fig, axes = plt.subplots(2,1, sharex=True, num=fignum, figsize=(5,6))
    fig.subplots_adjust(hspace=0)
    axes[0].tick_params(labelbottom='off')
    left, bottom, width, height = (0.15, 0.10, 0.75, 0.85)
    fraction = 0.8

    axes[0].set_position([left, bottom+(1-fraction)*height, width, fraction*height])
    axes[1].set_position([left, bottom, width, (1-fraction)*height])
    plot_counts(roi, axes=axes, outfile=None, **kwargs)
    
    plt.rcParams['axes.linewidth'] = oldlw

    axes[0].set_xlabel('') 
    axes[0].set_ylim(ymin=0.3)
    #axes[1].set_ylabel('fract. dev')
    if title is None:
        if hasattr(roi,'name'): fig.suptitle(roi.name)
    else: fig.suptitle(title)
    if counts_dir is not None:
        if os.path.isdir(counts_dir) and hasattr(roi,'name'):
            fout = os.path.join(counts_dir, ('%s_counts.png'%roi.name) )
        else:
            fout = counts_dir
        fig.savefig(fout)
        print 'saved counts plot to %s' % fout
    return axes
