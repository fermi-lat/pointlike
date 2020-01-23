"""
Code to generate an ROI counts plot 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/counts_plotter.py,v 1.6 2012/07/12 20:03:04 lande Exp $

Author M. Kerr, T. Burnett

"""
import os
from collections import deque
import numpy as np
import pylab as plt
from matplotlib import font_manager


def counts(r,integral=False):

    groupings = [deque() for x in xrange(len(r.bin_centers))]

    #group slw by energy
    for i,ei in enumerate(r.bin_centers):
        for band in r.bands:
            if band.e == ei:
                groupings[i].append(band)
        groupings[i] = list(groupings[i])

    #iso = np.asarray([ sum((band.bg_counts[1] for band in g)) for g in groupings]) * p
    #gal = np.asarray([ sum((band.bg_counts[0] for band in g)) for g in groupings]) * p
    dif = np.asarray([ np.asarray([band.phase_factor*band.bg_counts for band in g]).sum(axis=0) for g in groupings])
    obs = np.asarray([ sum((band.photons for band in g)) for g in groupings])
    src = np.asarray([ np.asarray([band.phase_factor*band.ps_counts*band.overlaps for band in g]).sum(axis=0) for g in groupings])
    
    if integral:
        for i in xrange(len(iso)):
            #iso[i] = iso[i:].sum()
            #gal[i] = gal[i:].sum()
            dif[i] = dif[i:].sum(axis=0)
            obs[i] = obs[i:].sum()
            src[i] = src[i:].sum(axis=0)

    return r.bin_edges,dif,src,obs,[b.name for b in r.bgm.bgmodels],[p.name for p in r.psm.point_sources]


def get_counts(roi, merge_non_free=True, merge_all=False, integral=False):
    """ reformat output from roi_plotting.counts as a dictionary
    
    """
    en,dif,src,obs,bg_names,ps_names = counts(roi,integral)
    en = (en[1:]*en[:-1])**0.5
    tot = src.sum(axis=1) + dif.sum(axis=1)

    # merge all sources, or ...
    if merge_all:
        new_src    = np.zeros([len(en),1])
        for i in range(src.shape[1]):
            new_src[:,0]+= src[:,i]
        src = new_src
        ps_names = ['All sources']
        
    # optionally merge all of the "frozen" sources for a cleaner plot
    elif merge_non_free:
        free_mask = np.asarray([np.any(m.free) for m in roi.psm.models]) 
        new_src    = np.zeros([len(en),free_mask.sum()+1])
        counter = 0
        for i in xrange(len(free_mask)):
            if free_mask[i] :
                new_src[:,counter] = src[:,i]
                counter += 1
            else:
                new_src[:,-1] += src[:,i]
        src        = new_src
        ps_names = [ps_names[i] for i in xrange(len(ps_names)) if free_mask[i] ]
        ps_names += ['Other Point Sources' ]
    models =  zip(bg_names+ps_names, np.hstack((dif,src)).T)
    return dict(energies=en, observed=obs, models=models, total=tot)


def plot_counts(roi,fignum=1, outfile=None,
        integral=False, max_label=10, merge_non_free=True,
        merge_all=False,
        axes = None,
        ):
    """Make counts and residual plots for the ROI roi
    keyword parameters
        fignum    : integer 
        outfile   : filename, or None,
        integral  : bool
            set True to show integral distributions
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
                obs_kw= dict(linestyle=' ', marker='o', color='black')
                ):
        ax.set_xscale('log')
        ax.set_yscale('log')
        en,obs,tot = count_data['energies'],count_data['observed'],count_data['total']
        for name, data in count_data['models']:
            assert len(en)==len(data), 'energy, data mismatch'
            if len(name)>20: name=name[:17]+'...'
            ax.loglog(en, data, label=name, **model_kw)
        ax.loglog( en, tot, label='Total Model', **total_kw)
        err = obs**0.5
        low_err = np.where(obs-err <= 0, 0.99*obs, err)
        ax.errorbar(en,obs,yerr=[low_err,err], label='Counts', **obs_kw )
        ax.set_ylabel('Counts > E' if integral else 'Counts per Bin')
        ax.set_xlabel('Energy (MeV)')
        prop = font_manager.FontProperties(size='x-small')
        ax.legend(loc=0,prop=prop)
        ax.grid(b=True)
        ax.set_ylim(ymin=0.3)
        
    def plot_residuals(ax, count_data, 
                ylabel='Fractional Residual',
                plot_kw=dict( linestyle=' ', marker='o', color='black',)
                ):
        en,obs,tot = count_data['energies'],count_data['observed'],count_data['total']
        ax.set_xscale('log')
        ax.errorbar(en,(obs-tot)/(tot), yerr=tot**-0.5, **plot_kw)
        ax.axhline(0, color = 'black')
        ybound = min( 0.5, np.abs(np.array(ax.get_ylim())).max() )
        ax.set_ylim((-ybound,ybound))
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Energy (MeV)')
        ax.grid(b=True)

    if axes is None:
        plt.close(fignum) # close it if exists
        fig, axes = plt.subplots(1,2, sharex=True, num=fignum, figsize=(12,6))

    count_data = get_counts(roi, integral=integral, merge_non_free=merge_non_free, merge_all=merge_all)

    plot_counts_and_models(axes[0], count_data)
    if len(axes>1): plot_residuals(axes[1], count_data)
    if outfile is not None: plt.savefig(outfile)


def roi_pipeline_counts_plot(roi, counts_dir=None, fignum=6, title=None, **kwargs):
    """ 
    Code used by the ROI pipeline to make a stacked plot
        
        roi : A ROIanalaysis object
            Uses the name as a title,
        counts_dir : None or the name of a folder
            In the folder case, makes a file name from the ROI name
            
        Creates the two Axes objects, and returns them
    """
    plt.close(fignum)
    oldlw = plt.rcParams['axes.linewidth']
    plt.rcParams['axes.linewidth'] = 2
    fig, axes = plt.subplots(2,1, sharex=True, num=fignum, figsize=(6,8))
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
    axes[1].set_ylabel('fract. dev')
    if title is not None:
        fig.suptitle(title)
    elif hasattr(roi,'name'): 
        fig.suptitle(roi.name)
    if counts_dir is not None:
        if os.path.isdir(counts_dir) and hasattr(roi,'name'):
            fout = os.path.join(counts_dir, ('%s_counts.png'%roi.name) )
        else:
            fout = counts_dir
        fig.savefig(fout)
        print ('saved counts plot to %s' % fout)
    return axes
