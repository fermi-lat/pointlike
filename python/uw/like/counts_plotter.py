"""
Code to generate an ROI counts plot 
$Header$

Author M. Kerr, T. Burnett

"""
import numpy as np
import pylab as plt
from matplotlib import font_manager
from uw.like import roi_plotting


def get_counts(roi, merge_non_free=True, merge_all=False, integral=False):
    """ reformat output from roi_plotting.counts as a dictionary
    
    """
    en,dif,src,obs,bg_names,ps_names = roi_plotting.counts(roi,integral)
    en = (en[1:]*en[:-1])**0.5
    tot = src.sum(axis=1) + dif.sum(axis=1)

    # optionally merge all of the "frozen" sources for a cleaner plot
    if merge_non_free:
        free_mask = np.asarray([np.any(m.free) for m in roi.psm.models])
        new_src    = np.zeros([len(en),free_mask.sum()+1])
        counter = 0
        for i in xrange(len(free_mask)):
            if free_mask[i]:
                new_src[:,counter] = src[:,i]
                counter += 1
            else:
                new_src[:,-1] += src[:,i]
        src        = new_src
        ps_names = [ps_names[i] for i in xrange(len(ps_names)) if free_mask[i]]
        ps_names += ['Other Point Sources']
    models =  zip(bg_names+ps_names, np.hstack((dif,src)).T)
    return dict(energies=en, observed=obs, models=models, total=tot)


def plot_counts(roi,fignum=1, outfile=None,
        integral=False, max_label=10, merge_non_free=True,
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
        ybound = min( 0.5, np.abs(np.array(ax.get_ylim()).max() ))
        ax.set_ylim((-ybound,ybound))
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Energy (MeV)')
        ax.grid(b=True)

    if axes is None:
        plt.close(fignum) # close it if exists
        fig, axes = plt.subplots(1,2, sharex=True, num=fignum, figsize=(12,6))

    count_data = get_counts(roi,merge_non_free, integral)

    plot_counts_and_models(axes[0], count_data)
    if len(axes>1): plot_residuals(axes[1], count_data)
    if outfile is not None: plt.savefig(outfile)

