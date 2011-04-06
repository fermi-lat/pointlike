"""
Plots to study residuals
$Header$
"""
import os, glob, pickle, types
import numpy as np
import pylab as plt
from ..pub import display_map # for skyplot
from ..pub import healpix_map # for HPskyfun
from uw.utilities import makerec, image
from skymaps import SkyDir

def get_outdir():
    version = int(open('version.txt').read())
    ret = 'uw%02d' % version
    print 'using folder %s' % ret
    return ret

def load_count_data(outdir=None, latexdir=r'W:\catalog\Extended_archive_v08' ):
    if outdir is None:
        outdir=get_outdir()
    if 'LATEXTDIR' not in os.environ: # this is now required to unpickle stuff
        os.environ['LATEXTDIR']=latexdir 
    pkfiles = glob.glob(os.path.join(outdir,'pickle', '*.pickle'))
    clist = [pickle.load(open(fname)) for fname in np.sort(pkfiles)]
    if len(clist)<1728: print 'warning: missing pickle files? found %d' %len(clist)
    return clist
    
def make_rec(clist=None):
    if clist is None: clist = load_count_data()
    if type(clist)==types.StringType: 
        clist = load_count_data(clist)
    rec = makerec.RecArray('hp12 ra dec glat glon ebin energy obs model'.split())
    for hp12,roi in enumerate(clist):
        sdir = roi['skydir']
        ra, dec, glat, glon = sdir.ra(), sdir.dec(), sdir.b(), sdir.l()
        counts = roi['counts']
        ebin = range(len(counts['energies']))
        for ebin, energy, obs, model in zip(ebin, counts['energies'], counts['observed'],counts['total']):
            rec.append(hp12, ra, dec, glat, glon, ebin, energy, obs, model)
    return rec()
 
def multiskydev(r=None, title=None, 
        cb_label='deviation relative to mean, sigma units', ait_kw={} ):
    if r is None:
        r = make_rec()

    energies =np.sort(list(set(r.energy)))
    nbins = len(energies)
    devs = np.array([(r.obs/r.model-1)[r.ebin==n] for n in range(nbins)])
    fig, axes = plt.subplots(4,4, figsize=(10,10))
    flatax = axes.flatten()
    for ax, dev, energy in zip(flatax[15-nbins:], devs, energies):
        mean, rms = dev.mean(), dev.std()
        print 'generating plot for energy bin %.0f MeV, mean, rms=%.2g,%.2g' % (energy,mean, rms)
        display_map.skyplot( dev/rms, vmin=-5, vmax=5, #(dev-dev.mean())/dev.std(), vmin=-4, vmax=4
            title='%.1f GeV' %(energy/1e3), title_kw=dict(fontsize=12),
            axes=ax,  nocolorbar=True,  ait_kw=ait_kw)
        #skydev(dev, energy, axes=ax, title_kw=dict(fontsize=12))
        ax.text( 0, -0.2, 'mean:%.1f%% rms=%.1f%%'% (100*mean,100*rms),fontsize=10, transform=ax.transAxes)
    flatax[0].set_axis_off()
    flatax[-1].set_axis_off()
    if title is not None: 
        fig.suptitle(title, fontsize=14)
    cax = fig.add_axes((0.2, 0.2, 0.6, 0.02))
    cb = fig.colorbar(flatax[1].get_images()[0], cax=cax, orientation='horizontal')
    fig.subplots_adjust(wspace=0.1, hspace=-0.5,top=1.05, bottom=0.15) #empirical
    cb.set_label(cb_label) #'deviation relative to mean, sigma units')
    fig.show()    

def multihpsky(r=None):
    """ like above, but just make a list of tables
    """
    if r is None: r=make_rec()
    energies =np.sort(list(set(r.energy)))
    nbins = len(energies)
    devs = np.array([(r.obs/r.model-1)[r.ebin==n] for n in range(nbins)])
    hparrays = []
    for dev, energy in zip(devs, energies):
        mean, rms = dev.mean(), dev.std()
        print 'generating array of deviations for energy bin %.0f MeV, mean, rms=%.2g,%.2g' % (energy,mean, rms)
        hparrays.append( healpix_map.HParray('%.1f GeV' %(energy/1e3),dev/rms) )
    return hparrays

def chisq(r=None):
    if r is None: r=make_rec()
    energies =np.sort(list(set(r.energy)))
    nbins = len(energies)
    sigs = np.array([((r.obs-r.model)**2/r.model)[r.ebin==n] for n in range(nbins)])
    return np.array([sigs[:,i].sum() for i in range(1728)])

        
def tofits(skymodelname, filename):
    r = make_rec(skymodelname)
    hpa = multihpsky(r)
    hpa.append(healpix_map.HParray('HEALPix index', np.arange(1728)))
    hpa.append(healpix_map.HParray('chi squared', chisq(r)))

    t = healpix_map.HEALPixFITS(hpa, nside=64)
    t.write(filename)

  