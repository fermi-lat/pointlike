"""
Residual maps

$Header$

"""

import os, pickle, glob
import numpy as np
import pylab as plt
import pandas as pd

from . import analysis_base
from . analysis_base import FloatFormat
from uw.like2.pub import healpix_map
from uw.like2 import tools, configuration
from skymaps import Band

class BandAnalysis(object):
    """ Manage data for a given band index
    """
    def __init__(self, dd,  band_index, nside=None):
        """
        dd : list of dicts with info for each ROI, in order
        band_index : integer
            the band index
        nside: integer, default None
            if specified, a smaller nside to sample down to
        """
       # these must be the same for all ROIs, get from the first
        i=band_index
        self.energy=energy = dd[0][i]['energy']
        event_type = dd[0][i]['event_type']
        self.nside = dd[0][i]['nside']
        nexpect=12*self.nside**2

        if 'inside' in dd[0][0]:
            # new format: all pixels for each ROI, so need to select those in the HEALPix dimond
            ids,model_s,data_s = [np.concatenate([d[i][q][d[i]['inside']] for d in dd]) 
                for q in ['ids','model','data']]
        else:
            # old format: preselected
            ids,model_s,data_s = [np.concatenate([d[i][q] for d in dd]) 
                for q in ['ids','model','data']]

        if nside is None:
            # No resampling
            if len(ids) != nexpect:
                print 'Warning: for nside={} expected {}, but missing {}'.format(self.nside, nexpect,
                    nexpect-len(ids))

            s = np.argsort(ids)
            self.model = model_s[s]
            self.data = data_s[s]
 
        else:
            # resample
            old_nside=self.nside
            self.nside = nside
            print 'Combining bins from nside={} to {}'.format(old_nside, nside)
            bnew=Band(nside).index
            bold=Band(old_nside).dir

            cc = map( lambda i:bnew(bold(i)), range(12*old_nside**2))
            self.data=np.zeros(12*nside**2)
            self.model=np.zeros(12*nside**2)
            for i,m,d in zip(ids, model_s, data_s):
                j = cc[i]
                self.data[j] +=d
                self.model[j]+=m

        self.hpdata = healpix_map.HParray('f{}_data'.format(energy), self.data)
        self.pulls = healpix_map.HParray('f{} pulls'.format(energy), (self.data-self.model)/np.sqrt(self.model))

        # Load galactic corrections from configuration info
        config=configuration.Configuration('.', quiet=True,)
        cfile=os.path.expandvars('$FERMI/diffuse/'+config.diffuse['ring']['correction'])
        print 'Loading corrections from file "{}"'.format(cfile)
        gc = pd.read_csv(cfile)
        self.galcorr = healpix_map.HParray('galcorr', gc['{}'.format(band_index/2)])

        self.label = '{:.0f} Mev {}'.format(energy, ['Front','Back'][event_type])
        #print self.label

    def __repr__(self):
        return 'BandAnalysis for {}, nside={}'.format(self.label, self.nside)

    
    def offset_map(self, **kwargs):
        u = healpix_map.HParray('', self.data/self.model-1) 
        u.plot(cmap=plt.get_cmap('coolwarm'), **kwargs);

    def offset_profile(self, bins=None, ylim=(-0.1,0.1), **kwargs):
        bins = bins if bins is not None else np.linspace(-90,90,91)
        f = Band(self.nside).dir
        sd = map(f, range(len(self.data)))
        glon = np.array([x.l() for x in sd])
        glat = np.array([x.b() for x in sd])
        center = (glon>60) & (glon<300)
        offset = (self.data/self.model-1)[center]
        prof = tools.Profile(glat[center], offset, bins)
        ax=prof.plot()
        plt.setp(ax,xlabel='glat',ylim=ylim, ylabel='<offset>', title=self.label);
        ax.figure.set_facecolor('white')

 
    def ait_pulls_plot(self):
        fig,ax=plt.subplots(figsize=(16,8))
        t=self.pulls.plot(axes=ax, vmin=-5, vmax=5, cmap=plt.get_cmap('coolwarm'),
            title='{} normalized residuals'.format(self.label))
        t.grid(color='white');
        return fig

    def residual_hist(self):
        """Normalized residual plots
        """
        nside=self.nside
        bdir=Band(nside).dir
        fig=plt.figure()
        glat= [bdir(i).b() for i in range(12*nside**2)]
        lolat,loname=np.array(np.abs(glat)<5, bool), '|b|<5'
        hilat,hiname=np.array(np.abs(glat)>10, bool), '|b|>10'
        hkw=dict(bins=np.linspace(-5,5,51), histtype='step', lw=2, log=True)
        def plotit(cut, name):
            ds = self.pulls.vec[cut]
            label='{} {:.2f} {:.2f}'.format(name, ds.mean(),ds.std())
            plt.hist(ds.clip(-5,5),label=label, **hkw )
        plotit(lolat, loname)
        plotit(hilat, hiname)
        plt.grid(alpha=0.5);
        g=lambda x: np.exp(-x**2/2.)
        x = np.linspace(-4,4,81)
        n=sum(lolat) 
        b =hkw['bins'];delta= b[1]-b[0]
        norm=sum(lolat)/delta/np.sqrt(np.pi)/4
        plt.plot(x, norm*g(x), '--g')
        plt.legend(title='select mean std')
        plt.ylim(ymin=0.8); plt.xlim(-5,5)
        plt.title('Normalized residuals for {}'.format(self.label));
        fig = plt.gcf()
        fig.set_facecolor('w')
        return fig
    
    def zea_plots(self,center=(0,0),size=90, size2=30):
        fig, axx = plt.subplots(2,1, figsize=(12,8), sharex=True)
        lookat=(0,0)
        self.hpdata.plot_ZEA(center=center,size=size, axes=axx[0] ,zea_kw=dict(size2=size2))
        axx[0].set_title('Data')
        self.pulls.plot_ZEA(center=center,size=size,  axes=axx[1], vmin=-5, vmax=5,
             cmap=plt.get_cmap('coolwarm'), zea_kw=dict(size2=size2)).grid(color='white');
        axx[1].set_title('pulls');
        return fig

    def correction_factor_plot(self, xlim=(0.8,1.2)):
        b12dir=Band(12).dir
        glat12= [b12dir(i).b() for i in range(12*12**2)]
        glon12= [b12dir(i).l() for i in range(12*12**2)]
        title='{:.0f} MeV Galactic correction factor'.format(self.energy)
        fig,ax = plt.subplots(figsize=(7,7))
        ax.plot(self.galcorr.vec.clip(*xlim), glat12, 'o');
        ax.axvline(1.0, color='green');
        ax.axhline(10, color='orange')
        ax.axhline(-10, color='orange');
        ax.set_title(title)
        ax.set_xlim(*xlim)
        return fig
    
    def ait_correction_factor_map(self, vmin=0.9, vmax=1.1):
        title='{:.0f} MeV Galactic correction factor'.format(self.energy)
        self.galcorr.plot(vmin=vmin, vmax=vmax, title=title);
    
    def zea_correction_factor_map(self, center=(0,0),size=90,size2=30, vmin=0.9, vmax=1.1):
        title='{:.0f} MeV Galactic correction factor'.format(self.energy)
        fig,ax = plt.subplots(figsize=(12,5))
        self.galcorr.plot_ZEA(center=center,size=size,  axes=ax, vmin=vmin, vmax=vmax, 
            zea_kw=dict(size2=size2),title=title,
            ).grid(color='white')
        return fig     



class ResidualMaps(analysis_base.AnalysisBase):
    """Residual Maps
    <br>Combine the maps of the count residuals made for each ROI and each energy band below 1 GeV.
    Only the pixels within the HEALPix tile are used from each ROI: about 1/3 of the total.
    """
    def setup(self, **kwargs):
        self.plotfolder='residual_maps'
        ff = sorted(glob.glob('residual_maps/*.pickle')); 
        assert len(ff)==1728, 'Did not find all files'
        self.dd = [pickle.load(open(f)) for f in ff]
        self.gc = GalacticCorrectionMaps()
        self.band_index,nside=0,None
        self.ba0 =self.band_analysis(band_index=0)
        plt.style.use('seaborn-bright')

    def band_analysis(self, band_index, nside=None):
        self.ba= BandAnalysis(self.dd, band_index, nside)
        return self.ba
    
    def extract(self, i=0):
        ids,model_s,data_s = [np.concatenate([d[0][q] for d in self.dd]) for q in ['ids','model','data']]
        s = np.argsort(ids)
        model = model_s[s]
        data = data_s[s]
        return model, data

    def offset_map(self, band_index=0, vmin=-0.1, vmax=0.1):
        """Residual maps
        """
        ba = BandAnalysis(self.dd,band_index)
        ba.offset_map(title=ba.label, vmin=vmin, vmax=vmax)
        return plt.gcf()

    def offset_profile(self, band_index=0, **kwargs):
        ba = BandAnalysis(self.dd,band_index)
        fig, ax = plt.subplots(figsize=(8,6))
        ba.offset_profile(ax=ax, **kwargs)
        return fig

    def offset_map_0(self):
        """Offset map for 133 MeV Front

        """
        return self.offset_map(0)
    def offset_map_1(self):
        """Offset map for 133 MeV Back

        """
        return self.offset_map(1)

    def offset_profile_0(self):
        """Offset profile for 133 MeV Front

        Selected latitudes 60 degrees about GC 
        """
        return self.offset_profile(0)

    def offset_profile_1(self):
        """Offset profile for 133 MeV Back

        Selected latitudes 60 degrees about GC 
        """
        return self.offset_profile(1)
    def gcmap(self):
        """Galactic correction map for 133 MeV
        """
        return self.gc.plot_map();
    def gcplots(self):
        """Galactic corrections
        Scan through the galactic plane.
        """
        return self.gc.gcplots()

    def residual_maps_ait(self):
        """All-sky Residual maps 

        """
        return self.ba0.ait_pulls_plot()

    def residual_maps_zea(self):
        """Inner Galaxy residual maps
        """
        return self.ba0.zea_plots();
    def residual_hist(self):
        """Residual histogram
        """
        return self.ba0.residual_hist();

    def all_plots(self, **kw):
        self.runfigures([
                 #self.offset_profile_0, self.offset_profile_1,
            self.gcmap,
            self.gcplots,
            self.residual_maps_ait,
            self.residual_maps_zea,
            self.residual_hist,
            self.offset_map_0, self.offset_map_1, 
       ]) 

class GalacticCorrectionMaps(object):
    def __init__(self):
        config=configuration.Configuration('.', quiet=True,)
        self.gc = gc=pd.read_csv(os.path.expandvars('$FERMI/diffuse/'+config.diffuse['ring']['correction']))
        sdirs = map(Band(12).dir, range(1728))
        gc['glat']= map(lambda s:s.b(), sdirs)
        gc['glon']= map(lambda s: s.l(),sdirs)
        gc.loc[gc.glon>180,'glon'] -= 360   

    def plot_map(self, energy_index=0 ):
        gc0 = self.gc['0']
        c0 = healpix_map.HParray('gc0', gc0/gc0.mean())
        c0.plot(vmin=0.6, vmax=1.4, cmap=plt.get_cmap('coolwarm'),
                title='133 MeV Galactic correction factor relative to {:.2f}'.format(gc0.mean()));
        return plt.gcf()

    def gcplots(self, ii=range(0,5), glon_cut=30,ylim=(0.5,1.6)):
        def gcplot(ct, ax=None, energy_index=0):
            if ax is None:
                fig, ax = plt.subplots(figsize=(8,5))
            energy = 10**(2.125 + 0.25*(energy_index))
            ax.plot(ct.glat, ct['{}'.format(energy_index)], 'd');
            ax.grid(alpha=0.5)
            ax.axhline(1, color='k')
            ax.text(0.05,0.9, '{:.0f} MeV '.format(energy),
                transform=ax.transAxes)
        ct =self.gc.query('abs(glon)<{}'.format(glon_cut))
        fig, axx = plt.subplots(len(ii),1, figsize=(8,2.*len(ii)),sharex=True, sharey=True)
        plt.subplots_adjust(top=0.95, hspace=0.05)
        for i,ax in zip(ii, axx):
            gcplot(ct, ax,i)
        axx[0].set_ylim(ylim)
        plt.setp(axx[-1], xlabel='glat', ylabel='correction factor');
        fig.suptitle('Galactic correction factors for ROIs with |glon|<{}'.format(glon_cut));
        return fig
