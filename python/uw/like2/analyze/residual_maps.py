"""
Residual maps

"""

import os, pickle, glob, healpy
import numpy as np
import pylab as plt
import pandas as pd

from . import analysis_base
from . analysis_base import FloatFormat
from uw.utilities import healpix_map
hpm = healpix_map
from uw.like2 import (diffuse, tools, configuration,)
from skymaps import Band

class BandAnalysis(object):
    """ Manage data for a given band index
    """
    def __init__(self, dd,  band_index, nside=None, outfile=None):
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

        keys = dd[0][0].keys()
        print 'keys: {}'.format(keys)
        cat = lambda q : np.concatenate([d[i][q][d[i]['inside']] for d in dd])
        if 'inside' in keys:
            # new format: all pixels for each ROI, so need to select those in the HEALPix dimond
            ids,model_s,data_s, source_s, diffuse_s = [cat(q) 
                for q in ['ids','model','data', 'source', 'galactic', ]]
            if 'isotropic' in keys:
                isotropic_s = cat('isotropic')
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
            sources = (source_s)[s]
            galactic = diffuse_s[s]
            if 'isotropic' in keys:
                isotropic = isotropic_s[s]
 
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

        self.hpdata = healpix_map.HParray('data', self.data)
        self.pulls = healpix_map.HParray('pulls', (self.data-self.model)/np.sqrt(self.model))
        self.resid = healpix_map.HParray('residuals', (self.data-self.model)/(self.model))
        self.sources = healpix_map.HParray('sources', sources)
        self.galactic = healpix_map.HParray('galactic', galactic)
        if 'isotropic' in keys:
            self.isotropic =  healpix_map.HParray('isotropic', isotropic)


        # Load galactic corrections from configuration info
        config=configuration.Configuration('.', quiet=True, postpone=True)
        try:
            cfile=os.path.expandvars('$FERMI/diffuse/'+config.diffuse['ring']['correction'])
            print 'Loading corrections from file "{}"'.format(cfile)
            gc = pd.read_csv(cfile)
            self.galcorr = healpix_map.HParray('galcorr', gc['{}'.format(band_index/2)])
        except:
            self.galcorr=None

        self.label = '{:.0f} Mev {}'.format(energy, ['Front','Back'][event_type])
        #print self.label

        if outfile:
            cols = [self.pulls, self.sources, self.galactic, self.hpdata,]
            if 'isotropic' in keys:
                cols = cols + [self.isotropic]
            t=healpix_map.HEALPixFITS( cols)
            t.write(outfile)
            print 'Wrote file {}'.format(outfile)

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

    def ait_pulls_plot(self, ax=None, vlim=5, pub=False):

        fig,ax=plt.subplots(figsize=(16,8)) if ax is None else (ax.figure, ax)
        t=self.pulls.plot(axes=ax, vmin=-vlim, vmax=vlim, cmap=plt.get_cmap('coolwarm'), ait_kw={},
            title='{} normalized residuals'.format(self.label) if not pub else '')
        
        t.grid(color='white');
        return fig

    # def residual_hist(self):
    #     """Normalized residual plots
    #     """
    #     nside=self.nside
    #     bdir=Band(nside).dir
    #     fig=plt.figure()
    #     glat= [bdir(i).b() for i in range(12*nside**2)]
    #     lolat,loname=np.array(np.abs(glat)<5, bool), '|b|< 5'
    #     hilat,hiname=np.array(np.abs(glat)>10, bool), '|b|>10'
    #     hkw=dict(bins=np.linspace(-5,5,51), histtype='step', lw=2, log=True)
    #     def plotit(cut, name):
    #         ds = self.pulls.vec[cut]
    #         label='{} {:4.2f} {:4.2f}'.format(name, ds.mean(),ds.std())
    #         plt.hist(ds.clip(-5,5),label=label, **hkw )
    #     plotit(lolat, loname)
    #     plotit(hilat, hiname)
    #     plt.grid(alpha=0.5);
    #     g=lambda x: np.exp(-x**2/2.)
    #     x = np.linspace(-4,4,81)
    #     n=sum(lolat) 
    #     b =hkw['bins'];delta= b[1]-b[0]
    #     norm=sum(lolat)/delta/np.sqrt(np.pi)/4
    #     plt.plot(x, norm*g(x), '--g')
    #     leg=plt.legend( title='   select mean std',prop=dict(size=10, family='monospace'))
    #     ltit = leg.get_title(); ltit.set_fontsize(10); ltit.set_family('monospace')
    #     plt.ylim(ymin=0.8); plt.xlim(-5,5)
    #     plt.title('Normalized residuals for {}'.format(self.label));
    #     fig = plt.gcf()
    #     fig.set_facecolor('w')
    #     return fig
    
    def residual_hist(self, ax=None, colors='orange green blue'.split(), pub=False):
        """Normalized residual plots
        """
        fig, ax = plt.subplots(figsize=(6,4)) if ax is None else (ax.figure, ax) 
        
        nside=self.pulls.nside
        pulls = self.pulls.vec
        
        bdir=Band(nside).dir
        glat= [bdir(i).b() for i in range(12*nside**2)]
        cuts = [np.array(np.abs(glat)<10, bool),np.array(np.abs(glat)>10,bool) ]
        sel  = [pulls[cut] for cut in cuts]
        labels = '|b|<10', '|b|>10'
        
        def plotit(s, label, color): 
            hkw=dict(bins=np.linspace(-5,5,51), histtype='step', lw=2, log=True)
            label='{:10} {:5.2f} {:5.2f}'.format(label, s.mean(),s.std())
            ax.hist(s.clip(-5,5),label=label,color=color, **hkw )
            #ovelay a gaussian with same total
            g=lambda x: np.exp(-x**2/2.)
            x = np.linspace(-4,4,81)
            b =hkw['bins']; delta= b[1]-b[0]
            norm = len(s) * delta/np.sqrt(2*np.pi)
            ax.plot(x, norm*g(x), '--', color=color) 
        
        for s,label, color in zip(sel, labels, colors):   
            plotit(s, label, color)
      
        leg=ax.legend(loc='lower center',
            title='      {:10} {:5} {:5}'.format('selection', 'mean', 'SD'),
                prop=dict(size=10, family='monospace'))
        ltit = leg.get_title()
        ltit.set_fontsize(10); ltit.set_family('monospace')
        ax.set(ylim=(0.8,None), xlim=(-5,5))
        ax.set(xlabel='Normalized Residual')
        if not pub:
            ax.set_title(' {}'.format(self.label));
            ax.grid(alpha=0.5)
        fig.set_facecolor('w')
        return fig

    def zea_plots(self,center=(0,0),size=90, size2=20):
        import matplotlib.colors as colors
        fig, axx = plt.subplots(3,1, figsize=(18,12), sharex=True)
        plt.subplots_adjust(hspace=0.04) 
        zea_kw, cmap = dict(size2=size2,), plt.get_cmap('coolwarm')
        plasma = plt.get_cmap('plasma')
        labelit = lambda ax, t: ax.text(0.01, 0.90, t, fontsize=16,
                        bbox=dict(facecolor='white', edgecolor='k',alpha=0.5),transform=ax.transAxes) 
        ax=axx[0]
        self.hpdata.plot_ZEA(center=center,size=size, axes=ax ,cmap=plasma, 
            norm=colors.LogNorm(),zea_kw=zea_kw).grid(color='grey')
        labelit(ax,'Data')
        
        self.pulls.plot_ZEA(center=center,size=size,  axes=axx[1], vmin=-5, vmax=5,
             cmap=cmap, zea_kw=zea_kw).grid(color='white');
        labelit(axx[1],'pulls: (-5,5)');

        self.resid.plot_ZEA(center=center,size=size,  axes=axx[2], vmin=-0.2, vmax=0.2,
             cmap=cmap, zea_kw=zea_kw).grid(color='white');
        labelit(axx[2],'residuals: (-0.2,0.2)');
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

class Components(object):
    """For interactive analysis of the all-sky components saved by a ResidualMaps run
    """
    def __init__(self, filename = 'plots/residual_maps/f131maps.fits'):
        from astropy.table import Table
        assert os.path.exists(filename)
        self.df=Table.read(filename, hdu=1).to_pandas()
        print 'Read file {} with columns {}'.format(filename, self.df.columns)


class ResidualMaps(analysis_base.AnalysisBase):
    """<h3>Residual Maps</h3>
    <br>Combine the maps of the count residuals made for each ROI and each energy band below 1 GeV,
    using the pipe run residualmaps.
    Only the pixels within the HEALPix tile are used from each ROI: about 1/3 of the total.
    <br>Most of the residuals here are for 133 MeV Front. 
    """
    def setup(self, **kwargs):
        self.plotfolder='residual_maps'
        ff = sorted(glob.glob('residual_maps/*.pickle')); 
        assert len(ff)==1728, 'Did not find all files'
        self.dd = [pickle.load(open(f)) for f in ff]

        self.band_index,nside=0,None
        try:
            self.ba0 =self.band_analysis(band_index=0, outfile=self.plotfolder+'/f131maps.fits')
            self.ba1 =self.band_analysis(band_index=1, outfile=self.plotfolder+'/f237maps.fits')

        except Exception, msg:
            print 'Fail to load maps: {}'.format(msg)

        plt.style.use('seaborn-bright')

        # generate a DF with ROI positions
        sd12 = map(Band(12).dir, range(1728))
        glon12 = np.array(map(lambda d: round(d.l(),1), sd12))
        glon12[glon12>180]-=360
        glat12 = np.array(map(lambda d: round(d.b(),1), sd12))
        self.df12 = pd.DataFrame(dict(glon=glon12, glat=glat12), index=range(1728))

        try:
            self.gc = GalacticCorrectionMaps()
        except Exception, msg:
            print "No corrections made"
            self.gc=None

    def band_analysis(self, band_index, nside=None, outfile=None):
        self.ba= BandAnalysis(self.dd, band_index, nside, outfile)
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
        ba = self.ba0 if band_index==0 else BandAnalysis(self.dd,band_index)
        ba.offset_map(title=ba.label, vmin=vmin, vmax=vmax)
        return plt.gcf()

    def offset_profile(self, band_index=0, **kwargs):
        ba = self.ba0 if band_index==0 else BandAnalysis(self.dd,band_index)
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
        
    def gcplots(self, glon_cut=60, ii=range(0,3), ylim=(0.75,1.25)):
        """Galactic corrections
        Scan through the galactic plane.
        """
        return self.gc.gcplots(glon_cut=glon_cut, ii=ii, ylim=ylim,)

    def residual_maps_ait(self):
        """All-sky Residual maps 

        """
        return self.ba0.ait_pulls_plot()

    def residual_maps_center(self):
        """Central Galaxy residual maps
        """
        return self.ba0.zea_plots();

    def residual_maps_cygnus(self):
        """Cygnus residual maps

            Brighest source is PSR J2021+4026.
        """
        return self.ba0.zea_plots(center=(90,0));
        
    def residual_maps_anti(self):
        """Anti-center residual maps
        """
        return self.ba0.zea_plots(center=(180,0));

    def residual_maps_vela(self):
        """Vela 133 MeV residual maps
        """
        return self.ba0.zea_plots(center=(-90,0));

    def residual_hist(self):
        """Residual histogram
        """
        return self.ba0.residual_hist();

    def residual_map_and_hist(self, pub=False):
        """Residual map and corresponding histogram

        For 131 MeV band.
        """
        fig,ax=plt.subplots(figsize=(16,8),gridspec_kw=dict(left=0, right=0.7))
        self.ba0.ait_pulls_plot(ax=ax, pub=pub);
        ax2=fig.add_axes([0.7,0.25, 0.3,0.5],label='hist')
        self.ba0.residual_hist( pub=pub, ax=ax2);
        return fig


    def residual_df(self, roi, i=0):
        """return a DF with the residuals for the given ROI number and energy band
        """
        d = self.dd[roi][i]
        nside = d['nside']
        model = d['model']
        index = d['ids']
        bdir = np.array(map(lambda i: Band(nside).dir(i).b(), index)).round(2)
        ldir = np.array(map(lambda i: Band(nside).dir(i).l(), index)).round(2)
        ldir[ldir>180] -=360
        return pd.DataFrame(dict(model=model.round(1), 
                                resid=d['data']-model, inside=d['inside'], 
                                glat=bdir, glon=ldir), index=index)
   
    def residual_scats(self, qstring, band_index=0,ax=None, scat=True, nocolorbar=True):
    
        fig,ax = plt.subplots(figsize=(12,5)) if ax is None else ax.figure, ax
            
        scatkw = dict( marker='d',vmin=-0.1, vmax=0.1, 
                    s=80, cmap=plt.get_cmap('coolwarm')) 
        
        plotkw = dict(marker='o', ls='')
        
        def residual_scat(roi_index):
            df=self.residual_df(roi_index, band_index)
            r = df.resid/df.model    
            return ax.scatter(df.glon, df.glat, c=r, **scatkw)

        def residual_plot(roi_index):
            df=self.residual_df(roi_index, band_index)
            r = df.resid/df.model   
            ax.plot(df.glat, r, color='blue', **plotkw)
            ax.plot(df.glat[df.inside], (r)[df.inside],color='red', **plotkw)
    
        roi_indeces = self.df12.query(qstring).index
        if scat:
            ax.set(xlim=(30,-30), ylim=(-15,15))

            scat =map(residual_scat, roi_indeces)[0]
            if not nocolorbar: plt.colorbar(scat)
        else:
            map(residual_plot, roi_indeces)
            ax.set(xlim=(-15,15), ylim=(-0.2,0.2))

        ax.axvline(0, color='lightgrey')
        ax.axhline(0, color='lightgrey')
        ax.set_title('band {}, {}'.format(band_index,qstring))
        return fig

    def ridge_systematics_scat(self):
        """Ridge systematic scatter plots

        scatter plots of pixels in subsets of ROIs in the ridge area
        color shows sytematic offset (range -20 to 20 pct.)
        """
        fig,axx = plt.subplots(3,1, figsize=(12,15), sharex=True)  
        self.residual_scats('abs(glat)==6.4 and abs(glon)<30', ax=axx[0])
        self.residual_scats('(abs(glat)==9.6 or glat==0) and abs(glon)<30', ax=axx[1]);
        self.residual_scats('(abs(glat)==3.2) and abs(glon)<30', ax=axx[2])
        return fig

    def ridge_systematics_plot(self, band_index=0):
        """Ridge systematic plots
        for selected ROIs along ridge, plots of fractional residual vs. galactic latitude.
            red points are pisles within the ROI active area, blue out to the 5 deg radius. 
        """
        kw = dict(scat=False, band_index=band_index)
        fig,axx = plt.subplots(3,1, figsize=(12,15), sharex=True)  
        self.residual_scats( 'abs(glat)==6.4 and abs(glon)<30', ax=axx[0], **kw);
        self.residual_scats( 'abs(glat)==3.2 and abs(glon)<30', ax=axx[1], **kw);
        self.residual_scats( '(abs(glat)==9.6 or glat==0) and abs(glon)<30', ax=axx[2], **kw);
        return fig

    def all_plots(self, **kw):

        self.runfigures([
                 #self.offset_profile_0, self.offset_profile_1,
            self.gcmap,
            self.gcplots,
            self.residual_maps_ait,
            self.residual_hist,
            self.residual_maps_center,
            self.residual_maps_cygnus,
            self.residual_maps_vela,
            self.residual_maps_anti,

            self.offset_map_0, self.offset_map_1,
            self.ridge_systematics_scat,
            self.ridge_systematics_plot,
             
       ]) 

class AllSkyMaps(object):
    """For interactive analysis of the all-sky maps saved by a ResidualMaps run
    """
    def __init__(self, filename = 'residual_maps/f131maps.fits'):
        from astropy.table import Table
        assert os.path.exists(filename)
        self.df=Table.read(filename, hdu=1).to_pandas()
        nside = int(np.sqrt(len(self.df)/12))
        print 'Read file "{}" nside={}, with columns {}'.format(filename, nside, list(self.df.columns))
        self.modelname = os.path.split(os.getcwd())[-1]

        # add locations to the DataFrame
        glon,glat=healpy.pix2ang(nside, range(12*nside**2), lonlat=True)
        glon[glon>180] -= 360
        self.df['glon']=glon
        self.df['glat']=glat
    
    def __repr__(self):
        cols = self.df.columns
        r = 'All-sky maps for model {}\n{:10} {:>12}\n'.format(self.modelname, 'column', 'sum')
        for col in self.df.columns:
            r += '{:10} {:12,}\n'.format(col, int(sum(self.df[col])))
        return r
                
    def __getitem__(self, name):
        """return HParray object """
        return healpix_map.HParray(name, self.df[name])
    
    def ait_plot(self, comp, grid_color='white', **kwargs):
        """ make an AIT map with grid and colorbar
        """
        fig,ax=plt.subplots(figsize=(12,6))
        kw= dict(cmap=None, title=comp)
        kw.update(kwargs)
        t=self[comp].plot(axes=ax, cbtext='counts/pixel',**kw)
        if grid_color: t.grid(color=grid_color)

    def pie(self, ax=None, query=None, colors=None):  
        """ make a pie chart with components
        """ 
        labels = 'galactic isotropic sources'.split()
        df = self.df if query is None else self.df.query(query)
        gal, src, data = [sum(df[c]) for c in 'galactic sources data'.split()  ]
        if not colors: colors = ['#1f77b4', '#8c564b', '#e377c2'] # hack to make same as gardian
        # extract by: colors=np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        sizes = np.array([gal, data-gal-src, src])*100. #iso is all but galactic and sources
        explode= [0,0, 0.1]
        fig, ax = plt.subplots(figsize=(6,6)) if ax is None else (ax.figure, ax)
        ax.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax.axis('equal'); # Equal aspect ratio ensures that pie is drawn as a circle.
        ax.set_title('{} component sizes for 133-MeV band'.format(self.modelname))

class GalacticCorrectionMaps(object):
    """Access to the galactic diffuse corrections
    """
    def __init__(self):
        config=configuration.Configuration('.', quiet=True,postpone=True)
        corr_key = config.diffuse['ring'].get('key', None)
        corr_file = config.diffuse['ring'].get('correction', None)
        if corr_key=='gal':
            ff,pp= analysis_base.load_pickles_from_zip()
            gcorr =  np.array([p['diffuse_normalization']['gal'] for p in pp])
            gc= pd.DataFrame(gcorr, columns=' 0 1 2 3 4 5 6 7'.split())

        elif corr_file is not None:
            gc=pd.read_csv(os.path.expandvars('$FERMI/diffuse/'+corr_file))
        else:
            raise Exception('Correction info not found')
        
        sdirs = map(Band(12).dir, range(1728))
        gc['glat']= map(lambda s:s.b(), sdirs)
        gc['glon']= map(lambda s: s.l(),sdirs)
        gc.loc[gc.glon>180,'glon'] -= 360   
        self.gc = gc
        self.energies = np.logspace(2.125,3.875,8) # band energies

    def hparray(self, nside=64, sigma=2.5):
        """ return a list ot HParray object for plots, etc
        Convert to the specified nside, smoothing by sigma (deg)
        """
        hpa = [(healpix_map.HParray(a, self.gc[a])) for a in '01234567']
        hpb = [healpix_map.HParray('{} MeV'.format(int(energy)), x.getcol(nside=64), sigma=sigma)
                 for x, energy in zip(hpa, self.energies)]   
        return hpb

    def plot_map(self, energy_index=0 ):
        gc0 = self.gc['{}'.format(energy_index)]
        c0 = healpix_map.HParray('gc0', gc0/gc0.mean())
        t=c0.plot(vmin=0.6, vmax=1.4, cmap=plt.get_cmap('coolwarm'),
                title='133 MeV Galactic correction factor relative to {:.2f}'.format(gc0.mean()))
        t.grid(color='grey')
        return plt.gcf()

    def gcplots(self, ii=range(0,5), glon_cut=30,ylim=(0.5,1.6)):
        def gcplot(ct, ax=None, energy_index=0):
            if ax is None:
                fig, ax = plt.subplots(figsize=(8,5))
            energy = 10**(2.125 + 0.25*(energy_index))
            ax.plot(ct.glat, ct['{}'.format(energy_index)], '.');
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

# create new diffuse map modified by 


def apply_patch(inpat, patchfile, outpat):

    def modify_by_patch(self, patch, newfile=None):
        modlayer=[]
        for i,layer in enumerate((self)):
            name = layer.name
            a = layer.vec
            if i < len(patch):
                p = patch[i]
            else:
                p = patch[-1]
            b = p.vec
            print layer.name, ',',
            modlayer.append(hpm.HParray(layer.name, hpm.multiply(a,b)))
        print
        outmap=hpm.HEALPixFITS(modlayer, energies=self.energies)
        if newfile is not None: 
            if newfile[0]!='/':
                newfile = os.path.expandvars('$FERMI/diffuse/')+newfile
            outmap.write(newfile)
        return outmap

    print 'Applying {} to {}'.format(patchfile, inpat)
    patch = diffuse.HealpixCube(patchfile); patch.load()
    for iet in 'front back'.split():
        print iet, ':',
        df = diffuse.HealpixCube(inpat.replace('*', iet)); df.load()
        modify_by_patch(df, patch, outpat.replace('*', iet))  