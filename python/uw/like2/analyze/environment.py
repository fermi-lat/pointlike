"""
Environment plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/environment.py,v 1.20 2018/01/27 15:39:29 burnett Exp $

"""

import os, pickle, types, glob
import numpy as np
import pylab as plt
import pandas as pd
from scipy import integrate, misc, optimize

from skymaps import BinnedPhotonData #,SkyDir

from . import roi_info
from .. import (configuration, diffuse, response)
from ..pipeline import stream
from ..pub import healpix_map as hpm;

class Environment(roi_info.ROIinfo):
    """ Environment plots"""
    
    def setup(self, **kw):
        super(Environment, self).setup(**kw)
        self.plotfolder='environment'
        self.config = configuration.Configuration('.',postpone=True, quiet=True)
        s=[]
        for sname in ['ring', 'isotrop']:
            z=[]
            for i,r in self.df.iterrows():
                j = r['diffuse_names'].index(sname)
                z.append(r['counts']['models'][j][1][:16])
            s.append(z)
        
        #s = [[ x[1]['counts']['models'][modelnumber][1][:16] for x in self.df.iterrows()] for modelnumber in range(2)]
        self.bdf = [pd.DataFrame(y, index=self.df.index) for y in s]
        
        # get diffuse conrrections
        
        self.sindec = np.sin(np.radians(np.asarray(self.df.dec,float)))
        self.singlat = np.sin(np.radians(np.array(self.df.glat, float)))
        try:
            t = [self.config.diffuse['isotrop']['correction'].replace('*',etn)
                    for etn in self.config.event_type_names]
        except:
            t = None
        self.isofiles = t
        u = self.config.diffuse['ring']
        if 'correction' in u.keys():
            self.galfile = os.path.expandvars('$FERMI/diffuse/') + u['correction']
        else:
            self.galfile = None
        try:    
            [self.isofiles_front, self.isofiles_back] = [os.path.expandvars('$FERMI/diffuse/')+ f for f in self.isofiles ]
        except:
            print ("no isofiles found")
            self.isofiles_front=self.isofiles_back=None

    def exposure_plots(self, energy=1000.):
        """ Exposure
        The ratio of the exposure to is mean  for the given energy and event type
        """
        cfg = configuration.Configuration(os.path.expandvars('.'), quiet=True);
        exp = cfg.irfs.exposure(0, energy)                                
        hf = hpm.HPskyfun('front-1000 exp', exp, 64);
        expf = hf.getcol()
        emeanf = expf.mean()
        euw=hpm.HParray('FRONT exposure @ {} MeV / {:.2e}'.format(energy, emeanf), expf/emeanf)
        fig,ax=plt.subplots(figsize=(12,6))
        euw.plot(axes=ax,vmin=0.80,vmax=1.20, title=euw.name,   
               cmap=plt.get_cmap('coolwarm')).grid(color='grey');

        return fig
     
    def get_psf(self, irfname=None, ):
        from uw.like import pypsf, pycaldb
        if irfname is None: irfname=self.config.irf
        cdm = pycaldb.CALDBManager(irf=irfname)
        self.psf_files=cdm.get_psf()
        return pypsf.CALDBPsf(cdm)

    def psf_plot(self, irfname=None, outfile='psf.csv', title=''):
        r"""PSF size
        
        <br>Plots of two science-driven measures of the size of the PSF, compared with the standard 68 percent containment.
        If  $P(\delta)$ is the normalized PSF as a function of deviation $\delta$ for a given energy and conversion type, 
        in principle averaged over incidence angle, but here for normal incidence, the two cases are:
        <ol>
        <li>Discrimination of a point source 
        signal in the presence of a uniform background; the <em>average</em> PSF is the inverse of the solid angle, with radius
        $$\begin{equation}
        1 / \sqrt{\pi \int_0^\pi P(\delta)^2 2\pi \sin\delta \ \mathrm{d}\delta}
        \end{equation}$$
        <li> Measurement of the position of a source. Assuming no background, this requires the expected value for the log of the
        PSF, which is the likelihood, as a function of $\delta$
        $$\begin{equation} 
        w(\delta) = \int_0^\pi  \sin\theta\ \mathrm{d} \theta\ P(\theta) \int_0^{2\pi} \mathrm{d}\phi \ln P\left(\sqrt{\delta^2 -2\delta\  \theta \cos\phi+\theta^2}\right)
        \end{equation}$$
        The resolution is the curvature, or second derivative of this function evaluated at $\delta=0$. The curvature at $\delta=0$ is
        $$\begin{equation} 
        \frac{\partial^2 w(\delta)}{\partial \delta^2} = \pi \int_0^\pi \sin\theta \ \mathrm{d}\theta \frac{P'(\theta)^2}{P(\theta)}
        \end{equation}$$
        where $P'(\theta) = \frac{\partial P(\theta)}{\partial \theta}$. This is consistent with equation (A2) in the 2FGL paper, 
        for the case with no background.
        </li>
        </ol>  
        <p>Of course the first measure is relevant in the background-dominated case, below a GeV or so, 
        the second when there is small background, above a few GeV.    
        <p>
        For a Gaussian PSF, $P(\delta)=\frac{1}{2\pi \sigma^2} \exp(-\frac{\delta^2}{2\sigma^2})$, these values are
        respectively $2\sigma$ and $\sigma$. (The 68 percent is in between, at 1.5 $\sigma$.)
        <br><br>PSF filenames: %(psf_files)s
        """
        psf = self.get_psf(irfname)
 
        def bkg_size(e, ct):
            f2 = lambda delta: psf(e,ct, delta)**2 * 2*np.pi*delta
            return np.degrees(1./np.sqrt(np.pi*integrate.quad(f2, 0, np.inf)[0]))
        
        def loc_size(e, ct):
            func = lambda x : psf(e,ct, x)
            fprime = lambda x : misc.derivative(func, x, dx=0.0001, order=5)
            integrand = lambda rp : rp * fprime(rp)**2/func(rp) * np.pi
            return np.degrees(1/np.sqrt(integrate.quad(integrand, 0, np.radians(5))[0]))
            
        
        egev = np.logspace(-1.+1/8., 2.5+1/8., 3.5*4+1)
        front, back = [[bkg_size(e*1e3,ct) for e in egev] for ct in range(2)]
        floc, bloc = [[loc_size(e*1e3,ct) for e in egev] for ct in range(2)]
        f68,b68  = [[psf.inverse_integral(e*1e3, ct) for e in egev] for ct in range(2)]
        fig,ax = plt.subplots(figsize=(6,6))
        for x, s, label in zip((front, back, floc, bloc, f68, b68),
                            ('-g', 'r', '--g', '--r', ':g', ':r'),
                            ('front bkg', 'back bkg','front loc', 'back loc', 'front 68', 'back 68')):
            ax.plot(egev, x, s, lw=2, label=label)
        
        plt.setp(ax, xlabel='Energy (GeV)', ylabel='PSF size (deg)', xscale='log', yscale='log',
            xlim=(0.1, 100), ylim=(0.02, 8), title=title)
        ax.legend(prop=dict(size=10)); ax.grid()
        #x.set_xticklabels('0.1 1 10 100'.split())
        #ax.set_yticklabels('0.01 0.1 1'.split())
        if outfile is None: return fig
        self.psf_df = pd.DataFrame(dict(front=front, floc=floc, back=back, bloc=bloc,f68=f68,b68=b68), 
                index=egev.round(3))
        self.psf_df.index.name='energy'
        self.psf_df.to_csv(os.path.join(self.plotfolder, outfile))
        print ('wrote file %s' % os.path.join(self.plotfolder, outfile))
        return fig
        
    def isotropic_spectrum(self, other=None):
        """ Isotropic Spectrum from template
        
        The spectrum used to define the isotropic diffuse component.
        <br>Files for front/back: %(idfiles)s
        <br>See also the corrections.
        """
        # look up filenames used to define the isotropic spectrum: either new or old diffuse spec; list or dict
        df=diffuse.diffuse_factory(self.config.diffuse['isotrop'])
        self.idfiles = [x.fullfilename for x in df]
        nf,nb = map(np.loadtxt, self.idfiles)
        df = pd.DataFrame([nf[:,0],nf[:,1],nb[:,1]],
            index='energy front back'.split()).T.query('900<energy<110000')

        fig, axs = plt.subplots(1,2, figsize=(12,5), dpi=50)
        def right(ax):
            ax.plot(df.energy, df.front/df.back, '-o');
            ax.axhline(1.0, color='k')
            plt.setp(ax, xscale='log', xlabel='Energy');ax.grid(True);
            ax.set_title('Isotropic flux front/back ratio', fontsize='small');
        def left(ax):
            ax.plot(df.energy, df.front*df.energy**2, '-g', label='front')
            ax.plot(df.energy, df.back*df.energy**2, '-r', label='back')
            plt.setp(ax, xlabel='Energy', ylabel='flux*e**2', xscale='log')
            ax.set_title('isotropic diffuse spectra', fontsize='small')
            ax.grid(True); ax.legend()
        for f,a in zip((left,right), axs.flatten()): f(a)
        # self.iso_df = pd.DataFrame(dict(front=front*energies**2, back=back*energies**2), index=(energies/1e3).round(3))
        # self.iso_df.index.name='energy'
        # self.iso_df.to_csv(os.path.join(self.plotfolder, 'isotropic.csv'))
        # print ('wrote file %s' % os.path.join(self.plotfolder, 'isotropic.csv'))
        return fig

    def limb_map(self, energy=100):
        """Limb plots
        """
        df=diffuse.diffuse_factory(config.diffuse['limb'])

        [dd.plot_map(energy, scale='linear', cbtext='flux', title='Limb %s, %d MeV'%(name,energy))\
            for dd,name in zip(df, self.config.event_type_names)]
    
    def limb_flux(self, energy=100, ra=0):
        """Limb flux 
        Note assume independentx7ZyIil9vaTFBDx7ZyIil9vaTFBD of RA: this is for RA=0.
        """
        from skymaps import SkyDir
        df=diffuse.diffuse_factory(self.config.diffuse['limb'])
        names = self.config.event_type_names
        sindec = np.linspace(-1,1,501)
        dec = np.degrees(np.arcsin(sindec));
        fig, ax = plt.subplots(figsize=(8,5))
        for i in range(len(df)):
            flux = np.array([df[i](SkyDir(ra,d), energy) for d in dec])
            ax.plot(sindec, flux, lw=2, label=names[i])
        ax.legend(); ax.grid();
        plt.setp(ax, xlabel='sin(Dec)', ylabel='flux', title='Limb flux @ %.0f MeV' %energy)
        return fig
    
    def get_background(self, roi):
        return [t.iloc[roi] for t in self.bdf]
    
    def diffuse_flux(self, rois=[0,888]):
        """Diffuse flux
        Predicted counts for the low latitude and high latitude ROIs.
        """
        fig, ax = plt.subplots(1,1, figsize=(6,6), dpi=150, sharey=True)
        egev = np.array(self.energy)/1e3
        if rois is None: rois = self.rois

        for r in rois:
            gal, iso = self.get_background(r)
            ax.plot(egev, gal, '-D', label='gal %d'%r)
            ax.plot(egev, iso, '--o', label='iso %d'%r)
        plt.setp(ax, xscale='log', xlim=(0.1,300), xlabel='Energy (GeV)',
                     yscale='log', ylim=(1e-1,1e6), ylabel='Diffuse counts/ROI')
        ax.legend(prop=dict(size=10)); ax.grid()
        return fig
        
    def ecliptic_coords(self):
        enp=SkyDir(270,90-23.439281) #ecliptic north pole
        gdir = [SkyDir(l,b, SkyDir.GALACTIC) for l,b in zip(self.df.glon, self.df.glat)]
        edir = np.array([ g.zenithCoords(enp) for g in gdir]); edir[0]
        sinlat = np.sin(np.radians(edir[:,1]))
        lon = edir[:,0]
        lon[lon>180] -= 360
        return lon, sinlat

    def equatorial_coords(self):
        gdir = [SkyDir(l,b, SkyDir.GALACTIC) for l,b in zip(self.df.glon, self.df.glat)]
        lon = np.array([x.ra() for x in gdir])
        lat = np.array([x.dec() for x in gdir])
        sinlat = np.sin(np.radians(lat))
        lon[lon>180] -= 360
        return lon, sinlat
    
    def cartesian_map_array(self, fn, vmin=None, vmax=None, bands=8, title='',cblabel='', 
            ecliptic=False, equatorial=False, nocolorbar=False, cmap=plt.get_cmap('coolwarm')):
        """
        Plot an array of cartesian maps
        
            fn : function object
                fn(iband) returns nside=12 HEALPix array
                has attributes vmin, vmax, title, cblabel
        """
        if vmin is None:vmin=fn.vmin
        if vmax is None: vmax=fn.vmax
        nrows, ncols = ((bands+1)//4, 4 ) if bands>=4 else (1, bands)
        
        fig, axx = plt.subplots(nrows, ncols, figsize=(3+3*ncols,1+3*nrows), sharex=True, sharey=True)
        plt.subplots_adjust(left=0.1, right=0.92, hspace=0.15, wspace=0.01, bottom=0.15)
        if ecliptic:
            lon, sinlat = self.ecliptic_coords()
        elif equatorial:
            lon, sinlat = self.equatorial_coords()
        else:
            lon = self.df.glon
            sinlat = self.singlat
        for iband,energy in enumerate(self.energy[:bands]):
            ax = axx.flatten()[iband] if bands>1 else axx
            scat=self.basic_skyplot(ax, lon, sinlat, fn(iband).clip(vmin,vmax),
                 title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, s=30, edgecolor='none', colorbar=False, labels=False, cmap=cmap)

        fig.text(0.5, 0.95, getattr(fn, 'title', title),  ha='center', size=14)
        if nocolorbar: return fig
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label(getattr(fn, 'cblabel', cblabel))
        fig.text(0.5, 0.025, 'longitude', ha='center', fontsize=14)
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center', fontsize=14)
        return fig

    class GalacticCorrection():
        def __init__(self, env):
            galdict = env.config.diffuse['ring']
            if galdict.get('key', None)=='gal':
                self.x = np.array([r['gal'] for r in env.df.diffuse_normalization])
                self.title = 'Galactic correction from {}'.format(env.skymodel)
            else:
                self.x = response.DiffuseCorrection(galdict['correction']).correction.iloc
                self.title = 'Galactic correction from {}'.format(galdict['correction'])
            self.cblabel='correction factor'
            self.vmin=0.9; self.vmax=1.1
            
        def __call__(self, energy_index,):
            return self.x[:,energy_index]
        def __getitem__(self,energy_index ):
            return self.x[:,energy_index]

    def anom(self):
        class Anom:
            # Functor to return array for plotting anomaly
            def __init__(self, gal):
                self.gal=gal
                self.title='ROI anomaly'
                self.clip = (-.1,.1)
                self.cblabel='Anomaly'
                self.vmin, self.vmax =self.clip
            def anomaly(self, index, energy_band=0):
                g = self.gal(energy_band)
                nbs=hpm.neighbor_pixels(index);
                return g[index]-g[nbs].mean()
            def __call__(self, iband):
                return np.array(map(lambda i: self.anomaly(i,iband), range(1728)))
        return Anom(self.GalacticCorrection(self))

    def galactic_correction_maps(self, vmin=0.8, vmax=1.2):
        """Galactic correction factor
        """
        return self.cartesian_map_array(self.GalacticCorrection(self),vmin=vmin,vmax=vmax)
    
    def galactic_correction_summary(self):
        """Galactic correction offsets
        
        Plot of the mean devition form 1.0, and RMS as the errorbar, of the galatic correction factors, per energy plane.
        """
        t = self.GalacticCorrection(self).x
        dd = dict()
        for i,c in enumerate(t.T):
            dd[i] = dict(offset=c.mean()-1,std=c.std())
        df = pd.DataFrame(dd).T
        fig, ax = plt.subplots()
        ax.errorbar(x=range(8), y=df.offset.values, yerr=df['std'].values, fmt='o')
        ax.axhline(0, color='grey')
        ax.set(xlabel='energy band', ylabel='offset', title='galactic correction offsets')
        return fig

    def gal_extra(self, sinbcut=0.4):
        """Special plots for galactic correction

        Each point is the normalization factor for an ROI in the range of $b$. 
        The shaded area is the range of the "patch" component of the glactic diffuse
        """
        def tplot(ax, glon, corr):
            ax.plot(glon, corr,  '.')
            ax.set(xlim=(180,-180),xlabel='longitude', ylabel='normalization factor');
            ax.axvspan(-90,60, color='orange' , alpha=0.2)
            ax.axhline(1.0, color='lightgrey');
            ax.axvline(0, color='lightgrey')
            ax.set_xticks(np.linspace(-180,180,9));

        t = self.GalacticCorrection(self)
        lon = self.df.glon
        sinlat = self.singlat
        scut = self.singlat>sinbcut; sum(scut)
        fig,axx=plt.subplots(4,2, figsize=(14,14), sharex=True,sharey=True,
                    gridspec_kw=dict(hspace=0.1, wspace=0.1, left=0.05, top=0.95))
        for i, ax in enumerate(axx.flatten()):
            tplot(ax,lon[scut], (t.x[:,i])[scut])
            ax.text(0.04, 0.9, '{:.0f} MeV'.format(self.energy[i]),transform=ax.transAxes )

        fig.suptitle('Normalization factors for sin(b)>{:.1f}'.format(sinbcut))
        return fig

    def galactic_correction_anomaly(self):
        """Galactic correction anomaly

            Maps of the difference of the normalization factor with respect to the average of the four
            nearest neighbors 
        """
        anom = self.anom()
        return self.cartesian_map_array(anom, vmin=-0.1, vmax=0.1,  bands=2);   

    class IsotropicCorrection(object):
        def __init__(self, residual, event_type_name):
            event_type_index = residual.config.event_type_names.index(event_type_name)
            self.x = response.DiffuseCorrection(residual.isofiles[event_type_index])
            self.title='Isotropic correction for %s'% (event_type_name,)
            self.cblabel = 'correction factor'
            self.vmin=0.5; self.vmax=1.5
            
        def __call__(self, energy_index,):
            return self.x[energy_index]

    def isotropic_correction_front(self):
        """Isotropic correction factor for front
        From file %(isofiles_front)s
        """
        return self.cartesian_map_array(self.IsotropicCorrection(self,'front'))

    def isotropic_correction_back(self):
        """Isotropic correction factor for back
        From file %(isofiles_back)s
        """
        return self.cartesian_map_array(self.IsotropicCorrection(self,'back'))

    def isotropic_correction(self):
        """Isotropic correction summary.
        From files  %(isofiles_back)s and  %(isofiles_front)s
        
        <br>While the files are 1728x8 arrays of corrections applied to each ROI and band, only the Back 
        varies for the first two energy bins.
        The first plot, for those back energy bins, I plot the average for |Dec|<30
         """
        isob = self.IsotropicCorrection(self,'back')
        isof = self.IsotropicCorrection(self,'front')
        sindec = np.sin(np.radians(np.array(self.df.dec.values,float)))
        fig, axx = plt.subplots(1,2, figsize=(12,5), sharey=True)
        ax=axx[1]
        for i in range(2):
            ax.plot(sindec, isob(i), '.', label='Energy Bin {}'.format(i));
        ax.set(xlabel='sin(Dec)',  title='Back correction vs. Dec.')

        ax=axx[0]
        for f, name in [(isof, 'Front'), (isob, 'Back')]:
            means = [f(i)[np.abs(sindec)<0.25].mean() for i in range(8)]
            ax.plot(means, 'o', label=name)
        ax.set_title('Correction factor vs Energy Bin')
        ax.set(xlabel='Energy Bin',ylabel='Correction Factor',)

        for ax in axx:
            ax.grid(alpha=0.5);
            ax.axhline(1.0, color='k', ls='--')
            ax.legend()

        return fig

    def load_isofits(self):
        if not os.path.exists('isotropic_fit'): return False
        # process isotrop
        files = sorted(glob.glob('isotropic_fit/*.pickle'))
        if len(files)>0:
            if len(files)<1728:
                msg= "found {} files, expected 1728".format(len(files))
                print (msg)
                raise Exception(msg)
            self.isofits = np.array([pickle.load(open(f)) for f in files]);
            model = '/'.join(os.getcwd().split('/')[-2:])
            streamdf= pd.DataFrame(stream.StreamInfo(model)).T
            snum=streamdf.query('stage=="fitisotropic"').index[-1]
            print ('loaded iso fits, generated by stream {} at {}'.format(snum,streamdf.loc[snum].date ))
            return True
  
    def dmap_info(self, out=None):
        """ formatted table of band contents """
        binfile = self.config.dataset.binfile
        dmap = BinnedPhotonData(binfile)
        print ('File: %s ' %binfile, file=out)
        print ('\n  index    emin      emax  type  nside     photons', file=out)
        total = 0
        def bignum(n):
            t = '%9d' % n
            return '  '+' '.join([t[0:3],t[3:6],t[6:]])
        for i,band in enumerate(dmap):
            fmt = '%5d'+2*'%10d'+2*'%6d'+'%12s'
            print (fmt % (i, round(band.emin()), round(band.emax()), 
                    band.event_class()&15, band.nside(), bignum(band.photons())))
            total += band.photons()
        print ('total%45s'% bignum(total), file=out)
        return dmap

    def correction_plots(self, cc, vmin=0.5, vmax=1.5, title=None, hist=False, start=0, cmap='coolwarm',
            cbtext='correction factor', **kwargs):
        from  matplotlib import patches

        nrows = cc.shape[1]/4
        #assert cc.shape[1]==8, 'Found shape {}'.format(cc.shape)

        if hist:
            hkw=dict(bins=np.linspace(vmin,vmax, 21), lw=1, histtype='step')
            fig,axx = plt.subplots(nrows,4, figsize=(14,3*nrows+1), sharex=True, sharey=False)
            plt.subplots_adjust(wspace=0.3, hspace=0.15)
        else:
            fig, axx = plt.subplots(nrows,4, figsize=(12,3*nrows), sharex=True, sharey=True)
            plt.subplots_adjust(left=0.10, wspace=0.1, hspace=0.15,right=0.92, top=0.90)
            
        for i,ax in enumerate(axx.flatten()):
            if i<start:
                ax.set_visible(False)
                continue
            if hist:
                h = np.array(cc[:,i],float)
                ax.hist(h.clip(vmin, vmax),  **hkw)
                ax.axvline(1.0, color='grey', ls='--')
                mypatch= patches.Patch(fill=False,lw=0, facecolor='none', 
                    label='{:4.1f} {:4.1f}'.format(100*(h.mean()-1),100*h.std()),)
                ax.legend(handles=[mypatch], facecolor='none', edgecolor='none')
            else:
                t,scat=self.skyplot(cc[:,i],ax=ax, vmin=vmin, vmax=vmax, 
                                    title='{:0f}'.format(self.energy[i]),
                        cmap=plt.get_cmap(cmap), colorbar=False, labels=False, **kwargs)
            ax.set_title('{:.0f} MeV'.format(self.energy[i]), fontsize=12)
            
        if not hist: 
            cbax = fig.add_axes((0.94, 0.15, 0.015, 0.7) )
            fig.colorbar(scat, cbax, orientation='vertical').set_label(cbtext, fontsize=12)
        fig.suptitle(title, fontsize=16)
        return fig

    def count_difference(self, vmin=-4, vmax=4, cmap='jet', get_data=False):
        """Count differences

        For each ROI and energy band, this is the differnce in counts implied by the galactic and isotropic factors,
        relative to the isotropic. (This latter helps normalize the different energies)
        """
        galcnt = np.array([ct['models'][0][1][:8] for ct in self.df.counts]) 
        isocnt = np.array([ct['models'][1][1][:8] for ct in self.df.counts]) 
        galcorr=self.GalacticCorrection(self).x
        x = [np.array([r['iso'][fb] for r in self.df.diffuse_normalization])for fb in 'front back'.split()]
        isocorr = 0.5*(x[0]+x[1])

        t  =(galcnt * (galcorr-1) + isocnt * (isocorr-1))/isocnt
        if get_data: #for interactive analysis
            return t
        return self.correction_plots(t, vmin=vmin, vmax=vmax, title ='count difference relative to isotropic',
                            cbtext='relative to isotropic', cmap=cmap);

    def all_plots(self, **kw):
        self.runfigures([
            #self.psf_plot, 
            self.exposure_plots, 
            self.isotropic_spectrum,self.diffuse_flux, #self.limb_flux,
            self.count_difference,
            self.galactic_correction_maps,
            self.galactic_correction_summary, 
            self.galactic_correction_anomaly,
            #self.gal_extra,
            #self.isotropic_correction,
            #self.isotropic_correction_front, self.isotropic_correction_back,
            ])
