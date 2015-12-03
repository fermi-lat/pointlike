"""
Residual plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/residuals.py,v 1.14 2015/07/24 17:56:02 burnett Exp $

"""

import os, glob, pickle
import numpy as np
import pylab as plt
import pandas as pd
import matplotlib.gridspec as gridspec
from skymaps import SkyDir 
from . import (roi_info,  analysis_base)
from .. import (tools, configuration, response)

class Residuals(roi_info.ROIinfo): 
    """<b>Residual plots</b>
 <p>
 Analysis of results from a measurement (using the uwpipeline task 'residuals') of the nomalization likelihood
 functions performed for every source and every energy band. 
 <p>
 That pipeline analysis makes individual likelihood fits to the normalization of all the components of the model for each ROI, 
 for each energy band, and separately for all event types. Event types are assumed to be front and back for now.
 <p>
 Each such fit is actually a measurement of the likelihood function, for which the maximum, plus and minus one sigma points, 
 TS, and pull are recorded. The pull is the deviation of the fit value from the nominal, usually 1.0, in sigma units. 
 <p>
 Since the galactic diffuse and isotropic, separately for front and back have ROI-specific correction factors, summary
 plots of these are shown.
 
    """ 
    require='residuals.zip'
    def setup(self, **kwargs):
        super(Residuals, self).setup(**kwargs)

        self.plotfolder = 'residuals'
        files, self.pkls = analysis_base.load_pickles_from_zip('residuals.zip')
        if len(files)!=1728:
            print 'Found %s files, expected 1727' % len(files)
            nn = set([int(x[15:19]) for x in files]); 
            print 'mising: %s' % set(range(1728)).difference(nn)
            raise
        self.energy = self.pkls[0]['ring']['all'].index
        self.sindec = np.sin(np.radians(np.asarray(self.df.dec,float)))
        self.singlat = np.sin(np.radians(np.array(self.df.glat, float)))
        # expect to be in a skymodel folder
        self.config = configuration.Configuration('.', postpone=True, quiet=True)
        try:
            t = [self.config.diffuse['isotrop']['correction'].replace('*',etn)
                for etn in self.config.event_type_names]
        except:
            t = None
        self.isofiles = t
        u = self.config.diffuse['ring']
        if 'correction' in u.keys():
            self.galfile = u['correction']
        else:
            self.galfile = None


    def resid_array(self, source_name, column_name, event_type='all'):
        """ extract a component from the array of residuals
        
        source_name : str
            for now, the name of a global source
        column_name :
            one of the column names
        returns a 2-d array, shape (1728,14)
        
        """
        empty = [np.nan]*len(self.energy)
        r= np.array([list(p[source_name][event_type][column_name]) if source_name in p else empty for p in self.pkls])
        assert r.shape==(1728,14), 'Failed shape requirement'
        return r
        
    class ResidualArray(object):
        """ manage access for plotting to component of the residual analysis"""
        
        def __init__(self, residual, column_name, source_name, event_type='all', vmin=-5, vmax=5):
            self.vmin=vmin; self.vmax = vmax
            self.cblabel='normalized residual' if column_name=='pull' else ''
            self.title='%s for %s, %s events'%(column_name,source_name, event_type)
            empty = [np.nan]*len(residual.energy)
            r= np.array([list(p[source_name][event_type][column_name]) if source_name in p else empty for p in residual.pkls])
            assert r.shape==(1728,14), 'Failed shape requirement'
            
            self.resid = r #residual.resid_array(source_name, colname, event_type=event_type)
            
        def __call__(self, iband):
            return self.resid[:,iband]
    
    
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
    
    def cartesian_map_array(self, fn, vmin=None, vmax=None, bands=8, 
            ecliptic=False, equatorial=False, nocolorbar=False):
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
        plt.subplots_adjust(right=0.9, hspace=0.15, wspace=0.1)
        if ecliptic:
            lon, sinlat = self.ecliptic_coords()
        elif equatorial:
            lon, sinlat = self.equatorial_coords()
        else:
            lon = self.df.glon
            sinlat = self.singlat
        for iband,energy in enumerate(self.energy[:bands]):
            ax = axx.flatten()[iband]
            scat=self.basic_skyplot(ax, lon, sinlat, fn(iband).clip(vmin,vmax),
                 title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, s=30, edgecolor='none', colorbar=False, labels=False)
        fig.text(0.5, 0.95, fn.title,  ha='center', size=12)
        if nocolorbar: return fig
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label(fn.cblabel)
        fig.text(0.5, 0.05, 'longitude', ha='center')
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center')
        return fig

    def update_correction(self, vin, vout=None):
        """ update the Galactic Diffuse correction factor array with new residuals
                NOTE: this is intended to be done by hand

        """
        if vout is None: vout=self.skymodel
        # get current residual array, replace any nan's with 1.0
        t = self.resid_array('ring', 'maxl')
        ra = t[:,:8] # only upto 10 GeV
        ra[~np.isfinite(ra)]=1.0
        # read in input correction array
        infile = os.path.expandvars('$FERMI/diffuse/galactic_correction_%s.csv'%vin)
        assert os.path.exists(infile), 'File %s not found' %infile
        cv_old = pd.read_csv(infile, index_col=0)
        # multiply input corrections by residuals ane write it out
        cv_new = ra * cv_old
        outfile = os.path.expandvars('$FERMI/diffuse/galactic_correction_%s.csv'%vout)
        cv_new.to_csv(outfile)
        print 'wrote new diffuse correction file %s' % outfile
    def update_iso_correction(self, vin, vout=None):
        """ Update the isotropic diffuse correction factors with current residuals 
        vin: string
            skymodel used to generate current correction.
        """
        if vout is None: 
            vout=self.skymodel
            print 'setting vout:', vout
        for et in ('front', 'back'):
            infile = os.path.expandvars('$FERMI/diffuse/isotropic_correction_%s_%s.csv'% (et, vin))
            assert os.path.exists(infile), 'File %s not found' %infile
            # get current residual array, replace any nan's with 1.0
            t = self.resid_array('isotrop', 'maxl', et)
            ra = t[:,:8] # only upto 10 GeV
            ra[ra>2.0]=2.0
            ra[ra<0.5]=0.5
            ra[~np.isfinite(ra)] = 1.0
            cv_old = pd.read_csv(infile, index_col=0)
            # multiply input corrections by residuals ane write it out
            cv_new = ra * cv_old
            outfile = os.path.expandvars('$FERMI/diffuse/isotropic_correction_%s_%s.csv'% (et,vout))
            cv_new.to_csv(outfile)
            print 'wrote new iso diffuse correction file %s' % outfile

    def norm_plot(self, name='isotrop', ax=None, ylim=(0.5,1.5)):
        """Isotropic Normalization vs Dec
        Only the isotropic component is allowed to vary; this is the resulting value.
        """
        lnorms =np.array([m[0] if m is not None else np.nan for m in self.diffuse_models(name)])
        high = np.array(np.abs(self.df.glat)>10, bool)
        if ax is None:
            fig,ax=plt.subplots(1,1, figsize=(10,5))
        else: fig=ax.figure
        ax.plot(self.sindec,      lnorms.clip(*ylim),  '.r' , label='|b|<10')
        ax.plot(self.sindec[high],lnorms[high].clip(*ylim),  'og', label='|b|>10')
        plt.setp(ax, ylim=ylim, xlim=(-1,1), xlabel='sin(Dec)',
            ylabel='normalization factor', title='%s normalization vs sin(Dec)'%name)
        ax.grid(); ax.legend(prop=dict(size=10))
        return fig

    def front_back_ridge(self, xlim=(0.8, 1.2)):
        """front/back galactic residual ratio on ridge
        Ridge is within 10 degrees in latitude, 60 degrees in longitude.
        <br>Top three rows: histograms
        <br>Bottom plot: average
        """
        fb = [self.resid_array('ring', 'maxl', et) for et in ['front', 'back']]
        fbr = fb[0]/fb[1]; 
        ridge = np.asarray((np.abs(self.df.glat)<10) & (np.abs(self.df.glon)<60)); 
        fbrr = fbr[ridge]

        fig =plt.figure(figsize=(12,9))
        gs1 = gridspec.GridSpec(3,4, bottom=0.35)
        axt = []; first=True
        for i in range(3):
            for j in range(4):
                axt.append(plt.subplot(gs1[i,j], sharex=None if first else axt[0] ))
                first=False
        
        gs2 = gridspec.GridSpec(1,1, top=0.3, right=0.6)
        axb = plt.subplot(gs2[0,0])
        means = []
        for i,ax in enumerate(axt):
            u = fbrr[:,i]
            uok = u[np.isfinite(u)]; means.append(uok.mean());
            ax.hist(u.clip(*xlim), np.linspace(xlim[0], xlim[1], 51));
            ax.axvline(1.0, color='b', ls='--')
            ax.text(0.1, 0.9, '%.0f MeV' % self.energy[i], size=10,  transform = ax.transAxes)
        plt.setp(axt[0] , xlim=xlim)
   
        ax=axb
        ax.semilogx(self.energy[:12], means, 'o-')
        ax.axhline(1.0, color='b', ls='--')
        ax.grid()
        plt.setp(ax, ylim=(0.9,1.10), xlabel='Energy [MeV]', ylabel='front/back ratio')
        return fig
        
    def front_back_strong(self):
        """Front/back ratio for strongest sources
        """
        sources = glob.glob('sources_*.csv')
        assert len(sources)>0, 'no sources_*.csv files fouind'
        filename = sources[0]
        assert os.path.exists(filename), 'sources csv file not found'
        sdf = pd.read_csv(filename, index_col=0)
        t =sdf.sort_index(by='ts')
        strong_names = t.ts[-4:].index
        print 'selected strong names: %s' % list(strong_names)
        roinames = [sdf.ix[sn]['roiname'] for sn in strong_names]
        rois = map(lambda s: int(s[-4:]), roinames)
        
        rdict = {}; edict={}
        for roi,name in zip(rois, strong_names):
            t = self.pkls[roi][name]
            u =np.array([np.array(t[et]['flux']) for et in ('front','back')]); 
            du = np.array([ np.array(t[et]['uflux']-t[et]['flux']) for et in ('front', 'back')]) 
            rdict[name] = u[0,:8]/u[1,:8]
            edict[name] = np.sqrt((du[0,:8]/u[0,:8])**2 +
                                  (du[1,:8]/u[1,:8])**2 )

        fig, axx = plt.subplots(2,2, figsize=(10,10), sharex=True, sharey=True)
        for (key, values), errors, ax  in zip(rdict.items(),edict.values(), axx.flatten()):
            ax.errorbar(x=self.energy[:8],y=values,yerr=errors, fmt='o')
            ax.grid()
            ax.axhline(1.0, color='k', ls ='--')
            plt.setp(ax, xlabel='Energy [MeV]', xscale='log', title=key, ylim=(0.85,1.15));
        return fig
        
    def isotropic_hists(self, bmin=10, xlim=(0.4,1.6), etnames = ('front', 'back', 'all') ):
        """Isotropic normalization hists
        Maximum likelihood values for normalization for front (green) and back (red), all (black)
        """
        maxl = [self.resid_array('isotrop', 'maxl', et) for et in etnames ]
        gcut = abs(self.df.glat)>bmin; 
        fig,axx = plt.subplots(2,4, figsize=(12,8), sharex=True)
        for i, ax in enumerate(axx.flatten()):
            z = [maxl[j][gcut,i].clip(*xlim) for j in range(len(etnames))]
            ax.hist(z, np.linspace( *xlim), histtype='step', color=('g','r', 'k'), label=etnames);
            ax.text(0.1,0.9, '%d'%self.energy[i],  transform = ax.transAxes)
            ax.axvline(1.0, ls = '--')
        fig.text(0.5, 0.05, 'normalization factor', ha='center')
        return fig

        
    def pull_maps_ring(self):
        """Pull plots for galactic diffuse
        """
        return self.cartesian_map_array( self.ResidualArray(self, 'pull', 'ring')); 
        
    def pull_maps_isotrop_front(self):
        """Pull plots for isotropic front
        """
        return self.cartesian_map_array( self.ResidualArray(self, 'pull', 'isotrop', 'front')); 
        
    def pull_maps_isotrop_back(self):
        """Pull plots for isotropic back
        """
        return self.cartesian_map_array( self.ResidualArray(self, 'pull', 'isotrop', 'back')); 

    def maxl_plots_isotrop_back(self):
        """Max Likelihood for Isotropic back
        """
        return self.cartesian_map_array( self.ResidualArray(self, 'maxl', 'isotrop', 'back', vmin=0.9, vmax=1.1), bands=4); 
        
    def maxl_plots_isotrop_front(self):
        """Max Likelihood for Isotropic front
        """
        return self.cartesian_map_array( self.ResidualArray(self, 'maxl', 'isotrop', 'front', vmin=0.9, vmax=1.1), bands=4); 
        
    def maxl_map_ring(self, vmin=0.9, vmax=1.1):
        """Max likelihood for ring"""
        return self.cartesian_map_array( self.ResidualArray(self, 'maxl', 'ring', vmin=vmin, vmax=vmax), ); 
    
    class GalacticCorrection():
        def __init__(self, residual):
            self.x = response.DiffuseCorrection(residual.config.diffuse['ring']['correction'])
            self.title = 'Galactic correction'
            self.cblabel='corection factor'
            self.vmin=0.9; self.vmax=1.1
            
        def __call__(self, energy_index,):
            return self.x[energy_index]

    def galactic_correction(self):
        """Galactic correction factor"""
        return self.cartesian_map_array(self.GalacticCorrection(self))
    
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
        """Isotropic correction factor for front"""
        return self.cartesian_map_array(self.IsotropicCorrection(self,'front'))

    def isotropic_correction_back(self):
        """Isotropic correction factor for back"""
        return self.cartesian_map_array(self.IsotropicCorrection(self,'back'))

    def plot_region_averages(self):
        """Special isotropic correction check
        This shows the average correction factor as a function of energy for two 
        high latitude regions, and all
        """
        lon, sinb = np.array(self.df.glon), self.singlat
        regiona=(-30>lon) & (lon>-80) & (np.abs(sinb-0.8)<0.1)
        regionb=(-100>lon) & (lon>-150) & (np.abs(sinb+0.8)<0.1)

        fig,axx = plt.subplots(1,2, figsize=(12,6), sharex=True, sharey=True)
        def plotit(ax, name):
            iso = self.IsotropicCorrection(self,name)
            fr = iso.x.correction
            e = self.energy[:8]
            total, ra, rb =fr.mean(), fr[regiona].mean(), fr[regionb].mean()
            ax.plot(e,total, 'o--b', label='all');
            ax.plot(e, ra, 'o--r', label='region A')
            ax.plot(e, rb, 'o--g', label='region B')
            ax.axhline(1.0, color='k')
            plt.setp(ax, title=name, xlabel='Energy {MeV]',xscale='log', ylabel='Factor')
            ax.legend(loc='lower left'); ax.grid()
        map( plotit, axx, ('front', 'back'))
        return fig
  
    def all_plots(self):
        self.runfigures([
            self.pull_maps_ring, self.pull_maps_isotrop_front, self.pull_maps_isotrop_back, 
            #self.norm_plot,
            #self.isotropic_hists,
            self.maxl_plots_isotrop_front, self.maxl_plots_isotrop_back,
            self.maxl_map_ring,
            #self.isotropic_correction_ait,
            self.galactic_correction, self.isotropic_correction_front, self.isotropic_correction_back,
            self.front_back_ridge, self.front_back_strong,
            self.plot_region_averages, 
            ])
