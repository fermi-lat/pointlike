"""
Diffuse fitting analysis

"""
import os, glob, pickle, healpy
import numpy as np
import pylab as plt
import pandas as pd
from  matplotlib import (patches, )
from skymaps import SkyDir, Band 
from . import (roi_info,  analysis_base)
from .. import ( diffuse, maps,)
from ..pipeline import stream

from uw.utilities import healpix_map
hpm=healpix_map

class DiffuseFits(roi_info.ROIinfo):
    """<b>Diffuse fit plots</b>
    <p>
    """
    def setup(self,**kwargs):
        super(DiffuseFits, self).setup()
        self.plotfolder = 'diffuse_fits'
        if not hasattr(self, 'energy'):
            self.energy = np.logspace(2.125, 5.875, 16)

    def summary(self):
        """Summary
        <pre>%(logstream)s</pre>
        """
        self.galfits=self.isofits=None
        galdict = self.config.diffuse['ring'] 
        self.galfile = galdict['filename'].split('/')[-1]
        if galdict.get('key', None) == 'gal' or 'correction' not in galdict:
            self.galcorr = None
        else:
            self.galcorr =galdict['correction']
        # Find streams 
        model = '/'.join(os.getcwd().split('/')[-2:])
        streamdf= pd.DataFrame(stream.StreamInfo(model)).T
        self.startlog()
        if os.path.exists('isotropic_fit'):
            # process isotrop
            files = sorted(glob.glob('isotropic_fit/*.pickle'))
            if len(files)>0:
                if len(files)<1728:
                    msg= "found {} files, expected 1728".format(len(files))
                    print (msg)
                    raise Exception(msg)
                self.isofits = np.array([pickle.load(open(f)) for f in files]);
                snum=streamdf.query('stage=="fitisotropic"').index[-1]
                print ('loaded iso fits, generated by stream {} at {}'.format(snum,streamdf.loc[snum].date ))
     
        if os.path.exists('galactic_fit'):
            files = sorted(glob.glob('galactic_fit/*.pickle'))
            if len(files)>0:
                if len(files)<1728:
                    msg= "found {} files, expected 1728".format(len(files))
                    print (msg)
                    ids= map(lambda f: int(f.split('.')[-2][-4:]), files);
                    print (np.array(sorted(list(set(range(1728)).difference(set(ids))))))
                    print ('trying to continue...')
                    #raise Exception(msg)

                self.galfits = np.array([pickle.load(open(f)) for f in files]); 
                snum=streamdf.query('stage=="fitgalactic" or stage=="postfitgalactic"').index[-1]
                print ('loaded gal fits, generated by stream {} at {}'.format(snum,streamdf.loc[snum].date ))
        
        # get current isotropic template and values at energies
        self.iso=diffuse.diffuse_factory(self.config.diffuse['isotrop'])
        print ('isotropic:\n {}'.format(self.iso))
        self.isoflux = np.array([np.array(map(lambda e: self.iso[i](None, e), self.energy)) for i in range(2)])
        
        # and current galactic 
        self.gal=diffuse.diffuse_factory(self.config.diffuse['ring'])
        print ('galactic:\n{}'.format(self.gal))
        self.logstream= self.stoplog()

    def correction_plots(self, cc, vmin=0.5, vmax=1.5, title=None, hist=False, start=0):
        if isinstance(cc, pd.DataFrame):
            cc = cc.as_matrix()
        nrows = cc.shape[1]/4
        #assert cc.shape[1]==8, 'Found shape {}'.format(cc.shape)
        if title is None:
            title = 'Galactic adjustments to: {}'.format(self.galcorr)
        if hist:
            hkw=dict(bins=np.linspace(vmin,vmax, 21), lw=1, histtype='step')
            fig,axx = plt.subplots(nrows,4, figsize=(14,3*nrows+1), sharex=True, sharey=False)
            plt.subplots_adjust(wspace=0.3)
        else:
            fig, axx = plt.subplots(nrows,4, figsize=(12,3*nrows), sharex=True, sharey=True)
            plt.subplots_adjust(left=0.10, wspace=0.1, hspace=0.1,right=0.92, top=0.92)
        for i,ax in enumerate(axx.flatten()):
            if i<start:
                ax.set_visible(False)
                continue
            if hist:
                h = cc[:,i]
                ax.hist(h.clip(vmin, vmax),  **hkw)
                ax.axvline(1.0, color='grey', ls='--')
                mypatch= patches.Patch(fill=False,lw=0, facecolor='none', 
                    label='{:4.1f} {:4.1f}'.format(100*(h.mean()-1),100*h.std()),)
                ax.legend(handles=[mypatch], facecolor='none', edgecolor='none')
            else:
                t,scat=self.skyplot(cc[:,i],ax=ax, vmin=vmin, vmax=vmax, title='{:0f}'.format(self.energy[i]),
                        cmap=plt.get_cmap('coolwarm'), colorbar=False,labels=False);
            ax.set_title('{:.0f} MeV'.format(self.energy[i]))
        if not hist: 
            cbax = fig.add_axes((0.94, 0.15, 0.015, 0.7) )
            fig.colorbar(scat, cbax, orientation='vertical').set_label('correction factor', fontsize=12)
        fig.suptitle(title, fontsize=14)
        return fig

    def corr_plot(self, c, ax=None, vmin=0.5, vmax=1.5, title=None, colorbar=True,cmap='coolwarm', **scatkw):
        """SkyPlot of fit or correction factors
        """
        assert c.shape==(1728,), 'Found shape {}'.format(c.shape)
        if ax is None:
            fig, ax = plt.subplots(figsize=(6,6))
        else: fig=ax.figure
        t,scat=self.skyplot(c,ax=ax, vmin=vmin, vmax=vmax,
                        cmap=cmap, colorbar=colorbar,labels=True, **scatkw)
        if title is not None:
            ax.set_title(title, fontsize=14)

    def galactic_fit_maps(self):
        """Galactic correction fits
        Results of normalization fits to adjust level of galactic diffuse flux
        """
        if self.galfits is None: return
        return self.correction_plots(self.galfits, title='Fit to {}'.format(self.galfile),vmin=0.98,vmax=1.02)
 
    def galactic_fit_hists(self):
        """Galactic correction fits
        """
        if self.galfits is None: return
        return self.correction_plots( self.galfits, title='Fit to {}'.format(self.galfile),vmin=0.98,vmax=1.02, hist=True)

    def write_spectral_cube(self):
        gf = self.galfits
        scube = [hpm.HParray('',gf[:,i] ).getcol(nside=128) for i in range(8)]
        sm = hpm.HEALPixSkymap(np.array(scube).T, self.energy[:8])
        sm.write(self.galfile.replace('.fits', '_corrections.fits'))

    
    # def all_plots(self):
    #      self.runfigures([
    #          self.summary,
    #          self.galactic_fit_maps,
    #          self.galactic_fit_hists,
    #     ])

class CombinedGalIso(roi_info.ROIinfo):
    """Combined galactic and isotropic fit analysis
        
        <p> %(logstream)s
        <p>This is an analysis of the 1728x8 separate fits to normalization factors for galactic and isotropic.
        Each fit optimizes the likelihood for the combined <i>front</i> and <i>back</i> data (except only
         <i>front</i> below 316 MeV) with respect to the two normalization factors. Such a fit can fail if the 
         components are too correlated; in this case the galactic is set to unity, and the isotropic modified.
        The same happens if the two-component fit is successful, but the correlation coeficient is less than -0.95.
        
        <p>The sources, from uw8608, are fixed.

    """
    def setup(self, **kwargs):
        """<b>Diffuse fit plots</b>
        <p>
        """
        super(CombinedGalIso, self).setup()
        self.plotfolder = 'diffuse_fits'
        # load files from "fitdiffuse run"
        files = sorted(glob.glob('diffuse_fit/*.pickle'))
        assert len(files)>0, 'no files found'
        if len(files)<1728:
            msg= "found {} files, expected 1728".format(len(files))
            print (msg)
            raise Exception(msg)
        self.fit_array = np.array([pd.read_pickle(f).values for f in files]);
        self.startlog()
        # Find streams 
        model = '/'.join(os.getcwd().split('/')[-2:])
        streamdf= pd.DataFrame(stream.StreamInfo(model)).T
        snum=streamdf.query('stage=="fitdiffuse"').index[-1]
        print ('Analysis of combined galactic and isotropic fits, generated by stream named "fitdiffuse" number {} at {}'.format(snum,streamdf.loc[snum].date ))
        print ('<br>galactic iemfile: {}'.format(self.config.diffuse['ring']['iemfile']))
        print ('<br>isotropic spectral files: {}'.format(self.config.diffuse['isotrop']['filename']))
        self.logstream= self.stoplog()

    def corr_scat(self, xlim=(0.5,1.5), ylim=(0.5,2.0), title='galactic vs isotropic'):
        """Correlation of galactic and isotropic fits

        """
        cc= self.fit_array
        nrows = cc.shape[1]/4
        fig,axx = plt.subplots(nrows,4, figsize=(14,3*nrows+1), sharex=True, sharey=True,
                            gridspec_kw=dict(wspace=0.1, hspace=0.15, left=0.07))
        lolat, label = np.abs(self.df.glon)<5 , '|b|<5'
        for i,ax in enumerate(axx.flatten()):
            ax.axhline(1.0, color='lightgrey')
            ax.axvline(1.0, color='lightgrey')
            x = cc[:,i,0].astype(float).clip(*xlim) 
            y = cc[:,i,1].astype(float).clip(*ylim)
            ax.plot(x ,       y,      '.')
            ax.plot(x[lolat], y[lolat],'.r', label=label)
            ax.set_title('{:.0f}'.format(self.energy[i]),fontsize=12)
            ax.legend()
        fig.suptitle(title)
        fig.text(0.5, 0.025, 'galactic factor',ha='center' )
        fig.text(0.01, 0.5, 'isotropic factor', va='center', rotation=90)
        return fig

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

    def galactic_fit_maps(self):
        """Galactic fit maps
        """
        return self.correction_plots( self.fit_array[:,:,0], title='galactic factor', vmin=0.75, vmax=1.25);

    def galactic_fit_hists(self):
        """Galactic fit histograms
        """
        return self.correction_plots( self.fit_array[:,:,0], title='galactic factor', hist=True, vmin=0.5, vmax=1.5);
    
    def isotropic_fit_maps(self):
        """Isotropic fit maps
        """
        return self.correction_plots( self.fit_array[:,:,1], title='isotropic factor',
            cmap='YlOrRd', vmin=0, vmax=5);

    def isotropic_fit_hists(self):
        """Isotropic fit histograms
        """
        return self.correction_plots( self.fit_array[:,:,1], title='isotropic factor', hist=True, 
            vmin=0.5, vmax=2.5);
    
    def count_difference(self, vmin=-4, vmax=4, cmap='jet'):
        """Count differences

        For each ROI and energy band, this is the differnce in counts implied by the galactic and isotropic factors,
        relative to the isotropic. (This latter helps normalize the different energies)
        """
        galcnt = np.array([ct['models'][0][1][:8] for ct in self.df.counts]) 
        isocnt = np.array([ct['models'][1][1][:8] for ct in self.df.counts]) 
        galcorr, isocorr = self.fit_array[:,:,0], self.fit_array[:,:,1]

        t  =(galcnt * (galcorr-1) + isocnt * (isocorr-1))/isocnt

        return self.correction_plots(t, vmin=vmin, vmax=vmax, title ='count difference relative to isotropic',
                            cbtext='relative to isotropic', cmap=cmap);
    def fit_quality(self):
        """Fit quality

        Quality, a chi-squared measure of the goodness of each fit. Expect this to be zero
        """
        return self.correction_plots( self.fit_array[:,:,3],cmap=plt.get_cmap('YlOrRd'), 
                title='fit quality', cbtext='quality', vmin=0, vmax=1)

    def adjustment_maps(self, sigma=1.5, vmin=0.95, vmax=1.05):
        """Adjustment maps

        Maps of the most recent fit, smoothed by 1.5 deg.
        """

        files = sorted(glob.glob('diffuse_fit_maps/*.pickle'))
        assert len(files)>0, 'no files found'
        if len(files)<1728:
            msg= "found {} files, expected 1728".format(len(files))
            print (msg)
            raise Exception(msg)

        nside=64
        cube = np.zeros([12*nside**2,8])
        index_table = healpix_map.make_index_table(12, nside)
        for i12,file in enumerate(files):
            indeces = index_table[i12]
            cube[indeces,:] = pickle.load(open(file)).T
        assert((cube==0).sum()==0), 'Fail to fill'
        hpm.HParray('test', cube[:,0]).plot(cmap='coolwarm', vmin=0.8, vmax=1.2)
        hpcube = [hpm.HParray('layer{}'.format(i), cube[:,i], sigma=1.5) for i in range(8)]
        self.hpfits = hpm.HEALPixFITS(hpcube)
        skymodel='uw9020'
        hpm.multi_ait(self.hpfits, vmin=0.95, vmax=1.05, cmap='coolwarm', 
            grid_color='white', cblabel='factor',
            title='Smoothed {} galactic diffuse normalization factors'.format(skymodel));

        return plt.gcf()


    def all_plots(self):
         self.runfigures([
            self.galactic_fit_maps,
            self.galactic_fit_hists,
            self.isotropic_fit_maps,
            self.isotropic_fit_hists,
            self.corr_scat,
            self.count_difference,
            self.fit_quality,
            self.adjustment_maps,
        ])

def update_correction(self):
    """
    """
    diffuse_dir = os.path.expandvars('$FERMI/diffuse/')
    if self.galcorr is None:
        outfile=diffuse_dir+self.galfile.replace('.fits', '_corr.csv')
    else:
        i = self.galcorr.find('_corr')
        assert i>0, 'Check galcorr: expected to find "_corr"'
        #corr_version=
    #pd.DataFrame(self.galfits).to_csv(outfile)
    print ('wrote file {}'.format(outfile))


def get_diffuse_norm():
    import pickle,glob
    pf = sorted(glob.glob('pickle/HP12*.pickle'))
    assert len(pf)==1728
    pd = [pickle.load(open(f)) for f in pf]
    r = np.array([p['diffuse_normalization']['gal'] for p in pd])
    return r

def clear_diffuse_norm():
    import pickle,glob
    pf = sorted(glob.glob('pickle/HP12*.pickle'))
    assert len(pf)==1728
    for f in pf:
        p = pickle.load(open(f))
        p['diffuse_normalization']['gal']= np.ones(8)
        pickle.dump(p, open(f, 'w'))

def check_bubble_maps(cubefile):
    if not os.path.exists(cubefile):

        bubble = [diffuse.HealpixCube(f) for f in 
                ['../uw8600/bubbles.fits',    '../uw8604/bubbles_v2.fits',
                '../uw8605/bubbles_v3.fits', '../uw8606/bubbles_v4.fits']]
        for b in bubble:
            b.load()

        # multiply them all together
        energies = bubble[0].energies
        pr = np.array([b.vec for b in  bubble[0]])

        for bb in bubble[1:]:
            pr *= np.array([b.vec for b  in bb])   

        t = [healpix_map.HParray('{:.0f} MeV'.format(energy), v) for energy,v in zip(energies, pr)]
        bcube = healpix_map.HEALPixFITS(t, energies=energies);
        healpix_map.multi_ait(t, vmin=0.8, vmax=2.0, cmap=plt.get_cmap('jet'), 
                    grid_color='grey', cblabel='ratio to diffuse', title='bubbles correction to diffuse');

        bcube.write(cubefile)
    else:
        print ('Using existing file {}'.format(cubefile))
    

class Polyfit(object):
    """ Manage a log parabola fit to every pixel"""
    def __init__(self, cubefile, sigsfile, start=0, stop=8, deg=2):
        
        check_bubble_maps(cubefile)
        m = diffuse.HealpixCube(cubefile); m.load()
        msig = diffuse.HealpixCube(sigsfile); msig.load()
        meas = np.array([m[i].vec for i in range(8)])
        sig  = np.array([msig[i].vec for i in range(8)])

        self.planes = np.array(range(start,stop)) # plane numbers
        self.values = meas[start:,:]
        weights = 100./sig[start:,:] #from percent
        self.wtmean = weights.mean(axis=1)

        self.fit, self.residuals, self.rank, self.svals, self.rcond =\
            np.polyfit(self.planes,self.values, deg=deg, full=True, w=self.wtmean)
            
        labels= 'intercept slope curvature'.split()   
        self.hpfit=[healpix_map.HParray(labels[deg-i], self.fit[i,:]) for i in range(deg,-1,-1)]
        
    def __getitem__(self, i):
        return self.fit[i]
    
    def ait_plots(self):
        healpix_map.multi_ait(self.hpfit, cmap=plt.get_cmap('jet'),  grid_color='grey')

    def __call__(self, x, n):
        if not hasattr(x, '__iter__'):
            x = np.array([x])
        fit= self.fit[:,n]; 
        t =fit.reshape(3,1)
        return ( t * np.vstack([x**2, x, np.ones(len(x))] )).sum(axis=0)
    
    def ang2pix(self, glon, glat):
        return healpy.ang2pix(64, glon, glat, lonlat=True)
        
    def get_fit(self, pix):
             
        y = self.values[:,pix]
        yerr = 1/self.wtmean
        fn = lambda xx : self(xx, pix)
        return y, yerr, fn
    
    def plot_fit(self, glon, glat, ax=None, axis_labels=True):
        pix = self.ang2pix(glon, glat)
        y, yerr, fn = self.get_fit(pix)

        fig, ax =plt.subplots() if ax is None else (ax.figure, ax)
        npl = len(self.planes)
        xx = np.linspace(self.planes[0]-0.5,self.planes[-1]+0.5,2*npl+1)

        ax.errorbar(self.planes, y, yerr=yerr, fmt='o', ms=8);
        ax.plot(xx, fn(xx), '-', lw=2);
        ax.text(0.05,0.9,'({:3.0f}, {:+2.0f})'.format(glon, glat), transform=ax.transAxes)
        if axis_labels:
            ax.set(ylabel='flux factor', xlabel='energy plane')
        ax.axhline(1.0, color='lightgrey')
        ax.grid(alpha=0.3)
        ax.set_xticks(self.planes[::2])

    def multiplot(self, glons, glats, grid_shape=(4,5), title=''):
 
        fig, axx = plt.subplots(grid_shape[0],grid_shape[1], figsize=(12,12), sharex=True, sharey=True,
                            gridspec_kw=dict(left=0.05, right = 0.95,top=0.95, wspace=0, hspace=0)  )
        for glon, glat, ax in zip(glons, glats, axx.flatten()):
            self.plot_fit( glon, glat, ax, axis_labels=False)
        fig.suptitle(title); 


class DiffuseFitsAnalysis(object):
    """Process diffuse analysis
    """
    def __init__(self, roi, roi_index, nside=64):
        """
        """
        self.roi = roi
        roi.setup_roi(roi_index) # Process object external for now
        self.ri = roi_index
        self.pdirs = map(Band(nside).dir, maps.make_index_table(12,nside)[roi_index])
        print ('Processing ROI index {}'.format(roi_index))
        # load diffuse fits for this model
        files = sorted(glob.glob('diffuse_fit/*.pickle'))
        assert len(files)>0, 'no files found'
        if len(files)<1728:
            msg= "found {} files, expected 1728".format(len(files))
            print (msg)
            raise Exception(msg)
        # return as an array 1728x8x2, last being (gal,iso)
        self.fa = (pd.read_pickle(files[roi_index]).values[:,:2]).astype(float)
        print ('Loaded diffuse fits for this RoI')
        
    def select_band(self, index):
        self.roi.select(index, event_type=None)
        energies = self.roi.energies
        assert len(energies)==1
        self.energy=energies[0]
        print ('Selected {:.0f} MeV'.format(self.energy))
        
        # get the counts in each pixel for gal,iso, and for [front] or [front,back] 
        self.dflux = np.array([[map(resp, self.pdirs)\
                                for resp in sb[:2]] for sb in self.roi.selected])
        
        m = self.dflux.mean(axis=2); s = self.dflux.std(axis=2)/m
        #for a,b in zip(m,s): print ('{:8.0f}, {:.3f}'.format(a,b))
        
    
    def process_band(self, iband, verbose=False):

        self.select_band(iband)

        dflux= self.dflux.astype(float)
        if verbose:
            m = dflux.mean(axis=2); s = dflux.std(axis=2)/m
            ids = 'Front Back'.split() if dflux.shape[0]==2 else ['Front']
            print ('mean\n', pd.DataFrame(m, index=ids, columns='gal iso'.split()))
            print ('rms/mean\n',pd.DataFrame(s, index=ids, columns='gal iso'.split()))

        x=self.fa[iband,:] 
        if verbose: print ('Fits\n', x)

        a64=Band(64).pixelArea()
        f=(dflux*a64); 

        # iso/gal ratio
        r = f[:,1,:]/f[:,0,:]; 
        if verbose: print ('ratio\n', r)

        # adjusted fits per (flux, pixel)
        af = x[0] + r*(x[1]-1); 
        if verbose: print ('adjusted fits\n',af)

        # isolate the galactic flux
        f_gal = f[:,0]; 
        if verbose: print ('galactic fit\n',f_gal)

        # weighted sum of the adjusted fits by the galactic flux
        waf = (f_gal*af).sum(axis=0) / f_gal.sum(axis=0);  
        if verbose: print ('weighted fit\n',self.waf)
        return waf
        
    def process_bands(self, ibands=range(8)):
        return np.array([self.process_band(iband) for iband in ibands])

                