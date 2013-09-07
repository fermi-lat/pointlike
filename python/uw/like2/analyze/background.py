"""
background analysis
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/background.py,v 1.8 2013/09/04 12:34:59 burnett Exp $

"""
import os, glob
import pandas as pd
import numpy as np
import pylab as plt
from scipy import integrate, optimize
from skymaps import SkyDir
from . import roi_info
from . analysis_base import FloatFormat, html_table

class Background(roi_info.ROIinfo):
    r"""Background analysis, including corrections and sensitivity
    <br>This is an analysis of the sensitivity of point source analysis to the presence of 
    diffuse background. The ingredients are the PSF and the galactic and isotropic background fluxes.
    <br>
    Detailed results are presented for weak and strong high- and low-latitude sources.
    """
    
    def introduction(self):
        r"""Introduction: Resolution for a point source
        The following is from a <a href="../../../../../documents/background_systematics.htm?skipDecoration">derivation</a> of the relative resolution for a fit to the flux from a point source, 
        with fixed spectral parameters:
        $$\begin{equation}  \frac{1}{{\sigma_s}^2} = \sum\limits_k    
          \int \mathrm{d}\Omega  \frac{P_k(\theta)^2}{ \alpha_k P_k(\theta)+\beta_k }
        \end{equation}$$
        where:
        <ul>
        <li>$\sigma_s$: Statistical uncertainty for the relative source flux; $s=1$ represents a multiplicative factor applied 
            to each $\alpha$ in the derivation 
        <li>$k$: Band index, including front or back
        <li>$\alpha_k$: predicted counts for the point source in band $k$, up to a common scale factor $s$
        <li>$\beta_k$: diffuse background count density at the source location, assumed constant within the PSF extent, 
        and known precisely <em>ab initio</em> or from a fit to the rest of the ROI
        <li>$P_k(\theta)$: the normalized ($\int\mathrm{d}\Omega P_k(\theta)=1$) PSF, as a function of angle $\theta$ from the source position
        <li>$\int\mathrm{d}\Omega$: Integral over solid angle of an ROI, typically a cone about the source position: 
            $\mathrm{d}\Omega=\sin(\theta)\mathrm{d}\theta \mathrm{d}\phi$. The derivation assumes that the PSF
            is contained within the ROI.
        </ul>
        <p>Consider a single band. For the background-dominated case, $\beta >> \alpha P(0)$, $\sigma_s^2=\beta/ \int\mathrm{d}\Omega P(\theta)^2$. 
        The integral is the average value of the PSF, equal to the inverse of a solid angle. This solid angle is referred to as the 
        PSF "footprint" below. The value for $\sigma_s^2$ is then the number of background counts in the footprint.
        If there is no background, then $\sigma_s^2=\alpha$, the number of source events, independent of the PSF.
        """
        return None
    
    def setup(self, **kw):
        super(Background, self).setup(**kw)
        self.plotfolder='background'
        
        s = [[ x[1]['counts']['models'][modelnumber][1][:16] for x in self.df.iterrows()] for modelnumber in range(2)]
        self.bdf = [pd.DataFrame(y, index=self.df.index) for y in s]
        self.roi_size=5 # wired-in for now
        self.rois=[888,0]
        self.example_sources='P7R41835 P7R42139 P7R42771 P7R41444'.split()
        # use to find count background for low, high latitude sources
        self.bgsources = 'P7R41835 P7R42771'.split()
        # get the current PSF
        self.config = eval(open('config.txt').read())
        self.psf = self.get_psf()
        
        # Get the full sedinfo for all sources
        files, pkls = self.load_pickles('sedinfo')
        # extract source names from file names
        def srcname(fname):
            i,j= fname.find('/'), fname.find('_sedinfo')
            return fname[i+1:j]
        srcnames = map(srcname, files)
        self.srcnames = srcnames
        self.sdict = dict(zip(srcnames, pkls))
        self.band_energy = np.zeros((28))
        self.source_type = np.zeros((28))
        e = np.array(np.sqrt(pkls[0]['elow']*pkls[0]['ehigh']),np.float32)
        self.band_energy[0::2] =e
        self.band_energy[1::2] =e
        self.source_type[0::2]=0
        self.source_type[1::2]=1
        self.solid_angle = 2*np.pi*(1-np.cos(np.radians(5)))
        
        #make source info available from the csv
        self.dfs = pd.read_csv(glob.glob('sources_uw*.csv')[-1], index_col=0)
        
        # make sedinfo into a DataFrame sedinfo
        self.sedinfo= pd.DataFrame(self.sdict).T
        self.sedinfo.index = map(lambda n: n.replace('_',' ').replace('p','+'), self.sedinfo.index)
        self.sedinfo.index.name = 'name'


    def get_background(self, roi):
        roiname='HP12_%04d' % roi
        return [t.ix[roiname] for t in self.bdf]
    
    def diffuse_flux(self, rois=None):
        """Diffuse flux
        Predicted counts per ROI, for low and high-latitude, front and back combined.
        """
        fig, ax = plt.subplots(1,1, figsize=(6,6), dpi=150, sharey=True)
        egev = np.array(self.energy)/1e3
        if rois is None: rois = self.rois

        for r in rois:
            gal, iso = self.get_background(r)
            ax.plot(egev, gal, '-d', label='gal %d'%r)
            ax.plot(egev, iso, '--o', label='iso %d'%r)
        plt.setp(ax, xscale='log', xlim=(0.1,300), xlabel='Energy (GeV)',
            yscale='log',  ylabel='Diffuse counts/ROI')
        ax.legend(prop=dict(size=10)); ax.grid()
        return fig

    def get_psf(self, irfname=None, ):
        from uw.like import pypsf, pycaldb
        if irfname is None: irfname=self.config['irf']
        cdm = pycaldb.CALDBManager(irf=irfname)
        self.psf_files=cdm.get_psf()
        return pypsf.CALDBPsf(cdm)

    def psf_background(self, bgsource_names=None, rois=None, outfile='psf_bkg.csv',):
        """Background counts in the PSF
        For galactic and isotropic backgrounds, and front and back, determine the number of counts 
        in the effective "footprint" of the PSF, see equation (1) in the introduction. These values are
        significant for two reasons:
        <ol><li>They define the limiting resolution for the number of source counts
        <li>The set a scale for a scheme to limit the resolution if background systematic uncertainties are an issue
        </ol>
        """
        # get effective PSF size, and area from environment analysis
        epsf = pd.read_csv('plots/environment/psf.csv', index_col=0)
        psa = [np.pi*np.radians(epsf[ct].values[:14])**2 for ct in ['front', 'back']]
        
        # assume 5 degrees for solid angle (should get the value actually used)
        solid_angle = 2*np.pi*(1-np.cos(np.radians(self.roi_size)))
        
        fig, axx = plt.subplots(2,2, figsize=(12,12), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.,top=0.9, left=0.1)
        energy = self.energy 
        
        if bgsource_names is None: bgsource_names = self.bgsources
        flist = []
        blist = []
        
        for k,source_name in enumerate(bgsource_names):
            res = self.resolution(source_name)
            gal, iso = [res[x] for x in ('gal', 'iso')]
            for ax, diffuse, what in zip(axx[k,:], (gal, iso), 'galactic isotropic'.split()):
                
                front, back = [ diffuse[i::2]*psa[i]/solid_angle for i in range(2)]
                flist.append(front.values)
                blist.append(back.values)
                
                ax.plot(energy, front, '-_', lw=2, label='front')
                ax.plot(energy, back, '-+r', lw=2, label='back')
                ax.legend(prop=dict(size=10)); ax.grid()
                ax.text(150,2e5, what+['-high lat','-low lat'][k], fontsize=12)
                plt.setp(ax, xscale='log', xlim=(100, 12000),
                    yscale='log', ylim=(1.0, 1e6), ylabel='counts in PSF' if what=='galactic' else '',)
            
        axx[0,1].set_xticklabels(['0.1', '1', '10'])
        fig.text(0.4, 0.05, 'Energy (GeV)')
        #make a DataFrame with the 8 plots
        self.psf_back_df = pd.DataFrame(dict(
            hgf=flist[0],  hgb=blist[0], 
            hif=flist[1],  hib=blist[1], 
            lgf=flist[2],  lgb=blist[2], 
            lif=flist[3],  lib=blist[3], 
             ), 
        index=(np.array(energy)*1e-3).round(3))
        self.psf_back_df.index.name='energy'
        self.psf_back_df.to_csv(os.path.join(self.plotfolder, outfile))
        print 'wrote file %s' % os.path.join(self.plotfolder, outfile)
    
        return fig
            
    def source_info(self,  source_name):
        return self.dfs.ix[source_name]
        
    def resolution(self, source_name, s=1):
        """ Return a data frame with resolution information
        s is a scale factor, to examine dependence on flux
        
        bratio : beta/alpha, where beta is average diffuse density for ROI, alpha the source counts
        gratio : beta/alpha, beta for galactic only. (So bratio-gratio is isotropic)
        pfun : value of the integral of P**2/(P+r), where P=PSR, r = bratio
        snratio: P(0)/bratio , signal/noise
        ivar : s * pfun, the expected inverse variance per bin
        """
        def Pfun( energy, ct, r):
            def f(x):
                P = self.psf(energy, int(ct), x)
                return x * P**2 /(r + P)
            return 2.*np.pi * integrate.quad(f, 0, np.pi/6) [0]
        try:
            p=self.sdict[source_name.replace(' ','_').replace('+','p')]
            data = np.array([self.band_energy.round(0), self.source_type, 
                p['counts'].round(0), p['bgcounts'][0].round(0), p['bgcounts'][1].round(0)])
        except:
            print '***Failed to find resolution info for source %s' % source_name
            raise
        df = pd.DataFrame(data, index='energy ct source gal iso'.split()).T
        alpha = s * df.source
        df['bratio']=(df.gal+df.iso)/self.solid_angle/alpha
        df['gratio']=(df.gal)/self.solid_angle/alpha
        df['pfun'] = [Pfun(e,ct,br) for e,ct, br in zip(df.energy, df.ct, df.bratio)]
        df['snratio'] = np.array([self.psf(e, int(ct), 0.)[0] for e,ct in zip(df.energy,df.ct)]) / df.bratio 
        df['ivar'] = alpha * df.pfun
        return df
        
    def plot_resolution(self, source_names=None):
        """contribution to the resolution for measuring the source strength
        These plots show the contributions from each band to the sum in equation (1).
        Representative sources are chosen for either high or low latitude and weak or strong.
        Each point is the inverse variance of the relative error, or $1/\sigma^2$.
        Combining front and back, the total relative error is shown in the table below.
        Note that it is evaluated at the nominal fit flux.
        
        %(resolution_summary)s
        """
        if source_names is None: source_names = self.example_sources
        src_info = map(self.source_info, source_names)
        fit_rel_sig = [ 100.*u['flux_unc']/u['flux'] for u in src_info]
        fig,axx = plt.subplots(2,len(source_names)/2, figsize=(10,10))
        plt.subplots_adjust(left=0.1)
        rel_sig=[]
        for ax, source_name in zip(axx.flatten(), source_names):
            try:
                df = self.resolution(source_name)
            except: 
                rel_sig.append(1)
                continue
            total_ivar = 0
            for ct, marker, label in zip((0,1), ('bo', 'rD'), 'front back'.split()):
                iv = df.ivar[ct::2]
                total_ivar += sum(iv)
                ax.plot(df.energy[ct::2]*1e-3,iv  , marker, label='%s, total=%.1f' % (label,sum(iv)))
            plt.setp(ax, xscale='log', xlabel='Energy (GeV)', xlim=(0.1, 500), title=source_name)
            ax.legend(prop=dict(size=10)); ax.grid()
            ax.set_xticklabels('0.1 1 10 100'.split())
            rel_sig.append( 100./np.sqrt(total_ivar))
        sdir = map(SkyDir, [s['ra'] for s in src_info], [s['dec'] for s in src_info])
        self.rtable = pd.DataFrame( [ [s['ts'] for s in src_info], [s.b() for s in sdir], [s['flux']*1e13 for s in src_info],[s['pindex'] for s in src_info],
                fit_rel_sig, rel_sig ],
            columns= source_names, index ='TS glat flux index meas_sig pred_sig'.split()).T 
        self.rtable.index.name = 'name'
        self.resolution_summary = html_table(self.rtable, dict(glat=',galactic latitude', index=',spectral index', 
            flux=',flux at pivot*1e13',
            meas_sig=',measured relative flux error',pred_sig=',predicted relative flux error'),float_format= FloatFormat(2))
        return fig
    
    def gal_corrections(self, bands=8):
        """Galactic diffuse corrections
        Analysis of the correction factors used for this model.<p>  %(correction_info)s
            <p>The diffuse correction factors were obtained by optimizing the IEM component separately
            for every ROI and the 8 bands below 10 GeV. Here we examine how well each of three alternate 
            strategies can account for these corrections:
            <ol><li>No correction
            <li>Adjust total flux only
            <li>Adjust the flux and spectral slope
            </ol> 
            <p>The histograms show the results, the residuals from comparing the corrected values with the prediction of each fit.
            The scale is percent, and the plot legends show the RMS for each case.
        """
        try:
            corr_csv = self.config['diffuse']['ring']['correction']
            self.correction_info='Correction info file: "%s"' % corr_csv
        except Exception, msg:
            print '*** no correction file found %s' % msg
            self.correction_info='No correction file foound'
            return
        t = self.get_background(888)[0].values[:bands]
        weights = t/t.sum()
        galcor = pd.read_csv(corr_csv, index_col=0)
        resid =  []
        pars1 = []
        pars2 = []
        for r, y in galcor.iterrows():
            ydata = y.values[:bands] # just to 1 GeV, even weighting
            func1 = lambda params: (ydata - params[0])*weights
            par1, code = optimize.leastsq(func1, (0,))
            func2 = lambda params: (ydata - params[0] - params[1]*np.arange(bands))*weights
            par2, code = optimize.leastsq(func2, (0,0))
            pars1.append( par1)
            pars2.append( par2 )
            resid.append([ 100.*u for u in (ydata-1, func1(par1), func2(par2))])
        resid = np.array(resid)

        fig, axx = plt.subplots(2,2, figsize=(10,10), sharex=True)
        plt.subplots_adjust(wspace=0.1)
        for j,ax in enumerate(axx.flatten()):
            for i,label in enumerate(('no corr', 'average', 'linear')):
                d = resid[:,i,j]
                ax.hist(d, np.linspace(-5,5,51), lw=2, histtype='step', 
                    label=label+' %.2f'%d.std())
            ax.grid()
            ax.legend(prop = dict(size=10))
            plt.setp(ax, xlabel='residual band %d'%j, xlim=(-5,5))
            
        return fig
 
    def gal_correction_maps(self, bands=4):
        """galactic correction maps
        The following maps show the distribution of the first %(nmaps)d corrections.
        %(correction_info)s
        %(gal_corr_maps)s
        
        """
        from . maps import Maps
        from uw.like2.pub import healpix_map as hm
        self.nmaps = bands
        try:
            corr_csv = self.config['diffuse']['ring']['correction']
            self.correction_info='Correction info file: "%s"' % os.path.join(os.getcwd(),corr_csv)
        except Exception, msg:
            print '*** no correction file found %s' % msg
            self.correction_info='No correction file found'
            self.gal_corr_maps=''
            return
        galcor = pd.read_csv(corr_csv, index_col=0)
        tt = [hm.HParray('band%d'%x,galcor['%d'%x]) for x in range(bands)]
        bb = [hm.HPresample(t) for t in tt]
        hpf = hm.HEALPixFITS(bb)
        fitsfile = os.path.join(self.plotfolder, 'diffuse_corrections.fits')
        hpf.write(fitsfile)
        mps = Maps()
        mps.plotfolder = self.plotfolder
        self.gal_corr_maps = mps.ait_plots(fitsfile, vmin=0.9, vmax=1.1)
   
    def read_covariance(self, iband=0):
        # requires that the covariance zip or folder, resulting from a special 'covariance' stage exists.
        files, pkls = self.load_pickles('covariance')
        print 'selecting band %d' % iband
        dp = dict()
        for pk in pkls:
            name = pk.get('name')
            try:
                cov =pk['covariance'][0][iband] #only lowest energy/conversion (second 0)
                sd = pk['skydir']
                sigs=cov['sigs'], #relative sigmas, order source, gal, iso
                if len(sigs)==1:
                    sigs = sigs[0]
                assert len(sigs)==3, 'problem with sigs: %s,len =%s' % (sigs, len(sigs))
                cc = cov['cc'] # correlation coefficients, order src-gal, src-iso, gal-iso
                assert  len(cc)==3, 'problem with cc: %s, len=%s' % (cc, len(cc) )
            except Exception, msg:
                print '*** fail to unpack source name %s: %s' %(name,msg)
                continue
            dp[name]= dict(ts=pk['ts'], skydir=pk['skydir'], 
                glat=sd.b(), cc=cc, sigs=sigs, unweight=cov['unweight'],
                galrat = cc[0]*sigs[0]/sigs[1], 
                isorat = cc[1]*sigs[0]/sigs[2])
        return pd.DataFrame(dp).T
        
    def covariance(self, iband=0):
        """Background sensitivity
        Evaluate the sensitivity of the measurement flux of point sources with
        respect to the covariance matrix of the source with respect to both.
        The plots are for the lowest energy band.
        
        """
        cv = self.read_covariance(iband=iband)
        
        def plot2(clipto=(-100,50)):
            fig, axx = plt.subplots(1,4, figsize=(14,5))
            aglat=np.abs(np.array(cv.glat.values, float))
            ax=axx[0]
            scat=ax.scatter(cv.galrat.clip(*clipto), cv.isorat.clip(*clipto), marker='.',s=80, c=aglat, )
            #plt.colorbar(scat,ax)
            plt.setp(ax, xlabel='galactic', xlim=clipto, ylabel='isotropic',ylim=clipto)
            ax.grid()
            ax=axx[1]
            ax.hist(cv.galrat.clip(*clipto), np.linspace(*clipto))
            plt.setp(ax, xlabel='galrat')
            ax=axx[2]
            ssig = np.array([s[0] for s in cv.sigs])
            ax.scatter(ssig.clip(0,1), cv.galrat.clip(*clipto), marker='.', s=80, c=aglat)
            plt.setp(ax, xlabel='source sigma', xlim=(0,1),ylim=(-40,0))
            ax=axx[3]
            cut = cv.ts>500
            ssig = np.array([s[0] for s in cv.sigs])
            ax.scatter(ssig[cut].clip(0,1), cv.galrat[cut].clip(*clipto), marker='.', s=80, c=aglat[cut])
            plt.setp(ax, xlabel='source sigma', xlim=(0,1),ylim=(-40,0))
            
            return fig
        return plot2()

    def all_plots(self):
        self.runfigures([self.introduction, self.plot_resolution,  self.psf_background,
            self.diffuse_flux, self.gal_corrections, self.gal_correction_maps,])
