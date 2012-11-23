"""
Make various diagnostic plots to include with a skymodel folder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/diagnostic_plots.py,v 1.9 2012/11/22 00:29:24 burnett Exp $

"""

import os, pickle, glob, zipfile, time, sys, types
import numpy as np
import pylab as plt
import pandas as pd
from uw.like2 import catrec
from skymaps import SkyDir, DiffuseFunction

class Diagnostics(object):
    """ basic class to handle data for diagnostics, collect code to make plots
    """
    def __init__(self, skymodel_dir='.'):
        """ skymodel_dir: string
            points to a directory containing a config.txt file, and perhaps other files
            
        """
        self.skymodel_dir = os.path.expandvars(skymodel_dir)
        assert os.path.exists(os.path.join(self.skymodel_dir, 'config.txt')), 'not a skymodel directory:%s'%skymodel_dir
        os.chdir(self.skymodel_dir)
        self.skymodel = os.path.split(os.getcwd())[-1]
        self.setup()
        if not os.path.exists('plots'):         os.mkdir('plots')
        self.plotfolder = os.path.join('plots', self.plotfolder)
        if not os.path.exists(self.plotfolder): os.makedirs(self.plotfolder)

        
    def setup(self):
        assert False, 'Base class not implemented'
    def describe(self):
        return 'no description'
    def set_plot(self, ax, fignum, figsize=(4,4)):
        if ax is None:
            plt.figure(fignum, figsize=figsize);
            ax = plt.gca()
        else:
            plt.sca(ax); 
        return ax
        
    def savefigure(self, name, caption=None, **kwargs):
        savefig_kw=dict(dpi=60, bbox_inches='tight', pad_inches=0.5) 
        savefig_kw.update(kwargs)
        savefile = os.path.join(self.plotfolder,'%s_%s.png'%(name,self.skymodel.replace('/','_')))
        plt.savefig(savefile, **savefig_kw)
        print 'saved plot to %s' % savefile
        if caption is not None:
            open(savefile.replace('.png','.txt'),'w').write(caption)
        return plt.gcf()

    def load_pickles(self,folder='pickle', offset=1):
        """
            load a set of pickles, return list from either zipfile or folder
        """
        pkls = []
        if os.path.exists(folder+'.zip'):
            print 'unpacking file %s.zip ...' % folder ,
            z=zipfile.ZipFile(folder+'.zip')
            files = sorted(z.namelist())[offset:] # skip  folder?
            print 'found %d files ' % len(files)
            opener = z.open
        else:
           files = sorted(glob.glob(os.path.join(folder,'*.pickle')))
           opener = open
        assert len(files)>0, 'no files found in %s' % folder 
        pkls = [pickle.load(opener(file)) for file in files]
        return files,pkls
        
    def multifig(self):
        fig,ax = plt.subplots(2,4, figsize=(14,8));
        fig.text(0.025,0.025, 'Asymmetry study %s' % time.asctime(),fontsize='small')
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.95)
        return ax.flatten()
    
    def multilabels(self, xtext, ytext, title=None):
        plt.subplots_adjust(bottom=0.2)
        plt.figtext(0.5,0.07, xtext, ha='center');
        plt.figtext(0.075, 0.5, ytext, rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
        
    def skyplot(self, values, ax=None, title='', vmin=None, vmax=None, 
                    labels=True, colorbar=True, cbtext=''):
        if ax is None:
            fig, ax = plt.subplots( 1,1, figsize=(4,4))
        else: fig = ax.figure
        scat =ax.scatter(self.rois.glon, self.rois.singlat, s=20, 
            c=values,  vmin=vmin, vmax=vmax,edgecolor='none')
        if title is not None:
            ax.set_title(title, fontsize='small')
        ax.axhline(0, color='k');ax.axvline(0,color='k')
        if labels: plt.setp(ax, xlabel='glon', ylabel='sin(glat)',)
        plt.setp(ax, xlim=(180,-180), ylim=(-1.02, 1.02),)
        ax.set_xticks([180,90,0,-90,-180])
        if colorbar:
            cb=plt.colorbar(scat)
            cb.set_label(cbtext)
        return fig


class CountPlots(Diagnostics):
    
    def setup(self):
        self.plotfolder = 'counts'

         # get the basic pickles with the model
        files, pkls = self.load_pickles()
        self.pkls = pkls # for development
        assert len(pkls)==1728, 'expect to find 1728 pickled roi files'
        sdirs = [r['skydir'] for r in pkls]
        glon = np.array([r['skydir'].l() for r in pkls]); 
        glon[glon>180]-=360
        glat = np.array([r['skydir'].b() for r in pkls])
        singlat = np.sin(np.radians(glat))
        roinames = [p['name'] for p in pkls]
        self.rois = pd.DataFrame(
            dict(glon=glon, glat=glat, singlat=singlat, 
                ra= [d.ra() for d in sdirs],
                dec = [d.dec() for d in sdirs],
                chisq=[r['counts']['chisq'] for r in pkls],
                bandts=[r['counts']['bandts'] for r in pkls],
                ),
                index = roinames )
        # dict of dataframes with count info. columns are energies
        self.energy = pkls[0]['counts']['energies'] # extract list from first pickle
        counts = [p['counts'] for p in pkls]
        self.counts=dict()
        for key in ['observed', 'total']:
            self.counts[key]= pd.DataFrame([x[key] for x in counts], index=roinames)
        # also model info
        for i,key in enumerate(['ring','isotrop', 'limb',]):
            t = []
            for j,p in enumerate(pkls):
                if key in p['diffuse_names']:
                    y=p['counts']['models'][i]
                    assert y[0]==key, 'wrong key, %s, %s'% (key, y[0])
                    t.append(y[1])
                else:
                    t.append(np.zeros(len(self.energy)))
            self.counts[key]= pd.DataFrame(t, index=roinames)

    def residual(self, ib):
        """ residual DF array for energy band ib 
        """
        obs   = self.counts['observed'].transpose().ix[ib]
        model = self.counts['total'].transpose().ix[ib]
        resid = (obs-model)/np.sqrt(model)
        return resid
     
    def residual_hists(self):
        fig,axx = plt.subplots(3,4, figsize=(12,12))
        for ib,ax in enumerate(axx.flatten()):
            resid = self.residual(ib)
            ax.hist(resid.clip(-5,5), np.linspace(-5,5,21))
            ax.set_title('%.0f MeV'% self.energy[ib], fontsize=10)
            ax.axvline(0, color='k')
            plt.setp(ax, xlim=(-5,5))
            ax.grid(True)
        self.savefigure('residual_hists')
    
    def residual_plot(self):
        res = [self.residual(ib) for ib in range(len(self.energy))]
        means = [x.mean() for x in res]
        stds  = [x.std() for x in res]
        fig,ax= plt.subplots(figsize=(4,4))
        ax.errorbar(self.energy, y=means, yerr=stds, fmt='o')
        plt.setp(ax, xscale='log', xlabel='energy', ylabel='average residual')
        ax.set_title('count residuals')
        ax.grid()
        ax.axhline(0, color='k')
        self.savefigure('average_residual_plot')
        
    def chisq_plots(self,  vmin=0, vmax=100):
        fig, axs = plt.subplots( 1,2, figsize=(6,3))
        ax = axs[0]
        scat =ax.scatter(self.rois.glon, self.rois.singlat, s=15, 
            c=self.rois.chisq,  vmin=vmin, vmax=vmax,edgecolor='none')
        ax.set_title('chisq', fontsize='small')
        ax.axhline(0, color='k');ax.axvline(0,color='k')
        plt.setp(ax, xlabel='glon', ylabel='sin(glat)',xlim=(180,-180), ylim=(-1.02, 1.02),)
        ax = axs[1]
        bins = np.linspace(0,100, 26)
        ax.hist(self.rois.chisq.clip(0,100), bins, label='all')
        ax.hist(self.rois.chisq.clip(0,100)[np.abs(self.rois.glat)<5], bins, color='red', label='|b|<5')
        ax.legend(loc='upper right', prop=dict(size=10)) 
        plt.setp(ax, xlabel='chisq')
        self.savefigure('chisq.png')
        return fig
        
    def residual_maps(self, vmin=-5, vmax=5):
        fig, axx = plt.subplots(3,4, figsize=(12,10))
        for ib,energy in enumerate(self.energy[:12]):
            ax = axx.flatten()[ib]
            self.skyplot(self.residual(ib).clip(vmin,vmax), ax=ax, title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, colorbar=False, labels=False)
        self.savefigure('residual_maps')
        return fig
        
    def resid_vs_dec(self, ib=0, ax=None, labels=True):
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        r = self.residual(ib)
        ax.plot(self.rois.dec, r, '.')
        galplane = np.abs(self.rois.glat)<5
        ax.plot(self.rois.dec[galplane], r[galplane], '+r', label='|b|<5')
        ax.axhline(0, color='k')
        plt.setp(ax, ylim=(-8,8), xlim=(-90,90))
        if labels: plt.setp(ax, xlabel='Dec', ylabel='normalized residual')
        ax.set_xticks((-90,-45,0,45,90))
        ax.set_title('%d MeV' % self.energy[ib], size=10)
        ax.legend(prop=dict(size=10))
        ax.grid()
        return fig
        
    def resid_vs_dec_multi(self):
        fig,axx = plt.subplots(2,4, figsize=(12,8))
        for ib, ax in enumerate(axx.flatten()):
            self.resid_vs_dec(ib, ax=ax, labels=None)
        self.multilabels('Dec', 'nomrmalized residual')
        self.savefigure('resid_vs_dec')
        return fig
    
    def all_plots(self):
        self.residual_hists()
        self.residual_plot()
        self.residual_maps()
        self.resid_vs_dec_multi()
        self.chisq_plots()


class FrontBackSedPlots(Diagnostics):
    """ in progress 
    """
    def setup(self):
        """
        Unpack the pickles, one per source, into convenient DataFrame objects
        """
        files, pkls = self.load_pickles('sedinfo')
        # get energies from first entry, assume all the same
        self.elow  = pkls[0]['elow']
        self.ehigh = pkls[0]['ehigh']
        self.energy= np.asarray(np.sqrt(self.elow*self.ehigh),int)
        
        # extract source names from file names
        def srcname(fname):
            i,j= fname.find('/'), fname.find('_sedinfo')
            return fname[i+1:j]
        srcnames = map(srcname, files)
        self.srcnames = srcnames

        # create DataFrame with basic source ifno
        makearray = lambda name : np.array([p[name] for p in pkls])
        glon  = makearray('glon'); 
        glon[glon>180] -= 360
        glat = makearray('glat')
        self.sourceinfo = pd.DataFrame( dict(
            ts=makearray('ts'), 
            glat=glat, glon=glon, singlat=np.sin(np.radians(glat)),
            ),
            index=srcnames)

        # create dictionary of data frames of TS, flux for front, back, both. Columns are energies
        self.flux = dict()
        for i,fkey in enumerate(['front','back', 'both']):
            self.flux[fkey]=dict()
            for key in [ 'bts', 'flux', 'uflux', 'lflux',]:
                self.flux[fkey][key]= pd.DataFrame( np.array([p[key][i, :] for p in pkls]),
                    index=srcnames)
        # derived for convenience
        a = self.flux['front']['flux']
        b = self.flux['back']['flux']
        self.asymmetry = (a-b)/(a+b)

        # dictionary of diffuse background: galactic and isotropic densities for front and back
        fgal,bgal = [pd.DataFrame(np.array([pkl['bgdensity'][0][i::2] for pkl in pkls]),\
                        index=srcnames) for i in range(2)]
        fiso,biso = [pd.DataFrame(np.array([pkl['bgdensity'][1][i::2] for pkl in pkls]),\
                        index=srcnames) for i in range(2)]
        self.diffuse = dict(fgal=fgal, bgal=bgal, fiso=fiso, biso = biso)
        
        self.plotfolder = 'front_back'
        
    def asym_plot(self, ib, axin=None, rcut=2,  fignum=31, size=15, **kwargs):
        """ ib: band
        """
        ax = self.set_plot( axin, fignum)
        fgal = self.diffuse['fgal'][ib] #gal[:,2*ib] #front, back diffuse density
        fiso = self.diffuse['fiso'][ib] #iso[:,2*ib]
        cut = fgal/fiso>rcut
        asym = self.asymmetry[ib]
        ref_flux = self.flux['both']['flux'][ib]
        if axin is None: size *= 2.0
        defaults =dict(edgecolors='none', s=size)
        defaults.update(kwargs)
        ax.scatter(ref_flux[cut], asym[cut],  c='r', label='gal/iso>%.1f'%rcut, **defaults); 
        ax.scatter(ref_flux[~cut],asym[~cut], c='g', label='gal/iso<%.1f'%rcut, **defaults); 
        plt.setp(ax, xscale='log', xticks=(10, 100), xticklabels=('10','100'), xlim=(4,1000), ylim=(-1.01, 1.01))
        ax.grid(True); ax.legend(prop=dict(size=10))
        ax.axhline(0, color='gray')
        ax.set_title('%0.f-%.0f MeV' % (self.elow[ib],self.ehigh[ib]), fontsize='small')
        if axin is None:
            plt.setp(ax, xlabel=' flux (eV/cm**2/s)', ylabel='front/back asymmery',
            )

    def asym_plots(self):
        map(self.asym_plot, range(8), self.multifig()); 
        self.multilabels('flux (eV/cm**2/s)','front/back asymmery','Asymmetries for all sources');
        self.savefigure('fb_asymmetry_test');
        return plt.gcf()
        
    def consistency_plot(self, ib, axin=None, fignum=13, vmin=-1, vmax=np.log10(60)):
        ax = self.set_plot( axin, fignum)
        ts_f   = self.flux['front']['bts']
        ts_b   = self.flux['back']['bts']
        ts_all = self.flux['both']['bts']
        signif = ts_f+ts_b-ts_all
        c = signif[ib]
        glon, singlat = self.sourceinfo.glon, self.sourceinfo.singlat
        scat = ax.scatter(glon, singlat, s=15 if axin is not None else 25, 
                      c=np.log10(c),  vmin=vmin, vmax=vmax,edgecolor='none')
        bad = c>60
        if sum(bad)>0:
            ax.scatter(glon[bad], singlat[bad], s=50, marker='s', c='k', 
                  edgecolor='none')
        ax.set_title('f-b check for %0.f-%.0f MeV' % (self.elow[ib],self.ehigh[ib]), fontsize='small')
        plt.setp(ax, xlim=(180,-180),  ylim=(-1.02, 1.02));
        ax.axhline(0, color='k');ax.axvline(0,color='k');
        if axin is None: plt.setp(ax,  xlabel='glon', ylabel='sin(glat)')
        return scat

    def consistency_plots(self):
        map(self.consistency_plot, range(8), self.multifig()); 
        self.multilabels('flux (eV/cm**2/s)','front/back asymmery','Asymmetries for all sources');
        self.savefigure('fb_consistency_test');
        return plt.gcf()
    
    def get_strongest(self):
        fluxes = self.flux['both']['flux'][0]
        cutat = sorted(fluxes)[-4]
        strong=fluxes>=cutat
        inds = np.arange(len(strong))[strong]
        #print 'Check strongest sources'
        #print 'fluxes: ', (fluxes[strong]).round()
        assert len(inds)==4, 'Must find four sources, maybe need to adjust cut on strength'
        return inds
        
    def ratio_fit(self, ib=0, axin=None, fignum=11):
        
        def checkflux( ind, ebin=0):
            fc = np.array([self.flux[x]['flux'][ebin][ind] for x in ['front','back']])
            fu = np.array([self.flux[x]['uflux'][ebin][ind] for x in ['front','back']])
            sig = fu-fc
            ratio = fc[0]/fc[1]; 
            rerr =ratio*np.sqrt((sig[0]/fc[0])**2+(sig[1]/fc[1])**2)
            return  ratio, rerr 
 
        ax = self.set_plot( axin, fignum)
        inds = self.get_strongest()
        
        name = [self.srcnames[ind] for ind in inds]
        realname = [{'P72Y3678':'3C454.3', 'PSR_J0835-4510':'Vela', 
                        'PSR_J0534p2200':'Crab', 'PSR_J0633p1746':'Geminga'}[n] for n in name]
        ratio = np.array([checkflux(ind,ib) for ind in inds])
        wts = 1/ratio[:,1]**2; sigma = 1/np.sqrt(np.sum(wts))
        mean  = np.sum( ratio[:,0]*wts)/np.sum(wts)
        #print '%.0f-%.0f: mean ratio = %.3f +/- %.3f' %(self.elow[ib],self.ehigh[ib],mean,sigma)
        ax.axhline(mean, color='g', lw=2, ls='--')
        ax.axhline(1.0, color='k')
        ax.errorbar( range(4), ratio[:,0],yerr=ratio[:,1],lw=2, fmt='', 
                 marker='o', linestyle='None',ms=10,capsize=5)
        ax.errorbar( 1.5, [mean], yerr=[sigma], elinewidth=4, fmt='', marker='x', ms=10,capsize=6, lw=2);
        plt.setp(ax, xlim=(-0.5, 3.5), ylim=(0.85,1.25));
        ax.yaxis.grid(True, linestyle='-', which='major', color='grey',alpha=0.5)
        ax.set_title('%.0f-%.0f MeV'%(self.elow[ib],self.ehigh[ib]), fontsize='medium')
        xticknames = plt.setp(ax, xticklabels=realname, xticks=range(4))
        if axin is None: ax.set_ylabel('front/back flux ratio')
        return (self.elow[ib],self.ehigh[ib],mean,sigma)
      
    def ratio_plots(self, fignum=12):
        vals = map(self.ratio_fit, range(8), self.multifig())
        plt.suptitle('Front/back flux ratios for strong sources')
        self.savefigure('flux_ratio_strong')
        
        ax = self.set_plot( None, fignum)
        
        y  = [v[2] for v in vals] 
        yerr = np.array([v[3] for v in vals])
        xmin = np.array([v[0] for v in vals])
        xmax = np.array([v[1] for v in vals])
        x = np.sqrt(xmin*xmax)
        xerr= (x-xmin, xmax-x)
        #print len(y),len(yerr)
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, marker='o', ms=12,fmt='', lw=2, linestyle='None')
        plt.setp(ax, xscale='log', ylim=(0.85,1.25), xlabel='Energy (MeV)', ylabel='front/back flux ratio',)
        ax.grid(True)
        ax.axhline(1.0, color='k')
        ax.set_title('Point source spectral fits', fontsize='medium')
        self.savefigure('fb_flux_vs_energy', dpi=60)
        
    def ts_hist(self, ib=0, fignum=201, space=np.logspace(1,3,21), **kwargs):
        plt.close(fignum); fig=plt.figure(fignum, figsize=(4.5,4.5), dpi=100)
        ax = plt.gca()
        defaults = dict( histtype='step', lw=2)
        defaults.update(kwargs)
        ts = self.sourceinfo.ts
        flux = [self.flux[x]['flux'][ib] for x in ['front','back','both']]
        
        ax.hist(ts ,             space,color='b', label='all', **defaults);
        ax.hist(ts[flux[2]==0], space,color='r', label='zero total flux', **defaults);
        ax.hist(ts[flux[1]==0], space,color='g', label='zero back flux',**defaults);
        ax.hist(ts[flux[0]==0], space,color='orange', label='zero front flux',**defaults);
        ax.grid();
        plt.setp(ax, xlabel='TS', xscale='log');
        ax.set_title('TS with zero flux, energy %.0f MeV'%self.energy[ib], fontsize='medium');
        ax.legend(prop=dict(size=10))  
        self.savefigure('ts_with_zero_flux')
 
    def all_plots(self):
        self.asym_plots()
        self.consistency_plots()
        self.ratio_plots()
        self.ts_hist()
        plt.close('all')
        

class SourceFitPlots(Diagnostics):
    def setup(self):
        assert os.path.exists('pickle.zip') or os.path.exists('pickle'), 'No pickled ROI data found'
        recfile = 'sources.rec'
        if not os.path.exists(recfile):
            print 'creating %s...' % recfile, ; sys.stdout.flush()
            catrec.create_catalog('.', save_local=True)
        sin = pickle.load(open(recfile))
        localized = -np.isnan(sin.a)
        self.use_localization =  sum(localized)>0
        if self.use_localization:
            cut = (sin.ts>10)*(localized)*(sin.pindex<3.5)+ sin.extended
        else:
            cut = (sin.ts>10)*(sin.pindex<3.5)+ sin.extended
        print np.sum(sin.ts>10), sum(-np.isnan(sin.a)),np.sum(sin.pindex<3.5) , np.sum(sin.extended)
        self.s = sin[cut]
        print 'found %d  sources, selecting %d for analysis %s' %\
            ( len(cut), sum(cut), ('' if self.use_localization else '(ignoring localization)'))
        self.plotfolder='sources'
        self.srcinfo = catrec.FitSource('.')
        
    def fitquality(self):
        fig, axs = plt.subplots(1,2, figsize=(7,3))
        plt.subplots_adjust(wspace=0.35)
        s = self.s
        fitqual = s.band_ts-s.ts
        from scipy import stats
        ndf=12
        chi2 = lambda x: stats.chi2.pdf(x,ndf)
        d = np.linspace(0,100,51); delta=d[1]-d[0]
        ax =axs[0]
        ax.hist(fitqual, d, log=False);
        ax.hist(fitqual[s.ts>500], d, label='TS>500');
        ax.plot(d, chi2(d)*len(fitqual)*delta/1.6, 'r', label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
        plt.setp(ax, xlabel='fit qual', ylim=(0,500))
        ax.grid(); ax.legend(prop=dict(size=10))
        ax = axs[1]
        ax.plot(s.ts, fitqual, '.'); 
        plt.setp(ax, xscale='log', xlabel='TS', xlim=(10,1e5),
             ylabel='fit qual',ylim=(1,1e3),yscale='log')
        ax.grid()
        self.savefigure('fitquality')#, 'Left: fit quality histogram; right: fit quality vs. TS'); 
        return fig #plt.close(fig)
        
    def _lowfluxplot(self, ax, cut, xmax=100., title='selected bad fits', energy=133):
        """ s is list of sources"""
        from skymaps import SkyDir
        s = self.s
        fitqual = s.band_ts-s.ts
        sd = map(SkyDir, s[cut].ra, s[cut].dec)
        glat = np.array([x.b() for x in sd])
        hilat = abs(glat)>5.0
        fss = [self.srcinfo(src) for src in s[cut]]
        fdata = np.array([fs['sedrec'].flux[0] for fs in fss])
        udata = np.array([fs['sedrec'].uflux[0] for fs in fss])
        ldata = np.array([fs['sedrec'].lflux[0] for fs in fss])
        fmodel = np.array([fs['model'](energy)*energy**2*1e6 for fs in fss])
        y = fdata/fmodel
        yerr=[(fdata-ldata)/fmodel,(udata-fdata)/fmodel]
        xhi,yhi,yerrhi = fmodel[hilat],  y[hilat],  [((fdata-ldata)/fmodel)[hilat],((udata-fdata)/fmodel)[hilat]]
        xlo,ylo,yerrlo = fmodel[~hilat], y[~hilat], [((fdata-ldata)/fmodel)[~hilat],((udata-fdata)/fmodel)[~hilat]]
        ax.errorbar(x=xhi, y=yhi, yerr=yerrhi, fmt='og', label='%d hilat sources'%sum(hilat))
        ax.errorbar(x=xlo, y=ylo, yerr=yerrlo, fmt='or', label='%d lowlat sources'%sum(~hilat))
        plt.setp(ax, xlabel=r'$\mathsf{model\ flux\ (eV\ cm^{-2} s^{-1}})$', xscale='log', 
            ylabel='data/model', ylim=(0,2.5),)
        ax.set_title( title, fontsize='medium')
        ax.set_xlim( (1,100) )

        ax.axhline(1.0,color='k')
        ax.legend(prop=dict(size=10))
        ax.grid()  
  
    def lowfluxplots(self):
        s = self.s
        fitqual = s.band_ts-s.ts
        bad = (s.ts>100)*(fitqual>50); 
        fig, ax = plt.subplots(1,2, figsize=(12,5))
        self._lowfluxplot(ax=ax[0], cut=bad, xmax=100, title='133 MeV flux ratio, selected bad fits')
        self._lowfluxplot(ax=ax[1], cut=s.ts>1000, xmax=2000, title='133 MeV flux ratio for TS>1000')
        self.savefigure('low_flux_plots')
        return fig
        
    def isotropic_plot(self):
        config = eval(open('config.txt').read())
        diffuse=config['diffuse']
        idfiles = [os.path.join(os.environ['FERMI'],'diffuse',diffuse[1][i]) for i in (0,1)]
        nf,nb = map(np.loadtxt, idfiles)
        energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
        fig, axs = plt.subplots(1,2, figsize=(7,3), dpi=50)
        ax= axs[1]
        ax.plot(energies, front/back, '-o');
        ax.axhline(1.0, color='k')
        plt.setp(ax, xscale='log', xlabel='Energy');ax.grid(True);
        ax.set_title('Isotropic flux front/back ratio', fontsize='small');
        ax = axs[0]
        ax.plot(energies, front*energies**2, '-g', label='front')
        ax.plot(energies, back*energies**2, '-r', label='back')
        plt.setp(ax, xlabel='Energy', ylabel='flux*e**2', xscale='log')
        ax.set_title('isotropic diffuse spectra', fontsize='small')
        ax.grid(True); ax.legend()
        self.savefigure('isodiffuse', 'isotropic spectra:%s'% list(diffuse[1]))
        
    def localization(self, bins=np.linspace(0,10,26)):
        fig, axx = plt.subplots(1,2,figsize=(7,3)); 
        plt.subplots_adjust(wspace=0.4)
        wp = self.s
        cut = wp.ts>10
        ax=axx[0]
        ax.hist(wp.delta_ts[cut].clip(0,10), bina)
        ax.hist(wp.delta_ts[wp.ts>100].clip(0,10), bins,label='TS>100')
        ax.legend(prop=dict(size=10))
        plt.setp(ax, xlabel='delta TS')
        ax=axx[1]
        ax.plot( wp.ts[cut],wp.delta_ts[cut].clip(0,10), '.')
        ax.grid()
        setp(ax, xscale='log', xlabel='TS', ylabel='delta TS')
        self.savefigure('localization')
    
    def cumulative_ts(self):
        s = self.s
        fig,ax = plt.subplots( figsize=(5,3))
        dom = np.logspace(1,5,1601)
        ax.axvline(25, color='k', lw=2) 
        ax.hist( s.ts ,dom, cumulative=-1, lw=2, color='g', histtype='step')
        plt.setp(ax, xscale='linear', ylabel='sources with greater TS', xlabel='TS',
             xlim=(10,40), ylim=(2000,4000) )

        cut = (s.ts>25)
        n25 = sum(cut); ax.plot([23,27], [n25,n25], '--r');
        n10 = len(s.ts); ax.plot([10,12], [n10,n10], '--r');
        ax.text(27,n25, '%d'%n25, fontsize='small', va='center')
        ax.text(12,n10, '%d'%n10, fontsize='small', va='center')
        ax.set_title('cumulative TS distribution', fontsize='medium')
        ax.grid()
        self.savefigure('cumulative_ts')
    
    def all_plots(self):
        self.fitquality()
        self.lowfluxplots()
        self.isotropic_plot()
        self.cumulative_ts()
        if self.use_localization:
            self.localization()
  
class GalDiffusePlots(Diagnostics):

    def diffuse_setup(self, which='gal'):
    
        self.which = which
        self.plotfolder = 'front_back_%s_diffuse' %which
        folder = '%sfits_all'%which
        if not os.path.exists(folder) and not os.path.exists(folder+'.zip'): folder = folder[:-4]
        files, pkls = self.load_pickles(folder, 0)
        assert len(files)==1728, 'Expect to find 1728 files in %s' %folder
        self.project_title='diffuse systematics for %s'%self.skymodel
        print self.project_title
        makearray = lambda name,eclass='both' :\
            np.array([p[eclass][name] if eclass is not None else p[name] for p in pkls])
        roinames = makearray('roiname',None)
        self.energy = makearray('energies')[0];
        # DataFrame of info per ROI
        glat   = makearray('glat',None); singlat = np.sin(np.radians(glat))
        glon   = makearray('glon',None); glon[glon>180] -= 360
        self.latcut, self.latcut_name = (abs(glat)<10, 'plane') if which=='gal' else (abs(glat)>20, 'high-lat')
        self.rois = pd.DataFrame(dict(glat=glat, glon=glon, singlat=singlat), 
            index=roinames)
        
        # create dictionary of data frames of flux for front, back, with values, energies, deltalike;
        #   Columns are energies
        self.flux = dict()
        for fkey in ['front','back','both']:
            print '\t',fkey+':',
            self.flux[fkey]=dict()
            for key in [ 'values', 'errors', ]: 
                print key,
                self.flux[fkey][key]= pd.DataFrame( np.array([[t.item() for t in p[fkey][key]] for p in pkls]),
                    index=roinames)
            self.flux[fkey]['deltalike'] = pd.DataFrame( np.array([ p[fkey]['loglike'] for p in pkls]), index=roinames)
            print 'deltalike'
                                                    

    def setup(self):
        self.diffuse_setup('gal')
        
    def like_scat(self, ib, axin=None,fignum=22, fb='both', vmin=0, vmax=2):
        ax=self.set_plot(axin,fignum); 
        ax.scatter(self.rois.glon, self.rois.singlat, 
                c=np.log10(self.flux[fb]['deltalike'].transpose().ix[ib]), 
                s=15 if axin is not None else 25,
                vmin=vmin, vmax=vmax, edgecolors='none')
        ax.set_title('%.0f MeV'%(self.energy[ib]),fontsize='small')
        plt.setp(ax, xlim=(180,-180), ylim=(-1.01, 1.01))
        ax.set_xticks([180,90,0,-90,-180])
        
    def like_scats(self):
        map(self.like_scat, range(8), self.multifig());
        self.multilabels('glon', 'sin(glat)', '%s diffuse significance'%self.which)
        self.savefigure('%s_significance'%self.which)
    
    def bfratio_hist(self, ib, axin=None, fignum=101, space = np.linspace(0.5, 1.5,26)):
        ax = self.set_plot( axin, fignum)
        f,b = [self.flux[fb]['values'].transpose().ix[ib] for fb in ['front', 'back'] ]
        bfratio = f/b
        ax.hist(bfratio, space, label='all', histtype='stepfilled',color='g')
        ax.hist(bfratio[self.latcut], space, histtype='stepfilled',
             label=self.latcut_name, color='r')
        ax.set_title('%.0f MeV' % self.energy[ib], fontsize='medium')
        plt.setp(ax, xlim=(space[0],space[-1]), )
        ax.axvline(1.0, color='k', lw=2)
        if axin is None:
            ax.set_xlabel('%s diffuse front/back fit ratio'%self.which);
            ax.legend(loc='upper left');
        ax.grid(True); 
        return (self.energy[ib],  bfratio[self.latcut].mean(), bfratio[self.latcut].std())
        
    def bfratio_hists(self):
        ax = self.multifig()
        self.bfratios = map(self.bfratio_hist, range(8),ax )
        self.multilabels('front/back fit', '', '%s Diffuse fit ratio'%self.which)
        ax[0].legend(loc='upper left',bbox_to_anchor=(-0.4,1.2));
        self.savefigure('bf_%s_fit.png'%self.which)
        
        html_rows = ['<tr><td>%.0f</td><td>%.2f</td></tr>' %v[:2] for v in self.bfratios]
        from IPython.core.display import HTML
        h=HTML('<table><tr><th>Energy</th><th>front/back ratio</th></tr>'+\
             ''.join(html_rows)+'</table>')
        open(os.path.join(self.plotfolder,'%s_fb_ratio.html'%self.which),'w').write(h.data)
  
    def diffuse_ratio_plot(self, fignum=121):
        ax = self.set_plot( None, fignum)
        vals = self.bfratios # must have generated
        x,y,yerr = [[v[i] for v in vals]  for i in range(3)]
        ax.errorbar(x, y, yerr=yerr, marker='o', ms=12,fmt='', lw=2, linestyle='None')
        plt.setp(ax, xscale='log',xlabel='Energy (MeV)', ylabel='front/back flux ratio',ylim=(0.55, 1.25))
        ax.grid(True)
        ax.axhline(1.0, color='k')
        ax.set_title('%s diffuse spectral fits'%self.which, fontsize='medium')
        self.savefigure('%s_diffuse_flux_ratio_vs_energy'%self.which, dpi=60)
        
    def diffuse_fit(self, axin, ind=0,  fignum=2, **kwargs):
        plane = abs(self.rois.glat)<10
        cut, cut_name = (plane, 'plane') if self.which=='gal' else (abs(self.rois.glat)>20, 'high-lat')
        space=np.linspace(0.5,1.5, 41)
        vals = self.flux['both']['values'].transpose().ix[ind] #values[:,ind]
        kw = dict( histtype='stepfilled')
        kw.update(kwargs)
        ax = self.set_plot(axin, fignum)
        ax.hist(vals, space,label='all', **kw)
        ax.hist(vals[cut], space ,color='r',label=cut_name,**kw);
        ax.grid(True);
        ax.axvline(1.0, color='grey');
        ax.set_xlim((0.5,1.5))
        ax.set_title('%.0f MeV'%self.energy[ind],fontsize='medium')
        ax.legend(prop=dict(size=10))
        
    def diffuse_fits(self):
        ax = self.multifig()
        self.multilabels('ratio', 'ROIs', '%s diffuse fit' %self.which)
        map(self.diffuse_fit, ax, range(8))
        self.savefigure('%sdiffuse_fits'%self.which)
    
    def all_plots(self):
        self.like_scats()
        self.bfratio_hists()
        self.diffuse_ratio_plot()
        self.diffuse_fits()
        
class IsoDiffusePlots(GalDiffusePlots):

    def setup(self):
        self.diffuse_setup('iso')

        
    def lowdiff_hist(self, ax, **kwargs):
        plane=abs(self.rois.glat)<10
        v0,v1 = [self.flux['both']['values'][i] for i in (0,1)]
        delta=v1-v0
        kw = dict(histtype='stepfilled'); kw.update(kwargs)
        space = np.linspace(-0.5, 0.5, 41)
        ax.hist(delta, space, label='all: mean %.2f'%delta.mean(), **kw)
        ax.hist(delta[plane], space, color='r', label='plane: mean %.2f'%delta[plane].mean(), **kw)
        ax.legend(); ax.grid();
        ax.axvline(0, color='grey')
        ax.set_title('%s diffuse fits %s'%(self.which,self.skymodel))
        ax.set_xlabel('bin1 - bin0 difference')


    def lowdiff_scat(self, ax, vmin=-0.2, vmax=0.2, **kwargs):
        v0,v1 = [self.flux['both']['values'][i] for i in (0,1)]
        delta=v1-v0
        kw = dict(edgecolor='none');kw.update(kwargs)
        t=ax.scatter(self.rois.glon, self.rois.singlat, c=delta, s=50,vmin=vmin, vmax=vmax, **kw)
        plt.setp(ax, xlabel='glon', ylabel='sin(glat)', xlim=(180,-180), ylim=(-1,1))
        ax.axhline(0, color='k')
        ax.axvline(0, color='k')
        ax.axhline(np.sin(np.radians(10.)), lw=2, color='grey')   
        ax.axhline(np.sin(np.radians(-10.)), lw=2, color='grey')
        # draw poles outline
        try:
            poles = pickle.load(open('../../polar_circle.pickle'))
        except:
            print 'could not find the polar_circle file'
            return
        for i in range(2):
            ax.plot(poles[0,:,0]-360, np.sin(np.radians(poles[0,:,1])), '-',lw=2, color='grey')
            ax.plot(poles[1,:,0], np.sin(np.radians(poles[1,:,1])), '-', lw=2,color='grey')

    def all_plots(self):
        self.like_scats()
        self.bfratio_hists()
        self.diffuse_ratio_plot()
        self.diffuse_fits()
        
        fig,ax=plt.subplots(1,2, figsize=(14,6))
        self.lowdiff_hist( ax[0])
        self.lowdiff_scat( ax[1])
        self.savefigure('isotropic_bin0_bin1_difference')
       
class LimbPlots(Diagnostics):
    def setup(self):
        self.plotfolder='limb'
        files, pkls = self.load_pickles()
        assert len(pkls)==1728, 'expect to find 1728 pickled roi files'
        
        cut = np.array(['limb' in p['diffuse_names'] for p in pkls])
        print 'Found %d rois with limb fits'% sum(cut)
        
        indexes = np.arange(1728)[cut]; 
        limb_pkls = [pkls[i] for i in indexes]

        limb_dirs = [p['skydir'] for p in limb_pkls]
        ra = np.array([d.ra() for d in limb_dirs])
        dec = np.array([d.dec() for d in limb_dirs])
        
        fitpars = np.array([p['diffuse'][2].get_all_parameters()  for p in limb_pkls])
        self.models = [  p['diffuse'][2]  for p in limb_pkls]
        front_flux=[]; back_flux=[]
        energy=133 #wire in first band for now
        for m in self.models:
            m.ct=0; front_flux.append( m(energy))
            m.ct=1; back_flux.append(  m(energy))      
        
        # get the diffuse function
        try:
            config = eval(open('config.txt').read())
            self.limb_file = os.path.join(os.path.expandvars('$FERMI/diffuse'),config['diffuse'][2])
            print 'loading diffuse definition %s' %self.limb_file
            df = DiffuseFunction(self.limb_file)
            flux = [df(sd, energy) for sd in limb_dirs]
        except Exception, msg:
            print 'Failed to load fluxes from diffuse definition: %s' %msg
            flux = np.zeros(len(indeces))
            
        self.df = pd.DataFrame(dict(ra=ra,dec=dec, 
                    glat = [d.b() for d in limb_dirs],
                    glon = [d.l() for d in limb_dirs],
                    fpar=fitpars[:,0], bpar=fitpars[:,1],
                    flux=flux,
                    ), 
                index=pd.Index(indexes, name='ROI index'))

    def describe(self):
        return self.df.describe()
        
    def polar_plots(self, values, title=None,
                vmin=None, vmax=None, vticks=5, vlabel=None):
        """
        values : array of float
            Must have same length as self.df
        Creates a Figure that must be cleared
        """
        fig, axes = plt.subplots(1,2, figsize=(8,4),subplot_kw=dict(polar=True))
        plt.subplots_adjust(bottom=0.12, top=0.90, wspace=0.3)
        plot_kw = dict(s=120, edgecolor='none', vmin=vmin, vmax=vmax)
        
        #galeq = pickle.load(open('analysis/galactic_equator.pickle'))
        galeq = [SkyDir(float(u),0, SkyDir.GALACTIC) for u in np.arange(0,360,1)]
        galra = np.array([sd.ra() for sd in galeq])
        galdec = np.array([sd.dec() for sd in galeq])
        def radius( dec, i): return [90-dec, 90+dec][i] 
        ra,dec =self.df.ra, self.df.dec
        for i, ax in enumerate(axes[:2]):
            cut = [dec>45, dec<-45][i]
            r = [90-dec, 90+dec][i]
            theta = np.radians(ra)
            sc =ax.scatter(theta[cut], radius(dec[cut],i), c=values[cut], **plot_kw);
            galtheta = galra; galrad = radius(galdec,i) 
            ax.plot(np.radians(galtheta[galrad<60]),galrad[galrad<60], '-', color='grey', lw=2)
            ax.set_ylim(ymax=50)
            ax.set_title(['North','South'][i]+' Pole', ha='right', fontsize='medium')

        cbax = fig.add_axes((0.25,0.08,0.5, 0.04))
        cb=plt.colorbar(sc, cbax, orientation='horizontal')
        if vlabel is not None: cb.set_label(vlabel)
        if vmin is not None and vmax is not None:
            cb.set_ticks(np.linspace(vmin, vmax, vticks))
        if title is not None: plt.suptitle(title)  
        return fig

    def ratio_hist(self):
        fig, ax = plt.subplots(figsize=(4,4))
        ratio = (self.df.fpar/self.df.bpar).clip(0,0.4)
        space = np.linspace(0, 0.4,21)
        ax.hist(ratio, space)
        ax.hist(ratio[abs(self.df.glat)>20], space, label='|b|>20')
        plt.setp(ax, xlabel='front/back ratio', xlim=(0,0.4))
        ax.set_xticks(np.linspace(0,0.4,5))
        ax.grid(); ax.legend(prop=dict(size=10))
        return fig
        
    def all_plots(self):
        self.polar_plots(self.df.bpar, vmin=0, vmax=2, title='Back normalization')
        self.savefigure('back_normalization_polar')
    
        self.polar_plots(self.df.fpar, vmin=0, vmax=0.5, title='Front normalization')
        self.savefigure('front_normalization_polar')
        self.polar_plots(self.df.flux, title='Nominal Flux at 133 MeV')
        self.savefigure('nominal_flux_polar', caption="""Nominal flux at 133 MeV, from diffuse file %s"""%self.limb_file)
        self.ratio_hist()
        self.savefigure('front_back_ratio_hist', caption="""
            Polar plots of the ratio of measured front to back""")
        
def main(args=None):
    opts = [('iso',    IsoDiffusePlots),
            ('gal',    GalDiffusePlots),
            ('sources',SourceFitPlots),
            ('fb',     FrontBackSedPlots),
            ('counts', CountPlots),
            ('limb',   LimbPlots),
            ]  
    keys,classes =  [[t[j] for t in opts] for j in (0,1)]
    if args is None:
        print 'require an argument: expect one, or list of %s' %keys
        return
    if type(args)==types.StringType:
        args = [args]
    for arg in args:
        i = keys.index(arg)
        if i<0: print 'found %s; expect one of %s' %(arg,keys)
        else:
            classes[i]('.').all_plots()
        
if __name__=='__main__':
    if len(sys.argv)<2:  main()
    main(sys.argv[1:])
