"""
Make various diagnostic plots to include with a skymodel folder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/diagnostic_plots.py,v 1.2 2012/11/15 15:33:36 burnett Exp $

"""

import os, pickle, glob, zipfile, time, sys
import numpy as np
import pylab as plt
import pandas as pd
from uw.like2 import catrec

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
    
         # get the basic pickles with the model
        files, pkls = self.load_pickles()
        self.rois =  pkls
        assert len(self.rois)==1728
        self.glon = np.array([r['skydir'].l() for r in self.rois]); self.glon[self.glon>180]-=360
        self.glat = np.array([r['skydir'].b() for r in self.rois])
        self.singlat = np.sin(np.radians(self.glat))
        self.plot_folder = 'chisq'
    
    def set_plot(self, ax, fignum, figsize=(4,4)):
        if ax is None:
            plt.figure(fignum, figsize=figsize);
            ax = plt.gca()
        else:
            plt.sca(ax); 
        return ax
        
    def savefigure(self, name=None, caption=None, **kwargs):
        if name is not None:
            savefig_kw=dict(dpi=60, bbox_inches='tight', pad_inches=0.5); savefig_kw.update(kwargs)
            plt.savefig(os.path.join(self.plotfolder,'%s_%s.png'%(name,self.skymodel.replace('/','_'))), **savefig_kw)
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
        
    def chisq_plots(self,  vmin=0, vmax=100):
        chisq = np.array([r['counts']['chisq'] for r in self.rois])
        fig, axs = plt.subplots( 1,2, figsize=(6,3))
        ax = axs[0]
        scat =ax.scatter(self.glon, self.singlat, s=15, c=chisq,  vmin=vmin, vmax=vmax,edgecolor='none')
        ax.set_title('chisq', fontsize='small')
        ax.axhline(0, color='k');ax.axvline(0,color='k')
        plt.setp(ax, xlabel='glon', ylabel='sin(glat)',xlim=(180,-180), ylim=(-1.02, 1.02),)
        ax = axs[1]
        bins = np.linspace(0,100, 26)
        ax.hist(chisq.clip(0,100), bins, label='all')
        ax.hist(chisq.clip(0,100)[np.abs(self.glat)<5], bins, color='red', label='|b|<5')
        ax.legend(loc='upper right', prop=dict(size=10)) 
        plt.setp(ax, xlabel='chisq')
        fig.savefig('plots/chisq.png')
        return fig


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
    def multifig(self):
        fig,ax = plt.subplots(2,4, figsize=(14,8));
        fig.text(0.025,0.025, 'Asymmetry study %s' % time.asctime(),fontsize='small')
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.95)
        return ax.flatten()
    def multilabels(self, xtext, ytext, title=None):
        plt.figtext(0.5,0.07, xtext, ha='center');
        plt.figtext(0.05, 0.5, ytext, rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)

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
        

class SourceFits(Diagnostics):
    def setup(self):
        assert os.path.exists('pickle.zip') or os.path.exists('pickle'), 'No pickled ROI data found'
        recfile = 'sources.rec'
        if not os.path.exists(recfile):
            print 'creating %s...' % recfile, ; sys.stdout.flush()
            catrec.create_catalog('.', save_local=True)
        sin = pickle.load(open(recfile))
        cut = (sin.ts>10)*(-np.isnan(sin.a))*(sin.pindex<3.5)+ sin.extended
        #cut = (sin.ts>10)*(sin.pindex<3.5)+ sin.extended
        print np.sum(sin.ts>10), sum(-np.isnan(sin.a)),np.sum(sin.pindex<3.5) , np.sum(sin.extended)
        self.s = sin[cut]
        print 'found %d sources, selecting %d for analysis (ignoring localization)' % (len(cut), sum(cut))
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
        
    def all_plots(self):
        self.lowfluxplots()
        self.isotropic_plot()
  
class GalDiffusePlots(Diagnostics):
    def diffuse_setup(self, which='gal'):
        self.plotfolder = 'front_back_%s_diffuse' %which
        folder = '%sfits_all'%which
        if not os.path.exists(folder) and not os.path.exists(folder+'.zip'): folder = folder[:-4]
        files, pkls = self.load_pickles(folder, 0)
        assert len(files)==1728, 'Expect to fine 1728 files in %s' %folder
        self.project_title='diffuse systematics for %s'%self.skymodel
        print self.project_title
        makearray = lambda name,eclass='both' :\
            np.array([p[eclass][name] if eclass is not None else p[name] for p in pkls])
        roinames = makearray('roiname',None)
        
        # DataFrame of info per ROI
        glat   = makearray('glat',None); singlat = np.sin(np.radians(glat))
        glon   = makearray('glon',None); glon[glon>180] -= 360

        self.rois = pd.DataFrame(dict(glat=glat, glon=glon, singlat=singlat), 
            index=roinames)
        
        # create dictionary of data frames of TS, flux for front, back,  Columns are energies
        self.flux = dict()
        for fkey in ['front','back',]:
            print fkey+':',
            self.flux[fkey]=dict()
            for key in [ 'values', 'errors', ]: 
                print key,
                self.flux[fkey][key]= pd.DataFrame( np.array([[t.item() for t in p[fkey][key]] for p in pkls]),
                    index=roinames)
            self.flux[fkey]['deltalike'] = pd.DataFrame( np.array([ p[fkey]['loglike'] for p in pkls]), index=roinames)
            print 'deltalike'
                                                    

    def setup(self):
        self.diffuse_setup('gal')
    def all_plots(self):
        print 'not implemented yet'
        
class IsoDiffusePlots(GalDiffusePlots):
    def setup(self):
        self.diffuse_setup('iso')
        
# temporary -- run in skymodel folder
if __name__=='__main__':
    dp = Diagnostics('.')
    dp.chisq_plots()
