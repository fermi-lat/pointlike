"""
Make various diagnostic plots to include with a skymodel folder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/diagnostic_plots.py,v 1.39 2012/12/23 13:33:10 burnett Exp $

"""

import os, pickle, glob, zipfile, time, sys, types, argparse
import numpy as np
import pylab as plt
import pandas as pd
from uw.like2 import catrec  #deprecated: remove when no longer use sources.rec
from skymaps import SkyDir, DiffuseFunction, Hep3Vector

class Diagnostics(object):
    """ basic class to handle data for diagnostics, collect code to make plots
    """
    def __init__(self, skymodel_dir='.', **kwargs):
        """ skymodel_dir: string
            points to a directory containing a config.txt file, and perhaps other files
            
        """
        self.skymodel_dir = os.path.expandvars(skymodel_dir)
        assert os.path.exists(os.path.join(self.skymodel_dir, 'config.txt')), 'not a skymodel directory:%s'%skymodel_dir
        os.chdir(self.skymodel_dir)
        self.skymodel = os.path.split(os.getcwd())[-1]
        self.setup(**kwargs)
        if not os.path.exists('plots'):         os.mkdir('plots')
        self.plotfolder = os.path.join('plots', self.plotfolder)
        if not os.path.exists(self.plotfolder): os.makedirs(self.plotfolder)

        
    def setup(self, **kwargs):
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
    def get_figure(self, ax=None, figsize=(5,4), **kwargs):
        if ax is not None:
            return ax.figure, ax
        return plt.subplots( figsize=figsize, **kwargs)
        
    def savefigure(self, name, title=None, caption=None, **kwargs):
        """ save a figure.
        name : string
            If name is the name of a function in the class, optionally define 
                the title as the first line, the caption the following lines
        """
        if hasattr(self, name):
            doclines = eval('self.%s' % name).__doc__.split('\n')
            doclines.append('')
            if caption is None:   caption = '\n'.join(doclines[1:])
            if title is None:     title = doclines[0]
        fig= plt.gcf()
        fig.text(0.02, 0.02, self.skymodel, fontsize=8)
        savefig_kw=dict(dpi=60, bbox_inches='tight', pad_inches=0.5) 
        savefig_kw.update(kwargs)
        localfile = '%s_%s.png'%(name, self.skymodel.replace('/','_'))
        savefile = os.path.join(self.plotfolder,localfile)
        
        plt.savefig(savefile, **savefig_kw)
        print 'saved plot to %s' % savefile
        # ceate a simple HTML file is there is a caption
        if caption is not None:
            if title is None: title = name.replace('_', ' ')
            html = '<h2>%s</h2> <img src="%s" /> <br> %s '% (title, localfile, caption)
            open(savefile.replace('.png','.html'),'w').write(html )
        return fig

    def runfigures(self, functions , **kwargs):
        """ functions: list of bound functions 
        """
        for function in functions:
            function(**kwargs)
            self.savefigure(function.__name__)
            
    def load_pickles(self,folder='pickle'):
        """
            load a set of pickles, return list from either zipfile or folder
            (need offset=1 if folder name zipped as well)
        """
        pkls = []
        if os.path.exists(folder+'.zip'):
            print 'unpacking file %s.zip ...' % folder ,
            z=zipfile.ZipFile(folder+'.zip')
            files = sorted(z.namelist()) # skip  folder?
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
        #fig.text(0.025,0.025, 'Asymmetry study %s' % time.asctime(),fontsize='small')
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.95)
        return ax.flatten()
    
    def multilabels(self, xtext, ytext, title=None):
        plt.subplots_adjust(bottom=0.2)
        plt.figtext(0.5,0.07, xtext, ha='center');
        plt.figtext(0.075, 0.5, ytext, rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
        
    def skyplot(self, values, ax=None, title='', vmin=None, vmax=None, 
                    labels=True, colorbar=True, cbtext=''):
        fig, ax = self.get_figure(ax)
        scat =ax.scatter(self.rois.glon, self.rois.singlat, s=20, 
            c=values,  vmin=vmin, vmax=vmax,edgecolor='none')
        if title is not None:
            ax.set_title(title, fontsize='small')
        ax.axhline(0, color='k');ax.axvline(0,color='k')
        if labels: 
            ax.set_xlabel('glon')
            ax.set_ylabel('sin(glat)', labelpad=-5) #note move label to right
        plt.setp(ax, xlim=(180,-180), ylim=(-1.02, 1.02),)
        ax.set_xticks([180,90,0,-90,-180])
        if colorbar:
            cb=plt.colorbar(scat)
            cb.set_label(cbtext)
        return fig, scat
     
    def ecliptic_angle(self, skydir):
        return np.degrees( SkyDir(270,90-23.439281).difference(skydir) ) -90.
        
    def draw_ecliptic(self, ax):
        """ draw the ecliptic path onto the axis ax
        """
        ecl_glon=[]
        ecl_singlat=[]
        zaxis = SkyDir(270,90-23.439281)
        xaxis = SkyDir()
        yaxis = zaxis.cross(xaxis)
        for phi in np.arange(0,2*np.pi, 0.05):
            t = np.sin(phi)*xaxis+np.cos(phi)*yaxis
            sd =SkyDir(Hep3Vector(*t))
            tglon=sd.l(); 
            if tglon>180:tglon-=360
            ecl_glon.append(tglon)
            ecl_singlat.append(np.sin(np.radians(sd.b())))
        ia = np.argsort(ecl_glon)
        ax.plot(np.array(ecl_glon)[ia], np.array(ecl_singlat)[ia], '-', color='gray') 


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
        self.roinames=roinames = [p['name'] for p in pkls]
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
        try:
            self.add_model_info()
        except:
            pass
    def add_model_info(self):
        for i,key in enumerate(['ring','isotrop', 'SunMoon', 'limb',]): # the expected order
            t = []
            for j,p in enumerate(self.pkls):
                if key in p['diffuse_names']:
                    y=p['counts']['models'][i]
                    assert y[0]==key, 'wrong key, roi %d: %s!=%s; list is %s'% (j,key, y[0], p['diffuse_names'])
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
        """ histograms of residuals """
        fig,axx = plt.subplots(3,4, figsize=(12,12))
        for ib,ax in enumerate(axx.flatten()):
            resid = self.residual(ib)
            ax.hist(resid.clip(-5,5), np.linspace(-5,5,21))
            ax.set_title('%.0f MeV'% self.energy[ib], fontsize=10)
            ax.axvline(0, color='k')
            plt.setp(ax, xlim=(-5,5))
            ax.grid(True)
    
    def residual_plot(self):
        """ plot of the average residual
        """
        res = [self.residual(ib) for ib in range(len(self.energy))]
        means = [x.mean() for x in res]
        stds  = [x.std() for x in res]
        fig,ax= plt.subplots(figsize=(4,4))
        ax.errorbar(self.energy, y=means, yerr=stds, fmt='o')
        plt.setp(ax, xscale='log', xlabel='energy', ylabel='average residual')
        ax.set_title('count residuals')
        ax.grid()
        ax.axhline(0, color='k')
        
    def chisq_plots(self,  vmin=0, vmax=100, bcut=10):
        """ chi squared plots
        chi squared distribution
        """
        fig, axs = plt.subplots( 1,2, figsize=(8,3))
        plt.subplots_adjust(wspace=0.3)
        ax = axs[1]
        self.skyplot(self.rois.chisq, ax=ax, vmin=vmin, vmax=vmax);
        ax = axs[0]
        bins = np.linspace(0,100, 26)
        ax.hist(self.rois.chisq.clip(0,100), bins, label='all')
        ax.hist(self.rois.chisq.clip(0,100)[np.abs(self.rois.glat)<bcut], bins, color='red', label='|b|<%d'%bcut)
        ax.legend(loc='upper right', prop=dict(size=10)) 
        plt.setp(ax, xlabel='chisq')
        ax.grid(True)
        return fig
        
    def residual_maps(self, vmin=-5, vmax=5):
        """ Maps of the residuals 
        """
        fig, axx = plt.subplots(3,4, figsize=(12,10))
        plt.subplots_adjust(right=0.9)
        for ib,energy in enumerate(self.energy[:12]):
            ax = axx.flatten()[ib]
            fig, scat=self.skyplot(self.residual(ib).clip(vmin,vmax), ax=ax, title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, colorbar=False, labels=False)
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('normalized residual')
        return fig
        
    def resid_vs_dec(self, ib=0, ax=None, ylim=(-8,8), labels=True):
        """ residual vs. decllination angle
        """
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        r = self.residual(ib).clip(*ylim)
        ax.plot(self.rois.dec, r, '.', color='gray')
        galplane = np.abs(self.rois.glat)<5
        ax.plot(self.rois.dec[galplane], r[galplane], 'or', label='|b|<5')
        ax.axhline(0, color='k')
        plt.setp(ax, ylim=ylim, xlim=(-90,90))
        if labels: plt.setp(ax, xlabel='Dec', ylabel='normalized residual')
        ax.set_xticks((-90,-45,0,45,90))
        ax.set_title('%d MeV' % self.energy[ib], size=10)
        ax.legend(prop=dict(size=10))
        ax.grid()
        return fig
        
    def all_plots(self):
        self.runfigures([self.residual_hists, 
            self.residual_plot, 
            self.residual_maps, 
            self.chisq_plots])
        

class FrontBackSedPlots(Diagnostics):
    """
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
        
    def asym_plot(self, ib, axin=None, rcut=2,  size=15, **kwargs):
        """ ib: band
        """
        fig, ax = self.get_figure( axin)
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
        return fig

    def asym_plots(self):
        map(self.asym_plot,  range(8), self.multifig()); 
        self.multilabels('flux (eV/cm**2/s)','front/back asymmery','Asymmetries for all sources');
        self.savefigure('fb_asymmetry_test');
        return plt.gcf()
        
    def consistency_plot(self, ib, axin=None, vmin=-1, vmax=np.log10(60)):
        fix, ax = self.get_figure( axin)
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
        
    def ratio_fit(self, ib=0, axin=None):
        
        def checkflux( ind, ebin=0):
            fc = np.array([self.flux[x]['flux'][ebin][ind] for x in ['front','back']])
            fu = np.array([self.flux[x]['uflux'][ebin][ind] for x in ['front','back']])
            sig = fu-fc
            ratio = fc[0]/fc[1]; 
            rerr =ratio*np.sqrt((sig[0]/fc[0])**2+(sig[1]/fc[1])**2)
            return  ratio, rerr 
 
        fig, ax = self.get_figure( axin)
        inds = self.get_strongest()
        
        name = [self.srcnames[ind] for ind in inds]
        try:
            realname = [{'P72Y3678':'3C454.3','P7R43539':'3C454.3', 'PSR_J0835-4510':'Vela', 
                            'PSR_J0534p2200':'Crab', 'PSR_J0633p1746':'Geminga'}[n] for n in name]
        except:
            print 'did not find new names: perhaps changed: looking for %s' %name
            realname=name
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
      
    def ratio_plots(self):
        vals = map(self.ratio_fit, range(8), self.multifig())
        plt.suptitle('Front/back flux ratios for strong sources')
        self.savefigure('flux_ratio_strong')
        
        fig, ax = self.get_figure( None)
        
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
        
    def ts_hist(self, ib=0,  space=np.logspace(1,3,21), **kwargs):
        fig,ax=self.get_figure(None)
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
 

class ROIinfo(Diagnostics):
    """ setup pickled DataFrame for ROI info.
    roi name is index
    columns as the individual ROI, except exclude name itself, and enter only list of source names for sources
    """
    def setup(self, **kwargs):
        self.plotfolder='rois'
        self.title='ROI summary'
      
        filename = 'rois.pickle'
        refresh = kwargs.pop('refresh', not os.path.exists('rois.pickle') 
                    or os.path.getmtime('rois.pickle')<os.path.getmtime('pickle.zip') )
        if (not os.path.exists(filename)) or (refresh):
            files, pkls = self.load_pickles('pickle')
            assert len(files)==1728, 'Expected to find 1728 files'
            rdict= dict()
            exclude = ('sources', 'name')
            for pkl in pkls:
                tdict =dict((key,item) for key,item in pkl.items() if key not in exclude )
                tdict.update(sources = pkl['sources'].keys())
                glon = pkl['skydir'].l() 
                if glon>180: glon-=360.
                glat = pkl['skydir'].b(); 
                tdict.update(glon = glon, glat=glat )
                rdict[pkl['name']] = tdict
            self.df = pd.DataFrame(rdict).transpose()
            self.df.save(filename)
            print 'saved %s' % filename
        else:
            print 'loading %s' % filename
            self.df = pd.load(filename)
        self.energy=self.df.ix[0]['counts']['energies']
        
    def skyplot(self, values, ax=None, title='', s=20, vmin=None, vmax=None, ecliptic=False,
                    labels=True, colorbar=True, cbtext=''):
        if ax is None:
            fig, ax = plt.subplots( 1,1, figsize=(5,4))
        else: fig = ax.figure
        singlat=np.sin(np.radians(list(self.df.glat)))
        scat =ax.scatter(self.df.glon, singlat, s=s, 
            c=values,  vmin=vmin, vmax=vmax,edgecolor='none')
        if title is not None:
            ax.set_title(title, fontsize='small')
        ax.axhline(0, color='k');ax.axvline(0,color='k')
        if labels: plt.setp(ax, xlabel='glon', ylabel='sin(glat)',)
        plt.setp(ax, xlim=(180,-180), ylim=(-1.02, 1.02),)
        ax.set_xticks([180,90,0,-90,-180])
        ax.get_yaxis().labelpad=-5 # clugy way to reduce space
        if ecliptic:
            self.draw_ecliptic(ax)
        
        if colorbar:
            cb=plt.colorbar(scat)
            cb.set_label(cbtext)
        return fig, scat # so can add colorbar later

    def model_counts(self, name, ib=0):
        """ list of counts per ROI for model name, energy bin ib
        """
        def cts(i):
            m = self.df.ix[i]['counts']['models']
            dn = self.df.ix[i]['diffuse_names']
            k = dn.index(name) if name in dn else -1
            return m[k][1][ib] if k>=0 else 0
        return np.array(map(cts, range(len(self.df))))

    def diffuse_models(self,  name):
        """ return list of referernces to the diffuse spectral models
        """
        def mdl(index):
            pkl = self.df.ix[index]
            m =pkl['diffuse_names']
            if name not in m: return None
            return pkl['diffuse'][m.index(name)]
        return map(mdl, range(1728))

    def counts_map(self, ib=0, title='', **kwargs):
        
        c = self.model_counts(self.source_name, ib)
        cbtext='counts'
        if kwargs.pop('log', True):
            c = np.log10(c)
            cbtext = 'log10(counts)'
        return self.skyplot(c, cbtext=cbtext,
            title='%s counts at %d MeV' % ( title, self.energy[ib]), **kwargs)
        
    def count_fraction(self, ib=0, title='', **kwargs):
        sm = self.model_counts(self.source_name, ib)
        tot = np.array([self.df.ix[i]['counts']['observed'][ib] for i in range(1728)])
        return self.skyplot(100*(sm/tot), title='%s count fraction at %d MeV' % (title, self.energy[ib]),
            cbtext='fraction (%)', **kwargs)

    def norm_map(self, vmin=0.8, vmax=1.2):
        """ Map of normalization factor for diffuse
        The normalization should be nominally 1.0.
        """ 
        models = self.diffuse_models(self.source_name)
        norms = [m.getp(0) for m in models]
        fig,ax = plt.subplots(figsize=(5,5))
        self.skyplot(norms, ax=ax, vmin=vmin, vmax=vmax, title='Normalization for %s'%self.title)
        return fig
        
    def all_plots(self, **kwargs):
        self.runfigures([self.norm_map], **kwargs)
        
        self.counts_map(title=self.title);
        self.savefigure('%s_counts'%self.source_name, title=self.title)
        self.count_fraction(title=self.title)
        self.savefigure('%s_count_fraction'%self.source_name, title=self.title)


class SunMoon(ROIinfo):
    def setup(self, **kwargs):
        super(SunMoon, self).setup(**kwargs)
        self.plotfolder='sunmoon'
        self.source_name='SunMoon'
        self.title='Sun/Moon'
    def all_plots(self, **kwargs):
        super(SunMoon, self).all_plots( ecliptic=True)


class Limb(ROIinfo):
    def setup(self, **kwargs):
        super(Limb, self).setup(**kwargs)
        self.plotfolder='limb'
        self.source_name='limb'
        self.title='Limb'
     
    def polar_plots(self, values, title=None,
                vmin=None, vmax=None, vticks=5, vlabel=None, thetamax=60):
        """
        values : array of float
        Creates a Figure that must be cleared
        """
        fig, axes = plt.subplots(1,2, figsize=(8,4),subplot_kw=dict(polar=True))
        plt.subplots_adjust(bottom=0.12, top=0.90, wspace=0.3)
        plot_kw = dict(s=120, edgecolor='none', vmin=vmin, vmax=vmax)
        
        galeq = [SkyDir(float(u),0, SkyDir.GALACTIC) for u in np.arange(0,360,1)]
        galra = np.array([sd.ra() for sd in galeq])
        galdec = np.array([sd.dec() for sd in galeq])
        def radius( dec, i): return [90-dec, 90+dec][i] 
        ra = np.array(map(lambda dir: dir.ra(), self.df.skydir))
        dec = np.array(map(lambda dir: dir.dec(), self.df.skydir))
        for i, ax in enumerate(axes[:2]):
            cut = [dec>(90-thetamax), dec<(thetamax-90)][i]
            r = [90-dec, 90+dec][i]
            theta = np.radians(ra)
            c = np.array(values)[cut]
            sc =ax.scatter(theta[cut], radius(dec[cut],i), c=c, **plot_kw);
            galtheta = galra; galrad = radius(galdec,i) 
            ax.plot(np.radians(galtheta[galrad<thetamax]),galrad[galrad<thetamax], '-', color='grey', lw=2)
            ax.set_ylim(ymax=thetamax)
            ax.set_title(['North','South'][i], ha='right', fontsize='small')

        cbax = fig.add_axes((0.25,0.08,0.5, 0.04))
        cb=plt.colorbar(sc, cbax, orientation='horizontal')
        if vlabel is not None: cb.set_label(vlabel)
        if vmin is not None and vmax is not None:
            cb.set_ticks(np.linspace(vmin, vmax, vticks))
        if title is not None: plt.suptitle(title)  
        return fig
        
    def ratio_hist(self, ratio, vmax=1.0, cut=None):
        """ Histogram of Front/Back limb flux ratio
        Histogram of the ratio of the front to the back flux. ROIs with |b|>20 degrees are shown in green. 
        """
        fig, ax = plt.subplots(figsize=(4,4))
        space = np.linspace(0, vmax,21)
        r = ratio.clip(0,vmax)
        ax.hist(r[cut], space)
        ax.hist(r[(abs(self.df.glat)>20) * cut], space, label='|b|>20')
        plt.setp(ax, xlabel='front/back ratio', xlim=(0,vmax))
        ax.set_xticks(np.linspace(0,vmax,5))
        ax.grid(); ax.legend(prop=dict(size=10))
        return fig
        
    def flux_vs_dec(self):
        """ front and back flux vs dec
        Plots of front and back flux normalizations, ploting ROIS with |b|>35
        """
        class PieceWise(object):
            """ functiod that is a piecewise set of straight lines"""
            def __init__(self, a,b):
                self.a, self.b =a,b
                self.n = len(a)
                self.s = [(b[i+1]-b[i])/(a[i+1]-a[i]) for i in range(self.n-1)]
            def __call__(self, x):
                if x<=self.a[0]: return self.b[0]
                for i in range(self.n-1):
                    if x<self.a[i+1]:
                        return self.b[i]+(x-self.a[i])*self.s[i]
                return self.b[-1]
        limbfun = dict(front =PieceWise([-1., -0.4, 0.4, 1.0],[0.75, 0, 0, 0.75]),
                back= PieceWise([-1., -0.7, -0.5, 0.5, 0.7, 0.85, 1.0],
                                [2.0,  0.5, 0,    0,   0.5,  1.2, 0.9])    )
        ra = np.array(map(lambda dir: dir.ra(), self.df.skydir))
        dec = np.array(map(lambda dir: dir.dec(), self.df.skydir))
        dom = np.linspace(-1,1,201) 
        dm = self.diffuse_models('limb')
        fpar,bpar = [np.array([m[i] if m else np.nan for m in dm] )for i in range(2)]
        fig, axx = plt.subplots(2,1, figsize=(8,6))
        plt.subplots_adjust(right=0.9)
        c=np.abs(self.df.glat)
        cut = c>35
        for ax, par, label  in zip(axx, [fpar,bpar], 'front back'.split()):
            scat=ax.scatter(np.sin(np.radians(dec))[cut], par[cut], c=c[cut], s=30,vmin=0, vmax=90)
            plt.setp(ax, xlim=(-1,1), xlabel='sin(dec)' if label=='back' else '',  ylim=(0,2))
            ax.plot(dom, map(limbfun[label],dom), '-', color='k', lw=2) 
            ax.grid()
            ax.text(-0.75, 1.6, label, fontsize=18)
        fig.text(0.05, 0.5, 'flux normalization factor', rotation='vertical', va='center')
        cax = fig.add_axes((0.94, 0.25, 0.02, 0.4))
        cb=plt.colorbar(scat, cax)
        cb.set_label('abs(glat)')
        plt.suptitle('Limb observed flux') 
    
    def all_plots(self, thetamax=60):
        super(Limb, self).all_plots()
        
        dm = self.diffuse_models('limb')
        fpar,bpar = [np.array([m[i] if m else np.nan for m in dm] )for i in range(2)]

        
        self.polar_plots(bpar, vmin=0, vmax=2, vticks=5, thetamax=thetamax, title='back normalization')
        self.savefigure('back_normalization_polar', title='Back normalization factor', caption="""\
        Polar plots of the back limb normalization\
        """)
      
        self.polar_plots(fpar, vmin=0, vmax=1.0, vticks=5, thetamax=thetamax, title='front normalization')
        self.savefigure('front_normalization_polar', title='Front normalization', caption="""\
        Polar plots of the front limb normalization\
        """)
        
        self.polar_plots(fpar/bpar, vmin=0, vmax=1.0, vticks=5, thetamax=thetamax, title='front/back ratio')
        self.savefigure('front_back_ratio_polar', title='Front/Back flux ratio', caption="""\
        Polar plots of the ratio of the front to back flux\
        """)
        
        #self.ratio_hist((fpar/bpar), cut=bpar>0.2)
        self.savefigure('ratio_hist')
        self.flux_vs_dec()
        self.savefigure('flux_vs_dec')


class Galactic(ROIinfo):
    def setup(self):
        super(Galactic, self).setup()
        self.plotfolder='gal'
        self.source_name='ring'
        self.title='Galactic'


class Isotropic(ROIinfo):

    def setup(self):
        super(Isotropic, self).setup()
        self.plotfolder='iso'
        self.source_name='isotrop'
        self.title='Isotropic'

    

class SourceInfo(Diagnostics):
    """ To be superclass for specific source plot stuff, creates or loads
        a DataFrame with all sources 
        (will do away with the recarray, which takes very long to generate, and is rather flaky)
        """
    def setup(self, **kwargs):
        self.plotfolder='sources' #needed by superclass
        filename = 'sources.pickle'
        refresh = kwargs.pop('refresh', not os.path.exists(filename) or os.path.getmtime(filename)<os.path.getmtime('pickle.zip'))
        if refresh:
            files, pkls = self.load_pickles('pickle')
            assert len(files)==1728, 'Expected to find 1728 files'
            sdict= dict()
            for pkl in pkls:
                for name, info in pkl['sources'].items():
                    model = info['model']
                    pars = np.empty(4); pars.fill(np.nan)
                    errs = np.empty(4); errs.fill(-2)
                    free = np.zeros(4, bool)
                    n = model.len()
                    pars[:n] = model.parameters
                    free[:n] = model.free
                    try:
                        errs[:n] = np.diag(model.get_cov_matrix())**0.5
                        errs[np.isnan(errs)]=-1
                        badfit = np.any(errs[model.free]<=0)
                    except Exception, msg:
                        print 'fail errors for %s:%s' % (name, msg)
                        badfit = True
                    sdict[name] = info
                    sdict[name].update(glat=info['skydir'].b(), glon=info['skydir'].l(),
                        roiname=pkl['name'], 
                        pars= pars, errs=errs, free=free, badfit=badfit,
                        e0 = model.e0,
                        modelname=model.name,
                        )
            self.df = pd.DataFrame(sdict).transpose()
            self.df.save(filename)
            print 'saved %s' % filename
        else:
            print 'loading %s' % filename
            self.df = pd.load(filename)
        self.energy = np.sqrt( self.df.ix[0]['sedrec'].elow * self.df.ix[0]['sedrec'].ehigh )
            
    def skyplot(self, values, proj=None, ax=None, s=20, vmin=None, vmax=None, ecliptic=False,
                labels=True, title='', colorbar=True, cbtext=''):
        """ 
        Make a sky plot of some quantity for a selected set of sources
        Parameters:
        ----------
            values: a DataFrame column, posibly a subset: 
                expect to have source name index to get position
            proj: None, or a function to map values to colors
            s : float
                size of dot to plot
        """
        assert hasattr(values, 'index'), 'skyplot: values arg must have index attribute'
        
        # generate arrays of glon and singlat using index 
        sd = self.df.ix[values.index, ['glat', 'glon']] # see page 101
        glon = sd.glon
        glon[glon>180]-=360
        singlat = np.sin(np.radians(list(sd.glat)))

        c = values if proj is None else map(proj, values)

        if ax is None:
            fig, ax = plt.subplots(figsize = (6,5))
        else: fig = ax.figure
        scat = ax.scatter(glon, singlat, s=s, c=c, 
                vmin=vmin, vmax=vmax, edgecolor='none')
        if title:
            ax.set_title(title, fontsize='small')
        
        plt.setp(ax, xlim=(180,-180),  ylim=(-1.02, 1.02));
        ax.axhline(0, color='k');ax.axvline(0,color='k');
        if labels: 
            ax.set_xlabel('glon')
            ax.set_ylabel('sin(glat)', labelpad=-5) #note move label to right

        plt.setp(ax, xlim=(180,-180), ylim=(-1.02, 1.02),)
        ax.set_xticks([180,90,0,-90,-180])
        if ecliptic:
            self.draw_ecliptic(ax)
        if colorbar:
            cb=plt.colorbar(scat)
            cb.set_label(cbtext)
        return fig
        
    def fluxinfo(self, ib=0, cut=None):
        """ extract flux info for energy bin ib, return as a DataFrame
        """
        if cut is None: cut=self.df.ts>25
        s = self.df[cut]
        energy = self.energy[ib]
        fdata = np.array([s.ix[i]['sedrec'].flux[0] for i in range(len(s))])
        udata = np.array([s.ix[i]['sedrec'].uflux[0] for i in range(len(s))])
        ldata = np.array([s.ix[i]['sedrec'].lflux[0] for i in range(len(s))])
        fmodel = np.array([s.ix[i]['model'](energy)*energy**2*1e6 for i in range(len(s))])
        return pd.DataFrame(dict(fdata=fdata, udata=udata, ldata=ldata, fmodel=fmodel, glat=s.glat, glon=s.glon),
            index=s.index)

    def cumulative_ts(self):
        df = self.df
        fig,ax = plt.subplots( figsize=(5,4))
        dom = np.logspace(1,5,1601)
        ax.axvline(25, color='gray', lw=1) 
        ax.hist( df.ts ,dom, cumulative=-1, lw=2, color='g', histtype='step')
        localized = ~np.array(pd.isnull(df.ellipse))
        extended = np.array(df.isextended, dtype=bool)
        unloc = ~ (localized | extended)
        ul = df[unloc * df.ts>10].sort_index(by='roiname')
        # this sorted, can print out if needed

        n = len(ul)
        if n>10:
            ax.hist(ul.ts ,dom, cumulative=-1, lw=2, color='r', histtype='step',
                label='no localization')
            ax.text(12,n, 'failed localization:%d'%n, fontsize=10, color='r')
        plt.setp(ax,  ylabel='sources with greater TS', xlabel='TS',
            xscale='log', yscale='log', xlim=(10, 1e4), ylim=(10,8000))
            
        # label the plot with number at given TS
        for t in (10,25):
            n = sum(df.ts>t) 
            ax.plot([t,2*t], [n,n], '-k');
            ax.text(2*t, n, '%d'%n, fontsize=10, va='center')
                
        ax.set_title('cumulative TS distribution', fontsize='medium')
        ax.grid()
        return fig

    def spectral_fit_consistency(self, ib=0, ax=None, minflux=2.,title=None, bcut=10, hist=False):
        if ax is None:
            fig, ax = plt.subplots(figsize=(5,4))
        else: fig = ax.figure
        
        fx = self.fluxinfo(ib)
        fxc = fx[fx.fmodel>minflux]
        
        x = fxc.fmodel
        y = fxc.fdata
        yerr = (fxc.udata-y)
        q = (y-x)/yerr
        if title is None:
            title = 'spectral fit consistency at %d MeV' % self.energy[ib]
        if not hist:
            self.skyplot(q, ax=ax, vmin=-5, vmax=5, title=title, cbtext='discrepancy in sigmas')
            return fig
        
        # make a hist here
        xlim =(-5,5)
        qlim = q.clip(*xlim)
        lowlat = (abs(fxc.glat)<bcut)
        dom = np.linspace(-5,5,46)
        hist_kw=dict(lw=2, histtype='stepfilled')

        ax.hist(qlim, dom, color='g',  label='%d all sources'%len(q),  **hist_kw)
        ax.hist(qlim[lowlat], dom, color='r',  label='%d |b|<%d' %(sum(lowlat), bcut),  **hist_kw)
        ax.set_xlabel('residual')
        ax.axvline(0, color='k')
        ax.set_xlim(xlim)
        ax.legend(prop=dict(size=10))
        ax.set_title( title, fontsize='medium')
        ax.grid()  
        return fig

    def lowenergyfluxratio(self, ax=None, cut=None, xmax=100., title='low energy fit consistency', energy=133, hist=False, minflux=2.0):
        if cut is None: cut=self.df.ts>25
        if ax is None:
            fig, ax = plt.subplots(figsize=(8,5))
        else: fig = ax.figure
        
        s = self.df[cut]
        fdata = np.array([s.ix[i]['sedrec'].flux[0] for i in range(len(s))])
        udata = np.array([s.ix[i]['sedrec'].uflux[0] for i in range(len(s))])
        ldata = np.array([s.ix[i]['sedrec'].lflux[0] for i in range(len(s))])
        fmodel = np.array([s.ix[i]['model'](energy)*energy**2*1e6 for i in range(len(s))])
        glat = np.array([x.b() for x in s.skydir])
        
        fluxcut = fmodel>minflux
        latcut  = abs(glat)>5.0
        hilat = fluxcut*(latcut)
        lolat = fluxcut*(~latcut)

        y = fdata/fmodel
        ylower, yupper =[(fdata-ldata)/fmodel,(udata-fdata)/fmodel]
        xhi,yhi,yerrhi = fmodel[hilat], y[hilat], [ylower[hilat],yupper[hilat]]
        xlo,ylo,yerrlo = fmodel[lolat], y[lolat], [ylower[lolat],yupper[lolat]]
        
        if not hist:
            ax.errorbar(x=xhi, y=yhi, yerr=yerrhi, fmt='og', label='%d hilat sources'%sum(hilat))
            ax.errorbar(x=xlo, y=ylo, yerr=yerrlo, fmt='or', label='%d lowlat sources'%sum(lolat))
            plt.setp(ax, xlabel=r'$\mathsf{model\ flux\ (eV\ cm^{-2} s^{-1}})$', xscale='log', 
                ylabel='data/model', ylim=(0,2.5), xlim=(minflux, 100) )
            ax.set_xticks([2,5,10,20,50,100])
            ax.set_title( title, fontsize='medium')

        else:
            dom = np.linspace(-4,4,51)
            hist_kw=dict(lw=2, histtype='step')
            q = ((y-1)/yupper).clip(-4,4)
            q[y==0]=-4
            
            ax.hist(q[hilat], dom, color='g',  label='%d hilat sources'%sum(hilat),  **hist_kw)
            ax.hist(q[lolat], dom, color='r',  label='%d lowlat sources'%sum(lolat), **hist_kw)
            ax.set_xlabel('residual')
            ax.axvline(0, color='k')
            ax.set_xlim((-4,4))
        ax.set_title( title, fontsize='medium')
        ax.legend(prop=dict(size=10))
        ax.grid()  

        return fig
    
    def ecliptic_hist(self, ax=None, title=''):
        ea = map(self.ecliptic_angle, self.df.skydir)
        fig, ax = self.get_figure(ax)
        ax.hist( ea, np.linspace(-90, 90, 91))
        plt.setp(ax, xlim=(-90,90), xlabel='ecliptic latitude')
        ax.set_xticks([-90, -45, 0, 45, 90])
        if title: ax.set_title(title)
        ax.grid(True)
        return fig
        
    def fit_quality(self):
        fig, ax = plt.subplots(figsize=(4,4))
        s = self.df
        s['beta'] = s.pars
        logparabola = s.modelname=='LogParabola'
        beta = np.array([pars[2] for pars in s.pars])
        cut=s.band_ts>20
        cut*=(logparabola)
        cut*=(beta<0.01)
        fitqual = s.band_ts-s.ts
        dom = np.linspace(0,50,26)
        ax.hist(fitqual.clip(0,50), dom, label='all non-PSF')
        ax.hist(fitqual[cut].clip(0,50), dom, label=' powerlaw')
        ax.grid(); ax.set_ylim(ymax=600); ax.set_xlabel('fit quality')
        ax.legend(prop=dict(size=10))
        return fig
        
    def pivot_vs_e0(self, xylim=(100, 4e4)):
        fig, ax = plt.subplots(figsize=(4,4))
        s = self.df
        cut = s.ts>10
        ax.plot(s.e0[cut].clip(*xylim), s.pivot_energy[cut].clip(*xylim), '.')
        plt.setp(ax, xscale='log',xlabel='e0', xlim=xylim, 
                    ylabel='pivot', yscale='log', ylim=xylim)
        ax.set_title('compare calculated pivot with e0', fontsize=10)
        ax.grid()
        return fig
    
    
    def all_plots(self):
        self.lowenergyfluxratio(hist=True)
        self.savefigure('low_energy_flux_ratio_hist')
        self.lowenergyfluxratio(hist=False)
        self.savefigure('low_energy_flux_ratio_scat')
        self.spectral_fit_consistency(hist=True)
        self.savefigure('spectral_fit_consistency_hist')
        self.spectral_fit_consistency()
        self.savefigure('spectral_fit_consistency_map')
        self.fit_quality()
        self.savefigure('fit_quality_powerlaw_check')
        self.pivot_vs_e0()
        self.savefigure('pivot_vs_e0')
        self.cumulative_ts()
        self.savefigure('cumulative_ts')
        plt.close('all')
  
class FluxCorr(SourceInfo):

    def setup(self, **kwargs):
        super(FluxCorr, self).setup(**kwargs)
        self.plotfolder='fluxcorr'
        self.source_name='fluxcorr'
        self.title='Source-diffuse flux dependence'
        
        # read in the flux correlation data, add new columns to the source DataFrame
        fs, ps = self.load_pickles('fluxcorr')
        ndf = None
        for x in ps:
            if x is not None:
                ndf = ndf.append(x) if ndf is not None else x
                
        self.ndf = ndf
        
    def flux_sensitivity(self):
        """ Galactic diffuse flux sensitivity
        
        Upper left: histogram of the <br>
        Upper right: scatter plot <br>
        Lower: Skymap for TS<100.
        """
        fig, axx = plt.subplots(2,2, figsize=(8,8))
        ax = axx[0,0]
        flux_ratio = self.ndf['average']
        bins = np.linspace(-30,5,36)
        ax.hist(flux_ratio.clip(bins[0],bins[-1]), bins, label='%d sources'%flux_ratio.count())
        ax.hist(flux_ratio[self.df.ts<100].clip(bins[0],bins[-1]), bins, label='TS<100')
        ax.hist(flux_ratio[self.df.ts<25].clip(bins[0],bins[-1]), bins, label='TS<25')
        plt.setp(ax, xlabel='flux sensitivity')
        ax.grid()
        ax.legend(loc='upper left', prop=dict(size=10))
        ax = axx[0,1]
        ax.plot( self.df['ts'], flux_ratio, '.')
        plt.setp(ax, xscale='log', xlim=(10,1e5), ylim=(-30,5) , xlabel='TS')
        ax.grid()
        ax = axx[1,1]
        tscut = (self.df.ts<100) & (self.df.ts>10)
        axx[1,0].set_axis_off()
        self.skyplot(-flux_ratio[tscut], ax=ax ,s=15, vmax=20, vmin=0,cbtext='abs(flux sensitivity)')
        return fig
        
    def all_plots(self):
        self.runfigures([self.flux_sensitivity])

        
class SourceFitPlots(Diagnostics):
    def setup(self):
        assert os.path.exists('pickle.zip') or os.path.exists('pickle'), 'No pickled ROI data found'
        self.plotfolder='sources'
        recfile = 'sources.rec'
        if not os.path.exists(recfile) or os.path.getmtime(recfile)> os.path.getmtime('pickle.zip'):
            print 'creating %s...' % recfile, ; sys.stdout.flush()
            catrec.create_catalog('.', save_local=True)
        
        sin = pickle.load(open(recfile))
        localized = -np.isnan(sin.a)
        self.use_localization =  sum(localized)>0
        if self.use_localization:
            cut = (sin.ts>10)*(localized)*(sin.pindex<3.5)+ sin.extended
            unloc = (sin.ts>10)*(~localized)*(sin.pindex<3.5)*(~sin.extended)
            self.unloc = pd.DataFrame( sin[unloc], index=pd.Index(sin[unloc].name, name='name'))
            if sum(unloc)>0:
                print 'Ignoring %d sources (%d with TS>25) due to localization failure'\
                %(sum(unloc), sum(unloc*(sin.ts>25)))
        else:
            cut = (sin.ts>10)*(sin.pindex<3.5)+ sin.extended
        print np.sum(sin.ts>10), sum(-np.isnan(sin.a)),np.sum(sin.pindex<3.5) , np.sum(sin.extended)
        self.s = sin[cut]
        print 'found %d  sources, selecting %d for analysis %s' %\
            ( len(cut), sum(cut), ('' if self.use_localization else '(ignoring localization)'))
        self.srcinfo = catrec.FitSource('.')
        self.df = pd.DataFrame(self.s, index=pd.Index(self.s.name, name='name'))
        
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
        self.savefigure('fitquality', title='Fit Quality', caption= 'Left: fit quality histogram; right: fit quality vs. TS'); 
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
        
    def localization(self, maxdelta=9, mints=10):
        bins=np.linspace(0,np.sqrt(maxdelta),26)
        fig, axx = plt.subplots(1,2,figsize=(8,4)); 
        plt.subplots_adjust(wspace=0.4)
        wp = self.s
        cut = wp.ts>mints
        ax=axx[0]
        ax.hist(np.sqrt(wp.delta_ts[cut].clip(0,maxdelta)), bins)
        ax.hist(np.sqrt(wp.delta_ts[wp.ts>100].clip(0,maxdelta)), bins,label='TS>100')
        ax.legend(prop=dict(size=10))
        ax.grid()
        plt.setp(ax, xlabel='sqrt(delta TS)')
        ax=axx[1]
        ax.plot( wp.ts[cut],np.sqrt(wp.delta_ts[cut].clip(0,maxdelta)), '.')
        ax.grid()
        plt.setp(ax, xscale='log', xlabel='TS', ylabel='sqrt(delta TS)')
        return fig
    def localization_quality(self, maxqual=10, mints=10):
        bins=np.linspace(0,maxqual,26)
        fig, axx = plt.subplots(1,2,figsize=(8,4)); 
        plt.subplots_adjust(wspace=0.4)
        wp = self.s
        cut = wp.ts>mints
        ax=axx[0]
        ax.hist(wp.qual[cut].clip(0,maxqual), bins)
        ax.hist(wp.qual[wp.ts>100].clip(0,maxqual), bins,label='TS>100')
        ax.legend(prop=dict(size=10))
        ax.grid()
        plt.setp(ax, xlabel='fit quality')
        ax=axx[1]
        ax.plot( wp.ts[cut],wp.qual[cut].clip(0,maxqual), '.')
        ax.grid()
        plt.setp(ax, xscale='log', xlim=(10,1e5), xlabel='TS', ylabel='fit quality')
        return fig   
    
    def all_plots(self):
        self.fitquality()
        self.lowfluxplots()
        self.isotropic_plot()

        #self.cumulative_ts()
        #self.savefigure('cumulative_ts', title='Cumulative, or logN-logTS plot', caption="""\
        #The cumulative histogram of the number of sources vs. TS""")

        if self.use_localization:
            self.localization()
            self.savefigure('localization', title='Localization plots', caption="""\
            Left: histogram of the square root of the TS difference from current position to
            the fit; corresponds the number of sigmas. Right: scatter plot of this vs. TS
            """)
            self.localization_quality()
            self.savefigure('localization_quality', title='Localization quality plots', caption="""\
            Left: histogram of the fit quality. This is a measure of the difference between the sampled
            TS map points and the prediction of the quadratic model. Right: scatter plot of the quality vs. TS.
            """)

  
class GalDiffusePlots(Diagnostics):

    def diffuse_setup(self, which='gal'):
    
        self.which = which
        self.plotfolder = which #'front_back_%s_diffuse' %which
        folder = '%sfits_all'%which
        if not os.path.exists(folder) and not os.path.exists(folder+'.zip'): folder = folder[:-4]
        files, pkls = self.load_pickles(folder)
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
        
    def like_scat(self, ib, axin=None, fb='both', vmin=0, vmax=2):
        fig, ax=self.get_figure(axin); 
        scat=ax.scatter(self.rois.glon, self.rois.singlat, 
                c=np.log10(self.flux[fb]['deltalike'].transpose().ix[ib]), 
                s=15 if axin is not None else 25,
                vmin=vmin, vmax=vmax, edgecolors='none')
        ax.set_title('%.0f MeV'%(self.energy[ib]),fontsize='small')
        plt.setp(ax, xlim=(180,-180), ylim=(-1.01, 1.01))
        ax.set_xticks([180,90,0,-90,-180])
        return scat
        
    def like_scats(self, title=None):
        fig,ax = plt.subplots(2,4, figsize=(14,8));
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.90, bottom=0.15)
        scats =map(self.like_scat, range(8), ax.flatten());
        plt.figtext(0.5,0.07, 'glon', ha='center');
        plt.figtext(0.05, 0.5, 'sin(glat)', rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scats[0], cbax, orientation='vertical')
        cb.set_label('log10(log likelihood difference)')
        return fig
    
    def bfratio_hist(self, ib, axin=None,  space = np.linspace(0.5, 1.5,26)):
        fig, ax = self.get_figure( axin)
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
        self.bfratios = map(self.bfratio_hist, range(8), ax )
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
        
    #def fit_map(self, axin=None, ib=0, fignum=2, vmin=0.75, vmax=1.25, **kwars):
    #    vals = self.flux['both']['values'].transpose().ix[ib] 
    #    ax = self.set_plot(axin, fignum)
    #    ax.scatter(self.rois.glon, self.rois.singlat, 
    #            c=vals,
    #            s=15 if axin is not None else 25,
    #            vmin=vmin, vmax=vmax, edgecolors='none')
    #    ax.set_title('%.0f MeV'%(self.energy[ib]),fontsize='small')
    #    plt.setp(ax, xlim=(180,-180), ylim=(-1.01, 1.01))
    #    ax.set_xticks([180,90,0,-90,-180])
    #    return ax.figure
    
        
    def diffuse_fits(self):
        ax = self.multifig()
        self.multilabels('ratio', 'ROIs', '%s diffuse fit' %self.which)
        map(self.diffuse_fit, ax, range(8))
        self.savefigure('%sdiffuse_fits'%self.which)
            
    def all_plots(self):
        self.like_scats()
        self.savefigure('%s_delta_log_likelihood_maps'%self.which)
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
       

opts = [('iso',    IsoDiffusePlots),
        ('gal',    GalDiffusePlots),
        ('sources',SourceFitPlots),
        ('fb',     FrontBackSedPlots),
        ('counts', CountPlots),
        ('limb',   Limb),
        ('sourceinfo',   SourceInfo),
        ('sunmoon',  SunMoon),
        ]  
        
def main(args):
    if type(args)==types.StringType: args = [args]
    np.seterr(invalid='warn', divide='warn')
    keys,classes =  [[t[j] for t in opts] for j in (0,1)]
    for arg in args:
        i = keys.index(arg)
        if i<0: print 'found %s; expect one of %s' %(arg,keys)
        else:
            try:
                classes[i]('.').all_plots()
            except FloatingPointError, msg:
                print 'Floating point error running %s: "%s"' % (keys[i], msg)
                print 'seterr:', np.seterr()
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='run a diagnostic output job; must be in skymodel folder')
    parser.add_argument('args', nargs='+', help='processsor identifier: must be one of %s' %[opt[0] for opt in opts])
    #parser.add_argument('-j','--joblist',  help='Optional list of jobs; assume local to $POINTLIKE_DIR', default='job_list')
    #parser.add_argument('--test', action='store_true', help='Do not run the pipeline createStream')
    args = parser.parse_args()
    main(args.args)
    
