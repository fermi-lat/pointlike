"""
Make various diagnostic plots to include with a skymodel folder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/diagnostic_plots.py,v 1.83 2013/04/09 21:24:09 burnett Exp $

"""

import os, pickle, glob, zipfile, time, sys, types, argparse, pyfits
import numpy as np
import pylab as plt
import pandas as pd
import pyfits
import mpl_toolkits.axes_grid.axes_size as Size
from mpl_toolkits.axes_grid import axes_grid, axes_size, Divider, make_axes_locatable
from matplotlib.colors import LogNorm

from skymaps import SkyDir, DiffuseFunction, Hep3Vector, Band
from uw.like2.pub import healpix_map
from uw.like2 import dataset

class Diagnostics(object):
    """ basic class to handle data for diagnostics, collect code to make plots
    """
    def __init__(self, skymodel_dir='.', **kwargs):
        """ skymodel_dir: string
            points to a directory containing a config.txt file, and perhaps other files
            
        """
        self.skymodel_dir = os.path.expandvars(skymodel_dir)
        assert os.path.exists(os.path.join(self.skymodel_dir, 'config.txt')), 'not a skymodel directory:%s'%skymodel_dir
        if skymodel_dir != '.': os.chdir(self.skymodel_dir)
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

    def get_figure(self, ax, figsize=(5,4), **kwargs):
        if ax is not None:
            return ax.figure, ax
        return plt.subplots( figsize=figsize, **kwargs)
 
    def subplot_array( self, hsize, vsize=(1.0,), figsize=(10,10)):
        """ Use the axes_divider module to make a single row of plots
        hsize : list of floats
            horizontal spacing: alternates Scaled for plot, Fixed for between plots
        vsize : list of floats
            vertical spacing
            
        ref:   http://matplotlib.org/mpl_toolkits/axes_grid/users/axes_divider.html
        """
        nx = (len(hsize)+1)/2
        ny = (len(vsize)+1)/2
        fig, axx = plt.subplots(ny,nx,squeeze=False, figsize=figsize) # just to make the axes, will move them
        sizer = lambda x,i: axes_size.Scaled(x) if i%2==0 else axes_size.Fixed(x)
        horiz = [ sizer(h,i) for i,h in enumerate(hsize) ]
        vert  = [ sizer(v,i) for i,v in enumerate(vsize) ]
        divider = Divider(fig, (0.1, 0.1, 0.8, 0.8), horiz, vert, aspect=False)
        for i,ax in enumerate(axx.flatten()):
            iy = i//nx; ix = i%nx
            ax.set_axes_locator(divider.new_locator(nx=2*ix, ny=2*iy))
        return fig, axx
        
    def savefigure(self, name, func=None, title=None, caption=None, **kwargs):
        """ save a figure.
        name : string
            If name is the name of a function in the class, optionally define 
                the title as the first line, the caption the following lines
        func : executable function, or None
            if not None, run the func, use it to get docs
            If func creates a figure, it must return it
        Note that the docstring may have %(xxx)s, which will be replaced by attribute xxx.
        """
        if func is not None:
            fig=func(**kwargs)
            fname = func.__name__
        else: fname = name
        if hasattr(self, fname):
            try:
                doclines = ((eval('self.%s' % fname).__doc__%self.__dict__).split('\n'))
                doclines.append('')
                if caption is None:   caption = '\n<p>'+'\n'.join(doclines[1:])+'</p>\n'
                if title is None:     title = doclines[0]
            except Exception, msg:
                print 'docstring processing problem: %s' % msg
        localfile = '%s_%s.png'%(name, self.skymodel.replace('/','_'))
        savefile = os.path.join(self.plotfolder,localfile)
        if title is None: title = name.replace('_', ' ')
        if fig is not None:
            fig.text(0.02, 0.02, self.skymodel, fontsize=8)
            savefig_kw=dict(dpi=60, bbox_inches='tight', pad_inches=0.5) 
            plt.savefig(savefile, **savefig_kw)
            print 'saved plot to %s' % savefile
            html = '<h3>%s</h3> <img src="%s" />\n <br> %s '% (title, localfile, caption if caption is not None else '')
        else:
            html = '<h3>%s</h3>\n <br>  %s' % (title, caption if caption is not None else '')
        open(savefile.replace('.png','.html'),'w').write(html )
        print 'saved html doc to %s' % os.path.join(os.getcwd(),savefile.replace('.png','.html'))
        return html

    def runfigures(self, functions, names=None,  **kwargs):
        """ 
        run the function, create web page containing them
        functions: list of bound functions 
        names: names to use instad of function names
        
        Expect to be called from all_plots, get a summary from its docstring
        """
        if names is None:
            names=[None]*len(functions)
        html = '<head>'+ HTMLindex.style + '\n<title>%s</title> </head>\n' % self.__class__.__name__
        html+='<body><h2>%(header)s</h2>'
        docstring = self.all_plots.__doc__
        if docstring is not None: html+=docstring
        for function, name in zip(functions,names):
            fname = name if name is not None else function.__name__
            fig = self.savefigure(fname, function, **kwargs)
            html+='\n'+ fig
        html+= '\n<hr>\nPage generated %4d-%02d-%02d %02d:%02d:%02d' % tuple(time.localtime()[:6])

        html+='\n</body>'
        t = os.path.split(os.getcwd())
        self.header='/'.join([t[-1], os.path.split(self.plotfolder)[-1]])
        try:
            text = html%self.__dict__
        except:
            print 'failed filling %s' % html
            raise
        open(os.path.join(self.plotfolder,'index.html'), 'w').write(text)
        print 'saved html doc to %s' %os.path.join(self.plotfolder,'index.html')
            
    def load_pickles(self,folder='pickle'):
        """
            load a set of pickles, return list from either zipfile or folder
            (need offset=1 if folder name zipped as well)
        """
        pkls = []
        if os.path.exists(folder+'.zip'):
            print 'unpacking file %s.zip ...' % (os.getcwd()+'/'+folder ,),
            z=zipfile.ZipFile(folder+'.zip')
            files = sorted(z.namelist()) # skip  folder?
            print 'found %d files ' % len(files)
            if folder=='pickle' and len(files)==1729:
                files = files[1:]
            opener = z.open
        else:
           files = sorted(glob.glob(os.path.join(folder,'*.pickle')))
           opener = open
        assert len(files)>0, 'no files found in %s' % folder 
        pkls = [pickle.load(opener(file)) for file in files]
        return files,pkls
        
    def multifig(self):
        fig,ax = plt.subplots(2,4, figsize=(14,8));
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.95)
        return ax.flatten()
    
    def multilabels(self, xtext, ytext, title=None):
        plt.subplots_adjust(bottom=0.2)
        plt.figtext(0.5,0.07, xtext, ha='center');
        plt.figtext(0.05, 0.5, ytext, rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
        
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
        
    def basic_skyplot(self, ax, glon, singlat, c,
                title=None, ecliptic=False, labels=True, colorbar=False, cbtext='', 
                aspect=180.,  **scatter_kw):
        """ basic formatting used for ROI and sources
            note that with aspect=180, the aspect ratio is 1:1 in angular space at the equator
        """
        cb_kw = scatter_kw.pop('cb_kw', {}) 
        ecliptic = scatter_kw.pop('ecliptic', False)
        scat = ax.scatter(glon, singlat, c=c, **scatter_kw)
        if title:
            ax.set_title(title, fontsize='small')
        
        plt.setp(ax, xlim=(180,-180),  ylim=(-1.02, 1.02));
        ax.axhline(0, color='k');ax.axvline(0,color='k');
        if labels: 
            ax.set_xlabel('glon')
            ax.set_ylabel('sin(glat)', labelpad=-5) #note move label to right

        plt.setp(ax, xlim=(180,-180), ylim=(-1.02, 1.02),aspect=aspect,)
        ax.set_xticks([180,90,0,-90,-180])
        ax.set_xticklabels([180,90,0,270, 180])
        if ecliptic:
            self.draw_ecliptic(ax)
        if colorbar:
            # supposed to be nice, didn't work with already-locatable?
            #http://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right", size="5%", pad=0.05)
            #cb=plt.colorbar(scat, cax=cax)
            cb=ax.figure.colorbar(scat, ax=ax, **cb_kw)
            cb.set_label(cbtext)    
        return scat       


class CountPlots(Diagnostics):
    require='pickle.zip'
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

    def counts_map(self):
        """ Sum, for E>100 Mev
        """
        obs = self.counts['observed']
        total = np.array([sum(x[1]) for x in obs.iterrows()])
        sy
    def residual(self, ib):
        """ residual DF array for energy band ib 
        """
        obs   = self.counts['observed'].transpose().ix[ib]
        model = self.counts['total'].transpose().ix[ib]
        resid = (obs-model)/np.sqrt(model)
        return resid
     
    def residual_hists(self):
        """ histograms of normalized residuals 
        subset for ridge (|b|<10, |l|<60) shown
        """
        fig,axx = plt.subplots(3,4, figsize=(12,12))
        ridge = ( np.abs(self.rois.glat)<10) * ( np.abs(self.rois.glon)<60 )

        for ib,ax in enumerate(axx.flatten()):
            resid = self.residual(ib)
            ax.hist(resid.clip(-5,5), np.linspace(-5,5,21))
            ax.hist(resid[ridge].clip(-5,5), np.linspace(-5,5,21))
            ax.set_title('%.0f MeV'% self.energy[ib], fontsize=10)
            ax.axvline(0, color='k')
            plt.setp(ax, xlim=(-5,5))
            ax.grid(True)
        return fig
    
    def residual_plot(self):
        """ plot of the average normalized residual
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
        
    def chisq_plots(self, hsize=(1.0, 0.8, 1.5, 0.5), vmin=0, vmax=50, bcut=10):
        """ chi squared plots
        chi squared distribution
        """
        #fig, axs = plt.subplots( 1,2, figsize=(8,3))
        #plt.subplots_adjust(wspace=0.3)
        fig, axs = self.subplot_array( hsize, figsize=(11,5))
        chisq = self.rois.chisq

        def chisky(ax):
            self.basic_skyplot(ax, self.rois.glon, self.rois.singlat, chisq, 
                s=60, vmin=vmin, vmax=vmax,  edgecolor='none', colorbar=True);
                
        def chihist(ax):
            bins = np.linspace(0, vmax, 26)
            lolat = np.abs(self.rois.glat)<bcut
            ax.hist(chisq.clip(0,vmax), bins, label='all: mean=%.1f'%chisq.mean())
            ax.hist(chisq.clip(0,vmax)[lolat], bins, color='red', label='|b|<%d (%.1f)'%(bcut, chisq[lolat].mean()))
            ax.legend(loc='upper right', prop=dict(size=10)) 
            plt.setp(ax, xlabel='chisq', xlim=(0,vmax))
            ax.grid(True)
            
        for f, ax in zip( (chihist, chisky), axs.flatten()): f(ax)
        return fig
        
    def residual_maps(self, vmin=-5, vmax=5):
        """ Maps of the residuals 
        """
        fig, axx = plt.subplots(3,4, figsize=(12,10), sharex=True, sharey=True)
        plt.subplots_adjust(right=0.9, hspace=0.15, wspace=0.1)
        for ib,energy in enumerate(self.energy[:12]):
            ax = axx.flatten()[ib]
            scat=self.basic_skyplot(ax, self.rois.glon, self.rois.singlat, self.residual(ib).clip(vmin,vmax),
                 title='%d MeV'%energy,
                vmin=vmin,vmax=vmax, s=15, edgecolor='none', colorbar=False, labels=False)
        #put colorbar at right        
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scat, cbax, orientation='vertical')
        cb.set_label('normalized residual')
        fig.text(0.5,0.05, 'galactic longitude', ha='center')
        fig.text(0.05, 0.5, 'sin(latitude)', rotation='vertical', va='center')
        return fig
        
    def resid_vs_dec(self, ib=0, ax=None, ylim=(-8,8), labels=True):
        """ residual vs. declination angle
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
        
    def ridge_spectral_residuals(self, ax=None, glat=(-10,10), glon=(-60,60), nolabels=False):
        """ Spectral residuals along the Galactic ridge
       
        Designed to match, except for different ROI definitions, the Saclay standard plot.
        """
        if ax is None:
            fig,ax = plt.subplots( figsize=(4,4))
        else: fig = ax.figure
        def cut( x, range):
            return (x>range[0])*(x<range[1])
        #ridge = ( np.abs(self.rois.glat)<10) * ( np.abs(self.rois.glon)<60 )
        ridge = cut(self.rois.glat, glat) * cut(self.rois.glon, glon)
        data =self.counts['observed'][ridge].sum()
        model = self.counts['total'][ridge].sum()
        x = self.energy
        y = data/model-1
        yerr = 1/np.sqrt(model) # needs correlation factor
        ax.errorbar(x, y, yerr=yerr, fmt='o')
        plt.setp(ax, xscale='log')
        if not nolabels:
            plt.setp(ax, xlabel='energy (MeV)', ylabel='(counts-model)/model')
        ax.axhline(0, color='gray')
        ax.grid()
        return fig

    def all_plots(self):
        """ Plots generated for any iteration, reflecting quality of counts histogram""" 
        self.runfigures([
            self.chisq_plots,
            self.residual_maps, 
            self.residual_plot, 
            self.residual_hists, 
            self.ridge_spectral_residuals,
            ])
        

class FrontBackSedPlots(Diagnostics):
    require = 'sedinfo.zip'
    """
    """
    def setup(self):
        """
        Unpack the pickles, one per source, into convenient DataFrame objects
        """
        try:
            files, pkls = self.load_pickles('sedinfo')
        except:
            raise Exception( 'No sedinfo files found: must run stage "sedinfo"')
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
        """ Asymmety maps
        Shows the front/back flux asymmetry for sources, vs. the energy flux, for each band.
        Red indicates large galactic diffuse flux, green small.
        """
        map(self.asym_plot,  range(8), self.multifig()); 
        self.multilabels('flux (eV/cm**2/s)','front/back asymmery','Asymmetries for all sources');
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
        """ Front-Back consistency
        Measure of the likelihood ratio test for front/back consistency.
        """
        map(self.consistency_plot, range(8), self.multifig()); 
        self.multilabels('flux (eV/cm**2/s)','front/back asymmery','Asymmetries for all sources');
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
      
    def fb_flux_vs_energy(self):
        """ Front-Back flux vs energy
        The front/back ratio for each of the four strongest soucces
        """
        self.vals = map(self.ratio_fit, range(8), self.multifig())
        plt.suptitle('Front/back flux ratios for strong sources')
        return plt.gcf()
        
    def fb_summary(self):
        """Front/Back flux vs energy summary.
        Weighted average of the four strongest sources.
        """
        fig, ax = self.get_figure( None)
        vals = self.vals # fb_flux_vs_energy first
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
        return fig
        
    def ts_hist(self, ib=0,  space=np.logspace(1,3,21), **kwargs):
        """TS histogram """
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
    
    def all_plots(self):
        """ Analysis of a special "sedinfo" run, which records SED information for all sources
        with fits to front and back only, as well as both.
        """
        self.runfigures([self.asym_plots, self.consistency_plots, self.fb_flux_vs_energy, self.fb_summary])


class ROIinfo(Diagnostics):
    """ setup pickled DataFrame for ROI info.
    roi name is index
    columns as the individual ROI, except exclude name itself, and enter only list of source names for sources
    """
    require=  'pickle.zip'
    def setup(self, **kwargs):
        self.plotfolder='rois'
        self.title='ROI summary'
        self.source_name='observed' # default for base class
        self.plots_kw={}
      
        filename = 'rois.pickle'
        refresh = kwargs.pop('refresh', not os.path.exists(filename) 
                    or os.path.getmtime(filename)<os.path.getmtime('pickle.zip') )
        if refresh:
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
                ra = pkl['skydir'].ra()
                dec= pkl['skydir'].dec()
                tdict.update(glon = glon, glat=glat, ra=ra, dec=dec )
                rdict[pkl['name']] = tdict
            self.df = pd.DataFrame(rdict).transpose()
            self.df.save(filename)
            print 'saved %s' % filename
        else:
            print 'loading %s' % filename
            self.df = pd.load(filename)
        # move this into refresh?
        rois = self.df
        rx = rois['ra dec glat glon'.split()] 
        rx['chisq'] = [r['chisq'] for r in rois['counts']]
        rx['npar'] = [len(p) for p in rois.parameters]
        rx['ring'] = [10**p[0] for p in rois.parameters]
        rx['iso']  = [10**p[1] for p in rois.parameters]
        rx.to_csv('rois.csv')
        print 'saved rois.csv'
        
        self.energy=self.df.ix[0]['counts']['energies']
        self.funcs = []
        self.fnames=[]
        
    def default_plots(self):
        # the set of defaults plots to generate: subclasses can add
        self.funcs = [self.counts_map, self.normalization, self.norm_unc, self.norm_vs_dec]
        self.fnames= map(lambda s: self.source_name+'_'+s, ['counts', 'normalization', 'norm_unc',  'norm_vs_dec'])

        
    def skyplot(self, values, ax=None, title='', ecliptic=False,
                    labels=True, colorbar=True, cbtext='', **scatter_kw):
        if ax is None:
            # note that different size might mean different marker size
            fig, ax = plt.subplots( 1,1, figsize=(6,5))
        else: fig = ax.figure
        singlat=np.sin(np.radians(list(self.df.glat)))
        scat_default = dict(s=60, marker='D', vmin=None, vmax=None, edgecolor='none')
        scat_default.update(scatter_kw)
        
        scat =self.basic_skyplot(ax, self.df.glon, singlat, c=values,
            labels=labels, colorbar=colorbar, cbtext=cbtext, **scat_default)
        
        return fig, scat # so can add colorbar later

    def model_counts(self, name, ib=0):
        """ 
        return an array, in ROI order, of the counts corresponding to the name
            name: string
                either 'observed', the name of a diffuse model, or 'sources' for the total
            ib : int or None
                the energy bin number, or return a sum if None
        """
        def select_band(x):
            return x[ib] if ib is not None else x.sum()
        if name=='observed':
            return np.array([ select_band(self.df.ix[i]['counts']['observed']) for i in range(1728)])
        def source_counts(i):
            m = self.df.ix[i]['counts']['models']
            dn = self.df.ix[i]['diffuse_names']
            r=0
            for nm, cnts in m:
                if nm in dn: continue
                r+= select_band(cnts)
            return r 
        if name=='sources':
            return np.array(map(source_counts, range(1728)))
            
        def cts(i):
            m = self.df.ix[i]['counts']['models']
            dn = self.df.ix[i]['diffuse_names']
            k = dn.index(name) if name in dn else -1
            return select_band(m[k][1]) if k>=0 else 0
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

    def counts_map(self,  hsize= (1.0, 1.7, 1.0, 1.4),  **kwargs):
        """ counts map for %(title)s
        Left:  each ROI, the total counts corresponding to the %(title)s component, 
        for %(energy_selection)s MeV. There are %(total_counts)s total.
        <br>Right: the fraction, %(total_fraction)s percent of the total.
        """
        ib = kwargs.pop('ib', None)
        if 'cb_kw' not in kwargs:
            kwargs['cb_kw'] = dict(shrink=0.70, anchor=(-1.0,0.5)) #ugly tweaking
        self.energy_selection= 'E=%.0f' %self.energy[ib] if ib is not None else 'E>100'
        sm = self.model_counts(self.source_name, ib)
        tot = self.model_counts('observed', ib)
        self.total_counts ='{:,d}'.format(int(np.sum(sm)))
        self.total_fraction = '%.1f' % (np.sum(sm)/np.sum(tot) * 100.)
        
        fig, ax = self.subplot_array(hsize, figsize=(12,6))
        def left(ax):
            if kwargs.pop('log', True):
                c,cbtext = np.log10(sm), 'log10(counts)'
            else: 
                c,cbtext = sm, 'counts'
            if ib is not None:
                return self.skyplot(c, ax=ax, cbtext=cbtext,  **Kwargs)
                 #   title='%s counts at %d MeV' % ( title, self.energy[ib]), **kwargs)
            else:
                return self.skyplot(c, ax=ax, cbtext=cbtext, **kwargs) #title=title,  **kwargs)
        def right(ax):
            if ib is not None:
                return self.skyplot(100*(sm/tot), ax=ax, cbtext='fraction (%)', **kwargs) #title='%s count fraction at %d MeV' % (title, self.energy[ib]),
            else:
                return self.skyplot(100*(sm/tot), ax=ax, cbtext='fraction (%)', **kwargs) #title=title,  cbtext='fraction (%)', **kwargs)
        for f, ax in zip((left,right), ax.flatten()): f(ax)
        return fig
            
#    def count_fraction(self,  title='', **kwargs):
#        """ Count Fraction for %(title)s
#        For each ROI, the fraction of %(title)s counts, 
#        for %(energy_selection)s MeV.
#        """
#        ib = kwargs.pop('ib', None)
#        self.energy_selection= 'E=%.0f' %self.energy[ib] if ib is not None else 'E>100'
#
#        sm = self.model_counts(self.source_name, ib)
#        tot = self.model_counts('observed', ib)
#        if ib is not None:
#            return self.skyplot(100*(sm/tot), title='%s count fraction at %d MeV' % (title, self.energy[ib]),
#                cbtext='fraction (%)', **kwargs)
#        else:
#            return self.skyplot(100*(sm/tot), title=title,  cbtext='fraction (%)', **kwargs)

    def skyplot_with_hist(self, values, xlabel, vmin, vmax, clip,  **kw):
    
        fig, ax = self.subplot_array( (1.0, 0.6, 1.5, 0.7), figsize=(10,5))
        stats = kw.pop('stats', True)
        def hist(ax):
            t = 'count %d\nmean  %.2f\nstd   %.2f'  %(len(values),values.mean(),values.std())
            ax.hist(values.clip(*clip), np.linspace(clip[0],clip[1], 51), label=t)
            ax.grid()
            ax.axvline(1.0, color='k')
            plt.setp(ax, xlabel=xlabel, xlim=clip)
            ax.legend( prop=dict(size=10))
                
        hist(ax[0,0])
        self.skyplot(values, ax=ax[0,1], vmin=vmin, vmax=vmax,  **kw)
        return fig

    
    def normalization(self, vmin=0.8, vmax=1.2, clip =(0.5,1.5), **kw):
        """ normalization factor for %(title)s
        The normalization should be nominally 1.0.
        Left: histogram: right: map
        """ 
        models = self.diffuse_models(self.source_name)
        norms = np.array([m.getp(0) if m is not None else np.nan for m in models])
        return self.skyplot_with_hist(norms, 'normalization', vmin, vmax, clip, **kw)
        
    def norm_unc(self, vmin=0, vmax=0.05, clip =(0, 0.05), **kw):
        """ normalization uncertainty for %(title)s
        The normalization should be nominally 1.0.
        Left: histogram: right: map
        """ 
        models = self.diffuse_models(self.source_name)
        unc = np.array([np.sqrt(m.get_cov_matrix()[0,0]) if m is not None else np.nan for m in models])
        return self.skyplot_with_hist(unc, 'normalization uncertainty', vmin, vmax, clip, **kw)

    def norm_vs_dec(self, vmin=0, vmax=90, size=15, ylim=(0.7,1.2), **kw):
        """ Normalization factor for %(title)s vs Dec
        The color represents the absolute Galactic latitude
        """
        models = self.diffuse_models(self.source_name)
        norms = [m.getp(0) if m is not None else np.nan for m in models]
        sindec = np.sin(np.radians(np.array(self.df.dec,float)))

        fig,ax = plt.subplots(figsize=(6,5))
        c=np.abs(self.df.glat.values)
        cut= c>vmin
        defaults =dict(edgecolors='none', s=size)
        scat=ax.scatter( sindec, norms, c=c, vmin=vmin, vmax=vmax, **defaults)
        plt.setp(ax, xlim=(-1,1), ylim=ylim, xlabel='sin(dec)', ylabel='normalization')
        ax.grid()
        cb =fig.colorbar(scat, ax=ax)
        cb.set_label('abs(b)')
        return fig
        
    def all_plots(self): #, other_html=None):
        """ ROI-based plots, for %(title)s diffuse component. These are based on the 
            individual ROI fits and examine only the normalization factor. See the spectral information, if present, for more
            information about the consistency of the model for this component.
        """
        self.runfigures(self.funcs, self.fnames, **self.plots_kw)


class Exposure(ROIinfo):
    
    def setup(self, **kw):
        super(Exposure, self).setup(**kw)
        self.plotfolder='exposure'
        # use the fact that the isotopic diffuse compoenent is isotropic, so that
        # the ratio of the computed counts, to the fit normalization, is proportional
        # to the exposure.
        iso = self.model_counts('isotrop')
        models = self.diffuse_models('isotrop')
        norms = np.array([m.getp(0) if m is not None else np.nan for m in models])
        self.relative_exp = iso/norms/(iso/norms).mean()
        
    def exposure_plots(self, hsize=(1.0,1.0,2.0,1.0, 2.0, 0.7),):
        """ exposure dependence
        Examine the relative exposure, per ROI. Express in terms of the mean. Note that
        ROIs are distributed uniformly over the sky.
        <br>Left: histogram, center: scatter plot vs. Declination; right: map on sky, in Galactic coordinates.
        
        """
        #fig, axx = plt.subplots(1,3, figsize=(15,4))
        fig, axx = self.subplot_array(hsize, figsize=(12,4))
        relative_exp = self.relative_exp
        label = 'exposure relative to mean'
        lim = (0.7, 1.6)
        def left(ax):
            ax.hist(relative_exp, np.linspace(*lim, num=25))
            plt.setp(ax, xlim=lim)# xlabel=label)
            ax.axvline(1.0, color='k')
            ax.grid()

        def center(ax):
            ax.plot(self.df.dec, relative_exp, '.')
            ax.grid()
            plt.setp(ax, xlim=(-90,90), xlabel='Dec (deg)',ylabel=label, ylim=lim)
            ax.set_xticks(range(-90,91,30))
            ax.axhline(1, color='k')
        def right(ax):
            self.skyplot(relative_exp, ax=ax, s=40)
        
        for f,ax in zip((left, center, right), axx.flatten()): f(ax)
        return fig
        
    def all_plots(self, **kw):
        """ Plots associated with the exposure"""
        self.runfigures([self.exposure_plots])
    

class SunMoon(ROIinfo):
    def setup(self, **kwargs):
        super(SunMoon, self).setup(**kwargs)
        self.plotfolder='sunmoon'
        self.source_name='SunMoon'
        self.title='Sun/Moon'
        t = np.any([x is not None for x in self.diffuse_models(self.source_name)])
        assert t, 'No sun-moon component in this sky model'
        self.default_plots()
        self.plots_kw=dict(ecliptic=True)
        

class Limb(ROIinfo):
    def setup(self, **kwargs):
        super(Limb, self).setup(**kwargs)
        self.plotfolder='limb'
        self.source_name='limb'
        self.title='Limb'
        dm = self.diffuse_models('limb')
        self.fpar,self.bpar = [np.array([m[i] if m else np.nan for m in dm] )for i in range(2)]

        self.default_plots()
        self.funcs += [self.bpar_plot, self.fpar_plot, self.flux_vs_dec, ]
        self.fnames+= ['limb_polar_back', 'limb_polar_front', 'limb_flux_vs_dec',]

     
    def polar_plots(self, values, title=None,
                vmin=0, vmax=2, vticks=5, vlabel=None, thetamax=60):
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
        
        #dm = self.diffuse_models('limb')
        #fpar,bpar = [np.array([m[i] if m else np.nan for m in dm] )for i in range(2)]
        
        
        fig, axx = plt.subplots(2,1, figsize=(8,6))
        plt.subplots_adjust(right=0.9)
        c=np.abs(self.df.glat)
        cut = c>35
        for ax, par, label  in zip(axx, [self.fpar,self.bpar], 'front back'.split()):
            scat=ax.scatter(np.sin(np.radians(dec))[cut], par[cut], c=c[cut], s=30,vmin=0, vmax=90, edgecolor='none')
            plt.setp(ax, xlim=(-1,1), xlabel='sin(dec)' if label=='back' else '',  ylim=(0,2))
            ax.plot(dom, map(limbfun[label],dom), '--', color='k', lw=1) 
            ax.grid()
            ax.text(-0.75, 1.6, label, fontsize=18)
        fig.text(0.05, 0.5, 'flux normalization factor', rotation='vertical', va='center')
        cax = fig.add_axes((0.94, 0.25, 0.02, 0.4))
        cb=plt.colorbar(scat, cax)
        cb.set_label('abs(glat)')
        plt.suptitle('Limb observed flux') 
        return fig
    
    def bpar_plot(self):
        """ Back normalization
        """
        return self.polar_plots(self.bpar)
    def fpar_plot(self):
        """ front normalization
        """
        return self.polar_plots(self.fpar)

class LimbRefit(Limb):
    def setup(self, **kw):
        self.plotfolder = 'limb_refit'
        self.source_name='limb'
        self.title='Limb refit'
        files, pickles = self.load_pickles('limb')
        rdict = dict()
        for f,p in zip(files,pickles):
            name = os.path.split(f)[-1][:9]
            front, back = p['model'].parameters if 'model' in p else (np.nan, np.nan)
            if 'model' not in p: p['model']= None
            ra, dec = p['ra'], p['dec']
            skydir = SkyDir(ra,dec)
            glat,glon = skydir.b(), skydir.l()
            rdict[name] = dict(ra=ra, dec=dec, glat=glat, glon=glon,skydir=skydir, front=front, back=back)
        self.df = pd.DataFrame(rdict).transpose()
        self.fpar = self.df.front
        self.bpar = self.df.back
        
    def all_plots(self):
        """ Special run to refit Limb normalization
        %(table)s
        """
        self.table = pd.DataFrame([self.df.front, self.df.back, ], 
                index=['front', 'back']).T.describe().to_html()
        self.runfigures([self.flux_vs_dec], ['limb_fit_norm_vs_dec'] )
        
class SunMoonRefit(ROIinfo):
    def setup(self, **kw):
        self.plotfolder = 'sunmoon_refit'
        self.title = 'SunMoon refit'
        files, pickles = self.load_pickles('sunmoon')
        rdict = dict()
        for f,p in zip(files,pickles):
            name = os.path.split(f)[-1][:9]
            model = p['model'] if 'model' in p else None
            norm = model.parameters[0] if model is not None else np.nan
            norm_unc = np.diag(model.get_cov_matrix())[0]**0.5 if model is not None else np.nan
            if 'model' not in p: p['model']= None
            ra, dec = p['ra'], p['dec']
            skydir = SkyDir(ra,dec)
            glat,glon = skydir.b(), skydir.l()
            if glon>180: glon-=360.
            rdict[name] = dict(ra=ra, dec=dec, glat=glat, glon=glon, skydir=skydir, 
                norm=norm, norm_unc=norm_unc,
                delta_likelihood = p['delta_likelihood'])
        self.df = pd.DataFrame(rdict).transpose()
    
    def sunmoon_normalization(self):
        """ Sun/Moon normalization 
        Note that it is defined for all directions
        """
        self.skyplot_with_hist(self.df.norm, 'norm', 0.5, 1.5, (0.5,1.5))
        
    def sunmoon_loglikelihood(self):
        """ Sun/Moon log likelihood change
        The improvement in the log likelihood when Sun/Moon normalization is freed
        """
        self.skyplot_with_hist(self.df.delta_likelihood, 'delta loglike', 0, 5, (0,5))

    def all_plots(self):
        """ Refit to SunMoon model, showing normalization, change in likelihood
        %(table)s"""
        self.table = pd.DataFrame([self.df.norm, self.df.norm_unc,self.df.delta_likelihood ], 
                index=['norm', 'norm_unc', 'delta_likelihood']).T.describe().to_html()
        self.runfigures([self.sunmoon_normalization, self.sunmoon_loglikelihood])
        

class Galactic(ROIinfo):
    def setup(self, **kw):
        super(Galactic, self).setup(**kw)
        self.plotfolder='gal'
        self.source_name='ring'
        self.title='Galactic'
        self.default_plots()


class Isotropic(Galactic):

    def setup(self, **kw):
        super(Isotropic, self).setup(**kw)
        self.plotfolder='iso'
        self.source_name='isotrop'
        self.title='Isotropic'
        self.default_plots()
        self.funcs += [self.isotropic_spectrum]
        self.fnames +=['isotropic_spectrum']
        # look up filenames used to define the isotorpic spectrum: either new or old diffuse spec; list or dict
        config = eval(open('config.txt').read())
        diffuse=config['diffuse']
        isokey = 'isotrop' if type(diffuse)==types.DictType else 1
        self.idfiles = [os.path.join(os.environ['FERMI'],'diffuse',diffuse[isokey][i]) for i in (0,1)]
        
    def isotropic_spectrum(self, other=None):
        """ Isotropic Spectrum from template
        
        The spectrum used to define the isotropic diffuse component.
        <br>Files for front/back: %(idfiles)s
        """
        nf,nb = map(np.loadtxt, self.idfiles)
        energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
        fig, axs = plt.subplots(1,2, figsize=(7,3), dpi=50)
        def right(ax):
            ax.plot(energies, front/back, '-o');
            ax.axhline(1.0, color='k')
            plt.setp(ax, xscale='log', xlabel='Energy');ax.grid(True);
            ax.set_title('Isotropic flux front/back ratio', fontsize='small');
        def left(ax):
            ax.plot(energies, front*energies**2, '-g', label='front')
            ax.plot(energies, back*energies**2, '-r', label='back')
            plt.setp(ax, xlabel='Energy', ylabel='flux*e**2', xscale='log')
            ax.set_title('isotropic diffuse spectra', fontsize='small')
            ax.grid(True); ax.legend()
        for f,a in zip((left,right), axs.flatten()): f(a)
        return fig
     
    def combined_spectra(self, other='isotrop_4years_P7_v9_repro_data_source.txt'):
        """ Special isotropic spectrum plot"""
        nf,nb = map(np.loadtxt, self.idfiles)
        fig, ax = plt.subplots(figsize=(5,5))
        energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
        ax.plot(energies, front*energies**2, '-g', label='1year_P76R_source_front')
        ax.plot(energies, back*energies**2, '-r', label='1year_P76R_source_back')
        plt.setp(ax, xlabel='Energy', ylabel='flux*e**2', xscale='log')
        ax.set_title('isotropic diffuse spectra', fontsize='small')
        f =os.path.join(os.path.expandvars('$FERMI/diffuse'),other)
        assert os.path.exists(f)
        iso = np.loadtxt(f)
        energies, fluxes =iso[:,0], iso[:,1] 
        ax.plot(energies, fluxes*energies**2, '-', label ='4years_P7-v9_source')
        plt.setp(ax, xlim=(1e2,1e4))
        ax.legend();ax.grid(True)
        return fig

    
class SourceTotal(ROIinfo):
    def setup(self, **kw):
        super(SourceTotal, self).setup(**kw)
        self.plotfolder='sourcetotal'
        self.source_name='sources'
        self.title='Sources'
        self.funcs = [self.counts_map]
        self.fnames=['source_counts']

    def all_plots(self, **kwargs):
        """ Counts for all sources, per RIO"""
    
        self.runfigures([self.counts_map], ['source_counts'], **kwargs)


class SourceInfo(Diagnostics):
    """ To be superclass for specific source plot stuff, creates or loads
        a DataFrame with all sources 
        """
    require='pickle.zip'
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
                        d = np.diag(model.get_cov_matrix())
                        d[d<0] =0
                        errs[:n] = np.sqrt(d)
                        errs[np.isnan(errs)]=-1
                        badfit = np.any(errs[model.free]<=0)
                    except Exception, msg:
                        print 'fail errors for %s:%s' % (name, msg)
                        badfit = True
                    ellipse = info.get('ellipse', None)
                    sdict[name] = info
                    pulsar = model.name.endswith('Cutoff')
                    betavalue = float(pars[2]) if not pulsar else np.nan
                    sdict[name].update(
                        glat=info['skydir'].b(), glon=info['skydir'].l(),
                        roiname=pkl['name'], 
                        pars= pars, errs=errs, free=free, badfit=badfit,
                        a = ellipse[2] if ellipse is not None else np.nan,
                        b = ellipse[3] if ellipse is not None else np.nan,
                        ang=ellipse[4] if ellipse is not None else np.nan,
                        locqual = round(ellipse[5],2) if ellipse is not None else np.nan,
                        delta_ts = ellipse[6] if ellipse is not None else np.nan,
                        freebits= np.sum( int(b)*2**i for i,b in enumerate(model.free)),
                        flux = pars[0],
                        flux_unc = errs[0],
                        pindex = pars[1],
                        pindex_unc = errs[1],
                        beta = betavalue,
                        beta_unc = errs[2] if not pulsar and pars[2]>0.002 else np.nan,
                        index2 = pars[3] if pulsar else pars[2],
                        index2_unc = errs[3] if pulsar and not np.isnan(pars[3]) else errs[2],
                        cutoff = pars[2] if pulsar else np.nan,
                        cutoff_unc = errs[2] if pulsar else np.nan,
                        e0 = model.e0,
                        modelname=model.name,
                        fitqual = round(sum(info['sedrec'].delta_ts),2),
                        eflux = pars[0]*model.e0**2*1e6,
                        psr = pulsar,
                        )
            df = pd.DataFrame(sdict).transpose()
            df.index.name='name'
            df['name'] = df.index
            ra = [x.ra() for x in df.skydir]
            dec = [x.dec() for x in df.skydir]
            df['ra'] = ra
            df['dec'] = dec
            self.df = df.sort_index(by='ra')
            self.df.save(filename)
            print 'saved %s' % filename

        else:
            print 'loading %s' % filename
            self.df = pd.load(filename)
        self.df['flux']    = [v[0] for v in self.df.pars.values]
        self.df['flux_unc']= [v[0] for v in self.df.errs.values]
        localized = ~np.array(pd.isnull(self.df.delta_ts))
        extended = np.array(self.df.isextended, bool)
        self.df['unloc'] = ~(localized | extended)
        unloc = self.df.unloc * (self.df.ts>10)
        self.df['poorloc'] = (self.df.a>0.2) + (self.df.locqual>8) + (self.df.delta_ts>2)
        if sum(unloc)>0:
            print '%d point sources (TS>10) without localization information' % sum(unloc)
            self.df.ix[unloc]['ra dec ts roiname'.split()].to_csv('unlocalized_sources.csv')
            print 'Wrote file "unlocalized_sources.csv"'
        self.df['flags'] = 0  #used to set bits below
        flags = self.df.flags
        flags[self.df.poorloc + self.df.unloc] += 4 ### bit 4
        print '%d sources flagged (4) as poorly or not localized' % sum(self.df.poorloc + self.df.unloc)

 
        self.energy = np.sqrt( self.df.ix[0]['sedrec'].elow * self.df.ix[0]['sedrec'].ehigh )
            
    def skyplot(self, values, proj=None, ax=None, ecliptic=False,
                labels=True, title='', colorbar=True, cbtext='', **scatter_kw):
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
        
        scatter_kw_default=dict(s=20, vmin=None, vmax=None, edgecolor='none')
        scatter_kw_default.update(scatter_kw)
        
        scat = self.basic_skyplot(ax, glon, singlat, c=c,
                title=title, ecliptic=ecliptic, colorbar=colorbar,cbtext=cbtext, **scatter_kw_default)
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
        return pd.DataFrame(dict(fdata=fdata, udata=udata, ldata=ldata, fmodel=fmodel, 
                glat=s.glat, glon=s.glon, roiname=s.roiname),
            index=s.index).sort_index(by='roiname')

    def cumulative_ts(self, ts=None, other_ts=None, tscut=(10,25), check_localized=True, label=None, otherlabel=None):
        """ Cumulative test statistic TS
        
        A logN-logS plot, but using TS. Important thresholds at TS=10 and 25 are shown.
        """
        usets = self.df.ts if ts is None else ts
        df = self.df
        fig,ax = plt.subplots( figsize=(8,6))
        dom = np.logspace(np.log10(9),5,1601)
        ax.axvline(25, color='gray', lw=1)
        ax.hist( usets ,dom, cumulative=-1, lw=2, color='g', histtype='step',label=label)
        if other_ts is not None:
            ax.hist( other_ts ,dom, cumulative=-1, lw=2, color='b', histtype='step',label=otherlabel)
        if check_localized:
            unloc = df.unloc
            ul = df[(unloc+df.poorloc) * usets>tscut[0]] 
            n = len(ul)
            if n>10:
                ax.hist(ul.ts ,dom, cumulative=-1, lw=2, color='r', histtype='step',
                    label='none or poor localization')
                ax.text(12, n, 'none or poor localization (TS>%d) :%d'%(tscut[0],n), fontsize=12, color='r')
        plt.setp(ax,  ylabel='# sources with greater TS', xlabel='TS',
            xscale='log', yscale='log', xlim=(9, 1e4), ylim=(9,8000))
        ax.set_xticklabels([' ', '10', '100', '1000'])
        ax.set_yticklabels(['', '10', '100', '1000'])
            
        # label the plot with number at given TS
        for t in tscut:
            n = sum(usets>t) 
            ax.plot([t,2*t], [n,n], '-k');
            ax.plot(t, n, 'og')
            ax.text(2*t, n, 'TS>%d: %d'%(t,n), fontsize=14, va='center')
                
        ax.grid()
        return fig

    def cumulative_counts(self):
        #assume rois.pickle available
        # extract the model counts for each source from it, and add to the DF
        if not hasattr(self,'roi_df'):
            self.roi_df = pickle.load(open('rois.pickle'))
        def counts(src_df,roi_df, name):
            """return fit counts for source name"""
            roiname = src_df.ix[name]['roiname']
            roi=roi_df.ix[roiname]
            names = list(roi['counts']['names'])
            i = names.index(name)
            if i<0: print name, i
            try: 
                t =np.sum(roi['counts']['models'][i][1])
                return t
            except:
                print 'fail to get counts for source', name, roiname, i
                return 0
        self.df['counts'] = [counts(self.df, self.roi_df, name) for  name in self.df.index]
        

    def non_psr_spectral_plots(self, index_min=1.0, index_max=3.5, beta_max=2.0):
        """ Plots showing spectral parameters for non-pulsar spectra
        Left: energy flux in eV/cm**2/s. This is the differential flux at the pivot energy
        <br> Center: the spectral index.
        <br> Right: the curvature index for the subset with log parabola fits.
        %(tail_check)s
        """
        fig, axx = plt.subplots( 1,3, figsize=(12,4))
        t = self.df.ix[(self.df.ts>10)*(self.df.modelname=='LogParabola')]['ts flux pindex beta e0 roiname'.split()]
        t['eflux'] = t.flux * t.e0**2 * 1e6
        ax = axx[0]
        [ax.hist(t.eflux[t.ts>tscut].clip(4e-2,1e2), np.logspace(-2,2,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xscale='log', xlabel='energy flux', xlim=(4e-2,1e2)); ax.grid(); ax.legend(prop=dict(size=10))
        ax = axx[1]
        [ax.hist(t.pindex[t.ts>tscut].clip(index_min,index_max), np.linspace(index_min,index_max,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xlabel='spectral index'); ax.grid(); ax.legend(prop=dict(size=10))
        ax = axx[2]
        [ax.hist(t.beta[(t.ts>tscut)*(t.beta>0.01)].clip(0,beta_max), np.linspace(0,beta_max,26), label='TS>%d' % tscut) for tscut in [10,25] ]
        plt.setp(ax, xlabel='beta'); ax.grid(); ax.legend(prop=dict(size=10))
        # get tails
        tail_cut = (t.eflux<5e-2)+((t.pindex<index_min)+(t.pindex>index_max))*t.beta.isnull()+(t.beta>beta_max)
        self.non_psr_tails=t[tail_cut]['ts eflux pindex beta roiname'.split()]
        if len(self.non_psr_tails)>0:
            self.tail_check = '<h4>Sources on tails:</h4>'\
                +self.non_psr_tails.sort_index(by='roiname').to_html()
            self.tail_check += '<p>Criteria: require index between 1 and 3.5 for powerlaw, beta<2.0 for log parabola'
        else:
            self.tail_check ='<p>No sources on tails'
        # flag sources
        flags = self.df.flags
        tails = self.non_psr_tails.index
        flags[tails] += 1 ### bit 1
        print '%d sources flagged (1) in tails of flux, index, or beta' % len(self.non_psr_tails)
        return fig
    
    def pulsar_spectra(self, index_min=0.5, index_max=2.5, cutoff_max=10000):
        """ Distributions for the LAT pulsars
        
        For each plot, the subset with a bad fit is shown.
        """
        fig, axx = plt.subplots( 1,4, figsize=(14,4))
        plt.subplots_adjust(wspace=0.3)
        t = self.df.ix[(self.df.ts>10)*(np.array(self.df.psr,bool))]['ts flux pindex cutoff e0 roiname fitqual'.split()]
        t['eflux'] = t.flux * t.e0**2 * 1e6
        badfit = t.fitqual>50

        def plot1(ax, efmin=1e-1,efmax=1e3):
            bins = np.logspace(np.log10(efmin),np.log10(efmax),26)
            vals = t.eflux.clip(efmin,efmax)
            ax.hist(vals, bins )
            ax.hist(vals[badfit], bins, color='red', label='bad fit')
            plt.setp(ax, xscale='log', xlabel='energy flux', xlim=(efmin,efmax)); ax.grid(); 
            ax.legend(prop=dict(size=10))

        def plot2(ax):
            bins = np.linspace(index_min,index_max,26)
            vals = t.pindex.clip(index_min,index_max)
            ax.hist(vals, bins)
            ax.hist(vals[badfit], bins, color='red', label='bad fit')
            plt.setp(ax, xlabel='spectral index'); ax.grid(); 
            ax.legend(prop=dict(size=10))
            
        def plot3(ax):
            bins = np.linspace(0,cutoff_max/1e3,26)
            vals = t.cutoff.clip(0,cutoff_max) /1e3
            ax.hist(vals, bins)
            ax.hist(vals[badfit], bins, color='red', label='bad fit')
            plt.setp(ax, xlabel='cutoff energy (GeV)'); ax.grid()
            ax.legend(prop=dict(size=10))
            
        def plot4(ax):
            xvals = t.cutoff.clip(0,cutoff_max) / 1e3
            yvals = t.pindex.clip(index_min,index_max)
            ax.plot(xvals, yvals, 'o')
            ax.plot(xvals[badfit], yvals[badfit], 'or', label='bad fit')
            plt.setp(ax, xlabel='cutoff energy', ylabel='spectral index')
            ax.grid(); ax.legend(prop=dict(size=10))

            
        for f,ax in zip((plot1,plot2,plot3,plot4,), axx.flatten()): f(ax)
        flags = self.df.flags
        tail_cut = (t.pindex<=index_min) + (t.pindex>index_max) + (t.cutoff>cutoff_max)

        tails = t.ix[tail_cut].index
        flags[tails] += 1 ### bit 1
        print '%d pulsar sources flagged (1) in tails of  index or cutoff' % sum(tail_cut)

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
        
    def fit_quality(self, xlim=(0,50), ndf=10, tsbandcut=20):
        """ Fit quality
        This is the difference between the TS from the fits in the individual energy bands, and that for the spectral fit.
        It should be distributed approximately as chi squared of at most 14-2 =12 degrees of freedom. 
        However, high energy bins usually do not contribute, so we compare with ndf=%(ndf)d.
        All sources with TS_bands>%(tsbandcut)d are shown.<br>
        Left: Power-law fits. Tails in this distribution perhaps could be 
        improved by converting to log parabola. 
        <br>Center: Log parabola fits.
        <br>Right: Fits for the pulsars, showing high latitude subset.
        <br> Averages: %(fit_quality_average)s
        %(badfit_check)s
        %(poorfit_table)s

        """
        from scipy import stats
        fig, axx = plt.subplots(1,3, figsize=(12,5))
        s = self.df
        psr = np.asarray(s.psr, bool)
        beta = s.beta
        logparabola = (~psr) * (beta>0.01)
        powerlaw = (~psr) * (beta.isnull() + (beta<0.01) )

        self.tsbandcut=tsbandcut
        cut=s.band_ts>tsbandcut
        
        dom = np.linspace(xlim[0],xlim[1],26)
        d = np.linspace(xlim[0],xlim[1],51); delta=dom[1]-dom[0]
        chi2 = lambda x: stats.chi2.pdf(x,ndf)
        fudge = 1.0 # to scale, not sure why
        hilat = np.abs(self.df.glat)>5
        self.average = [0]*3; i=0
        for ax, label in zip(axx[:2], ('powerlaw', 'logparabola',)):
            mycut=cut*eval(label)
            count = sum(mycut)
            ax.hist(s.fitqual[mycut].clip(*xlim), dom, label=label+' (%d)'%count)
            self.average[i]=s.fitqual[mycut].mean(); i+=1
            ax.plot(d, chi2(d)*count*delta/fudge, 'r', lw=2, label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
            ax.grid(); ax.set_xlabel('fit quality')
            ax.legend(prop=dict(size=10))
            
        def right(ax, label='PSR'):
            mycut = cut * (psr)
            count = sum(mycut)
            ax.hist(s.fitqual[mycut].clip(*xlim), dom, label=label+' (%d)' %count)
            ax.hist(s.fitqual[mycut*hilat].clip(*xlim), dom, label=label+' [|b|>5] (%d)' %sum(mycut*hilat))
            self.average[i]=s.fitqual[mycut].mean()
            ax.plot(d, chi2(d)*count*delta/fudge, 'r', lw=2, label=r'$\mathsf{\chi^2\ ndf=%d}$'%ndf)
            ax.grid();ax.set_xlabel('fit quality')
            ax.legend(loc='upper left', prop=dict(size=10))
        
        right(axx[2])
        self.df['badfit2'] =np.array(self.df.badfit.values, bool)
        t = self.df.ix[(self.df.badfit2)*(self.df.ts>10)].sort_index(by='roiname')
        print '%d sources with bad fits' %len(t)
        if len(t)>0:
            self.badfit = t[['ts', 'errs', 'roiname']]
            self.badfit_check = '<h4>Sources with missing errors:</h4>'+self.badfit.to_html()
        else: self.badfit_check = '<p>All sources fit ok.'
        self.fit_quality_average =  ', '.join( map(lambda x,n :'%s: %.1f' %(n,x) ,
                            self.average, 'powerlaw logparabola expcutoff'.split()) )
        self.ndf=ndf
        print 'fit quality averages:', self.fit_quality_average

        # Make tables (csv and html) of the poor fits

        t =s.ix[s.fitqual>50]['ra dec fitqual ts modelname roiname'.split()].sort_index(by='roiname')
        poorfit_csv = 'poor_spectral_fits.csv'
        t.to_csv(poorfit_csv)
        bs =sorted(list(set(t.roiname)))
        print 'Wrote out list of poor fits to %s, %d with fitqual>50, in %d ROIs' % (poorfit_csv, len(t), len(bs))
        # todo: make a function to do this nidcely
        poorfit_html = self.plotfolder+'/poorfits.html'
        open(poorfit_html,'w').write('<head>\n'+ HTMLindex.style + '</head>\n<body>'+t.to_html()+'\n</body>')
        self.poorfit_table = '<p> <a href="poorfits.html"> Table of %d poor fits, with fitqual>50 </a>' % (  len(t) )
        
        # flag sources that made it into the list
        self.df.flags[t.index] += 2
        print '%d sources flagged (2) as poor fits' %len(t)
        return fig
        
    def pivot_vs_e0(self, xylim=(100, 4e4)):
        """ pivot vs e0
        """
        fig, ax = plt.subplots(figsize=(4,4))
        s = self.df
        cut = s.ts>10
        ax.plot(s.e0[cut].clip(*xylim), s.pivot_energy[cut].clip(*xylim), '.')
        plt.setp(ax, xscale='log',xlabel='e0', xlim=xylim, 
                    ylabel='pivot', yscale='log', ylim=xylim)
        ax.set_title('compare calculated pivot with e0', fontsize=10)
        ax.grid()
        return fig
        
    def fitquality(self):
        """Fit Quality
        
        Left: fit quality histogram; right fit quality vs. TS'

        """
        fig, axs = plt.subplots(1,2, figsize=(7,3))
        plt.subplots_adjust(wspace=0.35)
        s = self.df
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
        return fig

    def flux_uncertainty(self):
        """ flux uncertainty compared with TS
        
        """
        fig, axx = plt.subplots(1,2, figsize=(10,4))
        plots=[]
        relflux_unc= self.df.flux_unc/self.df.flux
        ts = np.asarray(self.df.ts, float)
        ru= np.array(relflux_unc*100.,float)

        def plot1(ax):   
            dom = np.logspace(0,2,26)
            ax.hist(ru[ts>9], dom, label='%d sources'% sum(ts>9))
            for tsmin in (25,100,1000):
                ax.hist(ru[ts>tsmin], dom, label='TS<%d' % tsmin )
            plt.setp(ax, xscale='log', xlabel='relative flux uncertainty (%)', xlim=(1,100))
            ax.set_xticklabels([1,10,100])
            ax.grid()
            ax.legend(loc='upper left', prop=dict(size=10))
        plots.append(plot1)
            
        def plot2(ax):
            ax.plot(ts, ru*np.sqrt(ts)/100, '.')
            plt.setp(ax, xlabel='TS', xlim=(10,10000), xscale='log',
                 yscale='log',ylabel='ru*sqrt(ts)', ylim=(0.8,4))
            ax.plot([0.1,100], [0.1,100],'-g')
            ax.grid()
        plots.append(plot2)
            
        for plotf, ax in zip( (plots), axx.flatten(),):
            plotf(ax)
        return fig
        
    def spectral_fit_consistency_plots(self, energy=133., minflux=2.0, 
            title = 'low energy fit consistency'
        ):
        """ Spectral fit consistency for the lowest energy bin
        
        These plot show information the consistency of the lowest energy band with the spectrm
        defined by the full fit. <br>
        Left: distribution of the "pull" <br>
        Center: data/model ratio with errors, vs. the model flux.<br>
        Right: position in the sky of flagged sources <br>
        """
        cut=self.df.ts>25
        s = self.df[cut]
        fdata = np.array([s.ix[i]['sedrec'].flux[0] for i in range(len(s))])
        udata = np.array([s.ix[i]['sedrec'].uflux[0] for i in range(len(s))])
        ldata = np.array([s.ix[i]['sedrec'].lflux[0] for i in range(len(s))])
        pull = np.array([s.ix[i]['sedrec'].pull[0] for i in range(len(s))])
        fmodel = np.array([s.ix[i]['model'](energy)*energy**2*1e6 for i in range(len(s))])
        glat = np.array([x.b() for x in s.skydir])
        fluxcut = fmodel>minflux
        latcut  = abs(glat)>5.0
        hilat = fluxcut*(latcut)
        lolat = fluxcut*(~latcut)
        
        lowebad = np.abs(pull)>3
        self.df.flags[lowebad] += 4
        print 'Tagged %d sources with lowebad bit' % sum(lowebad)

        y = fdata/fmodel
        ylower, yupper =[(fdata-ldata)/fmodel,(udata-fdata)/fmodel]
        xhi,yhi,yerrhi = fmodel[hilat], y[hilat], [ylower[hilat],yupper[hilat]]
        xlo,ylo,yerrlo = fmodel[lolat], y[lolat], [ylower[lolat],yupper[lolat]]
        
        def error_bar(ax):
            ax.errorbar(x=xhi, y=yhi, yerr=yerrhi, fmt='og', label='%d hilat sources'%sum(hilat))
            ax.errorbar(x=xlo, y=ylo, yerr=yerrlo, fmt='or', label='%d lowlat sources'%sum(lolat))
            plt.setp(ax, xlabel=r'$\mathsf{model\ flux\ (eV\ cm^{-2} s^{-1}})$', xscale='log', 
                ylabel='data/model', ylim=(0,2.5), xlim=(minflux, 100) )
            ax.set_xticks([2,5,10,20,50,100])
            ax.set_title( title, fontsize='medium')
            ax.set_title( title, fontsize='medium')
            ax.legend(loc='upper left', prop=dict(size=10))
            ax.grid()  

        def hist(ax):
            dom = np.linspace(-3,3,26)
            hist_kw=dict(lw=2, histtype='step')
            q=pull.clip(-3,3) 
            
            ax.hist(q[hilat], dom, color='g',  label='%d hilat sources'%sum(hilat),  **hist_kw)
            ax.hist(q[lolat], dom, color='r',  label='%d lowlat sources'%sum(lolat), **hist_kw)
            ax.set_xlabel('pull')
            ax.axvline(0, color='k')
            ax.set_xlim((-3,3))
            ax.set_title( title, fontsize='medium')
            ax.legend(loc='upper left', prop=dict(size=10))
            ax.grid()  

        def skyplot(ax):
            pdf = pd.DataFrame(dict(pull=pull), index=s.index) # to atatch indx
            self.skyplot(pdf.pull[lowebad], ax=ax, vmin=-3, vmax=3, title=title, cbtext='pull')

        fig,ax = plt.subplots(1,3, figsize=(12,4))
        plt.subplots_adjust(wspace=0.3)
        for f, ax in zip( (hist, error_bar, skyplot), ax.flatten()):
            f(ax=ax)
        return fig
        
    def all_plots(self):
        """ Plots of source properties, from analysis of spectral fits. 
        See <a href="../localization/index.html"> localization </a> for localization plots.
        <h3>Census</h3>
        
        %(census_html)s
        <p>
        The columns are the number of sources with ts> than the header. 
        The rows labels are the first three characters of the source name, except 'ext' means extended.
        """
        df = self.df
        extended = np.asarray(df.isextended.values,bool)
        pointsource = ~extended

        def count(prefix, tsmin):
            if prefix=='ext':
                return sum(extended*(df.ts>tsmin))
            elif prefix=='total':
                return sum(df.ts>tsmin)
            names = df[pointsource*(df.ts>tsmin)]['name'].values    
            return sum([n.startswith(prefix) for n in names])
        prefixes = list(set( n[:3] for n in df[pointsource]['name'])) +['ext', 'total']

        census = dict()
        for x in (0, 10, 25):
            census[x] = [count(prefix, x) for prefix in prefixes]
        self.census=pd.DataFrame(census, index=prefixes)
        self.census_html = self.census.to_html()
       
        self.runfigures([self.cumulative_ts, self.fit_quality,self.spectral_fit_consistency_plots,
            self.non_psr_spectral_plots, self.pulsar_spectra, self.pivot_vs_e0,  ]
        )
        
        plt.close('all')
        csvfile='sources.csv'
        colstosave="""ra dec ts modelname freebits fitqual e0 flux flux_unc pindex pindex_unc index2 index2_unc
                 cutoff cutoff_unc locqual a b ang flags roiname""".split()
        self.df.ix[self.df.ts>10][colstosave].to_csv(csvfile)
        print 'saved truncated csv version to "%s"' %csvfile

        t =self.df[self.df.flags>0]['ra dec ts fitqual eflux pindex beta cutoff index2 flags roiname'.split()]
        t.to_csv('flagged_sources.csv')
        print 'wrote %d sources to flagged_sources.csv' % len(t)
        
    def flux_table(self, source_name):
        """ Return a DataFrame describing the sed info for the given source
        Columns, with energy fluxes in eV units, are:
            flux, lflux, uflux : measured flux, lower and upper uncertainty, or 0,0, 95% limit
            mflux : predicted flux for this bin
            TS : Test Statistic for the signal
            pull : signed square root of the 
        """
        si = self.df.ix[source_name]['sedrec']
        pull = np.sign(si.flux-si.mflux) * np.sqrt(si.delta_ts.clip(0,100))
        return pd.DataFrame(dict(flux=si.flux.round(1), TS=si.ts.round(1), lflux=si.lflux.round(1), 
            uflux=si.uflux.round(1), model=si.mflux.round(1), pull=pull.round(1) ),
                index=np.array(np.sqrt(si.elow*si.ehigh),int), columns='flux lflux uflux model TS pull'.split())
 
class Localization(SourceInfo):

    require='pickles.zip'
    def setup(self, **kw):
        super(Localization, self).setup(**kw)
        self.plotfolder = 'localization'
        # unpack the ellipse info into a new DataFrame
        self.ebox = pd.DataFrame([x if x is not None else [np.nan]*7 for x in self.df.ellipse], index=self.df.index)
        self.ebox.columns = 'fit_ra fit_dec a b ang locqual delta_ts'.split()
        self.ebox['roiname']=self.df.roiname
        self.ebox['ts'] = self.df.ts
        self.tscut = kw.get('tscut', 10.)
        self.acut =  kw.get('acut', 0.25)
        self.qualcut=kw.get('qualcut', 8.0)
        self.poorloc=self.ebox[ (self.ebox.locqual>self.qualcut) | (self.ebox.a>self.acut)]['ts a locqual delta_ts roiname'.split()].sort_index(by='ts',ascending=False)
        if len(self.poorloc)>0:
            print '%d poorly localized (locqual>%.1f or a>%.2f) '%\
                (len(self.poorloc), self.qualcut,self.acut)
            self.poorloc.to_csv('poorly_localized.csv')
            print 'wrote file "poorly_localized.csv"'

    def localization(self, maxdelta=9, mints=10):
        """Localization plots
            Left: histogram of the square root of the TS difference from current position to
            the fit; corresponds the number of sigmas. <br>
            Right: scatter plot of this vs. TS
            """
        bins=np.linspace(0,np.sqrt(maxdelta),26)
        fig, axx = plt.subplots(1,2,figsize=(10,5)); 
        plt.subplots_adjust(wspace=0.4)
        wp = self.ebox
        cut = self.df.ts>mints
        ax=axx[0]
        for tcut in (mints, 100):
            t = np.sqrt(wp.delta_ts[self.df.ts>tcut].clip(0,maxdelta))
            ax.hist(t, bins, label='ts>%d: mean=%.2f'%(tcut, t.mean()) )
        #ax.hist(np.sqrt(wp.delta_ts[self.df.ts>100].clip(0,maxdelta)), bins,label='TS>100\nmean:%f.1'%wp.delta)
        ax.legend(prop=dict(size=10))
        ax.grid()
        plt.setp(ax, xlabel='sqrt(delta TS)')
        ax=axx[1]
        ax.plot( self.df.ts[cut],np.sqrt(wp.delta_ts[cut].clip(0,maxdelta)), '.')
        ax.grid()
        plt.setp(ax, xscale='log', xlabel='TS', ylabel='sqrt(delta TS)')
        return fig
        
    def localization_quality(self, maxqual=10, mints=10, tscut=25):
        """Localization quality plots
            Left: histogram of the fit quality. This is a measure of the difference between the sampled
            TS map points and the prediction of the quadratic model. <br>
            Center: scatter plot of the quality vs. TS. <br>
            Right: locations of poorly-fit sources, see the <a href="poorly_localized_table.html">table</a>.
        """
        bins=np.linspace(0,maxqual,26)
        #ig, axx = plt.subplots(1,3,figsize=(13,5));
        fig, axxx = self.subplot_array( hsize=(1.0, 0.6, 1.0, 0.2, 2.0, 0.5), figsize=(13,5))
        axx = axxx[0]
        plt.subplots_adjust(wspace=0.4)
        wp = self.ebox
        cut = self.df.ts>mints
        ax=axx[0]
        for x in (mints, tscut):
            ax.hist(wp.locqual[self.df.ts>x].clip(0,maxqual), bins,label='TS>%d'%x)
        ax.legend(prop=dict(size=10))
        ax.grid()
        plt.setp(ax, xlabel='localization fit quality')
        ax=axx[1]
        ax.plot( self.df.ts[cut],wp.locqual[cut].clip(0,maxqual), '.')
        ax.grid()
        plt.setp(ax, xscale='log', xlim=(10,1e5), xlabel='TS', ylabel='localization fit quality')
        ax=axx[2]
        self.skyplot(self.poorloc.locqual, ax=ax, s=50, vmin=10, vmax=100, cbtext='localization quality')
        return fig 
        
    def r95(self, qualmax=10):
        """ Error circle radius
        R95 is the 95 %%%% containment radius. Here I show the semi-major axis.
        Applying cut quality < %(qualmax)d.
        """
        fig, ax = plt.subplots(1,2, figsize=(12,5))
        r95 = 60* 2.6 * self.ebox.a[self.ebox.locqual<qualmax]
        ts = self.df.ts[self.ebox.locqual<qualmax]
        self.qualmax = qualmax
        def hist(ax, rmax=30):
            bins = np.linspace(0,rmax,31)
            ax.hist(r95.clip(0,rmax), bins, label='all')
            for tsmin in (25,1000):
                ax.hist(r95[self.df.ts>tsmin].clip(0,rmax), bins, label='TS>%d' % tsmin)
            plt.setp(ax, xlabel='R95 (arcmin)')
            ax.grid(); ax.legend(prop=dict(size=10))
        def scat(ax):
            #ax.plot( ts, r95, '.')
            #plt.setp(ax, xscale='log', xlim=(10,1000), ylim=(0,30), xlabel='TS', ylabel='r95')
            ax.plot( ts**-0.5, r95.clip(0,30), '.')
            plt.setp(ax, xlim=(0,0.30),  ylim=(0,30.5), xlabel='1/sqrt(TS)', ylabel='R95 (arcmin)')
            ax.grid()
        for f,a in zip((hist,scat), ax): f(a)
        return fig
        
    def check_closeness(self, tol=0.15, bmin=0, tsmin=10):
        """ Closeness check
            <p>Table of pairs closer than %(close_tol).2f degrees
            %(close_table)s
        """
        cut = (np.abs(self.df.glat)>bmin) * (self.df.ts>tsmin)
        indeces = self.df.index[cut] 
        sdirs = self.df[ cut ]['skydir'].values
        name1=[]; name2=[]; distance=[]
        for i in range(len(sdirs)-1): 
            a = sdirs[i]
            for j in range(i+1, len(sdirs)):
                dist = np.degrees(a.difference(sdirs[j]))
                if dist< tol:
                    name1.append(indeces[i])
                    name2.append(indeces[j])
                    distance.append(dist)
                    print 'Closer than tolerance: sources %s, %s, %.2f deg' \
                        % (indeces[i], indeces[j], dist)
        self.close_tol = tol
        self.close_table = pd.DataFrame(dict(source1=name1, source2=name2, distance=distance)).to_html()
        return None
        

    def source_confusion(self, bmin=10, dtheta=0.1, nbins=50, deficit_angle=1.0, tsmin=10):
        """ Source Confusion
        Distribution of the distances to the nearest neighbors of all detected sources with |b|> %(bmin)s degrees.
        <br> Left:The number of entries per angular bin divided by the bin's solid angle. The overlaid curve is the expected
        distribution for a uniform distribution of sources with no confusion.
        <br> Right: ratio of measured to expected. 
        <br> Estimated loss: %(loss)s.
        """
        sdirs = self.df[ (np.abs(self.df.glat)>bmin) * (self.df.ts>tsmin) ]['skydir'].values
        closest = ([sorted(map(sdirs[i].difference, sdirs))[1] for i in range(len(sdirs))])
        z = np.degrees(closest)
        self.closest = z
        n = len(z)
        rho = n /(4*np.pi*np.degrees(1)**2) / (1-np.sin(np.radians(bmin)))
        f = lambda x : n*rho * np.exp( -np.pi*rho*x**2)
        
        def deficit(theta):
            n1, n1t =n-sum(z>theta), n*(1-np.exp(-np.pi * rho*theta**2))
            return n1t-n1, n, 100*(n1t-n1)/n
        loss = deficit(deficit_angle)
        self.loss ='%d/%d, or %.1f percent' % loss
        self.bmin = bmin
        print 'lost: ', self.loss
        n += loss[0]
        rho *=(1+loss[0]/n)
        fig, axx = plt.subplots(1,2, figsize=(12,5))
        plt.subplots_adjust(wspace=0.3)
        bins = np.arange(nbins+1)*dtheta
        cumarea=np.cos(np.radians(bins))*2*np.pi*np.degrees(1)**2
        dA = cumarea[:-1]-cumarea[1:]
        h = np.histogram(z,bins)[0]
        x = bins[:-1]+dtheta/2
        
        ax = axx[0]
        ax.errorbar(x, h/dA ,yerr=np.sqrt(h)/dA,  fmt='.', label='TS>%d: %d sources above |b|=%d'%(tsmin, len(sdirs),bmin))
        ax.plot(bins,f(bins), '--g')
        plt.setp(ax, yscale='log', ylim=(1,1000), 
            ylabel='Number of sources per square degree', xlabel='closest distance (deg)')
        ax.legend(prop=dict(size=10))
        ax.grid()
        
        ax=axx[1]
        ax.errorbar(x, h/dA/f(x), yerr=np.sqrt(h)/dA/f(x),  fmt='o', )
        ax.axhline(1.0, color='k')
        ax.grid()
        plt.setp(ax, xlim=(0,3), ylim=(0,1.5),  xlabel='closest distance (deg)',
            ylabel='ratio of detected to expected')
        return fig
    
    def poor_loc(self):
        """ Poorly localized sources
                %(poorly_localized_table_check)s
        """
        if len(self.poorloc)>0:
            poorly_localized_tablepath = os.path.join(self.plotfolder,'poorly_localized_table.html')
            open('poorly_localized_table.html','w').write(self.poorloc.to_html())
            print 'Wrote poorly_localized_table.html'
            open(os.path.join(poorly_localized_tablepath),'w').write(
                '<head>\n' + HTMLindex.style + '</head>\n<body>\n<h3>Poorly Localized Source Table</h3>'\
                            +  self.poorloc.to_html()+'\n</body>')
            print 'saved html doc%s' % os.path.join(poorly_localized_tablepath)
            self.poorly_localized_table_check =\
                        '<p><a href="%s"> Table of %d poorly localized (a>%.2f deg, or qual>%.1f with TS>%d) sources</a>'\
                        % ( 'poorly_localized_table.html',len(self.poorloc),self.acut,self.qualcut, self.tscut)
        else:
            self.poorly_localized_table_check ='<p>No poorly localized sources!'


    def all_plots(self):
        """Localization summary
        
        <p>Analysis of the localization information derived from the 'ellipse' field of the source records.
        """
        return self.runfigures([self.r95, self.localization,self.localization_quality,self.poor_loc,self.check_closeness,self.source_confusion])


class SourceComparison(SourceInfo):

    def setup(self, cat='gll_pscP72Y_v5r2_flags_assoc_v5r11p3.fit', #gll_psc_v06.fit', 
            catname='2FGL', **kw):
        super(SourceComparison, self).setup(**kw)
        self.catname=catname
        self.plotfolder='comparison'
        lat2 = os.path.expandvars('$FERMI/catalog/'+cat)
        assert os.path.exists(lat2), 'Did not find file %s' %cat
        ft = pyfits.open(lat2)[1].data
        name = ft.Source_Name
        ra = ft.RAJ2000
        dec= ft.DEJ2000
        id_prob = ft.ID_Probability[:,0]
        cat_skydirs = map (lambda x,y: SkyDir(float(x),float(y)), ra,dec)
        glat = [s.b() for s in cat_skydirs]
        glon = [s.l() for s in cat_skydirs]
        print 'generating closest distance to catalog "%s"' % cat
        closest = np.degrees(np.array([min(map(sdir.difference, cat_skydirs))for sdir in self.df.skydir.values]))
        self.df['closest'] = closest
        self.cat = pd.DataFrame(dict(ra=ft.RAJ2000,dec= ft.DEJ2000, ts=ft.Test_Statistic, 
                glat=glat, glon=glon, id_prob=id_prob), 
            columns = 'ra dec glat glon ts id_prob'.split(), # this to order them
            index=ft.Source_Name )
        self.cat.index.name='name'
        closest2 = np.degrees(np.array([min(map(sdir.difference, self.df.skydir.values)) for sdir in cat_skydirs]))
        self.cat['closest']= closest2
        
            
    def distance_to_cat(self, maxdist=0.5, tscuts=[10,50,500], nbins=26):
        """Associations of sources with 2FGL
        
        """
        fig,ax = plt.subplots( figsize=(4,4))
        for tscut in tscuts:
            ax.hist(self.df.closest[self.df.ts>tscut].clip(0,maxdist), np.linspace(0,maxdist,nbins), log=True,
             label='TS>%d'%tscut)
        ax.grid()
        ax.legend(prop=dict(size=10))
        plt.setp(ax, xlabel='closest distance to %s source'%self.catname)
        return fig
    
    def lost_plots(self, close_cut=0.25, minassocprob=0.8, maxts=250):
        """2FGL sources not present in new list
        Histogram of the 2FGL catalog TS and Galactic latitude for those sources more than %(close_cut).2f deg from a skymodel source. 
        The subset of sources with associations (prob>%(minassocprob)s) is shown. <br>
        Left: Distribution vs. TS.<br>
        Right: Distribution vs sine of Galactic latitude.
        """ 
        self.minassocprob=minassocprob
        self.close_cut = close_cut
        fig,axx = plt.subplots(1,2, figsize=(8,4))
        self.lost = self.cat.closest>close_cut
        print '%d sources from 2FGL further than %.2f deg: consider lost' % (sum(self.lost) , close_cut )
        self.cat.ix[self.lost].to_csv('2fgl_lost.csv')
        print '\twrite to file "2fgl_lost.csv"'
        lost_assoc = self.lost * self.cat.id_prob>0.8

        def left(ax):
            space = np.linspace(0,maxts,21)
            ax.hist(self.cat.ts[self.lost].clip(0,maxts), space, label='all (%d)'%sum(self.lost))
            ax.hist(self.cat.ts[lost_assoc].clip(0,maxts), space, label='associated(%d)' %sum(lost_assoc) )
            ax.legend(prop=dict(size=10))
            ax.grid()
            plt.setp(ax, xlabel='TS of %s source' %self.catname)

        def right(ax):
            space = np.linspace(-1,1,51)
            singlat = np.sin(np.radians(self.cat.glat))
            ax.hist(singlat[self.lost], space, label='all (%d)'%sum(self.lost))
            lost_assoc = self.lost * self.cat.id_prob>0.8
            ax.hist(singlat[lost_assoc], space, label='associated(%d)' %sum(lost_assoc) )
            ax.legend(prop=dict(size=10))
            ax.grid()
            plt.setp(ax, xlabel='sin(glat) of %s source' %self.catname, xlim=(-1,1))
            return fig
        for f, ax in zip((left,right), axx.flatten()): 
            f(ax)
        return fig
        
    def all_plots(self):
        """Results of comparison with 2FGL catalog
        """
        self.runfigures([ self.distance_to_cat, self.lost_plots])

    
class Associations(SourceInfo):

    def setup(self, **kw):
        super(Associations, self).setup(**kw)
        self.plotfolder='associations'
        probfun = lambda x: x['prob'][0] if x is not None else 0
        self.df['aprob'] = np.array([ probfun(assoc) for  assoc in self.df.associations])
        self.df['acat']  = np.array([ assoc['cat'][0] if assoc is not None else 'unid' for  assoc in self.df.associations])
        self.df10 = self.df.ix[self.df.ts>10]
        print 'associated: %d/%d' % (sum(self.df10.aprob>0.8), len(self.df10))
        
    def association_vs_ts(self):
        """ Association vs. TS
        
        Histogram of TS, showing the fraction that have associations.
        """
        fig, ax = plt.subplots( figsize=(5,5))
        bins = np.logspace(1,5,41)
        ax.hist(self.df10.ts, bins, label='all sources')
        ax.hist(self.df10.ts[self.df.aprob>0.8], bins, color='orange', label='associated')
        plt.setp(ax, xscale='log', xlabel='TS', xlim=(10,1e5))
        ax.legend(prop=dict(size=10)); ax.grid()
        return fig
    
    
    def all_plots(self):
        """ Analysis of associations
        %(atable)s
        """
        t = dict()
        for c in self.df10.acat:
            if c not in t: t[c]=1
            else: t[c]+=1
        html_rows = ['<tr><td>%s</td><td>%d</td></tr>' %x for x in t.items() ]
        self.atable = '<table><tr><th>Catalog</th><th>number</th></tr>'+\
             '\n'.join(html_rows)+'</table>'
             
        # compare with LAT pulsar catalog     
        tt = set(self.df.name[self.df.psr])
        pulsar_lat_catname = sorted(glob.glob(os.path.expandvars('$FERMI/catalog/srcid/cat/obj-pulsar-lat_*')))[-1]
        print 'opening LAT catalog file %s' %pulsar_lat_catname
        pp = pyfits.open(pulsar_lat_catname)[1].data
        lat = pd.DataFrame(pp, index=[n.strip() for n in pp.Source_Name])
        lat['ts'] = self.df[self.df.psr]['ts']
        lat['ROI_index'] = [Band(12).index(SkyDir(float(ra),float(dec))) for ra,dec in zip(lat.RAJ2000,lat.DEJ2000)]
        lat['skydir'] = map(SkyDir, lat.RAJ2000, lat.DEJ2000)
        lat['sourcedir'] = self.df.skydir[self.df.psr]
        lat['delta'] = [np.degrees(s.difference(t)) if not type(t)==float else np.nan for s,t in zip(lat.skydir,lat.sourcedir)]
        if sum(lat.delta>0.25)>0:
            print 'LAT pulsar catalog entries found more than 0.25 deg from catalog:'
            print lat[lat.delta>0.25]['ts delta'.split()]
        dc2names =set(pp.Source_Name)
        print 'sources with exp cutoff not in LAT catalog:', list(tt.difference(dc2names))
        print 'Catalog entries not found:', list(dc2names.difference(tt))
        missing = [ np.isnan(x) or x<10. for x in lat.ts]
        
        self.atable += '<h3>Compare with LAT pulsar catalog: %s</h3>' % os.path.split(pulsar_lat_catname)[-1]
        self.atable += '<p>Sources fit with exponential cutoff not in catalog %s' %list(tt.difference(dc2names))
        self.atable += '<p>%d LAT catalog entries not in the model (TS shown as NaN), or too weak.' % sum(missing)
        self.atable += lat[missing]['RAJ2000 DEJ2000 ts ROI_index'.split()].to_html()
        if sum(lat.delta>0.25)>0:
            self.atable += '<p>Pulsars located > 0.25 deg from nominal'\
                    + lat[lat.delta>0.25]['ts delta'.split()].to_html()

        
        self.runfigures([self.association_vs_ts,])
    
class Localization1K(Localization):
    """ load and analyze a special localization-only run"""
    require='localization.zip'
    
    def setup(self, zipname='localization', ecut=1000,  **kw):
        super(Localization1K, self).setup(**kw)
        try:
            f, pk = self.load_pickles(zipname)
        except:
            print 'failed to load %s, which should have zipped pickles of sources after localization attempt'%zipname
            raise
        d = dict( (x.name, x.ellipse if hasattr(x,'ellipse') else [np.nan]*7) for x in pk)
        self.ebox1K = pd.DataFrame( d ).T  
        self.ebox1K.columns = self.ebox.columns
        self.plotfolder='sources' #needed by superclass
        
    def ellipse_ratio(self):
        """ Ratio of error ellipse for E>1GeV
        Compare the error ellipse major axis for the full fit, with a fit using energies above 1 GeV 
        <br>Left: scatter plot of the ratio vs. the value, for TS>25
        <br>Right: Histogram of the ratio, showing subset with TS>100
        """
        fig, axx = plt.subplots(1,2, figsize=(12,5))
        aratio = self.ebox1K.a/self.ebox.a
        ratio_label = 'Ratio for E>1 GeV'
        cut = self.df.ts>25
        def plot1(ax):
            ax.plot( 2.5*60*self.ebox.a[cut],aratio[cut], '.')
            plt.setp(ax, xscale='log', xlim=(0.5,10), xlabel='error ellipse major axis (arc min)', 
                         yscale='linear', ylim=(0.95,1.1), ylabel=ratio_label)
            ax.axhline(1.0, color='k')
            ax.grid()
        def plot2(ax):
            xlim = (0.90, 1.25); bins=np.linspace(xlim[0],xlim[1],36)
            tscut=100
            cut2 = self.df.ts>tscut
            ax.hist(aratio[cut].clip(*xlim),bins , label='%d sources'%sum(cut.values))
            ax.hist(aratio[cut2].clip(*xlim), bins, label='TS>%d'%tscut)
            ax.legend(loc='upper right', prop=dict(size=10))
            plt.setp(ax, xlim=xlim, xlabel =ratio_label)
            ax.grid(); ax.axvline(1.0, color='k')
            
        plot1(axx[0])
        plot2(axx[1])
        
    def all_plots(self, **kw):
        self.runfigures([self.ellipse_ratio,])
    

class FluxCorr(SourceInfo):

    require='fluxcorr.zip'
    def setup(self, **kwargs):
        super(FluxCorr, self).setup(**kwargs)
        self.plotfolder='fluxcorr'
        self.source_name=kwargs.pop('source_name', 'fluxcorr')
        self.title='Source-galactic diffuse flux dependence'
        self.diffuse_name='Galactic'
        self.delta = kwargs.get('delta', 0.01) # default delta is 1%

        
        self._readdata()
        
    def _readdata(self):
        # read in the flux correlation data, in DataFrame, combine to DataFrame
        fs, ps = self.load_pickles(self.source_name)
        print 'Combining with Source info...'
        #first combine the roi DataFrames
        ndf = None 
        for x in ps:
            if len(x)>0:
                ndf = ndf.append(x) if ndf is not None else x
                
        self.df = self.df.join(ndf)
        self.emins=(100, 316, 1000) # assume 
        for emin in self.emins:
            self.df['flux_dependence_ratio_%d'%emin]= \
                (10**(ndf['par_p%d'%emin] - ndf['par_m%d'%emin])-1)/(2.*self.delta)
        
    def flux_sensitivity(self, axx=None, emin=100, **kwargs):
        """ %(diffuse_name)s diffuse flux sensitivity, emin=%(emin)s
        
        Let fs be the flux sensitivity, defined as the ratio of measured flux change to change in Galactic diffuse
        left: histogram of fs<br>
        center: scatter plot, flux sensitivity vs. TS<br>
        rignt: Skymap for fs<100.
        """
        self.emin = emin
        colorbar = kwargs.pop('colorbar', True)
        if axx is None:
            fig, axx = plt.subplots(1,3, figsize=(10,5))
        flux_ratio = self.df['flux_dependence_ratio_%d'%emin]
        def plot1(ax):
            bins = np.linspace(-30,5,36)
            ax.hist(flux_ratio.clip(bins[0],bins[-1]), bins, label='%d sources'%flux_ratio.count())
            ax.hist(flux_ratio[self.df.ts<100].clip(bins[0],bins[-1]), bins, label='TS<100')
            ax.hist(flux_ratio[self.df.ts<25].clip(bins[0],bins[-1]), bins, label='TS<25')
            plt.setp(ax, xlabel='flux sensitivity')
            ax.grid()
            ax.legend(loc='upper left', prop=dict(size=10))
        def plot2(ax):
            ax.plot( self.df['ts'], flux_ratio, '.')
            plt.setp(ax, xscale='log', xlim=(10,1e5), ylim=(-30,5) , xlabel='TS')
            ax.grid()
        def plot3(ax):
            tscut = (self.df.ts<100) & (self.df.ts>10)
            self.skyplot(-flux_ratio[tscut], ax=ax ,s=15, vmax=20, vmin=0, 
                colorbar=colorbar,cbtext='abs(flux sensitivity)')
        for ax, plot in zip(axx.flatten(), (plot1,plot2, plot3)):
            plot(ax)
   
    def flux_sensitivity_all(self, colorbar=False):
        """ %(diffuse_name)s diffuse flux sensitivity
        Rows are, from the bottom, for emin=%(emins)s MeV.<br>
        Let fs be the flux sensitivity, defined as the ratio of the measured flux change 
        to the change in %(diffuse_name)s diffuse flux.
        Columns are:
        left: histogram of fs;
        center: scatter plot, flux sensitivity vs. TS; 
        right: Skymap for TS<100, showing  where the highest sensitivity is located.
        """

        #fig, axx = plt.subplots(3,3, figsize=(9,12))
        fig, axx = self.subplot_array(hsize=(1.0, 0.75, 1.2, 0.5, 2.0,0.5 ),vsize=[1.0,0.75]*3, figsize=(10,10))
        for i, emin in enumerate(self.emins):
            self.flux_sensitivity(axx[i,:], emin, colorbar=colorbar)
        return fig
    
    def ratio_vs_stat(self):
        """Compare statistical error with systematic dependence
        The red diagonal line corresponds to the relative statistical error being equal to the flux change for a 1 percent change in the 
        %(diffuse_name)s flux. Sources below this line are at risk.
        The color represent the absolute value of the Galactic latitude.
        """
        fig, axx = plt.subplots(3,1, squeeze=False, figsize=(6,8))
        plt.subplots_adjust(left=0.2,right=0.85)
        def plot1(ax, emin=100):
            fdr = -self.df['flux_dependence_ratio_%d' % emin]
            relflux_unc = self.df['unc_z%d'%emin]*100.
            scat=ax.scatter( fdr, relflux_unc, c = np.abs(np.array(self.df.glat,float)),edgecolor='none', s=10)
            plt.setp(ax, xlabel='%s diffuse dependence ratio'%self.diffuse_name, xscale='log', xlim=(0.1,40),
                 yscale='log', ylim=(1,100));
            ax.plot([1.0, 100], (1,100.0),'r-')
            ax.text(8, 1.5, 'emin=%d'%emin)
            ax.text(25,30, '1%', color='r')
            ax.grid(True)
            return scat
        scats = map(plot1, axx.flatten(), self.emins)
        cbax = fig.add_axes((0.92, 0.3, 0.025, 0.4) )
        cb = plt.colorbar(scats[0], cbax, orientation='vertical')
        fig.text(0.1, 0.5, 'relative flux uncertainty (%)', va='center', rotation='vertical')
        cb.set_label('abs(b)')
        return fig
    
    def ts_ratio(self):
        """ Ratio of TS for fits at 100 and 1000 Mev
        Compare the ratio of the TS values for fits at 100 MeV, with that at 1 GeV, vs the 100 MeV TS
        """
        fig, ax = plt.subplots(1,1, figsize=(10,5))
        ax.plot(self.df.ts_z100, self.df.ts_z1000/self.df.ts_z100, '.')
        ax.axhline(1.0, color='k')
        plt.setp(ax, xscale='log',  xlim=(10,1e5), xlabel='TS', 
                    yscale='log', ylim=(0.1, 2), ylabel='TS_1000/TS_100')
        ax.grid()
        return fig
        
    def ts_dependence(self):
        """ TS dependence on diffuse
        """
        ts100 = self.df.ts_z100
        tsdiff = ts100 - self.df.ts_p100
        fig, axx = plt.subplots(2,1, figsize=(10,10))
        def plotit(ax, tsdiff):
            scat =ax.scatter(ts100, tsdiff/np.sqrt(ts100), s=20, edgecolor='none', c = np.abs(np.asarray(self.df.glat,float)))
            plt.setp(ax, xscale='log', xlabel='delta TS/sqrt(TS)', ylim=(-0.001, 2), xlim=(10,10000))
            ax.axhline(0, color='k')
            cbar=fig.colorbar(scat, ax=ax); cbar.set_label('|b|')
            ax.grid()
        plotit(axx[0], tsdiff)
        plotit(axx[1], self.df.ts_m100-ts100)
        return fig
    
    def all_plots(self):
        self.runfigures([self.flux_sensitivity_all, self.ratio_vs_stat, self.ts_ratio,])


class FluxCorrIso(FluxCorr):

    require = 'fluxcorriso.zip'
    def setup(self, **kw):
        super(FluxCorrIso,self).setup(source_name='fluxcorriso', **kw)
        self.title='Source-isotropic diffuse flux dependence'
        self.diffuse_name='Isotropic'
        
        self.plotfolder='fluxcorriso'

        
class GalacticSpectra(Diagnostics):

    require = 'galfits_all.zip'
    
    def diffuse_setup(self, which='gal'):
    
        self.which = which
        self.title=dict(gal='galactic', iso='isotropic')[which]
        self.plotfolder = self.title+'_spectra' 
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
                                                    
        self.plot_functions =[ [self.like_scats, self.diffuse_fits, self.bfratio_hists, self.diffuse_ratio_plot],
                    map( lambda s: self.which+'_'+s, ['likelihood_diff', 'diffuse_fits', 'bfratio', 'diffuse_ratio']),
                    ]

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
        """ Likelihood ratios for individual fits.
        These all-sky plots show, for each ROI and each energy band, the consistency of the %(title)s spectral fit 
        to a value determined for just that energy band. The distribution of the log likelihood should be approximately 
        the chi squared distribution of one degree of freedom. The lighter colors, especially red, indicate serious discrepancy.
        """
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
        """ Check the front/back consistency
        These histograms show the front/back ratio for all ROIs, and the %(latcut_name)s subset.
        """
        ax = self.multifig()
        self.bfratios = map(self.bfratio_hist, range(8), ax )
        self.multilabels('front/back fit', '', '%s Diffuse fit ratio'%self.which)
        ax[0].legend(loc='upper left',bbox_to_anchor=(-0.4,1.2));
        
        html_rows = ['<tr><td>%.0f</td><td>%.2f</td></tr>' %v[:2] for v in self.bfratios]
        from IPython.core.display import HTML
        h=HTML('<table><tr><th>Energy</th><th>front/back ratio</th></tr>'+\
             ''.join(html_rows)+'</table>')
        self.fb_ratio = h.data
        open(os.path.join(self.plotfolder,'%s_fb_ratio.html'%self.which),'w').write(h.data)
  
    def diffuse_ratio_plot(self, fignum=121):
        """ Front/back %(title)s diffuse ratio
        The front to back ratio, measured from the front/back fit ratios.
        
        <br>%(fb_ratio)s
        """
        ax = self.set_plot( None, fignum)
        vals = self.bfratios # must have generated
        x,y,yerr = [[v[i] for v in vals]  for i in range(3)]
        ax.errorbar(x, y, yerr=yerr, marker='o', ms=12,fmt='', lw=2, linestyle='None')
        plt.setp(ax, xscale='log',xlabel='Energy (MeV)', ylabel='front/back flux ratio',ylim=(0.55, 1.25))
        ax.grid(True)
        ax.axhline(1.0, color='k')
        ax.set_title('%s diffuse spectral fits'%self.which, fontsize='medium')
        
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
        """ %(title)s normalization
        For the eight lowest energy bands, the normalization factors, for both front and back. 
        
        """
        ax = self.multifig()
        self.multilabels('ratio', 'ROIs', '%s diffuse fit' %self.which)
        map(self.diffuse_fit, ax, range(8))
        return plt.gcf()
            
    def all_plots(self):
        """Set of plots to check consistency of %(title)s spectra. These result 
        from analysis of a special run that, for each ROI and each energy band, allows this diffuse component to be free.
        This is done three times: using only front, only back, and both.
        """
        self.runfigures(*self.plot_functions)
        

class IsotropicSpectra(GalacticSpectra):
    require = 'isofits.zip'

    def setup(self):
        self.diffuse_setup('iso')
        self.plot_functions[0] += [self.lowdiff_plots]
        self.plot_functions[1] += ['low_energy_difference']

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

    def lowdiff_plots(self):
        """ Isotropic bin0-bin1 differences 
        
        This is an indicator of the adequacy of the Limb contribution, sinc it only affects the lowest energy band.
        <br>Left: Histogram of the normalization difference between the two lowest energy bands.
        <br>Right: distribution of this over the plane, with the 45-degree Dec shown.
        """
        fig,ax=plt.subplots(1,2, figsize=(14,6))
        self.lowdiff_hist( ax[0])
        self.lowdiff_scat( ax[1])
        return fig
        

class SeedCheck(SourceInfo):
    require='seedcheck.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname='seedcheck'
        self.spectral_type='power law'
        self.load()
        
    def load(self):
        files, sources = self.load_pickles(self.seedname)
        sdict={}
        assoc={}
        for source in sources:
            name = source.name
            model = source.model
            pars = np.empty(4); pars.fill(np.nan)
            errs = np.empty(4); errs.fill(-2)
            free = np.zeros(4, bool)
            n = model.len()
            pars[:n] = model.parameters
            free[:n] = model.free
            try:
                diag = np.diag(model.get_cov_matrix())
                errs[:n] = [np.sqrt(x) if x>0 else -1 for x in diag[:n]]
                badfit = np.any(errs[model.free]<=0)
            except Exception, msg:
                print 'fail errors for %s:%s' % (name, msg)
                badfit = True
            has_adict = hasattr(source,'adict') and source.adict is not None
            sdict[name] = dict(
                ra =source.skydir.ra(), dec=source.skydir.dec(),
                ts=source.ts,
                delta_ts=source.ellipse[5] if hasattr(source, 'ellipse') else np.nan,
                r95 = 2.6*source.ellipse[2] if hasattr(source, 'ellipse') else np.nan,
                glat=source.skydir.b(), glon=source.skydir.l(),
                eflux=pars[0] * model.e0**2 *1e6,
                eflux_unc=errs[0] * model.e0**2 *1e6 if errs[0]>0 else np.nan,
                pindex = pars[1],
                pindex_unc = errs[1] if errs[1]>0 else np.nan,
                par2 = pars[2],
                par2_unc = errs[2] if errs[2]>0 else np.nan,
                e0 = model.e0,
                aprob = source.adict['prob'][0] if has_adict else 0,
                index = source.index,
                #gflux  = model.i_flux(), ## photon flux
                )
            assoc[name] = dict(
                acat = source.adict['cat'][0] if has_adict else None,
                aname= source.adict['name'][0] if has_adict else None,
                adelta_ts = source.adict['deltats'][0] if has_adict else None,
                aprob = source.adict['prob'][0] if has_adict else 0.,
                adict = source.adict if has_adict else None,
                )
        self.df = pd.DataFrame(pd.DataFrame(sdict).transpose(), 
            columns='ra dec glat glon ts delta_ts eflux eflux_unc pindex pindex_unc par2 par2_unc e0 r95 aprob index'.split() 
            )
        self.df.index.name='name'
        self.assoc = pd.DataFrame(assoc).transpose()
        self.assoc.index.name = 'name'
        # analyze associations, make summary
        acat=list(self.assoc.ix[self.assoc.aprob>0.8]['acat'].values)
        sa = list(set(acat))
        t = np.zeros(len(sa),int)
        for x in acat:
            t[sa.index(x)]+=1
        self.assoc_sum = zip(sa, t)
        #self.psr = [  s is not None and ('pulsar' in s or 'msp' in s ) for s in t] * associated
        #self.agn = [  s is not None and ('agn' in s or 'crates' in s or 'bzcat' in s or 'bllac' in s or 'qso' in s or 'cgrabs' in s) for s in t] * associated
    
    def select_candidates(self, tsmin):
        t=self.df[(self.df.ts>tsmin)*(self.df.r95<0.6)]['ra dec ts pindex r95 aprob'.split()]#.sort_index(by='ts')
        t['acat']=self.assoc.acat
        t['aname']=self.assoc.aname
        return t
    
    def seed_cumulative_ts(self, cut=None, label='all seeds'):
        """ Cumulative TS distribution for seeds 
        """
        v = self.df.ts
        if cut is not None: v=v[cut]
        fig = self.cumulative_ts(v, check_localized=False, label=label)
        ax = plt.gca()
        plt.setp(ax, ylim=(9,1000), xlim=(9,100))
        leg =ax.legend()
        pbox = leg.get_patches()[0]; pbox._height=0; pbox._y=5
        return fig
        
    def unassoc_seed_cumulative_ts(self):
        """ Cumulative TS distribution for seed sources that are not associated
        """
        return self.seed_cumulative_ts(cut=self.assoc.aprob<0.8, label='unassociated')
        
    def histo(self, ax, v, bins):
        ax.hist(v, bins)
        ax.hist(v[self.df.ts>10], bins, label='TS>10')
        ax.hist(v[self.df.ts>25], bins, label='TS>25')
        ax.legend(prop=dict(size=10))
        ax.grid()
    
    def localization(self):
        """ Localization results
        <br>Left: r95; right; delta TS
        """
        fig, ax = plt.subplots(1,2, figsize=(12,5))
        def r95(ax):
            v = 60.*self.df.r95; bins = np.linspace(0, 25,26)
            self.histo(ax, v[~pd.isnull(v)], bins)
            plt.setp(ax, xlabel='r95 (arcmin)')
        def delta_ts(ax):
            v = np.sqrt(list(self.df.delta_ts)); bins = np.linspace(0,10,26)
            self.histo(ax, v, bins)
            plt.setp(ax, xlabel='sqrt(delta_ts)', xlim=(0,10))

        for f, a in zip((r95, delta_ts), ax.flatten()):
            f(a)
        return fig

    
    def spectral_parameters(self, ax=None):
        """ Spectral fit parameters
        Flux vs. spectral index for %(spectral_type)s fit
        <br>histograms of sin(glat) and sqrt(delta_ts) for all, TS>10, and TS>25
        """
        fig, ax = plt.subplots(2,2, figsize=(12,12))
        good = self.df.ts>10
        super = self.df.ts>25
        def flux_index(ax):
            for cut, c,label in zip((good, super), ('.b', 'or'), ('TS>10', 'TS>25')):
                ax.plot(self.df.eflux[cut], self.df.pindex[cut], c, label=label)
            ax.grid()
            ax.legend(prop=dict(size=10))
            plt.setp(ax, ylim=(0.5,3.0), xlim=(0.1,10), xscale='log', ylabel='spectral index', xlabel='enegy flux (eV)')
        def singlat(ax):
            v = np.sin(np.radians(list(self.df.glat))); bins=np.linspace(-1,1,26)
            self.histo(ax, v, bins)
            plt.setp(ax, xlabel='sin(glat)')
        def skyplot(ax):
            glon = self.df.glon
            glon[glon>180]-=360
            ax.plot(glon, np.sin(np.radians(list(self.df.glat))), 'o')
            plt.setp(ax, xlim=(180,-180), xlabel='glon', ylabel='sin(glat)')
            #self.skyplot(self, self.df.ts, ax=ax)
        def index_vs_cutoff(ax):
            cutoff = self.df.par2
            for tsmin, marker in zip((10,25), ('.b', 'or')):
                cut = self.df.ts>tsmin
                ax.plot(cutoff[cut], self.df.pindex[cut],  marker, label='TS>%d'%tsmin)
            plt.setp(ax, ylabel='spectral index', xlabel='cutoff', ylim=(0.5,3.0), xlim=(0, 3000))
            ax.grid(); ax.legend(prop=dict(size=10))
        for f, a in zip((flux_index, index_vs_cutoff, singlat, skyplot), ax.flatten()):
            f(a)
            
        return fig

    def locations(self):
        """ Positions
        """
        return self.skyplot(self.df.ts)
        
    def all_plots(self):
        """ Results of analysis of seeds
        %(info)s
        """
        #t=self.df[(self.df.ts>6)*(-self.df.r95.isnull())]['ra dec ts pindex r95 aprob index'.split()]#.sort_index(by='ts')
        #t['acat']=self.assoc.acat
        #t['aname']=self.assoc.aname
        tsmin = 6
        t = self.select_candidates(tsmin)
        good_seeds = 'good_seeds.csv'
        t.to_csv(good_seeds)
        print 'wrote list that succeeded to %s' % good_seeds
        self.info = ' ' #self.df.describe().to_html()
        self.info += '<h3> Selected TS>%.0f and localized</h3>' %tsmin
        self.info += t.to_html()
        self.info += '<h3>Association summary:</h3>' #\n<pre>%s\n</pre>' %self.assoc_sum
        self.info += '<table border="1"><thead><tr><th>Catalog</th><th>Sources</th></tr></thead>\n<tbody>'
        for (c,n) in  self.assoc_sum:
            self.info += '<tr><td>%s</td><td>%5d</td></tr>' % (c,n)
        self.info += '</tbody></table>\n'
        self.runfigures([self.seed_cumulative_ts, self.spectral_parameters, self.localization,])
        
class PulsarSeedCheck(SeedCheck):
    require='pseedcheck.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname= 'pseedcheck'
        self.spectral_type = 'exponential cutoff'
        self.load()

    def all_plots(self):
        """ Results of analysis of pulsar seeds
        %(info)s
        """
        self.info = self.df.describe().to_html()
        self.runfigures([self.seed_cumulative_ts, self.unassoc_seed_cumulative_ts, self.spectral_parameters, self.localization],
                ('pulsar_cumulative_ts', 'pulsar_unassoc_cumulative_ts', 'pulsar_spectral_pars', 'pulsar_localization'))

class PGWSeedCheck(SeedCheck):
    require='seedcheck_PGW.zip'
    def setup(self, **kw):
        self.plotfolder = self.seedname= 'seedcheck_PGW'
        self.spectral_type = 'power law'
        self.load()

class HPtables(Diagnostics):
    """ process Healpix tables, inclucing TS residual map files generated by table, perhaps generating list of new seeds """
    require = 'ts_table' ## fix.
    def setup(self, **kw):
        fnames = glob.glob('hptables_ts*.fits')
        assert len(fnames)==1, 'expect one hptable*.fits file'
        self.fname=fnames[0]
        self.tables = pd.DataFrame(pyfits.open(self.fname)[1].data)
        self.plotfolder = 'hptables'
        self.tsname='ts'
        self.seedfile, self.seedroot, self.title= 'seeds.txt', 'SEED' ,'power-law'
        self.make_seeds(refresh=kw.pop('refresh', False))

     
    def make_seeds(self, refresh=False):
        """ may have to run the clustering application """
        if not os.path.exists(self.seedfile) or os.path.getmtime(self.seedfile)<os.path.getmtime(self.fname) or refresh:
            print 'reconstructing seeds: %s --> %s' % (self.fname, self.seedfile)
            cmd = 'python -m uw.like2.pipeline.check_ts %s %s --seedroot=%s --tsfield=%s' % (self.fname, self.seedfile, self.seedroot, self.tsname)
            print '-->',cmd
            os.system(cmd)
        self.seeds = pd.read_table(self.seedfile)
        self.n_seeds = len(self.seeds)
        print 'read in %d seeds from %s' % (self.n_seeds, self.seedfile)
    
    def kde_map(self, vmin=1e5, vmax=1e8, pixelsize=0.25):
        """Photon Density map
        All data, smoothed with a kernel density estimator using the PSF.
        """
        hpts = healpix_map.HParray('kde', self.tables.kde)
        hpts.plot(ait_kw=dict(pixelsize=pixelsize), norm=LogNorm(vmin, vmax))
        return plt.gcf()
     
    def ts_map(self, vmin=10, vmax=25, pixelsize=0.25):
        """ TS residual map 
        DIstribution of TS values for %(title)s residual TS study.
        """
        hpts = healpix_map.HParray(self.tsname, self.tables[self.tsname])
        hpts.plot(ait_kw=dict(pixelsize=pixelsize), vmin=vmin, vmax=vmax)
        return plt.gcf()
        
    def seed_plots(self):
        """ Seed plots
        
        Results of cluster analysis of the residual TS distribution. Analysis of %(n_seeds)d seeds from file %(seedfile)s. 
        <br>Left: size of cluster, in 0.15 degree pixels
        <br>Center: maximum TS in the cluster
        <br>Right: distribution in sin(|b|), showing cut.
        """
        z = self.seeds
        fig,axx= plt.subplots(1,3, figsize=(12,4))
        ax = axx[1]
        ax.hist(z.ts.clip(0,50), np.linspace(0,50));
        plt.setp(ax, xlabel='TS');
        ax.grid()
        ax =axx[0]
        ax.hist(z.size.clip(0,20), np.linspace(0,20,21)); 
        plt.setp(ax, xlabel='cluster size'); ax.grid()
        ax=axx[2]
        sinb = np.sin(np.radians(z.b))
        ax.hist(sinb, np.linspace(-1,1,26))
        ax.axvline(0, color='k')
        plt.setp(ax, xlabel='sin(b)')
        ax.grid()
        return fig
    
    def all_plots(self):
        """ Results of search for new sources, using a power-law spectrum"""
        self.runfigures([ self.seed_plots, self.ts_map, self.kde_map,])

class PTStable(HPtables):
    require = 'hptables_pts.fits'
    def setup(self, **kw):
        fnames = glob.glob('hptables_pts*.fits')
        assert len(fnames)==1, 'expect one hptable_pts*.fits file'
        self.fname=fnames[0]
        self.tables = pd.DataFrame(pyfits.open(self.fname)[1].data)

        self.tsname = 'pts'
        self.plotfolder = 'pulsar_ts'
        self.seedfile, self.seedroot, self.title = 'pseeds.txt', 'PSEED' ,'pulsar'
        self.make_seeds(kw.pop('refresh', False))
        
    def all_plots(self):
        """ Results of Pulsar seed search """
        self.runfigures([ self.seed_plots, self.ts_map,])


class Components(ROIinfo):
    def setup(self):
        self.plotfolder='components'
        self.title = 'All components'
        self.plots_kw={}
        self.components = [x() for x in (Galactic, Isotropic, Limb, SunMoon, SourceTotal)]
        self.funcs = np.hstack([x.funcs for x in self.components])
        self.fnames= np.hstack([x.fnames for x in self.components])

class Data(Diagnostics):
    require='config.txt'
    """ look at binned data """
    def setup(self):
        self.plotfolder = 'data'
        config = eval(open('config.txt').read())
        datadict = config['datadict']
        self.dataset = dataset.DataSet(datadict['dataname'], interval=datadict.get('interval',None),
            irf=config['irf'],)
        t =[]
        for band in self.dataset.dmap:
            t.append(dict(emin=band.emin(), emax=band.emax(), photons=band.photons(), 
                ec=band.event_class(), nside=band.nside(), pixels=band.size()))
        self.df = pd.DataFrame(t)
        sm= np.sum(self.df.photons[self.df.emin>=100])
        self.total_counts ='{:,d}'.format(int(sm))

            
    def plot_spectrum(self):
        """Spectrum of all data
        Total above 100 MeV: %(total100)s events.
        """
        fig, ax = plt.subplots(figsize=(4,4))
        photons = self.df.photons.values
        combine = photons[0::2]+photons[1::2]
        elow = self.df.emin.values[0::2]
        ehigh = self.df.emax.values[0::2]
        ax.plot(ehigh, combine, ls='steps', lw=2)
        ax.plot([elow[1],elow[1]], (100, combine[1]), '-b', lw=2) # vertical line for first bin
        plt.setp(ax, xscale='log', xlabel='Energy (MeV)', xlim=(10,1e6), yscale='log', ylabel='events/bin')
        ax.grid()
        return fig  

    def all_plots(self):
        """ plots involving the data """
        self.runfigures( [self.plot_spectrum,])
    
class HTMLindex():
    """ Manage the web browser pages
    """
    style="""
<style type="text/css">
body {	font-family:verdana,arial,sans-serif;
	font-size:10pt;
	margin:10px;
	background-color:white;
	}
p   { font-size:10pt; margin-left:25pt; }
pre { font-size:10pt; margin-left:25pt; 
    border-style:solid;
    border-width:thin;}
h5 {margin-left:25pt;}
table { margin-left:25pt; margin-top:15pt; font-size:8pt;
    border-style: solid; border-width: 1px;  border-collapse: collapse; }
table td { padding: 3px; }
a:link { text-decoration: none ; color:green}
a:hover { background-color:yellow; }
</style>"""

    menu_header="""<html> <head> <title>Plot index for model %(model)s </title> %(style)s 
    <script> function load(){ parent.content.location.href = '%(model_summary)s';} </script>
    </head>
<body onload="load()">
<h2><a href="%(model_summary)s", target="content">%(model)s</a></h2>"""

    top_nav="""<html> <head> <title>Top Nav</title> %(style)s 
    <script> function load(){ parent.menu.location.href = '%(last_model)s';} </script>
    </head>
<body onload="load()">
<h3>Skymodels</h3>""" 

    def __init__(self, folder='plots/*'):
        self.style = HTMLindex.style
    
        w= glob.glob(folder)
        if len(w)==0: 
            print 'Did not find any plot folders under %s' % folder
        z = dict( zip(w, [glob.glob(a+'/*.htm*') for a in w] ) )
        self.model = '/'.join(os.getcwd().split('/')[-2:])
        self.model_summary='plots/config.html' ## TODO: make a special page
        s= HTMLindex.menu_header % self.__dict__
        
        def parse_item(x):
            head, tail =os.path.split(x)
            name = os.path.splitext(tail)[0]
            n = name.find('_uw')
            return '<a href="%s" target="content">%s</a><br>' % (x,name[:n])

        for k in sorted(z.keys()):
            v = z[k]
            if len(v)==0: continue
            index = '%s/index.html'%k
            if index in v:
                v.remove(index)
                s += '\n<h4><a href="%s" target="content">%s</a></h4>'% (index, k.split('/')[-1])
            else:
                s += '\n<h4>%s</h4>'% k.split('/')[-1]
            s += '\n\t<p>' + '\n\t'.join(map(parse_item, v)) 
        self.ul = s + '</p>\n</body>'
        self.make_config_link()
        
    def _repr_html_(self):    return self.ul
    
    def make_config_link(self):
        html = '<head>%s</head><body><h3>%s - configuration and analysis history files</h3>' %(self.style,self.model)
        for filename in ('config.txt', 'dataset.txt', 'converge.txt', 'summary_log.txt'):
            if not os.path.exists(filename): continue
            html += '<h4>%s</h4>\n<pre>%s</pre>' % (filename, open(filename).read())
        html += '\n</body>'
        open('plots/config.html', 'w').write(html)
        print 'wrote plots/config.html'
    def create_menu(self, filename='plot_index.html'):
        ###summary = open(filename, 'w')
        open(filename, 'w').write(self.ul)
        print 'wrote menu %s' % os.path.join(os.getcwd(),filename)

    def update_top(self, filename='../../plot_browser/top_nav.html'):
        def parse_path(x): 
            'return relative path, model name'
            t = x.split('/')
            return  '/'.join(t[1:]) , '/'.join(t[2:4])
        def parse_model(x):
            return '<a href="%s" target="menu"> %s </a>' %(parse_path(x) )
        models = sorted(glob.glob('../../*/*/plot_index.html'))
        self.last_model = parse_path(models[-1])[0]
        s = HTMLindex.top_nav % self.__dict__
        s += '\n | '.join(map(parse_model, models))
        s += '\n</body></html>\n'
        open(filename, 'w').write(s)
        print 'wrote top menu %s' % os.path.join(os.getcwd(),filename)
    @staticmethod
    def head(title=''):
        return '<head><title>%s</title>\n'+HTMLindex.style+'</head>\n'

opts = dict(
        counts=  (CountPlots,),
        sources= (SourceInfo, Localization, SourceTotal,),
        diffuse= (Galactic, Isotropic, Limb, SunMoon),
        isotropic=(Isotropic,),
        galactic=(Galactic,),
        limb=    (Limb,),
        limb_refit=(LimbRefit,),
        sunmoon= (SunMoon,),
        sunmoon_refit = (SunMoonRefit,),
        isospect =  (IsotropicSpectra,),
        galspect =  (GalacticSpectra,),
        fb=      (FrontBackSedPlots,),
        fluxcorr=(FluxCorr,),
        fluxcorriso=(FluxCorrIso,),
        loc =    (Localization,),
        loc1K =  (Localization1K,),
        exp =    (Exposure,),
        exposure=(Exposure,),
        hptables = (HPtables,),
        tables = (HPtables,),
        sourcetotal=(SourceTotal,),
        seedcheck=(SeedCheck,),
        pseedcheck=(PulsarSeedCheck,),
        pts=     (PTStable,),
        data =   (Data,),
        comparison=(SourceComparison,),
        association=(Associations,),
        ) 
        
def main(args, update_top=False , raise_exception=False):
    np.seterr(invalid='warn', divide='warn')
    success=True
    if type(args)==types.StringType: args = args.split()
    for arg in args:
        if arg=='all':
            cs = set(np.hstack(opts.values()))
            for cls in cs:
                if os.path.exists(cls.require):
                    print 'running %s' % cls.__name__
                    try:
                        cls('.').all_plots()
                        plt.close('all')
                    except Exception, msg:
                        print '=====failed====\n %s\n=============='% msg
                else:
                    print 'skipped %s, missing %s' % (cls.__name__, cls.require)
            break
        if arg=='menu': 
            update_top=True
            continue

        if arg not in opts.keys():
            print 'found %s; expect one of %s' % (arg, opts.keys())
            continue
            success = False
        try:
            for cls in opts[arg]:
                cls('.').all_plots()
                plt.close('all')
        except FloatingPointError, msg:
            print 'Floating point error running %s: "%s"' % (arg, msg)
            print 'seterr:', np.seterr()
            success=False
        except Exception, msg:
            print 'Exception running %s: "%s"' % (arg, msg)
            if raise_exception: raise
            success = False
    if success: 
        HTMLindex().create_menu()
        if update_top: HTMLindex().update_top()
        
    return success    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='run a diagnostic output job; must be in skymodel folder')
    parser.add_argument('args', nargs='+', help='processsor identifier: must be one of %s' %opts.keys())
    parser.add_argument('--update_top', action='store_true', help='Update the top level Web  menu')
    args = parser.parse_args()
    if not main(args.args, update_top=args.update_top):
        raise Exception
    
