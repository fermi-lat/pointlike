"""
Isotropic fitting analysis

"""
import os, glob, pickle, zipfile
import numpy as np
import pylab as plt
import pandas as pd
from  matplotlib import (patches, gridspec)
from skymaps import SkyDir, Band 
from . import (diffuse_fits, analysis_base,)
from .. import (tools, configuration,  diffuse)
from ..pipeline import stream


class IsotropicModel(object):
    """Characterize the configuration for an isotroic diffuse model, 
       which includes effects of background from the Earth limb.
    """
    def __init__(self, config=None, sdbins=40):
        """config: a Configuration object, containing the IRF and isotropic files
               if None, load it from the configuration info in the local folder
        """ 
        if config is None:
            config = configuration.Configuration('.', quiet=True)
        self.config=config
        self.sdbins=sdbins
        self.isoinfo = self.config['diffuse']['isotrop']
        self.isomodel_file = self.isoinfo['filename']
        if self.isoinfo.get('key', None)=='iso':
            zipfilename='pickle'
            #z = zipfile.ZipFile(zipfilename)
            files = sorted( glob.glob('{}/*.pickle'.format(zipfilename)) ) 
            print 'Loading %d *.pickle files in folder %s' % (len(files), zipfilename)
            fnorm=[]; bnorm=[]
            for f in files:
                pk = pickle.load(open(f))
                t =pk['diffuse_normalization']['iso']
                fnorm.append(t['front'])
                bnorm.append(t['back'])

            self.isonorm= [pd.DataFrame(fnorm),pd.DataFrame(bnorm)]
            self.isonorm_file= '(local)'
        else:
            self.isonorm_file = self.isoinfo['correction']
            self.isofiles = [self.isoinfo['correction'].replace('*', etn)
                 for etn in self.config.event_type_names]
            try:    
                [self.isofiles_front, self.isofiles_back] = \
                [os.path.expandvars('$FERMI/diffuse/')+ f for f in self.isofiles ]
            except:
                raise Exception('no isotropic normalization files found in {}'.format(self.isoinfo))
            self.isonorm = [pd.read_csv(f, index_col=0) 
                for f in (self.isofiles_front,self.isofiles_back)]
            for i in range(2): #replace '0' with '0' etc.
                self.isonorm[i].columns=pd.RangeIndex(8)

        # determine the simplified model from the correction files 
        sdirs = map(lambda i: Band(12).dir(i), range(1728)) # list of SkyDirs for nside=12
        dec = np.array(map(lambda sd: sd.dec(), sdirs), float) # corresponding Dec values
        self.sindec = np.sin(np.radians(dec)) #list of sin(Dec)
        self.sdindex = np.digitize(self.sindec, np.linspace(-1,1,sdbins+1))-1 #dec bins for low energy back
        self.glat = np.array(map(lambda sd: sd.b(), sdirs),float)  # corresponding glat

        def decbin(norms, binned):
            """ return mean of norms, binned in nbins in sin(dec)
            """
            if not binned:
                return np.mean(norms)
            df=pd.DataFrame([self.sdindex,norms], index=['sdindex','norm']).T
            g = df.groupby('sdindex'); 
            m = g.mean(); 
            return m.norm.values

        self.back_model = [decbin(self.isonorm[1][i],True) for i in range(2)]\
            +             [decbin(self.isonorm[1][i],False) for i in range(2,8)]
        self.front_model =[decbin(self.isonorm[0][i],False) for i in range(8)]

    def __repr__(self):
        s = 'Isotropic Model\n\t     spectrum: {}\n\tnormalization: {}'.format(
            self.isomodel_file, self.isonorm_file)
        s += '\n  Normalization means per energy band'
        s += '\n\tFront: {}'.format(self.isonorm[0].values.mean(axis=0).round(2))
        s += '\n\tBack:  {}'.format(self.isonorm[1].values.mean(axis=0).round(2))
        return s

    def isonorms(self,  roi_index):
        """Return a dict with keys 'front', 'back', containing the normalization values for the given roi
        """
        r = dict()
        r['front'] = self.front_model
        ib = self.sdindex[roi_index]
        r['back']= [self.back_model[j][ib] for j in range(2)]+self.back_model[2:]
        return r

    def update_model_with_fit(self, isofits, roisel):
        """Update the current model with a set of normalization fits

            isofits: float array with shape (1728,8,2) with ROI fits
            roisel:  bool array with shape (1728) for selection
        """
        model = [self.front_model, self.back_model]
        print 'Adjusting model\n ie ib   model fit'
        for ib in range(2):
            for ie in range(8):
                fit = isofits[:,ie,ib][roisel].mean(axis=0)
                m = model[ib][ie]
                print '{:2d}{:2d}{:8.3f}{:8.3f}'.format(ie, ib, np.mean(m), fit)
                #print m, '->', m * fit
                model[ib][ie] = m * fit      
        return model
        files = sorted( glob.glob('{}/*.pickle'.format(zipfilename)) )

    def set_new_model(self, newiso, zipfilename='pickle/'):
        """ version with new iso for all"""
        files = sorted( glob.glob('{}/*.pickle'.format(zipfilename)) )
        assert len(files)==1728
        print 'updating with\n {}'.format(pd.DataFrame(newiso))
        def set_isonorm(f):
            pk = pickle.load(open(f))
            pk['diffuse_normalization']['iso']=newiso
            pickle.dump(pk, open(f, 'w') )
        map(set_isonorm, files)
        print 'Updated the normalizations: you must zip the files'


    def save_model(self, path='pickle/'):
        # Note, can save normalizations from current model into another folder's pickle
        ff = sorted(glob.glob('{}*.pickle'.format(path)))
        assert len(ff)==1728, 'Did not find 1728 pickle files in path {}'.format(path)
        for i,f in enumerate(ff):
            pk =pickle.load(open(f))
            pk['diffuse_normalization']['iso']=self.isonorms(i)
            pickle.dump(pk, open(f, 'w'))
        print 'Updated the normalizations: you must zip the files'

    def plot(self, ax=None):
        """Isotropic model summary
                
        <br>While the normalization factors are 1728x8 arrays of corrections applied to each ROI and band, only the Back 
        varies for the first two energy bins.
        <ul>
        <li>
        <b>Left</b>: The isotropic spectral model, based on %(isomodel_file)s</li>
        <li>
        <b>Right</b>: The normalization factor applied to all ROIs</li>
        

        </ul>
         """

        if ax is None:
            fig, axx = plt.subplots(1,2, figsize=(10,5))
            plt.subplots_adjust(wspace=0.3)
        else: fig=ax.figure

        def left(ax, emin=90,emax=1e4):
            ecut=(100, 1e4)
            fname = os.path.expandvars('$FERMI/diffuse/')+self.isoinfo['filename']
            idfiles = [fname.replace('**', fb.upper()).replace('*', fb) for fb in ['front','back']]
            nf,nb = map(np.loadtxt, idfiles)
            energies = nf[:,0]; front,back = nf[:,1],nb[:,1]
            ecut= (energies>=emin)& (energies<=emax)
            ax.plot(energies[ecut], (front*energies**2)[ecut], 'xg', label='front')

            ax.plot(energies[ecut], (back*energies**2)[ecut], '+r', label='back')
            plt.setp(ax, xlabel='Energy', ylabel=r'$E^2\ \times \ Flux$', xscale='log')
            ax.set_title('isotropic diffuse spectra')
            ax.grid(True); ax.legend()

        def center(ax):
            for f, name,start in [(self.front_model, 'Front',0), (self.back_model, 'Back',2)]:
                means = [np.mean(f[i]) for i in range(start,8)]
                ax.plot(range(start,8), means, '--o', label=name)
            ax.set_title('Normalization factor vs Energy Bin')
            ax.set(xlabel='Energy Bin',ylabel='Normalization Factor',)

        def right(ax):
            t = np.linspace(-1,1,self.sdbins+1)
            x = 0.5*(t[1:]+ t[:-1]) #bincenters
            for i in range(2):
                ax.plot(x, self.back_model[i], '.', label='Energy Bin {}'.format(i));
            ax.set(xlabel='sin(Dec)', ylabel='Normalization factor',  title='Back normalization vs. Dec.')

        left(axx[0]); center(axx[1])#; right(axx[2]); 
        for ax in axx[1:]:
            ax.grid(alpha=0.5);
            ax.axhline(1.0, color='k', ls='--')
            ax.legend()
        return fig

class Isotropic(diffuse_fits.DiffuseFits):
    """<b>Isotropic fit analysis</b>

    <p>This is a summary of a special check of the isotropic normalization. For each ROI,
    and then for each energy band, the run measures the best normalization of both
    front and back. 
    <p>Note that as of uw8500, Back is restricted to >300 MeV, so we assume it is truly isotropic, bypassing 
    lots of logic that allows every ROI to be different.

    """
    def setup(self,**kwargs):
        super(Isotropic, self).setup()
        self.plotfolder = 'isotropic'
        assert os.path.exists('isotropic_fit'), 'No isotropic analysis to summarize'
        self.fit_profile = {}
        self.isomodel = IsotropicModel(self.config)


    def summary(self, bcut=10):
        """Summary
        <pre>%(logstream)s</pre>
        """
        model = '/'.join(os.getcwd().split('/')[-2:])
        streamdf= pd.DataFrame(stream.StreamInfo(model)).T
        self.startlog()
        print self.isomodel 

        # process isotropic fit info
        files = sorted(glob.glob('isotropic_fit/*.pickle'))
        if len(files)>0:
            if len(files)<1728:
                msg= "found {} files, expected 1728".format(len(files))
                print msg
                raise Exception(msg)
            try:
                self.isofits = np.array([pickle.load(open(f))['val'] for f in files]);
            except: #old format
                self.isofits = np.array([pickle.load(open(f)) for f in files]);
                
            try:
                snum=streamdf.query('stage=="fitisotropic"').index[-1]
            except:
                print streamdf
                raise
            print 'loaded isotropic fits, generated by stream {} at {}'.format(snum,streamdf.loc[snum].date )
        
        # roi selection
        dec = np.array(self.df.dec,float)
        ra = np.array(self.df.ra, float)
        self.df['sindec'] = np.sin(np.radians(dec))
        glat = np.array(self.df.glat,float)
        glon = np.array(self.df.glon,float).copy()
        glon[glon>180]-=360
        self.df['highlat'] = highlat= np.abs(glat)>bcut

        # check fit data, collect means for |b|.bcut    
        print 'High latitude selection (|b|>{}): {}/{} ROIs'.format(bcut,sum(highlat), len(dec))
        self.fit_means = np.ones((8,2))

        def fit_sum(ib, start=0):
            print '\n{} normalization fit summary for |b|>{}\n\t{:4s}{:>8s}{:>8s}'.format(
                ['front','back'][ib], bcut,'band','mean','std')
            for i in range(start, self.isofits.shape[1]):
                t = self.isofits[:,i,ib][self.df.highlat]
                self.fit_means[i,ib]= np.mean(t)
                print '\t{:4d}{:8.3f}{:8.3f}'.format(i, np.mean(t), np.std(t))
        fit_sum( 0 )
        fit_sum( 1, 2 )  

        self.logstream= self.stoplog()


    def set_new_iso(self, update=False):
        t = np.array([iso.values.mean(axis=0) for iso in self.isomodel.isonorm]); 
        z = self.fit_means.T
        y = z*t
        newiso= dict(front=y[0], back=y[1])
        print 'new iso: {}'.format(newiso)
        if update:
            self.isomodel.set_new_model(newiso)
        else:
            print 'Set update to True to modify sky model'


    def multiplot(self, q, suptit, selection, ylim=(0.5,1.5),nbins=20, start=0):

        def dec_profile(t, ax, title, selection ):
            z=tools.Profile(self.df.sindec[selection],t[selection], np.linspace(-1,1,nbins+1))
            z.plot(ax=ax, legend=False);
            ax.set_xlim(-1,1); 
            ax.set_xlabel('sin(Dec)')
            ax.set_title('{}'.format(title), fontsize=12)
            ax.grid(alpha=0.5);
            ax.text(0.05,0.9, '<y>={:.3f}'.format(z.mean), transform=ax.transAxes)
            return z

        zz=[]
        fig,axx = plt.subplots(2,4, figsize=(12,6),sharey=True, sharex=True)
        for n,ax in enumerate(axx.flatten()):
            if n<start:
                ax.set_visible(False)
                continue
            zz.append(dec_profile(q[:,n],ax=ax,
                 title='{:.0f}-{:.0f} MeV'.format(10**(2+0.25*n),10**(2.25+0.25*n)), selection=selection))
            ax.axhline(1.0, color='grey', ls='--')
        axx[0,0].set_ylim(ylim)    
        fig.suptitle(suptit)
        self.fit_profile[suptit]=zz
        return fig


    @tools.decorate_with(IsotropicModel.plot)
    def model_plots(self):
        self.isomodel_file = self.isomodel.isomodel_file
        return self.isomodel.plot()

    def front_fit_map(self):
        """Front fit map
        """
        return self.correction_plots(self.isofits[:,:,0], title='front fit values', vmin=0.5, vmax=1.5)

    def back_fit_map(self):
        """Back fit map
        """
        return self.correction_plots(self.isofits[:,:,1],  title='back fit values', vmin=0.5, vmax=1.5, start=2)
    def front_fit_summary(self):
        """Front fit vs. Dec"""
        return self.multiplot(self.isofits[:,:,0], 'front', self.df.highlat, ylim=(0.8,1.2))

    def back_fit_summary(self):
        """Back fit vs. Dec"""
        return self.multiplot(self.isofits[:,:,1], 'back', self.df.highlat, ylim=(0.8,1.2), start=2)        

    def all_plots(self, **kw):
        self.runfigures([
            self.summary,
            self.model_plots,
            self.front_fit_map,
            self.back_fit_map,
            self.front_fit_summary,
            self.back_fit_summary,
        ] )

    def update_model(self):
        """Update the model using the results of the fit
        (intended for hand use: must call summary() first)
        """
        self.isomodel.update_model_with_fit(self.isofits,self.df.highlat )
        self.isomodel.save_model()

