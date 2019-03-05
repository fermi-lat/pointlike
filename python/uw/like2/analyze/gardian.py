"""gardian analysis code

This extracts files associated with a GarDian run to allow various plots
"""
import os, sys, glob, healpy
import matplotlib.pylab as plt
import numpy as np
import pandas as pd

from uw.utilities import healpix_map as hpm
from uw.like2 import diffuse
from skymaps import SkyDir, Band

class GComponent(object):
    """ Study components of Seth's model
    """
    def __init__(self, component, path, energy=131., nside=None):
        """
        component : string
            C0_n, HI_n, IC_n, where n=0..9
            DNMn, DNMp, 
            EGBfree,
            patch
            Moon, SolarDisk, SolarIC, 
            Sources, UNRESOLVED 
        """
        erange='E0_PSF3 E1_PSF23 E2_PSF123 E3_PSF0123'.split()
        if energy>=1000:
            ie=3
        elif energy>=300:
            ie=2
        else:
            ie=1
        def get_file(ie, pat):
            if pat !='': pat = '_'+pat
            t = glob.glob(path+erange[ie]+pat+'*.gz')
            assert len(t)>0, 'Failed to find pattern {}'.format( erange[ie]+pat+'*.gz')
            return t

        self.energy=energy
        self.cname=component
        files = get_file(ie, component)
        print 'Combining {:2d} file(s) for component {}'.format (len(files),  component)
        self.dmaps = [diffuse.HealpixCube(file) for file in files]
        for dm in self.dmaps: dm.setEnergy(self.energy)
        self.nside = nside
        if nside is not None and nside != self.dmaps[0].nside:
            pass #print 'Must reformat from {} to {}'.format(self.dmaps[0].nside, nside)
        else: self.nside = self.dmaps[0].nside
        self.column = None

    @property
    def map(self):
        """Return HEALPix map
        """
        if self.column is None:
            self.column = self.dmaps[0].column(self.energy).copy()
            for dm in self.dmaps[1:]:
                self.column+=  dm.column(self.energy)
        nside = self.dmaps[0].nside
        while self.nside < nside:
            print 'combine pixels in  {} from {} to {}'.format(self.cname, nside, nside/2)
            self.column = hpm.downsize(4.*self.column); 
            nside /=2  
        return self.column
    
    def __call__(self, l,b=0):
        ret = 0
        for dm in self.dmaps:
            ret+=dm(SkyDir(l,b,SkyDir.GALACTIC), self.energy)
        return ret
    
    def ait_plot(self, vmin=0, vmax=None, log=False, ax=None, cmap=plt.get_cmap('coolwarm'),
                 **kwargs):


        y = hpm.HParray('', self.map)
        y.plot(title='{} @ {:.2f} GeV'.format(self.cname, self.energy/1e3), ait_kw=dict(), 
               log=log, vmin=vmin, vmax=vmax,cmap=cmap, **kwargs).grid();
    
    def plot_along_plane(self, ax=None):
        x = np.linspace(-120, 120, 1000)
        y = map(self, x)
        fig, ax = plt.subplots(figsize=(8,4)) if ax is None else (ax.figure, ax)
        ax.plot(x,y, '-', label=self.cname);
        ax.grid(alpha=0.5)
        ax.set_xticks(np.linspace(-120,120,9))
        ax.set(xlim=(120,-120), xlabel='$l$', 
               title='{} along the plane @ {:.2f} GeV'.format(
                   self.cname, self.energy/1e3), ylabel='Counts/pixel');
    def plot_along_meridian(self, ax=None):
        x = np.linspace(-90, 90, 1000)
        y = map(lambda b: self(0,b), x)
        fig, ax = plt.subplots(figsize=(8,4)) if ax is None else (ax.figure, ax)
        ax.plot(x,y, '-', label=self.cname);
        ax.grid(alpha=0.5)
        ax.set_xticks(np.linspace(-90,90,13))
        ax.set(xlim=(-90,90), xlabel='$b$', 
               title='{} along l=0 @ {:.2f} GeV'.format(
                   self.cname, self.energy/1e3), ylabel='Counts/pixel');

class GarDian(dict):
    def __init__(self, energy=131., modelname='RH01',nside=64):
    
        root = '/nfs/farm/g/glast/g/diffuse/P8diffuse/'
        model_path = root+'results/gardian/8yr/InnerGalaxy{}_ModelComponentsCounts/'.format(modelname)
        final_path = root+'results/gardian/8yr/InnerGalaxy{}_FinalModels/'.format(modelname)

        self.energy=energy
        self.nside = nside
        
        print 'loading diffuse model {} from files with energy {:.2f} GeV from\n {}'.format(modelname,energy/1e3, model_path)
        clist = dict(patch='patch', dnm='DNM', h1='HI', co='CO', ic='IC', src='Sources', UN='UNRESOLVED', egb='EGBfree')
        for name, cname in clist.iteritems():
            self[name] = GComponent(cname, model_path, energy=energy, nside=nside)

        self['final'] = GComponent('', final_path, energy=energy, nside=nside)

        sdir = map(Band(nside).dir, range(12*nside**2))
        self.glat = map(lambda x: x.b(), sdir)
        self.glon = map(lambda x: x.l(), sdir)

        data_path =  root+'data/fermi_data/gardian_countMaps/SD_P305_4/'
        print 'loading counts map for energy {:.2f} GeV from \n {}'.format(energy/1e3, data_path)
        self['counts'] = counts = GComponent('countsMap', data_path, energy=energy, nside=self.nside)
        plt.rc('font', size=14)

    def write(self, filename):
        """generate a FITS file with all the (combined) maps for this energy
        """
        hh = [hpm.HParray(name, value.map) for name,value in self.iteritems()]
        t=hpm.HEALPixFITS(hh)
        t.write(filename)

    def along_plane(self, components, ylim=(0.1,None)):
        fig, ax = plt.subplots(figsize=(12,6))
        for name in components:
            self[name].plot_along_plane(ax=ax)
        ax.legend()
        ax.set_title('Diffuse Components along the plane at {:.2f} GeV'.format(self.energy/1e3))
        ax.set(yscale='log', ylim=ylim)
        print
    
    def along_meridian(self, components):
        fig, ax = plt.subplots(figsize=(12,6))
        for name in components:
            self[name].plot_along_meridian(ax=ax)
        ax.legend();
        ax.set_title('Diffuse Components along l=0 at {:.2f} GeV'.format(self.energy/1e3))
        ax.set(yscale='log', ylim=(0.1,None));
        print

def find_models(n=5):
    """Return dataframe with model names and creation times
    """
    import glob, time
    os.chdir('/nfs/farm/g/glast/g/diffuse/P8diffuse/results/gardian/8yr')
    s = glob.glob('InnerGalaxy*_ModelComponentsCounts')
    u = [(x[11:15], os.stat(x).st_mtime) for x in s]
    
    df = pd.DataFrame(u, columns='name mtime'.split()); 
    df['time']= [time.ctime(x) for x in df.mtime]
    df.index = df.name
    ret = df.sort_values(by='mtime', ascending=False)['time'.split()]
    print 'Most recent Gardian models\n{}'.format(ret.head(n))
    return ret
    

class GComp(object):
    """ A Galactic Diffuse components
    Loaded from file created by GarDian.

     """
    def __init__(self, modelname, energy=131, 
            root='/nfs/farm/g/glast/g/catalog/pointlike/fermi/gardian/' , reload=False):
        from astropy.table import Table
 
        cfile = '{}{}_{:d}.fits'.format(root, modelname, energy)
        if not os.path.exists(cfile) or reload:
            print 'Will create file {}'.format(cfile)
            GarDian(modelname=mname, energy=energy).write(cfile)

        assert os.path.exists(cfile), 'File "{}" not found'.format(cfile)
        self.df=Table.read(cfile, hdu=1).to_pandas(); 
        self.modelname, self.energy = modelname, energy
        self.nside = nside=int(np.sqrt(len(self.df)/12))
        print 'read file {}'.format(cfile)

        # add locations to the DataFrame
        glon,glat=healpy.pix2ang(nside, range(12*nside**2), lonlat=True)
        glon[glon>180] -= 360
        self.df['glon']=glon
        self.df['glat']=glat

        # a column with combined sources
        self.df['sources']=self.df['src']+self.df['UN']
    
        data,model = np.array(self.df['counts']), np.array(self.df['final'])
        resid = (data-model)/np.sqrt(model)
        #print resid.mean(), resid.std()

        self.pulls = hpm.HParray('f{} pulls'.format(energy), (data-model)/np.sqrt(model))
        self.df['pull']= self.pulls.vec
        self.label = '{} @ {:.0f} MeV'.format(self.modelname, self.energy)
        
    def __repr__(self):
        s = 'Gardian model {} for {} MeV\n\t{:10}{:>15}'.format(self.modelname,self.energy, 'column', 'sum')
        for col in self.df.columns:
            t = self.df[col]
            s += '\n\t{:10}{:15,}'.format(col, int(t.sum()))
        return s

    def __getitem__(self, name):
        """ returns an HParray object corresponding to the name
        """
        if name=='galactic':
            return hpm.HParray(name, self.df['final']-self.df['src'])
        if name=='isotropic': 
            return hpm.HParray(name, self.df['egb'])
        if name=='sources' : return hpm.HParray(name, self.df['src']+self.df['UN'])
        
        return hpm.HParray(name, self.df[name])

    
    def along_plane(self, components, ylim=(0.1,None)):
        fig, ax = plt.subplots(figsize=(12,6))
        for name in components:
            plot_along_plane(self[name], name, ax=ax)
        ax.legend()
        ax.set_title('{} Components along the plane at {:.2f} GeV'.format(self.modelname,self.energy/1e3))
        ax.set(yscale='log', ylim=ylim)
        print

    def pulls_hist(self, ax=None, cuts=['-5<glat<5', 'glat>10 | glat<-10'], 
                   colors='green orange blue'.split()):
        """Normalized residual plots
        """
        fig, ax = plt.subplots(figsize=(6,4)) if ax is None else (ax.figure, ax) 
        
        sel = [self.df.query(cut).pull for cut in cuts]
        
        def plotit(s, label, color): 
            hkw=dict(bins=np.linspace(-5,5,51), histtype='step', lw=2, log=True)
            label='{:20} {:5.2f} {:5.2f}'.format(label, s.mean(),s.std())
            ax.hist(s.clip(-5,5),label=label,color=color, **hkw )
            #ovelay a gaussian with same total
            g=lambda x: np.exp(-x**2/2.)
            x = np.linspace(-4,4,81)
            b =hkw['bins']; delta= b[1]-b[0]
            norm = len(s) * delta/np.sqrt(2*np.pi)
            ax.plot(x, norm*g(x), '--', color=color) 
        
        for s,label, color in zip(sel, cuts, colors):   
            plotit(s, label, color)

        ax.grid(alpha=0.5)
        
        leg=ax.legend(loc='lower center',
            title='      {:20} {:5} {:5}'.format('selection', 'mean', 'std'),
                prop=dict(size=10, family='monospace'))
        ltit = leg.get_title()
        ltit.set_fontsize(10); ltit.set_family('monospace')
        ax.set(ylim=(0.8,None), xlim=(-5,5))
        ax.set_title('Normalized residuals for {}'.format(self.label));
        fig.set_facecolor('w')
        return fig

    def ait_pulls_plot(self):
        fig,ax=plt.subplots(figsize=(16,8))
        t=self.pulls.plot(axes=ax, vmin=-5, vmax=5, cmap=plt.get_cmap('coolwarm'),
            title='{} normalized residuals'.format(self.label))
        t.grid(color='white');

    def ait_plot(self, comp, fraction=False, **kwargs):
        fig,ax=plt.subplots(figsize=(16,8))
        if not fraction:
            defaults = dict(log=True, vmin=1, vmax=1000., title=kwargs.get('title',comp),)
            defaults.update(**kwargs)
            t =self[comp].plot(axes=ax,  **defaults)
        else:
            defaults=dict(vmin=0, vmax=100.); defaults.update(**kwargs)
            a = self[comp].vec; b = self['final'].vec - self['src'].vec
            title = '{} / total (%)'.format(comp)
            t=hpm.HParray(title, a/b*100.).plot(axes=ax,title=title, cbtext='%',**defaults)
        t.grid(color='lightgrey')

    def pie(self, ax=None, query=None):   
        cols = 'h1 co patch dnm  ic egb src UN'.split()
        df = self.df if query is None else self.df.query(query)
        sizes = np.array([sum(df[col]) for col in cols  ])

        # combine UN and src
        sizes[-2] += sizes[-1]
        cols[-2] = 'src+UN'
        explode= [0]*6 + [0.1]*2
        fig, ax = plt.subplots(figsize=(6,6)) if ax is None else (ax.figure, ax)
        ax.pie(sizes[:-1], explode=explode[:-1], labels=cols[:-1], autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax.axis('equal'); # Equal aspect ratio ensures that pie is drawn as a circle.
        ax.set_title('XA01 component sizes for 131-MeV band')

class ExposureMap(object):
    """Manage an exposure map
    """
    
    def __init__(self,  event_type,
        pattern='expcube_car_P8_P305_8years_zmax90_P8R3_SOURCE_V2_*.fits',
        exproot = '/nfs/farm/g/glast/g/diffuse/P8diffuse/data/fermi_data/expcube10/' ):
            """
            """
            found = glob.glob(exproot+pattern)
            assert len(found)>0, 'No exposure files fit the pattern'
            filename = exproot+pattern.replace('*',event_type)
            assert os.path.exists(filename), 'Event type "{}" not found: possible were \n {}'.format(event_type,found)
            self.map = diffuse.FitsMapCube(filename)

    def __call__(self, energy, nside):
        """ return a HEALPix array for the energy, binned to nside"""
        self.map.setEnergy(energy)
        return hpm.HPskyfun('', self.map, nside).getcol()
    
    def plot(self, energy, ax=None, **kwargs):
        """ Make an AIT plot for the given energy"""
        self.map.setEnergy(energy)
        self.map.plot_map(axes=ax, **kwargs)
        