"""
Make various diagnostic plots to include with a skymodel folder

$Header$

"""

import os, pickle, glob
import numpy as np
import pylab as plt

class Diagnostics(object):
    """ basic class to handle data for diagnostics, collect code to make plots
    """
    def __init__(self, skymodel_dir):
        self.skymodel_dir = os.path.expandvars(skymodel_dir)
        assert os.path.exists(os.path.join(self.skymodel_dir, 'config.txt')), 'not a skymodel directory:%s'%skymodel_dir
        os.chdir(self.skymodel_dir)
        files, pkls = self.load_pickles()
        self.rois =  pkls
        assert len(self.rois)==1728
        self.glon = np.array([r['skydir'].l() for r in self.rois]); self.glon[self.glon>180]-=360
        self.glat = np.array([r['skydir'].b() for r in self.rois])
        self.singlat = np.sin(np.radians(self.glat))
        if not os.path.exists('plots'):
            os.mkdir('plots')
    
    def load_pickles(self,folder='pickle'):
        """
            load a set of pickles, return list from either zipfile or folder
        """
        pkls = []
        if os.path.exists(folder+'.zip'):
            print 'unpacking file %s.zip ...' % folder ,
            z=zipfile.ZipFile(folder+'.zip')
            files = sorted(z.namelist())#[1:] # skip  folder?
            print 'found %d files ' % len(files)
            opener = z.open
        else:
           files = sorted(glob.glob(os.path.join(folder,'*.pickle')))
           opener = open
        assert len(files)>0, 'no files found in %s' % folder 
        pkls = [pickle.load(opener(file)) for file in files]
        return files,pkls
        
    def chisq_plot(self,  vmin=0, vmax=100):
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