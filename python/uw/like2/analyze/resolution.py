"""
Resolution study

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/resolution.py,v 1.1 2017/11/17 22:44:09 burnett Exp $

author:  Toby Burnett
"""
import os
import numpy as np
from skymaps import SkyDir 
import pandas as pd
import matplotlib.pyplot as plt

from uw.like2 import process
from uw.like2.plotting import sed
from uw.irfs import irfman
from uw.utilities import image
etname=irfman.IrfManager.event_type_names

class Resolution(object):
    def __init__(self, source_name, modeldir='.'):
        self.sources = pd.read_pickle(os.path.join(modeldir,'sources.pickle'))
        self.sinfo = self.sources.ix[source_name]
        self.roi = int(self.sinfo['roiname'][-4:])
        self.r = process.Process('.', quiet=True)
        self.event_types =self.r.config.dataset.event_types
        print ('Loading ROI # {}'.format(self.roi))
        self.r.setup_roi(self.roi)
        self.source = self.r.get_source(source_name)
        plt.rcParams.update({'font.size': 16, 'axes.labelsize':16, 'axes.titlesize':16, 
          'figure.titlesize':20,'legend.fontsize':16, 'xtick.labelsize':16, 'ytick.labelsize':16,
                    'figure.facecolor':'white'}
                   )

    def sed_plots(self, ylim=(2,400)):
        s=self.source
        event_types = self.event_types
        srec = [self.r.get_sed(s.name, event_type=et, update=True) for et in event_types]
        fig,axx = plt.subplots(len(event_types)/2,2, figsize=(10, 2*len(event_types)), sharex=True, sharey=True)
        for ax, sedrec,et in zip(axx.flatten(), srec, event_types):
            sed.Plot(s, sedrec)(butterfly=False,axes=ax, axis=(1e2,1e6, 2, 80))#40, 400))
            ax.set_title(etname[et])
            ax.set_ylim(ylim)
        fig.suptitle(s.name, fontsize=14);

    def data_frame(self):
        s=self.source
        r=self.r
        event_types=self.event_types
        u = dict()
        pos = dict()
        def doit(et='all'):
            r.localize()
            t =s.ellipse
            pos[et]=SkyDir(t[0],t[1])
            r.fit(s.name)
            stats = np.array(s.model.statistical())[:,:2].flatten()
            if et=='all':
                t.append(0)
            else:
                t.append(np.degrees(pos['all'].difference(pos[et]))*3600)
            u[et if et=='all' else etname[et]]=[r.TS()]+list(stats)+t
        doit('all')    
        for et in event_types:
            r.select(event_type=et)
            doit(et)
        r.select()
        colnames='ts flux pindex flux_unc pindex_unc ra dec a b ang qual un delta'.split()
        v = pd.DataFrame(u, index=colnames, ).T
        v['sig']= np.sqrt(v.a*v.b)*3600
        
        if not os.path.exists('resolution'): os.mkdir('resolution')
        outfile = 'resolution/{}.pkl'.format(s.name.replace(' ','').replace('+','p'))
        v.to_pickle(outfile)
        print ('wrote resolution file {}'.format(outfile))
        self.v =v

class Plots(object):
    def __init__(self, *args):
        self.args=args

    def loc_pictures(self, size=0.02):
        sname, fb_res,psf_res=self.args
        fig,axx = plt.subplots(1,2, figsize=(10,5), sharey=True)
        pos = SkyDir(* [fb_res.ix['all'][x] for x in ('ra','dec')])
        for v, ax in zip((fb_res,psf_res), axx.flatten()):
            lw=4
            zea = image.ZEA(pos, axes=ax, size=size, pixelsize=size/15)
            zea.plot([pos])
            for name in v.index:
                ra,dec,a,b,ang = v.ix[name][5:10]
                zea.ellipse((ra,dec), (a,b,ang),lw=lw)
                lw=2
            zea.scale_bar(delta=10./3600., text='10 arcsec')
        fig.suptitle('Localization for {}'.format(sname), fontsize=16)

    def spectral_figures(self):
        sname, fb_res,psf_res=self.args
        fig,axx = plt.subplots(1,2, figsize=(12,5), sharex=True, sharey=True)
        def spectral_fig(v, ax, **kwargs):
            lw=4
            for name, w in v.iterrows():
                x,y = w['flux'], w['pindex']
                ax.errorbar(x,y, 
                        xerr=x*w['flux_unc'], yerr=y*w['pindex_unc'],
                        fmt='o', label=name, lw=lw, **kwargs);
                ax.text(x,y,name, fontsize=14)
                lw=2
            #ax.legend()
            ax.grid(alpha=0.5)
            ax.set_xlabel('Flux [cm^-2 s^-1]')
            ax.set_ylabel('Photon index')
        spectral_fig(psf_res,axx[0])
        spectral_fig(fb_res,axx[1])
        fig.suptitle('Spectral fit for {}'.format(sname), fontsize=14);

    def resolution_bar(self):
        sname, fb_res,psf_res=self.args
        x = [0,0.25]
        xlabel = ['PSF','FB']
        fig, ax= plt.subplots(figsize=(5,5))
        width=0.1
        def bar(df, x=0):
            bottom=0
            n = len(df)
            res = 1/df.sig**2
            for j in range(n-1,0,-1): #i in range(n-1):
                y,label = res[j], df.index[j]
                #print (y,bottom,  label)
                ax.bar(x-width/2, y, bottom=bottom, width=width, color='lightcyan',lw=2)
                ax.text(x, bottom+y/2, label, va='center', ha='center')
                bottom+=y

        bar(psf_res, x[0])
        bar(fb_res, x[1])
        #ax.legend()
        ax.set_xlim([x[0]-0.12,x[1]+0.12])
        ax.set_ylabel('Resolution [1/arcsec^2]')
        ax.set_xlabel('Event types used')
        plt.xticks(x, xlabel);
        ax.set_title('Localization resolution for {}'.format(sname))

    def resolution_piecharts(self):
        sname, fb_res,psf_res=self.args
        fig, axx = plt.subplots(2,2, figsize=(14,12))
        colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
        
        def pie(ax,labels, sig, title):
            z = 1/sig**2
            ax.pie(z[1:]/z[0], labels=labels, colors=colors, shadow=True, autopct='%1.0f%%',)
            ax.set_aspect('equal')
            ax.set_title(title)
            
        for v, axes in zip((fb_res, psf_res), (axx[0,:], axx[1,:])):
            for ax, t, title in zip(axes, 
                                (v.flux_unc, v.sig), 'Flux Localization'.split()):
                pie(ax,v.index[1:], t, title)
        fig.suptitle('Resolutions for {}'.format(sname))

 

