"""
Galactic diffuse refit spectra

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/galacticspectra.py,v 1.2 2013/08/20 19:50:23 burnett Exp $

"""

import os
import numpy as np
import pylab as plt
import pandas as pd

from . import roi_info
from . analysis_base import FloatFormat

class GalacticSpectra(roi_info.ROIinfo): #Diagnostics):
    """Galactic diffuse refits
      <br>Set of plots to check consistency of %(title)s spectra. These result 
        from analysis of a special run that, for each ROI and each energy band, allows this diffuse component to be free.
        This is done three times: using only front, only back, and both.
        <p>There two sets of plots: using both, how consistent is it with the expected unit normalization; 
        and is the front        consistent with the back?
        """

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
        self.rois = self.df = pd.DataFrame(dict(glat=glat, glon=glon, singlat=singlat), 
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
                                                    
        self.plot_functions =[ [self.like_scats, self.diffuse_fits, self.normalization_factor_scats, self.bfratio_hists, self.diffuse_ratio_plot, ],
                    map( lambda s: self.which+'_'+s, ['likelihood_diff', 'diffuse_fits', 'diffuse_normalization_factor', 'bfratio', 'diffuse_ratio', ]),
                    ]

    def setup(self, args=None):
        self.diffuse_setup('gal')
    #### these from old diagnostics: do I want them in AnalysisBase?
    def get_figure(self,axin):
        return plt.subplots( figsize=(5,5)) if axin is None else axin.figure, axin
    def multifig(self):
        fig,ax = plt.subplots(2,4, figsize=(14,8), sharex=True);
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.95)
        return ax.flatten()
    def multilabels(self, xtext, ytext, title=None):
        plt.subplots_adjust(bottom=0.2)
        plt.figtext(0.5,0.07, xtext, ha='center');
        plt.figtext(0.05, 0.5, ytext, rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
    def set_plot(self, ax, fignum, figsize=(4,4)):
        if ax is None:
            plt.close(fignum)
            plt.figure(fignum, figsize=figsize);
            ax = plt.gca()
        else:
            plt.sca(ax); 
        return ax
    ##################
        
    def like_scat(self, ib, axin=None, fb='both', vmin=0, vmax=1):
        fig, ax=self.get_figure(axin); 
        scat=ax.scatter(self.rois.glon, self.rois.singlat, 
                c=np.log10(self.flux[fb]['deltalike'].transpose().ix[ib]), 
                s=25 if axin is not None else 40,
                vmin=vmin, vmax=vmax, edgecolors='none')
        ax.set_title('%.0f MeV'%(self.energy[ib]),fontsize='small')
        plt.setp(ax, xlim=(180,-180), ylim=(-1.01, 1.01))
        ax.set_xticks([180,90,0,-90,-180])
        return scat
    
    def save_correction(self, filename='galactic_correction.csv', which='both'):
        """ a bit of code to create a csv file with the fit values, index by roiname
        """
        x = self.flux[which]['values']
        x.index.name='roiname'
        x.to_csv(filename)
        print 'wrote file %s' % filename
    
    def like_scats(self, title=None):
        """ Likelihood ratios for individual fits.
        These all-sky plots show, for each ROI and each energy band, the consistency of the %(title)s spectral fit 
        to a value determined for just that energy band. The distribution of the log likelihood should be approximately 
        the chi squared distribution of one degree of freedom. The lighter colors, especially red, indicate serious discrepancy.
        """
        fig,ax = plt.subplots(2,4, figsize=(14,8), sharex=True, sharey=True);
        plt.subplots_adjust(left=0.10, wspace=0.25, hspace=0.25,right=0.90, bottom=0.15,)
        scats =map(self.like_scat, range(8), ax.flatten());
        plt.figtext(0.5,0.07, 'glon', ha='center');
        plt.figtext(0.05, 0.5, 'sin(glat)', rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scats[0], cbax, orientation='vertical')
        cb.set_label('log10(log likelihood difference)')
        return fig
    def sky_scat(self, c, axin=None, vmin=0, vmax=1, title=None):
        fig, ax=self.get_figure(axin); 
        scat=ax.scatter(self.rois.glon, self.rois.singlat, 
                c=c, 
                s=25 if axin is not None else 40, 
                 marker='D',
                vmin=vmin, vmax=vmax, edgecolors='none')
        plt.setp(ax, xlim=(180,-180), ylim=(-1.01, 1.01))
        ax.set_xticks([120, 60 ,0,-60,-120])
        if title is not None: ax.set_title(title, size=12)
        return scat
        
    def sky_scats(self, v,  title=None, vmin=None, vmax=None, cb_label=None):
        fig,axx = plt.subplots(2,4, figsize=(14,8), sharex=True, sharey=True);
        plt.subplots_adjust(left=0.10, wspace=0.1, hspace=0.15,right=0.90, bottom=0.15)
        scats =[self.sky_scat( v[ib], axin=ax, vmin=vmin, vmax=vmax, title='%.0f MeV'%self.energy[ib]) for ib,ax in enumerate(axx.flatten())]
        plt.figtext(0.5,0.07, 'glon', ha='center');
        plt.figtext(0.05, 0.5, 'sin(glat)', rotation='vertical', va='center')
        if title is not None: plt.suptitle(title)
        cbax = fig.add_axes((0.92, 0.15, 0.02, 0.7) )
        cb=plt.colorbar(scats[0], cbax, orientation='vertical')
        if cb_label is not None: cb.set_label('fit value')
        return fig        

    def normalization_factor_scats(self, vmin=0.75, vmax=1.25):
        """Normalization factors
        The fit normalization factor for each ROI and the first eight energy bands
        """
        return self.sky_scats( self.flux['both']['values'], vmin=vmin, vmax=vmax, cb_label='fit value')
    
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
        return plt.gcf()
  
    def diffuse_ratio_plot(self):
        """ Front/back %(title)s diffuse ratio
        The front to back ratio, measured from the front/back fit ratios.
        
        <br>%(fb_ratio)s
        """
        fig, ax = plt.subplots( figsize=(4,4))
        plt.subplots_adjust(left=0.2, bottom=0.2) #not sure why this is necessary
        vals = self.bfratios # must have generated
        x,y,yerr = [[v[i] for v in vals]  for i in range(3)]
        ax.errorbar(x, y, yerr=yerr, marker='o', ms=12,fmt='', lw=2, linestyle='None')
        plt.setp(ax, xscale='log',xlabel='Energy (MeV)', ylabel='front/back flux ratio',ylim=(0.75, 1.25))
        ax.grid(True)
        ax.axhline(1.0, color='k')
        ax.set_title('%s diffuse spectral fits'%self.which, fontsize='medium')
        return fig
        
    def diffuse_fit(self, axin, ind=0,  fignum=2, xlim=(0.5,1.5), **kwargs):
    
        plane = abs(self.rois.glat)<10
        cut, cut_name = (plane, 'plane') if self.which=='gal' else (abs(self.rois.glat)>20, 'high-lat')
        space=np.linspace(xlim[0],xlim[1], 41)
        vals = self.flux['both']['values'].transpose().ix[ind] #values[:,ind]
        kw = dict( histtype='stepfilled')
        kw.update(kwargs)
        ax = self.set_plot(axin, fignum)
        ax.hist(vals.clip(*xlim), space,label='all', **kw)
        ax.hist(vals[cut].clip(*xlim), space ,color='r',label=cut_name,**kw);
        ax.grid(True);
        ax.axvline(1.0, color='grey');
        ax.set_xlim(xlim)
        ax.set_title('%.0f MeV'%self.energy[ind],fontsize='medium')
        #ax.legend(prop=dict(size=10))
        
    class DiffuseFit(object):
    
        def __init__(self, spectra,  xlim=(0.5,1.5), nbins=40, **kwargs):
            plane = abs(spectra.rois.glat)<10
            self.cut, self.cut_name = (plane, 'plane') if spectra.which=='gal' else (abs(spectra.rois.glat)>20, 'high-lat')
            self.space=np.linspace(xlim[0],xlim[1], nbins+1)
            self.xlim = xlim
            self.kw = dict( histtype='stepfilled')
            self.kw.update(kwargs)
            self.flux = spectra.flux
            self.energy = spectra.energy
            
        def __call__(self, ax, ind):
            vals = self.flux['both']['values'].transpose().ix[ind] #values[:,ind]
            ax.hist(vals.clip(self.xlim), self.space,label='all', **self.kw)
            ax.hist(vals[self.cut].clip(self.xlim), self.space ,color='r', label=self.cut_name,**self.kw);
            ax.grid(True);
            ax.axvline(1.0, color='grey');
            ax.set_xlim((self.space[0],self.space[1]))
            ax.set_title('%.0f MeV'%self.energy[ind],fontsize='medium')
            ax.legend(prop=dict(size=10))

            
    def fit_map(self, axin=None, ib=0, fignum=2, vmin=0.75, vmax=1.25, **kwars):
        vals = self.flux['both']['values'].transpose().ix[ib] 
        ax = self.set_plot(axin, fignum)
        ax.scatter(self.rois.glon, self.rois.singlat, 
                c=vals,
                s=15 if axin is not None else 25,
                vmin=vmin, vmax=vmax, edgecolors='none')
        ax.set_title('%.0f MeV'%(self.energy[ib]),fontsize='small')
        plt.setp(ax, xlim=(180,-180), ylim=(-1.01, 1.01))
        ax.set_xticks([180,90,0,-90,-180])
        return ax.figure
    
        
    def diffuse_fits(self, **kw):
        """ %(title)s normalization
        For the eight lowest energy bands, the normalization factors, for combined front and back.
        %(normalization_table)s
        """
        ax = self.multifig()
        self.multilabels('ratio', 'ROIs', '%s diffuse fit' %self.which)
        map(self.diffuse_fit, ax, range(8))
        ###### fix later
        #map(self.DiffuseFit(self), ax, range(8))
        ax[0].set_xticks(np.arange(0.5, 1.51, 0.25))
        ax[0].legend(loc='upper left',bbox_to_anchor=(-0.4,1.2));
        
        # add a table of the means and RMS values for all ROIS, and those within 5 deg.
        xx = self.flux['both']['values']
        plane = np.abs(self.df.glat)<5
        hilat = ~ plane
        av = xx[plane].mean()
        rms=xx[plane].std()
        av_hilat = xx[hilat].mean()
        rms_hilat= xx[hilat].std()
        av_all=xx.mean()
        rms_all=xx.std()
        z=pd.DataFrame(dict([
                     ('mean_plane',av.round(3)),      ('std_plane',rms.round(3)),
                     ('mean_hilat',av_hilat.round(3)),('std_hilat',rms_hilat.round(3)),
                     ('mean_all',av_all.round(3)),    ('std_all',rms_all.round(3)),
                     ]))
        z.index.name='band'
        zhtml = z.to_html(float_format=FloatFormat(3))
        self.normalization_table="""
        <p>Normalization statistics: 'plane' means |b|<5, 'hilat' is not plane <br>%s """ % zhtml
        open('normalization_stats.html','w').write(zhtml)
        print 'wrote HTML file to %s' % 'normalization_stats.html'
        return plt.gcf()
            
    def all_plots(self):
        self.runfigures(*self.plot_functions)