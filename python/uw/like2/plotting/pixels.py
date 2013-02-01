"""
not done yet
$Header$
"""
import os
import numpy as np
import pylab as plt
from matplotlib import ticker
from .. import roistat, bandlike
from skymaps import Band, SkyDir
from uw.utilities import image

class ROIplot(object):
    """ generate one or more plots of the pixels in an ROI
    """
    def __init__(self, b, roi_dir):
        """ b: an ROIBand object, loaded from the bands property of an ROIAnalysis
                only need the wsdl array
        """
        self.glat=[x.b() for x in b.wsdl]
        t=  np.asarray([x.l() for x in b.wsdl])
        self.glon = np.where(t>180, t-360, t)
        self.npix = len(b.wsdl)
        self.size = 9000/self.npix #default symbol size
        self.band = b
        self.roi_dir=roi_dir
        
    def __call__(self, a, val, classname='', dataname='', **kwargs):
        """ generate a color-coded scatter plot in galactic coordinates 
        a: the Axes object to draw in
        val : an array of values to associate with the pixels
        classname, dataname : strings, put at upper left and right
        **kwargs: passed to scatter, perhaps to set size, 's', or vmin,vmax
        """
        assert len(val)==self.npix, 'length of array, %d, not same as number of pixels, %d' \
                % (len(val), self.npix)
        zea = image.ZEA(self.roi_dir, size=20, axes=a, galactic=True)
        bsk = BandSkyFunction(self.band, val)
        zea.fill(bsk)
        zea.imshow( **kwargs)
        for i, stat in enumerate([
            ('min', val.min()), ('max', val.max()), 
            ('mean',val.mean()), ('total', sum(val)),  ]):
            a.text( 0.05, 0.12-i*0.035, '%-6s%7.1f' %stat, 
                transform = a.transAxes, fontsize=8, fontname='monospace')
        a.text( 0.05, 0.93, classname, transform = a.transAxes, fontsize=8)
        a.text( 0.95, 0.93, dataname,  transform = a.transAxes, ha='right',fontsize=8)

class BandSkyFunction(object):
    def __init__(self, band, values ):
        assert len(values)==len(band.wsdl), 'inconsistent pixel values'
        self.indexfun = Band(band.b.nside()).index
        self.indices = [self.indexfun(w) for w in band.wsdl]
        self.values = values
    def __call__(self,sd):  
        try:
            return self.values[self.indices.index(self.indexfun(sd))]
        except ValueError:
            return np.nan
        
def bplots(roistat, eband=0, saveto=None, title=None):
    """ Make set of (l,b) scatterplots of pixel-related quantities 
    roistat : an ROIstat instance
    """
    fig,ax = plt.subplots(2,6, sharex=True,sharey=True, figsize=(20,10) )
    
    for i in range(2*eband,2*eband+2):
        b = roistat.all_bands[i].band
        pixelarea = b.b.pixelArea()*(180/np.pi)**2 # pixel size in deg**2
        sp = ROIplot(b)
        classname= ['front','back'][i]
        for m in range(6):
            a = ax[i,m]
            val, dataname = (
                (b.bg_pix_counts[:,0], 'galactic'),
                (b.bg_pix_counts[:,1], 'isotropic'),
                (b.frozen_pix_counts, 'frozen sources'),
                (b.ps_all_pix_counts-b.frozen_pix_counts, 'variable sources'),
                (b.pix_counts,  'data'),
                (b.pix_counts-b.bg_all_pix_counts-b.ps_all_pix_counts, 'residual'),
               )[m]
            sp(a, val, classname, dataname,)
    if title  is not None: fig.suptitle(title, fontsize=18);
    if saveto is not None: plt.savefig(saveto)       

def bandmodelplots(bandmodels, roi_dir, saveto=None,title=None, combine_free=True):
    """ Make set of (l,b) scatterplots of pixel-related quantities
        bandmodels : a set of BandModelStat objects
    """
    nbm = len(bandmodels)
    #if combine_free:
    #    freepoint = [x.__class__.__name__=='BandPoint'  for x in bandmodels[0].free_sources]
    #nview =  len(bandmodels[0].free_sources)+3 if not conbine_free else 
    fig,ax = plt.subplots(nbm,nview, sharex=True,sharey=True, figsize=(20,10) ,squeeze=False)
    for i,bm in enumerate(bandmodels):
        sp = ROIplot(bm.band, roi_dir)
        classname= ['front','back'][bm.band.b.event_class() ]
#        if combine_free:
#            freepix = np.zeros(len(bm.data)
#            for freemod in bm.free_sources:
#                freepix += freemod.pix_counts
#            sp(ax[i,m], freemod.pix_counts, classname, freemod.source.name)
#
#            
#        else:
        for m,freemod in enumerate(bm.free_sources):
            sp(ax[i,m], freemod.pix_counts, classname, freemod.source.name)
        sp(ax[i,m+1], bm.fixed_pixels, classname, 'frozen')
        sp(ax[i,m+2], bm.data, classname, 'data')
        resids = (bm.data-bm.model_pixels)/np.sqrt(bm.model_pixels)
        sp(ax[i,m+3], resids, classname, 'residuals', vmin=-3,vmax=3)
        
    if title  is not None: fig.suptitle(title, fontsize=18);
    if saveto is not None: plt.savefig(saveto)
    return fig
    
def residual_plots(bandmodels, roi_dir, saveto=None, title=None):
    """ (l,b) scatterplots of the residuals
    makes a nband/2 by 2 array,
    """
    nbm = len(bandmodels)
    nview = len(bandmodels[0].free_sources)+3
    fig,ax = plt.subplots(2, (nbm+1)/2, sharex=True,sharey=True, figsize=(20,8) ,squeeze=False)
    for j,bm in enumerate(bandmodels):
        sp = ROIplot(bm.band, roi_dir)
        i = bm.band.b.event_class() 
        classname= ['front','back'][i]
        resids = (bm.data-bm.model_pixels)/np.sqrt(bm.model_pixels)
        sp(ax[i, j/2], resids, classname, 'energy %d'%(j/2), vmin=-3,vmax=3)
        
    if title  is not None: fig.suptitle(title, fontsize=18);
    if saveto is not None: plt.savefig(saveto)  
    return fig


class BandPlot(object):

    def __init__(self, roi , ib=8, size=5.):
        """ setup for plots of a band within an ROI. 
        roi : an ROIstat object
        ib :  integer
            index of the band, a BandLike objecrt in roi.selected_bands array
        """
        self.skydir = roi.roi_dir
        self.bl = roi.selected_bands[ib]
        self.classname = ['front','back'][ib%2]
        wsdl = self.bl.band.wsdl
        # use the ZEA object just for the ZEA transformation
        # define 1-degree 'pixels' transform
        zea = image.ZEA(self.skydir, pixelsize=1.0, size=size, galactic=True)
        r = map(zea.pixel, wsdl)
        self.x = [t[0]-size/2 for t in r]
        self.y = [t[1]-size/2 for t in r]
        self.data = [p.weight() for p in wsdl] 
        self.model = self.bl.model_pixels
        
    def _plot(self, vals, counts, valname, ax=None,vmin=None, vmax=None,):
        """ generate a plot, return the figure
        vals : array of values 
            corresponds to pixels with data
        counts : total predicted (for a model) in the ROI
        valname : descriptive name
        ax : Axes object, or None
        """

        if ax is None:
            fig, ax = plt.subplots( figsize=(7,5))
        else: fig =ax.figure
        scat=ax.scatter(self.x, self.y , c =vals, s=5, marker='D', edgecolor='none', vmin=vmin, vmax=vmax)
        ax.set_aspect(1.0)
        plt.setp(ax, aspect=1.0, xlim=(6,-6), ylim=(-8,8))
        #fig.colorbar(scat) #, cax=ax)

        for i, stat in enumerate([
            ('min', vals.min()), ('max', vals.max()), 
            ('mean',vals.mean()), ('total', sum(vals)), ('empty', counts-sum(vals)),  ]):
            ax.text( 0.04, 0.16-i*0.035, '%-6s%7.1f' %stat, 
                transform = ax.transAxes, fontsize=8, fontname='monospace')
        ax.text( 0.05, 0.93, self.classname, transform = ax.transAxes, fontsize=12)
        ax.text( 0.95, 0.93, valname,  transform = ax.transAxes, ha='right',fontsize=12
        )

        return fig

    def __call__(self, im=None, ax=None, vmin=None, vmax=None, ):
        """ generate a plot, return the figure
        im : integer or None
            the index of the source, a BandSource object in the BandLike list
            if None, plot the data
        ax : Axes object, or None
        """
        vals = np.array(self.data if im is None else self.bl[im].pix_counts)
        counts = sum(self.data) if im is None else self.bl[im].counts
        valname = 'data' if im is None else self.bl[im].source.name
        return self._plot(vals, counts, valname, ax=ax, vmin=vmin, vmax=vmax)

    def residuals(self, ax=None):
        return self._plot(self.data - self.model, 0,  'residuals', ax=ax)
        
    def point_sources(self, ax=None):
        vals = np.zeros(len(self.data))
        counts = 0
        for b in self.bl:
            if isinstance(b, bandlike.BandPoint):
                vals += b.pix_counts
                counts += b.counts
        return self._plot(vals, counts, 'all sources', ax=ax)
    def other_diffuse(self, f=2, ax=None):
        vals = np.zeros(len(self.data))
        counts = 0
        names = []
        for b in self.bl[f:]:
            if isinstance(b, bandlike.BandDiffuse):
                vals += b.pix_counts
                counts += b.counts
                names.append(b.source.name)
        return self._plot(vals, counts, ','.join(names), ax=ax)


    
        

#def residual_processor(roi, **kwargs):
#    print 'generating residual plots for', roi.name
#    outdir = kwargs['outdir']
#    rdir = os.path.join(outdir, 'residuals')
#    if not os.path.exists(rdir):   os.mkdir(rdir)
#    s = roistat.ROIstat(roi, lambda b: b.e<1000)
#    saveto = '%s/%s.png'%(rdir,roi.name)
#    residual_plots(s.all_bands, roi.roi_dir, title=roi.name, saveto=saveto)
#    print 'plots saved to ', saveto
#    
#def residual_process(setup, **kwargs):
#    """ a processor to create plots of residuals for all ROIs"""
#    setup['processor'] = 'bandplot.residual_processor'
#    pipeline.pipe.main(setup, **kwargs)