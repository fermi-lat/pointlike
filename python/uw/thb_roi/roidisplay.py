"""

$Header$
"""

import numpy as np
from skymaps import PySkyFunction, SkyDir
from uw.utilities import image, colormaps
import matplotlib
import pylab as plt

from uw.like.roi_plotting import ROIModelSkyFunction, DataSkyFunction, ppf

        


class ROIDisplay(object):
    """Manage the plotting of ROI info."""

    def init(self): 
        self.figsize = (12,8)
        self.fignum  = 3
        
        self.nside = 2**7
        self.pixelsize = 0.25
        self.size = 10
        self.nticks = 5

        self.log_counts_min = 1.0
        self.log_counts_max = 2.5

        self.label_sources = False

        self.galactic = True

        self.std_weight = True

        self.emin = 1e2
        self.emax = 1e5

        self.event_class = -1

        self.mm = self.cm = None

    def mathrm(self,st):
        return  r'$\mathrm{'+st.replace(' ','\/')+'}$'

    def __init__(self, roi_manager, **kwargs):
        
        self.init()
        self.__dict__.update(kwargs)
        self.rm = roi_manager

        from matplotlib import rcParams
        rcParams['xtick.major.size']=10 #size in points
        rcParams['xtick.minor.size']=6
        rcParams['ytick.major.size']=10 #size in points
        rcParams['ytick.minor.size']=6
        rcParams['xtick.labelsize']=12
        rcParams['ytick.labelsize']=12
        rcParams['font.family']='serif'
        #rcParams['text.usetex']=True
    
        import pylab as P
        from matplotlib import mpl,pyplot,ticker
        try:
            self.cmap_sls = colormaps.sls
            self.cmap_b    = colormaps.b
        except:
            self.cmap_sls = self.cmap_b = mpl.cm.jet

        P.ioff()
        fig = P.figure(self.fignum,self.figsize)
        P.clf()

        voff = -0.05
        imw = ( 0.55, 0.55, 0.55,  0.55,      0.25, 0.25) #0.38, 0.38)
        imh = ( 0.35, 0.35, 0.35,  0.35,      0.17, 0.17) # picture sizes
        imx = (-0.09, 0.23, 0.55, -0.09,      0.40, 0.40) #0.55, 0.55)
        imy = ( 0.60+voff, 0.60+voff, 0.60+voff,  0.15+voff,      0.32, 0.08)


        titles  = ['Modeled Counts', 'Observed Counts', 'Weighted Residuals', 'P-values','', '']
        xlabels = ['RA']*4 + ['p-values','weighted residuals'] 
        ylabels = ['Dec']*4 + ['Frequency','Frequency']

        for y in [titles,xlabels,ylabels]:
            for i in xrange(len(y)):
                y[i] = self.mathrm(y[i])

        for i,(t,xl,yl) in enumerate(zip(titles,xlabels,ylabels)):
            ax = self.__dict__['axes%d'%(i+1)] = fig.add_axes([imx[i],imy[i],imw[i],imh[i]])
            if i < 4: ax.set_aspect(1)
            ax.set_title(t)
            ax.set_xlabel(xl) 
            ax.set_ylabel(yl)

        self.norm1 = mpl.colors.Normalize(vmin=0,vmax=5)
        self.norm2 = mpl.colors.Normalize(vmin=self.log_counts_min,vmax=self.log_counts_max)
        self.norm3 = mpl.colors.Normalize(vmin=-5,vmax=5)

        #resid colorbar
        rcb_axes = fig.add_axes([imx[3]+imw[3]*0.72,imy[3],0.01,imh[3]])
        #rcb_norm = mpl.colors.Normalize(vmin=0,vmax=5)
        rcb        = mpl.colorbar.ColorbarBase(rcb_axes,
                                                            norm=self.norm1,
                                                            ticks = ticker.MaxNLocator(4),
                                                            cmap  = colormaps.b,
                                                            orientation='vertical')
        #rcb.set_label('')

        #counts colorbar
        ccb_axes1 = fig.add_axes([imx[0]+imw[0]*0.72,imy[0],0.01,imh[0]])
        ccb1        = mpl.colorbar.ColorbarBase(ccb_axes1,
                                                            norm=self.norm2,
                                                            ticks = ticker.MaxNLocator(4),
                                                            cmap  = colormaps.b,
                                                            orientation='vertical')
        ccb_axes2 = fig.add_axes([imx[1]+imw[1]*0.72,imy[1],0.01,imh[1]])
        ccb2        = mpl.colorbar.ColorbarBase(ccb_axes2,
                                                            norm=self.norm2,
                                                            ticks = ticker.MaxNLocator(4),
                                                            cmap  = colormaps.b,
                                                            orientation='vertical')

        ccb_axes3 = fig.add_axes([imx[2]+imw[2]*0.72,imy[2],0.01,imh[2]])
        ccb3        = mpl.colorbar.ColorbarBase(ccb_axes3,
                                                            norm=self.norm3,
                                                            ticks = ticker.MaxNLocator(4),
                                                            #cmap  = colormaps.b,
                                                            orientation='vertical')
        #rcb.set_label('')

        self.mm = ROIModelSkyFunction(roi_manager,mode=1,**self.__dict__) if self.mm is None else self.mm
        self.cm = DataSkyFunction(roi_manager.sa,mode=1,**self.__dict__) if self.cm is None else self.cm

    def model_plot(self):

        self.mm_zea = image.ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes1)
        self.mm_zea.fill(PySkyFunction(self.mm))
        self.axes1.imshow(np.log10(self.mm_zea.image),origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm2)
        self.mm_zea.grid()
        self.mm_zea.scale_bar(color='white')
        self.plot_sources(self.mm_zea,mc='k')
        
    def counts_plot(self):

        self.cm_zea = image.ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes2)
        self.cm_zea.fill(PySkyFunction(self.cm))
        self.axes2.imshow(np.log10(self.cm_zea.image),origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm2)
        self.cm_zea.grid()
        self.cm_zea.scale_bar(color='white')
        self.plot_sources(self.cm_zea,mc='k')

    def resids_plot(self):

        from scipy.stats import poisson

        self.resids_zea = image.ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes3)
        self.resids_zea.image = (self.mm_zea.image - self.cm_zea.image)/self.mm_zea.image**0.5
        self.axes3.imshow(self.resids_zea.image,origin='lower',interpolation='nearest',norm=self.norm3)
        self.resids_zea.grid()
        self.plot_sources(self.resids_zea,mc='k')
        
        pvals = 1 - poisson.cdf(self.cm_zea.image,self.mm_zea.image ) #0 problem?
        pvals = np.abs(np.where(pvals < 0.5, np.log10(pvals), np.log10(1-pvals)))
        #pvals = np.abs( np.tan ( pvals*np.pi - np.pi/2 ) )

        self.pvals_zea = image.ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes4)
        self.pvals_zea.image = pvals
        self.axes4.imshow(self.pvals_zea.image,origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm1)
        self.pvals_zea.grid()

        self.plot_sources(self.pvals_zea,mc='k')

    def hist_plot(self):

        import pylab as P
        from scipy.stats import poisson

        mc = np.asarray(self.mm.cache_pix.values())
        cc = np.asarray(self.cm.cache_pix.values())

        pvals = np.asarray([1-poisson.cdf(cc[i],mc[i]) for i in xrange(len(mc))])

        nb = 20
        av = float(len(pvals)) / nb
        self.axes5.hist(pvals,bins=np.linspace(0,1,20),histtype='step')
        self.axes5.axhline( av, color='red')
        lo,hi = ppf( (50.-95/2)/100., av), ppf( (50. + 95/2)/100.,av)
        self.axes5.axhspan( lo , hi , facecolor='red', alpha=0.3,label='95% Conf.')
        self.axes5.legend(loc='upper right')

        self.axes6.hist( (mc - cc)/mc**0.5, bins=np.linspace(-5,5,20), histtype='step')
        self.axes6.axvline(0,color='red')
        self.axes6.grid(True)
        self.axes6.set_xbound(lower=-5,upper=5)


    def show(self,to_screen=True,out_file=None):

        t =self.label_sources
        self.model_plot()
        self.label_sources=False #only label the model plot
        self.counts_plot()
        self.resids_plot()
        #self.hist_plot()
        self.label_sources = t

        import pylab as P
        if out_file is not None: P.savefig(out_file)
        if to_screen: P.show()

    def plot_sources(self, image, symbol='+',  fontsize=8, markersize=10, fontcolor='w', mc= 'green',**kwargs):
        nx = image.nx
        ny = image.ny
        ps = self.rm.psm.point_sources

        def allow(nx, ny, px, py, padx = 0.15, pady = 0.15):
            padx = padx * nx
            pady = pady * ny

            return (px > padx and px < nx - padx) and (py > pady and py < ny - pady)

        for p in ps:
            x,y = image.pixel(p.skydir)
            if not allow(nx,ny,x,y,0.02,0.02): continue
            image.axes.plot([x], [y], symbol, markersize=markersize, mec = mc, mfc = mc,**kwargs)
            image.axes.plot([x], [y], 'x', markersize=markersize, mec = 'white', mfc = 'white',**kwargs)
            if self.label_sources:
                dx, ha = (nx/50., 'right') if x>nx/2 else  (-nx/50., 'left')
                image.axes.text( x+dx, y, p.name, fontsize=fontsize, color=fontcolor, ha=ha)

