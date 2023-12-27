import numpy as np
import copy
from . roi_image import ModelImage,CountsImage,SmoothedCounts,SmoothedModel,SmoothedResidual,ROITSMapImage
from . pointspec_helpers import PointSource
from . Models import PowerLaw
from . SpatialModels import SpatialMap, PseudoSpatialModel
from uw.utilities import colormaps
from uw.utilities import region_writer
from uw.utilities import keyword_options
from skymaps import SkyDir

from scipy.stats import poisson,norm

import pylab as P
from matplotlib import colors,ticker,font_manager
from matplotlib.ticker import MaxNLocator
from matplotlib.axes import Axes

from . roi_plotting import tight_layout

def set_path_effects(object,**kwargs):
    """ Print a nice warning if withStroke doesn't exist. """
    try:
        from matplotlib.patheffects import withStroke
        object.set_path_effects([withStroke(**kwargs)])
    except ImportError as er:
        print ('Nicer plots require a newer version of matplotlib with patheffects.withStroke!')
        traceback.print_exc(file=sys.stdout)

class ROIMapPlotter(object):
    """ Base class for making a single plot. """

    defaults = (
        ('size',                        5, ),
        ('galactic',                 True, ),
        ('figsize',             (5.5,4.5), 'Size of figure in inches'),
        ('fignum',                   None, 'Matplotlib figure number'),
        ('title',                    None, 'Title of the plot'),
        ('fitsfile',                 None, ), 
        ('show_colorbar',            True, 'Show the colorbar'),
        ('extra_overlay',            None, 'Function which can be used to overlay stuff on the plot.'),
        ('overlay_kwargs',         dict(), 'kwargs passed into overlay_region'),
        ('interpolation',       'nearest', 'passed into imshow'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,**kwargs):

        self.roi=roi

        keyword_options.process(self, kwargs)

        self.pf = self.create_pyfits()

        self.data = self.pf[0].data
        self.header = self.pf[0].header

        if self.fitsfile is not None:
            self.pf.writeto(self.fitsfile,clobber=True)

    def show(self,filename=None, axes=None, cax=None):

        import pylab as P
        import pywcsgrid2
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

        d = self.pf[0].data

        if axes is None:
            fig = P.figure(self.fignum,self.figsize)
            P.clf()
            axes = pywcsgrid2.subplot(111, header=self.header)
        else:
            fig = axes.get_figure()

        self.axes = ax = axes

        im = ax.imshow(d, interpolation=self.interpolation, origin="lower", **self.imshow_kwargs())

        ax.axis[:].set_zorder(100)

        if self.show_colorbar:
            if cax is None:
                # add colorbar axes
                divider = make_axes_locatable(ax)
                self.cax = cax = divider.new_horizontal("5%", pad="2%", axes_class=Axes)
                fig.add_axes(cax)
                self.cbar = P.colorbar(im, cax=cax)
            else:
                # See comment for ROISmoothedSources's colobar code.
                self.cax = cax
                self.cbar = cax.colorbar(im)

        if self.title is not None:
            ax.set_title(self.title)

        ROIMapPlotter.overlay_region(self.roi,ax,self.header, **self.overlay_kwargs)

        self._additional_plotting()
        if self.extra_overlay is not None: self.extra_overlay(ax)

        tight_layout(fig)

        if filename is not None: P.savefig(filename)

    def imshow_kwargs(self):
        return dict()

    def _additional_plotting(self):
        pass

    @staticmethod
    def overlay_sources(roi, ax, 
                        exclude_sources=[], # don't plot these sources
                        override_kwargs = {}, # override kwargs for particular sources
                        **kwargs):
        """ Overlay sources in the ROI.

            exclude_sources is a list of sources or names of sources to not to overlay.

            override_kwargs is a dictionary where the keys are source
            names and the values are dictionaries of kwargs to pass
            into the overlay_source function. Useful if you want to
            format certain source markers specially.
        """
        for source in roi.get_sources():
            if source.name not in exclude_sources and source not in exclude_sources:
                pass_kwargs = kwargs.copy()

                if override_kwargs.has_key(source.name):
                    pass_kwargs.update(override_kwargs[source.name])
                elif override_kwargs.has_key(source):
                    pass_kwargs.update(override_kwargs[source])

                ROIMapPlotter.overlay_source(source, ax, **pass_kwargs)

    @staticmethod
    def overlay_source(source, ax, 
                       white_edge=True,
                       label_sources=False, 
                       text_color='black',
                       text_white_edge=False,
                       **kwargs
                      ):

        l = source.skydir.l()
        b = source.skydir.b()

        all_kwargs = dict(marker='x',color='black', markersize=12, zorder=2)
        all_kwargs.update(kwargs)
        
        # plot sources
        if white_edge and all_kwargs['marker'] in ['+', 'x']:
            white_kwargs = all_kwargs.copy()
            white_kwargs['color']='white'
            white_kwargs['markersize']=all_kwargs['markersize']+1
            white_kwargs['markeredgewidth']=2
            white_kwargs['zorder']=all_kwargs['zorder']-0.1

            ax["gal"].plot([l],[b],**white_kwargs)
        else:
            if white_edge: 
                all_kwargs['markeredgecolor'] = 'white'

        ax["gal"].plot([l],[b],**all_kwargs)

        if label_sources: 
            txt=ax["gal"].annotate(source.name, (l,b), 
                                   ha='center', va='top',
                                   color=text_color,
                                   xytext=(0,1.5*all_kwargs['markersize']), textcoords='offset points')
            if text_white_edge: 
                set_path_effects(txt,foreground="w", linewidth=2)

    @staticmethod
    def overlay_extensions(roi, **kwargs):
        for source in roi.get_extended_sources():
            ROIMapPlotter.overlay_extension(source, **kwargs)

    @staticmethod
    def overlay_extension(source, axes, header, white_edge=True, extension_color='black',
                          extension_zorder=None):
        import pyregion
        sm=source.spatial_model

        if isinstance(sm,PseudoSpatialModel) or isinstance(sm,SpatialMap):
            return

        region_string='\n'.join(region_writer.unparse_extension(sm,extension_color=extension_color))
        reg = pyregion.parse(region_string).as_imagecoord(header)

        # artist_list doesn't do anything
        patch_list, artist_list = reg.get_mpl_patches_texts()
        for p in patch_list: 
            if white_edge: set_path_effects(p,foreground="w", linewidth=2)
            if extension_zorder is not None:
                p.set_zorder(extension_zorder)
            axes.add_patch(p)

    @staticmethod
    def overlay_region(roi, ax, header, 
                       white_edge=True, # put a white edge around stuff
                       show_sources=True,
                       show_extensions=True, 
                       extension_color='black',
                       extension_zorder=None,
                       **kwargs):

        if show_extensions:
            # plot extended sources first so markers show up on top
            try:
                ROIMapPlotter.overlay_extensions(roi, axes=ax, header=header, 
                                                    white_edge=white_edge, 
                                                    extension_color=extension_color,
                                                    extension_zorder=extension_zorder)
            except ImportError as er:
                print ("To add extension information to plots, must install modules pyregion/pyparsing.")


        if show_sources:
            ROIMapPlotter.overlay_sources(roi, ax, 
                                              white_edge=white_edge, 
                                              **kwargs)



class ROITSMapPlotter(ROIMapPlotter):
    """ Create a residual TS map plot and overlay
        on it all of the sources in the ROI."""

    defaults = ROIMapPlotter.defaults + (
        ('pixelsize',               0.125, ),
    )

    def imshow_kwargs(self):
        # since this is a residual tsmap, never let the scale go below 25.
        norm=colors.Normalize(vmin=0, vmax=25) if np.max(self.pf['PRIMARY'].data) < 25 else None
        cmap=colormaps.b
        return dict(norm=norm, cmap=cmap)
    
    def create_pyfits(self):

        image=ROITSMapImage(self.roi,
                center=self.roi.roi_dir,
                pixelsize=self.pixelsize,
                size=self.size,
                galactic=self.galactic,
        )
        pf=image.get_pyfits()

        return pf

class ROISignificance(ROIMapPlotter):
    """ Make a plot of the poisson significance for the observed
        counts within a circual aperature of radius . """

    defaults = ROIMapPlotter.defaults + (
        ('kernel_rad',       0.25, 'Sum counts/model within radius degrees.'),
        ('conv_type',          -1,                        'Conversion type'),
    )
    defaults = keyword_options.change_defaults(defaults,'title','Poisson Significance')

    def create_pyfits(self):

        # Fit many pixels inside of the summing radius
        pixelsize=self.kernel_rad/10.0

        kwargs=dict(size=self.size,
                    pixelsize=pixelsize,
                    galactic=self.galactic,
                    conv_type=self.conv_type,
                    kerneltype='tophat',
                    kernel_rad=self.kernel_rad)

        counts=SmoothedCounts(self.roi,**kwargs)
        model=SmoothedModel(self.roi,**kwargs)


        pyfits = counts.get_pyfits()
        pyfits['PRIMARY'].data = ROISignificance.poisson_sigma(counts.image,model.image)

        return pyfits

    @staticmethod
    def poisson_sigma(obs,pred):
        """ Compute the sigma of the detection if you observe 
            a given number (obs) of counts and predict a 
            given number (pred) of counts.

            Note, it is numerically easier to deal with numbers close to 0
            then close to 1, so compute the sigma from the cdf or the sf
            depending upon which is smaller. """
        cdf=poisson.cdf(obs,pred)
        sf=poisson.sf(obs,pred)

        return np.where(cdf < 0.5, norm.ppf(cdf), norm.isf(sf))




class ROISmoothedSources(ROIMapPlotter):
    """ Make a smoothed residual plot which subtracts the diffuse emission. """



    defaults = ROIMapPlotter.defaults + (
        ('which',             None,    'Draw the smoothed point version of this source.'),
        ('conv_type',           -1,                                    'Conversion type'),
        ('overlay_psf',       True, 'Add a smoothed reference PSF on top of the counts.'),
        ('label_psf',         False,                        'Add a label on the PSF box.'),
        ('psf_size',             1,                         'Size of the PSF insert box'), 
        ('psf_loc',              4,                       'Location to put the psf box.'), 
        ('kerneltype',  'gaussian',                'Type of kernel to smooth image with'),
        ('kernel_rad',         0.1,            'Sum counts/model within radius degrees.'),
        ('override_center',   None,                               'Pick a better center'),
        ('cmap',              None,                                  'Show the colorbar'),
        ('colorbar_radius',   None,  """ If specified, calculate the intensity maximum
                                             pixel to be the maximum of all pixels within
                                             this radius. This is useful if you want to clip
                                             out a very bright nearby sources. By default,
                                             this radius will be the source size or the PSF
                                             size (if which!=None). If which==None, the
                                             default is that colorbar_radius=inf.        """),
        ('pixelsize_fraction', 5.0, 'The pixelsize is kernel_rad/pizelsize_fraction'),
    )

    defaults = keyword_options.change_defaults(defaults,'title','Smoothed Counts Map')

    def get_residual(self,**kwargs):
        """ Allow the particular method for getting the residual image to be overloaded. """

        residual = SmoothedResidual(self.roi,
                override_diffuse_sources=[i for i in self.roi.dsm.diffuse_sources if not hasattr(i,'skydir')],
                **kwargs)

        return residual

    def create_pyfits(self):

        # Fit many pixels inside of the summing radius
        self.pixelsize=self.kernel_rad/self.pixelsize_fraction

        roi = self.roi

        if self.which is not None:
            self.source = roi.get_source(self.which)
            self.center = self.source.skydir 
        else:
            self.source = None
            self.center = roi.roi_dir

        if self.override_center is not None: self.center = self.override_center

        self.smoothed_kwargs=dict(size=self.size,
                             pixelsize=self.pixelsize,
                             galactic=self.galactic,
                             conv_type=self.conv_type,
                             center=self.center,
                             per_solid_angle=True,
                             kerneltype=self.kerneltype,
                             kernel_rad=self.kernel_rad)

        residual = self.get_residual(**self.smoothed_kwargs)

        return residual.get_pyfits()

    def imshow_kwargs(self):
        self.max_intensity = ROISmoothedSources.get_max_intensity(self.source,
                                                                  self.pf,self.roi, 
                                                                  colorbar_radius=self.colorbar_radius)
        # sanity check
        self.max_intensity = max(self.max_intensity,1)

        if self.cmap is None: 
            self.cmap = colormaps.b 

        imshow_kwargs = dict(cmap=self.cmap, vmin=0, vmax=self.max_intensity)
        return imshow_kwargs 



    @staticmethod
    def get_max_intensity(source, pyfits, roi, colorbar_radius=None):
        """ Return the maximum value in the pyfits file either 
            within the extended source's size or otherwise
            within the 68% containment radius of the PSF (at
            the lowest energy). """
        try:
            import pyregion
        except:
            return pyfits[0].data.max()

        if source is None and colorbar_radius is None:
            return pyfits[0].data.max()

        if colorbar_radius is not None:
            ra,dec=roi.roi_dir.ra(),roi.roi_dir.dec()
            reg = pyregion.parse("fk5; circle(%.4f, %.4f, %.4f)" % (ra,dec,colorbar_radius))
            extensionmask = reg.get_mask(pyfits[0])

        elif hasattr(source,'spatial_model'):
            # For extended sources,
            # Get the maximum intensity value inside
            # the spatial model's extension
            extension_string='\n'.join(region_writer.unparse_extension(source.spatial_model,r68=True))
            reg = pyregion.parse(extension_string)
            extensionmask = reg.get_mask(pyfits[0])
        else:
            extensionmask = False # no mask

        # Get the maximum intensity inside the PSF (in lowest bin)
        emin=roi.bin_edges[0]
        ra,dec=source.skydir.ra(),source.skydir.dec()
        r68=roi.sa.psf.inverse_integral(emin,1,68) # 1=front entering events
        reg = pyregion.parse("fk5; circle(%.4f, %.4f, %.4f)" % (ra,dec,r68))
        psfmask = reg.get_mask(pyfits[0])

        # use whatever region mask is bigger
        return pyfits[0].data[psfmask | extensionmask].max()

    def _additional_plotting(self):
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        import pywcsgrid2

        if self.overlay_psf:

            if self.source is not None:
                model = self.source.model.copy()
            else:
                model = PowerLaw(index=2)

            # convert it to a point source placed at the origin
            point_version=PointSource(name='PSF',
                                      skydir=self.center,
                                      model=model)

            # create an image of the PSF (for our model).
            # Shrink down the image of the psf
            psf_kwargs=copy.deepcopy(self.smoothed_kwargs)
            psf_kwargs['size']=self.psf_size
            self.psf_model=SmoothedModel(self.roi,
                    override_point_sources=[point_version],
                    **psf_kwargs)
            self.psf_pyfits = self.psf_model.get_pyfits()

            # Normalize psf to have same maximum pixel scale
            # as residual image.
            self.psf_pyfits[0].data *= self.max_intensity/np.max(self.psf_pyfits[0].data)

            h_psf, d_psf = self.psf_pyfits[0].header, self.psf_pyfits[0].data
            self.axins = axins = zoomed_inset_axes(self.axes, zoom=1, loc=self.psf_loc,
                              axes_class=pywcsgrid2.Axes,
                              axes_kwargs=dict(wcs=h_psf))

            # Note, match color maps with parent.
            axins.imshow(d_psf, origin="lower", **self.imshow_kwargs())
            axins.axis[:].set_zorder(100)
            axins.axis[:].toggle(all=False)
            axins.axis[:].line.set_color('white')

            if self.label_psf:
                axins.add_inner_title("PSF", loc=3)

        self.cax.set_ylabel(r'$\mathrm{counts}/[\mathrm{deg}]^2$')



class ROISmoothedSource(ROISmoothedSources):
    """ Subclass ROISmoothedSources, but also subtract all background sources. """


    def get_residual(self,**kwargs):

        if self.source is None: raise Exception("Unable to subtract background sources unless which points at a real source.")

        self.roi.zero_source(which=self.source)

        residual = SmoothedResidual(self.roi,**kwargs)

        self.roi.unzero_source(which=self.source)

        return residual

class ROISmoothedResidual(ROISmoothedSources):
    """ Subclass ROISmoothedSources, but subtract all sources. """
    def get_residual(self,**kwargs):
        return SmoothedResidual(self.roi,**kwargs)



class ROISmoothedModel(ROISmoothedSources):
    """ Subclass ROISmoothedSources, but shows the model counts for all point+extended sources.
        This should be directly comparable to the ROISmoothedSources. """
    def get_residual(self,**kwargs):

        # Model counts for non-background sources.
        model = SmoothedModel(self.roi,
                              override_point_sources=self.roi.psm.point_sources,
                              override_diffuse_sources=[i for i in self.roi.dsm.diffuse_sources if hasattr(i,'skydir')],
                              **kwargs)

        return model





class ROIDisplay(object):
    """ Plotting a two dimensional map of the counts, model predicted,
        and residual counts in the ROI. Also plot a histogram of
        the weighted residuals and the p-values. """

    defaults = (
            ('figsize',       (7,6.5),                    'Size of the image'),
            ('fignum',           None,               'matplotlib figure number'),
            ('pixelsize',        0.25,               'size of each image pixel'),
            ('conv_type',          -1,                        'Conversion type'),
            ('size',               10,              'Size of the field of view'),
            ('nticks',              5,              'Number of axes tick marks'),
            ('galactic',         True,             'Coordinate system for plot'),
            ('countsfile',       None, 'Fits file to save the counts map data.'),
            ('modelfile',        None,  'Fits file to save the model map data.'),
            ('extra_overlay',    None, 'Function which can be used to overlay stuff on the plot.'),
            ('overlay_kwargs', dict(), 'kwargs passed into overlay_region'),
            ('title',            None),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        try:
            import pywcsgrid2
        except:
            raise Exception("You must install pywcsgrid2 to use this module.")
        
        self.roi = roi

        self.cm = CountsImage(self.roi,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)
        self.mm = ModelImage(self.roi,size=self.size,pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        self.cm_pf, self.mm_pf = self.cm.get_pyfits(), self.mm.get_pyfits()

        self.cm_d, self.mm_d = self.cm_pf[0].data, self.mm_pf[0].data
        self.res_d=(self.cm_d-self.mm_d)/self.mm_d**0.5

        self.h=self.cm_pf[0].header

        if self.countsfile is not None: 
            self.cm_pf.writeto(self.countsfile,clobber=True)
        if self.modelfile is not None: 
            self.mm_pf.writeto(self.modelfile,clobber=True)

        # Use same scale for the counts and model map
        self.counts_max = np.ceil(max(np.max(self.cm.image),np.max(self.mm.image)))
        self.counts_min = max(np.floor(min(np.min(self.cm.image),np.min(self.mm.image))),1.0)

    def add_cbar(self,im,ax):
        from matplotlib.axes import Axes
        import pylab as P
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal("5%", pad="2%", axes_class=Axes)
        self.fig.add_axes(cax)
        cbar = P.colorbar(im, cax=cax)

    def model_plot(self):

        from uw.utilities import colormaps
        im=self.ax_model.imshow(self.mm_d, norm=self.counts_norm, 
                            cmap=colormaps.b,
                            **self.imshow_args)

        #self.grid[0].cax.colorbar(im)
        self.ax_model.set_title('Model Counts')
        self.add_cbar(im,self.ax_model)
        
    def counts_plot(self):

        from uw.utilities import colormaps
        d=self.cm_d
        m=self.counts_min
        im=self.ax_counts.imshow(np.where(d>m,d,m), norm=self.counts_norm,
                            cmap=colormaps.b,
                            **self.imshow_args)
        self.ax_counts.set_title('Observed Counts')
        self.add_cbar(im,self.ax_counts)

    def resids_plot(self):

        im=self.ax_res.imshow(self.res_d,
                            norm=self.norm_res,
                            **self.imshow_args)
        self.ax_res.set_title('Weighted Residuals')
        self.add_cbar(im,self.ax_res)

    def hist_plot(self):
        self.ax_pvals.set_title('P-Values')

        mc = self.mm_d.flatten()
        cc = self.cm_d.flatten()
        rc = self.res_d.flatten()

        pvals = 1-poisson.cdf(cc,mc)

        nb = 20
        av = float(len(pvals)) / nb
        self.ax_pvals.hist(pvals,bins=np.linspace(0,1,20),histtype='step')
        self.ax_pvals.axhline( av, color='red')  
        lo,hi = ROIDisplay.ppf( (50.-95/2)/100., av), ROIDisplay.ppf( (50. + 95/2)/100.,av)
        self.ax_pvals.axhspan(lo , hi, facecolor='red', alpha=0.3)

        from matplotlib.offsetbox import AnchoredText
        from matplotlib.font_manager import FontProperties

        font=dict(size='small');
        at=AnchoredText('$95\%$ Conf.', loc=1, prop=font,frameon=False)
        self.ax_pvals.add_artist(at)
        set_path_effects(at.txt._text,foreground="w", linewidth=3)

        self.ax_resplot.set_title('Weighted Residuals')

        bins=np.linspace(-5,5,20)
        dx=bins[1]-bins[0]

        try:
            self.ax_resplot.hist(rc, bins=bins, histtype='step')
        except ValueError:
            pass
        b=np.linspace(-5,5,100)
        # overlay gaussian with same normalization
        self.ax_resplot.plot(b,(len(mc)*dx)*norm.pdf(b))
        self.ax_resplot.axvline(0,color='red')
        self.ax_resplot.set_xbound(lower=-5,upper=5)

    def show(self,filename=None):

        # taken from http://matplotlib.sourceforge.net/users/usetex.html
        import pywcsgrid2
        from matplotlib.gridspec import GridSpec

        self.imshow_args = dict(interpolation='nearest', origin='lower')


        self.counts_norm = colors.LogNorm(vmin=self.counts_min,vmax=self.counts_max)
        self.norm_res = colors.Normalize(vmin=-5,vmax=5)

        # first, divide the plot in 4x4
        self.fig = P.figure(self.fignum,self.figsize)
        P.clf()
        gs = GridSpec(4, 4)
        gs.update(wspace=1.5, hspace=1)
        self.ax_model = pywcsgrid2.subplot(gs[0:2, 0:2], header=self.h)
        self.ax_counts = pywcsgrid2.subplot(gs[0:2,2:4], header=self.h)
        self.ax_res = pywcsgrid2.subplot(gs[2:4, 0:2], header=self.h)

        # then divide the 4th space in two.
        self.ax_pvals = P.subplot(gs[2, 2:4])
        self.ax_pvals.yaxis.set_major_locator(MaxNLocator(3))

        self.ax_resplot = P.subplot(gs[3, 2:4])
        self.ax_resplot.yaxis.set_major_locator(MaxNLocator(3))

        self.model_plot()
        self.counts_plot()
        self.resids_plot()
        self.hist_plot()


        for ax in [self.ax_model, self.ax_counts, self.ax_res]:
            ax.axis[:].set_zorder(100)
            ROIMapPlotter.overlay_region(self.roi,ax,self.h, **self.overlay_kwargs)
            if self.extra_overlay is not None: self.extra_overlay(ax)

        if self.title is not None: self.fig.suptitle(self.title)
        if filename is not None: P.savefig(filename)


    @staticmethod
    def ppf(prob,mean):
        """Return the (approximate) Poisson percentage point function for given distribution.  Klugey."""
        if mean > 100: #normal approximation
            n = norm(mean,mean**0.5)
            return n.ppf(prob)        
        d = poisson(mean)
        prev = 0
        for i in range(1,200):        
            new = d.cdf(i)
            if new >= prob: break
            prev = new
        return (i-1)+(prob-prev)/(new-prev) #linear interpolation





class ROISmoothedDataModel(object):
    """ Plot (on the left) the diffuse subtracted smoothed counts and
        (on the right) the diffuse subtrcted smoothed model predicted
        counts. Useful to see if your model (qualitativly) looks like
        the right source. """

    defaults=ROISmoothedSources.defaults
    defaults = keyword_options.change_defaults(defaults,'figsize',(8,3.5))
    defaults = keyword_options.change_defaults(defaults,'title',None)

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi

        self.cmap = colormaps.b

        # Fit many pixels inside of the summing radius
        self.pixelsize=self.kernel_rad/5.0

        # Background subtracted counts
        counts_kwargs=keyword_options.defaults_to_kwargs(self,ROISmoothedSources)
        counts_kwargs['show_colorbar']=False
        counts_kwargs['overlay_psf']=False
        self.counts=ROISmoothedSources(roi,**counts_kwargs)

        # Model counts for non-background sources.
        model_kwargs=keyword_options.defaults_to_kwargs(self,ROISmoothedSources)
        model_kwargs['overlay_psf']=False
        self.model=ROISmoothedModel(roi,**model_kwargs)

    def show(self,filename=None):

        from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
        import pywcsgrid2

        self.fig = fig = P.figure(self.fignum,self.figsize)
        P.clf()

        header = self.counts.pf[0].header

        self.grid = grid = ImageGrid(fig, (1, 1, 1), nrows_ncols = (1, 2),
                              axes_pad=0.1, share_all=True,
                              cbar_mode="single", cbar_pad="2%",
                              cbar_location="right",
                              axes_class=(pywcsgrid2.Axes, 
                                          dict(header=header)))

        # This is kind of ugly, but use show() once to set the max_intensity, then show again
        # with matching colormaps. Kluge, but in an unimportant function.
        def show():
            self.counts.show(axes=self.grid[0])
            self.model.show(axes=self.grid[1],cax=grid.cbar_axes[0])
        
        show()

        # Same colorscale kluge as in ROISmoothedBeforeAfter
        max_intensity = max(self.counts.max_intensity, self.model.max_intensity)
        self.counts.max_intensity = self.model.max_intensity  = max_intensity
        show()

        self.grid[0].add_inner_title("Counts", loc=2)
        self.grid[1].add_inner_title("Model", loc=2)


        tight_layout(self.fig)

        self.header = self.h = [self.counts.header,self.model.header]

        if filename is not None: P.savefig(filename)



class ROISmoothedBeforeAfter(object):
    """ Make a 2x1 plot where the left plot
        has the diffuse subtracted and the right plot has the diffuse
        and all sources subtracted.

        This plot is nice for seeing how well the background source
        subtracting is doing. """
    defaults=ROISmoothedSources.defaults

    defaults = keyword_options.change_defaults(defaults,'figsize',(8,3.5))
    defaults = keyword_options.change_defaults(defaults,'title',None)

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):

        if kwargs.has_key('show_colorbar') or \
           kwargs.has_key('overlay_psf'):
            raise Exception("This feature doesn't work, for now...")

        keyword_options.process(self, kwargs)

        self.roi = roi

        sources_kwargs=keyword_options.defaults_to_kwargs(self,ROISmoothedSources)
        sources_kwargs['show_colorbar']=False
        self.smoothed_sources=ROISmoothedSources(roi,**sources_kwargs)

        source_kwargs=keyword_options.defaults_to_kwargs(self,ROISmoothedSources)
        source_kwargs['overlay_psf']=False
        self.smoothed_source=ROISmoothedSource(roi,**source_kwargs)

    def show(self,filename=None):

        from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
        import pywcsgrid2

        self.fig = fig = P.figure(self.fignum,self.figsize)
        P.clf()

        header = self.smoothed_source.pf[0].header

        self.grid = grid = ImageGrid(fig, (1, 1, 1), nrows_ncols = (1, 2),
                              axes_pad=0.1, share_all=True,
                              cbar_mode="single", cbar_pad="2%",
                              cbar_location="right",
                              axes_class=(pywcsgrid2.Axes, 
                                          dict(header=header)))

        def show():
            self.smoothed_sources.show(axes=self.grid[0])
            self.smoothed_source.show(axes=self.grid[1],cax=grid.cbar_axes[0])
        
        show()

        # use same color scale for each plot
        max_intensity = max(self.smoothed_sources.max_intensity, self.smoothed_source.max_intensity)
        self.smoothed_sources.max_intensity = self.smoothed_source.max_intensity  = max_intensity

        show()

        tight_layout(self.fig)

        self.header = self.h = [self.smoothed_sources.header,self.smoothed_source.header]

        if filename is not None: P.savefig(filename)
