"""
Plotting routines to display results of an ROI analysis.

     # Plot the counts in a vertical slice
     ROISlice(roi).show()

     # Plot the integral counts in a radius.
     ROIRadialIntegral(roi).show()


$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_plotting.py,v 1.93 2012/07/12 20:03:04 lande Exp $

author: Matthew Kerr, Joshua Lande
"""
import pprint
import numpy as np
from . roi_image import ModelImage,CountsImage,RadialCounts,RadialModel
from . roi_extended import ExtendedSource
from . pointspec_helpers import PointSource
from uw.utilities import keyword_options

import pylab as P
from matplotlib import ticker
from matplotlib.ticker import FuncFormatter


def tight_layout(fig):
    import traceback
    import sys
    try:
        fig.tight_layout()
    except:
        print ('To make pretty figures, you need a newer version of matplotlib which has tight_layout!')
        traceback.print_exc(file=sys.stdout)

def format_degrees(x,*args):
    r='%g' % x
    if '.' in r:
        return ('$%s$' % r).replace('.',r'.\!\!^\circ')
    else:
        return '$%s^\circ$' % r

DegreesFormatter=FuncFormatter(format_degrees)

        

class ROISlice(object):
    """ Create counts slice plot (integrating out one dimention and plotting
        the counts and model predicted counts in the other direction. """

    defaults = (
            ('which',           None,                       'Source to analyze'),
            ('figsize',        (7,6),                       'Size of the image'),
            ('fignum',          None,                'matplotlib figure number'),
            ('pixelsize',      0.125,                'size of each image pixel'),
            ('size',              10,               'Size of the field of view'),
            ('galactic',        True,              'Coordinate system for plot'),
            ('int_width',          2,            'Integration width for slice.'),
            ('conv_type',         -1,                         'Conversion type'),
            ('just_diffuse',    True, """Display the model predictions with 
                                               all point + extended sources 
                                             removed. The background is not 
                                                                     refit. """),
            ('aspoint',         True, """Display also the model predictions 
                                            for an extended source fit with 
                                           the point hypothesis. Only works 
                                          when which is an extended source. """),
            ('oversample_factor',  4, """ Calculate the model predictions 
                                          this many times more finely.  
                                          This will create a smoother 
                                          plot of model predictions. Set 
                                          to 1 if you want the model 
                                          predictions to 'look like' the 
                                          data.                             """),
            ('title',           None,                      'Title for the plot'),
            ('black_and_white',False, """ If True, make the plot black and 
                                          white (better for printing and
                                          publishing)                       """),
            ('legend',         True,  """ Show legend. """),
    )

    @staticmethod
    def set_color_cycle():
        P.gca().set_color_cycle(['k', 'b', 'g', 'r', 'm', 'y', 'k'])


    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        if not type(self.oversample_factor) == int:
            raise Exception("oversample_factor must be an integer.")

        self.roi   = roi

        manager,index=self.roi.mapper(self.which)
        if manager == roi.psm:
            self.source=manager.point_sources[index]
            self.pretty_name = 'Point'
        else:
            self.source=manager.diffuse_sources[index]
            self.pretty_name = self.source.spatial_model.pretty_name

        self.center=self.source.skydir

        self.get_counts()
        self.get_model()

    def get_counts(self):
        self.ci_x = CountsImage(self.roi,center=self.center,size=(self.size,self.int_width),
                                pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        self.counts_x=self.ci_x.image.sum(axis=0)
        self.counts_dx=ROISlice.range(self.counts_x,self.pixelsize,x_axis=True)

        if np.any(np.isnan(self.counts_x)):
            raise Exception("Error in calculating a slice, the x counts slice contains NaNs.")

        self.ci_y = CountsImage(self.roi,center=self.center,size=(self.int_width,self.size),
                                pixelsize=self.pixelsize,galactic=self.galactic,conv_type=self.conv_type)

        self.counts_y=self.ci_y.image.sum(axis=1)
        self.counts_dy=ROISlice.range(self.counts_y,self.pixelsize,x_axis=False)

        if np.any(np.isnan(self.counts_y)):
            raise Exception("Error in calculating a slice, the y counts slice contains NaNs.")

    @staticmethod
    def cache_roi(roi):
        roi.old_parameters = roi.parameters().copy() # save free parameters

    @staticmethod
    def uncache_roi(roi):
        roi.set_parameters(roi.old_parameters) # reset free parameters
        roi.__update_state__()

    def get_model(self):
        model_pixelsize=float(self.pixelsize)/self.oversample_factor

        kwargs=dict(center=self.center,pixelsize=model_pixelsize,
                    galactic=self.galactic,conv_type=self.conv_type)

        self.names = []
        self.mi_x = []
        self.mi_y = []

        self.names.append(self.pretty_name)
        self.mi_x.append(ModelImage(self.roi,size=(self.size,self.int_width),**kwargs))
        self.mi_y.append(ModelImage(self.roi,size=(self.int_width,self.size),**kwargs))

        if self.aspoint and isinstance(self.source,ExtendedSource):
            # replace with a point source

            ROISlice.cache_roi(self.roi)

            es=self.roi.get_source(self.which)
            
            ps=PointSource(name=es.name,model=es.model.copy(),skydir=es.spatial_model.center,leave_parameters=True)
            self.roi.add_source(ps)

            # only zero it after making a copy of the spectral part!
            self.roi.zero_source(es)

            self.roi.fit(estimate_errors=False)

            self.names.append('Point')
            self.mi_x.append(ModelImage(self.roi,size=(self.size,self.int_width),**kwargs))
            self.mi_y.append(ModelImage(self.roi,size=(self.int_width,self.size),**kwargs))

            self.roi.del_source(ps)
            self.roi.unzero_source(es)

            ROISlice.uncache_roi(self.roi)

        if self.just_diffuse:
            # hide all point + extended sources.

            sources = list(self.roi.psm.point_sources) + \
                    [ i for i in self.roi.dsm.diffuse_sources if isinstance(i,ExtendedSource) ]
            # don't zero already zeroed sources
            sources = [ i for i in sources if not i.model.iszero() ]

            ROISlice.cache_roi(self.roi)

            for source in sources: self.roi.zero_source(source)

            self.names.append('Diffuse')
            self.mi_x.append(ModelImage(self.roi,size=(self.size,self.int_width),**kwargs))
            self.mi_y.append(ModelImage(self.roi,size=(self.int_width,self.size),**kwargs))

            for source in sources: self.roi.unzero_source(source)

            ROISlice.uncache_roi(self.roi)

        self.models_x=[model.image.sum(axis=0)*self.oversample_factor for model in self.mi_x]
        self.models_y=[model.image.sum(axis=1)*self.oversample_factor for model in self.mi_y]

        self.model_dx = ROISlice.range(self.models_x[0],model_pixelsize,x_axis=True)
        self.model_dy = ROISlice.range(self.models_y[0],model_pixelsize,x_axis=False)

    @staticmethod
    def range(data,pixelsize,x_axis=False):
        if x_axis:
            return (len(data)/2.0-np.arange(len(data)))*pixelsize - pixelsize/2
        else:
            return (np.arange(len(data))-len(data)/2.0)*pixelsize + pixelsize/2

    @staticmethod
    def get_styles(black_and_white):
        return [
            dict(color='k' if black_and_white \
                 else 'red',linestyle='-'),
            dict(color='k' if black_and_white \
                 else 'blue',dashes=[5,3]),
            dict(color='k' if black_and_white \
                 else 'g',dashes=[5,3,1,3]),
        ]

    def plotx(self,axes,legend=False):

        ax = axes

        ROISlice.set_color_cycle()

        styles=ROISlice.get_styles(self.black_and_white)[0:len(self.names)]
        for name,model,style in zip(self.names,self.models_x,styles):
            ax.plot(self.model_dx,model,label=name,**style)

        ax.errorbar(self.counts_dx,self.counts_x,yerr=np.sqrt(self.counts_x),
                    label='Counts', fmt='.', color='black')

        ax.set_xlim(self.counts_dx[0],self.counts_dx[-1])
        ax.set_ylim(ymin=0)

        if legend: ax.legend(loc='upper right',numpoints=1)

        ax.xaxis.set_major_formatter(DegreesFormatter)


        ax.set_xlabel(r'$\Delta l$' if self.galactic else r'$\Delta \mathrm{RA}$')
        ax.set_ylabel('Counts')

    def ploty(self,axes,legend=False):

        ax = axes

        ROISlice.set_color_cycle()

        styles=ROISlice.get_styles(self.black_and_white)[0:len(self.names)]
        for name,model,style in zip(self.names,self.models_y,styles):
            ax.plot(self.model_dy,model,label=name,**style)

        ax.errorbar(self.counts_dy,self.counts_y,yerr=np.sqrt(self.counts_y),
                    label='Counts', fmt='.', color='black')

        ax.set_xlim(self.counts_dy[0],self.counts_dy[-1])
        ax.set_ylim(ymin=0)

        if legend: ax.legend(loc='upper right',numpoints=1)

        ax.xaxis.set_major_formatter(DegreesFormatter)

        ax.set_xlabel(r'$\Delta b$' if self.galactic else r'$\Delta \mathrm{Dec}$')
        ax.set_ylabel('Counts')


        # only need to show legend once

    def save_data(self,datafile):
        """ Note, shrink model predicted counts to be same size as regular counts,
            for an easier file format. """

        x,y=['l','b'] if self.galactic else ['ra','dec']
            
        results_dict = {
            x : {
                'Counts': [self.counts_dx.tolist(), self.counts_x.tolist()]
            },
            y : {
                'Counts': [self.counts_dy.tolist(), self.counts_y.tolist()]
            }
        }

        for name,modelx,modely in zip(self.names,self.models_x,self.models_y):
            results_dict[x][name]=[self.model_dx.tolist(), modelx.tolist()]
            results_dict[y][name]=[self.model_dy.tolist(), modely.tolist()]

        file = open(datafile,'w')
        try:
            import yaml
            file.write(yaml.dump(results_dict))
        except:
            import pprint
            file.write(pprint.pformat(results_dict))

        file.close()

    def show(self,filename=None, datafile=None, ax1=None, ax2=None):

        if ax1 is None and ax2 is None:
            fig = P.figure(self.fignum,self.figsize)
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
        elif ax1 is not None or ax2 is not None:
            raise Exception("Both ax1 and ax2 must be specified.")
        else:
            fig = ax1.get_figure()

        self.ax1, self.ax2 = ax1, ax2

        self.plotx(ax1, legend=self.legend)
        self.ploty(ax2)

        if self.title is None:
            self.title = 'Counts Slice'
            self.title += ' for %s' % self.source.name

        fig.suptitle(self.title)
        tight_layout(fig)

        if datafile is not None: self.save_data(datafile)

        if filename is not None: P.savefig(filename)


class ROIRadialIntegral(object):
    """ Create a radial integral plot which integrates radially
        the counts and model predicted counts and bins uniformly in theta^2. """

    defaults = (
            ('which',           None,                       'Source to analyze'),
            ('size',               2,                'Size of image in degrees'), 
            ('pixelsize',     0.0625, """ size of each image pixel. This is a misleading because the
                                          size of each pixel varies, since the image is uniform in theta^2. 
                                          This value is used to determine the total number of pixels using
                                          the formula npix=size/pixelsize and represents something
                                          like an average pixelsize."""),
            ('npix',            None, """ If specified, use this value instead of pixelsize. """),
            ('conv_type',         -1, 'Conversion type'),
            ('just_diffuse',    True, """ Display the model predictions with all point + extended
                                          sources removed. The background is not refit. """),
            ('aspoint',         True, """ Display also the model predictions for an extended source 
                                          fit with the point hypothesis. Only works when which is an
                                          extended source. """),
            ('oversample_factor',  4, """ Calculate the model predictions this many times more finely. 
                                          This will create a smoother plot of model predictions. Set 
                                          to 1 if you want the model predictions to 'look like' the 
                                          data."""),
            ('legend',          True, """Add a legend to the plot."""),
            ('title',           None,   'Title for the plot'),
            ('black_and_white',False, """If True, make the plot black and white (better for 
                                         printing/publishing)"""),
            ('figsize',        (7,6),   'Size of the image'),
            ('fignum',          None,   'matplotlib figure number'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        if not type(self.oversample_factor) == int:
            raise Exception("oversample_factor must be an integer.")


        self.roi   = roi

        manager,index=self.roi.mapper(self.which)
        if manager == roi.psm:
            self.source=manager.point_sources[index]
            self.pretty_name = 'Point'
        else:
            self.source=manager.diffuse_sources[index]
            self.pretty_name = self.source.spatial_model.pretty_name

        self.center=self.source.skydir

        self.get_counts()
        self.get_model()

    def get_counts(self):
        self.ci = RadialCounts(self.roi,center=self.center,size=self.size,pixelsize=self.pixelsize,conv_type=self.conv_type, npix=self.npix)
        self.counts = self.ci.image

        self.theta_sqr_co=self.ci.bin_centers_deg

        if np.any(np.isnan(self.counts)):
            raise Exception("Error in calculating a radial integral plot, the counts contains NaNs.")


    def get_model(self):
        kwargs=dict(center=self.center,size=self.size,conv_type=self.conv_type, 
                    pixelsize=self.pixelsize/self.oversample_factor,
                    npix=self.npix*self.oversample_factor if self.npix is not None else None)

        # Use lists to preserve legend order
        self.mi,self.names = [],[]

        self.mi.append(RadialModel(self.roi,**kwargs))
        self.names.append(self.pretty_name)

        if self.aspoint and isinstance(self.source,ExtendedSource):

            ROISlice.cache_roi(self.roi)

            es=self.roi.get_source(self.which)

            ps=PointSource(name=es.name,model=es.model.copy(),skydir=es.spatial_model.center,leave_parameters=True)
            self.roi.add_source(ps)

            # only zero it after making a copy of the spectral part!
            self.roi.zero_source(es)

            self.roi.fit(estimate_errors=False)

            self.mi.append(RadialModel(self.roi,**kwargs))
            self.names.append('Point')

            self.roi.del_source(ps)
            self.roi.unzero_source(es)

            ROISlice.uncache_roi(self.roi)

        if self.just_diffuse:

            sources = list(self.roi.psm.point_sources) + \
                    [ i for i in self.roi.dsm.diffuse_sources if isinstance(i,ExtendedSource) ]
            # don't zero already zeroed sources
            sources = [ i for i in sources if not i.model.iszero() ]

            ROISlice.cache_roi(self.roi)

            for source in sources: self.roi.zero_source(source)

            self.mi.append(RadialModel(self.roi,**kwargs))
            self.names.append('Diffuse')

            for source in sources: self.roi.unzero_source(source)

            ROISlice.uncache_roi(self.roi)

        self.theta_sqr_mo=self.mi[0].bin_centers_deg

        # scale the model to line up with the counts
        for i in self.mi: i.image*=self.oversample_factor
        self.models=[i.image for i in self.mi]

        for name,model in zip(self.names,self.models):
            if np.any(np.isnan(model)):
                raise Exception("Error in calculating a radial integral, model %s contains NaNs." % name)

    def save_data(self,datafile):

        results_dict = {}

        results_dict['Counts']=[ self.theta_sqr_co.tolist(), self.counts.tolist() ]
        for name,model in zip(self.names,self.models):
            results_dict[name]=[ self.theta_sqr_mo.tolist(), model.tolist() ]

        file = open(datafile,'w')
        try:
            import yaml
            file.write(yaml.dump(results_dict))
        except:
            import pprint
            file.write(pprint.pformat(results_dict))
        file.close()

    def show(self,filename=None,datafile=None, axes=None):

        if axes is None:
            fig = P.figure(self.fignum,self.figsize)
            axes = fig.add_subplot(111)
        else:
            fig = ax.get_figure()

        self.axes = ax = axes

        ROISlice.set_color_cycle()

        styles=ROISlice.get_styles(self.black_and_white)[0:len(self.names)]
        for name,model,style in zip(self.names,self.models,styles):
            ax.plot(self.theta_sqr_mo,model,label=name,**style)

        ax.errorbar(self.theta_sqr_co,self.counts,yerr=np.sqrt(self.counts),
                    label='Counts', fmt='.', color='black')

        if self.legend:
            ax.legend(loc='upper right',numpoints=1)

        ax.set_xlabel(r'$\Delta \theta^2 ([\mathrm{deg}]^2)$')
        ax.set_ylabel('Counts')
        ax.set_ylim(ymin=0)

        ax.set_ylim(ymin=0)

        if self.title is None:
            self.title = 'Radially Integrated Counts'
            if self.source is not None: self.title += ' for %s' % self.source.name

        ax.set_title(self.title)
        tight_layout(fig)

        if datafile is not None: self.save_data(datafile)

        if filename is not None: P.savefig(filename)


