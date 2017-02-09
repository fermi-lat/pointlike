""" image processing:
    class AIT for full sky
          ZEA for square region
          TSplot special adapter for ZEA
          
     author: T. Burnett tburnett@u.washington.edu

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/image.py,v 1.46 2014/07/11 21:38:17 burnett Exp $

"""
version = '$Revision: 1.46 $'.split()[1]

import sys, pylab, types, os
import math
import numpy as np
import pylab as pl
import pylab as plt
from matplotlib import pyplot, ticker
import matplotlib as mpl
from skymaps import SkyImage, SkyDir, double2, SkyProj,PySkyFunction,Hep3Vector
from math import exp
from numpy.fft import fft2,ifft2,fftshift
from scipy import optimize
import keyword_options

SkyImage.setNaN(np.nan) 

class Ellipse(object):
    def __init__(self, q):
        """ q: ellipical parameters
               a, b, phi  
        """
        self.q = q 

    def contour(self, r=1, count=50):
        """ return set of points in around closed figure"""
        s,c = math.sin(-self.q[2]), math.cos(self.q[2])
        a,b = self.q[0],self.q[1]
        x = []
        y = []
        for t in np.linspace(0, 2*math.pi, count):
            ct,st = math.cos(t), math.sin(t)
            x.append( r*(a*ct*s - b*st*c))
            y.append( r*(a*ct*c + b*st*s))
        return x,y      


class Rescale(object):

    def __init__(self, image, nticks=5, galactic=False):
        """ image: a SkyImage object
            nticks: suggested number of ticks for the ticker
            
            warning: fails if the north pole is in the image (TODO: figure out a sensible approach)
        """

        # get ra range from top, dec range along center of SkyImage
        nx,ny = image.nx, image.ny
        self.nx=nx
        self.ny=ny
        # convenient lat, lon functions for pixel coords
        lat = lambda x,y: image.skydir(x,y).l() if galactic else image.skydir(x,y).ra()
        lon = lambda x,y: image.skydir(x,y).b() if galactic else image.skydir(x,y).dec()
        xl,xr = lat(0,0), lat(nx,0)
        if xl<xr: # did it span the boundary?
            xr = xr-360 
        self.vmin, self.vmax = lon(0,0), lon(nx/2.,ny)
        ticklocator = ticker.MaxNLocator(nticks, steps=[1,2,5])
        self.uticks = [ix if ix>-1e-6 else ix+360\
              for ix in ticklocator.bin_boundaries(xr,xl)[::-1]] #reverse
        self.ul = xl
        self.ur = xr
        self.vticks = ticklocator.bin_boundaries(self.vmin,self.vmax)
        if len(self.vticks)==0: # protect against rare situatin
            self.vticks = [self.vmin,self.vmax]

        # extract positions in image coords,  text labels
        self.xticks = [image.pixel(SkyDir(x,self.vmin,SkyDir.GALACTIC if galactic else SkyDir.EQUATORIAL))[0]\
                       for x in self.uticks]
        #self.yticks = [image.pixel(SkyDir(xl,v))[1] for v in self.vticks]
        # proportional is usually good?
        try:
            yscale = ny/(lon(0,ny)-self.vmin)
            self.yticks = [ (v-self.vmin)*yscale for v in self.vticks]
            self.xticklabels = self.formatter(self.uticks)
            self.yticklabels = self.formatter(self.vticks)
        except Exception, msg:
            print 'formatting failure in image.py: {}'.format(msg)
            self.xticks=self.yticks=None
    def formatter(self, t):
        n=0
        s = np.abs(np.array(t))+1e-9
        for i in range(4):
            #print s, s-np.floor(s), (s-np.floor(s)).max()
            if (s-np.floor(s)).max()<1e-3: break
            s = s*10
            n+=1
        fmt = '%%5.%df'%n
        return [(fmt% x).strip() for x in t]

    def apply(self, axes):
        if self.xticks is None:
            return
        #note remove outer ones
        if len(self.xticks)>=3:
            axes.set_xticks(self.xticks[1:-1])
            axes.set_xticklabels(self.xticklabels[1:-1])
        axes.xaxis.set_ticks_position('bottom')
        #axes.set_xlim((0.5,self.nx+0.5)) # have to do again?

        if len(self.yticks)>=3:
            axes.set_yticks(self.yticks[1:-1])
            axes.set_yticklabels(self.yticklabels[1:-1])

        axes.yaxis.set_ticks_position('left')
        #axes.set_ylim((0.5,self.ny+0.5)) # have to do again?


class AITproj(object):
    def __init__(self, proj):
        self.proj = proj
        self.center = proj(0,0)
        self.scale = ((proj(180,0)[0]-self.center[0])/180., (proj(0,90)[1]-self.center[1])/90)
    def __call__(self,l,b):
        r = self.proj(l,b)
        return [(r[i]-self.center[i])/self.scale[i] for i in range(2)]


def draw_grid(self, labels=True, color='gray', pixelsize=0.5, textsize=8):
        label_offset = 5/pixelsize
        
        #my_axes = pylab.axes() #creates figure and axes if not set

        # # pylab.matplotlib.interactive(False)
        #my_axes.set_autoscale_on(False)
        #my_axes.set_xlim(0, 360/pixelsize)
        #my_axes.set_ylim(0, 180/pixelsize)
        #my_axes.set_axis_off()
        #my_axes.set_aspect('equal')
        ait = AITproj(self.projector.sph2pix) # the projector to use
        axes = self.axes
        
        bs = np.arange(-90, 91, 5)
        for l in np.hstack((np.arange(0, 360, 45),[180.01])):
            lstyle = '-' if int(l)==180 or int(l)==0 else '--' 
            # axes.plot([ait(l,b)[0] for b in bs], [ait(l,b)[1] for b in bs], lstyle, color=color)
            axes.plot([ait(l,b)[0] for b in bs], [ait(l,b)[1] for b in bs], lstyle, color=color)
            if labels:
                x,y = ait(l, 45) 
                axes.text(x,y, '%3.0f'%l ,size=textsize, ha='center')

        ls = np.hstack((np.arange(180, 0, -5), np.arange(355, 180,-5), [180.01]))
        for b in np.arange(-60, 61, 30):
            lstyle = '-' if int(b)==0 else '--'
            axes.plot([ait(l,b)[0] for l in ls], [ait(l,b)[1] for l in ls], lstyle, color=color)
            if labels:
                x,y = ait(180.1, b)                                               
                axes.text(x+label_offset,y+b/60*label_offset, '%+3.0f'%b, size=textsize, ha='center',va='center')
        if labels:
            for b in [90,-90]:
                x,y = ait(0,b)
                axes.text(x,y+b/90*label_offset,'%+3.0f'%b, size=textsize, ha='center',va='center') 


class AIT_grid():

    def __init__(self, fignum=20, axes=None, labels=True, color='gray', pixelsize=0.5, textsize=8, linestyle='-'):
        """Draws gridlines and labels for map.

        """

        if axes is None:
            fig=plt.figure(fignum, figsize=(12,6))
            fig.clf()
            self.axes = fig.gca()
        else: self.axes = axes
        self.pixelsize = pixelsize

        xsize,ysize = 325,162
        crpix = double2(xsize/pixelsize/2., ysize/pixelsize/2.)
        crval = double2(0,0)
        cdelt = double2(-pixelsize, pixelsize)
        self.proj = SkyProj('AIT', crpix, crval, cdelt, 0, True)
        
 
        self.axes.set_autoscale_on(False)
        self.axes.set_xlim(0, 360/self.pixelsize)
        self.axes.set_ylim(0, 180/self.pixelsize)
        self.axes.set_axis_off()
        self.axes.set_aspect('equal')
        self.extent= (self.ait(180,0)[0],self.ait(180.001,0)[0], self.ait(0,-90)[1], self.ait(0,90)[1])
        label_offset = 5/self.pixelsize
        bs = np.arange(-90, 91, 5)
        for l in np.hstack((np.arange(0, 360, 45),[180.01])):
            self.axes.plot([self.ait(l,b)[0] for b in bs], [self.ait(l,b)[1] for b in bs], linestyle, color=color)
            if labels:
                x,y = self.ait(l, 45) 
                self.axes.text(x,y, '%3.0f'%l ,size=textsize, ha='center')

        ls = np.hstack((np.arange(180, 0, -5), np.arange(355, 180,-5), [180.01]))
        for b in np.arange(-60, 61, 30):
            lstyle = '-' if int(b)==0 else linestyle
            self.axes.plot([self.ait(l,b)[0] for l in ls], [self.ait(l,b)[1] for l in ls], lstyle, color=color)
            if labels:
                x,y = self.ait(180.1, b)                                               
                self.axes.text(x+label_offset,y+b/60*label_offset, '%+3.0f'%b, size=textsize, ha='center',va='center')#, weight = 'bold')
        if labels:
            for b in [90,-90]:
                x,y = self.ait(0,b)
                self.axes.text(x,y+b/90*label_offset,'%+3.0f'%b, size=textsize, ha='center',va='center') 

    def skydir(self, x, y):
        """ from pixel coordinates to sky """
        return SkyDir(x+0.5, y+0.5, self.proj) 

    def pixel(self, sdir):
        """ return pixel coordinates for the skydir
        """
        x,y=self.proj.sph2pix(sdir.l(),sdir.b())
        return  (x-0.5,y-0.5)

    def ait(self, l, b):
        " convert lon, lat to car "
        # floats seem to be necessary: gcc SWIG oddity?
        z = self.proj.sph2pix(float(l), float(b))
        return (float(z[0]),float(z[1]))

    def plot(self, sources, marker='o', text=None, fontsize=8, colorbar=False, **kwargs):
        """ plot symbols at points 
        sources: list of SkyDir
        
        keywords:
            text: optional text strings (same length as sources if specified)
            fontsize: for text
            marker: symbol to use, see scatter doc.
            colorbar: set True to add a colorbar. (See method to attach label)
            
        kwargs: applied to scatter, use c as an array of floats, with optional
                cmap=None, norm=None, vmin=None, vmax=None
                to make color key for another value
                in that case, you can use Axes.colorbar
                set s to change from default size =10
        """
        X=[]
        Y=[]
        for i,s in enumerate(sources):
            x,y = self.ait(s.l(),s.b())
            X.append(x)
            Y.append(y)
            if text is not None:
                self.axes.text(x,y,text[i],fontsize=fontsize)
        #self.axes.plot(X,Y, symbol,  **kwargs)
        self.col=self.axes.scatter(X,Y, marker=marker,  **kwargs)
        if colorbar: self.colorbar()
        plt.draw_if_interactive()
        
    def colorbar(self, label=None, **kwargs):
        """ attach a colorbar to a plot, expect it was generated with the 'c' argument
        """
        if 'shrink' not in kwargs: kwargs['shrink'] = 0.7
        fig = self.axes.figure
        if 'col' not in self.__dict__:
            raise Exception('colorbar called with no mappable defined, say by plot(..., c=...)')
        cb =fig.colorbar(self.col, cax=None, ax=self.axes, **kwargs)
        fig.sca(self.axes) # restore selected axes objet
        if label: cb.set_label(label)
        plt.draw_if_interactive()




class AIT(object):
    """ Manage a full-sky image of a SkyProjection or SkyFunction, wrapping SkyImage
     """
    
    defaults= (
        ('pixelsize', 0.5, 'size, in degrees, of pixels'),
        ('galactic',  True, 'galactic or equatorial coordinates'),
        ('fitsfile',  '',   'if set, write the projection to a FITS file'),
        ('proj',      'AIT', 'could be ''CAR'' for carree or ''ZEA'': used by wcslib'),
        ('center',    None,   'if default center at (0,0) in coord system'),
        ('size',      180,   'make less for restricted size'),
        ('galactic',  True,  'use galactic coordinates'),
        ('earth',     False, 'if looking down at Earth'),
        ('axes',      None,   'set to use, otherwise create figure if necessary'),
        ('nocolorbar',False,  'set to turn off colorbar' ),
        ('cbtext',    None,   'text for colorbar'),
        ('background',None,   'if set, a value to apply to NaN: default is to not set pixels\n' 
                                'nb: do not set 0 if log scale'),
        )
    
    @keyword_options.decorate(defaults)
    def __init__(self, skyfun, **kwargs):
        """
        skyfun SkyProjection or SkyFunction object
        """
        keyword_options.process(self, kwargs)
        self.skyfun = skyfun
         # set up, then create a SkyImage object to perform the projection to a grid
        if self.center is None:
            self.center = SkyDir(0,0, SkyDir.GALACTIC if self.galactic else SkyDir.EQUATORIAL)
        self.skyimage = SkyImage(self.center, self.fitsfile, self.pixelsize, 
            self.size, 1, self.proj, self.galactic, self.earth)
        # we want access to the projection object, to allow interactive display via pix2sph function
        self.projector = self.skyimage.projector()
        self.x = self.y = 100 # initial def
        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()
        self.center = self.projector.sph2pix(0,0)
        self.scale = ((self.projector.sph2pix(180,0)[0]-self.center[0])/180., 
                      (self.projector.sph2pix( 0,90)[1]-self.center[1])/90)

        if skyfun is not None: 
            self.fill(skyfun)
            self.setup_image(self.earth)
        else:
            # special case: want to set pixels by hand
            pass
            
    def fill(self, skyfun):
        """ fill the image with the skyfunction"""
        if skyfun.__class__.__name__ !='PySkyFunction':
            def pyskyfun(v):
                return skyfun(SkyDir(Hep3Vector(v[0],v[1],v[2])))
            self.skyimage.fill(PySkyFunction(pyskyfun))
        else:
            self.skyimage.fill(skyfun)

            
    def setup_image(self, earth=False):
        # now extract stuff for the pylab image, creating a masked array to deal with the NaN values
        self.image = np.array(self.skyimage.image()).reshape((self.ny, self.nx))
        self.mask = np.isnan(self.image)
        if self.background is None:
            self.masked_image = np.ma.array( self.image, mask=self.mask)
        else:
            self.masked_image = np.ma.array( self.image)
            self.masked_image[self.mask]=self.background
        size = self.size
        if not earth:
            self.extent = (180,-180, -90, 90) if size==180 else (size, -size, -size, size)
        else:
            self.extent = (-180,180, -90, 90) if size==180 else (-size, size, -size, size)
       
        self.vmin ,self.vmax = self.skyimage.minimum(), self.skyimage.maximum()

    
    def __call__(self,l,b):
        """ return x, y plot coordinates given l,b """
        r = self.projector.sph2pix(l,b)
        return [(r[i]-self.center[i])/self.scale[i] for i in range(2)]

    def plot_coord(self, sdir):
        """ return x,y plot coordinates given a SkyDir, or (ra,dec) tuple"""
        if not isinstance(sdir, SkyDir): sdir=SkyDir(*sdir)
        if self.galactic:
            return self(sdir.l(),sdir.b())
        return self(sdir.ra(), sdir.dec())

    def grid(self, labels=False, color='gray', textsize=8, label_offset=10):
    	"""Draws gridlines and optional labels for map.  """
    
        bs = np.arange(-90, 91, 5)
        for l in np.hstack((np.arange(0, 360, 45),[180.01])):
            lstyle = '-' if int(l)==180 or int(l)==0 else '--' 
            #self.axes.plot([self(l,b)[0] for b in bs], [self(l,b)[1] for b in bs], lstyle, color=color)
            self.axes.plot(map(lambda b: self(l,b)[0], bs), 
                           map(lambda b: self(l,b)[1], bs), lstyle, color=color)
            if labels:
                x,y = self(l, 45) 
                self.axes.text(x,y, '%3.0f'%l ,size=textsize, ha='center')

        ls = np.hstack((np.arange(180, 0, -5), np.arange(355, 180,-5), [180.01]))
        for b in np.arange(-60, 61, 30):
            lstyle = '-' if int(b)==0 else '--'
            self.axes.plot(map(lambda l: self(l,b)[0], ls), 
                           map(lambda l: self(l,b)[1], ls), lstyle, color=color)
            if labels:
                x,y = self(180.1, b)                                               
                self.axes.text(x+label_offset,y+b/60*label_offset, '%+3.0f'%b, 
                    size=textsize, ha='center',va='center')
        if labels:
            for b in [90,-90]:
                x,y = self(0,b)
                self.axes.text(x,y+b/90*label_offset,'%+3.0f'%b, 
                    size=textsize, ha='center',va='center') 

    def plot(self, sources, marker='o', text=None, 
            colorbar=False, text_kw={}, **kwargs):
        """ plot symbols at points 
        sources: list of SkyDir
        
        keywords:
            text: optional text strings (same length as sources if specified)
            text_kw: dict
                keywords for the text function
            marker: symbol to use, see scatter doc.
            
            colorbar: set True to add a colorbar. (See method to attach label)
            
        kwargs: applied to scatter, use c as an array of floats, with optional
                cmap=None, norm=None, vmin=None, vmax=None
                to make color key for another value
                in that case, you can use Axes.colorbar
                set s to change from default size =10
        """
        if 'fontsize' not in text_kw: text_kw['fontsize']=8
        X=[]
        Y=[]
        for i,s in enumerate(sources):
            x,y = self.plot_coord(s)
            X.append(x)
            Y.append(y)
            if text is not None:
                self.axes.text(x,y,text[i], **text_kw )
        self.source_cb=self.axes.scatter(X,Y, marker=marker,  **kwargs)
        if colorbar: assert False, "not implemented" #self.colorbar()

    def on_move(self, event):
        """Reports mouse's position in galactic coordinates."""
        from numpy import fabs
        if event.xdata == None or event.ydata == None:
            pass 
        else:
            try:
                coords = self.proj.pix2sph(event.xdata, event.ydata)
                self.poslabel.set_text("long=%1.2f\n  lat=%1.2f" %(coords[0],coords[1]))
            except:
                self.poslabel.set_text("")
        self.figure.canvas.draw()
                  
    def imshow(self,  title=None, scale='linear', title_kw={}, **kwargs):
        """run imshow
        scale : string defining the translation or a ufunc
            the string can specify linear [default], log for log10, sqrt, or asinh
        """
        nocolorbar =kwargs.pop('nocolorbar', self.nocolorbar)
        cb_kw = kwargs.pop('colorbar_kw',
                dict(orientation='vertical', shrink=0.6 if self.size==180 else 1.0))
        from numpy import ma
        scale_fun = kwargs.pop('fun', lambda x : x)
        # set defaults
        imshow_kw =dict(origin='lower', interpolation='nearest', extent=self.extent) 
        imshow_kw.update( kwargs)
        if self.axes is None: 
            self.figure, self.axes = pylab.subplots(1,1, figsize=(10,5))
         
        if self.size==180: self.axes.set_axis_off()
        fun_dict = dict(linear=scale_fun, log=ma.log10, sqrt=ma.sqrt, asinh=ma.arcsinh)
        fun  = fun_dict.get(scale, None) if type(scale)==types.StringType else fun
        if fun is None:
            raise Exception('bad scale function: %s, must be one of %s'%(scale,fun_dict.keys()))
        m = self.axes.imshow(fun(self.masked_image), **imshow_kw)
                                        
        if not nocolorbar:
            self.colorbar =self.axes.figure.colorbar(m, ax=self.axes, **cb_kw)
            self.colorbar.set_label(self.cbtext)
            self.mappable = self.colorbar.mappable
        else:
            self.mappable=m
        self.title(title, **title_kw)

        # for interactive formatting of the coordinates when hovering
        ##pylab.gca().format_coord = self.format_coord # replace the function on the fly!

    def pcolor(self,  title=None, scale='linear',  **kwargs):
        'run pcolor'
        from numpy import ma, array
        import pylab
        #self.axes = pylab.axes()
        if self.galactic:
            xvalues=array([self.skydir(i,0).l() for i in range(self.nx+1)])
            yvalues=array([self.skydir(0,i).b() for i in range(self.ny+1)])
            pylab.xlabel('glon'); pylab.ylabel('glat')
        else:
             xvalues=array([self.skydir(i,0).ra() for i in range(self.nx+1)])
             yvalues=array([self.skydir(0,i).dec() for i in range(self.ny+1)])
             pylab.xlabel('ra'); pylab.ylabel('dec')

        if   scale=='linear':  pylab.pcolor(self.masked_image,   **kwargs)
        elif scale=='log':     pylab.pcolor(ma.log10(self.masked_image), **kwargs)
        else: raise Exception('bad scale: %s'%scale)
                                        
        self.colorbar=pylab.colorbar(orientation='horizontal', shrink=1.0 if self.size==180 else 1.0)

        self.title(title)
        if self.axes is None: self.axes = pylab.gca()

    def axislines(self, color='black',  **kwargs):
        ' overplot axis lines'
        import pylab
        self.axes.axvline(0, color=color, **kwargs)
        self.axes.axhline(0, color=color, **kwargs)
        pylab.axis(self.extent)
 
    def title(self, text=None, **kwargs):
        ' plot a title, default the name of the SkySpectrum'
        try:
            self.axes.set_title( text if text is not None else self.skyfun.name(), **kwargs)
        except AttributeError: #no name?
            pass

    def skydir(self, x, y):
        " from pixel coordinates to sky "
        xpixel = (180-x)*float(self.nx)/360.
        ypixel = (y+90)*float(self.ny)/180.
        if self.proj.testpix2sph(xpixel,ypixel) !=0: return None #outside valid region
        sdir = SkyDir(x, y, self.proj)
        return sdir

    def pixel(self, sdir):
        """ return pixel coordinates for the skydir"""
        if self.galactic: return  self.proj.sph2pix(sdir.l(),sdir.b())
        return  self.proj.sph2pix(sdir.ra(),sdir.dec())

    def format_coord(self, x, y):
        " replacement for Axes.format_coord"
        sdir = self.skydir(x,y)
        val  = self.skyfun(sdir)

        return 'ra,dec: (%7.2f,%6.2f); l,b: (%7.2f,%6.2f), value:%6.3g' %\
            ( sdir.ra(), sdir.dec(), sdir.l(), sdir.b(), val)
                
    def scale_bar(self,  delta=1,text='$1^o$', color='k'):
        """ draw a scale bar in lower left """
        xmin, xmax= self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()
        x1,y1 = 0.95*xmin + 0.05*xmax, 0.95*ymin+0.05*ymax
        sd = self.skydir(x1,y1)
        x2,y2 = self.pixel(SkyDir(sd.ra()-delta/math.cos(math.radians(sd.dec())), sd.dec())) 
        self.axes.plot([x1,x2],[y1,y1], linestyle='-', color=color, lw=2)
        self.axes.text( (x1+x2)/2, (y1+y2)/2+self.ny/200., text, ha='center', color=color, fontsize=10)

    def box(self, image, **kwargs):
        """ draw a box at the center, the outlines of the image """
        if 'lw' not in kwargs: kwargs['lw']=2
        nx,ny = image.nx, image.ny
        corners = [(0,0), (0,ny), (nx,ny), (nx,0), (0,0) ]
        dirs = [image.skydir(x,y) for x,y in corners]
        rp = [ self.pixel(sdir) for sdir in dirs]
        self.axes.plot( [r[0] for r in rp], [r[1] for r in rp], 'k', **kwargs)

def galactic_map(skydir, axes=None, pos=(0.77,0.88), width=0.2, 
            color='w', marker='s', markercolor='r', markersize=40):
    """ 
    insert a little map showing the galactic position
        skydir: sky coordinate for point
        axes: Axes object, use gca() if None
        pos: location within the map
        width: width, fraction of map size
        color: line color
        marker, markercolor, markersize ['s', 'r', 40] plot symbol,color,size (see scatter)
    returns the AIT_grid to allow plotting other points
    """
    # create new a Axes object positioned according to axes that we are using
    if axes is None: axes=pl.gca()
    b = axes.get_position()
    xsize, ysize = b.x1-b.x0, b.y1-b.y0
    axi = axes.figure.add_axes((b.x0+pos[0]*xsize, b.y0+pos[1]*ysize, width*xsize, 0.5*width*ysize))
    ait_insert=AIT_grid(axes=axi, labels=False, color=color)
    ait_insert.plot([skydir], s=markersize, marker=marker, c=markercolor, zorder=100)
    axes.figure.sca(axes) # restore previous axes
    return ait_insert 

        
class ZEA(object):
    """ Manage a square image SkyImage
    """
    defaults = (
        ('size',     2,     'size of image in degrees'), 
        ('pixelsize',0.1,   'size, in degrees, of pixels'), 
        ('galactic', False, 'galactic or equatorial coordinates'), 
        ('fitsfile', '',    'set non-empty to write out as a FITS file'), 
        ('axes',     None,  'Axes object to use: \nif None, ...'), 
        ('nticks',   5,     'number of tick marks to attempt'), 
        ('nocolorbar', False, 'set to suppress the colorbar'),
        50*'-',
        ('proj',     'ZEA', 'projection name: can change if desired'),
        ('size2',    -1,    'vertical size, if specified' ),

    )

    @keyword_options.decorate(defaults)
    def __init__(self, center, **kwargs):
        """
        center SkyDir specifying center of image or tuple
            tuple interpreted as (l,b) or (ra,dec) depending on self.galactic

        """
        keyword_options.process(self, kwargs)
 
        if not isinstance(center, SkyDir):
            self.center=SkyDir(center[0],center[1], SkyDir.GALACTIC if self.galactic else SkyDir.EQUATORIAL)
        else:
            self.center = center
        # set up, then create a SkyImage object to perform the projection to a grid and manage an image
        self.skyimage = SkyImage(self.center, self.fitsfile, self.pixelsize, self.size, 1, self.proj, self.galactic, False, self.size2)
        
        # now extract stuff for the pylab image
        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()

        # we want access to the projection object, to allow interactive display via pix2sph function
        self.projector = self.skyimage.projector()
        self.set_axes()
        self.cid = None #callback id

    # note the 1/2 pixel offset: WCS convention is that pixel centers are integers starting from 1.
    # here we use 0 to nx, 0 to ny standard for image.

    def skydir(self, x, y):
        """ from pixel coordinates to sky """
        return SkyDir(x+0.5, y+0.5, self.projector) 

    def pixel(self, sdir):
        """ return  coordinates for the skydir
        """
        x,y = self.projector.sph2pix(sdir.ra(),sdir.dec()) \
            if not self.galactic else self.projector.sph2pix(sdir.l(),sdir.b())
        return  (x-0.5,y-0.5)

    def inside(self, sdir):
        """ is the direction sdir inside the boundary """
        x,y =self.pixel(sdir)
        return x> 0 and y>0 and x<self.nx and  y<self.ny
    def fill(self, skyfun):
        """ fill the image from a SkyFunction
            sets self.image with numpy array appropriate for imshow
        """
        if skyfun.__class__.__name__ !='PySkyFunction':
            def pyskyfun(v):
                return skyfun(SkyDir(Hep3Vector(v[0],v[1],v[2])))
            self.skyimage.fill(PySkyFunction(pyskyfun))
        else:
            self.skyimage.fill(skyfun)
        self.image = np.array(self.skyimage.image()).reshape((self.ny, self.nx))
        self.vmin ,self.vmax = self.skyimage.minimum(), self.skyimage.maximum()
        return self.image

    def set_axes(self):
        """ configure the axes object
           if self.axes is None, simply use gca()
          (set coordinate scale offset by 0.5 from WCS standard)
        """
        if self.axes is None:
            figure = pyplot.gcf()
            if len(figure.get_axes())==0:
                # no axes in the current figure: add one that has equal aspect ratio
                h,w = figure.get_figheight(), figure.get_figwidth()
                if w>h:
                    figure.add_axes((0.18, 0.15, h/w*0.75, 0.75))
                else:
                    figure.add_axes((0.18, 0.15, 0.75, w/h*0.75))
            
            self.axes=pyplot.gca()
        self.axes.set_aspect(1)
        self.axes.set_xlim((0.0,self.nx))
        self.axes.set_ylim((0.0,self.ny))
        self.axes.set_autoscale_on(False) 
        r =Rescale(self,self.nticks, galactic = self.galactic)

        r.apply(self.axes)

        labels = ['l','b'] if self.galactic else ['RA','Dec'] 
        self.axes.set_xlabel(labels[0]);self.axes.set_ylabel(labels[1])
        
    def grid(self, nticks=None, **kwargs):
        """ draw a grid
        """
        label_offset = 5 # units degrees?

        if nticks is None: nticks=self.nticks

        r = Rescale(self, nticks, galactic = self.galactic)
        r.apply(self.axes)
        self.axes.xaxis.set_ticks_position('none')
        try:
            # need this for 1.0.0 due to bug that turns off labels
            self.axes.yaxis.set_tick_params(which='both', right=False, left=False)
        except:
            self.axes.yaxis.set_ticks_position('none')
        
        uticks, vticks = r.uticks, r.vticks
        cs = SkyDir.GALACTIC if self.galactic else SkyDir.EQUATORIAL
        for u in uticks:
            w = [self.pixel(SkyDir(u,v,cs)) for v in  np.linspace(r.vmin,r.vmax, 2*nticks)]
            self.axes.plot([q[0] for q in w], [q[1] for q in w], '-k', **kwargs)
        for v in vticks:
            w = [self.pixel(SkyDir(u,v,cs)) for u in np.linspace(r.ul, r.ur,2*nticks)]
            self.axes.plot([q[0] for q in w], [q[1] for q in w], '-k', **kwargs)
        return r


    def scale_bar(self,  delta=1,text='$1^o$', color='k'):
        """ draw a scale bar in lower left """
        xmin, xmax= self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()
        x1,y1 = 0.95*xmin + 0.05*xmax, 0.95*ymin+0.05*ymax
        sd = self.skydir(x1,y1)
        if self.galactic:
            x2,y2 = self.pixel(SkyDir(sd.l()-delta/math.cos(math.radians(sd.b())), sd.b(),SkyDir.GALACTIC)) 
        else:
            x2,y2 = self.pixel(SkyDir(sd.ra()-delta/math.cos(math.radians(sd.dec())), sd.dec())) 
        self.axes.plot([x1,x2],[y1,y1], linestyle='-', color=color, lw=3)
        self.axes.text( (x1+x2)/2, (y1+y2)/2+self.ny/80., text, ha='center', color=color)

    def imshow(self, **kwargs):
        """ run imshow on the image, presumably set by a fill: set up for colorbar.
        
        """
        if 'image' not in self.__dict__: raise Exception, 'no data to show: must run fill first'
        # set up kw 
        self.imshow_kw = dict(cmap=None, norm=None, interpolation='nearest')
        self.imshow_kw.update(kwargs)
        fun = self.imshow_kw.pop('fun', lambda x: x)
        self.cax = self.axes.imshow(fun(self.image), **self.imshow_kw)
        
    def colorbar(self, label=None, **kwargs):
        """ 
        draw a color bar using the pylab colorbar facility
        note that the 'shrink' parameter needs to be adjusted if not a full figure
        Must have called imshow, which will will be used for default cmap, norm
        
        returns the colorbar object
        """
        if self.nocolorbar: return
        if 'cax' not in self.__dict__: raise Exception('You must call imshow first')
        cb_kw=dict(orientation= 'vertical',
                pad    = 0.01,
                ticks  = ticker.MaxNLocator(4),   
                fraction=0.10,
                shrink  = 1.0,
                cmap    = self.imshow_kw['cmap'],
                norm    = self.imshow_kw['norm'],
                )
        cb_kw.update(kwargs)
                
        self.cb=self.axes.figure.colorbar(self.cax,  **cb_kw)
        if label is not None: self.cb.set_label(label)
        return self.cb

    def box(self, image, **kwargs):
        """ draw a box at the center, the outlines of the image
            +image An object of this class, or implementing the skydir function
        """
        if 'lw' not in kwargs: kwargs['lw']=2
        nx,ny = image.nx, image.ny
        corners = [(0,0), (0,ny), (nx,ny), (nx,0), (0,0) ]
        dirs = [image.skydir(x,y) for x,y in corners]
        rp = [ self.pixel(sdir) for sdir in dirs]
        self.axes.plot( [r[0] for r in rp], [r[1] for r in rp], 'k', **kwargs)


    def plot_source(self, name, source, symbol='+', fontsize=10, **kwargs):
        """ plot symbols at points
            name: text string 
            source: a SkyDir
        """
        x,y = self.pixel(source)
        if x<0 or x> self.nx or y<0 or y>self.ny: return False
        self.axes.plot([x],[y], symbol, 
            markersize=kwargs.pop('markersize',12), **kwargs)
        #self.axes.text(x,y, name, fontsize=fontsize, **kwargs)
        if name is not None:
            self.axes.text( x+self.nx/100., y+self.nx/100., name, 
                fontsize=fontsize, **kwargs)
        return True

    def plot(self, sources, marker='o', text=None, 
            colorbar=False, text_offset=(0,0),text_kw={},
            c=None, s=20, **kwargs):
        """ plot symbols at points 
        see AIT.plot
        """
        if 'fontsize' not in text_kw: text_kw['fontsize']=8
        X=[]
        Y=[]
        C=[] if c is not None else None
        dx,dy=text_offset
        for i,t in enumerate(sources):
            if not self.inside(t): continue
            x,y = self.pixel(t)
            X.append(x)
            Y.append(y)
            if c is not None: C.append(c[i])
            if text is not None:
                self.axes.text(x+dx,y+dy,text[i], **text_kw )
        self.source_cb=self.axes.scatter(X,Y,c=C, marker=marker, s=s, **kwargs)
        #if colorbar: assert False, "not implemented" #self.colorbar()
        if colorbar: 
            plt.colorbar(self.source_cb) #self.colorbar()

    def cross(self, sdir, size, text=None, **kwargs):
        """ draw a cross at sdir,
            size: half-length of each arm, in deg.
        
        """    
        x,y = self.pixel(sdir)
        if x<0 or x> self.nx or y<0 or y>self.ny: return False
        pixelsize = self.pixelsize
        delta = size/pixelsize
        axes = self.axes
        axes.plot([x-delta, x+delta], [y,y], '-k', **kwargs)
        axes.plot([x,x], [y-delta, y+delta], '-k', **kwargs)
        if text is not None:
            if 'lw' in kwargs: kwargs.pop('lw') # allow lw for the lines. 
            axes.text(x,y, text, **kwargs)
        return True
        
    def ellipse(self, sdir, par, symbol='-', **kwargs ):
        """ sdir: SkyDir or (ra, dec) or (l,b)
            ellipse parameters: a, b, ang (all deg) 
        """
        if not isinstance(sdir, SkyDir):
            sdir = SkyDir(sdir[0],sdir[1], SkyDir.GALACTIC if self.galactic else SkyDir.EQUATORIAL)

        x0,y0 = self.pixel(sdir)
        a,b,ang = par
        if self.galactic:
            # adjust angle for local coordinate
            up = SkyDir(sdir.ra()+a, sdir.dec())
            xup,yup =self.pixel(up)
            ang += np.degrees(np.arctan2( yup-y0,xup-x0))
        ellipse = Ellipse([a,b,np.radians(ang)])
        x,y = ellipse.contour(1.0)
        pixelsize=self.pixelsize #scale for plot
        self.axes.plot(x0+np.asarray(x)/pixelsize, y0+np.asarray(y)/pixelsize, symbol, **kwargs)
        
    def clicker(self, onclick=None):
        """ enable click processing: default callback prints the location
            callback example:
            def default_onclick(event):
                print 'button %d, %s' % (event.button, zea.skydir(event.xdata,event.ydata))
        """
        def default_onclick(event):
            print 'button %d, %s' % (event.button, self.skydir(event.xdata,event.ydata))
        if onclick==None: onclick=default_onclick
        if self.cid is not None: self.noclicker()
        self.cid=self.axes.figure.canvas.mpl_connect('button_press_event', onclick)
    
    def noclicker(self):
        if self.cid is not None: self.axes.figure.canvas.mpl_disconnect(self.cid)
        self.cid=None
    
    def smooth(self,scale=0.1,smoother=None):
        """ smooth the image using a Gaussian kernel.  Reverse process by calling ZEA.unsmooth.

            NB -- if more than one image with the same dimension is to be smoothed, it is
            more computationally efficient to create a single GaussSmoothZEA object and use it
            to make smoothed images.  This can be done manually, or using the smoother returned
            by this message and passing it as the argument for future calls to smooth for
            other ZEA objects.

            scale: the smoothing scale (std. dev.) in deg

            returns: the GaussSmoothZEA object for use in smoothing additional images
        """
        if 'image' not in self.__dict__.keys(): return
        gsz = smoother if smoother is not None else GaussSmoothZEA(self,scale)
        self.original = self.image.copy()
        self.image    = gsz(self.image)
        return gsz

    def unsmooth(self):
        """ replace smoothed image with original unsmoothed image."""
        if 'original' not in self.__dict__.keys(): return
        self.image = self.original
        self.__dict__.pop('original')

    def circle(self, sdir, size, fill=False, **kwargs):
        """draw a cirle
        """
        self.axes.add_artist(plt.Circle(self.pixel(sdir), size/self.pixelsize, fill=fill, **kwargs))

def ZEA_test(ra=90, dec=80, size=5, nticks=8, galactic=False, **kwargs):
    """ exercise (most) everything """
    pyplot.clf()
    q = ZEA(SkyDir(ra,dec), size=size, nticks=nticks, galactic=galactic, **kwargs)
    q.grid(color='gray')
    q.scale_bar(1, '$1^0$')
    q.axes.set_title('test of ZEA region plot')
    t=q.cross( SkyDir(ra,dec), 1, 'a red cross, arms +/- 1 deg', color='r', lw=2)
    if not t: print 'failed to plot the cross'
    q.plot_source('(80,76)', SkyDir(80,76), 'd')
    q.plot_source('(110,74)', SkyDir(110,74), 'x')
    q.ellipse( SkyDir(ra,dec), (1, 0.25, 0))
    for dec in np.arange(-90, 91, 2):
        q.plot_source( '(%d,%d)'%(ra,dec), SkyDir(ra,dec), 'x')
    #def myfun(v): return v[0] # x-component of skydir to make a pattern
    #q.fill(PySkyFunction(myfun))
    q.fill(lambda sd: sd.ra())
    q.imshow()
    q.colorbar()
    q.axes.figure.show()
    return q

class ZEA_from_fits(ZEA):
    """ subclass of ZEA with input FITS file, produced by ZEA"""
    defaults = ZEA.defaults

    @keyword_options.decorate(defaults)
    def __init__(self, filename, **kwargs):
        """
        filename: string | skymaps.SkyImage object
        """
        keyword_options.process(self,kwargs)
        if not isinstance(filename, SkyImage):
            try:
                self.skyimage =SkyImage(filename)
            except Exception, msg:
                raise Exception('failed to load file %s: %s)' % (filename, msg))
        else:
            self.skyimage=filename
        s = self.skyimage
        image_array = np.array(s.image())
        self.nx,self.ny = s.naxis1(), s.naxis2()
        self.projector = s.projector()
        image_array.resize(self.nx,self.ny)
        self.image=image_array  
        self.center = SkyDir(*self.projector.pix2sph(self.nx/2,self.ny/2.))
        left = SkyDir(*self.projector.pix2sph(0,self.ny/2.))
        edge = SkyDir(*self.projector.pix2sph(self.nx,self.ny/2.))
        self.size = 2.*np.degrees(self.center.difference(edge))


        self.set_axes()

class TSplot(object):
    """
    Create a "TS plot" 
    Uses the ZEA class for display

    """
    defaults = dict(
        pixelsize  = None,  # passed to ZEA
        axes       = None, # passed to ZEA
        nticks     = 4,    # passed to ZEA
        fitsfile   = '', 
        galmap     = True, # determine if overlay a galactic map
        galpos     = (0.77,0.88), # position if galmap set
        scalebar   = True,
            )
    def __init__(self, tsmap, center, size,  **kwargs):
        """
        parameters:
        *tsmap*   a SkyFunction, that takes a SkyDir argument and returns a value
        *center*  SkyDir direction to center the plot
        *size*    (degrees)  width of plot
        *pixelsize* [None] size (degrees) of a pixel: if not specified, will be size/10
        *axes* [None] Axes object to use: if None, use current
        *nticks* [4] Suggestion for labeling
        *fitsfile*[''] 
        *galmap* [True] overplot a little map in galactic coordinates showing the position
        *scalebar" [True] overplot a scalebar in lower left
        **kwargs  additional args for ZEA, like galactic
        """
        self.__dict__.update(TSplot.defaults) 
        for key in self.__dict__.keys():
            if key in kwargs: self.__dict__[key] = kwargs.pop(key)
        self.tsmap = tsmap
        self.size=size
        if self.pixelsize is None: self.pixelsize=size/10.
        npix = round(size/self.pixelsize)
        self.pixelsize=size/npix
        self.zea= ZEA(center, size=size, pixelsize=self.pixelsize, axes=self.axes, 
                nticks=self.nticks,fitsfile=self.fitsfile, **kwargs)
        print 'TSplot: filling %d pixels (size=%.2f, npix=%d)...'%( (size/self.pixelsize)**2, size, npix)
        sys.stdout.flush()
        self.zea.fill(tsmap)
        # create new image that is the significance in sigma with respect to local max
        self.tsmaxpos=tsmaxpos = find_local_maximum(tsmap, center) # get local maximum, then check that is in the image
        x,y = self.zea.pixel(tsmaxpos)
        if x>=0 and x < self.zea.nx and y>=0 and y<self.zea.ny:
            tsmaxval = tsmap(tsmaxpos)
        else: # not in image: use maximm actually in the image
            tsmaxval = self.zea.image.max()
        tmap = tsmaxval-self.zea.image
        tmap[tmap<0] = 0
        self.image =np.sqrt(tmap)
        self.cb=None
        # np.sqrt(-2* np.log(1-np.array([0.68,0.95, 0.99]))
        self.clevels = np.array([1.51, 2.45, 3.03])

         
    def show(self, colorbar=True):
        """
        Generate the basic plot, with contours, scale bar, color bar, and grid
        """
        norm2 = mpl.colors.Normalize(vmin=0, vmax=5)
        cmap2 = mpl.cm.hot_r

        self.nx,self.ny = self.image.shape
        axes = self.zea.axes
        t = axes.imshow(self.image, origin='lower', 
            extent=(0,self.nx,0,self.ny),
            cmap=cmap2, norm=norm2,
            interpolation='bilinear')
        if colorbar:
            if self.cb is None:
                self.cb = pyplot.colorbar(t, ax=axes, 
                    cmap=cmap2, norm=norm2,
                    ticks=ticker.MultipleLocator(),
                    orientation='vertical',
                    #shrink=1.0,
                    )
            self.cb.set_label('$\mathrm{sqrt(TS difference)}$')

        ct=axes.contour(  np.arange(0.5, self.nx,1), np.arange(0.5, self.ny, 1), self.image,
            self.clevels, 
            colors='k', linestyles='-' )
        if axes.get_xlim()[0] !=0:
            print 'Warning: coutour modified: limits', axes.get_xlim(), axes.get_ylim()
        cfmt={} 
        for key,t in zip(self.clevels,['68%','95%', '99%']): cfmt[key]=t
        plt.clabel(ct, fmt=cfmt, fontsize=8)
        #axes.set_xlim((0,nx)); axes.set_ylim((0,ny))
        #print 'after reset', axes.get_xlim(), axes.get_ylim()
        if self.scalebar:
            t = 3
            if self.size< 0.03*t:
                self.zea.scale_bar(1/60.,  "1'", color='w')
            elif self.size<0.6*t:
                self.zea.scale_bar(0.1, "$0.1^o$", color='w')
            elif self.size<1.1*t:
                self.zea.scale_bar(0.5, "$0.5^o$", color='w')
            else:
                self.zea.scale_bar(1.0, '$1^o$', color='w')
        self.zea.grid(color='gray')
        
        if self.galmap:
            galactic_map(self.zea.center, axes=self.axes, 
                color='w', marker='s', markercolor='r', pos=self.galpos);
 

    def overplot(self, quadfit, sigma=1.0,contours=None, **kwargs):
        """
        OVerplot contours from a fit to surface
        
        quadfit: either: a uw.like.quadform.Localize object used for the fit, 
                        or: an array [ra,dec, a,b, phi, qual] 
        sigma:  scale factor
        contours: [None] default is the 68,95,99% assuming the radius is 1 sigma
                if specified, must be a list, eg [1.0]
        """
        axes = self.zea.axes
        if contours is None: contours = self.clevels
        if getattr(quadfit,'__iter__', False):
            ra,dec,a,b,phi=quadfit[:5]
            qual = None if len(quadfit)==5 else quadfit[5]
            ellipse = Ellipse([a,b,np.radians(phi)])
            x,y = ellipse.contour(1.0)
        else:
            x,y = quadfit.ellipse.contour(quadfit.fit_radius)
            ra,dec = quadfit.ra, quadfit.dec
            qual = quadfit.quality()
            
        pixelsize=self.zea.pixelsize #scale for plot
        x0,y0 = self.zea.pixel(SkyDir(ra,dec))
        f=sigma/pixelsize #scale factor
        xa = f*np.array(x)
        ya = f*np.array(y)
        axes.plot([x0],[y0], '+g')
        for r in contours:
            axes.plot(r*xa+x0,r*ya+y0, '--g', **kwargs);
        if qual is not None:
            axes.text(0.95, 0.07,'fit quality=%.1f'%qual, color='g',
                horizontalalignment='right', fontsize=10,
                transform = axes.transAxes)

    def cross(self, sdir, size, label=None,  fontsize=12, markersize=10,  **kwargs):
        """ make a cross at the position, size defined in celestial coordinats
        """
      
        image=self.zea
        x,y = image.pixel(sdir)
        pixelsize = image.pixelsize
        delta = size/pixelsize
        axes = self.zea.axes
        axes.plot([x-delta, x+delta], [y,y], '-k', **kwargs)
        axes.plot([x,x], [y-delta, y+delta], '-k', **kwargs)
        if label is not None:
            nx = image.nx 
            if 'lw' in kwargs: kwargs.pop('lw')
            image.axes.text( x+nx/100., y+nx/100., label, fontsize=fontsize, **kwargs)

    def plot(self, loc, label=None, symbol='+',  fontsize=12, markersize=10,  **kwargs):
        """ plot a single point at the celestial location
            return the tsmap value there
        """
        image = self.zea
        nx = image.nx 
        x,y = image.pixel(loc)
        image.axes.plot([x], [y], symbol, markersize=markersize, **kwargs)
        if label is not None:
            textargs = {}
            if 'color' in kwargs: textargs['color']=kwargs['color']
            image.axes.text( x+nx/100., y+nx/100., label, fontsize=fontsize, **textargs)
        return self.tsmap(loc) if hasattr(self, 'tsmap') else None
        
    def clicker(self, onclick=None):
        return self.zea.clicker(onclick)
    def noclicker(self):
        self.zea.noclicker()
        
    def saveto(self, filename):
        """ (in progress)
        more reliable way to save the image as is to a FITS file
        """
        self.skyimage.reimage(self.gcdir, filename, self.pixelsize, self.size, self.proj, True)
    
class TSplotFromFITS(TSplot):
    """subclass of TSplot, load data from FITS image file"""

    def __init__(self, filename, **kwargs):
        self.__dict__.update(TSplot.defaults)
        try:
            self.zea = ZEA_from_fits(filename)
        except Exception, msg:
            raise Exception('failed to load file %s: %s)' % (filename,msg) )
        if not isinstance(filename, SkyImage):
            t =os.path.split(os.path.splitext(filename)[0])[-1].split('_')
            self.sourcename= ' '.join(t[:-1]).replace('p', '+')
        else:
            self.sourcename = kwargs.pop('sourcename', '(not specified)')
        zi = self.zea.image
        self.image=np.sqrt(zi.max()-zi)  
        self.cb = None
        self.clevels = np.array([1.51, 2.45, 3.03])
        self.scalebar=True
        self.size = self.zea.size


class GaussKernel(object):

    def __init__(self,center,scale):
        """Center -- the center of an image to be Gaussian-smoothed
           scale  -- the standard deviation of the Gaussian in degrees
        """
        self.center = center
        self.scale  = np.radians(scale)
        self.norm   = 1./(self.scale*(2*np.pi)**0.5)

    def __call__(self,v,skydir = None):
        sd = skydir or SkyDir(Hep3Vector(v[0],v[1],v[2]))
        diff = self.center.difference(sd)/self.scale
        if diff > 5: return 0
        return self.norm*exp(-0.5*diff**2)        

    def pyskyfun(self): return PySkyFunction(self)


class GaussSmoothZEA(object):
    """ Support Gaussian smoothing for ZEA images.  Create a ZEA object with the desired
        properties and then use it to instantiate this class.  The instance of this class
        will store the Fourier transform of the gaussian kernel.  Then, any ZEA image
        with the same dimension can be quickly smoothed using the same GaussSmoothZEA
        instance.
    """

    def __init__(self,zea,scale):
        """zea is an instance of ZEA.
           scale is the smoothing scale in degrees."""
        gc = GaussKernel(zea.center,scale)
        im_backup = zea.image.copy()
        zea.fill(gc.pyskyfun())
        zea.image /= zea.image.sum()
        self.fft_kernel = fft2(zea.image)
        zea.image = im_backup

    def __call__(self,image):
        """image is a 2d numpy array which must have same dimensions as ZEA used to create this object."""
        fft_image = fft2(image)
        return np.real(fftshift(ifft2(self.fft_kernel*fft_image)))
        
def find_local_maximum( mapfun, startpos):
    """ 
        looks for local maximum for a SkyFunction, starting at startpos
    """
    class LocalMax(object):
        """ helper class """
        def __init__(self, mapfun, center):
            self.tsf=mapfun
            self.sdir = center
            self.ra,self.dec = self.sdir.ra(), self.sdir.dec()
            self.cdec= math.cos(math.degrees(self.dec))
        def __call__(self,par):
            ra = self.ra+par[0]/self.cdec
            dec= self.dec+par[1]
            return -self.tsf(SkyDir(ra,dec))
        def find(self):
            dx,dy = optimize.fmin(self, (0,0),disp=0)
            return SkyDir(self.ra+dx/self.cdec, self.dec+dy)

    return LocalMax(mapfun, startpos).find()
    
if __name__=='__main__':
    ZEA_test(dec=0, nticks=6)
    pass

