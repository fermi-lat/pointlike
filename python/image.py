""" image processing

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/image.py,v 1.3 2007/11/24 19:50:27 burnett Exp $

"""

import matplotlib

def date_tag():
    import datetime, pylab
    pylab.figtext(0.04, 0.02, str(datetime.datetime.today())[:16], size=8)

#---------------------------------------------------------------------------
def energy_range(level):
    " energy range associated with Healpix level: note wired in"
    import math
    emin = 100.*math.pow( 2.35, level-6)
    emax =2.35*emin
    return emin,emax

#---------------------------------------------------------------------------
class Image(object):
    def __init__(self, fun, s, level=9, scale = 0.5, step = 0.02):
        """fun is a SkySpectrum object, s either has a dir() or ra(),dec() method"""

        import math
        grid=arange(-scale, scale+step/2, step)
        emin, emax = energy_range(level)
        if 'ra' in dir(s): ra, dec = s.ra(), s.dec()
        else: ra,dec = s.dir().ra(), s.dir().dec()
        cosdec = math.cos(math.radians(dec))
        self.image = array([fun.integral(SkyDir(ra-dx/cosdec, dec-dy), emin, emax)\
                   for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
        self.scale, self.level, self.ra, self.dec = scale, level, ra, dec
    def show(self, **kwargs):
        import pylab
        scale=self.scale
        pylab.imshow(self.image, extent=[-scale, scale, -scale, scale],interpolation='nearest', **kwargs)
        pylab.axvline(0, color='white')
        pylab.axhline(0, color='white')
        pylab.colorbar()
        emin, emax = energy_range(self.level)
        pylab.title('level %d (%d-%d)' %(self.level, emin, emax), size=10)

#---------------------------------------------------------------------------
def make_image( fun, sdir, level=9, scale = 0.5, step = 0.02):
    """ Return an image array ready to be plotted by imshow
        fun: a SkySpectrum object
        sdir: a SkyDir to be the center
        level: integrate the SkySpectrum over the energy range
        scale: half-size, in degrees
        step:  pixel size, degrees
    """
    import math
    grid=arange(-scale, scale+step/2, step)
    emin, emax = energy_range(level)
    ra, dec = sdir.ra(), sdir.dec()
    cosdec = math.cos(math.radians(dec))
    # note reversal of both coordinates
    image = array([fun.integral(SkyDir(ra-dx/cosdec, dec-dy), emin, emax)\
                   for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
    return image

#---------------------------------------------------------------------------
def show_image( fun,sdir, level=9, scale=0.5, step=0.01):
    import pylab
    pylab.imshow(make_image(fun, sdir, level, scale, step), extent=[-scale, scale, -scale, scale])
    pylab.axvline(0, color='white')
    pylab.axhline(0, color='white')
    pylab.colorbar()
    emin, emax = energy_range(level)
    pylab.title('level %d (%d-%d)' %(level, emin, emax), size=10)

#---------------------------------------------------------------------------
class LogNorm(matplotlib.colors.normalize):
    ' replacement for the non-existent colors.LogNorm'
    def __init__(self, vmin=None, vmax=None):
        self.vmin, self.vmax= vmin, vmax
    def __call__(self, val):
        ' expect val to be a masked arrray'
        from matplotlib.numerix import ma, minimum, maximum
        logval = ma.log10(val)
        #self.autoscale(logval) #calls base class
        vmin, vmax = ma.minimum(logval), ma.maximum(logval) #self.vmax
        print 'autoscale: vmin =%f, vmax=%f' % (vmin, vmax)
        return (logval-vmin)/(vmax-vmin)

#---------------------------------------------------------------------------
class AIT(object):
    """ Manage a full-sky image of a SkyProjection or SkyFunction, wrapping SkyImage
     """
    
    def __init__(self, skyfun, pixelsize=0.5, galactic=True, fitsfile='', proj='AIT'):
        """
        skyfun SkyProjection or SkyFunction object
        pixelsize [0.5] size, in degrees, of pixels
        galactic [True] galactic or equatorial coordinates
        fitsfile [''] if set, write the projection to a FITS file
        proj ['AIT'] could be 'CAR' for carree: used by wcslib

        """
        from pointlike import SkyImage, SkyDir
        from numpy import array, isnan, ma
        from pylab import normalize
        
        self.skyfun = skyfun
        self.galactic = galactic
        self.pixelsize = pixelsize
        # set up, then create a SkyImage object to perform the projection to a grid
        center = SkyDir(0,0, SkyDir.GALACTIC if galactic else SkyDir.EQUATORIAL)
        self.skyimage = SkyImage(center, fitsfile, pixelsize, 180, 1, proj, galactic)
        self.skyimage.fill(skyfun)
        
        # now extract stuff for the pylab image, creating a masked array to deal with the NaN values
        self.nx, self.ny = self.skyimage.naxis1(), self.skyimage.naxis2()
        self.image = array(self.skyimage.image()).reshape((self.ny, self.nx))
        self.mask = isnan(self.image)
        self.masked_image = ma.array( self.image, mask=self.mask)
        self.extent = (180,-180, -90, 90)
        self.vmin ,self.vmax = self.skyimage.minimum(), self.skyimage.maximum()

        # we want access to the projection object, to allow interactive display via pix2sph function
        self.proj = self.skyimage.projector()
        self.x = self.y = 100 # initial def
                  
    def show(self,  title=None, scale='linear',  **kwargs):
        'run imshow'
        from numpy import ma
        import pylab
        # change defaults
        if 'origin'        not in kwargs: kwargs['origin']='lower'
        if 'interpolation' not in kwargs: kwargs['interpolation']='nearest'
        if 'extent'        not in kwargs: kwargs['extent']=self.extent
        
        if   scale=='linear':  pylab.imshow(self.masked_image,   **kwargs)
        elif scale=='log':     pylab.imshow(ma.log10(self.masked_image), **kwargs)
        else: raise Exception('bad scale: %s'%scale)
                                        
        pylab.colorbar()
        if self.galactic:
            pylab.xlabel('glon'); pylab.ylabel('glat')
        else:
            pylab.xlabel('ra'); pylab.ylabel('dec')
        self.title(title)

        # for interactive formatting of the coordinates when hovering
        pylab.gca().format_coord = self.format_coord # replace the function on the fly!


    def axes(self, color='black',  **kwargs):
        ' overplot axis lines'
        import pylab
        pylab.axvline(0, color=color, **kwargs)
        pylab.axhline(0, color=color, **kwargs)
        pylab.axis(self.extent)
 
    def title(self, text=None, **kwargs):
        ' plot a title, default the name of the SkySpectrum'
        import pylab
        try:
            pylab.title( text if text else self.skyfun.name(), **kwargs)
        except AttributeError: #no name?
            pass

    def format_coord(self, x, y):
        " replacement for Axes.format_coord"
        from pointlike import SkyDir
        xpixel = (180-x)*float(self.nx)/360.
        ypixel = (y+90)*float(self.ny)/180.
        if self.proj.testpix2sph(xpixel,ypixel) !=0: return '' #outside valid region
        l, b = self.proj.pix2sph(xpixel,ypixel)
        sdir = SkyDir(l, b, SkyDir.GALACTIC if self.galactic else SkyDir.EQUATORIAL)
        val  = self.skyfun(sdir)

        return 'ra,dec: (%7.2f,%6.2f); l,b: (%7.2f,%6.2f), value:%6.3g' %\
            ( sdir.ra(), sdir.dec(), sdir.l(), sdir.b(), val)
                
           
