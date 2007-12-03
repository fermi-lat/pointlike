"""
Test addition of diffuse, nearby sources

"""
from numpy import arange, array, log10
from pointlike import *
galdiffuse_file = r'F:\glast\extlib\extFiles\v0r7\galdiffuse\GP_gamma_v0r0p1.fits'
galdiffuse_file = r'F:\glast\extlib\extFiles\v0r7\galdiffuse\GP_gamma.fits'
galdiffuse = DiffuseFunction(galdiffuse_file)
PointSourceLikelihood.set_levels(8)  # always use this minimum level

def pixeldata(pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'):
    return PhotonMap(pixelfile, 'PHOTONMAP')
                          
def date_tag():
    import datetime, pylab
    pylab.figtext(0.04, 0.02, str(datetime.datetime.today())[:16], size=8)

def energy_range(level):
    " energy range associated with Healpix level: note wired in"
    import math
    emin = 100.*math.pow( 2.35, level-6)
    emax =2.35*emin
    return emin,emax


def plot_diffuse():
    import pylab as pl
    energies = 10**(arange(1.5,5,0.05))
    pole = SkyDir(0,90, SkyDir.GALACTIC)
    center = SkyDir(0,0, SkyDir.GALACTIC)
    pole_flux = array([galdiffuse(pole, e) for e in energies])
    center_flux = array([galdiffuse(center, e) for e in energies])
    extra_gal = array([galdiffuse.extraGal(e) for e in energies])
    pl.cla();
    pl.loglog( energies, (center_flux), '-', lw=2, label = 'galactic center')
    pl.loglog( energies, (pole_flux), '-r', lw=2, label = 'North pole')
    pl.loglog( energies, (extra_gal),  '-g', lw=2, label = 'extra galactic')
    pl.legend()
    pl.axis(( 30, 4e4, 1e-10, 1e-4))
    pl.xlabel('energy (MeV)')
    pl.ylabel(r'$\rm{flux\ (ph\ cm^{-2}\ s^{-1}\ MeV^{-1})}$')
    pl.grid()
    pl.figtext(0.1,0.95, 'Diffuse file: %s'%galdiffuse_file, size=10)
    date_tag()
    

def diffuse_data(layer=0):
    import pyfits
    return pyfits.open(galdiffuse_file)[layer].data

def sc2map():
    pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
    return PhotonMap(pixelfile, 'PHOTONMAP')

    
def vela_fit(photonmap):
    vela_dir = SkyDir(128.73, -45.2 )
    like = PointSourceLikelihood(photonmap, 'vela', vela_dir)
    like.printSpectrum()
    return like

def strong_like(photonmap):
    like = PointSourceLikelihood(photonmap, 'DC2_3EGJ1635m1751', SkyDir(248.788, -17.861))
    like.set_verbose(1)
    like.printSpectrum()
    return like

def weak_like(photonmap):
    like = PointSourceLikelihood(photonmap, 'HLCloud_SC1_05', SkyDir(248.4804, -18.294))
    like.set_verbose(1)
    like.printSpectrum()
    return like

def make_fits(photonmap):
    return (strong_like(photonmap), weak_like(photonmap))

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
    def show(self):
        import pylab
        scale=self.scale
        pylab.imshow(self.image, extent=[-scale, scale, -scale, scale],interpolation='nearest')
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

def show_image( fun,sdir, level=9, scale=0.5, step=0.01):
    import pylab
    pylab.imshow(make_image(fun, sdir, level, scale, step), extent=[-scale, scale, -scale, scale])
    pylab.axvline(0, color='white')
    pylab.axhline(0, color='white')
    pylab.colorbar()
    emin, emax = energy_range(level)
    pylab.title('level %d (%d-%d)' %(level, emin, emax), size=10)
    
                

def multi_image( fun, sdir,  levels=range(8,12), scale=1, step=0.02):
    import pylab as pl
    nplots = len(levels)
    pl.figtext(0.05, 0.95, '%s at (%5.2f,%5.2f)'% (fun.name(), sdir.ra(), sdir.dec()), size=10)
    if nplots>1:
        for i in range(nplots):
            pl.subplot(2,(nplots+1)//2,i+1)
            show_image(fun, sdir, levels[i], scale, step)
    else:
        show_image(fun, sdir, levels[0], scale, step)
        
    date_tag()


def pixelArea(level):
    """ area, in sr, for a pixel at the given level"""
    return Healpix(2**level).pixelArea()


def rescale_fit( like, level=8):
    """ return a SkySpectrum from a fit scaled to photons/pixel with the level
    """
    ret = CompositeSkySpectrum()
    area = Healpix(2**level).pixelArea()
    ret.add( like, area)


def multiback( gal, strong, weak,  factors=[1,2,3]):
    for f in factors:
        back = CompositeSkySpectrum(gal, f*1e10)
        back.add(strong)
        PointSourceLikelihood.set_diffuse(back)
        weak.setDir(weak.dir())
        weak.maximize()
        w = [weak[level].TS() for level in range(8,12)]
        print '%5.1f %6.1f %6.1f %6.1f %6.1f %6.1f' % tuple([f]+w+ [sum(w)])

class Candidate(object):
    def __init__(self,photonmap):
        """ refit candidates from 5-degree circle around the 3EG blazar"""    
        candinfo=\
        """UW_J2531m179  253.056 -17.8938 
            UW_J2533m173  253.316 -17.3096 
            UW_J2529m193  252.9   -19.3173
            UW_J2519m192  251.938 -19.2024 
            UW_J2524m161  252.418 -16.1355 
            UW_J2513m166  251.271 -16.5517 
            UW_J2517m145  251.683 -14.5489 
            UW_J2503m150  250.349 -14.9796 
            UW_J2488m179  248.782 -17.8592 
            UW_J2485m183  248.499 -18.3096 
            UW_J2483m187  248.333 -18.7189 
            UW_J2485m179  248.514 -17.8798 
            UW_J2483m181  248.273 -18.1156 
            UW_J2473m206  247.314 -20.6061 
            UW_J2462m187  246.214 -18.7433 
            UW_J2460m185  245.95   -18.503
            UW_J2471m176  247.118 -17.5794 
            UW_J2481m153  248.149 -15.3457 
            UW_J2463m167  246.338 -16.7067 
            UW_J2452m177  245.229 -17.72""".split('\n')
        self.name=[line.split()[0] for line in candinfo]
        self.sdir=[SkyDir(float(line.split()[1]),float(line.split()[2])) for line in candinfo]
        self.fitdir=[]
        self.like=[]
        self.TS=[]
        for i in range(len(self.name)):
            print 'fitting %s at (%5.2f,%5.2f)' % (self.name[i] , self.sdir[i].ra(), self.sdir[i].dec()),
            l = PointSourceLikelihood(photonmap, self.name[i], self.sdir[i])
            self.like.append(l)
            l.set_verbose(0)
            Ts = l.maximize()
            print 'TS=%6.0f' % Ts
            self.TS.append(Ts)
            l.localize(1,5)
            self.fitdir.append(l.dir())
    def __getitem__(self, n):  return self.like[n]
    def __len__(self):        return len(self.like)

    def plot(self):
        import pylab
        ra = [ d.ra() for d in self.fitdir]
        dec= [d.dec() for d in self.fitdir]
        pylab.plot( ra, dec, '+')
        pylab.axvline(248.788, color='red'); pylab.axhline(-17.861, color='red')
        for i in range(len(self.TS)):
            pylab.text( ra[i], dec[i], '%4.0f' % self.TS[i], size=10)
        pylab.plot([248.4804], [ -18.294], 'ob')
        
                        
    def refit(self, center, gal=galdiffuse, exposure=3e10):
        import math
        cdir = center.dir()
        print 'refitting candidates within 1 deg of source at (%5.2f,%5.2f) with TS=%5.0f'\
              % (cdir.ra(), cdir.dec(), center.TS() )
        olddiff = PointSourceLikelihood.get_diffuse()
        back = CompositeSkySpectrum(gal, exposure)
        back.add(center)
        PointSourceLikelihood.set_diffuse(back)

        for i in range(len(self.TS)):
            like = self.like[i]
            rdir = like.dir()
            diff = math.degrees(cdir.difference(rdir))
            if diff>0.05 and diff <1 :
                print 'refitting  source at (%5.2f,%5.2f): TS= %5.0f'% (rdir.ra(), rdir.dec(),self.TS[i] ),
                like.setDir(rdir) #does a refit
                like.maximize()
                print "-->%5.0f" % like.TS()
                
        PointSourceLikelihood.set_diffuse(olddiff)
                

print 'On first time, run pixeldata(), diffuse '  
    
    
    
