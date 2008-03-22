"""
 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit.py,v 1.2 2008/03/21 02:48:22 burnett Exp $
"""
try: import uw.pointlike
except: pass

from pointlike import SkyDir, PointSourceLikelihood, Data, DiffuseFunction, CompositeSkySpectrum
import os, sys, math

class Source(object):
    def __init__(self, *args):
        if len(args)==1:
            tokens = args[0].split()
            self.name = tokens[0]
            self.ra, self.dec = [float(tok) for tok in tokens[1:3]]
        else:
            self.name = args[0]
            self.ra, self.dec = [float(a) for a in args[1:3]]

        self.sdir = SkyDir(self.ra, self.dec)

class Fitter(object):
    
    def __init__(self, source, data,  background=None, skip=3, verbose=0):
        psl = PointSourceLikelihood(data.map(), source.name, source.sdir)
        if background is not None: psl.set_diffuse(background)
        psl.set_verbose(verbose)
        psl.maximize()
        self.sigma =psl.localize(skip)
        if verbose>0: psl.printSpectrum()
        self.TS = psl.TS()
        self.ra = psl.dir().ra()
        self.dec= psl.dir().dec()
        self.delta = math.degrees(psl.dir().difference(source.sdir))

def photonmap(filename, pixeloutput=None):
    """ return a Data object, determined one of 3 ways:
        * the name of a file containing a list of photon files, preceded by an @
        * the name of a photon file, that will be expanded by glob
        * the name of a photon map file.
        if pixeloutput is set, write the photonmap to the file
    """
    data = None
    if filename[0]=='@':
        # it is a list of data files
        filelist = [line for line in file(filename[1:]) if len(line)>0]
        data =  Data(filelist, -1, 0, 0)
    elif filename[-5:]=='.fits':
        # a fits file: either data to read, or a photonmap
        import pyfits, glob
        files = glob.glob(filename)
        if len(files)==0:
            raise Exception('no such file(s): %s' %filename)
        hd = pyfits.open(files[0])
        if len(hd)==1:
            raise Exception('Invalid data file, apparent image file with primary only')
        if hd[1].name=='PHOTONMAP':
            hd.close()
            if len(files)>1: print 'Warning: more than one photonmap file not supported'
            data = Data(files[0], 'PHOTONMAP')
        else:
            hd.close()
            data = Data(files, -1, 0 , 0)
    else:
        raise Exception('filename %s not a valid list of files or fits file' % filename)
    if pixeloutput is not None:
        data.map().write(pixeloutput)
        print 'created a photonmap file: %s' % pixeloutput
    return data
    
def sourcelist(filename):
    if filename[-5:]=='.fits':
        import pyfits
        data = pyfits.open(filename)[1].data
        names = data.field('NAME')
        ras = data.field('RA')
        decs= data.field('DEC')
        return [Source(names[i], ras[i], decs[i]) for i in range(len(names))]
           
    f = open(filename)
    return [Source(line) for line in f]
#----------------------------------------------------------------------------------------
        
class Background(object):
    """ manage the background
    """
    def __init__(self, diffusefilename, exposure=3e10):
        self.background=None
        if diffusefilename is not None:
            print 'setting up background from file %s, exposure=%.2g' % (diffusefilename, exposure)
            if 'GLAST_EXT' in os.environ and diffusefilename=='galdiffuse':
                diffusefilename = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'GP_gamma.fits')

        if diffusefilename !='':
            self._df = DiffuseFunction(diffusefilename) # need to keep this reference open!
            self.background = CompositeSkySpectrum(self._df, exposure)
    def __call__(self):
        return self.background


if __name__=='__main__':
    t = Background('galdiffuse')
    pass
