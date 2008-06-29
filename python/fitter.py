"""
 Utility classes or functions to implement pointlike 
 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/fitter.py,v 1.5 2008/06/10 17:57:01 burnett Exp $
"""
try: import uw.pointlike
except: pass

from pointlike import SkyDir, PointSourceLikelihood, Data, DiffuseFunction, CompositeSkySpectrum
import os, sys, math

#----------------------------------------------------------------------------------------
class Source(object):
    def __init__(self, *args):
        if len(args)==1:
            tokens = args[0].split()
            self.name = tokens[0]
            try:
                self.ra, self.dec = [float(tok) for tok in tokens[1:3]]
            except:
                print 'invalid source line: "%s"' %args[0]
                raise
        else:
            self.name = args[0]
            self.ra, self.dec = [float(a) for a in args[1:3]]

        self.sdir = SkyDir(self.ra, self.dec)
#----------------------------------------------------------------------------------------

class Fitter(object):
    
    def __init__(self, source, data,  background=None, skip=2, verbose=0, localize=True):
        psl = PointSourceLikelihood(data.map(), source.name, source.sdir)
        if background is not None: psl.set_diffuse(background)
        psl.set_verbose(verbose)
        self.TS = psl.maximize()
        self.name = source.name 
        if localize:
            self.sigma = psl.localize(skip) 
            self.delta = math.degrees(psl.dir().difference(source.sdir))
        else:
            self.sigma = self.delta=-1
        self.ra = psl.dir().ra()
        self.dec= psl.dir().dec()
        self.sdir = psl.dir()
        self.like = psl # access to the likelihood functions
        if self.sigma< 1. : #test for convergence of localization, or not done
            if verbose>0: psl.printSpectrum()
            self.photons = [psl[i].photons() for i in range(len(psl))]
            self.alpha   = [psl[i].alpha() for i in range(len(psl))]

#----------------------------------------------------------------------------------------

def photonmap(filename, eventtype=-1, pixeloutput=None, tstart=0, tstop=0):
    """ return a Data object, determined one of 3 ways:
        * the name of a file containing a list of photon files, preceded by an @
        * the name of a photon file, that will be expanded by glob
        * the name of a photon map file.
        if pixeloutput is set, write the photonmap to the file
        eventtype is -1 for all events, 0/1 for front back
    """
    data = None
    if filename[0]=='@':
        # it is a list of data files
        filelist = [line.strip() for line in file(filename[1:]) if len(line)>0 and line[0]!='#']
        print filelist
        data =  Data(filelist, eventtype, tstart, tstop)
    elif filename[-5:]=='.fits' or filename[-4:]=='.fit' :
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
            data = Data(files, eventtype, tstart , tstop)
    else:
        raise Exception('filename %s not a valid list of files or fits file' % filename)
    if pixeloutput is not None:
        data.map().write(pixeloutput)
        print 'created a photonmap file: %s' % pixeloutput
    return data
#----------------------------------------------------------------------------------------
    
def sourcelist(filename):
    if filename[-5:]=='.fits':
        import pyfits
        data = pyfits.open(filename)[1].data
        names = data.field('NAME')
        ras = data.field('RA')
        decs= data.field('DEC')
        return [Source(names[i], ras[i], decs[i]) for i in range(len(names))]
           
    f = open(filename)
    return [Source(line) for line in f if line[0]!='#']

#----------------------------------------------------------------------------------------

if __name__=='__main__':
    print 'testing Background with galdiffuse option'
    t = Background('galdiffuse')
    pass
