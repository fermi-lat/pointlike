"""
format:

    pointfit [parameters] <datafile> <sourcelist> [<outputfile>]

where:
    <datafile>: Name of a file, including glob-style wild cards, of photon data in FITS
        format or a merit ROOT tuple, or a FITS file with pixel data. If preceded
        by an at-sign, the remainder is the name of an ascii file containing a list of
        files to open
    <sourcelist>: a file, ascii or FITS, containing a list of seed points for fit.
        Expected fields name, ra, dec. There should be a field (name TBD) with TS
        values allowing a ranking of significance or flux. pointfit will use this
        to sort nearby pairs. If not found, default behavior is to not sort.
        [Not currently working, will be fixed]
    <outputfile>: Optional parameter specifying file name to write results, in
        format: name, input_ra, input_dec, fit_ra, fit_dec TS, sigma, pull, [name of background source(s)]
        where sigma is the fit error circle radius in degrees: >90 indicates a convergence error
              pull is the difference between fit and seed position in units of sigma

Optional parameters:
    --exposure=: Either the name of a gtexpcube-generated file, or a value, in cm^2 s,
        to apply to all energies and directions. (3e10 is appropriate for a year.)
        If not specified, use the GTI information from the datafile.
        [This is TODO: Is there a function member to sum the GTI intervals?]
        Note that this is only actually needed for overlapped sources in the Galactic
        plane, if spectral information is not required.
    --diffuse=:  Define a diffuse file (see the flag --galdiffuse to define it in the context of the science tools)
    --galdiffuse: Flag to use the galprop-generated galactic diffuse file which is distributed with the science tools)
    --eventtype= [-1] Event selection if datafile is event data. -1 means front and back,
        0/1 for front/back specifically.
    --write=<output>: if set, and the datafile is event data, write a pixelfile for
        subsequent input
    -v or --verbose [0] set verbosity
    --binsperdecade [0] default is 2.35 ratio. otherwise energy binning is set.
    --emin [500]  minimum energy, used to select bands


 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit.py,v 1.11 2008/08/21 18:59:30 burnett Exp $
"""
import os, sys, types
from numpy import arange

from pointlike import DiffuseFunction
from pointlike import SourceList, Source, PointSourceLikelihood, Background, Data

#----------------------------------------------------------------------------------------
class Fitter(object):
    """ this is an attempt to be consistent with an older interface that is used by ASP
    """
    def __init__(self, source, data,  background=None, skip=2, verbose=0, localize=True):
        psl = PointSourceLikelihood(data.map(), source.name(), source.dir())
        if background is not None: PointSourceLikelihood.set_diffuse(background)
        psl.set_verbose(verbose)
        self.TS = psl.TS()
        self.name = source.name() 
        if localize:
            self.sigma = psl.localize(skip) 
            self.delta = math.degrees(psl.dir().difference(source.dir()))
        else:
            self.sigma = self.delta=-1
        self.ra = psl.dir().ra()
        self.dec= psl.dir().dec()
        self.sdir = psl.dir()
        self.like = psl # access to the likelihood functions

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
        #print filelist
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
        if hd[1].name=='PHOTONMAP' or hd[1].name=='BANDS':
            hd.close()
            if len(files)>1: print 'Warning: more than one photonmap file not supported'
            data = Data(files[0], hd[1].name)
        else:
            hd.close()
            data = Data(files, eventtype, tstart , tstop)
    else:
        raise Exception('filename %s not a valid list of files or fits file' % filename)
    if pixeloutput is not None:
        data.map().write(pixeloutput)
        print 'created a photonmap file: %s' % pixeloutput
    return data
#--------------------------------------------------------
    
def set_diffuse(diffusefilename='galdiffuse', exposure=1e10):
    if 'GLAST_EXT' in os.environ and diffusefilename=='galdiffuse':
        diffusefilename = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'GP_gamma.fits')
    diffuse = DiffuseFunction(diffusefilename)
    background = Background(diffuse, exposure)
    PointSourceLikelihood.set_diffuse(background)
    return (diffuse, background) # return to keep references
    
#--------------------------------------------------------

def main():

    def help(msg=None):
        if msg: print '\nError: %s' % msg
        print __doc__
        sys.exit(0)
    
    from getopt import getopt, GetoptError
    import sys


    options = 'b:w:v'
    long_options= [ 'diffuse=','write=', 'verbose', 'galdiffuse', 
                    'eventtype=', 'exposure=', 'binsperdecade=', 'emin=']

    try:    
        (opts, args) = getopt(sys.argv[1:], options, long_options )
    except GetoptError, msg:
        help(msg)

    outputpixelfile= background=None
    diffusefilename='galdiffuse' # wire in for now
    verbose=0
    exposure=3e10 # this is appropriate for 1 year. 
    eventtype=-1  # all events
    binsperdecade=0 # default binning
    emin = 500
    TSmin= 10  # minimum TS for final list
                                    
    for (opt,val) in opts:
        if   opt=='-b' or opt=='--diffuse'  : diffusefilename = val
        elif opt=='-w' or opt=='--write'    : outputpixelfile = val
        elif opt=='-v' or opt=='--verbose'  : verbose =1
        elif opt=='--galdiffuse'            : diffusefilename='galdiffuse' #flag
        elif opt=='--eventtype'             : eventtype= int(val)
        elif opt=='--emin'                  : emin = float(val)
        elif opt=='--binsperdecade'         : binsperdecade=float(val)
        elif opt=='--exposure'              :
            try: exposure= float(val)
            except: exposure = val

    if len(args)>0: eventfilename=args[0]
    else: help('No event file name specified')
    if binsperdecade>0:
        bins = 10**arange(2,5,1./binsperdecade)
        Data.setEnergyBins(bins)
    
    data = photonmap(eventfilename, pixeloutput=outputpixelfile, eventtype=eventtype)
    data.info()

    if len(args)>1: sourcefilename= args[1]
    else: help('No source file list')

    if diffusefilename is not None:
        print 'setting up background from file %s' % diffusefilename
        (t1,t2)=set_diffuse(diffusefilename, exposure)
    PointSourceLikelihood.set_energy_range(emin)
    PointSourceLikelihood.set_verbose(verbose)
    SourceList.set_data(data.map())
    sourcelist = SourceList(sourcefilename)
    sourcelist.sort_TS()
    sourcelist.refit()
    sourcelist.filter_TS(TSmin)
    sourcelist.sort_ra()
    sourcelist.dump()
    if len(args)>2:
        sourcelist.dump(args[2])


#--------------------------------------------------------
    
if __name__=='__main__':
    main()
    
