"""
format:

    specfit [parameters] <datafile> <sourcelist> [<outputfile>]

where:
    <datafile>: Name of a file, including glob-style wild cards, of photon data in FITS
        format or a merit ROOT tuple, or a FITS file with pixel data. If preceded
        by an at-sign, the remainder is the name of an ascii file containing a list of
        files to open
    <sourcelist>: a file, ascii or FITS, containing a list of seed points for fit.
        Expected fields name, ra, dec, [TS]. There should be a field (name TBD) with TS
        values allowing a ranking of significance or flux. pointfit will use this
        to sort nearby pairs. If not found, default behavior is to not sort.
    <outputfile>: Optional parameter specifying file name to write results, in
        format: name, input_ra, input_dec, fit_ra, fit_dec TS, sigma, diff,
                [name of background source(s)], list of values of E*dN/dE at energy band centers
        where sigma is the fit error circle radius: >90 indicates a convergence error
              diff is the difference between fit and initial position
        The model parameters are listed as ModelName:: p1 (err) p2 (err) etc.

Optional parameters:
    --exposure=: Either the name of a gtexpmap-generated file, or a value, in cm^2 s,
        to apply to all energies and directions. (3e10 is appropriate for a year.)
        The name is interpreted as a root to which the suffices '_front.fits' and '_back.fits'
        will be appended.  Both files need not exist.
        If not specified, use the GTI information from the datafile.
        [This is TODO: Is there a function member to sum the GTI intervals?]
        Note that this is only actually needed for overlapped sources in the Galactic
        plane, if spectral information is not required.
    --galdiffuse: Flag to use the galprop-generated galactic diffuse file
    --eventtype= [-1] Event selection if datafile is event data. -1 means front and back,
        0/1 for front/back specifically.
    --write=<output>: if set, and the datafile is event data, write a pixelfile for
        subsequent input
    -v or --verbose [0] set verbosity
    --emin - left edge of minimum energy bin (default 100 MeV)
    --emax - right edge of maximum energy bin (default 2.51e5 MeV)
    --enumbins - (default 17)
    --plotpath - if set, generate plots with the spectral energy density and power law
        spectral fit, of form plotpath+name+_sed.png.
    --model - specify a model type, e.g. PowerLaw (default), ExpCutoff... see Models.py
    --printspec - print out spectral fit results
    --fitter - which spectral fitter to use -- default Marginal Poisson


 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/specfit.py,v 1.2 2008/07/25 18:27:43 kerrm Exp $
"""
# setup to import pointlike
try: #Try block only for UW environment
   import uw.pointlike
except: pass

import pointlike as pl

import os, sys, types
from numpy import arange

import pointlike as pl
from SourceLib import *
from Response import *

#from pointlike import DiffuseFunction

#from pointlike import SourceList, Source, PointSourceLikelihood, Background, Data
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
        data =  pl.Data(filelist, eventtype, tstart, tstop)
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
            data = pl.Data(files[0], hd[1].name)
        else:
            hd.close()
            data = pl.Data(files, eventtype, tstart , tstop)
    else:
        raise Exception('filename %s not a valid list of files or fits file' % filename)
    if pixeloutput is not None:
        data.map().write(pixeloutput)
        print 'created a photonmap file: %s' % pixeloutput
    return data

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
                    'eventtype=', 'exposure=', 'binsperdecade=','emin='
                    'enumbins=','plotpath=','fitter=','printspec=','model=']

    try:    
        (opts, args) = getopt(sys.argv[1:], options, long_options )
    except GetoptError, msg:
        help(msg)

    outputpixelfile= background=plopath=None
    diffusefilename='galdiffuse' # wire in for now
    verbose=0
    exposure=3e10 # this is appropriate for 1 year. 
    eventtype=-1  # all events
    binsperdecade=5. # default binning
    emin=300.
    enumbins=11. #5 per decade with defaults (300-30000 MeV)
    fitter='MP'
    emap=ExposureMap(const_val=3e10)
    plotpath=None
    printspec=True
    model='PowerLaw'
                                    
    for (opt,val) in opts:
        if   opt=='-b' or opt=='--diffuse'  : diffusefilename = val
        elif opt=='-w' or opt=='--write'    : outputpixelfile = val
        elif opt=='-v' or opt=='--verbose'  : verbose =1
        elif opt=='--galdiffuse'            : diffusefilename='galdiffuse' #flag
        elif opt=='--eventtype'             : eventtype= int(val)
        elif opt=='--binsperdecade'         : binsperdecade=float(val)
        elif opt=='--emin'                  : emin = float(val)
        elif opt=='--enumbins'              : enumbins = float(val)
        elif opt=='--plotpath'              : plotpath = str(val)
        elif opt=='--fitter'                : fitter = str(val)
        elif opt=='--printspec'             : printspec = True #flag
        elif opt=='--model'                 : model = str(model)
        elif opt=='--exposure'              :
            try:
               exposure= float(val)
               emap=ExposureMap(const_val=exposure)
            except: #file name provided?
               exposure=str(val)
               print '\n\n\n'+exposure
               emap=ExposureMap(emap_file=exposure+'_front.fits',emap_file2=exposure+'_back.fits')

    if len(args)>0: eventfilename=args[0]
    else: help('No event file name specified')
    if binsperdecade>0: #Custom binning
       bins = emin*10**(N.arange(enumbins+1)/binsperdecade)
       pl.Data.setEnergyBins(bins)
    else:
      bins=100*2.35**N.arange(11)      
    bands=EnergyBands(bins[:-1],[bins[-1]])
    
    data = photonmap(eventfilename, pixeloutput=outputpixelfile, eventtype=eventtype)
    data.info()

    if len(args)>1: sourcefilename= args[1]
    else: help('No source file list')

    if diffusefilename is not None:
        print 'setting up background from file %s' % diffusefilename
        if 'GLAST_EXT' in os.environ and diffusefilename=='galdiffuse':
            diffusefilename = os.path.join(os.environ['GLAST_EXT'],'extFiles','v0r7','galdiffuse', 'GP_gamma.fits')

        diffuse = pl.DiffuseFunction(diffusefilename)
        if (not emap.const_fe) and (not emap.const_be): #Both front and back exposure set
            background = pl.Background(diffuse, emap.fe, emap.be)
        else:
            background = pl.Background(diffuse, exposure)
        pl.PointSourceLikelihood.set_diffuse(background)

    pl.PointSourceLikelihood.set_energy_range(bins[0]*1.01)
    pl.SourceList.set_data(data.map())
    sourcelist = pl.SourceList(sourcefilename)
    sourcelist.sort_TS()
    sourcelist.refit()
    sourcelist.dump()
    if len(args)>2:
        sourcelist.dump(args[2])
        if plotpath is None and not printspec: return
    print 'Finished processing C++ SourceList.'
    
    
    #Now use the SourceList to construct Source objects
    sl=SourceLibrary()
    response=ModelResponse(bands,emap)
    global_data=GlobalData(emap,response,method=fitter)
    sl.add_global_data(global_data)
    sl.add_sourcelist(sourcelist)
    for s in sl:
        sl.fit(s().name(),model,savepath=plotpath,printfit=printspec)

    
   
#run ../../spectrum_dev4/specfit --exposure=P6_V1_DIFFUSE_100b_30_300000 --printspec=True ft1_first_diff_v2.fits d:/common/first_light/associated_sources.txt

#--------------------------------------------------------
    
if __name__=='__main__':
    main()
    
