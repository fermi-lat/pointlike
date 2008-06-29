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
    <outputfile>: Optional parameter specifying file name to write results, in
        format: name, input_ra, input_dec, fit_ra, fit_dec TS, sigma, diff, [name of background source(s)]
        where sigma is the fit error circle radius: >90 indicates a convergence error
              diff is the difference between fit and initial position

Optional parameters:
    --exposure=: Either the name of a gtexpmap-generated file, or a value, in cm^2 s,
        to apply to all energies and directions. (3e10 is appropriate for a year.)
        If not specified, use the GTI information from the datafile.
        [This is TODO: Is there a function member to sum the GTI intervals?]
        Note that this is only actually needed for overlapped sources in the Galactic
        plane, if spectral information is not required.
    --galdiffuse: Flag to use the galprop-generated galactic diffuse file
    --minlevel= [8] Minimum healpix level for fit [question: skip or not?]
    --eventtype= [-1] Event selection if datafile is event data. -1 means front and back,
        0/1 for front/back specifically.
    --write=<output>: if set, and the datafile is event data, write a pixelfile for
        subsequent input
    -v or --verbose [0] set verbosity


 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit.py,v 1.5 2008/04/09 13:02:18 burnett Exp $
"""
# setup to import pointlike
try: import uw.pointlike
except: pass

import os, sys, types

from fitter import Fitter,  photonmap, sourcelist
#from pointlike import DiffuseFunction, CompositeSkySpectrum
from background import Background



#--------------------------------------------------------

def main():

    def help(msg=None):
        if msg: print '\nError: %s' % msg
        print __doc__
        sys.exit(0)
    
    from getopt import getopt, GetoptError
    import sys


    options = 'b:w:v'
    long_options= [ 'diffuse=','write=', 'verbose', 'galdiffuse', 'minlevel=',
                    'eventtype=', 'exposure=']

    try:    
        (opts, args) = getopt(sys.argv[1:], options, long_options )
    except GetoptError, msg:
        help(msg)

    outputpixelfile= diffusefilename=background=None
    verbose=0
    exposure=3e10 # this is appropriate for 1 year. 
    minlevel=8
    eventtype=-1  # all events
    radius = 8.7  
                                    
    for (opt,val) in opts:
        if   opt=='-b' or opt=='--diffuse'  : diffusefilename = val
        elif opt=='-w' or opt=='--write'    : outputpixelfile = val
        elif opt=='-v' or opt=='--verbose'  : verbose =1
        elif opt=='--galdiffuse'            : diffusefilename='galdiffuse' #flag
        elif opt=='--minlevel'              : minlevel = int(val)
        elif opt=='--eventtype'             : eventtype= int(val)
        elif opt=='--exposure'              :
            try: exposure= float(val)
            except: exposure = val

    if len(args)>0: eventfilename=args[0]
    else: help('No event file name specified')
    
    data = photonmap(eventfilename, pixeloutput=outputpixelfile, eventtype=eventtype)

    if len(args)>1: sourcefilename= args[1]
    else: help('No source file list')

    background = Background(diffusefilename, exposure)

    sources = sourcelist(sourcefilename)

    out = None if len(args)==2 else file(args[2], 'w')
    for source in sources:
        fit = Fitter(source, data , background=background(), verbose=verbose)
        print >>out, ('%-20s'+'%10.4f'*4+                             '%10.2f'+'%10.4f'*2)\
              %( source.name, source.ra, source.dec, fit.ra, fit.dec, fit.TS, fit.sigma, fit.delta)

#--------------------------------------------------------
    
if __name__=='__main__':
    print 'running pointfit'
    main()
    
