"""
 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit.py,v 1.2 2008/03/21 02:48:22 burnett Exp $
"""
try: import uw.pointlike
except: pass

import os, sys

from fitter import Fitter, Background, photonmap, sourcelist
        

def main():
    
    from getopt import getopt, GetoptError
    import sys

    def help(msg=None):
        if msg: print '\nError: %s' % msg
        sys.exit(0)

    options = 'b:w:v'
    long_options= [ 'diffuse=','write=', 'verbose', 'galdiffuse']

    try:    
        (opts, args) = getopt(sys.argv[1:], options, long_options )
    except GetoptError, msg:
        help(msg)

    outputpixelfile= diffusefilename=background=None
    verbose=0
    exposure=3e10 # this is appropriate for 1 year. Todo: make it settable, from an exposure file
                                    
    for (opt,val) in opts:
        if   opt=='-b' or opt=='--diffuse' : diffusefilename = val
        if   opt=='-w' or opt=='--write'   : outputpixelfile = val
        if   opt=='-v' or opt=='--verbose' : verbose =1
        if   opt=='--galdiffuse'           : diffusefilename='galdiffuse' #flag 

    if len(args)>0: eventfilename=args[0]
    else: help('No event file name specified')
    
    data = photonmap(eventfilename, pixeloutput=outputpixelfile)

    if len(args)>1: sourcefilename= args[1]
    else: help('No source file list')

    background = Background(diffusefilename, exposure)

    sources = sourcelist(sourcefilename)

    out = None if len(args)==2 else file(args[2], 'w')
    for source in sources:
        fit = Fitter(source, data , background=background(), verbose=verbose)
        print >>out, ('%-20s'+'%10.4f'*4+                             '%10.2f'+'%10.4f'*2)\
              %( source.name, source.ra, source.dec, fit.ra, fit.dec, fit.TS, fit.sigma, fit.delta)
    

if __name__=='__main__':
    main()
