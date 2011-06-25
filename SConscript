# -*- python -*-
# @file SConscript
# @brief scons build specifications
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/SConscript,v 1.142 2011/06/22 19:06:41 lande Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: pointlike-07-07-28

import os

#specify package name, applications
package= 'pointlike'
libname = package+'Lib'
testname = 'test_'+package
apps   =['pointfit', 'pointfind', 'alignment']

# this part is standard: assume includes, a shareable lib, zero or more applications, a test program
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package=package, toBuild='shared')

progEnv.Tool(libname)

lib = libEnv.SharedLibrary(package, listFiles(['src/*.cxx']))

swigEnv = progEnv.Clone()
pyLib = swigEnv.SwigLibrary('_'+package,'src/swig_setup.i')

progEnv.Tool('registerTargets', 
             package   = package, 
             includes  = listFiles([package+'/*.h']),
             libraryCxts = [[lib, libEnv]],
             swigLibraryCxts = [[pyLib, swigEnv]],
             binaryCxts  = [[progEnv.Program(name, listFiles(['src/%s/*.cxx'%name])), progEnv] for name in apps], 
             testAppCxts  = [[progEnv.Program(testname, listFiles(['src/test/*.cxx'])), progEnv]],
	     data = ['data/test_events.root'],
             python = (['src/pointlike.py','python/pointlike_defaults.py',
			'python/pointfit.py', 'python/pointfit_setup.py', 
                        'python/test_pointlike_setup.py']+listFiles(['python/uw'],recursive=True)))
