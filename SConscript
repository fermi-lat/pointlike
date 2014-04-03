# -*- python -*-
# @file SConscript
# @brief scons build specifications
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/SConscript,v 1.185 2014/01/27 03:04:43 burnett Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: pointlike-08-00-04

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
             ###testAppCxts  = [[progEnv.Program(testname, listFiles(['src/test/*.cxx'])), progEnv]],
	     data = ['data/test_events.root'],
             python = (['src/pointlike.py','python/pointlike_defaults.py',
			'python/pointfit.py', 'python/pointfit_setup.py', 
                        'python/test_pointlike_setup.py']+listFiles(['python/uw'],recursive=True)))
