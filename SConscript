# -*- python -*-
# @file SConscript
# @brief scons build specifications
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/SConscript,v 1.100 2010/04/12 18:46:31 burnett Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: pointlike-07-04-01

import os

#specify package name, applications
package= 'pointlike'
apps   =['pointfit', 'pointfind', 'alignment']

# this part is standard: assume includes, a shareable lib, zero or more applications, a test program
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package=package, toBuild='shared')

progEnv.Tool(package+'Lib')

lib = libEnv.SharedLibrary(package, listFiles(['src/*.cxx']))

swigEnv = progEnv.Clone()
pyLib = swigEnv.SwigLibrary('_pointlike','src/swig_setup.i')

progEnv.Tool('registerTargets', 
             package   = package, 
             includes  = listFiles([package+'/*.h']),
             libraryCxts = [[lib, libEnv]],
             swigLibraryCxts = [[pyLib, swigEnv]],
             binaryCxts  = [[progEnv.Program(name, listFiles(['src/%s/*.cxx'%name])), progEnv] for name in apps], 
             testAppCxts  = [[progEnv.Program('test_'+package, listFiles(['src/test/*.cxx'])), progEnv]],
             python = (['src/pointlike.py','python/pointlike_defaults.py',
                       'python/test_pointlike_setup.py']+listFiles(['python/uw'],recursive=True)))


