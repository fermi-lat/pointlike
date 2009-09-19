# -*- python -*-
# @file SConscript
# @brief scons build specifications
#
# $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/SConscript,v 1.75 2009/09/11 19:26:43 lande Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: pointlike-06-23-04

import os

#specify package name, applications
package= 'pointlike'
apps   =['pointfit', 'pointfind', 'alignment']

# this part is standard: assume includes, a shareable lib, zero or more applications, a test program
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

progEnv.Tool(package+'Lib')
libEnv.Tool(package+'Lib', depsOnly = 1)
lib = libEnv.SharedLibrary(package, listFiles(['src/*.cxx']))

swigEnv = progEnv.Clone()

swigEnv.Append(RPATH = [swigEnv['LIBDIR'], 
                        os.path.join(os.environ['CLHEPROOT'],'lib')])
pyLib = swigEnv.LoadableModule('_pointlike','python/swig_setup.i')

progEnv.Tool('registerTargets', 
             package   = package, 
             includes  = listFiles([package+'/*.h']),
             libraryCxts = [[lib, libEnv]],
             swigLibraryCxts = [[pyLib, swigEnv]],
             binaryCxts  = [[progEnv.Program(name, listFiles(['src/%s/*.cxx'%name])), progEnv] for name in apps], 
             testAppCxts  = [[progEnv.Program('test_'+package, listFiles(['src/test/*.cxx'])), progEnv]],
             python = ['python/pointlike.py','python/pointlike_defaults.py',
                       'python/test_pointlike_setup.py'])


