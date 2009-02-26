# -*- python -*-
# @file SConscript
# @brief scons build specifications
#
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/SConscript,v 1.50 2008/12/30 20:31:24 glastrm Exp $
# Authors: Toby Burnett <tburnett@u.washington.edu>
# Version: pointlike-06-16-00

#specify package name, applications
package= 'pointlike'
apps   =['pointfit', 'pointfind', 'alignment']

# this part is standard: assume includes, a shareable lib, zero or more applications, a test program
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

progEnv.Tool(package+'Lib')
libEnv.Tool(package+'Lib', depsOnly = 1)
lib = libEnv.SharedLibrary(package, listFiles(['src/*.cxx']))

swigEnv.Replace(SHLIBPREFIX = '_')
#swigEnv.Replace(SHLIBSUFFIX = '.pyd')
swigEnv.Append(RPATH = swigEnv['LIBDIR'])
pyLib = swigEnv.SharedLibrary(package,'python/swig_setup.i')

progEnv.Tool('registerObjects', 
    package  = package, 
    includes = listFiles([package+'/*.h']),
    libraries= [lib], 
    binaries = [progEnv.Program(name, listFiles(['src/%s/*.cxx'%name])) for name in apps], 
    testApps = [progEnv.Program('test_'+package, listFiles(['src/test/*.cxx']))],
    python = ['python/%s.py'%package, pyLib],
 )
