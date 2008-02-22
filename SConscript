#$Id$

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

pointlikeStaticLib = libEnv.StaticLibrary('pointlike', listFiles(['src/*.cxx']))
pointlikeSharedLib = libEnv.SharedLibrary('pointlike', listFiles(['src/*.cxx']))

progEnv.Tool('pointlikeLib')
pointfit = progEnv.Program('pointfit', listFiles(['src/pointfit/*.cxx']))
pointfind = progEnv.Program('pointfind', listFiles(['src/pointfind/*.cxx']))
test_pointlike = progEnv.Program('test_pointlike', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'pointlike', libraries = [pointlikeStaticLib, pointlikeSharedLib], testApps = [test_pointlike],
             binaries = [pointfit, pointfind], includes = listFiles(['facilities/*.h']))

