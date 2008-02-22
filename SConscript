# @file SConscript
# @brief build info
#
# $Header$

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

progEnv.Tool('pointlikeLib')

progEnv.Tool('registerObjects', 
    package  = 'pointlike', 
    includes = listFiles(['pointlike/*.h']),
    libraries= [
        libEnv.SharedLibrary('pointlike', listFiles(['src/*.cxx']))], 
    binaries = [
        progEnv.Program('pointfit',  listFiles(['src/pointfit/*.cxx'])),
        progEnv.Program('pointfind', listFiles(['src/pointfind/*.cxx'])),
	progEnv.Program('alignment', listFiles(['src/alignment/*.cxx'])),
    ], 
    testApps = [progEnv.Program('test_pointlike', listFiles(['src/test/*.cxx']))],
 )

