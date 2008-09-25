# @file pointlikeLib.py
# @brief scons package dependencies
#
#$Id: pointlikeLib.py,v 1.4 2008/02/26 05:25:20 glastrm Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['pointlike'])
    depends = 'astro healpix skymaps embed_python'.split()
    for pack in depends: env.Tool(pack+'Lib')
    env.Tool('addLibrary', library=env['minuitLibs'])
    
def exists(env):
    return 1
