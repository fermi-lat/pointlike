# @file pointlikeLib.py
# @brief scons package dependencies
#
#$Id: pointlikeLib.py,v 1.5 2008/09/25 02:59:10 glastrm Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['pointlike'])
    depends = 'astro healpix skymaps embed_python'.split()
    for pack in depends: env.Tool(pack+'Lib')
    env.Tool('addLibrary', library=env['clhepLibs'])
    env.Tool('addLibrary', library=env['minuitLibs'])
    
def exists(env):
    return 1
