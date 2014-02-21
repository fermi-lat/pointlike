# @file pointlikeLib.py
# @brief scons package dependencies
#
#$Id: pointlikeLib.py,v 1.6 2009/09/19 00:39:14 wallacee Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['pointlike'])
    env.Tool('addLibrary', library=env['rootLibs'])
    depends = 'astro healpix skymaps embed_python'.split()
    for pack in depends: env.Tool(pack+'Lib')
    env.Tool('addLibrary', library=env['clhepLibs'])
    
def exists(env):
    return 1
