# @file pointlikeLib.py
# @brief scons package dependencies
#
#$Id: pointlikeLib.py,v 1.2 2008/02/22 15:07:08 burnett Exp $
def generate(env, **kw):
    env.Tool('addLibrary', library = ['pointlike'])
    depends = 'astro healpix skymaps embed_python'.split()
    for pack in depends: env.Tool(pack+'Lib')
    
def exists(env):
    return 1
