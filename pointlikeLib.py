#$Id$
def generate(env, **kw):
    env.Tool('addLibrary', library = ['pointlike'])
    env.Tool('astroLib')
    env.Tool('healpixLib')
    env.Tool('skymapsLib')
    env.Tool('embed_pythonLib')
    
def exists(env):
    return 1
