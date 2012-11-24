"""
#Python script to create a stream, using current folder
"""
import os, sys

def createStream(stage, sdir):
    pipeline='/afs/slac/g/glast/ground/bin/pipeline -m PROD createStream '
    cmd=pipeline+' -D "stage=%s, SKYMODEL_SUBDIR=%s" UWpipeline' %(stage,sdir)
    return cmd

assert os.path.exists('config.txt'), 'expect this folder to have a file config.txt'
cwd = os.getcwd()
m = cwd.find('skymodels')
assert m>0, 'did not find "skymodels" in path to cwd, which is %s' %cwd
sdir = cwd[m:]
assert len(sys.argv)==2, 'expect a single argument'
stage = sys.argv[1]
cmd =createStream(stage, sdir)
print cmd
os.system(cmd)