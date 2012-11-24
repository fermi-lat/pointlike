"""
Python module to create a stream, using current folder
$Header$

"""
import os, sys, argparse

def createStream(stage, sdir, job_list):
    pipeline='/afs/slac/g/glast/ground/bin/pipeline -m PROD createStream '
    cmd=pipeline+' -D "stage=%s, SKYMODEL_SUBDIR=%s, job_list=%s" UWpipeline' %(stage,sdir,job_list)
    return cmd

def main( args ):
    assert os.path.exists('config.txt'), 'expect this folder to have a file config.txt'
    cwd = os.getcwd()
    m = cwd.find('skymodels')
    assert m>0, 'did not find "skymodels" in path to cwd, which is %s' %cwd
    sdir = cwd[m:]
    for stage in args.stage:
        cmd =createStream(stage, sdir, args.joblist)
        print cmd
        if not args.test:
            os.system(cmd)
   
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Submit one or more UWpipeline streams; must be in skymodel folder')
    parser.add_argument('stage', nargs='+', help='the stage indentifier(s)')
    parser.add_argument('-j','--joblist',  help='Optional list of jobs; assume local to $POINTLIKE_DIR', default='job_list')
    parser.add_argument('--test', action='store_true', help='Do not run the pipeline createStream')
    args = parser.parse_args()
    main(args)
