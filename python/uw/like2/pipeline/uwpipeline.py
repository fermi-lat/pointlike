"""
task UWpipeline Interface to the ISOC PipelineII

$Header$
"""
import os, argparse
import numpy as np
#from . import pipe
from uw.like2.pipeline import check_data
from uw.like2.pipeline import pipeline_job
from uw.like2.pipeline import summary_plots
from uw.like2.pipeline import check_converge

class Stage(dict):
    def __init__(self, help='', **kwargs):
        super(Stage,self).__init__(help=help, **kwargs)
    pass

# List of possible stages
# not coordinated with redundant lists in pipeline_job, summary_plots, check_converge
    
stagenames = dict(
    create=Stage( help='Create a new skymodel'),
    update_full=Stage(help='perform update'),
    update = Stage(help='perform update'),
    update_beta = Stage(help='update beta'),
    update_pivot = Stage(help='update pivot'),
    finish = Stage(help='perform localization'),
    tables = Stage(help='create tables'),
    sedinfo =Stage(help=''),
    diffuse=Stage(help=''),
    isodiffuse=Stage(help=''),
    limb=Stage(help=''),
    fluxcorr=Stage(help=''),
    pulsar_table=Stage(help=''),
)

keys = stagenames.keys()
stage_help = 'stage name, or sequential stages separaged by : must be one of %s' %keys

class StartStream(object):
    """ setup, start a stream """
    def main(self, args):
        pipeline='/afs/slac/g/glast/ground/bin/pipeline -m PROD createStream '
        for stage in args.stage:
            cmd=pipeline+' -D "stage=%s, SKYMODEL_SUBDIR=%s, job_list=%s" UWpipeline' \
                %(stage, args.skymodel,args.job_list)
            print '-->' , cmd
            if not args.test:
                os.system(cmd)

startup = StartStream()

class Proc(dict):
    def __init__(self, run=None, help='', **kwargs):
        super(Proc, self).__init__(self,  help=help, **kwargs)
        self.run = run
    def __call__(self, args):
        if self.run is not None:
            self.run.main(args)
        else:
            raise Exception( 'no code to run for proc %s'% self)
 
procnames = dict(
    start = Proc(startup, help='start a stream'),
    check_data = Proc(check_data, help='check that required data files are present'),
    job   = Proc(pipeline_job, help='run a parallel pipeline job'),
    check_converge = Proc(check_converge, help='check for convergence, combine results, possibly submit new stream'),
    summary =Proc(summary_plots, help='Process summaries, need stage'),
    )
    
def check_names(stage, proc):
    if len(stage)==0:
        if proc is  None:
            raise Exception('No proc or stage argement specified')
        if proc not in procnames:
            raise Exception('proc name "%s" not in list %s' % (proc, procnames,keys()))
        return
    if stage[0] is None: 
        raise Exception('no stage specified')
    for s in stage:
        for t in s.split(':'):
            if t not in keys:
                raise Exception('"%s" not found in possible stage names, %s' %(t, keys))

def check_environment(args):
    assert os.path.exists('config.txt'), 'expect this folder to have a file config.txt'
    cwd = os.getcwd()
    m = cwd.find('skymodels')
    assert m>0, 'did not find "skymodels" in path to cwd, which is %s' %cwd
    if args.stage[0] is None :
        pass #    raise Exception( 'No stage specified: either command line or os.environ')
    else:
        os.environ['stage']=args.stage[0]

    if 'SKYMODEL_SUBDIR' not in os.environ:
        os.environ['SKYMODEL_SUBDIR'] = cwd
    # add these to the Namespace object for convenience
    args.__dict__.update(skymodel=cwd, pointlike_dir=cwd[:m], stream=-1)


def main( args ):
    check_environment(args)
    check_names(args.stage, args.proc)
    proc = args.proc
    print '-->', proc
    procnames[proc](args)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=' start a UWpipeline stream, or run a proc; ')
    parser.add_argument('stage', nargs='*', default=[os.environ.get('stage', None)], help=stage_help)
    parser.add_argument('-p', '--proc', default=os.environ.get('PIPELINE_PROCESS', 'start'),
        help='proc name as defined by the UWpipeline xml, one of: %s'%procnames.keys())
    parser.add_argument('--job_list', default=os.environ.get('job_list', 'joblist.txt'), help='file used to allocate jobs')
    parser.add_argument('--test', action='store_true', help='Do not run' )
    args = parser.parse_args()
    main(args)
    
 
