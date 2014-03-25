"""
Run after a successful UWpipeline/job_task

Summarize the execution,

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/summary_plots.py,v 1.3 2012/12/23 20:19:11 burnett Exp $
"""
import os, sys,  argparse

from uw.like2.pipeline import pipe, processor, diagnostic_plots


def main(args):
    """ 
    """
    pointlike_dir=args.pointlike_dir 
    skymodel  =args.skymodel 
    stream    =args.stream
    stagelist = args.stage[0] 
   
    absskymodel = os.path.join(pointlike_dir, skymodel)

    os.chdir(absskymodel) 
    t = stagelist.split(':',1)
    if len(t)==2:
        stage, nextstage = t 
    else: stage,nextstage = t[0], None

    if stage.split('_')[0]=='update' or stage=='counts':
        diagnostic_plots.main('counts');

    elif stage=='sedinfo':
        diagnostic_plots.main('fb')

    elif stage=='create' or stage=='create_reloc':
        pass
    elif stage=='diffuse':
        diagnostic_plots.main('gal')

    elif stage=='isodiffuse':
        diagnostic_plots.main('iso')

    elif stage=='limb':
        pass
    elif stage=='finish':
        diagnostic_plots.main('sources')
        diagnostic_plots.main('sourceinfo') # join thise 

    elif stage=='tables':
        pass
    else:
        print 'stage %s not recognized for summary'%stage 

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run after a set of pipeline jobs, check status, accumulate info')
    parser.add_argument('--stage', default=os.environ.get('stage', '?'), help='the stage indentifier(s)')
    parser.add_argument('--pointlike_dir', default=os.environ.get('POINTLIKE_DIR', '.'),
            help='top level folder with pointlike')
    parser.add_argument('--skymodel', default= os.environ.get('SKYMODEL_SUBDIR', '.'),
            help='folder, from pointlike_dir, to the skymodel. Default $SKYMODEL_SUBDIR, set by pipeline')
    parser.add_argument('--stream', default = os.environ.get('PIPELINE_STREAM', '0'),
            help='stream number')
    parser.add_argument('--test', action='store_true', help='test mode')
    args = parser.parse_args()
    main(args)

