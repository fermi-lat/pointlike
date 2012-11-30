"""
Run after a successful UWpipeline/job_task

Summarize the execution, and perhaps perform specific analysis or checks,
mostly zipping the files generated by the multiple jobs

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/check_converge.py,v 1.2 2012/11/26 16:06:33 burnett Exp $
"""
import os, sys, glob, zipfile, logging, datetime, argparse
import numpy as np
import pandas as pd

from uw.like2.pipeline import pipe, processor, diagnostic_plots


def streamInfo( stream, path='.'):
    streamlogs=sorted(glob.glob(os.path.join(path,'streamlogs','stream%s.*.log'%stream)))
    substream = [ int(f.split('.')[-2]) for f in streamlogs]
    print '\t found %d streamlog files...' %(len(streamlogs),), 
    etimes =np.array([ float(open(f).read().split('\n')[-2].split()[-1][:-1]) for f in streamlogs])
    print '\texecution times: mean=%.1f, max=%.1f'% (etimes.mean(), etimes.max())
    times = pd.DataFrame(dict(etimes=etimes), index=pd.Index(substream, name='substream'))#,name='stream%d'%stream)
    bad = times[times.etimes>2*times.etimes.mean()]
    if len(bad)>0:
        print '\tsubstreams with times>2*mean:'
        print bad
    return times



def create_stream(newstage):
    cmd = 'cd %s;  /afs/slac/g/glast/ground/bin/pipeline  createStream -D "stage=%s, SKYMODEL_SUBDIR=%s" UWpipeline '\
        %(pointlike_dir, newstage, skymodel)
    rc=os.system(cmd)
    if rc==0:
        print '\n----> started new stream with stage %s'% newstage
    else: print '\n***** Failed to create new stream: tried %s'%cmd

def main(args):
    """ 
    """
    pointlike_dir=args.pointlike_dir # = os.environ.get('POINTLIKE_DIR', '.')
    skymodel=args.skymodel # = os.environ.get('SKYMODEL_SUBDIR', sys.argv[1] )
    stream=args.stream #= os.environ.get('PIPELINE_STREAM', '26')
    stagelist = args.stage #os.environ.get('stage', 'update' if len(sys.argv)==2 else sys.argv[2])
   
    absskymodel = os.path.join(pointlike_dir, skymodel)

    def make_zip(fname,  ext='pickle'):
        ff = glob.glob(os.path.join(absskymodel, fname, '*.'+ext))
        assert len(ff)>0, 'expected to find files in folder %s' %fname
        print 'found %d *.%s in folder %s ...' % ( len(ff),ext, fname,) ,
        with zipfile.ZipFile(os.path.join(absskymodel, fname+'.zip'), 'w') as z:
            for filename in ff:
                z.write( filename, os.path.join(fname,os.path.split(filename)[-1]))
        print ' zipped into file %s.zip' %fname


    if not args.test:
        tee = processor.OutputTee(os.path.join(absskymodel, 'summary_log.txt'))

    streamInfo(stream, absskymodel)

    os.chdir(absskymodel) # useful for diagnostics below
    s
    current = str(datetime.datetime.today())[:16]
    print '\n%s stage %s stream %s model %s ' % (current, stagelist, stream,  absskymodel)

    t = stagelist.split(':',1)
    if len(t)==2:
        stage, nextstage = t 
    else: stage,nextstage = t[0], None

    if stage=='update' or stage=='update_full':
        make_zip('pickle') # always update; latter diagnostic will use it, not the files
        diagnostic_plots.main('counts');
        logto = open(os.path.join(absskymodel,'converge.txt'), 'a')
        pipe.check_converge(absskymodel, log=logto)
        r = pipe.roirec(absskymodel)
        q = pipe.check_converge(absskymodel, add_neighbors=False)
        if max(r.niter)<8 and len(q)>1:
            create_stream('update')
        else:
            create_stream('sedinfo')
    elif stage=='sedinfo':
        make_zip('sedinfo')
        make_zip('sedfig','png')
        diagnostic_plots.main('fb')

    elif stage=='create' or stage=='create_reloc':
        ff = glob.glob(os.path.join(absskymodel, 'pickle', '*.pickle'))
        print 'found %d pickled ROI files' % len(ff)
        
        if nextstage is None:
            create_stream('update_full') # always start an update

    elif stage=='diffuse':
        make_zip('galfit_plots', 'png')
        make_zip('galfits_all')
        diagnostic_plots.main('gal')

    elif stage=='isodiffuse':
        make_zip('isofit_plots', 'png')
        make_zip('isofits')
        diagnostic_plots.main('iso')

    elif stage=='limb':
        #make_zip('limb_plots', 'png')
        make_zip('limb')

    elif stage=='finish':
        make_zip('pickle')
        make_zip('tsmap', 'png')
        make_zip('sedfig', 'png')
        diagnostic_plots.main('sources')

    elif stage=='tables':
        make_zip('ts_table')
        make_zip('kde_table')
        make_zip('counts_table')
    else:
        print 'stage %s not recognized for summary'%stage 
    if not args.test:
        if nextstage:
            create_stream(nextstage)
        tee.close()

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

