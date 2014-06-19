"""
Run after a successful UWpipeline/job_task

Summarize the execution, 
mostly zipping the files generated by the multiple jobs and submitting follow-up streams
diagnostic plots are now done by the summary task


$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/check_converge.py,v 1.30 2014/06/19 13:39:16 burnett Exp $
"""
import os, sys, glob, zipfile, logging, datetime, argparse, subprocess
import numpy as np
import pandas as pd

from uw.like2 import (tools, )
from uw.like2.pipeline import pipe, stream 
from uw.like2.pub import healpix_map


def streamInfo( stream, path='.'):
    streamlogs=sorted(glob.glob(os.path.join(path,'streamlogs','stream%s.*.log'%stream)))
    if len(streamlogs)==0:
        print 'Did not find stream %d logfiles'% stream
        return
    etimes = []
    substream = []
    nex = []
    for f in streamlogs:
        text=open(f).read()
        etimes.append(float(text.split('\n')[-2].split()[-1][:-1]))
        nex.append( sum(np.array(['setup' in line for line in text.split('\n') ] )))
        substream.append( int(f.split('.')[-2]) )
    print '\t found %d streamlog files for stream %s ...' %(len(streamlogs), stream), 
    times = pd.DataFrame(dict(etimes=etimes,nex=nex), index=pd.Index(substream, name='substream'))#,name='stream%d'%stream)
    etimes = times.etimes
    print '\texecution times: mean=%.1f, max=%.1f'% (etimes.mean(), etimes.max())
    print '\ttotal execution time: %.1f hours' % (etimes.sum()/3600.)
    print '\tnumber of executions: %s' % np.histogram(times.nex, range(7))[0][1:]
    bad = times[times.etimes>2*times.etimes.mean()]
    if len(bad)>0:
        print '\tsubstreams with times>2*mean:'
        print bad
    return times

def main(args):
    """ 
    """
    def fixpath(s):
        if os.path.exists(s): return s
        r= s.replace('/a/wain025/g.glast.u55','/afs/slac/g/glast/groups')
        assert os.path.exists(r), 'paths do not exist: %s or %s' % (s,r)
        return r
    pointlike_dir=fixpath(args.pointlike_dir) # = os.environ.get('POINTLIKE_DIR', '.')
    skymodel     =fixpath(args.skymodel) # = os.environ.get('SKYMODEL_SUBDIR', sys.argv[1] )
    
    stream = args.stream 
    stagelist = args.stage
    if hasattr(stagelist, '__iter__'): stagelist=stagelist[0] #handle local or from uwpipeline.py
   
    absskymodel = os.path.join(pointlike_dir, skymodel)

    def make_zip(fname,  ext='pickle', select=None):
        if select is not None:
            ff = glob.glob(os.path.join(absskymodel, select))
        else:
            ff = glob.glob(os.path.join(absskymodel, fname, '*.'+ext))
        if len(ff)==0:
            print 'no files found to zip in folder %s' %fname
            return
        if len(ff)!=1728 and ext=='pickle':
            raise Exception('Found %d pickle files, expected 1728'%len(ff))
        print 'found %d *.%s in folder %s ...' % ( len(ff),ext, fname,) ,
        with zipfile.ZipFile(os.path.join(absskymodel, fname+'.zip'), 'w') as z:
            for filename in ff:
                z.write( filename, os.path.join(fname,os.path.split(filename)[-1]))
        print ' zipped into file %s.zip' %fname
        
    def create_stream(newstage, **kw):
        print 'Starting stage %s' % newstage
        args = dict(stage=newstage, POINTLIKE_DIR=pointlike_dir, SCRIPT=pointlike_dir, SKYMODEL_SUBDIR=skymodel)
        args.update(**kw)
        ps = stream.PipelineStream()
        
        #cmd = ['/afs/slac/g/glast/ground/bin/pipeline','createStream']
        #cmd += [' -D '] +  ['"' + ', '.join(['%s=%s' %item for item in args.items()]) + '"']
        vars = ', '.join(['%s=%s' %item for item in args.items()]) 
        #cmd += ['UWpipeline']
        ps(newstage, vars)
        #print '\n-->' + ' '.join(cmd)
        
        try:
            os.system(' '.join(cmd))
            #print subprocess.check_output(cmd)
        except Exception, msg:
            print '***Failed: "%s"' % msg
            raise


    if not args.test:
        tee = tools.OutputTee(os.path.join(absskymodel, 'summary_log.txt'))

    streamInfo(stream, absskymodel)

    os.chdir(absskymodel) # useful for diagnostics below
    current = str(datetime.datetime.today())[:16]
    print '\n%s stage %s stream %s model %s ' % (current, stagelist, stream,  absskymodel)
    if os.path.exists('failed_rois.txt'):
       
       failed = sorted(map(int, open('failed_rois.txt').read().split()))
       print 'failed rois %s' % failed
       raise Exception('failed rois %s' % failed)

    t = stagelist.split(':',1)
    if len(t)==2:
        stage, nextstage = t 
    else: stage,nextstage = t[0], None

       
    if stage.split('_')[0]=='update':
        make_zip('pickle') # always update; latter diagnostic will use it, not the files
        logto = open(os.path.join(absskymodel,'converge.txt'), 'a')
        qq=pipe.check_converge(absskymodel, log=logto)
        r = pipe.roirec(absskymodel)
        q = pipe.check_converge(absskymodel, add_neighbors=False)
        #open(os.path.join(pointlike_dir,'update_roi_list.txt'), 'w').write('\n'.join(map(str, sorted(qq))))
        open('update_roi_list.txt', 'w').write('\n'.join(map(str, sorted(qq))))
        if stage!='update_only':
            if  len(q)>1:
                if len(qq)> 200:
                    create_stream('update')
                else:
                    create_stream('update', job_list='$SKYMODEL_SUBDIR/update_roi_list.txt')
            else:
                create_stream('finish')
            
    elif stage=='sedinfo':
        make_zip('sedinfo')
        make_zip('sedfig','png')

    elif stage=='create' or stage=='create_reloc':
        # always do one more stream after initial
        make_zip('pickle')
        if nextstage is None:
            create_stream('update_full') # always start an update

    elif stage=='diffuse':
        make_zip('galfit_plots', 'png')
        make_zip('galfits_all')

    elif stage=='isodiffuse':
        make_zip('isofit_plots', 'png')
        make_zip('isofits')

    elif stage=='limb':
        make_zip('limb')

    elif stage=='finish' or stage=='counts':
        make_zip('pickle')

    elif stage=='tables':
        names = 'ts kde counts'.split() 
        healpix_map.assemble_tables(names)

    elif stage=='pulsar_detection':
        healpix_map.assemble_tables(['pts'])
        
    elif stage=='seedcheck':
        make_zip('seedcheck', select='seedcheck/SEED*') #needs implementation
        
    elif stage=='pseedcheck':
        make_zip('pseedcheck', select='seedcheck/PSEED*')
        
    elif stage=='seedcheck_PGW':
        make_zip('seedcheck_PGW', select='seedcheck/PGW*')
        
    else: # catch fluxcorr, any others like
        if os.path.exists(stage):
            make_zip(stage)
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
    parser.add_argument('--test', action='store_false', help='test mode') ######################
    args = parser.parse_args()
    main(args)

