"""
setup and run pointlike all-sky analysis for subset of ROIs

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipeline_job.py,v 1.21 2016/10/28 21:11:15 burnett Exp $
"""
import os, sys, logging, time, random
from collections import OrderedDict

import numpy as np

### Set variables corresponding to the environment
# these can be overridden by setting environment variables
# note that all are treated as strings
# they are converted to global variables in main

defaults=OrderedDict([
    ('HOSTNAME', '?'),
    ('PIPELINE_STREAMPATH', '-1'),
    ('PIPELINE_STREAM', '-1'),
    ('POINTLIKE_DIR','.'),
    ('SKYMODEL_SUBDIR','.'),
    ('stage', '?'),
    ('roi_list', '0-2'), # comma-delimited, like 1,5,6, or two numbers separated by a dash 9-30
    ])

def main( factory=None, **args):

    global roi_list, stage, HOSTNAME, PIPELINE_STREAMPATH, POINTLIKE_DIR,SKYMODEL_SUBDIR
    print ('\npointlike configuration')
    defaults['roi_list']=args.pop('roi_list', defaults['roi_list'])
    for key,default_value in defaults.items():
        # set local variables from default names, values from environment if set
        value = os.environ.get(key, default_value).strip()
        globals()[key]=value
        print ('\t%-20s: %s' % (key, value))
    np.seterr(invalid='warn', divide='warn')

    #### Check to see it abort flag is set
    ####
    try:
        stream = PIPELINE_STREAMPATH.split('.')[0]
        if os.path.exists(stream):
            print ('**** Found abort file for stream {}: stoping execution!'.format(stream))
            sys.stdout.flush()
            return -1
    except Exception as msg:
        print ('Failed to get stream: {}: continuing'.format(msg))
    
    # parse the list of ROIs: either a range or comma-delimited list of integers in range 0-1727
    t = roi_list.split('-')
    if len(t)>1:
        roi_list = range(int(t[0]), int(t[1])+1)
    else:
        u = filter(lambda x:len(x)>0, roi_list.split(','))
        assert len(u)>0, 'No entries found in list of rois: %s' % roi_list
        roi_list = map(int, u)
    first_roi, last_roi = roi_list[0], roi_list[-1]
    assert set(roi_list).issubset(range(1728)), 'ROI indeces, "%s", not in range 0-1727' % roi_list
   

    skymodeldir =SKYMODEL_SUBDIR #.replace('/a/wain025/g.glast.u55/','/afs/slac/g/glast/groups/') 
        
    streamlogdir = os.path.join(POINTLIKE_DIR,skymodeldir,'streamlogs')
    streamlogfile=os.path.join(streamlogdir,'stream%s.%04d.log' % ( PIPELINE_STREAMPATH.split('.')[0], int(PIPELINE_STREAM)) )
    if not os.path.exists(streamlogdir): os.mkdir(streamlogdir)
    print ('Logging execution progress to %s' % streamlogfile)
    if os.path.exists(streamlogfile):
       print ('...appending to existing ')
       with open(streamlogfile) as f:
            lines= f.read().split('\n')
            if len(lines)<2:
                print ('file is empty.')
            else:
                lastline = lines[-1] if len(lines[-1])>0 else lines[-2]
                if 'Start roi' in lastline: 
                   first_roi = int(lastline.split()[-1])
                   print ('Resuming execution with roi at %s' %first_roi)
                   assert first_roi in roi_list, 'Logic error roi %d not found in %d' %(first_roi, roi_list)
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO, filename=streamlogfile   )

    ### set up object with skymodel and data
    tzero = tnow = logging.time.time()
    logging.info('Start setup: rois %s ... %s on host %s' % (int(first_roi), int(last_roi), 
        os.environ.get('HOSTNAME','unknown')))
    
    ### Random delay to avoid possible collision with other jobs starting at the same time in the same machine
    time.sleep(0.5*random.random())
    
    if factory is not None:
        print ('--> executing process_roi in %s' % factory)
        g = factory(**args)
    else: 
        print ('-->Test execution')
        def g(n):
            print ('Would process roi %d' % n)
    mstage = stage
    stage=stage.split(':')[0] # allows multiple stages, separated by colons
        
    
    tprev, tnow= tnow, logging.time.time()
    logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))

    ### process eash ROI

    for s in roi_list[roi_list.index(first_roi):] :
        logging.info('Start roi %d' % s )
        g(s)
        tprev, tnow= tnow, logging.time.time()
        logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))
        print ('Elapsed time for ROI {}: {:.1f}'.format(s, tnow-tprev))
        sys.stdout.flush()
    print ('\n####################\nTotal elapsed time: {:.1f}'.format(tnow -tzero); sys.stdout.flush())
