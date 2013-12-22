"""
setup and run pointlike all-sky analysis for subset of ROIs

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipeline_job.py,v 1.15 2013/12/16 16:14:41 burnett Exp $
"""
import os, sys, logging
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
    print '\npointlike configuration'
    defaults['roi_list']=args.pop('roi_list', defaults['roi_list'])
    for key,default_value in defaults.items():
        # set local variables from default names, values from environment if set
        value = os.environ.get(key, default_value).strip()
        globals()[key]=value
        print '\t%-20s: %s' % (key, value)
    np.seterr(invalid='warn', divide='warn')
    
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
    

    skymodeldir =SKYMODEL_SUBDIR.replace('/a/wain025/g.glast.u55/','/afs/slac/g/glast/groups/')  
    streamlogdir = os.path.join(POINTLIKE_DIR,skymodeldir,'streamlogs')
    streamlogfile=os.path.join(streamlogdir,'stream%s.%04d.log' % ( PIPELINE_STREAMPATH.split('.')[0], int(PIPELINE_STREAM)) )
    if not os.path.exists(streamlogdir): os.mkdir(streamlogdir)
    print 'Logging execution progress to %s' % streamlogfile
    if os.path.exists(streamlogfile):
       print '...appending to existing '
       with open(streamlogfile) as f:
            lines= f.read().split('\n')
            lastline = lines[-1] if len(lines[-1])>0 else lines[-2]
            if 'Start roi' in lastline: 
               first_roi = int(lastline.split()[-1])
               print 'Resuming execution with roi at %s' %first_roi
               assert first_roi in roi_list, 'Logic error roi %d not found in %d' %(first_roi, roi_list)
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO, filename=streamlogfile   )

    ### set up object with skymodel and data
    tzero = tnow = logging.time.time()
    logging.info('Start setup: rois %s ... %s on host %s' % (int(first_roi), int(last_roi), 
        os.environ.get('HOSTNAME','unknown')))
    if factory is not None:
        print '--> executing process_roi in %s' % factory
        g = factory(**args)
    else: 
        print '-->Test execution'
        def g(n):
            print 'Would process roi %d' % n
    mstage = stage
    stage=stage.split(':')[0] # allows multiple stages, separated by colons
        
    
    tprev, tnow= tnow, logging.time.time()
    logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))

    ### process eash ROI

    for s in roi_list[roi_list.index(first_roi):] :
        logging.info('Start roi %d' % s )
        try:
            g(s)
        except Exception, msg:
            print '***Exception raised for roi %d, "%s'  %(s,msg)
            with open('failed_rois.txt', 'wa') as bad:
                bad.write('s')
        tprev, tnow= tnow, logging.time.time()
        logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))
        print 'Elapsed time: %.1f' % (tnow-tprev)
        sys.stdout.flush()

#if __name__=='__main__':
#    main()