"""
setup and run pointlike all-sky analysis for subset of ROIs

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipeline_job.py,v 1.13 2013/05/14 15:35:15 burnett Exp $
"""
import os, sys, logging
from collections import OrderedDict

import numpy as np
#from uw.like2 import main


### Set variables corresponding to the environment
# these can be overridden by setting environment variables
# note that all are treated as strings
# if None, must be set

defaults=OrderedDict([
    ('HOSTNAME', '?'),
    ('PIPELINE_STREAMPATH', '-1'),
    ('PIPELINE_STREAM', '-1'),
    ('POINTLIKE_DIR','.'),
    ('SKYMODEL_SUBDIR','.'),
    ('stage', '?'),
    ('begin_roi','2'), # defaults for test: #2 has one source
    ('end_roi', '4'),
    ('roi_list', ''), # comma-delimited list, to replace begin_roi, end_roi if present

    ])

def main( update=None):

    print '\npointlike skymodel configuration'
    for key,default_value in defaults.items():
        item = (key, os.environ.get(key, default_value).strip())
        exec('%s="%s"'%item)
        print '\t%-20s: %s' % item
    if update is not None:
        print '--> executing process_roi in %s' % update.__class__.__name__
        g = update
    else: 
        print 'Test execution'
        def g(n):
            print 'Would process roi %d' % n
    np.seterr(invalid='warn', divide='warn')
    
    if roi_list =='':
        roi_list = range(int(begin_roi), int(end_roi))
    else:
        if roi_list.endswith(','):
            roi_list=roi_list[:-1]
        roi_list = map(int, roi_list.split(','))
    first_roi, last_roi = roi_list[0], roi_list[-1]
    

    skymodeldir =SKYMODEL_SUBDIR.replace('/a/wain025/g.glast.u55/','/afs/slac/g/glast/groups/')  
    streamlogdir = os.path.join(POINTLIKE_DIR,skymodeldir,'streamlogs')
    streamlogfile=os.path.join(streamlogdir,'stream%s.%04d.log' % ( PIPELINE_STREAMPATH.split('.')[0], int(PIPELINE_STREAM)) )
    if not os.path.exists(streamlogdir): os.mkdir(streamlogdir)
    print 'Logging execution progress to %s' % streamlogfile
    if os.path.exists(streamlogfile):
       print 'found existing stream log file:%s' % streamlogfile
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
        tprev, tnow= tnow, logging.time.time()
        logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))
        sys.stdout.flush()

#if __name__=='__main__':
#    main()