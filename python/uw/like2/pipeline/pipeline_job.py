"""
setup and run pointlike all-sky analysis for subset of ROIs

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/pipeline_job.py,v 1.9 2012/12/27 15:49:02 burnett Exp $
"""
import os, sys, logging
from collections import OrderedDict

import numpy as np
from uw.like2.pipeline import pipe


### Set variables corresponding to the environment
# these can be overridden by setting environment variables
# note that all are treated as strings
# if None, must be set

defaults=OrderedDict([
	('HOSTNAME', '?'),
	('PIPELINE_STREAMPATH', '-1'),
	('PIPELINE_STREAM', '-1'),
    
    ('POINTLIKE_DIR','.'),
    ('SKYMODEL_SUBDIR','?'),
    ('stage', '?'),
    ('begin_roi','2'), # defaults for test: #2 has one source
    ('end_roi', '3'),
    ])
 
def main( update):

    print '\npointlike skymodel configuration'
    for key,default_value in defaults.items():
        item = (key, os.environ.get(key, default_value))
        exec('%s="%s"'%item)
        print '\t%-20s: %s' % item
    print '--> executing:\n', update.setup()
    np.seterr(invalid='warn', divide='warn')

    streamlogdir = os.path.join(POINTLIKE_DIR,SKYMODEL_SUBDIR,'streamlogs')
    streamlogfile=os.path.join(streamlogdir,'stream%s.%04d.log' % ( PIPELINE_STREAMPATH.split('.')[0], int(PIPELINE_STREAM)) )
    if not os.path.exists(streamlogdir): os.mkdir(streamlogdir)
    print 'Logging execution progress to %s' % streamlogfile
    if os.path.exists(streamlogfile):
       print 'found existing steam log file:%s' % streamlogfile
       with open(streamlogfile) as f:
            lines= f.read().split('\n')
            lastline = lines[-1] if len(lines[-1])>0 else lines[-2]
            if 'Start roi' in lastline: 
               begin_roi = lastline.split()[-1]
               print 'Resuming execution with roi at %s' %begin_roi
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO, filename=streamlogfile   )

    ### set up object with skymodel and data
    tzero = tnow = logging.time.time()
    logging.info('Start setup: rois %s-%s on host %s' % (int(begin_roi), int(end_roi)-1, 
        os.environ.get('HOSTNAME','unknown')))

    mstage = stage
    stage=stage.split(':')[0] # allows multiple stages, separated by colons
        
    g = update.g()
    tprev, tnow= tnow, logging.time.time()
    logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))

    ### process eash ROI

    for s in range(int(begin_roi), int(end_roi)):
        logging.info('Start roi %d' % s )
        g(s)
        tprev, tnow= tnow, logging.time.time()
        logging.info('Finish: elapsed= %.1f (total %.1f)' % ( tnow-tprev, tnow-tzero ))
        sys.stdout.flush()

#if __name__=='__main__':
#    main()