"""
Check that the data specification for this stream is valid, perhaps creating the intermediate files
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/check_data.py,v 1.2 2012/12/23 13:32:10 burnett Exp $
"""
import os, sys, glob, zipfile, logging, datetime
import numpy as np

from uw.like2.pipeline import pipe, processor
from uw.like2 import dataset

def main(args=None):

    if args is not None:
        stage = args.stage[0]
    else:
        stage = os.environ.get('stage', 'update' )

    pointlike_dir = os.environ.get('POINTLIKE_DIR', '.')
    skymodel = os.environ.get('SKYMODEL_SUBDIR', '.')
    stream = os.environ.get('PIPELINE_STREAM', '-1')
    absskymodel = os.path.join(pointlike_dir, skymodel)

    tee = processor.OutputTee(os.path.join(absskymodel, 'summary_log.txt'))

    current = str(datetime.datetime.today())[:16]
    print '\n%s stage %s stream %s model %s ' % (current, stage, stream,  absskymodel)

    rc = dataset.validate(absskymodel, nocreate=True)
    print 'Validated' if rc else 'NOT validated'
    tee.close()

    if not rc: raise Exception('Failed to validate data')
 
if __name__=='__main__':
    main()