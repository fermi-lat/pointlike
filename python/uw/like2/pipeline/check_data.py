"""
Check that the data specification for this stream is valid, perhaps creating the intermediate files
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/check_data.py,v 1.6 2013/03/05 19:49:25 burnett Exp $
"""
import os, sys, glob, zipfile, logging, datetime
import numpy as np

from uw.like2 import (dataset, tools)

def main(args=None):

    if args is not None:
        stage = args.stage[0]
        nocreate = True
    else:
        # if called directly, may create
        stage = os.environ.get('stage', 'create' )
        nocreate = False
    if stage!='create':
        print 'assume validated'
        return
        

    pointlike_dir = os.environ.get('POINTLIKE_DIR', '.')
    skymodel = os.environ.get('SKYMODEL_SUBDIR', '.')
    stream = os.environ.get('PIPELINE_STREAM', '-1')
    absskymodel = os.path.join(pointlike_dir, skymodel)

    if args is not None:
        tee = tools.OutputTee(os.path.join(absskymodel, 'summary_log.txt'))

    current = str(datetime.datetime.today())[:16]
    print '\n%s stage %s stream %s model %s ' % (current, stage, stream,  absskymodel)

    rc = dataset.validate(absskymodel, nocreate=nocreate)
    print 'Validated' if rc else 'NOT validated'
    if args is not None:
        tee.close()

    if not rc: raise Exception('Failed to validate data')
 
if __name__=='__main__':
    main()