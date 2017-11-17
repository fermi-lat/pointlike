"""
Check that the data specification for this stream is valid, perhaps creating the intermediate files
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/check_data.py,v 1.10 2015/12/03 17:33:07 burnett Exp $
"""
import os, sys, glob, zipfile, logging, datetime
import numpy as np

from uw.like2 import (dataset, tools)

def main(args=None):
    if os.path.exists('kill'):
        raise Exception('Kill this job')
    if args is not None:
        stage = args.stage[0]
        nocreate = False # Was True, but only important to run the creation for large data sets
    else:
        # if called directly, may create
        stage = os.environ.get('stage', 'create' )
        nocreate = False
    stage=stage.split('_')[0]
    if stage!='create' and stage!='monthly':
        print 'Not creating a model: assume data validated'
        return
        

    pointlike_dir = os.environ.get('POINTLIKE_DIR', '.')
    skymodel = os.environ.get('SKYMODEL_SUBDIR', '.')
    stream = os.environ.get('PIPELINE_STREAM', '-1')
    absskymodel = os.path.join(pointlike_dir, skymodel) if args is not None else '.'

    if args is not None:
        tee = tools.OutputTee(os.path.join(absskymodel, 'summary_log.txt'))

    current = str(datetime.datetime.today())[:16]
    print '\n%s stage %s stream %s model %s ' % (current, stage, stream,  absskymodel)

    if 'CUSTOM_IRF_DIR' not in os.environ and os.path.exists(os.path.expandvars('$FERMI/custom_irfs')):
        os.environ['CUSTOM_IRF_DIR'] = os.path.expandvars('$FERMI/custom_irfs')

    rc = dataset.validate(absskymodel, nocreate=nocreate, quiet=False)
    print 'Data is validated' if rc else 'NOT validated'
    if args is not None:
        tee.close()

    if not rc: raise Exception('Failed to validate data')
 
if __name__=='__main__':
    main()
