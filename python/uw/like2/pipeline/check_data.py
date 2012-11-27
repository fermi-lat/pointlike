"""
Check that the data specification for this stream is valid, perhaps creating the intermediate files
$Header$
"""
import os, sys, glob, zipfile, logging, datetime
import numpy as np

from uw.like2.pipeline import pipe, processor
from uw.like2 import dataset

pointlike_dir = os.environ.get('POINTLIKE_DIR', '.')
skymodel = os.environ.get('SKYMODEL_SUBDIR', sys.argv[1] if len(sys.argv)>1 else '' )
stream = os.environ.get('PIPELINE_STREAM', '0')
absskymodel = os.path.join(pointlike_dir, skymodel)

tee = processor.OutputTee(os.path.join(absskymodel, 'summary_log.txt'))

stagelist = os.environ.get('stage', 'update' if len(sys.argv)<3 else sys.argv[2])
current = str(datetime.datetime.today())[:16]
print '\n%s stage %s stream %s model %s ' % (current, stagelist, stream,  absskymodel)

rc = dataset.validate(absskymodel, nocreate=True)
tee.close()

if not rc: raise Exception('Failed to validate data')