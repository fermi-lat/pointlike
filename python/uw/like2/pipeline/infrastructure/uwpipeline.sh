## Interface to uw/like2/pipeline/uwpipeline.py
#
# This script is submitted to the batch system by the pipeline 
# Expect that POINTLIKE_DIR and SKYMODEL_SUBDIR are set, pointing to the model
# 
ulimit -c 0

echo source $POINTLIKE_DIR/configure.sh
source $POINTLIKE_DIR/configure.sh

# this has to point to a writeable directory
export MPLCONFIGDIR=$POINTLIKE_DIR/.matplotlib

echo python -m uw/like2/pipeline/uwpipeline $*

python -m uw/like2/pipeline/uwpipeline $*
