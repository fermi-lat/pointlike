## Interface to uw/like2/pipeline/uwpipeline.py
#
# Expect that POINTLIKE_DIR and SKYMODEL_SUBDIR are set, pointing to the model
# 

echo source $POINTLIKE_DIR/configure.sh
source $POINTLIKE_DIR/configure.sh

# this has to point to a writeable directory
export MPLCONFIGDIR=$POINTLIKE_DIR/.matplotlib

echo python -m uw/like2/pipeline/uwpipeline $*

python -m uw/like2/pipeline/uwpipeline $*
