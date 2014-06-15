#
# Setup to run pointlike in the catalog setup
#
release=09-33-03
#release=09-31-01
build=redhat6-x86_64-64bit-gcc44

catalog_dir=/afs/slac/g/glast/groups/catalog/
pointlike=$catalog_dir/pointlike
#mypython=/afs/slac/g/glast/users/burnett/python
mypython=$pointlike/python

# pointlike uses $FERMI to find default paths for data, diffuse, catalog files
export FERMI=$pointlike/fermi

# location of custom IRF files not found under CALDB structure
# (this is now automatic with like2)
# (oops, not yet in create stage)
# export CUSTOM_IRF_DIR=$pointlike/custom_irfs

# Define GLAST_EXT to externals linking
export GLAST_EXT=/afs/slac/g/glast/ground/GLAST_EXT/$build

# BASE_DIR should be set to root of base installation, INST_DIR
builds_dir=/nfs/farm/g/glast/u35/ReleaseManagerBuild
export BASE_DIR=$builds_dir/$build/Optimized/ScienceTools/$release

# If local packages overriding, set INST_DIR to $FERMI/inst_dir
#export INST_DIR=$FERMI/inst_dir
export INST_DIR=$BASE_DIR

#run the scons-generated setup script, then update paths and env vars 
source $INST_DIR/bin/$build-Optimized/_setup.sh
export CALDB=${BASE_DIR}/irfs/caldb/CALDB/data/glast/lat
export TIMING_DIR=${GLAST_EXT}/extFiles/v0r9/jplephem
export PATH=$catalog_dir/python/anaconda/bin:${PATH}
export PYTHONPATH=$mypython:${PYTHONPATH}

# IPython configuration in $pointlike/.ipython
# No, this does not support multiple users: each must have own
# There is a reason to set this env var to nfs, since notebook server can die if token expires
#export IPYTHONDIR=$pointlike/.ipython

# use a local matplotlib configuration. Important for batch to use agg at least
export MPLCONFIGDIR=$pointlike/.matplotlib

# convenient aliases,to be run in a skymodel folder 
alias uwpipeline='python -m uw/like2/pipeline/uwpipeline'
# for building superceding packages, in INST_DIR by default
alias scons='scons --with-GLAST-EXT=$GLAST_EXT'
alias supersede='(cd $BASE_DIR; scons --supersede=$INST_DIR --compile-opt)'
alias diagnostic_plots='python -m uw/like2/pipeline/diagnostic_plots'
alias plot_browser='firefox $FERMI/skymodels/plot_browser/index.html'
# short forms for starting stream, generating plots
alias pipe='python -m uw/like2/pipeline/uwpipeline'
alias plots='python -m uw/like2/analyze/app'
#start a notebook server
alias notebook='ipython notebook --notebook-dir=$mypython/uw/like2/notebooks'

# should prevent core dumps when run in batch
ulimit -c 0

if [ "$PS1" ]; then
  echo "MYPYTHON :" $mypython
  echo "FERMI    :" $FERMI
  echo "INST_DIR :" $INST_DIR
  echo "BASE_DIR :" $BASE_DIR
fi
