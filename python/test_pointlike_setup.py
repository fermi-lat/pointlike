#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.8 2008/01/25 22:36:11 burnett Exp $
from  pointlike_defaults import *

# the test data set does not have a diffuse background
Diffuse.file=""

# override these for stability, to compare with the output from previous versions
PointSourceLikelihood.minlevel=6
PointSourceLikelihood.gamma_list=[0,0,0,0,0,2.25,
        2.27,2.22,2.25,2.25,2.29,2.14,2.02,1.87]
PointSourceLikelihood.sigma_list= [0,0,0,0,0,0.343,
        0.335,0.319,  0.431, 0.449, 0.499, 0.566, 0.698, 0.818]

