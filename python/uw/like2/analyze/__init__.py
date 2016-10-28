"""
Package containing source analysis code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/__init__.py,v 1.19 2016/04/27 02:22:06 burnett Exp $

"""
# this is a list of modules in this package that have a class implementing an 'all_plots()' method.
# see app.py, which has an "import *" to include all of them
# (comment out those not relevant)
__all__=[
 'associations',
# 'components',
 'background', 
 'counts',
 'pulsars',
 #'collection',
 'config',
 'cputime',
 'data',
 #'diffuse_info',
 'diffuse_correction',
 'environment',
 'export',
 'find_peak',
 #'fluxcorr',
# 'fluxcorriso',
 #'frontback',
 #'galactic',
 #'galacticspectra',
 #'gtlikecomparison',
 'hptables',
 #'isotropic',
 #'isotropicspectra',
 #'limb',
 #'limbrefit',
 'localization',
 #'localization1k',
# 'pgwseedcheck',
 'maps',
 'pgwave',
 #'ptstable',
 #'pulsarseedcheck',
 #'roi_info',
 'residuals',
 'seedcheck',
 'sourcecomparison',
 'sourceinfo',
 #'sourcetotal',
 #'sunmoon',
 #'sunmoonrefit',
 'uwsourcecomparison',
 'transientinfo',
 'transient_catalog',
 'residual_maps',
 'gc_comparison',
 ]

