"""
Package containing source analysis code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/__init__.py,v 1.17 2015/08/16 01:11:36 burnett Exp $

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
 ]

