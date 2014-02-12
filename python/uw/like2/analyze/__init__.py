"""
Package containing source analysis code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/__init__.py,v 1.12 2014/01/16 19:36:43 burnett Exp $

"""
# this is a list of modules in this package that have a class implementing an 'all_plots()' method.
# see app.py, which has an "import *" to include all of them
__all__=[
 'associations',
# 'components',
 'background', 
 'counts',
 'collection',
 'config',
 'cputime',
 'data',
 'diffuse_info',
 'environment',
 'export',
 'find_peak',
 'fluxcorr',
# 'fluxcorriso',
 'frontback',
 'galactic',
 'galacticspectra',
 'gtlikecomparison',
 'hptables',
 'isotropic',
 'isotropicspectra',
 'limb',
 'limbrefit',
 'localization',
 'localization1k',
# 'pgwseedcheck',
 'maps',
 'ptstable',
 'pulsarseedcheck',
 #'roi_info',
 'residuals',
 'seedcheck',
 'sourcecomparison',
 'sourceinfo',
 'sourcetotal',
 'sunmoon',
 'sunmoonrefit',
 'uwsourcecomparison',
 ]

