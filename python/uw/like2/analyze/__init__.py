"""
Package containing source analysis code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/__init__.py,v 1.6 2013/08/03 19:39:48 burnett Exp $

"""
# this is a list of modules in this package that have a class implementing an 'all_plots()' method.
# see app.py, which has an "import *" to include all of them
__all__=[
 'associations',
# 'components',
 'background', 
 'countplots',
 'collection',
 'config',
 'data',
 'environment',
 'export',
 'find_peak',
 'fluxcorr',
# 'fluxcorriso',
 'frontbacksedplots',
 'galactic',
 'galacticspectra',
 'gtlikecomparison',
 'hptables',
 'isotropic',
# 'isotropicspectra',
 'limb',
 'limbrefit',
 'localization',
 'localization1k',
# 'pgwseedcheck',
 'maps',
 'ptstable',
 'pulsarseedcheck',
 'roi_info',
 'seedcheck',
 'sourcecomparison',
 'sourceinfo',
 'sourcetotal',
 'sunmoon',
 'sunmoonrefit',
 'uwsourcecomparison',
 ]

