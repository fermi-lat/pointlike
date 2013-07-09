"""
Package containing source analysis code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/__init__.py,v 1.3 2013/06/21 20:51:41 burnett Exp $

"""
# this is a list of modules in this package that have a class implementing an 'all_plots()' method.
# see app.py, which has an "import *" to include all of them
__all__=[
 'associations',
# 'components',
 'background', 
 'countplots',
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

