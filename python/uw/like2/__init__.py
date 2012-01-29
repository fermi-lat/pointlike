"""
Package containing likelihood analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/__init__.py,v 1.3 2012/01/27 15:10:59 burnett Exp $
Authors:  T. Burnett, M. Kerr, E. Wallace, M. Roth, J. Lande

Modules in this package, and basic dependencies

main 
    roisetup 
        diffuse 
            models
        dataset
        skymodel
            sources
    roistat
        sourcelist
        bandlike
    localization 
        plotting.tsmaps
    sedfuns
    printing
    
contained packages
    pipeline
    plotting
    pub
"""
# for interactive convenience
from .main import factory
from skymaps import SkyDir
