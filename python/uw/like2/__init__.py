"""
Package containing likelihood analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/__init__,v 1.2 2011/09/28 17:35:52 burnett Exp $
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

from .main import factory
