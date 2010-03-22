"""
THe UW system for performing spectral analysis 

at this top level define a simple "factory" function giving a simple user access to the system

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/__init__.py,v 1.3 2010/03/18 20:36:18 burnett Exp $
Authors:  T. Burnett, M. Kerr, E. Wallace, M. Roth
"""

def factory(**kwargs):
    """
    return a like.thb_roi.myroi.ROIfactory object, iself a uw.like.pointspec.SpectralAnalysis subclass.
    It implements a __call__ function that returns a uw.thb_roi.myroi.MyROI, 
    a subclass of ROIAnalysis roi_analysis.ROIAnalysis

    keyword arguments
        ==========   =============
        keyword      description
        ==========   =============
        dataset      None      desciption of data set to use; if None, get "all data"
        emin         [100]     Default emin for fits
        free_radius  [1.0]     fit all sources within this radius
        prune_radius [0.1]     Do not include catalog sources within this radius of specified direction
        fit_bg_first [False]   When doing a spectral fit, first optimize the background
        use_gradient [True]    When doing a spectral fit, set the "use_gradient" option
        ==========   =============
    see uw.like.pointspec.SpectralAnalysis docstring for more keywords
    """
    from thb_roi import roi_factory

    defaults = {
        'dataset'     : None,
        'emin'        : 100,
        'free_radius' : 1.0,
        'prune_radius': 0.1,
        'fit_bg_first': False,
        'irf'         : 'P6_v8_diff',
        'use_gradient': True,
    }
    defaults.update(kwargs)
    return roi_factory.ROIfactory( **defaults )
    
