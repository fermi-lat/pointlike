"""
$Header$
author: T.Burnett <tburnett@u.washington.edu>
"""

import copy
from uw.like import pointspec, Models, pointspec_helpers
from skymaps import SkyDir
# in this package
import data, roi_setup, catalog, myroi

class ROIfactory(pointspec.SpectralAnalysis):
    """ subclass of SpectralAnalysis, coded as a ROI "factory"
    """

    def __init__(self,  **kwargs):
        
        defaults = {
            'fit_bg_first': False,
            'use_gradient': True,
            }
        self.log = None
        defaults.update(kwargs)
        analysis_environment =  data.MyAnalysisEnvironment( **defaults)
        super(ROIfactory,self).__init__( analysis_environment, **defaults)

        self.cb = roi_setup.ConsistentBackground(self.ae, self.background, quiet=self.quiet)
        aux = self.__dict__.pop('aux_cat', None)
        if aux:
            print 'adding sources from catalog'
            self.cb.append(aux)
            
        if not self.quiet: print >>self.log, self

    def __str__(self):
        s = 'ROIfactory configuration:\n'
        ignore = ('psf', 'exposure', 'cb', 'mc_energy', 'mc_src_id')
        for key in sorted(self.__dict__.keys()):
            if key in self.ae.__dict__ or key in ignore: continue # avoid duplication internal functions
            s += '\t%-20s: %s\n' %(key, self.__dict__[key])
        return s


    def roi(self, point_sources = None, bgmodels = None, previous_fit = None, **kwargs):

        """
        return an ROIAnalysis object with default settings.

        point_sources    [None] a list of PointSource objects to merge with a Catalog list
                         (if None, the nearest catalog source will be fit)

        bgmodels         a list of ROIBackgroundModels with which to override the default
                         isotropic and Galactic backgrounds (optional)

        previous_fit     [None] a file containing the results of an earlier spectral fit;
                         if not None, set spectral values to this fit
                         ***WARNING*** not tested!

        Optional Keyword Arguments:
            ==========   =============
            keyword      description
            ==========   =============
            nocat        [False] if True, do not add additional sources from a catalog 
            bg_smodels   [None]  a list of spectral models to replace the default ones in ConsistentBackground
                                 i.e., a custom set of spectral scaling models
            glat         [None]  the Galactic latitude of the source; sets default free parameters in diffuse
            fit_emin     [100,100] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            ==========   =============
        """
        ps_manager, bg_manager = super(ROIfactory, self).roi(point_sources, bgmodels, previous_fit, no_roi=True, **kwargs)
        
        return myroi.MyROI(ps_manager,bg_manager, self, quiet=self.quiet, **kwargs)


    def __call__(self, sources=None, max_roi=None, min_roi=None, roi_dir=None, **kwargs):
        """ 
        return a MyROI object

        sources -- a list (name, skydir, model) of sources, or a text string for one source
                   with:
                      name ra dec 
                      name     (look up in catalog if so)
                    Assume PowerLaw
        roi_dir -- center of ROI: if None, use first source
        min_roi, max_roi -- if specified, change default settings
    
        Optional Keyword Arguments:
            ==========   =============
            keyword      description
            ==========   =============
            nocat        [False] if True, do not add additional sources from a catalog 
            bg_smodels   [None]  a list of spectral models to replace the default ones in ConsistentBackground
                                 i.e., a custom set of spectral scaling models
            glat         [None]  the Galactic latitude of the source; sets default free parameters in diffuse
            fit_emin     [100,100] minimum energies (separate for front and back) to use in spectral fitting.
            fit_emax     [1e5,1e5] maximum energies (separate for front and back) to use in spectral fitting.
            free_radius  [see factory] fit all sources within this radius
            prune_radius [see factory] Do not include catalog sources within this radius of specified direction
            use_gradient [True]    When doing a spectral fit, set the "use_gradient" option
            model        'powerlaw' One of: 'powerlaw', 'expcuoff', 'logparabola'
            ==========   =============
        """

        ps = None
        if sources is not None:
            if type(sources)==type(' '):
                try:
                    name, ra, dec = catalog.find_source(sources)
                except:
                    print 'expected string with name ra dec'
                    raise
                modelname = kwargs.pop('model', 'powerlaw')
                models ={'powerlaw':    Models.PowerLaw(), 
                         'expcutoff':   Models.ExpCutoff(), 
                         'logparabola': Models.LogParabola(),
                        }
                if modelname not in models:
                    raise Exception('model "%s" not recognized: %s expect one of: '% (modelname, models.keys()))
                model_object = copy.copy(models[modelname])
                ps = [pointspec_helpers.PointSource(SkyDir(float(ra), float(dec)), name, model_object)]
            else:
                ps = []
                for s in sources:
                    name,dir = s[:2]
                    model = s[2] if len(s)==3 else Models.PowerLaw();
                    ps.append(pointspec_helpers.PointSource(dir,name,model) )
        self.roi_dir = roi_dir or ps[0].skydir
        if min_roi: self.minROI=min_roi  
        if max_roi: self.maxROI=max_roi 
        # pass on default values for optional args
        passon_list = ('free_radius', 'prune_radius', 'fit_bg_first', 'use_gradient')
        for x in passon_list:
                if x not in kwargs: kwargs[x] = self.__dict__[x]

        emin,emax = self.emin, self.emax
        r = self.roi(point_sources = ps, fit_emin=[emin,emin],fit_emax=[emax,emax], **kwargs)

        # if a different direction, we need to disable the original, most likely the nearest
        # a newer
        if r.psm.point_sources[1].name.strip() == r.name:
            r.psm.models[1].p[0]=-20

        return r
        
    def source_list(self):
        """ return a recarray of all the sources in the current model
        """
        return self.cb.cm.source_recarray()
 