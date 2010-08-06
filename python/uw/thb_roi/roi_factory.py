"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/roi_factory.py,v 1.6 2010/06/20 13:53:17 burnett Exp $
author: T.Burnett <tburnett@u.washington.edu>
"""

import os, copy, types
import numpy as np
from uw.like import pointspec, Models, pointspec_helpers, roi_managers, roi_diffuse
from skymaps import SkyDir
# in this package
import roi_setup, myroi, config, findsource

## this should go to Models 
def make_model( name='PowerLaw', pars=(1e-12/1.53, 2.3), quiet=True):
    """ make a fit model:  
        name: ['PowerLaw']
        pars: [('1e-12/1.53, 2.3)] If more default parameters, will fill in from defaults
    """
    try:
        description = Models.DefaultModelValues.simple_models[name] # look it up
        factory = eval('Models.'+name)  # the class 
        defaults = list(description['p'])
        if pars is not None:
            for i in range(len(pars)): defaults[i] = pars[i]
        if not quiet: print 'Set up model %s, inital parameters %s' % (name, np.asarray(defaults))
        
    except KeyError:
        raise Exception('model "%s" not recognized: %s expect one of: '\
                        % (name, Models.DefaultModelValues.keys()))
    except:
        print 'Failed to make model %s with parameters %s' % (name ,pars) 
        raise
    return factory(p=defaults)
    

class ROIfactory(pointspec.SpectralAnalysis):
    """ subclass of SpectralAnalysis, coded as a ROI "factory"
    """

    def __init__(self, **kwargs): 
        
        self.log = kwargs.pop('log', None)
        analysis_environment = config.AE(**kwargs)
        super(ROIfactory,self).__init__( analysis_environment)

        ae = analysis_environment
        catpars = dict(  # this is ugly, need general solution: want to pass in ae, extract only the ones catman wants.
                pulsar_dict = ae.__dict__.pop('pulsar_dict', None), 
                quiet = ae.quiet,
                verbose = ae.verbose,
                )
        self.catman= roi_setup.CatalogManager(os.path.join(analysis_environment.catdir, analysis_environment.catalog),
            **catpars)

        aux = self.__dict__.pop('aux_cat', None)
        if aux:
            if not self.quiet: print 'adding sources from catalog'
            self.cb.append(aux)
            
        if not self.quiet: 
            print >>self.log, self
            if self.log is not None: self.log.close()
        


    def __str__(self):
        s = 'ROIfactory analysis environment:\n'
        s += self.ae.__str__()
        #ignore = ('ae', 'psf', 'exposure', 'cb', 'catman')
        #for key in sorted(self.__dict__.keys()):
        #    if key in ignore: continue # avoid duplication internal functions and copy in ae
        #    s += '\t%-20s: %s\n' %(key, self.__dict__[key])
        return s

    def __call__(self, sources, max_roi=None, min_roi=None, roi_dir=None, **kwargs):
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
            model        'PowerLaw' One of: 'PowerLaw', 'ExpCuoff', 'LogParabola', ... (see uw.like.Models)
            model_par    [1e-12, 2.3]    a list of initial model parameters
            ==========   =============
        """
        ps = None
        model_specified = True
        modelname = kwargs.pop('model', None)
        model_par = kwargs.pop('model_par', (1e-12, 2.3))
        if modelname is None:
            model_specified = False
            modelname = 'PowerLaw'
         
        if type(sources)==types.StringType:
            #try:
            name, ra, dec = findsource.find_source(sources, catalog=self.catalog)
            if not self.quiet: print '%s --> %s at (%s, %s)' % (sources, name, ra, dec)
            #except Exception, arg:
            #    raise Exception, 'source name "%s" unrecognized, reason %s' % (sources, arg)
            themodel = make_model(modelname, model_par)
            ps = [pointspec_helpers.PointSource(SkyDir(float(ra), float(dec)), name, themodel)]
        else:
            ps = []
            for s in sources:
                name,dir = s[:2]
                themodel = s[2] if len(s)==3 else make_model(modelname, model_par);
                ps.append(pointspec_helpers.PointSource(dir,name, themodel) )
         
        skydir =self.roi_dir = roi_dir or ps[0].skydir
        if min_roi: self.minROI=min_roi  
        if max_roi: self.maxROI=max_roi 
        
        # pass on default values for optional args
        passon_list = ('free_radius', 'prune_radius', 'fit_bg_first', 'use_gradient')
        for x in passon_list:
                if x not in kwargs: kwargs[x] = self.__dict__[x]

        # add background point sources from the CatalogManager
        ps = ps + self.catman(skydir, self.maxROI)
        excluded = self.catman.exclude
        if excluded is not None:
            if not model_specified:
                ps[0].model = excluded.model # if replacing a cat source, use its model
            elif modelname=='ExpCutoff':
                # specified ExpCufoff: copy catalog values for flux, index
                ps[0].model.p[:2] = excluded.model.p[:2]
        
        ps_manager = roi_managers.ROIPointSourceManager(ps, skydir,quiet=self.quiet)
        
        # following is deprecated!
        #bg_manager = roi_managers.ROIBackgroundManager(self, self.bgmodels(skydir), self.roi_dir,quiet=self.quiet)
        # new way
        diffuse_mapper = lambda x: roi_diffuse.ROIDiffuseModel_OTF(self, x, skydir)
        diffuse_sources = pointspec_helpers.get_default_diffuse( *self.diffuse)

        diffuse_models = [diffuse_mapper(ds) for ds in diffuse_sources]
        bg_manager = roi_managers.ROIDiffuseManager(diffuse_models,skydir,quiet=self.quiet)

        emin,emax = self.fit_emin, self.fit_emax
        r = myroi.MyROI(skydir, ps_manager, bg_manager, self, 
                        point_sources = ps,
                        fit_emin=[emin,emin], fit_emax=[emax,emax],
                        quiet=self.quiet, **kwargs)

        # if a different direction, we need to disable the original, most likely the nearest
        if len(r.psm.point_sources)>1 and r.psm.point_sources[1].name.strip() == r.name:
            r.psm.models[1].p[0]=-20

        return r
    
        
    def source_list(self):
        """ return a recarray of all the sources in the current model
        """
        return self.catman.source_recarray()
 