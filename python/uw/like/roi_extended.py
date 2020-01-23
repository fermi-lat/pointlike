""" Provides classes to encapsulate and manipulate extended sources. 

    This code all derives from objects in roi_diffuse.py

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_extended.py,v 1.78 2012/10/01 18:12:54 lande Exp $

    author: Joshua Lande
"""

import sys
from SpatialModels import RadiallySymmetricModel,Disk,SpatialModel,SpatialMap
from uw.utilities.convolution import ExtendedSourceConvolution,AnalyticConvolutionCache
from roi_diffuse import DiffuseSource,ROIDiffuseModel,SmallBand
from textwrap import dedent
from skymaps import SkyDir,Background
from scipy.optimize import fmin 
import numpy as np
import numbers
from uw.like.Models import PowerLaw
from uw.utilities import keyword_options
from . pypsf import PsfOverlap

class ExtendedSource(DiffuseSource):
    """ Class inherting from DiffuseSource but implementing a spatial source. 
        The main difference is the requirement of a spatial model to accomany 
        a spectral model. """

    defaults = (
        ('name',None,'The name of the extended source.'),
        ('model',None,'a Model object.'),
        ('spatial_model',None,"""The spatial model to use. 
                                       This is a SpatialModel object."""),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,**kwargs):
        """ Make the naming consistent with the PointSource object so that
            extended sources 'feel' like point sources.

            """
        keyword_options.process(self, kwargs)

        if self.model == None: self.model = PowerLaw()
        if self.spatial_model == None: self.spatial_model = Disk()

        if not isinstance(self.spatial_model,SpatialModel):
            raise Exception("The diffuse_model passed to an Extended Source must inherit from SpatialModel.")

        super(ExtendedSource,self).__init__(
            diffuse_model = self.spatial_model,
            scaling_model = self.model,
            name          = self.name)

        self.model.background = False

    @property
    def skydir(self): return self.spatial_model.center

    @property
    def smodel(self): 
        """ No reason to keep a model & smodel. """
        return self.model

    @smodel.setter
    def smodel(self, value): self.model = value


    def __str__(self,indent=''):
        return indent+('\n'+indent).join(['\n',
                          '='*60,
                          'Name:\t\t%s'%(self.name),
                          'R.A. (J2000):\t\t%.5f'%(self.spatial_model.center.ra()),
                          'Dec. (J2000):\t\t%.5f'%(self.spatial_model.center.dec()),
                          'Model:\t\t%s'%(self.model.full_name()),
                          '\t'+self.model.__str__(indent='\t'), 
                          'SpatialModel:\t%s'%(self.spatial_model.full_name()),
                          '\t'+self.spatial_model.__str__(indent='\t')
                         ])

    def copy(self):
        """ Create a deep copy of an extended source. """
        return ExtendedSource(name=self.name,
                              spatial_model=self.spatial_model.copy(),
                              model=self.model.copy())


###=========================================================================###


class ROIExtendedModel(ROIDiffuseModel):
    """ Implements the ROIDiffuseModel interface for the
        representation of a spatial source which can be represented as
        an analytic function inherting from the SpatialModel class."""

    @staticmethod 
    def factory(spectral_analysis,extended_source,*args,**kwargs):
        """ This static method acts just like the constructor to
            ROIExtendedModel and ROIExtendedModelAnalytic but decides on
            the fly whether or not the spatial model in the ExtendedSource
            is radially symmetric and returns the analytic convolution
            extended model for the radially symmetric source and the
            numeric convolution extended model for non-radially symmetric
            sources. """

        #if not isinstance(extended_source,ExtendedSource):
        #    raise Exception("The extended_source option passed to ROIExtendedModel.factory() must inherit from ExtendedSource.")

        spatial_model = extended_source.spatial_model

        if isinstance(spatial_model,RadiallySymmetricModel):
            return ROIExtendedModelAnalytic(spectral_analysis,extended_source,*args,**kwargs)
        elif isinstance(spatial_model,SpatialModel) or isinstance(spatial_model,SpatialMap):
            return ROIExtendedModel(spectral_analysis,extended_source,*args,**kwargs)
        else:
            raise Exception("The extended_source.dmodel option passed to ROIExtendedModel.factory must inherit from SpatialModel.")

    def setup(self):
        """ Unlike background models, always do the convolution around 
            the spatial model's center. """
        self.exp = self.sa.exposure.exposure; psf = self.sa.psf
        self.active_bgc = ExtendedSourceConvolution(self.extended_source,psf)

    def set_state(self,band):
        self.active_bgc.do_convolution(band)
        self.current_energy = energy=band.psf.eopt if band.psf.__dict__.has_key('eopt') else band.e

    def initialize_counts(self,bands,roi_dir=None):
        rd = self.roi_dir if roi_dir is None else roi_dir
        self.bands = [SmallBand() for i in xrange(len(bands))]

        es = self.extended_source
        sm = es.model

        for iband,(myband,band) in enumerate(zip(self.bands,bands)):
            if not self.quiet: 
                status_string = '...convolving band %2d/%2d'%(iband+1,len(self.bands))
                print (status_string,)
                sys.stdout.flush()

            self.set_state(band)

            exposure=band.exp.value

            myband.er = exposure(es.spatial_model.center,self.current_energy)/exposure(rd,self.current_energy)

            myband.es_pix_counts = self._pix_value(band.wsdl)*band.b.pixelArea()

            myband.overlaps = self._overlaps(rd,band)

            if not self.quiet: 
                print ('\b'*(2+len(status_string)))
                sys.stdout.flush()

        if not self.quiet: print()

    def update_counts(self,bands,model_index):
        """Update models with free parameters."""
        sm = self.extended_source.model
        mi = model_index

        for myband,band in zip(self.bands,bands):
            myband.es_counts = band.expected(sm)*myband.er

            band.bg_counts[mi] = myband.overlaps*myband.es_counts
            if band.has_pixels:
                band.bg_pix_counts[:,mi] = myband.es_pix_counts * myband.es_counts

    def gradient(self,bands,model_index):
        """ Return the gradient of the log likelihood with respect to the spectral parameters of
            this model. Note that this calculation is essentially identical to that of point
            sources since extended sources decouple the spectral and spatial parts, as is done
            for point soruces. """
        sm  = self.extended_source.model
        np  = sm.len() 
        nfp = sm.free.sum()

        # special case -- no free parameters
        if nfp == 0: return []

        gradient = [0]*nfp

        for myband,band in zip(self.bands,bands):

            grad    = band.gradient(sm)[sm.free]*myband.er # correct for exposure
            apterm =  band.phase_factor*myband.overlaps
            if band.has_pixels:
                pixterm = (band.pix_weights*myband.es_pix_counts).sum()
            else:
                pixterm = 0
            gradient += grad * (apterm - pixterm)

        return gradient

    def _pix_value(self,pixlist):
        return self.active_bgc(pixlist,self.active_bgc.cvals)

    def _overlaps(self,center,band=None):
        return self.active_bgc.overlap(center,band.radius_in_rad)

    def __init__(self,spectral_analysis,extended_source,roi_dir,name=None,*args,**kwargs):

        self.extended_source = self.diffuse_source = extended_source

        super(ROIExtendedModel,self).__init__(
            spectral_analysis, extended_source,
            roi_dir, extended_source.name,*args,**kwargs)

    def __str__(self):
        es = self.extended_source
        sm = es.spatial_model

        return '%s fitted with %s\n%s\n%s fitted with %s\n%s' % \
                (es.name,sm.full_name(),
                 '\t'+sm.__str__(indent='\t'),
                 es.name,es.model.full_name(),
                 '\t'+es.model.__str__(indent='\t'))

    def fit_extension(self,roi,tolerance=0.05, maxcalls=500,
                      bandfits=False, error="UMINOS",init_grid=None, estimate_errors=True, 
                      verbose=False,
                      **kwargs):
        """ Fit the extension of this extended source by fitting all non-fixed spatial paraameters of 
            self.extended_source. The likelihood at the best position is returned.

Arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  roi                   An ROIAnalysis object.
  tolerance             This is passed into Minuit while fitting the
                        extension parameters.
  bandfits      [False] Whether to use a spectral independent model when
                        fitting extension.
  error       ["HESSE", "UMINOS"] 
                        The fitting algorithm to use when calculating
                        errors.  UMINOS will call uw.utilitites.minuit
                        function Minuit.uncorrelated_minos_error function.
                        Specify None to not calculate errors (or set
                        estimate_errors=False).
  init_grid      [None] A list of spatial parameters. The likelihood
                        for each of the spatial parametesr is tested and
                        the minimum one is used as the initial guess
                        when running minuit. Useful for doing initial
                        grid search.

                        - The spatial part of init_grid
                          should be measured relative to the spatial model's center: (0,0)
                        - The init_grid values should be absolute, none of them should
                          be given in log space.
  **kwargs             Will be passed into the ROIAnalysis.fit() function.
  =========   =======================================================

        Implemenation notes, do not directly fit RA and Dec. Instead,
        rotate the coordiante to (0,0) and fit x & y deviation (measured
        in degrees) in this rotated coordinate system. For each function
        iteration, rotate back to the original space and set the direction
        to there.  The error on RA (Dec) is thus the distance error
        in the direction of increasing RA (Dec). This makes the method
        robust anywhere in the sky. """

        from uw.utilities.minuit import Minuit
        from . roi_state import PointlikeState

        es = self.extended_source

        sm = es.spatial_model
        cs = sm.coordsystem

        origin=SkyDir(0,0,cs)

        # keep the roi quiet during fit, so lots of print is supressed during the iteration.
        roi.old_quiet = roi.quiet
        roi.quiet = True

        if roi.TS(which=self.name,quick=True,bandfits=bandfits) < 1:
            print ("Warning: initial TS<1 (in point hypothesis) so TS_ext will likely not be trustworthy")

        init_state = PointlikeState(roi)
        init_spatial  = sm.get_parameters(absolute=False)

        # Fix scoping issue in likelihood_wrapper
        # In python 3, could use 'nonlocal' parameter.
        # http://eli.thegreenplace.net/2011/05/15/understanding-unboundlocalerror-in-python/
        d=dict(best_spectral = PointlikeState(roi), ll_best = -np.inf, ll_0=-np.inf)

        if not np.any(sm.free):
            raise Exception("Unable to fit the diffuse source %s's extension. No parameters to fit. Perhaps when you wrote your xml file, you forgot to set some of the spatial parameters to be free?" % es.name)

        # Remember thet get_parameters only returns the free parameters.
        # here we need to get the first two parameters (position) so
        # we can fit relative to it.
        init_lon,init_lat = sm.p[0:2]
        init_dir          = SkyDir(init_lon,init_lat,cs)


        # Define an axis and angle that takes the equator point (of the
        # same longitude) to the original position.
        # Fit values around at the equator around this new point.
        axis=SkyDir(init_lon-90,0,cs)
        theta=np.radians(init_lat)

        # Fit in the rotated coodinate system.
        init_spatial[0:2]=[init_lon,0]

        def likelihood_wrapper(p):
            """ Helper function which takes in the spatial parameters
                for the extended source and returns the logLikelihood.

                Note: Sometimes the fitter gets confused. For example,
                if the last iteration moved to a spatial value really far
                from the initial position, the model will very pooorly
                converge and during the next fill be unable to move
                back to reasonable spectral values.  To get around this,
                the code remembers the initial and best spectral values
                and tries those out whenever the overall likelihood is
                worse than it was before.  """
            p=p.copy()

            lon,lat=p[0:2]
            # New direction in rotated coordiante system
            new_dir = SkyDir(lon,lat,cs)


            # This rotation takes (0,0) back to the initial position
            new_dir().rotate(axis(),theta)

            # update spatial model, then modify inside the ROI.
            sm.set_parameters(p=p[2:],absolute=False, center=new_dir)

            if verbose:
                print ('Extension fitting for spatial model with p(absolute=False)=',sm.get_parameters(absolute=False).tolist())

            temp=self.quiet;self.quiet=True
            self.initialize_counts(roi.bands)
            self.quiet=temp
            roi.update_counts()

            if bandfits:
                ll=roi.bandFit(es)
            else:
                if verbose:
                    print ('About to do first fit:')
                    roi.print_summary(indent='    ')

                ll=roi.fit(estimate_errors=False,**kwargs)

                if verbose:
                    print ('After first fit:')
                    roi.print_summary(indent='    ')

                if ll < d['ll_best']:

                    prev_state=PointlikeState(roi)
                    d['best_state'].restore(just_spectra=True)

                    if verbose:
                        print ('likelihood=%s worse than best likelihood=%s, starting from best spectral parameters' % (ll, d['ll_best']))
                        roi.print_summary(indent='    ')

                    ll_alt=roi.fit(estimate_errors=False,**kwargs)

                    if verbose:
                        print ('After fit with best parameters:')
                        roi.print_summary(indent='    ')

                    if ll_alt > ll: 
                        if verbose: 
                            print ('Likelihood better than previous likelihood, so keeping it')
                        ll=ll_alt
                    else: 
                        if verbose:
                            print ('Likelihood worse than previous likelihood, so discarding it')
                        prev_state.restore(just_spectra=True)
                else:
                    if verbose:
                        print ('likelihood=%s is better than best likelihood=%s' % (ll, d['ll_best']))

                if ll < d['ll_0']:

                    prev_state=PointlikeState(roi)
                    init_state.restore(just_spectra=True)

                    if verbose:
                        print ('Likelihood=%s worse than initial likelihood=%s, starting from initial spectral parameters' % (ll, d['ll_0']))
                        roi.print_summary(indent='    ')

                    ll_alt=roi.fit(estimate_errors=False,**kwargs)

                    if verbose:
                        if verbose:
                            print ('After fit with initial parameters:')
                        roi.print_summary(indent='    ')

                    if ll_alt > ll: 
                        if verbose:
                            print ('Likelihood better than previous likelihood, so keeping it')
                        ll=ll_alt
                    else: 
                        if verbose:
                            print ('Likelihood worse than previous likelihood, so discarding it')
                        prev_state.restore(just_spectra=True)
                else:
                    if verbose:
                        print ('likelihood=%s is better than initial likelihood=%s' % (ll, d['ll_0']))

                if ll > d['ll_best']: 
                    d['ll_best'] = ll
                    d['best_state'] = PointlikeState(roi)

            if not self.quiet: print ('%s, logL = %.3f, dlogL = %.3f' % (sm.pretty_string(),ll,ll-d['ll_0']))

            if verbose: 
                print();print()

            return -ll # return negative log likelihood, suitable for minimizing

        f=likelihood_wrapper


        # shut up print just the first time
        old_quiet = self.quiet; self.quiet = True; 
        d['ll_0'] = d['ll_best'] = -f(init_spatial); 
        self.quiet = old_quiet

        if not self.quiet: 
            print ('Localizing %s source %s Using %s' % (sm.pretty_name,es.name,
                                                        'BandFits' if bandfits else 'Spectral Fits'))

        if init_grid is not None:
            print ('Testing initial grid values')
            ll = [] 
            transformed = []
            for _ in init_grid:
                # convert into log space.
                sm.set_parameters([0,0,_] if isinstance(_,numbers.Real) else _,absolute=True)
                param=sm.get_parameters(absolute=False)
                param[0] += init_lon # rotate longitude to fit area
                # easy way to turn parameters into their non-absolute value.
                ll.append(likelihood_wrapper(param))
                transformed.append(param)

            index=np.argmin(ll)
            start_spatial=transformed[index]
            if not self.quiet: print ('Now starting with the best initial parameters')
        else:
            start_spatial = init_spatial

        # Display which parameters are in log space.
        relative_names = sm.get_param_names(absolute=False)

        limits=sm.get_limits(absolute=False)
        limits[0]+= init_lon

        m = Minuit(f,start_spatial,
                   up=0.5,
                   maxcalls=maxcalls,
                   tolerance=tolerance,
                   printMode = -1 if self.quiet else 1,
                   param_names=relative_names,
                   limits    = limits,
                   fixed     = ~sm.free,
                   steps     = sm.get_steps())

        best_spatial,fval = m.minimize(method="SIMPLEX")

        if estimate_errors is True and error is not None:
            if not self.quiet: print ('Calculating Covariance Matrix')
            if error == 'UMINOS':
                cov_matrix = m.uncorrelated_minos_error()
            else:
                cov_matrix = m.errors(method=error)
        else:
            cov_matrix = np.zeros([len(sm.p),len(sm.p)])

        sm.set_cov_matrix(cov_matrix)

        if not self.quiet: print ('setting source to best fit parameters')
        likelihood_wrapper(best_spatial)

        # Fit once more at the end to get the right errors.
        try:
            roi.fit(estimate_errors=True,**kwargs)
        except Exception, err:
            print (err)

        roi.quiet = roi.old_quiet

        # return log likelihood from fitting extension.
        final_dir=sm.center
        delt = np.degrees(final_dir.difference(init_dir))
        return final_dir,0,delt,-2*(d['ll_0']+fval)

    def modify_loc(self,bands,center):
        self.extended_source.spatial_model.modify_loc(center)
        self.initialize_counts(bands)
    
    def TS_ext(self,roi,refit=True,**kwargs):
        """ Calculates TS_ext for an extended source, which is defined as 2*(ll_ext-ll_pt).
            This code does so by shinking the extended source down to be very small and
            relocalizing the shrunk source to calculate the likelihood in the null hypothesis.
            After calculating TS_ext, the ROI is reset to the state before the function is called.
        
            If refit is false, the null hypothesis (shrunken source) is not relocalized.
            Generally, this is a great idea since it gives a biased estimate of significance,
            but may be useful to quickly get approximate values

            Any additional argument passed into this function is passed into fit_extension when
            the null hypothesis is relocalized. """
        if kwargs.has_key('error') or kwargs.has_key('update'): 
            raise Exception("argument 'error' cannot be passed into function TS_ext")

        es = self.extended_source
        sm = es.spatial_model

        if not sm.can_shrink():
            raise Exception("Unable to calculate ts_ext for %s source %s. No way to shrink to null hypothesis." % (sm.pretty_name,es.name))

        old_roi_p   = roi.get_parameters().copy()

        def l():
            if kwargs.has_key('bandfits') and kwargs['bandfits']:
                return roi.bandFit(es)
            else:
                return -roi.logLikelihood(roi.get_parameters())

        def f():
            d={'use_gradient':kwargs['use_gradient']} if kwargs.has_key('use_gradient') else {}
            try:
                roi.fit(estimate_errors=True,**d)
            except Exception, err:
                print (err)

        ll_disk = l()

        sm.shrink()
        self.initialize_counts(roi.bands)

        if not roi.quiet: print ('Refitting position for the null hypothesis')
        f() # have to fit with shrunk spatial model to get a reasonable starting spectrum for extension fit.
        if refit: self.fit_extension(roi,estimate_errors=False,**kwargs)

        if not roi.quiet: print ('Redoing spectral fit in the null hypothesis')

        f() # have to refit in case bandfits was used during the localization.

        ll_point = l()

        ts_ext = 2*(ll_disk - ll_point)

        sm.unshrink()
        self.initialize_counts(roi.bands)

        roi.set_parameters(old_roi_p) 
        roi.__set_error__()

        # reset model predictions
        roi.__update_state__()

        return ts_ext

###=========================================================================###


class ROIExtendedModelAnalytic(ROIExtendedModel):
    """ Implements the ROIDiffuseModel interface for a radially
        symmetric extended source. Utilize teh semi-analytic
        convolution for a more efficient pdf calculation.

        Any of the optional keyword arguments to uw.utilities.convolution's 
        AnalyticConvolutionCache class will be passed on to that class.  """

    defaults = ROIExtendedModel.defaults

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        super(ROIExtendedModelAnalytic,self).__init__(*args,**kwargs)
        self.overlap = PsfOverlap()

    def setup(self):
        self.exp = self.sa.exposure.exposure; 
        psf = self.sa.psf

        d={'num_points':self.num_points} if self.__dict__.has_key('num_points') else {}
        self.active_bgc = AnalyticConvolutionCache(self.extended_source,psf,**d)

    def set_state(self,band):
        self.active_bgc.do_convolution(band)
        self.current_energy = energy=band.psf.eopt if band.psf.__dict__.has_key('eopt') else band.e

    def _pix_value(self,pixlist):
        return self.active_bgc(pixlist)

    def _overlaps(self,center,band,radius=None):
        """ Calculate the fraction of the PDF not contained within the ROI
            using the PsfOverlap object (but override the pdf and integral
            function to use instead the extended source pdf. """
        return self.overlap(band=band,
                            roi_dir=center,
                            ps_dir=self.extended_source.spatial_model.center,
                            radius_in_rad=radius if radius is not None else band.radius_in_rad,
                            override_pdf=self.active_bgc,
                            override_integral=self.active_bgc.integral)

class BandFitExtended(object):

    def __init__(self,which,energy_band,roi):
        """ extendedsource
            which             index of source
            energy_band       ROIEnergyBand object to fit.
            bands             all energy bands in roi. """

        self.energy_band = energy_band
        self.bands       = self.energy_band.bands
        self.all_bands   = roi.bands

        self.which = which

        self.all_mybands = roi.dsm.bgmodels[self.which].bands

        # list of the mybands corresponding to energy_bands
        self.mybands = []

        # Create a lsit of myband object corresponding to the bands 
        # in the energy_band.
        for eb_band in self.bands:
            for band,myband in zip(self.all_bands,self.all_mybands):
                if eb_band == band:
                    self.mybands.append(myband)
                    break

    def bandLikelihoodExtended(self,parameters, band, myband):

        new_counts = parameters[0]*myband.er

        old_counts = band.bg_counts[self.which]

        tot_term = (band.bg_all_counts + band.ps_all_counts + myband.overlaps*(new_counts - old_counts))*band.phase_factor

        pix_term = (band.pix_counts * 
                            np.log(
                                band.bg_all_pix_counts + band.ps_all_pix_counts + myband.es_pix_counts*(new_counts - old_counts)
                            )
                      ).sum() if band.has_pixels else 0.

        return tot_term - pix_term

    def energyBandLikelihoodExtended(self,parameters,m):
        m.set_parameters(parameters)
        return sum(self.bandLikelihoodExtended([b.expected(m)],b,mb) for b,mb in 
                   zip(self.bands,self.mybands))

    def normUncertaintyExtended(self):
        tot = 0
        for b,mb in zip(self.bands,self.mybands):
            if not b.has_pixels: continue
            my_pix_counts = mb.es_pix_counts*b.expected(self.m)*mb.er
            all_pix_counts= b.bg_all_pix_counts + b.ps_all_pix_counts - mb.es_pix_counts*b.bg_counts[self.which] + my_pix_counts
            tot += (b.pix_counts * (my_pix_counts/all_pix_counts)**2).sum()
        return tot**-0.5

    def fit(self,saveto=None):

        bad_fit = False
        self.m = PowerLaw(free=[True,False],e0=(self.energy_band.emin*self.energy_band.emax)**0.5) # fix index to 2
        f = self.energyBandLikelihoodExtended

        self.fit = fmin(f,self.m.get_parameters(),disp=0,full_output=1,args=(self.m,))

        def upper_limit():

            flux_copy = self.m[0]
            zp          = self.energyBandLikelihoodExtended(np.asarray([-20]),self.m)

            # NB -- the 95% upper limit is calculated by assuming the likelihood is peaked at
            # 0 flux and finding the flux at which it has fallen by 1.35; this is a two-sided
            # 90% limit, or a one-sided 95% limit -- that's how it works, right?
            def f95(parameters):
                return abs(self.energyBandLikelihoodExtended(parameters,self.m) - zp - 1.35)
            
            # for some reason, can't get fsolve to work here.  good ol' fmin to the rescue
            self.energy_band.uflux = 10**fmin(f95,np.asarray([-11.75]),disp=0)[0]
            self.energy_band.lflux = None
            self.energy_band.flux  = None

            self.m[0] = flux_copy

        # if flux below a certain level, set an upper limit
        if self.m[0] < 1e-20:
            bad_fit = True
            upper_limit()

        else:
            try:
                err = self.normUncertaintyExtended()
            except:
                bad_fit = True
                err = 0 

            self.energy_band.flux  = self.m[0] 
            self.energy_band.uflux = self.energy_band.flux*(1 + err)
            self.energy_band.lflux = max(self.energy_band.flux*(1 - err),1e-30)

        if saveto is not None:
            for b,mb in zip(self.bands,self.mybands): 
                b.__dict__[saveto] = (b.expected(self.m)*mb.er if not bad_fit else -1)

        if bad_fit:
            self.energy_band.ts = 0
        else:
            null_ll = sum(self.bandLikelihoodExtended([0],b,mb)
                for b,mb in zip(self.bands,self.mybands))
            alt_ll  = sum(self.bandLikelihoodExtended([b.expected(self.m)*mb.er],b,mb)
                    for b,mb in zip(self.bands,self.mybands))
            self.energy_band.ts = 2*(null_ll - alt_ll)

