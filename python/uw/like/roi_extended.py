""" Provides classes to encapsulate and manipulate extended sources. 

    This code all derives from objects in roi_diffuse.py

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_extended.py,v 1.49 2011/04/04 22:56:25 kerrm Exp $

    author: Joshua Lande
"""

from SpatialModels import RadiallySymmetricModel,Disk,SpatialModel,SpatialMap
from uw.utilities.convolution import ExtendedSourceConvolutionCache,AnalyticConvolutionCache
from roi_diffuse import DiffuseSource,ROIDiffuseModel,SmallBand
from textwrap import dedent
from skymaps import SkyDir,Background
from scipy.optimize import fmin 
import numpy as N
from uw.like.Models import PowerLaw
from uw.utilities import keyword_options
from uw.like.pypsf import PsfOverlap

class ExtendedSource(DiffuseSource):
    """ Class inherting from DiffuseSource but implementing a spatial source. 
        The main difference is the requirement of a spatial model to accomany 
        a spectral model. """

    defaults = (
        ('name',None,'The name of the extended source.'),
        ('model',None,'a Model object.'),
        ('spatial_model',None,"""The spatial model to use. 
                                       This is a SpatialModel object."""),
        ('free_parameters',True,"""Set all spectral model and 
                                   spatial_model parameter's free 
                                   status to this value"""),
        ('leave_parameters',False,"""Don't set any of the spectral 
                                     or spatial parameter's free status. 
                                     Leave them at whatever their default 
                                     values are.""")
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

        if not self.leave_parameters:
            for i in xrange(len(self.model.free)): self.model.free[i] = self.free_parameters
            for i in xrange(len(self.spatial_model.free)): self.spatial_model.free[i] = self.free_parameters

    @property
    def skydir(self): return self.spatial_model.center

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
        self.active_bgc = ExtendedSourceConvolutionCache(self.extended_source,psf)

    def set_state(self,band):
        self.active_bgc.do_convolution(band)
        self.current_energy = energy=band.psf.eopt if band.psf.__dict__.has_key('eopt') else band.e

    def initialize_counts(self,bands,roi_dir=None):
        rd = self.roi_dir if roi_dir is None else roi_dir
        self.bands = [SmallBand() for i in xrange(len(bands))]

        es = self.extended_source
        sm = es.model

        for myband,band in zip(self.bands,bands):

            self.set_state(band)

            exposure=band.exp.value

            myband.er = exposure(es.spatial_model.center,self.current_energy)/exposure(rd,self.current_energy)

            myband.es_pix_counts = self._pix_value(band.wsdl)*band.b.pixelArea()

            myband.overlaps = self._overlaps(rd,band)

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

    def fit_extension(self,roi,tolerance=0.05, bandfits=False, error="HESSE",init_grid=None, use_gradient=True, estimate_errors=True):
        """ Fit the extension of this extended source by fitting all non-fixed spatial paraameters of 
            self.extended_source. The likelihood at the best position is returned.

Arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  roi                   An ROIAnalysis object.
  tolerance             This is passed into Minuit while fitting the extension parameters.
  bandfits      [False] Whether to use a spectral independent model when fitting extension.
  use_gradient  [False] If bandfits is False, whether to use the analytic gradient when fitting 
                        spectral parameters (during an iteration for a given extension).
  error       ["HESSE"] The fitting algorithm to use when calculating errors.
                        Specify None to not calculate errors.
  init_grid      [None] A list of spatial parameters. The likelihood
                        for each of the spatial parametesr is tested and the minimum
                        one is used as the initial guess when running minuit. Useful
                        for doing initial grid search. 

                        - The spatial part of init_grid
                          should be measured relative to the spatial model's center: (0,0)
                        - The init_grid values should be absolute, none of them should
                          be given in log space.
  =========   =======================================================

        Implemenation notes, do not directly fit RA and Dec. Instead,
        rotate the coordiante to (0,0) and fit x & y deviation (measured
        in degrees) in this rotated coordinate system. For each function
        iteration, rotate back to the original space and set the direction
        to there.  The error on RA (Dec) is thus the distance error
        in the direction of increasing RA (Dec). This makes the method
        robust anywhere in the sky. """

        from uw.utilities.minuit import Minuit

        es = self.extended_source

        sm = es.spatial_model
        cs = sm.coordsystem

        origin=SkyDir(0,0,cs)

        # keep the roi quiet during fit, so lots of print is supressed during the iteration.
        quiet = roi.quiet
        roi.quiet = True

        if roi.TS(which=self.name,quick=True,bandfits=bandfits) < 1:
            raise Exception("Unable to localize a source with initial TS<1")

        init_spectral = es.model.get_parameters()
        init_spatial = sm.get_parameters(absolute=False)

        if not N.any(sm.free):
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
        theta=N.radians(init_lat)

        # Fit in the rotated coodinate system.
        init_spatial[0:2]=[init_lon,0]

        def likelihood_wrapper(p):
            """ Helper function which takes in the spatial parameters
                for the extended source and returns the logLikelihood.

                Note: Sometimes the fitter gets confused. For example,
                if the last iteration moved to a spatial value really far
                from the initial position, the source flux will get fit
                to 0 (as it should), but during the next iteration at a
                more reasonable spatial value, the fit will be unable
                to reconverge on the best flux and stays at 0. To get
                around this, I cache the initial fit value (which is
                presumably pretty good), and whenever the current fit
                logLikelihood is less then the initial logLikelihood,
                I rest the loglikelihood to
                
                N.B. roi.dsm.update_counts handled by fit() calling logLikelihood()
                """
            p=p.copy()

            lon,lat=p[0:2]
            # New direction in rotated coordiante system
            new_dir = SkyDir(lon,lat,cs)


            # This rotation takes (0,0) back to the initial position
            new_dir().rotate(axis(),theta)

            # Now add the rotated spatial part back to the list.
            if cs == SkyDir.GALACTIC:
                p[0:2]=new_dir.l(),new_dir.b()
            elif cs == SkyDir.EQUATORIAL:
                p[0:2]=new_dir.ra(),new_dir.dec()

            # Do the convolution here.
            sm.set_parameters(p=p,absolute=False)
            self.initialize_counts(roi.bands)
            roi.update_counts()
            # Note: roi.dsm.update_counts called by the fit function.


            if bandfits:
                ll=roi.bandFit(es)
            else:
                ll=roi.fit(estimate_errors=False,use_gradient=use_gradient)

                if ll < ll_0:
                    prev_fit=es.model.get_parameters()
                    es.model.set_parameters(init_spectral)
                    ll_alt=roi.fit(estimate_errors=False,use_gradient=use_gradient)

                    if ll_alt > ll: ll=ll_alt
                    else: es.model.set_parameters(prev_fit)

            if not quiet: print '%s, logL = %.3f, dlogL = %.3f' % (sm.pretty_string(),ll,ll-ll_0)
            return -ll

        f=likelihood_wrapper

        ll_0 = 0
        old_quiet = quiet; quiet = True; 
        ll_0 = -f(init_spatial); 
        quiet = old_quiet

        if quiet: print 'Localizing %s source %s Using %s' % (sm.pretty_name,es.name,
                                                    'BandFits' if bandfits else 'Spectral Fits')

        if init_grid is not None:
            print 'Testing initial grid values'
            ll = [] 
            transformed = []
            for _ in init_grid:
                # convert into log space.
                sm.set_parameters(_,absolute=True)
                param=sm.get_parameters(absolute=False)
                param[0] += init_lon # rotate longitude to fit area
                # easy way to turn parameters into their non-absolute value.
                ll.append(likelihood_wrapper(param))
                transformed.append(param)

            index=N.argmin(ll)
            start_spatial=transformed[index]
            if quiet: print 'Now starting with the best initial parameters'
        else:
            start_spatial = init_spatial

        # Display which parameters are in log space.
        relative_names = sm.get_param_names(absolute=False)

        limits=sm.get_limits(absolute=False)
        limits[0]+= init_lon

        m = Minuit(f,start_spatial,
                   up=0.5,
                   maxcalls=500,
                   tolerance=tolerance,
                   printMode = -1 if quiet else 1,
                   param_names=relative_names,
                   limits    = limits,
                   fixed     = ~sm.free,
                   steps     = sm.get_steps())

        best_spatial,fval = m.minimize(method="SIMPLEX")

        if estimate_errors is True and error is not None:
            if not quiet: print 'Calculating Covariance Matrix'
            cov_matrix = m.errors(method=error)
        else:
            cov_matrix = N.zeros([len(sm.p),len(sm.p)])

        sm.set_cov_matrix(cov_matrix)

        if not quiet: print 'setting source to best fit parameters'
        likelihood_wrapper(best_spatial)

        # Fit once more at the end to get the right errors.
        roi.fit(estimate_errors=True,use_gradient=use_gradient)

        roi.quiet = quiet

        # return log likelihood from fitting extension.
        final_dir=sm.center
        delt = N.degrees(final_dir.difference(init_dir))
        return final_dir,0,delt,-2*(ll_0+fval)

    def modify_loc(self,bands,center):
        self.extended_source.spatial_model.modify_loc(center)
        self.initialize_counts(bands)
    
    def TS_ext(self,roi,refit=True,**kwargs):
        """ Refit refers to wheter the localization of the extended source should
            be relocalized for the null hypothesis. Generally, this is a great idea
            so the default is to do so.

            Any additional argument passed into this function is passed into fit_extension when
            the null hypothesis is localized. 
            
            Neither 'error' or 'update' are allowed to be passed into fit_extension,
            because it is of no use to calculate for the null hypothesis. """
        if kwargs.has_key('error') or kwargs.has_key('update'): 
            raise Exception("argument 'error' cannot be passed into function TS_ext")

        es = self.extended_source
        sm = es.spatial_model

        old_sm_p       = sm.p.copy()
        old_sm_cov     = sm.cov_matrix.copy()
        old_sm_free    = sm.free.copy()
        old_roi_p   = roi.get_parameters().copy()

        def l():
            if kwargs.has_key('bandfits') and kwargs['bandfits']:
                return roi.bandFit(es)
            else:
                return -roi.logLikelihood(roi.get_parameters())

        def f():
            d={'use_gradient':kwargs['use_gradient']} if kwargs.has_key('use_gradient') else {}
            roi.fit(estimate_errors=True,**d)

        ll_disk = l()

        try:
            sm.shrink()
        except:
            raise Exception("Unable to calculate ts_ext for %s source %s. No way to shrink to null hypothesis." % (sm.pretty_name,es.name))
        self.initialize_counts(roi.bands)

        if not roi.quiet: print 'Refitting position for the null hypothesis'
        f() # have to fit with shrunk spatial model to get a reasonable starting spectrum for extension fit.
        if refit: self.fit_extension(roi,error=None,**kwargs)

        if not roi.quiet: print 'Redoing spectral fit in the null hypothesis'

        f() # have to refit in case fitpsf or bandfits was used during the localization.

        ll_point = l()

        ts_ext = 2*(ll_disk - ll_point)

        sm.set_parameters(p=old_sm_p,absolute=False)
        sm.cov_matrix = old_sm_cov
        sm.free = old_sm_free

        roi.set_parameters(old_roi_p) 
        roi.__set_error__()

        self.initialize_counts(roi.bands)
        # reset point source

        roi.__update_state__()

        return ts_ext

###=========================================================================###


class ROIExtendedModelAnalytic(ROIExtendedModel):
    """ Implements the ROIDiffuseModel interface for a radially
        symmetric extended source. Utilize teh semi-analytic
        convolution for a more efficient pdf calculation.

        Any of the optional keyword arguments to uw.utilities.convolution's 
        AnalyticConvolutionCache class will be passed on to that class.  """

    defaults = [
        ['fitpsf',False,'Method where the psf is calculated by fitting a single king function to the actual psf.'],
        ['fastpsf',False,'Method where the psf is apprximated by weighting the sigmas and gammas before evaluating the psf.']
    ]

    @staticmethod
    def set_fastpsf(fastpsf):
        was = ROIExtendedModelAnalytic.defaults[1][1]
        ROIExtendedModelAnalytic.defaults[1][1] = fastpsf
        return was

    @keyword_options.decorate(defaults)
    def __init__(self,*args,**kwargs):
        super(ROIExtendedModelAnalytic,self).__init__(*args,**kwargs)
        self.already_fit=False
        self.overlap = PsfOverlap()

    def setup(self):
        self.exp = self.sa.exposure.exposure; 
        psf = self.sa.psf

        d={'num_points':self.num_points} if self.__dict__.has_key('num_points') else {}
        self.active_bgc = AnalyticConvolutionCache(self.extended_source,psf,**d)

    def set_state(self,band):
        self.active_bgc.do_convolution(band,self.fitpsf, self.fastpsf)
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

    def initialize_counts(self,bands,roi_dir=None):
        # probably there is a cleaner way to do this, but since you should only 
        # use this feature in the process of localizing one source, it is 
        # probably good enough.

        if self.fitpsf and not self.already_fit: 
            # Note that this must be done before calculating the pdf.
            fitter=BandFitter(bands)
            fitter.fit()
            self.already_fit=True

        super(ROIExtendedModelAnalytic,self).initialize_counts(bands,roi_dir)

    def fit_extension(self,*args,**kwargs):
        """ Localiztion of radially symmetric extended sources is exactly the same as
            regular extended sources. The only difference is radially symmetric
            sources can be fit with the addition of an optional flag 'fitpsf'.

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  fitpsf        [True]  Use an approximation to the PSF in each energy bin
                        where PSF is weighted in each energy bin by an assumed
                        spectral index and is fit by a single king function
                        which is then convovled with the extended source shape.
                        This leads to a significantly faster source localization.
  =========   =======================================================
        """
        if kwargs.has_key('roi'):
            roi = kwargs['roi']
        else:
            roi=args[0]

        fitpsf=True if self.fastpsf!=True else False
        if kwargs.has_key('fitpsf'):
            fitpsf = kwargs.pop('fitpsf')

        if fitpsf:
            if self.fastpsf: raise Exception("Cannot specify fitpsf for localization since fastpsf is already selected.")

            if not roi.quiet: print 'Changing to fitpsf for localization step.'
            self.fitpsf,old_fitpsf=True,self.fitpsf

        stuff = super(ROIExtendedModelAnalytic,self).fit_extension(*args,**kwargs)

        if fitpsf:
            if not roi.quiet: print 'Setting back to original PDF accuracy.'
            self.fitpsf=old_fitpsf
            self.initialize_counts(roi.bands)
            roi.update_counts()
        return stuff

class BandFitter(object):
    """ This class has a somewhat weird purpose. Basically the PSF is
        a king function, (or for the newpsf two king functions). But
        this psf shape is a function of energy and the incident theta
        angle of the photos. When studying extended sources, we have a
        semianalytic formula to convolve a king function with the shape
        we are interested in studying. On the other hand, pointlike
        does a binned analysis binning in course energy bins and not
        in photon theta angle. So when calculating an extended source
        model counts, we must convolve our spatial model with the PSF
        for different theta and energy values in a particular energy
        bin and then integrate over energy & theta weighting by the LAT's
        exposure and the spectral model's spectrum. For point sources,
        the PSF is integrated over 16 sub energy values and 8 theta
        values for each energy bin to get the model predicted counts.
        But to study extended sources, for each band doing 16x8 of these
        convolutions becomes untenable in the process of fitting the
        spatial extension of the source.

        To work around this, we are interested in 'faking' by calculating
        what the psf truly should by correctly weighting the different
        energies and theta values in an energy bin with an assumed
        spectrum (defaulting to a powerlaw index 2) and then fitting a
        single king function to this psf, which is considered the energy
        bin average of the psf and should be a reasonable representation
        of the width of the source, for the sake of fitting the extension
        of a source.

        N.B. For the newstyle psf which is represtened by the sum of
        two king functions, it is best to weight each king function
        independently and end up with an average core and average tail
        PSF for the particular energy bin.  """

    def __init__(self,bands):
        """ Pass in the roi.bands object. """
        self.bands=bands

    def fit(self,weightspectrum=PowerLaw()):
        """ Fit the psf. the fit values will be stored inside of
            each band in the variables band.fit_gamma and band.fit_sigma,
            or for the newstyle psf as band.fit_gamma_core
            and band.fit_sigma_core & band.fit_gamma_tail and
            band.fit_sigma_tail. """

        print 'Fitting PSF to weighted average'
        for band in self.bands: 
            self._fit_band(band,weightspectrum)

    def psf_base(self,g,s,delta):
        # Fancy vectorized code taken from BandCALDBPsf. Maybe one
        # day these functions can share a common base.
        u = 0.5 * N.outer(delta,1./s)**2
        y = (1-1./g)*(1+u/g)**(-g)/(2*N.pi*s**2)
        return y

    def _fit_band(self,band,weightspectrum):
        """ Fit a psf with a default spectral index of 2. """

        psf=band.sa.psf

        newstyle = psf.newstyle

        # Get a rough estimate of the gamma & sigma
        # in the middle by weighting sigmas & gammas in the 
        # geomeric mean energy over cos theta.
        # This serves as the starting fit value.
        scale = psf.scale_func[band.ct](band.e)
        if newstyle:
            nc,nt,gc,gt,sc,st,w = psf.get_p(band.e,band.ct)
            gamma_middle=sum(w*(nc*gc+nt*gt))
            sigma_middle=sum(w*(nc*sc+nt*st))
        else:
            g,s,w=psf.get_p(band.e,band.ct)
            gamma_middle=sum(g*w)
            sigma_middle=sum(s*w)

        sigma_middle=sigma_middle*scale

        rlist = N.linspace(0,20*sigma_middle,10000)
        pdf = N.zeros_like(rlist)
        pdf_weight = N.zeros_like(rlist)

        weight_sum,s_list=0,[]

        # Use the same simpson integral used to calculate
        # predicted counts for a point source.
        for e,sp in zip(band.sp_points,band.sp_vector):
            scale = psf.scale_func[band.ct](e)
            if newstyle:
                nc,nt,gc,gt,sc,st,w = psf.get_p(e,band.ct)
                sc *= scale
                st *= scale
                pdf += (w*(nc*self.psf_base(gc,sc,rlist)+
                           nt*self.psf_base(gt,st,rlist))).sum(axis=1)*sp*weightspectrum(e)
                s_list += list(sc) + list(st)
            else:
                g,s,w = psf.get_p(e,band.ct).copy()
                s *= scale 
                pdf += (w*self.psf_base(g,s,rlist)).sum(axis=1)*sp*weightspectrum(e)
                s_list += list(s)

            weight_sum += sp*weightspectrum(e)
        
        # Normalize
        pdf /= weight_sum 

        # Function returns the difference between the true pdf and 
        # the psf for a given gamma & sigma. Note that the most
        # accurate thing is the probability per unit radius summed over
        # all of the radii.
        f=lambda p: N.std(rlist*(self.psf_base(p[0],p[1],rlist).sum(axis=1)-pdf))

        fit = fmin(f,[gamma_middle,sigma_middle],full_output=False,disp=False)
        fit_gamma,fit_sigma = fit

        # the fit sigma should be somewhere between the two extremes.
        # Not so sure about the gamma tail parameter being a good
        # measure.
        if fit_sigma < min(s_list) or fit_sigma > max(s_list): 
            print dedent("""\
                Error: Sigma fit outside of a 
                reasonable range. fit sigma is %g,
                minimum sigma is %g. maximum sigma
                is %g.""" % \
                (fit_sigma,min(s_list),max(s_list)))

        band.fit_gamma,band.fit_sigma=fit_gamma,fit_sigma


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
                            N.log(
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
            zp          = self.energyBandLikelihoodExtended(N.asarray([-20]),self.m)

            # NB -- the 95% upper limit is calculated by assuming the likelihood is peaked at
            # 0 flux and finding the flux at which it has fallen by 1.35; this is a two-sided
            # 90% limit, or a one-sided 95% limit -- that's how it works, right?
            def f95(parameters):
                return abs(self.energyBandLikelihoodExtended(parameters,self.m) - zp - 1.35)
            
            # for some reason, can't get fsolve to work here.  good ol' fmin to the rescue
            self.energy_band.uflux = 10**fmin(f95,N.asarray([-11.75]),disp=0)[0]
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
