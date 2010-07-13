""" Provides classes to encapsulate and manipulate extended sources. 

    This code all derives from objects in roi_diffuse.py

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_extended.py,v 1.9 2010/07/06 23:01:05 lande Exp $

    author: Joshua Lande
"""

from SpatialModels import RadiallySymmetricModel,Gaussian,SpatialModel
from uw.utilities.convolution import BackgroundConvolution,BackgroundConvolutionNorm,AnalyticConvolution
from roi_diffuse import DiffuseSource,ROIDiffuseModel_OTF
from textwrap import dedent
from skymaps import SkyDir,Background
from scipy.optimize import fmin 
import numpy as N
from uw.like.Models import PowerLaw

class ExtendedSource(DiffuseSource):
    """ Class inherting from DiffuseSource but implementing a spatial source. 
        The main difference is the requirement of a spatial model to accomany 
        a spectral model. """

    def __init__(self,name=None,model=None,spatial_model=None,free_parameters=True,leave_parameters=False):
        """ Make the naming consistent with the PointSource object so that
            extended sources 'feel' like point sources.

            spatial_model should inherit from SpatialModel. """
        if spatial_model is None: spatial_model = Gaussian()

        self.spatial_model = spatial_model

        if model is None: model = PowerLaw()

        if not isinstance(self.spatial_model,SpatialModel):
            raise Exception("The diffuse_model passed to an Extended Source must inherit from SpatialModel.")

        super(ExtendedSource,self).__init__(
            diffuse_model = spatial_model.get_PySkySpectrum(),
            scaling_model = model,
            name          = name)

        self.smodel.background = False

        if not leave_parameters:
            for i in xrange(len(self.smodel.free)): self.smodel.free[i] = free_parameters
            for i in xrange(len(self.spatial_model.free)): self.spatial_model.free[i] = free_parameters

    def __str__(self):
        return '\n'.join(['\n',
                          ''.join(['=']*60),
                          'Name:\t\t%s'%(self.name),
                          'R.A. (J2000):\t\t%.5f'%(self.spatial_model.center.ra()),
                          'Dec. (J2000):\t\t%.5f'%(self.spatial_model.center.dec()),
                          'Model:\t\t%s'%(self.smodel.name),
                          'SpatialModel:\t%s'%(self.spatial_model.pretty_name)])


###=========================================================================###


class ROIExtendedModel(ROIDiffuseModel_OTF):
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

        if not isinstance(extended_source,ExtendedSource):
            raise Exception("The extended_source option passed to ROIExtendedModel.factory() must inherit from ExtendedSource.")

        spatial_model = extended_source.spatial_model

        if isinstance(spatial_model,RadiallySymmetricModel):
            return ROIExtendedModelAnalytic(spectral_analysis,extended_source,*args,**kwargs)
        elif isinstance(spatial_model,SpatialModel):
            return ROIExtendedModel(spectral_analysis,extended_source,*args,**kwargs)
        else:
            raise Exception("The extended_source.dmodel option passed to ROIExtendedModel.factory must inherit from SpatialModel.")

    def init(self,*args,**kwargs):
        self.pixelsize = 0.05
        self.npix      = 512
        self.nsimps    = 16 # for consistency with point sources.

    def setup(self):
        """ Use the Normalized convolution object and always do the
            convolution around the spatial model's center. """
        exp = self.sa.exposure.exposure; psf = self.sa.psf
        self.bg  = [Background(self.dmodel[0],exp[0],exp[1])]
        self.bgc = [BackgroundConvolution(self.extended_source.spatial_model.center,self.bg[0],psf,
                    npix=self.npix,pixelsize=self.pixelsize)]


    def __init__(self,spectral_analysis,extended_source,roi_dir,name=None,*args,**kwargs):

        self.extended_source = self.diffuse_source = extended_source

        super(ROIExtendedModel,self).__init__(
            spectral_analysis, extended_source,
            roi_dir, extended_source.name,*args,**kwargs)

    def __str__(self):
        es = self.extended_source
        sm = es.spatial_model

        return '%s fitted with %s\n%s\n%s fitted with %s\n%s' % \
                (es.name,sm.pretty_name,
                 sm.__str__(),
                 es.name,es.smodel.pretty_name,
                 es.smodel.__str__())

    def localize(self,roi,which,tolerance=1e-3,update=False,verbose=False,
                 bandfits=True, error="HESSE",init_grid=None,fitpsf=False):
        """ Localize this extended source by fitting all non-fixed spatial paraameters of 
            self.extended_source. The likelihood at the best position is returned.

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  bandfits      [True]  Spectrum independent fitting.
  tolerance     [1e-3]  Tolerance set when integrating hypergeometirc
                        function to calcualte PDF.
  update        [False] Number of points to calculate the PDF at.
                        interpolation is done in between.
  verbose       [False] Make noise when fitting

  error         [HESSE] The fitting algorithm to use when calculating errors.

  init_grid     [None]  A list of spatial parameters. The likelihood
                        for each of the spatial parametesr is tested and the minimum
                        one is used as the initial guess when running minuit. Useful
                        for doing initial grid search. 

                        - The spatial part of init_grid
                          should be measured relative to the spatial model's center: (0,0)
                        - The init_grid values should be absolute, none of them should
                          be given in log space.

  fitpsf        [False] Use an approximation to the PSF in each energy bin
                        where PSF is weighted in each energy bin by an assumed
                        spectral index and is fit by a single king function
                        which is then convovled with the extended source shape.
                        This leads to a significantly faster source localization.
  =========   =======================================================

        Implemenation notes, do not directly fit RA and
              Dec. Instead, rotate the coordiante to (0,0) and
              fit x & y deviation (measured in degrees) in this 
              rotated coordinate system. For each function iteration, rotate back to the
              original space and set the direction to there.  This is a
              more robust fitting algorithm for high latitude and allows
              for errors that are physical distance. 
              
        N.B This code is ugly, mainly because you don't want to directly
        fit the spatial parameters of a source but fit the relative
        deviation of the source away from its starting value. This makes
        the method robust anywhere in the sky. This gets bad because you
        may not be directly fitting the spatail parameters (by fixing the
        spatial parameters), so the code has to intellegently """

        from uw.utilities.minuit import Minuit

        if fitpsf:
            if verbose: print 'Changing to fitpsf accuracy for localization step.'
            self.nsimps,old_nsimps=0,self.nsimps
            self.fitpsf,old_fitpsf=True,self.fitpsf

        es = self.extended_source
        sm = es.spatial_model
        pn = sm.get_param_names(absolute=True,all=False)
        cs = sm.coordsystem

        origin=SkyDir(0,0,cs)

        old_dsm_quiet = roi.dsm.quiet
        old_quiet = roi.quiet
        roi.quiet = True

        init_spectral = self.smodel.get_parameters()
        init_spatial = sm.get_parameters(absolute=False)

        if len(init_spatial) < 1:
            print 'Unable to localize diffuse source %s. No parameters to fit.' % es.name
            return

        # Remember thet get_parameters only returns the free parameters.
        # here we need to get the first two parameters (position) so
        # we can fit relative to it.
        init_lon,init_lat = sm.p[0:2]
        init_dir          = SkyDir(init_lon,init_lat,cs)

        # Fit in coordinate system rotated so y=z=0.
        if cs == SkyDir.GALACTIC:
            init_spatial[(pn=='l')|(pn=='b')] = 0 
        elif cs == SkyDir.EQUATORIAL:
            init_spatial[(pn=='RA')|(pn=='Dec')] = 0 

        def loglike():
            """ Wrapper to call the standard likelihood function and,
                if it fails to converge, retry the likelihood from 
                the initial spectral values. This is useful in case
                the fitting algorithm moves way too far away and
                the source flux gets fit to 0, the fitter may be
                too far away from the minimum to converge. """
            ll=roi.fit(estimate_errors=False)

            if ll < ll_0:
                prev_fit=self.smodel.get_parameters()
                self.smodel.set_parameters(init_spectral)
                ll_alt=roi.fit(estimate_errors=False)

                if ll_alt > ll:
                    ll=ll_alt
                else:
                    self.smodel.set_parameters(prev_fit)
            return ll

        def bandlike():
            """ Helper function which takes in the spatial parameters
                for the extended source and returns the logLikelihood.
            """

            roi.bgm.update_counts(roi.bands)
            roi.psm.update_counts(roi.bands)

            ll = 0
            if 'energy_bands' not in roi.__dict__.keys(): roi.setup_energy_bands()

            # Note, negative sign to be consistent with value returned
            # by roi.fit()
            for eb in roi.energy_bands: 
                ll -= eb.bandFitDiffuse(which=which)

            return ll

        def likelihood_wrapper(p):
            """ Helper function which takes in the spatial parameters
                for the extended source and returns the logLikelihood.

                Implemenation note: Sometimes the fitter gets
                confused. For example, if the last iteration moved to
                a spatial value really far from the initial position,
                the source flux will get fit to 0 (as it should), but
                during the next iteration at a more reasonable spatial
                value, the fit will be unable to reconverge on the best
                flux and stays at 0. To get around this, I cache the
                initial fit value (which is presumably pretty good),
                and whenever the current fit logLikelihood is less then
                the initial logLikelihood, I rest the loglikelihood to 
                
                N.B. roi.dsm.update_counts handled by fit() calling logLikelihood()
                """

            # Only pull ra & dec/l & b out of parameters if they are being fit.
            # otherwise they are, by definition, equal to 0 (in rotated 
            # coordinates)
            if cs == SkyDir.GALACTIC:
                lon     = float(p[pn=='l'] or 0)
                lat     = float(p[pn=='b'] or 0)
                p=p[(pn!='l')&(pn!='b')]
            elif cs == SkyDir.EQUATORIAL:
                lon     = float(p[pn=='RA'] or 0)
                lat     = float(p[pn=='Dec'] or 0)
                p=p[(pn!='RA')&(pn!='Dec')]

            # New direction in rotated coordiante system
            new_dir = SkyDir(lon,lat,cs)

            # find vector perpendicular to both origin & initial dir,
            # rotate around that vector by the distance from the initial
            # angle to the origin. This rotation takes the origin to the
            # initial position and small deviations around the origin to
            # corresponding small deviations.
            cross=origin.cross(init_dir)()
            theta=origin.difference(init_dir)

            # This rotation takes (0,0) back to the initial position
            new_dir().rotate(cross,theta)

            # Now add the rotated spatial part back to the list.
            if cs == SkyDir.GALACTIC:
                if N.any(pn=='b'): p=N.append(new_dir.b(),p)
                if N.any(pn=='l'): p=N.append(new_dir.l(),p)
            elif cs == SkyDir.EQUATORIAL:
                if N.any(pn=='Dec'): p=N.append(new_dir.dec(),p)
                if N.any(pn=='RA'): p=N.append(new_dir.ra(),p)

            # Do the convolution here.
            sm.set_parameters(p=p,absolute=False)
            self.initialize_counts(roi.bands)

            ll = bandlike() if bandfits else loglike()

            if verbose: print '%s, logL = %.2f, dlogL = %.2f' % (sm.pretty_string(),ll,ll-ll_0)
            return -ll

        f=likelihood_wrapper

        ll_0 = 0
        old_verbose = verbose; verbose = False; ll_0 = -f(init_spatial); verbose = old_verbose

        if verbose: print 'Localizing %s source %s Using %s' % (sm.pretty_name,es.name,
                                                    'BandFits' if bandfits else 'Spectral Fits')

        if init_grid is not None:
            print 'Testing initial grid values'
            ll = [] 
            transformed = []
            for _ in init_grid:
                # easy way to turn parameters into their non-absolute value.
                sm.set_parameters(_,absolute=True)
                param=sm.get_parameters(absolute=False)
                ll.append(likelihood_wrapper(param))
                transformed.append(param)

            index=N.argmin(ll)
            start_spatial=transformed[index]
            if verbose: print 'Now starting with the best initial parameters'
        else:
            start_spatial = init_spatial

        # Display which parameters are in log space.
        relative_names = sm.get_param_names(absolute=False,all=False)

        m = Minuit(f,start_spatial,up=.5,maxcalls=20000,tolerance=tolerance,
                   printMode=verbose,param_names=relative_names,
                   limits=sm.get_limits(absolute=False),
                   steps  =sm.get_steps())
        best_spatial,fval = m.minimize(method="SIMPLEX")

        if verbose: print 'Calculating Covariance Matrix'

        cov_matrix = m.errors(method=error)

        sm.set_cov_matrix(cov_matrix)

        if update:
            if verbose: print 'setting source to best fit parameters'
            likelihood_wrapper(best_spatial)
        else:
            if verbose: print 'setting source back to original parameters'
            likelihood_wrapper(init_spatial)

        roi.quiet = old_quiet

        if fitpsf:
            if verbose: print 'Setting back to original PDF accuracy.'
            self.nsimps=old_nsimps
            self.fitpsf=old_fitpsf
            self.initialize_counts(roi.bands)

        # return log likelihood from fitting extension.
        return -fval

###=========================================================================###


class ROIExtendedModelAnalytic(ROIExtendedModel):
    """ Implements the ROIDiffuseModel interface for a radially
        symmetric extended source. Utilize teh semi-analytic
        convolution for a more efficient pdf calculation.

        Any of the optional keyword arguments to uw.utilities.convolution's 
        AnalyticConvolution class will be passed on to that class.  """

    def init(self):
        super(ROIExtendedModelAnalytic,self).init()
        self.fitpsf=False
        self.already_fit=False


    def setup(self):
        self.exp = self.sa.exposure.exposure; 
        psf = self.sa.psf

        if self.fitpsf and self.nsimps > 0:
            raise Exception("ROIExtendedModelAnalytic Error: fitpsf and nsimps>0 are incompatable.")

        self.bgc = AnalyticConvolution(self.extended_source.spatial_model,psf,**self.__dict__)

    def set_state(self,energy,conversion_type,band):
        self.current_energy = energy
        self.current_exposure = self.exp[conversion_type].value(self.extended_source.spatial_model.center,energy)
        self.bgc.do_convolution(energy,conversion_type,band)

    def _ap_value(self,center,radius):
        return self.current_exposure*self.bgc.ap_average(center,radius)

    def _pix_value(self,pixlist):
        return self.current_exposure*self.bgc(pixlist)

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


class BandFitter(object):
    """ This class has a somewhat wierd function. Basically the PSF is
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

        scale = psf.scale_func[band.ct](band.e)
        p=psf.get_p(band.e,band.ct)
        gamma_middle=sum(p[0]*p[2])/sum(p[2])
        sigma_middle=sum(p[1]*p[2])/sum(p[2])
        sigma_middle=sigma_middle*scale

        rlist = N.linspace(0,20*sigma_middle,10000)
        pdf = N.zeros_like(rlist)
        pdf_weight = N.zeros_like(rlist)
        if newstyle:
            raise Exception("PSF Fitting is not yet implemented for the newstyle PSF. Bug Josh to add this")
        else:
            weight_sum,s_list= 0,[]

            # Use the same simpson integral used to calculate
            # predicted counts for a point source.
            for e,sp in zip(band.sp_points,band.sp_vector):
                g,s,w = psf.get_p(e,band.ct).copy()
                scale = psf.scale_func[band.ct](e)
                s *= scale 
                pdf += (w*self.psf_base(g,s,rlist)).sum(axis=1)*sp*weightspectrum(e)
                s_list += list(s)

                weight_sum += sp*weightspectrum(e)
            
            # Normalize
            pdf /= weight_sum 

            # Function returns the difference between the true pdf and 
            # the psf for a given gamma & sigma.
            f=lambda p: N.std(self.psf_base(p[0],p[1],rlist).sum(axis=1)-pdf)

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

