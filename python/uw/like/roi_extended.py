from SpatialModels import RadiallySymmetricModel,Gaussian,SpatialModel
from uw.utilities.minuit import Minuit
from uw.utilities.convolution import BackgroundConvolution,BackgroundConvolutionNorm,AnalyticConvolution
from roi_diffuse import DiffuseSource,ROIDiffuseModel_OTF
from scipy.optimize import fmin
from skymaps import SkyDir,Background
import numpy as N

class ExtendedSource(DiffuseSource):
    """ Class inherting from DiffuseSource but implementing a spatial source. 
        The main difference is the requirement of a spatial model to accomany 
        a spectral model. """

    def __init__(self,name=None,model=None,spatial_model=None,free_parameters=True):
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

        for i in xrange(len(self.smodel.free)): self.smodel.free[i] = free_parameters

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
        self.nsimps    = 0

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

    def localize(self,roi,which,bandfits=True,tolerance=1e-3,update=False,verbose=False,
                 error="HESSE",init_grid=None):
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
  =========   =======================================================

        Implemenation notes, do not directly fit RA and
              Dec. Instead, rotate the coordiante to (0,0) and
              fit longitude and latitude in this rotated coordiante
              system. For each function iteration, rotate back to the
              original space and set the direction to there.  This is a
              more robust fitting algorithm for high latitude and allows
              for errors that are physical distance. 
              
        N.B This code is ugly, mainly because you don't want to directly
        fit the spatial parameters of a source but fit the relative
        deviation of the source away from its starting value. This makes
        the method robust anywhere in the sky. This gets bad because you
        may not be directly fitting the spatail parameters (by fixing the
        spatial parameters), so the code has to intellegently """

        es = self.extended_source
        sm = es.spatial_model
        pn = sm.get_param_names()
        cs = sm.coordsystem

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
                    ll = ll_alt
                else:
                    self.smodel.set_parameters(prev_fit)
            return ll

        def bandlike():
            """ Helper function which takes in the spatial parameters
                for the extended source and returns the logLikelihood.
            """
            ll = 0
            if 'energy_bands' not in roi.__dict__.keys(): roi.setup_energy_bands()

            for eb in roi.energy_bands: 
                eb.bandFit(which=which,saveto='bandfits',diffuse=True)

            for band in roi.bands:

                # Only include likelihood from good bands.
                if band.bandfit < 0: 
                    continue
                ll  += b.bandLikelihood(band.bandfit,which,diffuse=True)
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

            # in rotated coordiante system
            new_dir = SkyDir(lon,lat,cs)

            # This rotation takes (0,0) back to the initial position
            new_dir().rotateY(N.radians(-init_lat)).rotateZ(N.radians(init_lon))

            if cs == SkyDir.GALACTIC:
                if N.any(pn=='b'): p=N.append(new_dir().b(),p)
                if N.any(pn=='l'): p=N.append(new_dir().l(),p)
            elif cs == SkyDir.EQUATORIAL:
                if N.any(pn=='Dec'): p=N.append(new_dir().dec(),p)
                if N.any(pn=='RA'): p=N.append(new_dir().ra(),p)

            sm.set_parameters(p=p,absolute=False)
            self.initialize_counts(roi.bands)

            if not bandfits:
                ll = loglike()
            else:
                ll = bandlike()


            if verbose: print '%s, logL = %.2f, dlogL = %.2f' % (sm.pretty_string(),ll,ll-ll_0)
            return -ll

        f=likelihood_wrapper

        ll_0 = 0
        ll_0 = -f(init_spatial)

        print 'Localizing %s source %s' % (sm.pretty_name,es.name)

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
            print 'Starting with the best initial parameters: ',init_grid[index]
        else:
            start_spatial = init_spatial

        
        # these guys should eventually be prompoted to input parameters

        m = Minuit(f,start_spatial,up=.5,maxcalls=20000,tolerance=tolerance,
                   printMode=verbose,param_names=pn,
                   limits=sm.get_limits(absolute=False))
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

        # return log likelihood from fitting extension.
        return -fval

###=========================================================================###


class ROIExtendedModelAnalytic(ROIExtendedModel):
    """ Implements the ROIDiffuseModel interface for a radially
        symmetric extended source. Utilize teh semi-analytic
        convolution for a more efficient pdf calculation.

        Any of the optional keyword arguments to uw.utilities.convolution's 
        AnalyticConvolution class will be passed on to that class.  """

    def setup(self):
        self.exp = self.sa.exposure.exposure; 
        psf = self.sa.psf
        self.bgc = AnalyticConvolution(self.extended_source.spatial_model,psf,**self.__dict__)

    def set_state(self,energy,conversion_type):
        self.current_energy = energy
        self.current_exposure = self.exp[conversion_type].value(self.extended_source.spatial_model.center,energy)
        self.bgc.do_convolution(energy,conversion_type)

    def _ap_value(self,center,radius):
        return self.current_exposure*self.bgc.ap_average(center,radius)

    def _pix_value(self,pixlist):
        return self.current_exposure*self.bgc(pixlist)


