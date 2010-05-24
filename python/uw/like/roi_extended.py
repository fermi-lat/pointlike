from SpatialModels import RadiallySymmetricModel,Gaussian,SpatialModel
from uw.utilities.convolution import BackgroundConvolution,BackgroundConvolutionNorm,AnalyticConvolution
from roi_diffuse import DiffuseSource,ROIDiffuseModel_OTF
from scipy.optimize import fmin
from skymaps import SkyDir,Background
import numpy as N

class ExtendedSource(DiffuseSource):

    def __init__(self,name=None,model=None,spatial_model=None,free_parameters=True):
        """ Make the naming consistent with the PointSource object so that
            extended sources 'feel' like point sources. 

            spatial_model should inherit from SpatialModel
            
            """
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

    def pretty_string(self):
        """ Return a one line string of the spatial parameters of the extended source. """
        return 'dir = (%.3f,%.3f), ext = %s' % (self.spatial_model.center.ra(),
                                              self.spatial_model.center.dec(),
                                              self.spatial_model.pretty_string())


###=========================================================================###


class ROIExtendedModel(ROIDiffuseModel_OTF):
    """ Implements the ROIDiffuseModel interface for the
        representation of a spatial source which can be
        represented as an analytic function inherting from
        the SpatialModel class."""

    @staticmethod 
    def factory(spectral_analysis,extended_source,*args,**kwargs):

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
            spectral_analysis = spectral_analysis,
            diffuse_source =    extended_source,
            roi_dir =           roi_dir,
            scaling_model =     extended_source.smodel,
            name =              extended_source.name)

    def __str__(self):
        es = self.extended_source

        return '%s fitted with %s\n%s\n%s fitted with %s\n%s' % \
                (es.name,es.smodel.pretty_name,
                 es.smodel.__str__(),
                 ' '*len(es.name),es.spatial_model.pretty_name,
                 es.spatial_model.__str__())

    def localize(self,roi,which):

        es = self.extended_source

        old_dsm_quiet = roi.dsm.quiet
        old_quiet = roi.quiet
        roi.quiet = True

        ll_0 = -roi.logLikelihood(roi.parameters())

        init_vals = es.smodel.p.copy()

        def likelihood_wrapper(params):
            """ Helper function which takes in the spatial parameters for the extended
                source and returns the logLikelihood. 
                
                Implemenation note: Sometimes the fitter gets
                confused. For example, if the last iteration moved to
                a spatial value really far from the initial position,
                the source flux will get fit to 0 (as it should), but
                during the next iteration at a more reasonable spatial
                value, the fit will be unable to reconverge on the best
                flux and stays at 0. To get around this, I cache the
                initial fit value (which is presumably pretty good), and
                whenever the current fit logLikelihood is less then the
                initial logLikelihood, I rest the loglikelihood to
                
                """
            p=params[2:]
            center= SkyDir(*params[0:2])
            es.spatial_model.update(p=p, center=center)
            self.initialize_counts(roi.bands)
            # roi.dsm.update_counts handled by fit() calling logLikelihood()
            roi.fit()

            ll=float(-roi.logLikelihood(roi.parameters()))

            if ll < ll_0:
                prev_fit =self.smodel.p.copy()
                self.smodel.p = init_vals.copy()
                roi.fit()
                ll_alt=float(roi.logLikelihood(roi.parameters()))

                if ll_alt > ll:
                    ll = ll_alt
                else:
                    self.smodel.p = prev_fit

                ll=max(ll_alt,ll)

            if not old_quiet: print '%s, logL = %.2f, dlogL = %.2f' % (es.pretty_string(),ll,ll-ll_0)
            return -ll

        print 'Localizing %s source %s' % (es.spatial_model.pretty_name,es.name)

        init_dir=es.spatial_model.center
        init_params=[init_dir.ra(),init_dir.dec()]+es.spatial_model.p
        f=fmin(likelihood_wrapper,init_params,full_output=1,
               maxiter=10000,maxfun=20000, ftol=1e-2, disp=0 if old_quiet else 1)

        roi.quiet = old_quiet

###=========================================================================###


class ROIExtendedModelAnalytic(ROIExtendedModel):
    """ Implements the ROIDiffuseModel interface for a radially
        symmetric extended source. Utilize teh semi-analytic
        convolution for a more efficient pdf calculation."""

    def setup(self):
        self.exp = self.sa.exposure.exposure; 
        psf = self.sa.psf
        self.bgc = AnalyticConvolution(self.extended_source.spatial_model,psf)

    def set_state(self,energy,conversion_type):
        self.current_energy = energy
        self.current_exposure = self.exp[conversion_type].value(self.extended_source.spatial_model.center,energy)
        self.bgc.do_convolution(energy,conversion_type)

    def _ap_value(self,center,radius):
        solid_angle=2*N.pi*(1-N.cos(radius))
        return self.current_exposure/solid_angle

    def _pix_value(self,pixlist):
        return self.current_exposure*self.bgc(pixlist)


