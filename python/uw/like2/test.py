"""
All like2 testing code goes here, using unittest
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/test.py,v 1.1 2013/11/07 17:21:28 burnett Exp $
"""
import os, sys, unittest
import numpy as np

from uw.like2 import ( configuration, 
    diffusedict as diffuse,
    sources,
    exposure,
    response,
    extended,
    )
import skymaps
from skymaps import SkyDir, Band
# to define the configuration
config = None
ecat = None

class Bandlite(object):
    """ behaves like ROIBand, but does not require data
    Default initialization is for back, first energy bin above 100 MeV
    """
    def __init__(self, roi_dir, config, radius=5, event_type=1, emin=10**2, emax=10**2.25):#energy=10**2.125):
        self.skydir=roi_dir
        self.radius =5
        self.event_type = event_type
        self.emin, self.emax = emin, emax
        self.energy = energy =  np.sqrt(emin*emax)
        self.psf=config.psfman(event_type, energy)
        self.exposure = config.exposureman(event_type, energy)
        
    def set_energy(self, energy):
        self.psf.setEnergy(energy)
        self.exposure.setEnergy(energy)
        self.energy=energy
    @property
    def radius_in_rad(self): return np.radians(self.radius)
    @property #alias, for compatibilty, but deprecated
    def sd(self): return self.skydir
    @property
    def solid_angle(self):
        return np.pi*self.radius_in_rad**2 

def setup_config_dir(skymodelname='P202/uw29'):
    return os.path.join(os.path.expandvars('$HOME/skymodels'),skymodelname)
    
class TestSetup(unittest.TestCase):
    def setUp(self, force=False):
        """Configuration assuming P202_uw29, back 133 MeV"""
        global config
        if config is None:
            self.config_dir=setup_config_dir()
            config =  configuration.Configuration(self.config_dir, quiet=True, postpone=True)
        self.config = config
        self.skydir = Band(12).dir(2)#SkyDir()
        
class TestConfig(TestSetup):
    
    def test_psf(self):
        """check that the PSF is set up
        """
        psf = self.config.psfman(1, 1000)
        
        self.assertDictContainsSubset(dict(event_type=1, energy=1000), psf.__dict__)
        psf.setEnergy(133)
        self.assertEquals(133, psf.energy)
        self.assertAlmostEqual(81.0904920, psf(0)[0])
        
    def test_exposure(self):
        """check exposure"""
        exposure = self.config.exposureman(1, 1000)
        self.assertDictContainsSubset(dict(et=1, energy=1000), exposure.__dict__, 'full dict: %s'%exposure.__dict__)
        exposure.setEnergy(133)
        self.assertEquals(133, exposure.energy)
        self.assertAlmostEqual(2.1604e10 , exposure(self.skydir), delta=1e8)
        
    def test_exposure_integral(self, expect=921.5):
        """test the exposure integral at a point
        """
        emin, e, emax = np.logspace(2, 2.25, 3)
        exp = self.config.exposureman(1, e)
        model = sources.PowerLaw(1e-11,2)
        f1 = exposure.ExposureIntegral(exp,self.skydir, emin, emax)(model)
        f = lambda model : exp.model_integral(self.skydir, model, emin, emax)
        f2 = f(model) #exp.model_integral(self.skydir, model, emin, emax)
        self.assertAlmostEquals(f1,f2, delta=1e-2)
        self.assertAlmostEquals(expect, f1, delta=0.1)
        # need to check value print f2, f(model.gradient), (f(sources.PowerLaw(1.1e-11,2))-f(model))
        
    def test_bandlite(self):
        """ check bandlite"""
        band = Bandlite(self.skydir, self.config)
        self.assertDictContainsSubset(dict(radius=5, event_type=1), band.__dict__, str(band.__dict__))
        
    ## move this to fitting, factor out loading data perhaps
    def xtest_get_bands(self):
        """check loading bands from data"""
        bands = self.config.get_bands(self.skydir, minROI=5)
        self.assertEquals(32, len(bands))
        self.assertDictContainsSubset(dict(event_type=0), bands[0].__dict__, 'first band not event type 0')

class TestDiffuse(TestSetup):
    
    def setUp(self, **kwargs):
        super(TestDiffuse,self).setUp(**kwargs)
        self.band = Bandlite(self.skydir,self.config)
        
        
    def test_isotrop(self):
        """istropic source with constant model"""
        glsrc = sources.GlobalSource(name='isotrop', skydir=None,
            model=sources.Constant(1.0),
            dmodel=diffuse.diffuse_factory(['isotrop_4years_P7_V15_repro_v2_source_%s.txt'%s 
                                            for s in self.config.event_class_names]))
        self.resp =resp= response.IsotropicResponse(glsrc, self.band)
        self.assertAlmostEquals(4766, resp.counts, delta=1) # warning: seems to be 4850 in old version
        self.assertAlmostEquals(199772, resp(resp.roicenter), delta=10)

    def test_cached_galactic(self):
        """cached galactic source"""
        glsrc = sources.GlobalSource(name='ring', skydir=None,
                model=sources.Constant(1.0),
                dmodel = diffuse.diffuse_factory(dict(filename='template_4years_P7_v15_repro_v2_4bpd.zip'))
                )
        self.glsrc=glsrc
        self.resp =resp= response.CachedDiffuseResponse(glsrc, self.band)
        self.assertAlmostEquals(747, resp.ap_average, delta=1)
        self.assertAlmostEquals(1391, resp.counts, delta=1) # expected  1535
        self.assertAlmostEquals(57312, resp(resp.roicenter), delta = 100)
        
        
    def test_map_cube(self):
        """a MapCube source"""
        glsrc = sources.GlobalSource(name='ring1', skydir=None,
                model=sources.Constant(1.0),
                dmodel = diffuse.diffuse_factory(dict(filename='template_4years_P7_v15_repro_v2.fits'))
                )
        self.glsrc=glsrc
        self.resp =resp= response.DiffuseResponse(glsrc, self.band)
        self.assertAlmostEquals(742, resp.ap_average, delta=1)
        self.assertAlmostEquals(1381, resp.counts, delta=1) # expected  1535
        self.assertAlmostEquals(56876, resp(resp.roicenter), delta = 100)

    def test_healpix_diffuse(self):
        """a MapCube source"""
        glsrc = sources.GlobalSource(name='ring2', skydir=None,
                model=sources.Constant(1.0),
                dmodel = diffuse.diffuse_factory(dict(filename='template_4years_P7_v15_repro_v2_nside256_bpd4.fits',
                    type='Healpix',  ))
                )
        self.glsrc=glsrc
        self.resp =resp= response.DiffuseResponse(glsrc, self.band)
        self.assertAlmostEquals(722, resp.ap_average, delta=1)
        self.assertAlmostEquals(1344, resp.counts, delta=1) # expected  1535
        self.assertAlmostEquals(55354, resp(resp.roicenter), delta = 100)
        
    
 
    
class TestPoint(TestSetup):
    def setUp(self, **kwargs):
        super(TestPoint,self).setUp(**kwargs)
        self.band = Bandlite(self.skydir,self.config)
  
    
    def test_point(self):
        """point source at the center"""
        ptsrc =sources.PointSource(name='test', skydir=self.skydir, 
                                           model=sources.PowerLaw(1e-11, 2.0))
        self.resp =resp = response.PointResponse(ptsrc, self.band)
        self.assertAlmostEquals(0.65, resp.overlap, delta=0.01)
        self.assertAlmostEquals(1.0, resp._exposure_ratio)
        self.assertAlmostEquals(600, resp.counts, delta=1.)
        self.assertAlmostEquals(75042, resp(resp.source.skydir), delta=10)
        
    def make_test_source(self, offset, expected_overlap):
        model = sources.PowerLaw(1e-11,2.0)
        source = sources.PointSource(name='test source', 
            skydir=SkyDir(self.skydir.ra(), self.skydir.dec()+offset), model=model)
        conv = source.response(self.band) 
        overlap = conv.overlap
        self.assertAlmostEqual(expected_overlap, overlap, delta=0.001)
        #a,b = conv.evaluate_at([source.skydir, self.band.sd])/1e12
        #c = conv(source.skydir)/1e12
        #self.assertAlmostEqual(a,c)
        
    def test_create_on_center(self):
        self.make_test_source(0, 0.650)
    def test_create_2deg(self):
        self.make_test_source(2, 0.604)
    def test_create_6deg(self):
        self.make_test_source(4, 0.450)
   
class TestExtended(TestSetup):

    def setUp(self):
        super(TestExtended, self).setUp()
        global ecat
        if ecat is None:
            ecat = extended.ExtendedCatalog(self.config.extended)

    def extended(self, source_name='W28', roi=None, expect=(0,), psf_check=True, model=None, quiet=True):
        source = ecat.lookup(source_name)
        if model is not None:
            source.model = model
        self.assertIsNotNone(source, 'Source %s not in catalog'%source_name)
        b12 = skymaps.Band(12);
        roi_index=roi if roi is not None else b12.index(source.skydir)
        roi_dir = b12.dir(roi_index) 
        difference = np.degrees(roi_dir.difference(source.skydir))
        band = Bandlite(roi_dir, self.config)
        if not quiet:
            print 'Using ROI #%d, distance=%.2f deg' %( roi_index, difference)
            print 'Testing source "%s at %s" with band parameters' % (source, source.skydir)
            for item in band.__dict__.items():
                print '\t%-10s %s' % item
        self.resp = conv = source.response(band)
        if not quiet:
            print 'overlap: %.3f,  exposure_ratio: %.3f' %( conv.overlap,conv.exposure_ratio)
            print 'PSF overlap: %.3f'% conv.psf_overlap
        self.assertAlmostEqual(expect[0], conv.overlap, delta=1e-2)
        if psf_check:
            self.assertAlmostEqual(conv.overlap, conv.psf_overlap, delta=1e-2)
        if len(expect)>1:
            self.assertAlmostEqual(expect[1], conv.counts, delta=1)
        
    def test_W28(self):
        self.extended('W28', expect=(0.590,3141), 
            model=sources.LogParabola(3.38e-11, 2.27, 0.127, 1370))
    def test_W30_in840(self):
        self.extended('W30', 840, expect=(0.377,1038,), 
            model=sources.LogParabola(1.31e-11,2.15, 0.036, 1430) )
    def test_LMC(self):
        self.extended('LMC', expect=(0.593,), psf_check=False)

    def test_Cygnus_Cocoon(self):
        self.extended('Cygnus Cocoon', expect=(0.551,), psf_check=False)


        
def test_suite(t):
    return unittest.TestLoader().loadTestsFromTestCase(t)
    
def run(t='all'): 
    if t=='all':
        suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    else:
        suite = unittest.TestLoader().loadTestsFromTestCase(t)
    print 'running %d tests' % suite.countTestCases()
    unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)