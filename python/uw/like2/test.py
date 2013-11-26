"""
All like2 testing code goes here, using unittest
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/test.py,v 1.18 2013/11/25 23:48:43 burnett Exp $
"""
import os, sys, unittest
import numpy as np
import skymaps
from skymaps import SkyDir, Band

from uw.like2 import ( configuration, 
    diffusedict as diffuse,
    sources,
    bands,
    exposure,
    extended,
    roimodel,
    dataset,
    bandlike,
    views,
    sedfuns,
    associate,
    main,
    )

# globals: references set by setUp methods in classes as needed
config_dir = os.path.expandvars('$HOME/skymodels/P202/uw29')
config = None
ecat = None
roi_index = 840
roi_sources = None
roi_bands = None
blike = None
likeviews = None
user = None

def setup(name):
    global config, ecat, roi_sources, roi_bands, blike, likeviews, user
    gnames = 'config ecat roi_sources roi_bands blike likeviews user'.split()
    assert name in globals() and name in gnames
    if config is None :
        config =  configuration.Configuration(config_dir, quiet=True, postpone=True)
    if ecat is None:
        ecat = extended.ExtendedCatalog(config.extended)
    if (name=='roi_sources' or name=='blike' or name=='likeviews') and roi_sources is None:
        roi_sources = roimodel.ROImodel(config, ecat=ecat, roi_index=roi_index)
    if (name=='roi_bands' or name=='blike' or name=='likeviews') and roi_bands is None:
        roi_bands = bands.BandSet(config, roi_index)
    if name=='blike' and blike is None:
        assert roi_bands is not None and roi_sources is not None
        roi_bands.load_data()
        blike = bandlike.BandLikeList(roi_bands, roi_sources)
    if name=='likeviews' and likeviews is None:
        assert roi_bands is not None and roi_sources is not None
        if roi_bands.pixels==0:
            roi_bands.load_data()
        likeviews = views.LikelihoodViews(roi_bands, roi_sources)
    if name=='user' and user is None:
        user = main.ROI_user(config_dir, roi_index)
    return eval(name)
        
class TestSetup(unittest.TestCase):
    def setUp(self, force=False):
        """Configuration assuming P202_uw29, back 133 MeV"""
        self.config = setup('config')
        # use ROI 2 for some simple tests
        self.skydir = Band(12).dir(2)#SkyDir()
        
class TestConfig(TestSetup):
    """Test aspects of the configuration    
    """
    def test_psf(self):
        """check that the PSF is set up
        """
        psf = self.config.psfman(1, 1000)
        self.assertDictContainsSubset(dict(event_type=1, energy=1000), psf.__dict__)
        psf.setEnergy(133)
        self.assertEquals(133, psf.energy)
        self.assertAlmostEqual(81.0904920, psf(0)[0], msg='value expected with IRF for pass 7')
        
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
        band = bands.Bandlite(self.skydir, self.config)
        self.assertDictContainsSubset(dict(radius=5, event_type=1), band.__dict__, str(band.__dict__))

class TestDiffuse(TestSetup):
    
    def setUp(self, **kwargs):
        super(TestDiffuse,self).setUp(**kwargs)
        self.back_band = bands.Bandlite(self.skydir,self.config, event_type=1)
        self.front_band = bands.Bandlite(self.skydir,self.config, event_type=0)
        
    def test_factory(self):
        """test the factory setup """
        for t in ['junk.txt', ('junk.txt','t'),'tst_PowerLaw(1e-11, c )', 
                      dict(file='template_4years_P7_v15_repro_v2_4bpd.zip'),
                  ]:
            self.assertRaises(diffuse.DiffuseException, diffuse.diffuse_factory, t )
            
        test_values = [ 
                ('isotrop_4years_P7_V15_repro_v2_source_front.txt', 
                    'isotrop_4years_P7_V15_repro_v2_source_back.txt'),
                dict(filename='template_4years_P7_v15_repro_v2_4bpd.zip',
                        correction='../../P202/uw29/galactic_correction_uw26a_v2.csv', 
                        systematic=0.0316),
                'template_4years_P7_v15_repro_v2.fits',
                'limb_PowerLaw(1e-11, 4.0)',
                ]
        for t in test_values:
            diffuse.diffuse_factory(t)

        ## test that evaluting returns the same object id
        ids = map(lambda f: id(diffuse.diffuse_factory(f)[0]), [test_values[i] for i in (0,1,0)])
                
        self.assertEquals(ids[0],ids[2], msg='expect the same object')
        self.assertNotEquals(ids[0],ids[1], msg='expect different objects')
        
        
    def test_isotrop(self):
        """istropic source with constant model"""
        source = sources.GlobalSource(name='isotrop', skydir=None,
            model=sources.Constant(1.0),
            dmodel=diffuse.diffuse_factory(['isotrop_4years_P7_V15_repro_v2_source_%s.txt'%s 
                                            for s in self.config.event_class_names]))
        self.resp =resp=source.response(self.back_band)
        self.assertAlmostEquals(4626, resp.counts, delta=1) # warning: seems to be 4850 in old version
        self.assertAlmostEquals(193864, resp(resp.roicenter), delta=10)

    def test_cached_galactic(self):
        """cached galactic source"""
        source = sources.GlobalSource(name='ring', skydir=None,
                model=sources.Constant(1.0),
                dmodel = diffuse.diffuse_factory(dict(filename='template_4years_P7_v15_repro_v2_4bpd.zip'))
                )
        resp= source.response(self.back_band)
        self.response_check(resp, (747, 1391, 57312))
        
    def load_map_cube(self):
        return sources.GlobalSource(name='ring1', skydir=None,
                model=sources.Constant(1.0),
                dmodel = diffuse.diffuse_factory(dict(filename='template_4years_P7_v15_repro_v2.fits'))
                )
        
    def test_map_cube(self):
        """a MapCube source"""
        source = self.load_map_cube()
        resp= source.response( self.back_band)
        self.response_check(resp, (742, 1381, 56876))

    def test_map_cube_front(self):
        """a MapCube source"""
        source = self.load_map_cube()
        resp= source.response( self.front_band)
        self.response_check(resp, (618, 1151, 46894))

    def load_healpix(self):
        source = sources.GlobalSource(name='ring2', skydir=None,
                model=sources.Constant(1.0),
                dmodel = diffuse.diffuse_factory(dict(
                    filename='template_4years_P7_v15_repro_v2_nside256_bpd4.fits',
                    type='Healpix',  ))
                )
        return source

    def response_check(self, resp, expect):
        self.assertAlmostEquals(expect[0], resp.ap_average, delta=1)
        if len(expect)==1: return
        self.assertAlmostEquals(expect[1], resp.counts, delta=1) 
        if len(expect)==2: return
        self.assertAlmostEquals(expect[2], resp(resp.roicenter), delta = 100)

    def test_healpix(self):
        """a Healpix source"""
        source = self.load_healpix()
        self.resp =resp= source.response(self.back_band)
        self.response_check(resp, (722, 1344, 55354))

    def test_healpix_front(self):
        """a Healpix source"""
        source = self.load_healpix()
        self.resp =resp= source.response(self.front_band)
        self.response_check(resp, (602, 1120,  45657))
        
    def test_limb(self):
        """The PowerLaw limb -- not understood yet"""
        source = sources.GlobalSource(name='limb', skydir=None,
            model = sources.FBconstant(2.0, 1.0),
            dmodel=diffuse.diffuse_factory('limb_PowerLaw(1e-11, 4.0)'))
        self.resp_back = source.response(self.back_band)
        self.assertAlmostEquals(2545, self.resp_back.counts, delta=10)
        self.resp_front = source.response(self.front_band)
        self.assertAlmostEquals(2136, self.resp_front.counts, delta=10)
    
class TestPoint(TestSetup):
    def setUp(self, **kwargs):
        super(TestPoint,self).setUp(**kwargs)
        self.back_band = bands.Bandlite(self.skydir,self.config)
  
    
    def test_point(self):
        """response of point source at the center"""
        ptsrc =sources.PointSource(name='test', skydir=self.skydir, 
                                           model=sources.PowerLaw(1e-11, 2.0))
        self.resp =resp = ptsrc.response(self.back_band)
        self.assertAlmostEquals(0.65, resp.overlap, delta=0.01)
        self.assertAlmostEquals(1.0, resp._exposure_ratio)
        self.assertAlmostEquals(600, resp.counts, delta=1.)
        self.assertAlmostEquals(75042, resp(resp.source.skydir), delta=10)
        
    def make_test_source(self, offset, expected_overlap):
        model = sources.PowerLaw(1e-11,2.0)
        source = sources.PointSource(name='test source', 
            skydir=SkyDir(self.skydir.ra(), self.skydir.dec()+offset), model=model)
        conv = source.response(self.back_band) 
        overlap = conv.overlap
        self.assertAlmostEqual(expected_overlap, overlap, delta=0.001)
        #a,b = conv.evaluate_at([source.skydir, self.back_band.sd])/1e12
        #c = conv(source.skydir)/1e12
        #self.assertAlmostEqual(a,c)
        
    def test_create_2deg(self):
        self.make_test_source(2, 0.604)
    def test_create_6deg(self):
        self.make_test_source(4, 0.450)
   
class TestExtended(TestSetup):

    def setUp(self):
        super(TestExtended, self).setUp()

    def extended(self, source_name='W28', roi=None, expect=(0,), psf_check=True, model=None, quiet=True):
        source = ecat.lookup(source_name)
        if model is not None:
            source.model = model
        self.assertIsNotNone(source, 'Source %s not in catalog'%source_name)
        b12 = skymaps.Band(12);
        roi_index=roi if roi is not None else b12.index(source.skydir)
        roi_dir = b12.dir(roi_index) 
        difference = np.degrees(roi_dir.difference(source.skydir))
        band1 = bands.Bandlite(roi_dir, self.config)
        if not quiet:
            print 'Using ROI #%d, distance=%.2f deg' %( roi_index, difference)
            print 'Testing source "%s at %s" with band parameters' % (source, source.skydir)
            for item in band1.__dict__.items():
                print '\t%-10s %s' % item
        self.resp = conv = source.response(band1)
        if not quiet:
            print 'overlap: %.3f,  exposure_ratio: %.3f' %( conv.overlap,conv.exposure_ratio)
            print 'PSF overlap: %.3f'% conv.psf_overlap
        if psf_check:
            self.assertAlmostEqual(conv.overlap, conv.psf_overlap, delta=1e-2)
        self.assertAlmostEqual(expect[0], conv.overlap, delta=1e-2)
        if len(expect)>1:
            self.assertAlmostEqual(expect[1], conv.counts, delta=10)
        
    def test_W28(self):
        self.extended('W28', expect=(0.65,3483), 
            model=sources.LogParabola(3.38e-11, 2.27, 0.127, 1370))
    def test_W30_in840(self):
        self.extended('W30', 840, expect=(0.377,1038,), 
            model=sources.LogParabola(1.31e-11,2.15, 0.036, 1430) )
    def test_LMC(self):
        self.extended('LMC', expect=(0.593,), psf_check=False)

    def test_Cygnus_Cocoon(self):
        self.extended('Cygnus Cocoon', expect=(0.551,), psf_check=False)

class TestROImodel(TestSetup):

    def setUp(self):
        setup('roi_sources')
        self.pars = roi_sources.parameters.get_parameters()
    def tearDown(self):
        roi_sources.parameters.set_parameters(self.pars)

    def test_properties(self):
        rs = roi_sources
        self.assertEquals(14, sum(rs.free))
        self.assertEquals(82, len(rs.free))
        self.assertEquals(82, len(rs.source_names))
        self.assertEquals(35, len(rs.parameter_names))
        self.assertEquals(35, len(rs.bounds))
        
    def test_source_access(self):
        rs = roi_sources
        #finding a source
        self.assertRaises(roimodel.ROImodelException, rs.find_source, 'junk')
        self.assertEquals('PSR J1801-2451', rs.find_source('PSR*').name)
        self.assertEquals('P7R42735',rs.find_source('*2735').name)
        self.assertEquals('P7R42735',rs.find_source(None).name)
        # adding and removing sources
        self.assertRaises(roimodel.ROImodelException, rs.add_source, rs[0])
        s = rs[0].copy()
        s.name = 'ring2'
        rs.add_source(s)
        self.assertEquals('ring2', rs[-1].name)
        self.assertEquals('ring2', rs.find_source('ring2').name)
        rs.del_source('ring2')
        self.assertRaises(roimodel.ROImodelException, rs.find_source, 'ring2')
    
    def test_parameters(self):
        rs = roi_sources
        k = 3
        parz = rs.parameters.get_parameters()
        self.assertEquals(0, sum(rs.parameters.dirty))
        rs.parameters[k]=0.3
        self.assertEquals(1, sum(rs.parameters.dirty))
        self.assertEquals(0.3, rs.parameters[k])
        rs.parameters.set_parameters(parz)
        
        # subset tests
        ps = rs.parsubset()
        pars =ps.get_parameters()
        pz = pars.copy()
        pz[2]+=0.1
        ps.set_parameters(pz)
        self.assertTrue( np.all( pz== ps.get_parameters()))
        ps.select(8)
        self.assertEquals(pz[8], ps[0])
        # should check other features ...

        
    
class TestBands(TestSetup):

    def setUp(self):
        setup('roi_bands')
            
    def test_load_data(self):
        print roi_bands
        roi_bands.load_data()
        print roi_bands
        self.assertEquals( 88484, roi_bands.pixels)

        
class TestLikelihood(TestSetup):
    def setUp(self):
        self.bl = setup('blike')
        self.init = blike.log_like()
        print 'initial loglike: %.1f ...' % self.init ,
    def tearDown(self):
        self.assertAlmostEquals(self.init, blike.log_like(), 1)
        
    def test_unweight(self):
        self.assertAlmostEquals(0.092, self.bl[0].make_unweight().round(3))
        self.assertAlmostEquals(0.033, self.bl[1].make_unweight().round(3))
    def test_weights(self):
        'check weights for band 1'
        b1 = self.bl[1]
        weights = b1.data / b1.model_pixels
        self.assertAlmostEquals(0.98, weights.mean(), delta=0.01)
        self.assertAlmostEquals(0.032, weights.std(), delta=0.002)
    def test_hessian(self):
        bl = self.bl
        hess = bl.hessian()
        self.assertTrue( np.all(hess.diagonal()>0), msg='diagonal: %s' % hess.diagonal())
        s = np.sqrt(hess.diagonal())
        corr = hess / np.outer(s,s)
        t = np.array(corr.T - corr).flatten()
        self.assertTrue( np.abs(t).max()<0.02)

    def test_bandsubset(self):
        bl = self.bl
        bl.selected = bl
        self.assertTrue(bl.selected == bl)
        bl.selected = bl[1]
        bl.selected = bl[:4]
        parta = bl.log_like()
        bl.selected = bl[4:]
        partb = bl.log_like()
        bl.selected= bl
        total = bl.log_like()
        self.assertAlmostEquals(total, parta+partb)
        
class TestFitterView(TestSetup):
    def setUp(self):
        setup('likeviews')
        
    def test_ts(self, sourcename='P7R42722', expect=800):
        """--> set up a subset fitter view, use it to checka TS value"""
        with likeviews.fitter_view(sourcename) as fv:
            self.assertAlmostEquals(expect, fv.ts(), delta=1)
            
    def test_fitting(self):
        """-->generate a fitter_view, use it to maximize the likelihood"""
        with likeviews.fitter_view() as t:
            a = t()
            self.assertEquals(t.log_like(), -a)
            b, g, sig = t.maximize()
            self.assertAlmostEquals(-697607., a, delta=1)
            self.assertAlmostEquals(-697645., b, delta=1)
        
class TestSED(TestSetup):
    def setUp(self):
        setup('likeviews')
        self.init = likeviews.log_like()
        print 'initial loglike: %.1f ...' % self.init ,
    def tearDown(self):
        self.assertAlmostEquals(self.init, likeviews.log_like())

    def test_sourceflux(self, sourcename='W28', checks=(61.695, 63.816, 4889, 5025)):
        """-->create and check the SED object"""
        with sedfuns.SED(likeviews, sourcename) as sf:
            poiss = sf.full_poiss
            errors =poiss.errors
            pp = sf.all_poiss()
            bandts = np.array([x.ts for x in pp]).sum()
            print 'errors, TS, bandts: %.3f, %.3f %.3f %.3f' % (tuple(errors)+(poiss.ts,bandts)),
            self.assertAlmostEquals(checks[0], errors[0], delta=1e-2)
            self.assertAlmostEquals(checks[1], errors[1], delta=1e-2)
            self.assertAlmostEquals(checks[2], poiss.ts, delta=100.)
            self.assertAlmostEquals(checks[3], bandts, delta=1.0)

class TestLocalization(TestSetup):
    def setUp(self):
        setup('likeviews')
        print 'initial loglike: %.1f ...' % likeviews.log_like() ,

    def test(self):
        """--> generate a view, check values, and restore"""
        init = likeviews.log_like()
        with likeviews.tsmap_view('P7R42722') as tsm:
            self.assertEquals(0, tsm())
            self.assertAlmostEquals(-0.297, tsm((266.6, -28.86)), delta=0.1)
        self.assertAlmostEquals(init, likeviews.log_like(), delta=0.1)
        
class TestUser(TestSetup):
    def setUp(self):
        setup('user')
        self.init = user.log_like()

    def tearDown(self):
        self.assertAlmostEquals(self.init, user.log_like())

    def test_fit(self, selects=(0, '_Norm', None), 
            expects=(697608.2, 697617.7 ,697645.4)):
        for select, expect in zip(selects, expects):
            wfit, pfit = user.fit(select, summarize=False)
            self.assertAlmostEquals(expect, wfit, delta=0.1)

    def test_localization(self, source_name='P7R42722'):
        t = user.localize(source_name, quiet=True)
        self.assertAlmostEquals(0.0062, t['a'], delta=1e-3)
        self.assertAlmostEquals(0.210, t['qual'], delta=1e-3)
 
    def testTS(self, source_name='P7R42722', expect=800):
        """-->compute a Test Statistic"""
        ts = user.TS(source_name)
        self.assertAlmostEquals(expect, ts, delta=1)
        
    def testSED(self, source_name='P7R42722'):
        """-->measure a full SED (not yet)"""
        pass
        
    
    
class TestAssociations(TestSetup):
    def test(self):
        assoc = associate.SrcId()
        t = assoc('test', (266.5980,  -28.8680), 0.01)
        self.assertTrue(set(['ra', 'deltats', 'ang', 'name', 'prior', 'density', 'dec', 'prob',
                'dir', 'cat']).issubset(t.keys()))
        self.assertAlmostEquals(2.614, t['deltats'][0], delta=0.001)
     
    
test_cases = (TestConfig, 
    TestPoint, TestDiffuse, 
   # TestExtended, 
    TestROImodel, 
    TestBands, 
    TestLikelihood,
    TestFitterView,
    TestSED,
    TestLocalization,
    TestAssociations,
    #### no memory to do this at the same time since it creates duplicate large objects
    #TestUser,
    )
    
def run(t='all', loader=unittest.TestLoader(), debug=False): 
    if t=='all':
        suite = unittest.TestSuite()
        for test_class in test_cases:
            tests = loader.loadTestsFromTestCase(test_class)
            suite.addTests(tests)
    else:
        suite = loader.loadTestsFromTestCase(t)
    print 'running %d tests %s' % (suite.countTestCases(), 'in debug mode' if debug else '') 
    if debug:
        suite.debug()
    else:
        unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)
    
    
if __name__=='__main__':
    run()