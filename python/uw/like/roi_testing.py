"""
Module to perfrom routine testing of pointlike's many features.'

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_testing.py,v 1.8 2011/07/24 02:53:35 lande Exp $

author: Matthew Kerr, Toby Burnett, Joshua Lande
"""
import os
import unittest

import numpy as np

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification,SpectralAnalysis
from uw.like.pointspec_helpers import PointSource,get_diffuse_source
from uw.like.Models import PowerLaw,LogParabola,SumModel,ProductModel
from uw.like.SpatialModels import Disk,Gaussian,SpatialMap
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_monte_carlo import SpectralAnalysisMC



class PointlikeTest(unittest.TestCase):


    def setUp(self):

        # Create/store files in $SIMDIR
        self.assertTrue(os.environ.has_key('SIMDIR'),'$SIMDIR must be defiend.')
        self.assertTrue(os.path.exists(os.environ['SIMDIR']),'$SIMDIR must exist.')

        self.MAX_ALLOWED_PULL = 3

    def compare_model(self,fit,true):
        norm_pull=(fit.model['Norm']-true.model['Norm'])/fit.model.error('Norm')
        #print 'norm_pull = ',norm_pull
        self.assertTrue(abs(norm_pull)<self.MAX_ALLOWED_PULL,'norm pull=%.1f is bad.' % norm_pull)

        index_pull=(fit.model['Index']-true.model['Index'])/fit.model.error('Index')
        #print 'index_pull = ',index_pull
        self.assertTrue(abs(index_pull)<self.MAX_ALLOWED_PULL,'index pull=%.1f is bad.' % index_pull)

    def compare_spatial_model(self,fit,true,lsigma):
        """ Compare a source 'fit's spatial model to the source 'true's spatial model. """

        if hasattr(true,'spatial_model') and hasattr(fit,'spatial_model'):
            sigma_pull = (fit.spatial_model['Sigma'] - true.spatial_model['Sigma'])/fit.spatial_model.error('Sigma')
            #print 'sigma_pull = ',sigma_pull
            self.assertTrue(abs(sigma_pull)<self.MAX_ALLOWED_PULL,'sigma pull=%.1f is bad.' % sigma_pull)

        dist_pull = np.degrees(true.skydir.difference(fit.skydir))/lsigma
        #print 'dist_pull = ',dist_pull
        self.assertTrue(abs(dist_pull)<self.MAX_ALLOWED_PULL,'dist pull=%.1f is bad.' % dist_pull)

    @staticmethod
    def get_roi(name,center,point_sources,diffuse_sources,emin=1e2,emax=1e5):

        ft1='$SIMDIR/%s_ft1.fits' % name
        ft2='$SIMDIR/%s_ft2.fits' % name

        if os.path.exists(os.path.expandvars(ft1)) and \
           os.path.exists(os.path.expandvars(ft2)):
            sa_object=SpectralAnalysis
        else:
            sa_object=SpectralAnalysisMC

        ds=DataSpecification(ft1files=ft1,
                             ft2files=ft2,
                             ltcube='$SIMDIR/%s_ltcube.fits' % name, 
                             binfile='$SIMDIR/%s_binfile.fits' % name
                            )

        sa=sa_object(ds,
                     irf='P7SOURCE_V6',
                     binsperdec = 2,
                     mc_energy=True,
                     tstart=0,
                     tstop=604800, # 7 days
                     quiet=True,
                     maxROI =5, minROI = 5,
                    )


        roi=sa.roi(roi_dir=center,
                   point_sources=point_sources,
                   diffuse_sources=diffuse_sources,
                   fit_emin = emin, fit_emax = emax,
                  )

        return roi

    #@unittest.skip("skip")
    def test_models(self):

        print '\nTesting the uw.like.Models.Model object.\n'

        model=PowerLaw(p=[1e-7,2])

        self.assertAlmostEqual(1e-7,model['Norm'])
        self.assertAlmostEqual(1e-7,model['norm'],msg='case insensitive')
        model['Norm']=1e-6
        self.assertAlmostEqual(1e-6,model['Norm'])

        model=PowerLaw(index=1)
        self.assertAlmostEqual(1,model['index'],msg="test kwargs passing to object")
        self.assertAlmostEqual(1e-11,model['norm'], msg="test default value")

        copy_model=model.copy()
        self.assertAlmostEqual(1e-11,copy_model['norm'], msg="test default value")
        model['norm'] = 1e-10
        self.assertAlmostEqual(1e-11,copy_model['norm'], msg="unchanged by copy")

        model=PowerLaw(random_input=3)
        self.assertEqual(model.random_input,3,'test random kwarg')

    #@unittest.skip("skip")
    def test_sum_product_models(self):

        print '\nTesting the uw.like.Models.Model object.\n'

        m1=PowerLaw(index=2);m1.set_flux(1)
        m2=LogParabola(beta=2); m1.set_flux(1)

        sum_model=SumModel(m1,m2)
        prod_model=ProductModel(m1,m2)

        for energy in [1e2,1e3,1e4]:
            self.assertAlmostEqual(sum_model(energy),m1(energy)+m2(energy))
            self.assertAlmostEqual(prod_model(energy),m1(energy)*m2(energy))

    #@unittest.skip("skip")
    def test_spatial_models(self):

        print '\nTesting the uw.like.SpatialModels.SpatialModel object.\n'

        # test default values
        spatial_model = Disk()
        self.assertAlmostEqual(0.1,spatial_model['Sigma']) 
        self.assertAlmostEqual(0,spatial_model['ra']) 
        self.assertAlmostEqual(0,spatial_model['dec']) 

        spatial_model = Disk(p=[1.5],center=SkyDir(22,22,SkyDir.GALACTIC))
        self.assertAlmostEqual(1.5,spatial_model['Sigma']) # test getitem/setitem 
        
        spatial_model['Sigma'] = 0.5
        self.assertAlmostEqual(0.5,spatial_model['Sigma']) 
        spatial_model['sigma'] = 0.5
        self.assertAlmostEqual(0.5,spatial_model['sigma']) # case insensitive
        self.assertAlmostEqual(22,spatial_model.center.l())
        self.assertAlmostEqual(22,spatial_model.center.b())

        spatial_model = Disk(p=[1.5],center=SkyDir(22,22,SkyDir.GALACTIC),coordsystem=SkyDir.GALACTIC)
        self.assertAlmostEqual(22,spatial_model['l'])
        self.assertAlmostEqual(22,spatial_model['b'])

        # test what happens when only center is specified
        spatial_model = Disk(center=SkyDir(22,22))
        self.assertAlmostEqual(0.1,spatial_model['sigma'])
        self.assertAlmostEqual(22,spatial_model['ra'])
        self.assertAlmostEqual(22,spatial_model['dec'])

        spatial_model = Disk(sigma=1.5, l=22, b=22)
        self.assertAlmostEqual(1.5,spatial_model['Sigma'])
        self.assertAlmostEqual(22,spatial_model['l'])
        self.assertAlmostEqual(22,spatial_model['b'])

        self.assertRaises(Exception,spatial_model.__setitem__,['Sigma',-1])

        spatial_model['l'] = -100
        self.assertAlmostEqual(-100, spatial_model['l'])
        for k in ['ra','dec']:
            self.assertRaises(Exception,spatial_model.__getitem__,k)
            self.assertRaises(Exception,spatial_model.__setitem__,[k,0])

        spatial_model = Disk(p=1.5, ra=22, dec=22)
        self.assertAlmostEqual(1.5,spatial_model['Sigma'])
        self.assertAlmostEqual(22,spatial_model['ra'])
        self.assertAlmostEqual(22,spatial_model['dec'])
        for k in ['l','b']:
            self.assertRaises(Exception,spatial_model.__getitem__,k)
            self.assertRaises(Exception,spatial_model.__setitem__,[k,0])

        spatial_model = Disk(p=[1.5], l=22, b=22)
        self.assertAlmostEqual(1.5,spatial_model['Sigma'])

        copy_model=spatial_model.copy()
        self.assertAlmostEqual(1.5,copy_model['Sigma'])
        spatial_model['Sigma'] = 0.5
        self.assertAlmostEqual(1.5,copy_model['Sigma']) # unchanged by copy

        model=Disk(random_input=3)
        self.assertEqual(model.random_input,3,'test random kwarg')

    #@unittest.skip("skip")
    def test_extended_source(self):

        print '\nAnalyze a simulated extended source against an isotropic background (E>10GeV)\n'

        center=SkyDir(0,0,SkyDir.GALACTIC)

        # Sreekumar-like isotropic
        point_sources=[]
        diffuse_sources=[
            get_diffuse_source('ConstantValue',None,'PowerLaw',None,'Isotropic Diffuse')
        ]

        model = PowerLaw(p=[1,2])
        model.set_flux(1e-4)

        spatial_model = Gaussian(p=[1],center=center)

        es_mc = ExtendedSource(name='source',spatial_model=spatial_model,model=model)
        es_fit = es_mc.copy()
        diffuse_sources.append(es_fit)

        roi = PointlikeTest.get_roi('extended_test',center,point_sources,diffuse_sources, emin=1e4)
        global roi_ext;roi_ext=roi # helps with debugging

        #roi.modify(which='source',sigma=0.3)
        roi.modify(which='source',spatial_model=Gaussian(0.3))

        roi.fit(use_gradient=True)
        roi.fit_extension(which='source')
        roi.localize(update=True)
        roi.fit(use_gradient=True)

        self.compare_model(es_fit,es_mc)
        self.compare_spatial_model(es_fit,es_mc,roi.lsigma)

        self.assertTrue(roi.TS(which='source')>25,'The source should be significant')
        self.assertTrue(roi.TS_ext(which='source')>25,'And significantly extended')

        es_mc.spatial_model.save_template('$SIMDIR/extended_template.fits')

        template_source=ExtendedSource(
            name='template_source',
            model=es_mc.model,
            spatial_model=SpatialMap(file='$SIMDIR/extended_template.fits')
        )

        roi.del_source(which='source')
        roi.add_source(template_source)

        roi.fit()

        self.compare_model(template_source,es_mc)

        self.assertTrue(roi.TS(which='template_source')>25,'Make sure these functions work similary with spatial_map')


    #@unittest.skip("skip")
    def test_point_source(self):

        print '\nAnalyze a simulated point source against an isotropic background\n'

        center=SkyDir(0,0)

        # Sreekumar-like isotropic
        diffuse_sources=[
            get_diffuse_source('ConstantValue',None,'PowerLaw',None,'Isotropic Diffuse')
        ]

        model = PowerLaw(p=[1,2])
        model.set_flux(1e-6)
        ps_mc = PointSource(name='source',skydir=center,model=model)
        ps_fit = ps_mc.copy()
        point_sources=[ps_fit]

        roi = PointlikeTest.get_roi('point_test',center,point_sources,diffuse_sources)
        global roi_pt;roi_pt=roi # helps with debugging

        roi.fit(use_gradient=True)
        roi.localize(update=True)
        roi.fit(use_gradient=True)

        self.compare_model(ps_fit,ps_mc)
        self.compare_spatial_model(ps_fit,ps_mc,roi.lsigma)

if __name__ == '__main__':

    print '\nPerforming Automated tests of Pointlike\n'

    import numpy as np
    np.seterr(all='warn')

    unittest.main()
