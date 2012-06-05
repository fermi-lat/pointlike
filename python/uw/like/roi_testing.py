"""
Module to perfrom routine testing of pointlike's many features.'

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_testing.py,v 1.10 2012/01/14 00:13:57 lande Exp $

author: Matthew Kerr, Toby Burnett, Joshua Lande
"""
import os
import unittest

import numpy as np

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification,SpectralAnalysis
from uw.like.pointspec_helpers import PointSource,get_diffuse_source,get_default_diffuse
from uw.like.Models import PowerLaw,LogParabola,SumModel,ProductModel
from uw.like.SpatialModels import Disk,Gaussian,SpatialMap
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_monte_carlo import SpectralAnalysisMC



class PointlikeTest(unittest.TestCase):


    MAX_ALLOWED_PULL = 3
    USE_GRADIENT = False
    VERBOSE = False
    
    @staticmethod
    def p(str):
        print '\n\n%s %s %s\n\n' % ('#'*5,str,'#'*5)

    def setUp(self):

        # Create/store files in $SIMDIR
        self.assertTrue(os.environ.has_key('SIMDIR'),'$SIMDIR must be defiend.')
        self.assertTrue(os.path.exists(os.environ['SIMDIR']),'$SIMDIR must exist.')

    def compare_model(self,fit,true):
        norm_pull=(fit.model['Norm']-true.model['Norm'])/fit.model.error('Norm')
        if PointlikeTest.VERBOSE:
            print 'norm_pull = ',norm_pull
        self.assertTrue(abs(norm_pull)<self.MAX_ALLOWED_PULL,'norm pull=%.1f is bad.' % norm_pull)

        index_pull=(fit.model['Index']-true.model['Index'])/fit.model.error('Index')
        if PointlikeTest.VERBOSE:
            print 'index_pull = ',index_pull
        self.assertTrue(abs(index_pull)<self.MAX_ALLOWED_PULL,'index pull=%.1f is bad.' % index_pull)

    def compare_spatial_model(self,fit,true,lsigma):
        """ Compare a source 'fit's spatial model to the source 'true's spatial model. """

        if hasattr(true,'spatial_model') and hasattr(fit,'spatial_model'):
            sigma_pull = (fit.spatial_model['Sigma'] - true.spatial_model['Sigma'])/fit.spatial_model.error('Sigma')
            if PointlikeTest.VERBOSE:
                print 'sigma_pull = ',sigma_pull
            self.assertTrue(abs(sigma_pull)<self.MAX_ALLOWED_PULL,'sigma pull=%.1f is bad.' % sigma_pull)

        dist_pull = np.degrees(true.skydir.difference(fit.skydir))/lsigma
        if PointlikeTest.VERBOSE:
            print 'dist_pull = ',dist_pull
        self.assertTrue(abs(dist_pull)<self.MAX_ALLOWED_PULL,'dist pull=%.1f is bad.' % dist_pull)

    @staticmethod
    def get_roi(name,center,point_sources,diffuse_sources,emin=1e2,emax=1e5):

        ft1='$SIMDIR/%s_ft1.fits' % name
        ft2='$SIMDIR/%s_ft2.fits' % name

        ds=DataSpecification(ft1files=ft1,
                             ft2files=ft2,
                             ltcube='$SIMDIR/%s_ltcube.fits' % name, 
                             binfile='$SIMDIR/%s_binfile.fits' % name
                            )

        sa=SpectralAnalysisMC(ds,
                              seed=0,
                              emin=emin,
                              emax=emax,
                              irf='P7SOURCE_V6',
                              binsperdec = 2,
                              mc_energy=True,
                              tstart=0,
                              tstop=604800, # 7 days
                              quiet=not PointlikeTest.VERBOSE,
                              roi_dir=center,
                              maxROI=5, minROI=5,
                             )


        roi=sa.roi(roi_dir=center,
                   point_sources=point_sources,
                   diffuse_sources=diffuse_sources)

        return roi


    #@unittest.skip("skip")
    def test_extended_source(self):

        PointlikeTest.p('USE_GRADIENT=%s' % PointlikeTest.USE_GRADIENT)

        if PointlikeTest.VERBOSE:
            PointlikeTest.p('Analyze a simulated extended source against an isotropic background (E>10GeV)')

        center=SkyDir(0,0,SkyDir.GALACTIC)

        # Sreekumar-like isotropic
        point_sources=[]
        diffuse_sources=[
            get_diffuse_source('ConstantValue',None,'PowerLaw',None,'Isotropic Diffuse')
        ]

        model = PowerLaw(p=[1,2])
        model.set_flux(1e-4)

        if PointlikeTest.VERBOSE:
            PointlikeTest.p('Simulating gaussian source with sigma=1 degrees')
        spatial_model = Gaussian(p=[1],center=center)

        es_mc = ExtendedSource(name='source',spatial_model=spatial_model,model=model)
        es_fit = es_mc.copy()
        diffuse_sources.append(es_fit)

        roi = PointlikeTest.get_roi('extended_test',center,point_sources,diffuse_sources, emin=1e4)
        global roi_ext;roi_ext=roi # helps with debugging

        if PointlikeTest.VERBOSE:
            print roi

        if PointlikeTest.VERBOSE:
            PointlikeTest.p('Setting initial spatial model to 0.3 degrees')
        roi.modify(which='source',spatial_model=Gaussian(0.3))

        if PointlikeTest.VERBOSE: print roi
        roi.fit(use_gradient=PointlikeTest.USE_GRADIENT)
        if PointlikeTest.VERBOSE: print roi
        roi.fit_extension(which='source', use_gradient=PointlikeTest.USE_GRADIENT)
        roi.localize(update=True)
        roi.fit(use_gradient=PointlikeTest.USE_GRADIENT)

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

        if PointlikeTest.VERBOSE:
            print '\nAnalyze a simulated point source against the galactic + isotropic diffuse\n'

        center=SkyDir(0,0)

        diffuse_sources=get_default_diffuse(
            diffdir='$GLAST_EXT/diffuseModels/v2r0p1/',
            gfile='ring_2year_P76_v0.fits',
            ifile='isotrop_2year_P76_source_v1.txt')

        model = PowerLaw(p=[1,2])
        model.set_flux(1e-6)
        ps_mc = PointSource(name='source',skydir=center,model=model)
        ps_fit = ps_mc.copy()
        point_sources=[ps_fit]

        roi = PointlikeTest.get_roi('point_test',center,point_sources,diffuse_sources)
        global roi_pt;roi_pt=roi # helps with debugging

        if PointlikeTest.VERBOSE:
            print roi

        roi.fit(use_gradient=PointlikeTest.USE_GRADIENT)
        if PointlikeTest.VERBOSE: print roi
        roi.localize(update=True)
        roi.fit(use_gradient=PointlikeTest.USE_GRADIENT)
        if PointlikeTest.VERBOSE: 
            roi.print_summary()
            print roi

        self.compare_model(ps_fit,ps_mc)
        self.compare_spatial_model(ps_fit,ps_mc,roi.lsigma)

if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", default=False,help="Output more verbosely")
    parser.add_argument("--use-gradient", default=False, action='store_true')
    args=parser.parse_args()
    sys.argv = sys.argv[0:1]

    PointlikeTest.USE_GRADIENT = args.use_gradient
    PointlikeTest.VERBOSE = args.verbose

    if PointlikeTest.VERBOSE: 
        PointlikeTest.p('Performing Automated tests of Pointlike')

    import numpy as np
    np.seterr(all='ignore')

    unittest.main()
