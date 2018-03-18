"""
Module to perfrom routine testing of pointlike's many features.'

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_testing.py,v 1.17 2012/07/16 16:44:09 lande Exp $

author: Matthew Kerr, Toby Burnett, Joshua Lande
"""
import os
from os.path import join,abspath
import sys
import unittest

import numpy as np

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification,SpectralAnalysis
from uw.like.pointspec_helpers import PointSource,get_diffuse_source,get_default_diffuse
from uw.like.Models import PowerLaw,PowerLawFlux,LogParabola,SumModel,ProductModel,FileFunction
from uw.like.SpatialModels import Disk,Gaussian,SpatialMap
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_monte_carlo import SpectralAnalysisMC

from uw.utilities import path


class PointlikeTest(unittest.TestCase):


    MAX_ALLOWED_PULL = 3
    USE_GRADIENT = False
    VERBOSE = False
    
    @staticmethod
    def p(str):
        print '\n\n%s %s %s\n\n' % ('#'*5,str,'#'*5)

    def setUp(self):

        # Create/store files in $SIMDIR
        self.assertTrue(os.environ.has_key('SIMDIR'),"""$SIMDIR must be defiend. If it does not exist, please make a new folder and set $SIMDIR to point to it. All simulated data will be put into it.""")
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
    def get_roi(name,center,point_sources,diffuse_sources,emin=1e2,emax=1e5,binsperdec=2):

        simdir=path.expand('$SIMDIR/%s' % name)

        if not os.path.exists(simdir):
            os.makedirs(simdir)


        ft1='%s/ft1.fits' % simdir
        ft2='%s/ft2.fits' % simdir

        ds=DataSpecification(ft1files=ft1,
                             ft2files=ft2,
                             ltcube='%s/ltcube.fits' % simdir, 
                             binfile='%s/binfile.fits' % simdir
                            )

        sa=SpectralAnalysisMC(ds,
                              seed=0,
                              emin=emin,
                              emax=emax,
                              irf='P7SOURCE_V6',
                              binsperdec = binsperdec,
                              mc_energy=True,
                              tstart=0,
                              tstop=604800, # 7 days
                              quiet=not PointlikeTest.VERBOSE,
                              roi_dir=center,
                              maxROI=5, minROI=5,
                              savedir='%s/gtobssim' % simdir)


        roi=sa.roi(roi_dir=center,
                   point_sources=point_sources,
                   diffuse_sources=diffuse_sources)

        return roi


    @unittest.skipIf("--skip-extended" in sys.argv,'Skip time consuming extended source test')
    def test_extended_source(self):

        PointlikeTest.p('USE_GRADIENT=%s' % PointlikeTest.USE_GRADIENT)

        if PointlikeTest.VERBOSE:
            PointlikeTest.p('Analyze a simulated extended source against an isotropic background (E>10GeV)')

        center=SkyDir(0,0)

        # Sreekumar-like isotropic
        point_sources=[]
        diffuse_sources=[
            get_diffuse_source('ConstantValue',None,'PowerLaw',None,'Isotropic Diffuse')
        ]

        model = PowerLaw(index=2)
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

        if PointlikeTest.VERBOSE:
                PointlikeTest.p('Now, switching from Disk soruce to template source.')

        roi.del_source(which='source')
        template_source=ExtendedSource(
            name='template_source',
            model=es_mc.model,
            spatial_model=SpatialMap(file='$SIMDIR/extended_template.fits')
        )

        roi.add_source(template_source)

        roi.fit(use_gradient=PointlikeTest.USE_GRADIENT)

        self.compare_model(template_source,es_mc)

        self.assertTrue(roi.TS(which='template_source')>25,'Make sure these functions work similary with spatial_map')


    @unittest.skipIf("--skip-ps1" in sys.argv,'skip')
    def test_ps1(self):

        if PointlikeTest.VERBOSE:
            print '\nAnalyze a simulated point source against the galactic + isotropic diffuse\n'

        center=SkyDir(0,0)

        diffuse_sources=get_default_diffuse(
            diffdir='$GLAST_EXT/diffuseModels/v2r0p1/',
            gfile='ring_2year_P76_v0.fits',
            ifile='isotrop_2year_P76_source_v1.txt')

        model = PowerLaw(index=2)
        model.set_flux(1e-6)
        ps_mc = PointSource(name='source',skydir=center,model=model)
        ps_fit = ps_mc.copy()
        point_sources=[ps_fit]

        roi = PointlikeTest.get_roi('ps1',center,point_sources,diffuse_sources)
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


    @unittest.skipIf("--skip-ps2" in sys.argv,'skip')
    def test_ps2(self):

        if PointlikeTest.VERBOSE:
            PointlikeTest.p('\nAnalyze a simulated point source\n')

        center=SkyDir(0,0)

        model = PowerLawFlux(int_flux=1e-5, index=3, emin=1e2, emax=1e5)
        ps = PointSource(name='source',skydir=center,model=model)
        point_sources=[ps]
        diffuse_sources=None

        roi = PointlikeTest.get_roi('ps2',center,point_sources,diffuse_sources, emin=1e2, emax=1e5, binsperdec=4)
        global roi_pt;roi_pt=roi # helps with debugging

        if PointlikeTest.VERBOSE:
            print roi

        if PointlikeTest.VERBOSE:
            PointlikeTest.p('\nTesting if SED points are within errors\n')

        bf=roi.plot_sed(which='source', filename='sed.pdf', merge=False)

        elow=bf.rec.elow
        ehigh=bf.rec.ehigh
        e = np.sqrt(elow*ehigh) # geometric mean

        flux=bf.rec.flux
        significant=bf.rec.flux > 0 # pointlike convention for upper limit
        uflux = bf.rec.uflux
        lflux = bf.rec.lflux
        true_flux=e**2*model(e)
        pull = np.where(flux>true_flux, (flux-true_flux)/uflux, (true_flux-flux)/lflux)

        self.assertTrue(np.all(pull[significant] < PointlikeTest.MAX_ALLOWED_PULL), 'All SED points within errors of true value')
        if PointlikeTest.VERBOSE:
            PointlikeTest.p('\nPull for SED points is %s\n' % pull[significant])

        # I am not sure this strictly has to always be true. Maybe there is a better test - J.L.
        self.assertTrue(np.all(uflux[~significant] > true_flux[~significant]), 'All SED upper limits must be above spectra')

    @unittest.skipIf("--skip-ff" in sys.argv,'skip')
    def test_ff(self):
        """ Simulate from a filefunction object and test that the best
        fit flux
            is consistent with the simulated flux. """
        name='ff'

        model = PowerLaw(index=2)
        model.set_flux(1e-6)
        simdir=path.expand('$SIMDIR/%s' % name)
        if not os.path.exists(simdir):
            os.makedirs(simdir)

        filename = abspath(join(simdir,'file_function.txt'))
        model.save_profile(filename,10,1e6)
        ff = FileFunction(file=filename)

        center=SkyDir(0,0)
        ps = PointSource(name='source',skydir=center,model=ff)
        point_sources=[ps]
        diffuse_sources=None
        roi = PointlikeTest.get_roi(name,center,point_sources,diffuse_sources, emin=1e2, emax=1e5, binsperdec=4)

        if PointlikeTest.VERBOSE:
            roi.print_summary()
            print roi

        roi.fit(use_gradient=PointlikeTest.USE_GRADIENT)

        if PointlikeTest.VERBOSE:
            roi.print_summary()
            print roi

        fit,error = ff.i_flux(1e2,1e5,error=True)
        true= model.i_flux(1e2,1e5,error=False)
        self.assertPull(fit,true,error,'flux')

    def assertPull(self,fit,true,error,message):
        pull = abs((fit-true)/error)
        if PointlikeTest.VERBOSE:
            PointlikeTest.p('%s, pull=%s' % (message,pull))
        self.assertTrue(pull<self.MAX_ALLOWED_PULL,'pull=%.1f is bad.' % pull)


if __name__ == '__main__':
    from argparse import ArgumentParser

    # Make sure not to open plots to screen
    from matplotlib import rcParams
    rcParams['backend'] = 'Agg'

    parser = ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", default=False,help="Output more verbosely")
    parser.add_argument("--use-gradient", default=False, action='store_true')
    parser.add_argument("--skip-extended", default=False, action='store_true')
    parser.add_argument("--skip-ps1", default=False, action='store_true')
    parser.add_argument("--skip-ps2", default=False, action='store_true')
    args=parser.parse_args()
    sys.argv = sys.argv[0:1]

    PointlikeTest.USE_GRADIENT = args.use_gradient
    PointlikeTest.VERBOSE = args.verbose

    if PointlikeTest.VERBOSE: 
        PointlikeTest.p('Performing Automated tests of Pointlike')

    import numpy as np
    np.seterr(all='ignore')

    unittest.main(verbosity=2 if args.verbose else 0)
