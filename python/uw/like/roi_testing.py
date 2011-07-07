"""
Module to perfrom routine testing of pointlike's many features.'

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_testing.py,v 1.2 2011/07/06 05:03:19 lande Exp $

author: Matthew Kerr, Toby Burnett, Joshua Lande
"""
import datetime
import os

import numpy as np

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification,SpectralAnalysis
from uw.like.pointspec_helpers import PointSource,get_diffuse_source
from uw.like.Models import PowerLaw
from uw.like.SpatialModels import Disk,Gaussian
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_monte_carlo import SpectralAnalysisMC


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
                 irf='P6_V3_DIFFUSE',
                 binsperdec = 4,
                 mc_energy=True,
                 tstart=0,
                 tstop=datetime.timedelta(days=7).total_seconds(),
                 quiet=True
                )


    global roi # helps with debugging
    roi=sa.roi(roi_dir=center,
               point_sources=point_sources,
               diffuse_sources=diffuse_sources,
               fit_emin = emin, fit_emax = emax,
              )

    return roi

def compare_model(fit,true):
    norm_mc=true.model['Norm']
    index_mc=true.model['Index']
    [norm,index],[norm_err,index_err]=fit.model.statistical(absolute=True)
    print ' > True Norm = %.2e, Fit Norm = %.2e +/- %.2e (pull=%.1f)' % (norm_mc,norm,norm_err,(norm-norm_mc)/norm_err)
    print ' > True Index = %.2f, Fit Index = %.2f +/- %.2f (pull=%.1f)' % (index_mc,index,index_err,(index-index_mc)/index_err)

def compare_spatial_model(fit,true,lsigma):
    """ Compare a source 'fit's spatial model to the source 'true's spatial model. """

    if hasattr(true,'spatial_model') and hasattr(fit,'spatial_model'):
        sigma_mc=true.spatial_model['Sigma']
        sigma=fit.spatial_model['Sigma']
        sigma_err=fit.spatial_model.statistical(absolute=True)[1][2]
        print ' > True Ext = %.2f, Fit Ext = %.2f +/- %.1g%% (pull=%.1f)' % (
            sigma_mc,sigma,sigma_err,(sigma-sigma_mc)/sigma_err)

    print ' > True Pos = (%.3f,%.3f), Fit Pos = (%.3f,%.3f), dist=%.3f, err=%.3f (pull=%.1f)' % \
            (true.skydir.ra(),true.skydir.dec(),
             fit.skydir.ra(),fit.skydir.dec(),
             np.degrees(true.skydir.difference(fit.skydir)),lsigma,
             np.degrees(true.skydir.difference(fit.skydir))/lsigma)


def test_models():

    print 'Testing the Model + SpatialModel objects.'
    print

    model=PowerLaw(p=[1e-7,2])
    print '> Norm = 1e-7 = %s' % model['Norm']

    model['Norm']=1e-6
    print '> Norm = 1e-6 = %s' % model['Norm']

    print '> Index = 2 = %s' % model['Index']

    spatial_model = Disk(p=[1.5],center=SkyDir(22,22,SkyDir.GALACTIC))
    print '> Sigma = 1.5 = %s' % spatial_model['Sigma']
    
    spatial_model['Sigma'] = 0.5
    print '> Sigma = 0.5 = %s' % spatial_model['Sigma']

def test_extended_source():

    print 'Analyze a simulated extended source against an isotropic background (E>10GeV)'

    center=SkyDir(0,0,SkyDir.GALACTIC)

    # Sreekumar-like isotropic
    point_sources=[]
    diffuse_sources=[
        get_diffuse_source('ConstantValue',None,'PowerLaw',None,'Isotropic Diffuse')
    ]

    model = PowerLaw(p=[1,2])
    model.set_flux(1e-5)

    spatial_model = Gaussian(p=[1],center=center)

    es_mc = ExtendedSource(name='source',spatial_model=spatial_model,model=model)
    es_fit = es_mc.copy()
    diffuse_sources.append(es_fit)

    roi = get_roi('extended_test',center,point_sources,diffuse_sources, emin=1e4)

    roi.fit(use_gradient=True)
    roi.fit_extension(which='source')
    roi.localize(update=True)
    roi.fit(use_gradient=True)

    compare_model(es_fit,es_mc)
    compare_spatial_model(es_fit,es_mc,roi.lsigma)


def test_point_source():

    print 'Analyze a simulated point source against an isotropic background'

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

    roi = get_roi('point_test',center,point_sources,diffuse_sources)

    roi.fit(use_gradient=True)
    roi.localize(update=True)
    roi.fit(use_gradient=True)

    compare_model(ps_fit,ps_mc)
    compare_spatial_model(ps_fit,ps_mc,roi.lsigma)

if __name__ == '__main__':

    print 'Performing Automated tests of Pointlike'
    print

    # Create/store files in $SIMDIR
    if not os.environ.has_key('SIMDIR'):
        raise Exception('$SIMDIR must be defiend.')
    if not os.path.exists(os.path.expandvars("$SIMDIR")):
        raise Exception('$SIMDIR must exist.')

    #np.seterr(all='raise')
    np.seterr(all='warn')

    test_models()
    test_point_source()
    test_extended_source()
