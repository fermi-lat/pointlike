"""
Module to perfrom routine testing of pointlike's many features.'

$Header: $

author: Matthew Kerr, Toby Burnett, Joshua Lande
"""
import datetime
import os

import numpy as np

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification,SpectralAnalysis
from uw.like.pointspec_helpers import PointSource,get_diffuse_source
from uw.like.Models import PowerLaw
from uw.like.roi_monte_carlo import SpectralAnalysisMC


def test_point_source_fit():

    print 'Analyze a simulated point source agains an isotropic background'

    center=SkyDir(0,0)

    ft1='$SIMDIR/test_point_source_ft1.fits'
    ft2='$SIMDIR/test_point_source_ft2.fits'

    if os.path.exists(os.path.expandvars(ft1)) and \
       os.path.exists(os.path.expandvars(ft2)):
        sa_object=SpectralAnalysis
    else:
        sa_object=SpectralAnalysisMC

    ds=DataSpecification(ft1files=ft1,
                         ft2files=ft2,
                         ltcube='$SIMDIR/test_point_source_ltcube.fits',
                         binfile='$SIMDIR/test_point_source_binfile.fits'
                        )

    sa=sa_object(ds,
                 irf='P6_V3_DIFFUSE',
                 emin = 1e2,
                 emax = 1e5,
                 binsperdec = 4,
                 mc_energy=True,
                 tstart=0,
                 tstop=datetime.timedelta(days=7).total_seconds(),
                 quiet=True
                )

    # Sreekumar-like isotropic
    diffuse_sources=[
        get_diffuse_source('ConstantValue',None,'PowerLaw',None,'Isotropic Diffuse')
    ]

    model = PowerLaw(p=[1,2])
    model.set_flux(1e-6)
    ps_mc = PointSource(name='source',skydir=center,model=model)
    ps_fit = ps_mc.copy()
    point_sources=[ps_fit]

    roi=sa.roi(roi_dir=center,
               point_sources=point_sources,
               diffuse_sources=diffuse_sources,
              )

    roi.fit(use_gradient=True)
    roi.localize(update=True)
    roi.fit(use_gradient=True)

    [norm_mc,index_mc]=ps_mc.model.statistical()[0]
    [norm,index],[norm_err,index_err]=ps_fit.model.statistical(absolute=False)

    print ' > True Norm = %.2e, Fit Norm = %.2e +/- %.1g%%' % (norm_mc,norm,norm_err)
    print ' > True Index = %.2f, Fit Index = %.2f +/- %.1g%%' % (index_mc,index,index_err)
    print ' > True Pos = (%.3f,%.3f), Fit Pos =  (%.3f,%.3f), dist=%.3f, err=%.3f' % \
            (ps_mc.skydir.ra(),ps_mc.skydir.dec(),
             ps_fit.skydir.ra(),ps_fit.skydir.dec(),
             np.degrees(ps_mc.skydir.difference(ps_fit.skydir)),roi.lsigma)

if __name__ == '__main__':

    print 'Performing Automated tests of Pointlike'
    print

    # Create/store files in $SIMDIR
    if not os.environ.has_key('SIMDIR'):
        raise Exception('$SIMDIR must be defiend.')
    if not os.path.exists(os.path.expandvars("$SIMDIR")):
        raise Exception('$SIMDIR must exist.')

    test_point_source_fit()
