"""
Module implements a wrapper around gtobssim to allow
less painful simulation of data.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_monte_carlo.py,v 1.47 2012/03/17 00:10:17 lande Exp $

author: Joshua Lande
"""
import os
import re
import types
from textwrap import dedent
import shutil
import collections
from tempfile import mkdtemp
from GtApp import GtApp

import pyfits
import numpy as np
import pywcs

from skymaps import IsotropicSpectrum,IsotropicPowerLaw,DiffuseFunction,\
        PySkyFunction,Hep3Vector,SkyImage,SkyDir,PythonUtilities,IsotropicConstant
from . SpatialModels import Gaussian,EllipticalGaussian,RadiallySymmetricModel
from . Models import PowerLaw,PowerLawFlux,Constant,FileFunction
from . pointspec import DataSpecification,SpectralAnalysis
from . pointspec_helpers import PointSource,PointSourceCatalog
from . roi_diffuse import DiffuseSource
from . roi_extended import ExtendedSource
from . SpatialModels import Gaussian,EllipticalGaussian,SpatialModel,RadiallySymmetricModel,SpatialMap

from uw.utilities import keyword_options
from uw.utilities.keyword_options import decorate, process, change_defaults, get_default, get_row
from uw.utilities.fitstools import rad_mask

import pyLikelihood

# I think this is the smallest possible float for 32 bit numbers
# http://www.psc.edu/general/software/packages/ieee/ieee.php
SMALL_NUMBER = 2**-149

class FitsShrinker(object):

    def __init__(self, diffuse_model, skydir, radius):

        self.diffuse_model = diffuse_model
        self.skydir = skydir
        self.radius = radius

        self.test_good_file()


    def test_good_file(self):

        diffuse_model = self.diffuse_model
        assert isinstance(diffuse_model, pyfits.core.HDUList)

        h = diffuse_model['PRIMARY'].header
        # make sure 2 or 3 dimesional data
        assert h['NAXIS'] in [2, 3]

        # Here, make sure diffuse function is allsky
        assert abs(abs(h['NAXIS1']*h['CDELT1'])-360) < 1
        assert abs(abs(h['NAXIS2']*h['CDELT2']) -180) < 1

        # make sure diffuse file has same binning
        assert abs(h['CDELT1']) == abs(h['CDELT2'])

    @staticmethod
    def rad_cut(lons,lats,galactic,skydir, radius):
        """ Return a boolean array that is true when the distance from the 
            skydir is less than the radius (in degrees). """
        if galactic:
            lon0,lat0 = skydir.l(),skydir.b()
        else:
            lon0,lat0 = skydir.ra(),skydir.dec()

        lon0,lat0=np.radians(lon0),np.radians(lat0)
        lons,lats = np.radians(lons),np.radians(lats)

        cos_diffs = np.sin(lats)*np.sin(lat0)+np.cos(lats)*np.cos(lat0)*np.cos(lons-lon0)

        # remember, increasing angle, decreasing cos(angle)
        mask = -cos_diffs < -np.cos(np.radians(radius))
        return mask


    @staticmethod
    def pyfits_to_sky(data, header):
        """ Return two arrays with the same dimensions as data but
            containing either RA & Dec or L & B for the data. 
        """
        wcs = pywcs.WCS(header)

        dataT = data.transpose()
        if len(data.shape) == 3 and header['naxis'] == 3:
            xlen,ylen,zlen = dataT.shape
            x,y,z = np.mgrid[0:xlen, 0:ylen, 0:zlen]

            length = len(x.flatten())
            indices = np.hstack((x.reshape((length,1)),
                                 y.reshape((length,1)),
                                 z.reshape((length,1))))

            lons,lats,energy=wcs.all_pix2sky(indices, 1).transpose()

            lons=lons.reshape(dataT.shape).transpose()
            lats=lats.reshape(dataT.shape).transpose()
            energy=energy.reshape(dataT.shape).transpose()

            if data.shape[0] == 1:
                # Sometimes the energy axis is empty, in which case
                # dont' return energies
                return lons,lats
            else:
                return lons,lats,energy

        elif len(data.shape) == 2 and header['naxis'] == 2:

            dataT = data.transpose()

            xlen,ylen = dataT.shape
            x,y = np.mgrid[0:xlen, 0:ylen]

            length = len(x.flatten())
            indices = np.hstack((x.reshape((length,1)),
                                 y.reshape((length,1))))

            lons,lats=wcs.all_pix2sky(indices, 1).transpose()

            lons=lons.reshape(dataT.shape).transpose()
            lats=lats.reshape(dataT.shape).transpose()

            return lons,lats
        else:
            raise Exception("Algorithm only defiend for length 2 or 3 arrays.")


    @staticmethod
    def smaller_range(values, min, max):
        """ values should be a sorted list of numbers.

            This function will return a boolean array of the values
            between min and max, but the list will be inclusive so that it will
            include the first number smaller than min and the first number
            larger than max.

                >>> print FitsShrinker.smaller_range(np.asarray([1,2,3,5]), 2, 3)
                [False  True  True False]

                >>> print FitsShrinker.smaller_range(np.asarray([1,2,3,5]), 1.9, 3.1)
                [ True  True  True  True]

            If there is no number smaller than min (or larger than max), it just includes
            the whole list

                >>> print FitsShrinker.smaller_range(np.asarray([1,2,3,5]), 0.9, 5.1)
                [ True  True  True  True]
        """
        good = np.ones_like(values).astype(bool)
        g = min < values
        if sum(~g) > 1:
            index = np.where(~g)[0][:-1]
            good[index] = False
        g = values < max
        if sum(~g) > 1:
            index = np.where(~g)[0][1:]
            good[index] = False
        return good


class MapShrinker(FitsShrinker):

    def shrink(self, **kwargs):

        x = self.diffuse_model

        lons,lats = FitsShrinker.pyfits_to_sky(x['PRIMARY'].data, x['PRIMARY'].header)

        galactic=True
        inside_cut = FitsShrinker.rad_cut(lons,lats,galactic,self.skydir, self.radius)

        x['PRIMARY'].data[~inside_cut] = SMALL_NUMBER


class DiffuseShrinker(FitsShrinker):
    def __init__(self, diffuse_model, skydir, radius, emin, emax):
        """ Take in an allsky diffuse model and create a diffuse model which predicts
            emission only inside the given radius and energy range. 
            

            A simple test that creates a shrunk file and shows that the pixes are 0 far away
            from the center
            
                >>> from os.path import expandvars
                >>> galname=expandvars('$GLAST_EXT/diffuseModels/v0r0/gll_iem_v02.fit')
                >>> diff=pyfits.open(galname)
                >>> shrinker = DiffuseShrinker(diff, skydir=SkyDir(55,60, SkyDir.GALACTIC), radius=10, emin=1e3, emax=1e5)
                >>> shrinker.shrink()
                >>> from tempfile import NamedTemporaryFile
                >>> tempfile = NamedTemporaryFile()
                >>> diff.writeto(tempfile.name)
                >>> cut = DiffuseFunction(tempfile.name)
                >>> uncut = DiffuseFunction(galname)

            First, test at the skydir:

                >>> sd = SkyDir(55,60, SkyDir.GALACTIC)
                >>> np.allclose(cut(sd, 1e4), uncut(sd, 1e4))
                True

            Note, we appleid a 1 GeV energy cut, so in teh cut file we cannot access very small eneriges:

                >>> not np.allclose(uncut(sd, 1e2),0)
                True
                >>> cut(sd, 1e2)
                Traceback (most recent call last):
                    ...
                RuntimeError: Diffuse function: energy out of range: 100

            Then test at a point far away, it should be 0:

                >>> sd = SkyDir(100,60, SkyDir.GALACTIC)
                >>> np.allclose(cut(sd, 1e4), 0)
                True
        """

        super(DiffuseShrinker,self).__init__(diffuse_model, skydir, radius)
        self.emin = emin
        self.emax = emax

    def shrink(self):

        d = self.diffuse_model

        self.all_energies = np.asarray(d['energies'].data.field('energy'))

        # create a list of the indices in the self.all_energy array for energies
        # that are allowed by the energy cut.
        self.good_indices = FitsShrinker.smaller_range(self.all_energies, self.emin, self.emax)
        self.good_energies = self.all_energies[self.good_indices]

        layers = sum(self.good_indices)

        # first, strip out unwanted energies
        d['PRIMARY'].data = data = d['PRIMARY'].data[self.good_indices]
        d['ENERGIES'].data = d['ENERGIES'].data[self.good_indices]

        d['ENERGIES'].header['NAXIS2'] = layers

        header = d['PRIMARY'].header
        header['NAXIS3'] = layers
        header['CRVAL3'] = d['ENERGIES'].data.field('ENERGY')[0]

        # Create a 2 dimensinoal header from the 3d header
        header_2d = header.copy()
        header_2d['NAXIS'] = 2
        for i in ['NAXIS3', 'CRVAL3', 'CDELT3', 'CRPIX3', 'CTYPE3', 'CUNIT3']:
            del header_2d[i]

        for i in range(data.shape[0]):
            # Cut each energy level separately, less memory intensive
            # for large galactic diffuse models
            lons,lats = FitsShrinker.pyfits_to_sky(data[i], header_2d)

            galactic=True
            inside_cut = FitsShrinker.rad_cut(lons,lats,galactic,self.skydir, self.radius)

            d['PRIMARY'].data[i][~inside_cut] = SMALL_NUMBER


class NoSimulatedPhotons(Exception):
    pass

class MonteCarlo(object):
    """ This object allowes for the simulation of data. Simply pass in
        a list of point sources and a list of diffuse sources and the
        object will simulate the region.

        The simulated data will get stored in the file 'ft1'. If you pass
        in a valid ft2 file, the data will be simulted using that pointing
        history. Otherwise, a default rocking profile will be used.  """

    defaults = (
            ('seed',               0, " Random number generator seen when running obssim. This should be varied for simulations."),
            ('savedir',         None, " If specified, save temporary files to this directory."),
            ('tempbase', '/scratch/', " Directory to put temporary files in. Can be cleaned up with the cleanup function."),
            ('tstart',          None, " Required if ft2 does not point to a real file."),
            ('tstop',           None, " Required if ft2 does not point to a real file."),
            ('ft2',             None, " If exists, simulate with this ft2 file. If not, create."),
            ('ltfrac',          None, " If set, pass value into gtobssim."),
            ('emin',             100, " Minimum energy"),
            ('emax',          100000, " Maximum energy"),
            ('conv_type',         -1, " Conversion Type"),
            ('roi_dir',         None, " Center of ROI. Gtobssim will use the use_ac flag if this is specified."),
            ('maxROI',          None, " Maximum ROI Size. Gtobssim will use the use_ac flag if this is specified."),
            ('quiet',          False, " Surpress output."),
            ('mc_energy',      False, " Use MC Energy"),
            ('gtifile',         None, " Make the simulated data have only the same GTIs as the GTIs in this file (typically an ft1 file or ltcube)."),
            ('diffuse_pad',     10.0,  " How many area outside ROI should the diffuse emission be simulated to account for energy dispersion."),
            ('roi_pad',          5.0,  " Include counts from this far ouside maxROI in the ft1 file to account for ragged edge effect in pointlike."),
            ('energy_pad',       2.0,  """ Lower energy of simluated photons is emin/energy_bad and 
                                          upper energy of simulated photos is energy_pad*emax.
                                          This allows for energy dispersion effects to be 
                                          naturally accounted for in the monte carlos simulation. """),
    )

    @staticmethod
    def strip(name):
        """ Create source names that can be used in files and gtobssim's xml.
            strip out periods and parenthesis.
            Replace sapces with underbars, pluses and minumes with p & m,
            and put an underbar in the front of the name if it begins
            with a digit. """
        name = name.replace(' ', '_').replace('+','p').replace('-','m')
        name = re.sub('[\.()]','',name)
        if name[0].isdigit(): name = '_' + name
        return name

    @decorate(defaults)
    def __init__(self,ft1,irf, sources, **kwargs):
        """ Constructor does not require a data_specification. """
        process(self, kwargs)

        self.ft1=ft1
        self.irf=irf
        self.sources=sources

        if self.savedir is not None:
            self.savedir=os.path.abspath(self.savedir)
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            self.save_output=True
        else:
            if not os.path.isdir(self.tempbase):
                raise Exception("Argument tempbase must point to a directory where temporary files can be placed.")
            self.savedir=mkdtemp(prefix=self.tempbase)
            self.save_output=False


        if not isinstance(self.ft1,types.StringType):
            if len(self.ft1) != 1: raise Exception(dedent("""\
                                                          Exactly one ft1 file may be specified by roi_monte_carlo script
                                                          (%s were acutally specified)""" % len(self.ft1)))
            self.ft1=self.ft1[0]

        if not isinstance(self.ft2,types.StringType):
            if len(self.ft2) != 1: raise Exception(dedent("""\
                                                   Exactly one ft2 file may be specified by roi_monte_carlo script
                                                   (%s were acutally specified)""" % len(self.ft2)))
            self.ft2=self.ft2[0]

        if os.path.exists(self.ft1):
            raise Exception("Unable to run MC simulation because file %s already exists." % self.ft1)

        if self.ft2 is not None and os.path.exists(self.ft2):
            self.use_existing_ft2 = True
        else:
            self.use_existing_ft2 = False

        if self.use_existing_ft2:
            if self.tstart is not None or self.tstop is not None:
                raise Exception("Since simulating from exisitng ft2 file, tstart and tstop can not be set.")
        else:
            if self.tstart is None or self.tstop is None:
                raise Exception("tstart and tstop must be specified since ft2 does not exist.")

            if self.tstop - self.tstart < 1: raise Exception("tstart and tstop must describe a real range.")

        if self.use_existing_ft2: self.set_time_from_ft2()

        if self.use_existing_ft2:
            if not self.gtifile: 
                print "Since an already existing ft2 file is set, it is strongly recommended that you specify a file with the desired GTIs (typically the ft1 or ltcube file."
        else:
            if self.gtifile:
                print "gtifile can only be set for existing ft2 files."

    def set_time_from_ft2(self):
        # Note, get the start & stop times from the actual
        # data instead of the fits header
        # See https://jira.slac.stanford.edu/browse/OBS-18
        ft2 = pyfits.open(self.ft2)
        self.tstart = ft2['SC_DATA'].data.field('START')[0]
        self.tstop = ft2['SC_DATA'].data.field('STOP')[-1]
        ft2.close()

    def larger_energy_range(self):
        """ Get an energy range larger then the desired simulated points
            (to correct for energy dispersion). """
        return self.emin/self.energy_pad,self.energy_pad*self.emax

    def _make_ps(self,ps,mc_emin,mc_emax,indent):

        assert isinstance(ps,PointSource)

        if not self.quiet: print 'Processing source %s' % ps.name

        model = ps.model

        if isinstance(model,PowerLaw) or isinstance(model,PowerLawFlux):
            xml=[
                '<source name="%s" flux="%s">' % (MonteCarlo.strip(ps.name),model.i_flux(mc_emin,mc_emax,cgs=True)*1e4),
                '  <spectrum escale="MeV">',
                '  <particle name="gamma">',
                '    <power_law emin="%s" emax="%s" gamma="%s"/>' % (mc_emin,mc_emax,model.getp(1)),
                '  </particle>',
                '  <celestial_dir ra="%s" dec="%s"/>' % (ps.skydir.ra(),ps.skydir.dec()),
                '  </spectrum>',
                '</source>'
            ]
            return indent+('\n'+indent).join(xml)
        else:
            if isinstance(model,FileFunction):
                flux=model.i_flux(self.energy[0],self.energy[-1],cgs=True)*1e4
                spectral_filename=model.file
            else:
                flux=model.i_flux(mc_emin,mc_emax,cgs=True)*1e4

                spectral_filename = '%s_spectra_%s.txt' % (MonteCarlo.strip(ps.name),model.name)
                model.save_profile(filename=spectral_filename, emin=mc_emin, emax=mc_emax)

            xml=[
                '<source name="%s">' % MonteCarlo.strip(ps.name),
                '  <spectrum escale="MeV">',
                '    <SpectrumClass name="FileSpectrum"',
                '      params="flux=%g,specFile=%s"/>' % (flux,spectral_filename),
                '    <celestial_dir ra="%s" dec="%s"/>' % (ps.skydir.ra(),ps.skydir.dec()),
                '  </spectrum>',
                '</source>'
            ]
            return indent+('\n'+indent).join(xml)


    @staticmethod
    def _make_profile(name,spatial_model,numpoints=200):
        temp='%s_extension_profile_%s.txt' % (MonteCarlo.strip(name),spatial_model.name)
        radius,pdf = spatial_model.approximate_profile()
        open(temp,'w').write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(radius,pdf)]))
        return temp

    def _make_radially_symmetric(self,es,mc_emin, mc_emax, indent):

        name=es.name
        sm=es.spatial_model
        model = es.model

        assert isinstance(sm,RadiallySymmetricModel)
 
        flux=model.i_flux(mc_emin,mc_emax,cgs=True)*1e4

        ra,dec=sm.center.ra(),sm.center.dec()

        spatial_filename='%s_extension_profile_%s.txt' % (MonteCarlo.strip(name),sm.name)
        sm.save_profile(spatial_filename)

        spectral_filename = '%s_spectra_%s.txt' % (MonteCarlo.strip(es.name),model.name)
        model.save_profile(filename=spectral_filename, emin=mc_emin, emax=mc_emax)

        xml=[
            '<source name="%s">' % MonteCarlo.strip(es.name),
            '  <spectrum escale="MeV">',
            '    <SpectrumClass name="RadialSource"',
            '      params="flux=%s, profileFile=%s, specFile=%s, ra=%s, dec=%s"/>' % \
                    (flux,spatial_filename,spectral_filename,ra,dec),
            '    <use_spectrum frame="galaxy"/>',
            '  </spectrum>',
            '</source>',
        ]
        return indent+('\n'+indent).join(xml)

    def _make_extended_source(self,es,mc_emin,mc_emax,indent):

        sm=es.spatial_model
        model = es.model

        assert isinstance(sm,SpatialModel)

        flux=model.i_flux(mc_emin,mc_emax,cgs=True)*1e4

        ra,dec=sm.center.ra(),sm.center.dec()

        if isinstance(sm,SpatialMap):
            print 'WARNING: gtobssim can only use plate-carree projection fits files!'
            spatial_filename=sm.file
        else:
            spatial_filename='%s_spatial_template_%s.fits' % (MonteCarlo.strip(es.name),sm.name)
            # Allegedly simulated templates must only be in the plate-carree projection
            # http://www.slac.stanford.edu/exp/glast/wb/prod/pages/sciTools_observationSimTutorial/obsSimTutorial.htm
            sm.save_template(spatial_filename, proj='CAR')

        if isinstance(model,PowerLaw) or isinstance(model,PowerLawFlux):

            index = model['index']

            xml=[
                '<source name="%s">' % MonteCarlo.strip(es.name),
                '  <spectrum escale="MeV">',
                '    <SpectrumClass name="MapSource"',
                '      params="%s,%s,%s,%s,%s"/>' % \
                        (flux,index,spatial_filename,mc_emin,mc_emax),
                '    <use_spectrum frame="galaxy"/>',
                '  </spectrum>',
                '</source>',
            ]

        else:

            spectral_filename='%s_spectra_%s.txt' % (MonteCarlo.strip(es.name),model.name)
            model.save_profile(filename=spectral_filename, emin=mc_emin, emax=mc_emax)

            xml=[
                '<source name="%s">' % MonteCarlo.strip(es.name),
                '  <spectrum escale="MeV">',
                '    <SpectrumClass name="FileSpectrumMap"',
                '      params="flux=%s, fitsFile=%s, specFile=%s, emin=%s, emax=%s"/>' % \
                        (flux,spatial_filename,spectral_filename,mc_emin,mc_emax),
                '    <use_spectrum frame="galaxy"/>',
                '  </spectrum>',
                '</source>',
            ]
        return indent+('\n'+indent).join(xml)

    def _make_powerlaw_gaussian(self,es,mc_emin,mc_emax,indent):
        """ Make an extended source. """
        model = es.model

        assert isinstance(es.spatial_model,Gaussian) or \
               isinstance(es.spatial_model,EllipticalGaussian)

        if isinstance(model,PowerLaw) or isinstance(model,PowerLawFlux):

            flux=model.i_flux(mc_emin,mc_emax,cgs=True)*1e4
            gamma=model.getp('Index')
            center=es.spatial_model.center
            ra,dec=center.ra(),center.dec()
            

            if isinstance(es.spatial_model,Gaussian):
                major=es.spatial_model.get_parameters(absolute=True)[2]
                minor=major
                angle=0
            elif isinstance(es.spatial_model,EllipticalGaussian):
                major,minor,angle=es.spatial_model.get_parameters(absolute=True)[2:]
                angle *= -1 # For some reason, gtobssim measures the angle anti-east of north...

            xml=[
                '<source name="%s">' % MonteCarlo.strip(es.name),
                '   <spectrum escale="MeV">',
                '      <SpectrumClass name="GaussianSource"',
                '          params="%s, %s, %s, %s, %s, %s, %s, %s, %s"/>' % \
                                        (flux,gamma,ra,dec,major,minor,angle,mc_emin,mc_emax),
                '      <use_spectrum frame="galaxy"/>',
                '   </spectrum>',
                '</source>'
            ]
            return indent+('\n'+indent).join(xml)

        else:
            raise Exception("Can only parse PowerLaw gaussian sources.")

    @staticmethod
    def make_isotropic_fits(radius, skydir=None, **kwargs):
        """ Note, if there is an ROI cut, we can make this
            isotropic file not allsky. """

        allsky = MonteCarlo.make_allsky_isotropic_fits(**kwargs)
        if radius < 180:
            shrinker = MapShrinker(allsky, skydir, radius)
            shrinker.shrink()
        return allsky

    @staticmethod
    def make_allsky_isotropic_fits(pixelsize=1):
        """ Create an allsky pyfits file with 1s in it.  """

        assert 180 % pixelsize == 0
        import pywcs

        wcs = pywcs.WCS(naxis=2)

        nxpix = 360/pixelsize
        nypix = 180/pixelsize

        wcs.wcs.crpix = [nxpix/2.0 + 0.5, nypix/2.0 + 0.5]
        wcs.wcs.cdelt = [-pixelsize, pixelsize]
        wcs.wcs.crval = [0,0]
        wcs.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]

        header = wcs.to_header()

        data = np.ones((nypix,nxpix), dtype=float)
        hdu = pyfits.PrimaryHDU(header=header, data=data)

        hdulist = pyfits.HDUList([hdu])
        return hdulist


    @staticmethod
    def isotropic_spectral_integrator(filename):
        """ Assume a powerlaw connecting each point. 

            This code is adopated from the function pl_integral in 
            genericSources/src/FileSpectrum.cxx which is used by gtobssim
            to integrate the isotropic spectrum.
        
            Returns the flux differential in solid angle: ph/m^2/sr 
            
            N.B. multiply by 10^4 to convert from ph/cm^2/sr to ph/m^2/sr
            """
        file=np.genfromtxt(filename,unpack=True)
        energy,flux=file[0],file[1]
        emin,emax = energy[0], energy[-1]

        e1,e2=energy[:-1],energy[1:]
        f1,f2=flux[:-1],flux[1:]

        gamma=np.log(f2/f1)/np.log(e2/e1)

        n0 = f1/e1**gamma

        flux=np.sum(
            np.where(gamma != -1,
                    n0/(gamma+1)*(e2**(gamma+1)-e1**(gamma+1)),
                    n0*np.log(e2/e1)
            )
        )*10**4
        return flux,emin,emax

    @staticmethod
    def spatial_integrator_2d(filename):
        """ Computes the spatial integral of a 2D fits file.


                >>> from tempfile import NamedTemporaryFile
                >>> from skymaps import SkyDir
                >>> import numpy as np

            If the file is allsky and has inside of it the value '1', then the integral should return
            the total solid angle of a sphere (4*pi).

                >>> allsky=MonteCarlo.make_isotropic_fits(radius=180, pixelsize=1)
                >>> tempfile = NamedTemporaryFile()
                >>> allsky.writeto(tempfile.name, clobber=True)
                >>> map_area = MonteCarlo.spatial_integrator_2d(tempfile.name)
                >>> np.allclose(map_area, 4*np.pi)
                True

            When the fits file contains the value '1' in it only within a given radius, the map 
            integral should should be 4*pi*(1-cos(radius)). Note, we need smaller
            pixel size to get a good agreement due to ragged edge effect at corner of circle.

                >>> def map_area(radius):
                ...     tempfile = NamedTemporaryFile()
                ...     allsky=MonteCarlo.make_isotropic_fits(skydir=SkyDir(134,83,SkyDir.GALACTIC), radius=radius, pixelsize=0.25)
                ...     allsky.writeto(tempfile.name, clobber=True)
                ...     map_area = MonteCarlo.spatial_integrator_2d(tempfile.name)
                ...     return map_area
                >>> map_area = np.vectorize(map_area)
                >>> true_area = lambda radius: 2*np.pi*(1-np.cos(np.radians(radius)))
                >>> np.allclose(map_area([10,30,50]), true_area([10,30,50]), atol=1e-3, rtol=1e-3)
                True
        """
        fits=pyfits.open(filename)
        data=fits[0].data
        assert len(data.shape) == 2 or data.shape[0] == 1
        if data.shape == 3: data = data[1]

        # get the solid angles from gtlike. Probably
        # could be rewritten in pointlike, but not necessary
        mapcube=pyLikelihood.MapCubeFunction2(filename)
        sa=np.asarray(mapcube.wcsmap().solidAngles()).transpose()
        map_integral = np.sum(sa*data)
        fits.close()
        return map_integral

    @staticmethod
    def isone(model):
        """ Return 1 if model predicts 1 everywhere. """
        if isinstance(model,Constant) and model['scale'] == 1:
            return 1
        if isinstance(model,PowerLaw) and model['norm'] == 1 and model['index'] == 1 and \
           hasattr(model,'index_offset') and model.index_offset == 1:
            return 1
        return 0

    def _make_isotropic_diffuse(self,ds,mc_emin,mc_emax,indent):

        dm=ds.dmodel[0]
        sm=ds.smodel

        if not MonteCarlo.isone(sm):
            raise Exception("Can only run gtobssim with IsotropicSpectrum diffuse models if model predicts 1.")

        spectral_file=dm.name()
        energies,spectra=np.genfromtxt(spectral_file,unpack=True)[0:2]

        smaller_range = FitsShrinker.smaller_range(energies, mc_emin, mc_emax)
        if np.any(smaller_range == False):
            # If any energies are outside the range, cut down spectral file.
            cut_spectral_file = os.path.basename(spectral_file).replace('.txt','_cut.txt')
            np.savetxt(cut_spectral_file, zip(energies[smaller_range], spectra[smaller_range]))
        else:
            cut_spectral_file = spectral_file

        spatial_file=os.path.basename(spectral_file).replace('.txt','_spatial.fits')
        if self.roi_dir is not None and self.maxROI is not None:
            radius=self.maxROI + self.diffuse_pad
        else:
            radius=180

        if not self.quiet: print '.. Making isotropic model for %s' % ds.name
        allsky=MonteCarlo.make_isotropic_fits(skydir=self.roi_dir, radius=radius)
        allsky.writeto(spatial_file)


        # flux is ph/cm^2/sr to ph/m^2
        flux, emin, emax = MonteCarlo.isotropic_spectral_integrator(cut_spectral_file)

        if not self.quiet: print '.. Integrating isotropic model for %s' % ds.name

        # multiply by solid angle to convert to ph/m^2
        flux*=MonteCarlo.spatial_integrator_2d(spatial_file)


        ds = [ 
            '<source name="%s">' % MonteCarlo.strip(ds.name),
            '  <spectrum escale="MeV">',
            '    <SpectrumClass name="FileSpectrumMap"',
            '       params="flux=%s,fitsFile=%s,' % (flux,spatial_file),
            '       specFile=%s,emin=%s,emax=%s"/>' % (cut_spectral_file,emin,emax),
            '    <use_spectrum frame="galaxy"/>',
            '  </spectrum>',
            '</source>'
        ]
        return indent+('\n'+indent).join(ds)

    def _make_powerlaw_diffuse(self,ds,mc_emin,mc_emax,indent):
        dm=ds.dmodel[0]
        sm=ds.smodel

        if not MonteCarlo.isone(sm):
            raise Exception("Can only run gtobssim with DiffuseFunction diffuse models if model predicts 1.")

        # flux in ph/cm^2/s/sr b/n 100MeV & infinity
        flux_pointlike=dm.flux()
        index=dm.index()

        x=PowerLawFlux(p=[flux_pointlike,index],emin=100,emax=np.inf)

        if self.roi_dir is not None and self.maxROI is not None:
            ra,dec=self.roi_dir.ra(),self.roi_dir.dec()
            # just to be safe from the large low energy psf
            radius=self.maxROI + self.diffuse_pad

            # gtobssim wants ph/m^2/s, integrate over solid angle
            flux=x.i_flux(mc_emin,mc_emax)*(2*np.pi*(1-np.cos(np.radians(radius))))*10**4
        else:
            ra,dec,radius=0,0,180

            # gtobssim wants ph/m^2/s
            flux=x.i_flux(mc_emin,mc_emax)*4*np.pi*10**4


        ds = [
            '<source name="%s">' % MonteCarlo.strip(ds.name),
            '   <spectrum escale="MeV">',
            '      <SpectrumClass name="Isotropic"',
            '                     params="flux=%g,gamma=%g,emin=%g,emax=%g,ra=%g,dec=%g,radius=%g"/>' % (flux,index,mc_emin,mc_emax,ra,dec,radius),
            '      <use_spectrum frame="galaxy"/>',
            '   </spectrum>',
            '</source>',
        ]
        return indent+('\n'+indent).join(ds)


    @staticmethod
    def diffuse_integrator(filename):
        """ A simpel test, this file $GLAST_EXT/diffuseModels/v1r0/obssim_v02_P6_V11_DIFFUSE.xml
            claims that the normaliztiano of gll_iem_v02_P6_V11_DIFFUSE.fit
            should be 8.423. The value we find is a little bit different, but
            not by very much (8.409):

                >>> from os.path import expandvars
                >>> filename = expandvars('$GLAST_EXT/diffuseModels/v1r0/gll_iem_v02_P6_V11_DIFFUSE.fit')
                >>> np.allclose(MonteCarlo.diffuse_integrator(filename), 8.423, rtol=1e-2)
                True
        """
        flux = MonteCarlo.mapIntegral(filename)*10**4
        return flux

    @staticmethod
    def mapIntegral(filename):
        """ No idea why, but the function was implemented directly
            in MapCubeFunction but broke in MapCubeFunction2. Now there
            is no swig interface to MapCubeFunction.

            This code is shamelessly stolen from MapCubeFunction::mapIntegral
            in MapCubeFunction.cxx

            http://www-glast.stanford.edu/cgi-bin/viewcvs/Likelihood/src/MapCubeFunction.cxx

            This function is an ugly mess and could probably be properly rewritten
            as a fully vecotrized numpy function, but this is good enough, for now.

        """

        fits=pyfits.open(filename)
        energies=fits[1].data.field('Energy')
        data=fits[0].data

        nxpix,nypix=data.shape[1],data.shape[2]

        mapcube=pyLikelihood.MapCubeFunction2(filename)
        sa=np.asarray(mapcube.wcsmap().solidAngles()).transpose()

        # First, integrate the spatial part
        d = (sa*data).sum(axis=2).sum(axis=1)

        # now, free up the otherwise very costly memory
        del(data)
        del(mapcube)
        del(sa)
        fits.close()

        # Then integrate the spectral part connecting each point
        # with a powerlaw
        p = MonteCarlo.powerLawIntegral
        map_integral = np.sum(p(energies[0:-1],energies[1:],
                                d[0:-1],d[1:]).sum(axis=0))

        return map_integral

    @staticmethod
    def powerLawIntegral(x1, x2, y1, y2):
        """ This function is blatently stolen from
            MapCubeFunction::powerLawIntegral
            in MapCubeFunction.cxx:
                http://www-glast.stanford.edu/cgi-bin/viewcvs/Likelihood/src/MapCubeFunction.cxx
        """
        gamma = np.log(y2/y1)/np.log(x2/x1)
        n0 = y1/x1**gamma
        gp1 = gamma + 1.;
        integral = np.where(gamma != 1., n0/gp1*(x2**gp1 - x1**gp1), n0*np.log(x2/x1))
        return integral

    def _make_isotropic_constant(self,ds,mc_emin,mc_emax,indent):

        dm=ds.dmodel[0]
        sm=ds.smodel

        if dm.constant() != 1:
            raise Exception("When simulation IsotropicConstant source, the constant must be 1")

        if isinstance(sm,PowerLaw) or isinstance(sm,PowerLawFlux):
            index=sm['index']

            if self.roi_dir is not None and self.maxROI is not None:
                ra,dec=self.roi_dir.ra(),self.roi_dir.dec()
                # just to be safe from the large low energy psf
                radius=self.maxROI + self.diffuse_pad

                # gtobssim wants ph/m^2/s, integrate over solid angle
                flux=sm.i_flux(mc_emin,mc_emax)*(2*np.pi*(1-np.cos(np.radians(radius))))*10**4
            else:
                ra,dec,radius=0,0,180

                # gtobssim wants ph/m^2/s
                flux=sm.i_flux(mc_emin,mc_emax)*4*np.pi*10**4

            ds = [
                '<source name="%s">' % MonteCarlo.strip(ds.name),
                '   <spectrum escale="MeV">',
                '      <SpectrumClass name="Isotropic"',
                '                     params="flux=%g,gamma=%g,emin=%g,emax=%g,ra=%g,dec=%g,radius=%g"/>' % (flux,index,mc_emin,mc_emax,ra,dec,radius),
                '      <use_spectrum frame="galaxy"/>',
                '   </spectrum>',
                '</source>',
            ]
        else:
            raise Exception("Unable to create XML for source %s" % ds.name)

        return indent+('\n'+indent).join(ds)





    def _make_diffuse(self,ds,mc_emin,mc_emax,indent):

        dm=ds.dmodel[0]
        sm=ds.smodel

        # galactic diffuse
        if not MonteCarlo.isone(sm):
            raise Exception("Can only run gtobssim with DiffuseFunction diffuse models where the spectral model is a PowerLaw with norm and index 1.")

        allsky_filename=dm.name()

        if self.roi_dir is not None and self.maxROI is not None:
            print '.. Shrinking diffuse model %s' % ds.name
            radius=self.maxROI + self.diffuse_pad

            allsky = pyfits.open(allsky_filename)
            shrinker = DiffuseShrinker(allsky, skydir=self.roi_dir, radius=radius, emin=mc_emin, emax=mc_emax)
            shrinker.shrink()
            filename = os.path.basename(allsky_filename).replace('.fits','_cut.fits')
            allsky.writeto(filename)


        else:
            filename = allsky_filename


        print '.. Integrating diffuse model %s' % ds.name
        flux = MonteCarlo.diffuse_integrator(filename)


        ds = [
            '<source name="%s">' % MonteCarlo.strip(ds.name),
            '   <spectrum escale="MeV">',
            '      <SpectrumClass name="MapCube" params="%s,%s"/>' % (flux,filename),
            '      <use_spectrum frame="galaxy"/>',
            '   </spectrum>',
            '</source>'
        ]
        return indent+('\n'+indent).join(ds)

    def _make_ds(self,ds,*args,**kwargs):
        if not self.quiet: print 'Processing source %s' % ds.name

        if len(ds.dmodel) > 1:
            raise Exception("Can only run gtobssim for diffuse models with a combined front/back diffuse model.")

        if isinstance(ds.dmodel[0],IsotropicConstant):
            return self._make_isotropic_constant(ds,*args,**kwargs)

        elif isinstance(ds.dmodel[0],IsotropicSpectrum):
            return self._make_isotropic_diffuse(ds,*args,**kwargs)

        elif isinstance(ds.dmodel[0],IsotropicPowerLaw):
            return self._make_powerlaw_diffuse(ds,*args,**kwargs)

        elif isinstance(ds.dmodel[0],DiffuseFunction):
            return self._make_diffuse(ds,*args,**kwargs)

        elif isinstance(ds,ExtendedSource):
            if (isinstance(ds.spatial_model,Gaussian) or \
                    isinstance(ds.spatial_model,EllipticalGaussian)) and \
                    (isinstance(ds.model,PowerLaw) or isinstance(ds.model,PowerLawFlux)):
                return self._make_powerlaw_gaussian(ds,*args,**kwargs)
            elif isinstance(ds.spatial_model,RadiallySymmetricModel):
                return self._make_radially_symmetric(ds,*args,**kwargs)
            else:
                return self._make_extended_source(ds,*args,**kwargs)
        else:
            raise Exception("Can not simulate diffuse source %s. Unknown diffuse source type %s." % (ds.name,type(ds.dmodel[0])))
            
    def _make_model(self,indent='  '):

        xml = ['<source_library title="Library">']
        src = []

        self.point_sources = []
        self.diffuse_sources = []

        if not isinstance(self.sources,collections.Iterable): 
            self.sources=[self.sources]

        for source in self.sources:
            assert isinstance(source,PointSource) or isinstance(source,DiffuseSource)
            if isinstance(source,PointSource):
                self.point_sources.append(source)
            elif isinstance(source,DiffuseSource):
                self.diffuse_sources.append(source)

        mc_emin,mc_emax=self.larger_energy_range()

        for ps in self.point_sources:
            xml.append(self._make_ps(ps,mc_emin, mc_emax,indent))
            src.append(MonteCarlo.strip(ps.name))

        for ds in self.diffuse_sources:
            xml.append(self._make_ds(ds,mc_emin,mc_emax,indent))
            src.append(MonteCarlo.strip(ds.name))

        xml.append('</source_library>')
    
        self.xmlfile=os.path.join(self.savedir,'source_library.xml')
        temp=open(self.xmlfile,'w')
        temp.write('\n'.join(xml))
        temp.close()

        self.srclist=os.path.join(self.savedir,'source_list.txt')
        temp=open(self.srclist,'w')
        temp.write('\n'.join(src))
        temp.close()

    def gtmktime_from_file(self, evfile, scfile, outfile, gtifile):
        """ Apply the GTIs from gtifile to the file evfile. """

        if not self.quiet: print 'Applying GTIs from file %s to simulated data' % gtifile

        # Note, add on gtis to 'evfile'. This is a big distructive,
        # but should cause no real harm.
        e = pyfits.open(evfile, mode='update')

        # Temporarily address issue https://jira.slac.stanford.edu/browse/OBS-20
        if len(e) == 4 and e[2].name == 'GTI' and e[3].name == 'GTI' and e[3].data == None:
            del(e[3])
        
        g = pyfits.open(gtifile)
        e['GTI'] = g['GTI']
        e.flush()

        app=GtApp('gtmktime')
        app.run(evfile=evfile,
                outfile=outfile,
                scfile=scfile,
                # Here, apply no additional gtmktime cuts
                filter='T', 
                roicut='no',
                # This will cut on previously added GTI cuts
                apply_filter='yes',
               )

    def simulate(self, dry_run=False):
        """ dry_run: make all the tempfiles and print out gtobssim command, but 
                     don't acutally run it.
        """

        if not self.quiet: print 'Simulating in energy range from %g MeV to %g MeV' % (self.emin, self.emax)

        old_dir=os.getcwd()

        if not self.quiet: print 'working in directory',self.savedir
        os.chdir(self.savedir)
        self._make_model()

        irfs=self.irf

        if self.conv_type == 0:
            irfs+="::FRONT"
        elif self.conv_type == 1:
            irfs+="::BACK"

        if self.roi_dir is not None and self.maxROI is not None:
            use_ac="yes"
            ra=self.roi_dir.ra()
            dec=self.roi_dir.dec()
            # Add an extra 5 degrees due to ragged edge healpix effect in pointlike
            radius=self.maxROI + self.roi_pad
        else:
            use_ac="no"
            ra,dec,radius=0,0,180
        
        os.environ['PFILES']=self.savedir+';'+os.environ['PFILES'].split(';')[-1]; # to set the writing pfiles to the savedir
        if not self.quiet: print 'current pfile settings',os.getenv('PFILES'); #should add the output for debugging reasons
        app=GtApp('gtobssim');
        if self.ltfrac is not None: app['ltfrac']=self.ltfrac
        app.run(infile=self.xmlfile,
                srclist=self.srclist,
                scfile=self.ft2 if self.ft2 is not None and os.path.exists(self.ft2) else 'none',
                evroot='sim',
                startdate="2001-01-01 00:00:00",
                tstart=self.tstart,
                simtime=self.tstop-self.tstart,
                use_ac=use_ac, ra=ra, dec=dec, radius=radius,
                irfs=irfs,
                edisp='no' if self.mc_energy else 'yes',
                seed=self.seed,
                emin=self.emin,
                emax=self.emax,
                maxrows=200000000, # make sure all photons are in one file
                dry_run=dry_run, # to pass all keywords from simulate through to app.run();
       )

        os.chdir(old_dir)

        if dry_run: return

        ft1=os.path.join(self.savedir,'sim_events_0000.fits')

        if not os.path.exists(ft1): raise NoSimulatedPhotons()

        if self.use_existing_ft2:

            if self.gtifile is not None:
                self.gtmktime_from_file(evfile=ft1,scfile=self.ft2,outfile=self.ft1, gtifile=self.gtifile)
            else:
                shutil.move(ft1,self.ft1)

        else:
            ft2=os.path.join(self.savedir,'sim_scData_0000.fits')
            shutil.move(ft2,self.ft2)

            shutil.move(ft1,self.ft1)

    def __del__(self):
        """ Remove folder with simulation stuff. """
        if not self.save_output:
            if not self.quiet: print 'Removing directory',self.savedir
            shutil.rmtree(self.savedir)


class SpectralAnalysisMC(SpectralAnalysis):
    """ The intention of this object is to act just like SpectralAnalysis
        but with the difference that it will go ahead and simulate the
        data for you.

        When you create this object, you can pass in a DataSpecification
        object with an existing ft2 file and a not existing ft1 file and
        the ft1 file will be created with the specified pointing history.

        Otherwise, specify tstart and tstop and the data will be simulated
        using a default rocking profile for the desired time range.

        The ROI will be simulated self consistently using the
        point_sources and diffuse_sources passed to the roi() function
        and the resulting ROI with simulated data is returned. """

    defaults = SpectralAnalysis.defaults + (
        'keywords for gtobssim simulation',
    )

    defaults += tuple(get_row(MonteCarlo.defaults,i) for i in ['seed', 'ltfrac', 'tempbase', 'savedir'])

    for i in ['tstart', 'tstop']:
        defaults=change_defaults(defaults, i, get_default(MonteCarlo.defaults,i))

    defaults=change_defaults(defaults,'event_class',0)

    @decorate(defaults)
    def __init__(self, data_specification, **kwargs):
        """ Don't do anything here. """
        self.dataspec = data_specification

        process(self, kwargs)

        if self.event_class != 0:
            raise Exception("event_class must be set to 0 for MC data.")

    def roi(self, roi_dir, point_sources = [], diffuse_sources = [], **kwargs):
        """ First, simulate the requested data. 
            Then properly create a new spectral analysis object.
            And use it to return an ROI. 
            
            This is not a great design structure. 
            But since this object doesn't know what the point + diffuse
            soruces are until the roi() fucntion is called, I couldn't find
            a better way to do it.
            
            Also note that roi_from_xml will correctly be inherited from
            SpectralAnalysis and call this function. """

        assert len(self.dataspec.ft1files) == 1

        if not os.path.exists(self.dataspec.ft1files[0]):
            monte_carlo=MonteCarlo(
                ft1=self.dataspec.ft1files, 
                gtifile = self.dataspec.ltcube if os.path.exists(self.dataspec.ltcube) else None,
                sources = point_sources + diffuse_sources,
                seed=self.seed, 
                tempbase=self.tempbase,
                tstart=self.tstart,
                tstop=self.tstop, 
                ft2=self.dataspec.ft2files,
                ltfrac=self.ltfrac, 
                emin=self.emin,
                emax=self.emax,
                conv_type=self.conv_type,
                irf=self.irf,
                roi_dir=roi_dir,
                maxROI=self.maxROI,
                mc_energy=self.mc_energy,
                quiet=self.quiet,
                savedir=self.savedir)

            monte_carlo.simulate()
            del(monte_carlo)

        # Create a new SpectralAnalysis object with
        # the now existing ft1/ft2 files (Yo Dawg!)
        if self.tstart is None: self.tstart = 0
        if self.tstop is None: self.tstop = 0
        sa=SpectralAnalysis(self.dataspec,**keyword_options.defaults_to_kwargs(self,SpectralAnalysis))
        return sa.roi(roi_dir=roi_dir,
                      point_sources=point_sources, 
                      diffuse_sources=diffuse_sources,**kwargs)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
