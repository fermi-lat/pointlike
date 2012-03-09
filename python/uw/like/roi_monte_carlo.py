"""
Module implements a wrapper around gtobssim to allow
less painful simulation of data.

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/roi_monte_carlo.py,v 1.41 2012/03/07 03:47:20 lande Exp $

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

import pyLikelihood

# I think this is the smallest possible float for 32 bit numbers
# http://www.psc.edu/general/software/packages/ieee/ieee.php
SMALL_NUMBER = 2**-149

class DiffuseShrinker(object):
    defaults = (
        ('galactic', False, "Save resulting fits file in galactic coordiantes"),
        ('proj', 'CAR', "Projection"),
    )

    @staticmethod
    def smaller_range(values, min, max):
        """ values should be a sorted list of numbers.

            This function will return a boolean array of the values
            between min and max, but the list will be inclusive so that it will
            include the first number smaller than min and the first number
            larger than max.

                >>> print DiffuseShrinker.smaller_range(np.asarray([1,2,3,5]), 2, 3)
                [False  True  True False]

                >>> print DiffuseShrinker.smaller_range(np.asarray([1,2,3,5]), 1.9, 3.1)
                [ True  True  True  True]

            If there is no number smaller than min (or larger than max), it just includes
            the whole list

                >>> print DiffuseShrinker.smaller_range(np.asarray([1,2,3,5]), 0.9, 5.1)
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

    @decorate(defaults)
    def __init__(self, diffuse_model, skydir, radius, emin, emax, **kwargs):
        """ Take in an allsky diffuse model and create a diffuse model which predicts
            emission only inside the given radius and energy range. """

        process(self, kwargs)

        assert type(diffuse_model) == str and os.path.exists(diffuse_model)

        # I am not sure if this algorithm makes sense for very-larger areas of the sky
        # For now, just crash
        assert radius < 90

        self.diffuse_model = diffuse_model
        self.skydir = skydir
        self.radius = radius
        self.emin = emin
        self.emax = emax

        self.allsky_fits = pyfits.open(self.diffuse_model)

        h = self.allsky_fits['PRIMARY'].header
        # Here, make sure diffuse function is allsky
        assert abs(abs(h['NAXIS1']*h['CDELT1'])-360) < 1
        assert abs(abs(h['NAXIS2']*h['CDELT2']) -180) < 1

        # make sure diffuse file has same
        assert h['CDELT1'] == h['CDELT2']
        self.pixelsize = h['CDELT1']

        self.allsky_image = SkyImage(self.diffuse_model)

        self.all_energies = np.asarray(self.allsky_fits['energies'].data.field('energy'))

        # create a list of the indices in the self.all_energy array for energies
        # that are allowed by the energy cut.
        self.good_indices = DiffuseShrinker.smaller_range(self.all_energies, self.emin, self.emax)
        self.good_energies = self.all_energies[self.good_indices]

    def shrink(self, smaller_file):

        shutil.copyfile(self.diffuse_model, smaller_file)

        layers = sum(self.good_indices)

        # first, strip out unwanted energies
        x = pyfits.open(smaller_file, mode='update')
        x['PRIMARY'].data = x['PRIMARY'].data[self.good_indices]
        x['ENERGIES'].data = x['ENERGIES'].data[self.good_indices]

        x['ENERGIES'].header['NAXIS2'] = layers

        x['PRIMARY'].header['NAXIS3'] = layers
        x['PRIMARY'].header['CRVAL3'] = x['ENERGIES'].data.field('ENERGY')[0]

        x.flush()
        del(x)

        # second, set to 0 pixels outside of desired region
        img=SkyImage(smaller_file)
        for i in range(layers):
            wsdl = img.get_wsdl(i)

            weights = np.zeros(len(wsdl),dtype=float)
            arclength = np.zeros(len(wsdl),dtype=float)

            PythonUtilities.get_wsdl_weights(weights, wsdl)
            PythonUtilities.arclength(arclength,wsdl, self.skydir)

            rad_cut = np.degrees(arclength) > self.radius
            weights[rad_cut] = SMALL_NUMBER

            PythonUtilities.set_wsdl_weights(weights, wsdl)

            img.setLayer(i)
            img.set_wsdl(wsdl)
            
        img.save()
        del(img)



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
            ('diffuse_pad',      10,  " How many area outside ROI should the diffuse emission be simulated to account for energy dispersion."),
            ('roi_pad',           5,  " Include counts from this far ouside maxROI in the ft1 file to account for ragged edge effect in pointlike."),
            ('energy_pad',      0.5,  """ Lower energy of simluated photons is energy_pad*emin and 
                                          upper energy of simulated photos is (1+energy_pad)*emax.
                                          This allows for energy dispersion effects to be 
                                          naturally accounted for in the monte carlos simulation. """),
    )

    @staticmethod
    def strip(name):
        return re.sub('[ \.()]','',name)

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
        ft2 = pyfits.open(self.ft2)[1]
        self.tstart = ft2.data.field('START')[0]
        self.tstop = ft2.data.field('STOP')[-1]

    def larger_energy_range(self):
        """ Get an energy range larger then the desired simulated points
            (to correct for energy dispersion). """
        return self.energy_pad*self.emin,(1+self.energy_pad)*self.emax

    def _make_ps(self,ps,mc_emin,mc_emax,indent):

        assert isinstance(ps,PointSource)

        if isinstance(ps.model,PowerLaw) or isinstance(ps.model,PowerLawFlux):
            xml=[
                '<source name="%s" flux="%s">' % (MonteCarlo.strip(ps.name),ps.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4),
                '  <spectrum escale="MeV">',
                '  <particle name="gamma">',
                '    <power_law emin="%s" emax="%s" gamma="%s"/>' % (mc_emin,mc_emax,ps.model.getp(1)),
                '  </particle>',
                '  <celestial_dir ra="%s" dec="%s"/>' % (ps.skydir.ra(),ps.skydir.dec()),
                '  </spectrum>',
                '</source>'
            ]
            return indent+('\n'+indent).join(xml)
        else:
            if isinstance(ps.model,FileFunction):
                flux=ps.model.i_flux(self.energy[0],self.energy[-1],cgs=True)*1e4
                specfile=ps.model.file
            else:
                flux=ps.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4
                specfile=MonteCarlo._make_specfile(ps.name,ps.model,mc_emin,mc_emax)

            xml=[
                '<source name="%s">' % MonteCarlo.strip(ps.name),
                '  <spectrum escale="MeV">',
                '    <SpectrumClass name="FileSpectrum"',
                '      params="flux=%g,specFile=%s"/>' % (flux,specfile),
                '    <celestial_dir ra="%s" dec="%s"/>' % (ps.skydir.ra(),ps.skydir.dec()),
                '  </spectrum>',
                '</source>'
            ]
            return indent+('\n'+indent).join(xml)

    @staticmethod
    def _make_specfile(name,model,emin,emax,numpoints=200):
        temp='%s_spectra_%s.txt' % (MonteCarlo.strip(name),model.name)
        energies=np.logspace(np.log10(emin),np.log10(emax),numpoints)
        fluxes=model(energies)
        # josh's fix to the clipping.../ thanks!
        while fluxes[0] == 0: energies,fluxes=energies[1:],fluxes[1:]
        while fluxes[-1] == 0: energies,fluxes=energies[:-1],fluxes[:-1]
        # okay that's to check for local divergences
        if np.any(fluxes==0): raise Exception("Error: 0s found in differential flux")
        open(temp,'w').write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(energies,fluxes)]))
        return temp

    @staticmethod
    def _make_profile(name,spatial_model,numpoints=200):
        temp='%s_extension_profile_%s.txt' % (MonteCarlo.strip(name),spatial_model.name)
        radius,pdf = spatial_model.approximate_profile()
        open(temp,'w').write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(radius,pdf)]))
        return temp

    def _make_radially_symmetric(self,es,mc_emin, mc_emax, indent):

        name=es.name
        sm=es.spatial_model

        assert isinstance(sm,RadiallySymmetricModel)
 
        flux=es.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4

        ra,dec=sm.center.ra(),sm.center.dec()

        profile_name='%s_extension_profile_%s.txt' % (MonteCarlo.strip(name),sm.name)
        sm.save_profile(profile_name)

        specfile=MonteCarlo._make_specfile(name,es.model,mc_emin,mc_emax)

        xml=[
            '<source name="%s">' % MonteCarlo.strip(es.name),
            '  <spectrum escale="MeV">',
            '    <SpectrumClass name="RadialSource"',
            '      params="flux=%s, profileFile=%s, specFile=%s, ra=%s, dec=%s"/>' % \
                    (flux,profile_name,specfile,ra,dec),
            '    <use_spectrum frame="galaxy"/>',
            '  </spectrum>',
            '</source>',
        ]
        return indent+('\n'+indent).join(xml)

    def _make_extended_source(self,es,mc_emin,mc_emax,indent):

        sm=es.spatial_model

        assert isinstance(sm,SpatialModel)

        flux=es.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4

        ra,dec=sm.center.ra(),sm.center.dec()

        if isinstance(sm,SpatialMap):
            print 'WARNING: gtobssim can only use plate-carree projection fits files!'
            spatialfile=sm.file
        else:
            spatialfile='%s_spatial_template_%s.fits' % (MonteCarlo.strip(es.name),sm.name)
            # Allegedly simulated templates must only be in the plate-carree projection
            # http://www.slac.stanford.edu/exp/glast/wb/prod/pages/sciTools_observationSimTutorial/obsSimTutorial.htm
            sm.save_template(spatialfile, proj='CAR')

        if isinstance(es.model,PowerLaw) or isinstance(es.model,PowerLawFlux):

            index = es.model['index']

            xml=[
                '<source name="%s">' % MonteCarlo.strip(es.name),
                '  <spectrum escale="MeV">',
                '    <SpectrumClass name="MapSource"',
                '      params="%s,%s,%s,%s,%s"/>' % \
                        (flux,index,spatialfile,mc_emin,mc_emax),
                '    <use_spectrum frame="galaxy"/>',
                '  </spectrum>',
                '</source>',
            ]

        else:

            specfile=MonteCarlo._make_specfile(es.name,es.model,mc_emin,mc_emax)

            xml=[
                '<source name="%s">' % MonteCarlo.strip(es.name),
                '  <spectrum escale="MeV">',
                '    <SpectrumClass name="FileSpectrumMap"',
                '      params="flux=%s, fitsFile=%s, specFile=%s, emin=%s, emax=%s"/>' % \
                        (flux,spatialfile,specfile,mc_emin,mc_emax),
                '    <use_spectrum frame="galaxy"/>',
                '  </spectrum>',
                '</source>',
            ]
        return indent+('\n'+indent).join(xml)

    def _make_powerlaw_gaussian(self,es,mc_emin,mc_emax,indent):
        """ Make an extended source. """

        assert isinstance(es.spatial_model,Gaussian) or \
               isinstance(es.spatial_model,EllipticalGaussian)

        if isinstance(es.model,PowerLaw) or isinstance(es.model,PowerLawFlux):

            flux=es.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4
            gamma=es.model.getp('Index')
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
    def make_isotropic_fits(filename, radius, proj="CAR",  pixelsize=1, skydir=None):
        """ Note, if there is an ROI cut, we can make this
            isotropic file not allsky. """
        galactic=True
        gal_center=SkyDir(0,0,SkyDir.GALACTIC)
        img=SkyImage(gal_center,filename,pixelsize,180,1,proj,galactic)

        if radius >= 180:
            fill=lambda v: 1
        else:

            assert skydir is not None

            def fill(v):
                sd = SkyDir(Hep3Vector(*v))
                if np.degrees(sd.difference(skydir)) <= radius:
                    return 1
                else:
                    return SMALL_NUMBER

        skyfun=PySkyFunction(fill)
        img.fill(skyfun)
        del(img)

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

                >>> tempfile = NamedTemporaryFile()
                >>> MonteCarlo.make_isotropic_fits(tempfile.name, radius=180, pixelsize=1)
                >>> map_area = MonteCarlo.spatial_integrator_2d(tempfile.name)
                >>> np.allclose(map_area, 4*np.pi)
                True

            When the fits file contains the value '1' in it only within a given radius, the map 
            integral should should be 4*pi*(1-cos(radius)). Note, we need smaller
            pixel size to get a good agreement due to ragged edge effect at corner of circle.

                >>> for radius in [10, 30, 50]:
                ...     tempfile = NamedTemporaryFile()
                ...     MonteCarlo.make_isotropic_fits(tempfile.name, skydir=SkyDir(134,83,SkyDir.GALACTIC), radius=radius, pixelsize=0.25)
                ...     map_area = MonteCarlo.spatial_integrator_2d(tempfile.name)
                ...     true_area = 2*np.pi*(1-np.cos(np.radians(radius)))
                ...     np.allclose(map_area, true_area, atol=1e-2, rtol=1e-2)
                True
                True
                True
        """
        fits=pyfits.open(filename)
        data=fits[0].data
        assert len(data.shape) == 2 or data.shape[0] == 1
        if data.shape == 3: data = data[1]

        nxpix,nypix=data.shape[1],data.shape[2]

        # get the solid angles from gtlike. Probably
        # could be rewritten in pointlike, but not necessary
        mapcube=pyLikelihood.MapCubeFunction2(filename)
        sa=np.asarray(mapcube.wcsmap().solidAngles()).transpose()
        map_integral = np.sum(sa*data)

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

        smaller_range = DiffuseShrinker.smaller_range(energies, mc_emin, mc_emax)
        if np.any(smaller_range == False):
            # If any energies are outside the range, cut down spectral file.
            cut_spectral_file = os.path.basename(spectral_file).replace('.txt','_cut.txt')
            np.savetxt(cut_spectral_file, zip(energies[smaller_range], spectra[smaller_range]))
        else:
            cut_spectral_file = spectral_file

        spatial_file=os.path.basename(spectral_file).replace('.txt','_spatial.fits')
        if self.roi_dir is not None and self.maxROI is not None:
            radius=self.maxROI + self.diffuse_pad
            MonteCarlo.make_isotropic_fits(spatial_file, skydir=self.roi_dir, radius=radius)
        else:
            # multiply by solid angle
            MonteCarlo.make_isotropic_fits(spatial_file, radius=180)


        # flux is ph/cm^2/sr to ph/m^2
        flux, emin, emax = MonteCarlo.isotropic_spectral_integrator(cut_spectral_file)

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

        map_integral = \
                np.sum(sa*\
                    MonteCarlo.powerLawIntegral(energies[0:-1],energies[1:],
                    data[0:-1,:,:],data[1:,:,:]).sum(axis=0))

        return map_integral

    @staticmethod
    def powerLawIntegral(x1, x2, y1, y2):
        """ This function is blatently stolen from
            MapCubeFunction::powerLawIntegral
            in MapCubeFunction.cxx

            http://www-glast.stanford.edu/cgi-bin/viewcvs/Likelihood/src/MapCubeFunction.cxx
        """
        # transpose to make energy the third dimension so that vector oprations will
        # work seamlessly with the 1 dimesnional energy arrays
        y1=y1.transpose()
        y2=y2.transpose()
        #return np.transpose((y1+y2)/2 * (x2-x1))

        gamma = np.log(y2/y1)/np.log(x2/x1)
        n0 = y1/x1**gamma

        gp1 = gamma + 1.;
        integral = np.where(gamma != 1., n0/gp1*(x2**gp1 - x1**gp1), n0*np.log(x2/x1))
        return np.transpose(integral)

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
            radius=self.maxROI + self.diffuse_pad

            filename = os.path.basename(allsky_filename).replace('.fits','_cut.fits')

            shrinker = DiffuseShrinker(allsky_filename, self.roi_dir, radius, mc_emin, mc_emax)
            shrinker.shrink(filename)

        else:
            filename = allsky_filename


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
            if not self.quiet: print 'Processing source %s' % ps.name
            xml.append(self._make_ps(ps,mc_emin, mc_emax,indent))
            src.append(MonteCarlo.strip(ps.name))

        for ds in self.diffuse_sources:
            if not self.quiet: print 'Processing source %s' % ds.name
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

    def simulate(self,**kwargs):
        ''' understands all keywords that GtApp.run() can handle, especially dry_run=True and verbosity=X '''

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
                **kwargs # to pass all keywords from simulate through to app.run();
       )

        os.chdir(old_dir)

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

    defaults += tuple(get_row(MonteCarlo.defaults,i) for i in ['seed', 'ltfrac', 'tempbase'])

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
            quiet=self.quiet)

        monte_carlo.simulate()

        ds=DataSpecification(ft1files=monte_carlo.ft1,
                             ft2files=monte_carlo.ft2,
                             binfile=self.dataspec.binfile,
                             ltcube=self.dataspec.ltcube)

        # Create a new SpectralAnalysis object with
        # the now existing ft1/ft2 files (Yo Dawg!)
        if self.tstart is None: self.tstart = 0
        if self.tstop is None: self.tstop = 0
        sa=SpectralAnalysis(ds,**keyword_options.defaults_to_kwargs(self,SpectralAnalysis))
        return sa.roi(roi_dir=roi_dir,
                      point_sources=point_sources, 
                      diffuse_sources=diffuse_sources,**kwargs)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
