"""
Module implements a wrapper around gtobssim to allow
less painful simulation of data.

$Header:

author: Joshua Lande
"""
import os
import re
import shutil
import collections
from tempfile import mkdtemp,NamedTemporaryFile
from GtApp import GtApp

import numpy as N
from skymaps import IsotropicSpectrum,IsotropicPowerLaw,DiffuseFunction,PySkyFunction
from . SpatialModels import Gaussian,EllipticalGaussian,RadiallySymmetricModel
from . Models import PowerLaw,PowerLawFlux,Constant
from . pointspec import DataSpecification,SpectralAnalysis
from . pointspec_helpers import PointSource
from . roi_extended import ExtendedSource
from . SpatialModels import Gaussian,EllipticalGaussian,RadiallySymmetricModel

from uw.utilities import keyword_options 

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
            ('point_sources',     [], "List of Point Sources to use in the simulation.."),
            ('diffuse_sources',   [], "List of Diffuse Sources to use in the simulation."),
            ('seed',               0, "Random number generator seen when running obssim. This should be varied for simulations."),
            ('tempbase', '/scratch/', "Directory to put temporary files in. Can be cleaned up with the cleanup function."),
            ('tstart',          None, "Required if ft2 does not point to a real file."),
            ('tstop',           None, "Required if ft2 does not point to a real file."),
            ('ft2',             None, "If exists, simulate with this ft2 file. If not, create."),
            ('ltfrac',          None, "If set, pass value into gtobssim."),
            ('emin',             100, "Minimum energy"),
            ('emax',          100000, "Maximum energy"),
            ('conv_type',         -1, "Conversion Type"),
            ('irf',             None, "Instrument Response Function."),
            ('roi_dir',         None, "Center of ROI. Gtobssim will use the use_ac flag if this is specified."),
            ('maxROI',          None, "Maximum ROI Size. Gtobssim will use the use_ac flag if this is specified."),
            ('quiet',          False, "Surpress output."),
            ('mc_energy',      False, "Use MC Energy"),
    )

    @staticmethod
    def strip(name):
        return re.sub('[ \.()]','',name)

    @keyword_options.decorate(defaults)
    def __init__(self,ft1,**kwargs):
        """ Constructor does not require a data_specification. """
        keyword_options.process(self, kwargs)

        self.ft1=ft1

        if isinstance(self.ft1,collections.Iterable):
            if len(self.ft1) != 1: raise Exception("...")
            self.ft1=self.ft1[0]

        if isinstance(self.ft2,collections.Iterable):
            if len(self.ft2) != 1: raise Exception("...")
            self.ft2=self.ft2[0]

        if os.path.exists(self.ft1):
            raise Exception("Unable to run MC simulation because file %s already exists." % self.ft1)

        if not os.path.isdir(self.tempbase):
            raise Exception("Argument tempbase must point to a directory where temporary files can be placed.")

        if self.ft2 is None:
            if self.tstart is None or self.tstop is None:
                raise Exception("tstart and tstop must be specified if ft2 is not set.")
        else:
            if os.path.exists(self.ft2):
                if self.tstart is not None or self.tstop is not None:
                    raise Exception("Since ft2 points to an existing file, tstart and tstop can not be set.")
            else:
                if self.tstart is None or self.tstop is None:
                    raise Exception("Since ft2 does not point to an already existing file, tstart and tstop must be set.")

        if self.tstop is not None and self.tstart is not None and (self.tstop - self.tstart < 1):
            raise Exception("tstart and tstop must describe a real range.")

        if self.ft2 is not None and os.path.exists(self.ft2):
            self.set_time_from_ft2()

        if self.irf is None:
            raise Exception("...")

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
        return 0.5*self.emin,1.5*self.emax

    def _make_ps(self,ps,indent='  '):

        assert isinstance(ps,PointSource)

        # account for energy dispersion by conservative 50% oversimulation of events.
        mc_emin,mc_emax=self.larger_energy_range()

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
            flux=ps.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4
            specfile=self._make_specfile(ps.model,mc_emin,mc_emax)
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

    def _make_specfile(self,model,emin,emax,numpoints=200):
        temp=NamedTemporaryFile(dir='.',delete=False)
        energies=N.logspace(N.log10(emin),N.log10(emax),numpoints)
        fluxes=model(energies)
        temp.write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(energies,fluxes)]))
        temp.close()
        return temp.name

    def _make_profile(self,spatial_model,numpoints=200):
        temp=NamedTemporaryFile(dir='.',delete=False)
        edge=spatial_model.effective_edge()
        radius=N.linspace(0,edge,numpoints)
        pdf=spatial_model.at_r_in_deg(radius)
        temp.write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(radius,pdf)]))
        temp.close()
        return temp.name

    def _make_radially_symmetric(self,es,indent):

        sm=es.spatial_model

        assert isinstance(sm,RadiallySymmetricModel)
 
        mc_emin,mc_emax=self.larger_energy_range()

        flux=es.model.i_flux(mc_emin,mc_emax,cgs=True)*1e4

        ra,dec=sm.center.ra(),sm.center.dec()

        profile=self._make_profile(sm)
        specfile=self._make_specfile(es.model,mc_emin,mc_emax)

        xml=[
            '<source name="%s">' % MonteCarlo.strip(es.name),
            '  <spectrum escale="MeV">',
            '    <SpectrumClass name="RadialSource"',
            '      params="flux=%s, profileFile=%s, specFile=%s, ra=%s, dec=%s"/>' % \
                    (flux,profile,specfile,ra,dec),
            '    <use_spectrum frame="galaxy"/>',
            '  </spectrum>',
            '</source>',
        ]
        return indent+('\n'+indent).join(xml)

    def _make_gaussian(self,es,indent='  '):
        """ Make an extended source. """

        assert isinstance(es.spatial_model,Gaussian) or \
               isinstance(es.spatial_model,EllipticalGaussian)

        mc_emin,mc_emax=self.larger_energy_range()

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
    def make_isotropic_allsky(filename):
        img=SkyImage(SkyDir(0,0,SkyDir.GALACTIC),filename,6,180,1,"CAR",True)
        
        one=lambda x: 1
        skyfun=PySkyFunction(one)
        img.fill(skyfun)
        del(img)

    @staticmethod
    def isotropic_integrator(filename):
        """ Assume a powerlaw connecting each point. 

            This code is adopated from the function pl_integral in 
            genericSources/src/FileSpectrum.cxx which is used by gtobssim
            to integrate the isotropic spectrum.
        
            Returns a value suitable for gtobssim in units of ph/m^2 """
        file=N.genfromtxt(filename,unpack=True)
        energy,flux=file[0],file[1]

        e1,e2=energy[:-1],energy[1:]
        f1,f2=flux[:-1],flux[1:]

        gamma=N.log(f2/f1)/N.log(e2/e1)

        n0 = f1/e1**gamma

        return N.sum(
            N.where(gamma != -1,
                    n0/(gamma+1)*(e2**(gamma+1)-e1**(gamma+1)),
                    n0*N.log(e2/e1)
            )
        )*4*N.pi*10**4

    def _make_isotropic(self,ds,savedir,indent):

        dm=ds.dmodel[0]
        sm=ds.smodel

        isotropic_spectrum=dm.name()

        if not (isinstance(sm,Constant) and sm.getp(0) == 1):
            raise Exception("Can only run gtobssim with IsotropicSpectrum diffuse models where the spectral model is a Constant with norm 1.")

        energies=N.genfromtxt(isotropic_spectrum,unpack=True)[0]

        emin,emax=energies[0],energies[-1]

        # multiply by 4pi * 10^4 to convert from ph/cm^2/sr to ph/m^2
        flux=MonteCarlo.isotropic_integrator(isotropic_spectrum)

        isotropic_allsky=os.path.join(savedir,'isotropic_allsky.fits')
        MonteCarlo.make_isotropic_allsky(isotropic_allsky)

        ds = [ 
            '<source name="%s">' % MonteCarlo.strip(ds.name),
            '  <spectrum escale="MeV">',
            '    <SpectrumClass name="FileSpectrumMap"',
            '       params="flux=%s,fitsFile=%s,' % (flux,isotropic_allsky),
            '       specFile=%s,emin=%s,emax=%s"/>' % (isotropic_spectrum,emin,emax),
            '    <use_spectrum frame="galaxy"/>',
            '  </spectrum>',
            '</source>'
        ]
        return indent+('\n'+indent).join(ds)

    def _make_powerlaw(self,ds,savedir,indent):
        dm=ds.dmodel[0]
        sm=ds.smodel

        if not (isinstance(sm,PowerLaw) and N.all(sm._p==0) and sm.index_offset==1): 
            raise Exception("Can only run gtobssim with DiffuseFunction diffuse models where the spectral model is a PowerLaw with norm and index 1.")

        # flux in ph/cm^2/s/sr b/n 100MeV & infinity
        flux_pointlike=dm.flux()
        index=dm.index()

        x=PowerLawFlux(p=[flux_pointlike,index],emin=100,emax=N.inf)

        mc_emin,mc_emax=self.larger_energy_range()


        if self.roi_dir is not None and self.maxROI is not None:
            ra,dec=self.roi_dir.ra(),self.roi_dir.dec()
            # just to be safe from the large low energy psf
            radius=self.maxROI+20 

            # gtobssim wants ph/m^2/s, integrate over solid angle
            flux=x.i_flux(mc_emin,mc_emax)*(2*N.pi*(1-N.cos(N.radians(radius))))*10**4
        else:
            ra,dec,radius=0,0,180

            # gtobssim wants ph/m^2/s
            flux=x.i_flux(mc_emin,mc_emax)*4*N.pi*10**4


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

        # get the solid angles from gtlike. Probably
        # could be rewritten in pointlike, but not necessary
        import pyLikelihood
        mapcube=pyLikelihood.MapCubeFunction2(filename)
        wcsmap=mapcube.wcsmap()
        solid_angles=N.empty((nxpix,nypix))
        sa=wcsmap.solidAngles()
        for i in xrange(nypix): solid_angles[:,i]=sa[i]

        map_integral = \
                N.sum(solid_angles*\
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
        #return N.transpose((y1+y2)/2 * (x2-x1))

        gamma = N.log(y2/y1)/N.log(x2/x1)
        n0 = y1/x1**gamma

        gp1 = gamma + 1.;
        integral = N.where(gamma != 1., n0/gp1*(x2**gp1 - x1**gp1), n0*N.log(x2/x1))
        return N.transpose(integral)

    def _make_diffuse(self,ds,indent):

        dm=ds.dmodel[0]
        sm=ds.smodel

        # galactic diffuse
        if not (isinstance(sm,PowerLaw) and N.all(sm._p==0) and sm.index_offset==1): 
            raise Exception("Can only run gtobssim with DiffuseFunction diffuse models where the spectral model is a PowerLaw with norm and index 1.")

        filename=dm.name()

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

    def _make_ds(self,ds,savedir,indent='  '):
        if len(ds.dmodel) > 1:
            raise Exception("Can only run gtobssim for diffuse models with a combined front/back diffuse model.")

        if isinstance(ds.dmodel[0],IsotropicSpectrum):
            return self._make_isotropic(ds,savedir,indent)

        elif isinstance(ds.dmodel[0],IsotropicPowerLaw):
            return self._make_powerlaw(ds,savedir,indent)

        elif isinstance(ds.dmodel[0],DiffuseFunction):
            return self._make_diffuse(ds,indent)

        elif isinstance(ds,ExtendedSource):
            if isinstance(ds.spatial_model,Gaussian) or isinstance(ds.spatial_model,EllipticalGaussian):
                return self._make_gaussian(ds,indent)
            elif isinstance(ds.spatial_model,RadiallySymmetricModel):
                return self._make_radially_symmetric(ds,indent)
            else:
                raise Exception("Can simulate diffuse source %s. No code to simulate general non-radially symmetric sources." % ds.name)
        else:
            raise Exception("Can simulae diffuse source %s to xml suitable for gtobssim." % ds.name)
            
    def _make_model(self,savedir,indent='  '):

        xml = ['<source_library title="Library">']
        src = []

        if not isinstance(self.point_sources,collections.Iterable):
            self.point_sources=[self.point_sources]

        if not isinstance(self.diffuse_sources,collections.Iterable):
            self.diffuse_sources=[self.diffuse_sources]

        for ps in self.point_sources:
            if not self.quiet: print 'Processing source %s' % ps.name
            xml.append(self._make_ps(ps,indent=indent))
            src.append(MonteCarlo.strip(ps.name))

        for ds in self.diffuse_sources:
            if not self.quiet: print 'Processing source %s' % ds.name
            xml.append(self._make_ds(ds,savedir,indent=indent))
            src.append(MonteCarlo.strip(ds.name))

        xml.append('</source_library>')
    
        self.xmlfile=os.path.join(savedir,'source_library.xml')
        temp=open(self.xmlfile,'w')
        temp.write('\n'.join(xml))
        temp.close()

        self.srclist=os.path.join(savedir,'source_list.txt')
        temp=open(self.srclist,'w')
        temp.write('\n'.join(src))
        temp.close()

    def simulate(self):

        old_dir=os.getcwd()
        self.tempdir=mkdtemp(prefix=self.tempbase)
        if not self.quiet: print 'working in tempdir',self.tempdir
        os.chdir(self.tempdir)
        self._make_model(savedir=self.tempdir)

        irfs=self.irf

        if self.conv_type == 0:
            irfs+="::FRONT"
        elif self.conv_type == 1:
            irfs+="::BACK"

        if self.roi_dir is not None and self.maxROI is not None:
            use_ac="yes"
            ra=self.roi_dir.ra()
            dec=self.roi_dir.dec()
            radius=self.maxROI
        else:
            use_ac="no"
            ra,dec,radius=0,0,180

        app=GtApp('gtobssim')
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
                maxrows=200000000 # make sure all photons are in one file
        )

        os.chdir(old_dir)

        ft1=os.path.join(self.tempdir,'sim_events_0000.fits')
        if os.path.exists(ft1):
            shutil.copy(ft1,self.ft1)
        else:
            raise NoSimulatedPhotons()


        ft2=os.path.join(self.tempdir,'sim_scData_0000.fits')
        if self.ft2 is not None and not os.path.exists(self.ft2): 
            shutil.copy(ft2,self.ft2)

    def __del__(self):
        """ Remove folder with simulation stuff. """
        shutil.rmtree(self.tempdir)


class SpectralAnalysisMC(SpectralAnalysis):
    """ The intention of this object is to act just like SpectralAnalysis
        but with the difference that it will go ahead and simulate the
        data for you.

        When you create this object, you can pass in a DataSpecifciation
        object with an existing ft2 file and a not existing ft1 file and
        the ft1 file will be created with the specified pointing history.

        Otherwise, specify tstart and tstop and the data will be simulated
        using a default rocking profile for the desired time range.

        The ROI will be simulated self consistently using the
        point_sources and diffuse_sources passed to the roi() function
        and the resulting ROI with simulated data is returned. """

    defaults = SpectralAnalysis.defaults + (
            ('seed',               0, "See MonteCarlo."),
            ('ltfrac',          None, "See MonteCarlo."),
            ('tempbase', '/scratch/', "See MonteCarlo.")
    )
    keyword_options.change_defaults(defaults,'event_class',0)
    keyword_options.change_defaults(defaults,'tstart',None)
    keyword_options.change_defaults(defaults,'tstop',None)

    @keyword_options.decorate(defaults)
    def __init__(self, data_specification, **kwargs):
        """ Don't do anything here. """
        self.dataspec = data_specification

        keyword_options.process(self, kwargs)

        if self.event_class != 0:
            raise Exception("event_class must be set to 0 for MC data.")

    def roi(self, roi_dir,
            point_sources = [], catalogs = [], catalog_mapper = None,
            diffuse_sources = [], catalog_include_radius = None,
            **kwargs):
        """ First, simulate the requested data. 
            Then properly create a new spectral analysis object.
            And use it to return an ROI. 
            
            This is not a great design structure. 
            But since this object doesn't know what the point + diffuse
            soruces are until the roi() fucntion is called, I couldn't find
            a better way to do it.
            
            Also note that roi_from_xml will correctly be inherited from
            SpectralAnalysis and call this function. """

        if not isinstance(catalogs,collections.Iterable): catalogs = [catalogs]
        for cat in catalogs:
            if not isinstance(cat,PointSourceCatalog):
                cat = catalog_mapper(cat)
            point_sources,diffuse_sources = cat.merge_lists(roi_dir,
                    self.maxROI+5 if catalog_include_radius is None else catalog_include_radius,
                    point_sources,diffuse_sources)

        monte_carlo=MonteCarlo(
                ft1=self.dataspec.ft1files, 
                point_sources=point_sources, 
                diffuse_sources=diffuse_sources,
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
                mc_energy=self.mc_energy)

        monte_carlo.simulate()

        ds=DataSpecification(
                ft1files=monte_carlo.ft1,
                ft2files=monte_carlo.ft2,
                binfile=self.dataspec.binfile,
                ltcube=self.dataspec.ltcube
        )
    
        # Create a new SpectralAnalysis object with
        # the now existing ft1/ft2 files (Yo Dawg!)
        sa=SpectralAnalysis(ds,**keyword_options.defaults_to_kwargs(self,SpectralAnalysis))
        return sa.roi(roi_dir=roi_dir,
                      point_sources=point_sources, 
                      diffuse_sources=diffuse_sources,**kwargs)
