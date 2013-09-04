"""
Set up an ROI factory object

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roisetup.py,v 1.31 2013/05/18 23:31:17 burnett Exp $

"""
import os, sys, types
import numpy as np
import pandas as pd
import skymaps
from . import dataset, skymodel, diffuse , roi_bands
from .. utilities import keyword_options, convolution
from .. like import roi_extended, pypsf, pycaldb, pointspec2 #, roi_bands
from .. data import dataman

        
class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime.
    
    Also handles an ad-hoc exposure correction
    """

    def __init__(self, dataset, **datadict): 
        """
        Parameters
        ----------
        dataset :  DataSet object
            for CALDB, aeff, some parameters
            
        datadict['exposure-correction'] : list of strings defining functions of energy
            the correction factors to apply to front, back 
        """

        def make_exposure():
            if dataset.exposure_cube is not None:
                ## use pregenerated gtexpcube2 cube; turn off interpolation
                return [skymaps.DiffuseFunction(f,1000.,False) for f in dataset.exposure_cube]
                
            skymaps.EffectiveArea.set_CALDB(dataset.CALDBManager.CALDB)
            skymaps.Exposure.set_cutoff(np.cos(np.radians(dataset.thetacut)))
            inst = ['front', 'back']
            aeff_files = dataset.CALDBManager.get_aeff()
            ok = [os.path.exists(file) for file in aeff_files]
            if not all(ok):
                raise DataSetError('one of CALDB aeff files not found: %s' %aeff_files)
            self.ea  = [skymaps.EffectiveArea('', file) for file in aeff_files]
            if dataset.verbose: print ' -->effective areas at 1 GeV: ', \
                    ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
            
            if dataset.use_weighted_livetime and hasattr(dataset, 'weighted_lt'):
                return [skymaps.Exposure(dataset.lt,dataset.weighted_lt,ea) for ea in self.ea]
            else:
                return  [skymaps.Exposure(dataset.lt,ea) for ea in self.ea]
                
        self.exposure = make_exposure()

        correction = datadict.pop('exposure_correction', None)
        if correction is not None:
            self.correction = map(eval, correction)
            energies = [100, 1000, 10000]
            print 'Exposure correction: for energies %s ' % energies
            for i,f in enumerate(self.correction):
                print ('\tfront:','\tback: ','\tdfront:', '\tdback')[i], map( f , energies)
        else:
            self.correction = lambda x: 1.0, lambda x: 1.0
        iem_correction = datadict.pop('iem_correction', None)
        if iem_correction is not None:
            self.correction.append( eval(iem_correction))
            self.correction.append( eval(iem_correction))
            energies = [133, 237, 421, 749, 1.333]
            print 'Diffuse correction: for energies %s ' % energies
            for i,f in enumerate(self.correction[3:]):
                print ('\tfront: ','\tback: ',)[i], map( f , energies)
        self.systematic = datadict.pop('systematic', None)
        if self.systematic is not None:
            print 'Galactic diffuse systematic uncertainty: %f' % self.systematic
            
    def value(self, sdir, energy, event_class):
        return self.exposure[event_class].value(sdir, energy)*self.correction[event_class](energy)
        
class ExposureCorrection(object):
    """ logarithmic interpolation function
    """
    def __init__(self, a,b, ea=100, eb=300):
        self.c = (b-a)/np.log(eb/ea)
        self.d =  a -self.c*np.log(ea)
        self.a, self.b = a,b
        self.ea,self.eb = ea,eb
    def __call__(self, e):
        if e>self.eb: return self.b
        if e<self.ea: return self.a
        return self.c*np.log(e) + self.d
    def plot(self, ax=None, **kwargs):
        import pylab as plt
        if ax is None: 
            ax = plt.gca()
        dom = np.logspace(1.5, 3.5, 501) 
        ax.plot(dom, map(self, dom), lw=2, color='r', **kwargs)
        plt.setp(ax, xscale='log', xlim=(dom[0], dom[-1]))
 
class DiffuseCorrection(ExposureCorrection):
    """ binned low-energy diffuse correction functor
    """
    def __init__(self, vals, binsperdec=4):
        """ vals: array-type of float
                values to correct at the energies from 100 MeV defined by binsperdec
        """
        self.vals, self.binsperdec=vals, binsperdec
        
    def __call__(self, e):
        loc = np.log10(e/100.)*self.binsperdec
        if loc<0 or loc>=len(self.vals): return 1.0
        return self.vals[int(loc)]
        

class ROIfactory(object):
    """
    combine the dataset and skymodel for an ROI
    
    """
    defaults =(
        ('analysis_kw', dict(irf=None,minROI=5,maxROI=5, emin=100, emax=316277, quiet=False),'roi analysis keywords'),
        ('skymodel_kw', {}, 'skymodel keywords'),
        ('convolve_kw', dict( resolution=0.125, # applied to OTF convolution: if zero, skip convolution
                            pixelsize=0.05, # ExtendedSourceConvolution
                            num_points=25), # AnalyticConvolution
                                    'convolution parameters'),
        ('irf', None,  'Set to override saved value with the skymodel: expect to find in custom_irf_dir'),
        ('diffuse', None, 'Set to override saved value with skymodel'),
        ('extended', None, 'Set to override saved value with skymodel'),
        ('selector', skymodel.HEALPixSourceSelector,' factory of SourceSelector objects'),
        ('data_interval', 0, 'Data interval (e.g., month) to use'),
        ('nocreate', True, 'Do not allow creation of a binned photon file'),
        ('galactic_correction', None, 'Name of file with diffuse correction factors'), 
        ('galactic_systematic', None, 'Systematic uncertainty for the galactic counts'), 
        ('quiet', False, 'set to suppress most output'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, modeldir, dataspec=None, **kwargs):
        """ 
        parameters
        ----------
        modeldir: folder containing skymodel definition
        dataspec : string or dict or None
                used to look up data specification
                if string, equivalent to dict(dataname=dataspec); otherwise the dict must have
                a dataname element
                if None, use the data used to generate the skymodel
        """
        keyword_options.process(self, kwargs)
        print 'ROIfactory setup: \n\tskymodel: ', modeldir
        # extract parameters used by skymodel for defaults
        input_config = eval(open(os.path.expandvars(modeldir+'/config.txt')).read())
        for key in 'extended diffuse irf galactic_correction galactic_systematic'.split():
            if self.__dict__[key] is None or self.__dict__[key]=='None': 
                #print '%s: %s replace from skymodel: "%s"' %(key, kwargs.get(key,None), input_config.get(key,None))
                self.__dict__[key]=input_config.get(key, None)

        # check for skymodel_kw in the config,txt file, use if found
        config_skymodel = input_config.get('skymodel_config', None)
        if config_skymodel is not None:
            self.skymodel_kw.update(config_skymodel)
            print 'Using skymodel config: %s' % self.skymodel_kw
            
        self.skymodel = skymodel.SkyModel(modeldir, diffuse=self.diffuse,  **self.skymodel_kw)

        if isinstance(dataspec,dataman.DataSet):
            self.dataset = dataspec
            self.dataset.CALDBManager = pycaldb.CALDBManager(irf = self.analysis_kw.get('irf',None),
                                                             psf_irf = self.analysis_kw.get('psf_irf',None),
                                                             CALDB = self.analysis_kw.get('CALDB',None),
                                                             custom_irf_dir=self.analysis_kw.get("irfdir",None))
            self.data_manager = self.dataset(self.data_interval) #need to save a reference to avoid a segfault
            self.exposure = pointspec2.ExposureManager(self.data_manager,self.dataset.CALDBManager)

            self.exposure.correction = [lambda e: 1,lambda e : 1] #TODO
             
        else: # not a DataSet or dict
            if dataspec is None:
                print 'dataspec is None: loading datadict from skymodel.config'; sys.stdout.flush()
                datadict = self.skymodel.config['datadict']
                if type(datadict)==types.StringType: datadict=eval(datadict)
                if isinstance(datadict, dataman.DataSet): 
                    interval = self.skymodel.config.get('interval', None)
                    if interval is None: interval = self.skymodel.config.get('data_interval', None)
                    assert interval is not None, 'did not find interval or data_interval in skymodel.config'
                    dset = datadict[interval]
                    assert hasattr(dset, 'binfile'), 'Not a DataSet? %s' % dataset
                    datadict = dict(dataname=dset)
                else:
                    assert type(datadict)==types.DictType, 'expected a dict'
            elif isinstance(dataspec, dict):
                # a dictionary
                datadict = dataspec
            else:
                print 'looking up dataspec %s' % dataspec
                datadict = dict(dataname=dataspec)\
                        if type(dataspec)!=types.DictType else dataspec
            print '\tdatadict: ', datadict
            if True: #self.analysis_kw.get('irf',None) is None:
                t = self.__dict__['irf']
                if t[0] in ('"',"'"): t = eval(t)
                self.analysis_kw['irf'] = t
            print '\tirf:\t%s' % self.analysis_kw['irf'] ; sys.stdout.flush()
            #datadict = dict(dataname=dataspec, ) \
            #        if type(dataspec)!=types.DictType else dataspec
            exposure_correction=datadict.pop('exposure_correction', None)        
            self.dataset = dataset.DataSet(datadict['dataname'], interval=datadict.get('interval',None),
                    nocreate = self.nocreate,
                    **self.analysis_kw)
            print self.dataset
            self.exposure = ExposureManager(self.dataset, exposure_correction=exposure_correction, systematic=self.galactic_systematic)
        
        if 'CUSTOM_IRF_DIR' not in os.environ and os.path.exists(os.path.expandvars('$FERMI/custom_irfs')):
            os.environ['CUSTOM_IRF_DIR'] = os.path.expandvars('$FERMI/custom_irfs')
        self.psf = pypsf.CALDBPsf(self.dataset.CALDBManager)
 
        # set up glactic correction if requested, and file exists
        if self.galactic_correction is not None:
           if os.path.exists(self.galactic_correction):
               self.dcorr = pd.read_csv(self.galactic_correction, index_col=0)
           else:
               print 'Galactic correction file "%s" not found, no correction applied' % self.galactic_correction
               self.galactic_correction=None
        convolution.AnalyticConvolution.set_points(self.convolve_kw['num_points'])
        convolution.ExtendedSourceConvolution.set_pixelsize(self.convolve_kw['pixelsize'])

    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        show = """analysis_kw selector""".split()
        for key in show:
            s += '\t%-20s: %s\n' %(key,
                self.__dict__[key] if key in self.__dict__.keys() else 'not in self.__dict__!')
        return s

    def _diffuse_sources(self, src_sel):
        """ return  the diffuse, global and extended sources
        """
        skydir = src_sel.skydir()
        assert skydir is not None, 'should use the ROI skydir'
        # get all diffuse models appropriate for this ROI
        globals, extended = self.skymodel.get_diffuse_sources(src_sel)

        try:
            if hasattr(self,'data_manager'):
                bpd = self.data_manager.dataspec.binsperdec
            else:
                bpd = self.dataset.binsperdec
            global_models = [diffuse.mapper(self, src_sel.name(), skydir, source, binsperdec = bpd) for source in globals]
        except Exception, msg:
            print self.dataset, msg
            raise

        def extended_mapper( source):
            if not self.quiet:
                print 'constructing extended model for "%s", spatial model: %s' \
                    %(source, source.spatial_model.__class__.__name__)
            return roi_extended.ROIExtendedModel.factory(self,source,skydir)
        extended_models = map(extended_mapper, extended)
        return global_models, extended_models

    def _local_sources(self, src_sel):
        """ return the local sources with significant overlap with the ROI
        """
        ps = self.skymodel.get_point_sources(src_sel)
        return np.asarray(ps)

    def roi(self, *pars, **kwargs):
        """ return an object based on the selector, with attributes for creating  roi analysis:
            list of ROIBand objects
            list of models: point sources and diffuse
        pars, kwargs : pass to the selector
        """
        roi_kw = kwargs.pop('roi_kw',dict())
        # allow parameter to be a name or a direction
        sel = pars[0]
        source_name=None
        if type(sel)==types.IntType:
            index = int(sel) # needs to be int if int type
        elif type(sel)==skymaps.SkyDir:
            index = self.skymodel.hpindex(sel)
        elif type(sel)==types.StringType:
            index = self.skymodel.hpindex(self.skymodel.find_source(sel).skydir)
            source_name=sel
        elif type(sel)==tuple and len(sel)==2: # interpret a tuple of length 2 as (ra,dec)
            index = self.skymodel.hpindex(skymaps.SkyDir(*sel))
        else:
            raise Exception( 'factory argument "%s" not recognized.' %sel)
        ## preselect the given source after setting up the ROI
        ## (not implemented here)
        #
        src_sel = self.selector(index, **kwargs)

        class ROIdef(object):
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)
            def __str__(self):
                return 'ROIdef for %s' %self.name
        skydir = src_sel.skydir()  
        global_sources, extended_sources = self._diffuse_sources(src_sel)
        if isinstance(self.dataset,dataman.DataSet):
            bandsel = BandSelector(self.data_manager,
                                   self.psf,
                                   self.exposure,
                                   skydir)
            bands = bandsel(minROI = self.analysis_kw['minROI'],
                            maxROI = self.analysis_kw['maxROI'],
                            emin = self.analysis_kw['emin'],
                            emax = self.analysis_kw['emax'])
        else:
            bands = self.dataset(self.psf,self.exposure,skydir)
            
        # setup galactic diffuse correction and/or counts to pass into ROI generation
        if self.galactic_correction is not None:
            # look up the diffuse correction factors for this ROI, by name
            self.exposure.dcorr = self.dcorr.ix[src_sel.name()].values

        return ROIdef( name=src_sel.name() ,
                    roi_dir=skydir, 
                    bands=bands,
                    global_sources=global_sources,
                    extended_sources = extended_sources,
                    point_sources=self._local_sources(src_sel), 
                    exposure=self.exposure,
                    **roi_kw)
                    
    def __call__(self, *pars, **kwargs):
        """ alias for roi() """
        return self.roi(*pars, **kwargs)
        
    def reload_model(self):
        """ Reload the sources in the model """
        self.skymodel._load_sources()

class BandSelector(object):
    """Class to handle selection of bands with new data management code.
    
    Should be moved somewhere more sensible."""

    def __init__(self,data_manager,psf,exposure,roi_dir):
        class SA(object):
            def __init__(self,**kw):self.__dict__.update(kw)
        #make sure exposure is the right kind of ExposureManager and that
        #data_manager is the same one use for the exposure
        assert(hasattr(exposure,'data_manager'))
        data_manager.dataspec.check_consistency(exposure.data_manager.dataspec)
        self.sa = SA(psf=psf,exposure=exposure)
        self.data_manager = data_manager
        self.roi_dir = roi_dir
    
    def __call__(self,minROI=7,maxROI=7,emin=100,emax=1e6,band_kw=dict()):
        self.sa.minROI = minROI
        self.sa.maxROI = maxROI
        bands = []
        for i,band in enumerate(self.data_manager.bpd):
            if (band.emin() + 1) >= emin and (band.emax() - 1) < emax:
                bands.append(roi_bands.ROIBand(band, self.sa, self.roi_dir,**band_kw))
        return np.asarray(bands)

        

def main(modeldir='3years/uw10', skymodel_kw={}):
    rf = ROIfactory(modeldir, **skymodel_kw)
    return rf
