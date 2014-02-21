"""
Manage the analysis configuration

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/configuration.py,v 1.18 2014/02/11 04:17:42 burnett Exp $

"""
import os, sys, types
import numpy as np
from . import ( dataset, exposure, psf)
import skymaps
from uw.utilities import keyword_options
        

class Configuration(object):
    """
    Manage an analysis configuration: the IRF, exposure and data
    This is used by the band class EnergyBand to set up the PSF and exposure for the energy 
    range and sky position.
    """
    defaults =(
        ('event_type_names', ('front', 'back'), 'names to describe event classes'),
        ('bins_per_decade', 4, 'number of bins per decade'),
        ('irf', None,  'Set to override value in config.txt : may find in custom_irf_dir'),
        ('extended', None, 'Set to override value in config.txt'),
        ('nocreate', True, 'Do not allow creation of a binned photon file'),
        ('quiet', False, 'set to suppress most output'),
        ('postpone', False, 'set true to postpone loading data until requested'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, configdir='.', **kwargs):
        """ 
        parameters
        ----------
        configdir: folder containing configuration definition, file config.txt
        It is a Python dictionary, and must have the keys:
            irf : a text string name for the IRF
            diffuse: a dict defining the diffuse components
                see diffuse.DIffuseDict
            datadict: a dict with at least a key 'dataname'
            extended: string which is the path to the extended folder
                see extended.ExtendedCatalog
            
        If there is a key 'input_model' and the file 'pickle.zip' does not exist in the directory
        will look in input_model['path'] folder
        """
        keyword_options.process(self, kwargs)
        # make absolute path
        if configdir.startswith('$'):
            self.configdir = os.path.expandvars(configdir)
        elif not configdir.startswith('/'):
            self.configdir = os.path.join(os.getcwd(), configdir)
        else:
            self.configdir = configdir
        
        if not os.path.exists(self.configdir):
            raise Exception( 'Configuration folder %s not found' % self.configdir)
        if not os.path.exists(os.path.join(self.configdir, 'config.txt')):
            raise Exception('Configuration file "config.txt" not found in %s' % self.configdir)
        if not self.quiet: print 'Using configuration file "config.txt" in folder: %s' % self.configdir
        
        # extract parameters config
        try:
            config = eval(open(self.configdir+'/config.txt').read())
        except Exception, msg:
            raise Exception('Failed to evaluate config file: %s' % msg)
        for key in 'extended irf'.split():
            if self.__dict__.get(key, None) is None: 
                if not self.quiet and kwargs.get(key) is not None:
                    print '%s : override config.txt, "%s" => "%s"' %(key, kwargs.get(key,None), config.get(key,None))
                self.__dict__[key]=config.get(key, None)

        # set up IRF
        
        irf = config['irf']
        if 'CUSTOM_IRF_DIR' not in os.environ and os.path.exists(os.path.expandvars('$FERMI/custom_irfs')):
            os.environ['CUSTOM_IRF_DIR'] = os.path.expandvars('$FERMI/custom_irfs')
        custom_irf_dir = os.environ['CUSTOM_IRF_DIR']
        
        # define diffuse dictionary for constucting spatial models
        
        self.diffuse = config['diffuse']
        
        # set up dataset
        
        datadict = config.get('datadict',None)
        assert isinstance(datadict, dict), 'Did not find datadict dictionary in config.txt'
        self.dataset = dataset.DataSet(datadict['dataname'], interval=datadict.get('interval',None),
                nocreate = self.nocreate,
                irf = irf,
                quiet = self.quiet,
                postpone=self.postpone,
                )
        if not self.quiet:  print self.dataset
        
        # use dataset to extract psf and exposure, set up respective managers
        
        exposure_correction=datadict.pop('exposure_correction', None)        
        self.exposureman = exposure.ExposureManager(self.dataset, exposure_correction=exposure_correction)
        
        self.psfman = psf.PSFmanager(self.dataset)
        
        # check location of model
        
        input_model = config.get('input_model', None)
        self.modeldir = self.configdir
        if not os.path.exists(os.path.join(self.configdir, 'pickle.zip')) and input_model is not None:
            self.modeldir = os.path.expandvars(input_model['path'])
            try:
                self.modelname = input_model['file']
            except KeyError:
                self.modelname = 'pickle.zip'
            if not os.path.exists(self.modeldir):
                t = os.path.expandvars(os.path.join('$FERMI', self.modeldir))
                if not os.path.exists(t):
                    raise Exception('No source model file found in %s or %s' %(self.modeldir, t) )
                self.modeldir=t
        if not os.path.exists(os.path.join(self.modeldir, self.modelname)):
            print 'WARNING: pickle.zip not found in %s: no model to load' % self.modeldir
        elif not self.quiet:
            print 'Will load healpix sources from %s/%s' % (self.modeldir,self.modelname)
            
    def __repr__(self):
        return '%s.%s: %s' %(self.__module__, self.__class__.__name__, self.configdir)
        
    def get_bands_for_ROI(self, index):
        return self.get_bands(skymaps.Band(12).dir(index), minROI=5)
    
    def get_bands(self, roi_dir, **kwargs):
        """ return a list of dataset.ROIband objexts for an ROI
        """
        dset = self.dataset
        dset.load()
        band_kwargs = dict(emin=dset.emin, emax=dset.emax, minROI=dset.minROI, maxROI=dset.maxROI)
        band_kwargs.update(kwargs)
        radius = band_kwargs['minROI'] # fixed radius now
        bandlist = []
        for band in dset.dmap:
            emin,emax, event_type =  band.emin(), band.emax(), band.event_class()&5
            if (emin + 1) < band_kwargs['emin'] or (emax - 1) >band_kwargs['emax']: continue
            #print int(emin), event_class
            energy= np.sqrt(emin*emax)
            bandlist.append( ROIBand(band, self.psfman(event_type,energy), self.exposureman(event_type,energy), 
                roi_dir, radius))
        return np.asarray(bandlist)
        
    def select_event_type(self, which):
        """convenience function for selecting an event type index
        which : None, int, or string
            if None, return None to indicate all
        
        returns None or the index of the selected event type
        """
        if which is None or which=='all': return None
        etnames = self.event_type_names
        try:
            if type(which)==str:
                which = which.lower()
                return etnames.index(which)
            t = etnames[which]
            return which
        except Exception, msg:
            print 'Bad event type, "%s": %s\nMust be one of %s or a valid index' % (which, msg, etnames)
            raise

 
 
