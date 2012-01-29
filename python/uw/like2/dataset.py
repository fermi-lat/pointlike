"""  
 Setup the ROIband objects for an ROI
 
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/dataset.py,v 1.9 2012/01/11 14:09:17 burnett Exp $

    authors: T Burnett, M Kerr, J. Lande
"""
version='$Revision: 1.9 $'.split()[1]
import os, glob, types 
import cPickle as pickle
import numpy as np
import skymaps
from ..like import pycaldb, roi_bands
from ..data import dataman
from ..data.dataman import DataSpec

from ..utilities import keyword_options

class DataSetError(Exception):pass

class DataSpecification(object):
    """ 
    """
    def __init__(self, folder, **data):
        """
        folder : string
            the path to the folder where the dictionary was found
        data : dict
            dictionary that contains entries sufficient to describe the data, and how to generate it
        """
        if folder=='.': folder = os.getcwd()
        if 'pickle' in data:
            print 'found pickle %s' % data['pickle']
            t = pickle.load(open(data['pickle']))
            print t.__dict__
            self.__dict__.update(data)
        else:
            for key in 'ft1files ft2files binfile ltcube'.split():
                if key in data and data[key] is not None:
                    data[key]=os.path.expandvars(data[key]) 
                    if not os.path.isabs(key):
                        data[key] = os.path.join(folder,data[key])
                    # need a check, but will fail if need to glob
                    #assert os.path.exists(data[key]), 'DataSpec: file %s not found' % data[key]
            self.__dict__.update(data)
        #print 'data spec:\n', str(self.__dict__)

    def __str__(self):
        return self.data_name
        

        
class DataSet(dataman.DataSpec):
    """
    subclass of data.dataman.DataSpec that uses dicts to specify data, 
        and sets up instrument response and  ROI setup
        
    Note the __call__ function
    """
    defaults = dataman.DataSpec.defaults + (
        ' legacy key words that should be in DSS',
        ('zenithcut',105,'Maximum spacecraft pointing angle with respect to zenith to allow'),
        ('thetacut',66.4,'Cut on photon incidence angle'),
        ('use_weighted_livetime',False,'Use the weighted livetime'),
        
        ' new feature',
        ('pickle', None, 'pickled dataspec'),
        ('legacy', True, 'set True to read old files'),
        
        'keywords controlling instrument response',
        ('irf',None,'Which IRF to use'),
        ('psf_irf',None,'specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!'),
        ('CALDB',None,'override the CALDB specified by $CALDB.'),
        ('custom_irf_dir',None,'override the CUSTOM_IRF_DIR specified by the env. variable'),
        
        'keywords defining actual ROI setup',
        ('emin',100,'Minimum energy for selected bands'),
        ('emax',1e6,'Maximum energy for selected bands'),
        ('minROI', 7, 'minimum ROI radius'),
        ('maxROI', 7, 'maximum ROI radius'),
        ('phase_factor', 1.0, 'phase factor for pulsar phase analysis'),

        ('verbose', False, 'more output'),)
    
    @keyword_options.decorate(defaults)
    def __init__(self, dataset_name, **kwargs):
        """
        Create a new DataSet object.

        dataset_name: an instance of DataSpecification with links to the FT1/FT2,
                            and/or binned data / livetime cube needed for analysis
                            (see docstring for that class) 
        """

        dataspec = self._process_dataset(dataset_name).__dict__
        dataspec.update(
                ft1=dataspec.pop('ft1files',None), 
                ft2=dataspec.pop('ft2files',None),
                #binsperdec=4,
                )
        dataspec.update(kwargs)
        super(DataSet,self).__init__( **dataspec)
        assert self.irf is not None, 'irf was not specifed!'
        self.CALDBManager = pycaldb.CALDBManager(
                irf=self.irf, 
                psf_irf=self.psf_irf,
                CALDB=self.CALDB,
                custom_irf_dir=self.custom_irf_dir)
        self.lt = skymaps.LivetimeCube(self.ltcube,weighted=False) ###<< ok?
        self._load_binfile()

    def _process_dataset(self,dataset,month=None):
        """ Parse the dataset as either a DataSpecification object, a dict, or a string lookup key.
            month: sub spec.
        """
        if hasattr(dataset,'binfile'): # dataset is DataSpecification instance
            return dataset
        if hasattr(dataset,'pop'): # dataset is a dict
            if 'data_name' not in dataset.keys():
                dataset['data_name'] = 'Custom Dataset %d'%id(dataset)
            dataman.DataSpec.datasets[id(dataset)] = dataset
            return dataman.DataSpec(id(dataset),month=month)
        # it is a string, check dictionary in ., then $FERMI/data
        folders = ['.'] + glob.glob( os.path.join(os.path.expandvars('$FERMI'),'data'))
        for folder in folders :
            dict_file=os.path.join(folder, 'dataspec.py')
            if os.path.exists(dict_file):
                try:
                    ldict = eval(open(dict_file).read())
                except Exception, msg:
                    raise DataSetError( 'Data dictionary file %s not valid: %s' % (dict_file, msg))
                if dataset in ldict: 
                    print 'found dataset %s in %s' % (dataset, folder)
                    return DataSpecification(folder, **ldict[dataset])
        # not found: this is deprecated, leave for backwards consisency
        raise DataSetError('dataset name %s not found in %s' % (dataset, folders))

    def __call__(self, psf, exposure, roi_dir, **kwargs):
        """ return array of ROIBand objects for given direction, radius, emin, emax
        
        Parameters
        ----------
        sa : object of a class
            must have psf and 
        """
        self.bands = []
        band_kwargs= dict()
        band_kwargs.update(kwargs)
        class ForBand(object):
            """ mix stuff that ROIBand wants """
            def __init__(self, minROI, maxROI):
                self.__dict__.update(psf=psf, exposure=exposure, 
                        minROI=minROI, maxROI=maxROI)
        forband = ForBand(self.minROI,self.maxROI)
        for band in self.dmap:
            if (band.emin() + 1) >= self.emin and (band.emax() - 1) < self.emax:
                self.bands.append(roi_bands.ROIBand(band, forband, roi_dir, **band_kwargs))
        return np.asarray(self.bands)

    def _load_binfile(self):
        if not self.quiet: print 'loading binfile %s ...' % self.binfile ,
        self.dmap = skymaps.BinnedPhotonData(self.binfile)
        if not self.quiet: print 'found %d bands, energies %.0f-%.0f MeV'\
                % (len(self.dmap), self.dmap[1].emin(), self.dmap[len(self.dmap)-1].emax())

        if self.verbose:
            self.datainfo()
            print '---------------------'
            
    def info(self, out=None):
        """ formatted table of band contents """
        print >>out, 'File: %s ' %self.binfile
        print >>out, '\n  index    emin      emax  class nside     photons'
        total = 0
        def bignum(n):
            t = '%9d' % n
            return '  '+' '.join([t[0:3],t[3:6],t[6:]])
        for i,band in enumerate(self.dmap):
            fmt = '%5d'+2*'%10d'+2*'%6d'+'%12s'
            print fmt % (i, round(band.emin()), round(band.emax()), 
                    band.event_class()&3, band.nside(), bignum(band.photons()))
            total += band.photons()
        print >>out, 'total%45s'% bignum(total)

        #self.dmap.info()


def main(datadict=dict(dataname='P7_V4_SOURCE_4bpd'), analysis_kw=dict(irf='P7SOURCE_V6')):
    """ for testing: must define a dataset, and and at least an IRF"""
    dataset = DataSet(datadict['dataname'], **analysis_kw)
    return dataset
