"""  
 Setup the ROIband objects for an ROI
 
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/dataset.py,v 1.5 2011/10/20 21:41:27 burnett Exp $

    authors: T Burnett, M Kerr, J. Lande
"""
version='$Revision: 1.5 $'.split()[1]
import os, glob, types
import numpy as np
from uw.like import pixeldata, pypsf, pycaldb, pointspec_helpers, roi_bands
from uw.utilities import keyword_options
from .. pipeline import dataspec

class DataSet(object):
    """ 
    Manage the data, producing a set of ROIBand objects with __call__
    """

    defaults = (
        'keywords controlling data binning and livetime calculation',
        ('exp_radius',20,'radius (deg) to use if calculate exposure or a ROI. (180 for full sky)'),
        ('zenithcut',105,'Maximum spacecraft pointing angle with respect to zenith to allow'),
        ('thetacut',66.4,'Cut on photon incidence angle'),
        ('event_class',3,'select class level (3 - diffuse; 2 - source; 1 - transient; 0 - Monte Carlo)'),
        ('conv_type',-1,'select conversion type (0 - front; 1 - back; -1 = front + back)'),
        ('tstart',0,'Default no cut on time; otherwise, cut on MET > tstart'),
        ('tstop',0,'Default no cut on time; otherwise, cut on MET < tstop'),
        ('recalcgti',False,'if True, try to get GTI from GT1 files; otherwise, try from livetime cube or binned data file'),
        ('binsperdec',4,'energy binning granularity when binning FT1'),
        ('emin',100,'Minimum energy'),
        ('emax',1e6,'Maximum energy'),
        ('use_weighted_livetime',False,'Use the weighted livetime'),
        
        'keywords for monte carlo data',
        ('mc_src_id',-1,'set to select on MC_SRC_ID column in FT1'),
        ('mc_energy',False,'set True to use MC_ENERGY instead of ENERGY'),
        
        'keywords controlling instrument response',
        ('irf',None,'Which IRF to use'),
        ('psf_irf',None,'specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!'),
        ('CALDB',None,'override the CALDB specified by $CALDB.'),
        ('custom_irf_dir',None,'override the CUSTOM_IRF_DIR specified by the env. variable'),
        
        'keywords defining actual ROI setup',
        ('minROI', 7, 'minimum ROI radius'),
        ('maxROI', 7, 'maximum ROI radius'),
        ('phase_factor', 1.0, 'phase factor for pulsar phase analysis'),
        
        'miscelleneous',
        ('quiet',False,'Set True to suppress (some) output'),
        ('verbose',False,'More output'),
    )
    @keyword_options.decorate(defaults)
    def __init__(self, dataset_name, **kwargs):
        """
        Create a new DataSet object.

        data_specification: an instance of DataSpecification with links to the FT1/FT2,
                            and/or binned data / livetime cube needed for analysis
                            (see docstring for that class) """

        self.dataspec = self._process_dataset(dataset_name)

        self.__dict__.update(self.dataspec.__dict__)
        keyword_options.process(self, kwargs)
        assert self.irf is not None, 'Must provide the name of an IRF to use'
        self.CALDBManager = pycaldb.CALDBManager(irf=self.irf,psf_irf=self.psf_irf,
            CALDB=self.CALDB,custom_irf_dir=self.custom_irf_dir)

         #TODO -- sanity check that BinnedPhotonData agrees with analysis parameters
        self.pixeldata = pixeldata.PixelData(self.__dict__)
        self.exposure  = pointspec_helpers.ExposureManager(self)
        self.psf = pypsf.CALDBPsf(self.CALDBManager)
    
    def __str__(self):
        #s = 'Data selection, class %s\n' %self.__class__.__name__
        show = """data_name irf minROI maxROI emin emax""".split()
        d = {}
        for key in show:
            d[key] = self.__dict__.get(key, 'not in self.__dict__!')
        return str(d)
            #s += '\t%-20s: %s\n' %(key,
            #    self.__dict__[key] if key in self.__dict__.keys() else 'not in self.__dict__!')
        #return s

    def _process_dataset(self,dataset,month=None):
        """ Parse the dataset as either a DataSpecification object, a dict, or a string lookup key.
            month: sub spec.
        """
        if hasattr(dataset,'binfile'): # dataset is DataSpecification instance
            return dataset
        if hasattr(dataset,'pop'): # dataset is a dict
            if 'data_name' not in dataset.keys():
                dataset['data_name'] = 'Custom Dataset %d'%id(dataset)
            dataspec.DataSpec.datasets[id(dataset)] = dataset
            return dataspec.DataSpec(id(dataset),month=month)
        # it is a string, check dictionary in ., then $FERMI/data
        folders = ['.'] + glob.glob( os.path.join(os.path.expandvars('$FERMI'),'data'))
        for folder in folders :
            dict_file=os.path.join(folder, 'dataspec.py')
            if os.path.exists(dict_file):
                try:
                    ldict = eval(open(dict_file).read())
                except:
                    print 'Data dictionary file %s not valid' % ldict
                    raise
                if dataset in ldict: 
                    print 'found dataset %s in %s' % (dataset, folder)
                    return dataspec.DataSpecification(folder, **ldict[dataset])
        # not found: this is deprecated, leave for backwards consisency
        raise RuntimeError('dataset name %s not found in %s' % (dataset, folders))
        return dataspec.DataSpec(dataset,month=month)

    def __call__(self, roi_dir, **kwargs):
        """ return array of ROIBand objects for given direction, radius, emin, emax
        """
        self.bands = []
        band_kwargs= dict()
        band_kwargs.update(kwargs)
        for band in self.pixeldata.dmap:
            if (band.emin() + 1) >= self.emin and (band.emax() - 1) < self.emax:
                # note: pass self to ctor for minROI, maxROI, exposure, psf.band_psf
                self.bands.append(roi_bands.ROIBand(band, self, roi_dir, **band_kwargs))

        return np.asarray(self.bands)

    def info(self):
        self.pixeldata.dmap.info()
        

def main(dataname='2years_z100', irf='P7SOURCE_V6',**kwargs):
    """ for testing: must define a dataset, and and at least an IRF"""
    from uw.pipeline import dataspec
    ds = dataspec.DataSpec(dataname)
    analysis_kw = dict(irf=irf, minROI=7, maxROI=7)
    analysis_kw.update(kwargs)
    ds = DataSet(ds, **analysis_kw)
    print ds
    return ds
