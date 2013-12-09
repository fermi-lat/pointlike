"""  
 Setup the ROIband objects for an ROI
 
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/dataset.py,v 1.28 2013/11/08 04:30:05 burnett Exp $

    authors: T Burnett, M Kerr, J. Lande
"""
version='$Revision: 1.28 $'.split()[1]
import os, glob, types 
import cPickle as pickle
import numpy as np
import skymaps, pointlike #from the science tools
from uw.data import dataman, dssman
from uw.utilities import keyword_options
from uw.like import pycaldb

class DataSetError(Exception):pass

class DataSpecification(object):
    """ 
    """
    def __init__(self, folder, quiet=True, **data):
        """
        folder : string
            the path to the folder where the dictionary was found
        data : dict
            dictionary that contains entries sufficient to describe the data, and how to generate it
            
        Note: if interval found, modify the binfile and ltcube name to include
        """
        if folder=='.': folder = os.getcwd()
        if 'pickle' in data:
            print 'found pickle %s' % data['pickle']
            t = pickle.load(open(data['pickle']))
            #print t.__dict__
            self.__dict__.update(t.__dict__)
        else:
            for key in 'ft1files ft2files binfile ltcube'.split():
                if key in data and data[key] is not None:
                    data[key]=os.path.expandvars(data[key]) 
                    if not os.path.isabs(key):
                        data[key] = os.path.join(folder,data[key])
                    # need a check, but will fail if need to glob
                    #assert os.path.exists(data[key]), 'DataSpec: file %s not found' % data[key]
            interval = data.get('interval', None)
            if interval is not None:
                for name in ('binfile', 'ltcube'):
                    data[name]=data[name].replace('.fit', '_%s.fit'%interval)
            self.__dict__.update(data)
        if not quiet: print 'data spec:\n', str(self.__dict__)

    def __str__(self):
        return self.data_name
        
class Interval(dict):
    """ dictionary for looking up an interval
        Use if no intervals.py file found
        returns an MET tuple, e.g.
        Interval()['day01']
        
    """
    start = 239557417
    day = 24*60*60
    length= [day, 7*day, 365.25*day/12., 365.25*day]
    def __getitem__(self, x):
        for i, kind in enumerate(['day', 'week', 'month', 'year']):
            if not x.startswith(kind): continue
            n = int(x[len(kind):])
            s,w = Interval.start, Interval.length[i]
            return (s+(n-1)*w, s+n*w)
        raise KeyError(x)
        
class DataSet(dataman.DataSpec):
    """
    subclass of data.dataman.DataSpec that uses dicts to specify data, 
        and sets up instrument response and  ROI setup
        
    Note the __call__ function
    """
    defaults = dataman.DataSpec.defaults + (
        ' legacy key words that should be in DSS',
        ('zenithcut',100,'Maximum spacecraft pointing angle with respect to zenith to allow'),
        ('thetacut',66.4,'Cut on photon incidence angle'),
        ('use_weighted_livetime',False,'Use the weighted livetime'),
        ('event_class', 'source', 'either source or clean'),

        ('gti',    None, 'good time interval'),
        ('dss',    None, 'DSS keyword object'),
        ('version', None, ''),
        ('binner', None, ''),
        ('bins',  None,  ''),
        ('interval', None, 'Name for a specific time period: key in a dictionary to (start, stop)'),

        ' new feature',
        ('pickle', None, 'pickled dataspec'),
        ('legacy', False, 'set True to read old files'),
        
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
        #('phase_factor', 1.0, 'phase factor for pulsar phase analysis'),

        ('nocreate', False, 'Set True to suppress file creation'),
        ('postpone', False, 'Set true to postpone opening data until needed'),
        ('quiet', True, 'set False for some output'),
        ('verbose', False, 'more output'),)
        
    
    @keyword_options.decorate(defaults)
    def __init__(self, dataset_name, **kwargs):
        """
        Create a new DataSet object.

        dataset_name: an instance of DataSpecification with links to the FT1/FT2,
                            and/or binned data / livetime cube needed for analysis
                            (see docstring for that class) 
        """

        dataspec = self._process_dataset(dataset_name, quiet=kwargs.get('quiet',False),
            interval=kwargs.pop('interval',None)).__dict__
        dataspec.update(
                ft1=dataspec.pop('ft1files',None), 
                ft2=dataspec.pop('ft2files',None),
                )
                
        dataspec.update(kwargs)
        # Now invoke the superclass to actually load the data, which may involve creating the binfile and livetime cube
        super(DataSet,self).__init__(  **dataspec)
        assert self.irf is not None, 'irf was not specifed!'
        self.CALDBManager = pycaldb.CALDBManager(
                irf=self.irf, 
                psf_irf=self.psf_irf,
                CALDB=self.CALDB,
                custom_irf_dir=self.custom_irf_dir)
        if self.exposure_cube is None:
            self.lt = skymaps.LivetimeCube(self.ltcube,weighted=False) ###<< ok?
        else:
            self.lt = None
        if not self.postpone:
            self._load_binfile()
    

    def _process_dataset(self, dataset, interval=None, month=None, quiet=False):
        """ Parse the dataset as either a DataSpecification object, a dict, or a string lookup key.
            interval: string
                name of a 
            month: sub spec.
        """
        self.name= ''
        if hasattr(dataset,'binfile'): # dataset is DataSpecification instance
            return dataset
        if hasattr(dataset,'pop'): # dataset is a dict
            if 'data_name' not in dataset.keys():
                dataset['data_name'] = 'Custom Dataset %d'%id(dataset)
            dataman.DataSpec.datasets[id(dataset)] = dataset
            return dataman.DataSpec(id(dataset),month=month)
        # it is a string, check dictionary in ., then $FERMI/data
        self.name=dataset
        folders = ['.'] + glob.glob( os.path.join(os.path.expandvars('$FERMI'),'data'))
        for folder in folders :
            dict_file=os.path.join(folder, 'dataspec.py')
            if not os.path.exists(dict_file):
                continue
            try:
                ldict = eval(open(dict_file).read())
            except Exception, msg:
                raise DataSetError( 'Data dictionary file %s not valid: %s' % (dict_file, msg))
            if interval is not None:
                try:
                    pyfile = os.path.join(folder, 'intervals.py')
                    idict = eval(open(pyfile).read()) if os.path.exists(pyfile) else Interval()
                    gr = idict[interval]
                    gti_mask = skymaps.Gti([gr[0]], [gr[1]])
                    if True: #self.verbose: 
                        if not quiet: print 'apply gti mask %s, %s' %(interval, gti_mask)
                except Exception, msg:
                    raise DataSetError('Interval dictionary file %s, key %s, problem: %s'\
                                % (pyfile, interval, msg))
            else: gti_mask=None
            if dataset in ldict: 
                if not quiet: print 'Opening dataset %s from key in %s' % (dataset, dict_file)
                # translate event class name to appropriate bit
                ddict=ldict[dataset]
                #ddict['event_class_bit']= 4 ######FOR PASS8 now 
                ddict['event_class_bit']=dict(source=2, clean=3, extraclean=4)[ddict.get('event_class','source').lower()]
                return DataSpecification(folder,  interval=interval, gti_mask=gti_mask, **ddict)
        raise DataSetError('dataset name "%s" not found in %s' % (dataset, folders))

    def load(self):
        if self.postpone:
            self._load_binfile()
            self.postpone=False
 
    def _load_binfile(self):
        if not self.quiet: print 'loading binfile %s ...' % self.binfile ,
        self.dmap = skymaps.BinnedPhotonData(self.binfile)
        if not self.quiet: print 'found %d photons in %d bands, energies %.0f-%.0f MeV'\
                % (self.dmap.photonCount(),len(self.dmap), self.dmap[1].emin(), self.dmap[len(self.dmap)-1].emax())

        if self.verbose:
            self.dmap.info()
            print '---------------------'
            
    def info(self, out=None):
        """ formatted table of band contents """
        self.load()
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

    def dss_info(self, indent=''):   
        """ return a formatted table of the DSS keywords
        """
        s = indent+'DSS: %-15s  %-10s%-10s%-10s\n'% tuple('name value units ref'.split() )
        s+= indent
        s+= indent.join(['       %-15s  %-10s%-10s %s\n' %\
            (dss['TYP'],dss['VAL'],dss['UNI'],dss['REF'] ) for dss in self.dss])
        return s
    
    def __str__(self):
        """ Pretty print of cuts/data."""
        s = [] #collections.deque()
        s.append('Bins per decade: {0}'.format(self.binsperdec))
        def process_ft(label, files):
            if files is None: 
                s.append(label + '\tNone')
                return
            s.append(label)
            s.append(label)
            if len(files) < 10:
                s.append('\n\t'.join(files))
            else:
                s.append('\n\t'.join(files[:5]))
                s.append('...')
                s.append('\n\t'.join(files[-5:]))
        process_ft('FT1 files: ',self.ft1files)
        process_ft('FT2 files: ',self.ft2files)
        s.append('Binned data: {0}'.format(self.binfile))
        if not self.postpone:
            s.append('           :  %d photons, %d energy bands from %d to %d MeV'\
                    % (self.dmap.photonCount(), len(self.dmap), self.dmap[1].emin(), self.dmap[len(self.dmap)-1].emax()))
        else:
            s.append('            (load postponed)')
        s.append('Livetime cube: {0}'.format(self.ltcube))
        s.append(self.gti.__str__())
        s.append(self.dss_info())
        return 'dataset "%s":\n  '%self.name+'\n  '.join(s)
    def __repr__(self):
        return self.__str__()


def main(datadict=dict(dataname='P7_V4_SOURCE_4bpd'), analysis_kw=dict(irf='P7SOURCE_V6')):
    """ for testing: must define a dataset, and and at least an IRF"""
    dataset = DataSet(datadict['dataname'], **analysis_kw)
    return dataset
    
def validate(model_path='.', interval=None, nocreate=False, logfile='dataset.txt'):
    """
    validate the dataset for a model, as defined by the pipeline architecture
    
    model_path: string
        full path to the folder containing the model info, expect to find a file config.txt
        can contain environment variables with $ prefix
    
    interval: [string | None]
        name for an interval to use: override any specified in the config.txt
        
    nocreate: bool
        if False (default), try to create the binned photon data and livetime cubes
        if True, will return False if either file has to be generated
        
    logfile: string
        if not None, name of a file to write info
        
    returns True if the files exist or have been created OK
    """
    actual_path =  os.path.expandvars(model_path)
    assert os.path.exists(actual_path), 'Model path %s not found' % actual_path
    try:
        modelspec = eval(open(os.path.join(actual_path, 'config.txt')).read())
    except Exception, msg:
        raise DataSetError('could not evaluate config.txt: '+msg)
    
    datadict = modelspec['datadict']
    dataname = datadict['dataname']
    if interval is None:
        interval= datadict.get('interval', None)
    try:
        dset = DataSet(dataname, interval=interval, irf=modelspec['irf'], nocreate=nocreate)
        if logfile is not None: open(logfile, 'w').write(dset.__str__())
    except Exception, msg:
        print 'Failed: %s ' % msg
        raise 
        return False
    return True
    
