"""
Define an analysis environment for the UW pointlike ROI analysis
"""

import os, glob, types 
import skymaps 

fermi_root = None
data_path = None
galprop_path = None
catalog_path = None
default_catalog = 'gll_psc_v03.fit' # 1FGL release
default_catalog = 'gll_psc18month_uw11b.fits' # later

def setup():
    global fermi_root, data_path, galprop_path, catalog_path
    # UW local configuration
    if os.name=='posix':
        fermi_root = '/phys/groups/tev/scratch1/users/Fermi'
    elif os.environ['computername']=='GLAST-TS':
        fermi_root = 'f:\\glast'
    elif os.environ['computername']=='FERMI-TS':
        fermi_root = 'T:\\'
    else: # all others
        fermi_root = 'D:\\fermi'

    # expect to find folders for ft1/ft2 data, galprop, and catalog info
    data_path    = os.path.join(fermi_root, 'data')
    galprop_path = os.path.join(data_path, 'galprop')
    catalog_path = os.path.join(fermi_root, 'catalog')
    for t in (data_path, galprop_path, catalog_path):
        if not os.path.exists(t):
            raise Exception('path does not exist: "%s"' % t)
            

def data_join(*pars):
    if data_path is None: setup()
    return os.path.join(data_path, *pars)
    
def data_glob(*pars):
    """ run glob and sort"""
    t = glob.glob(data_join(*pars))
    t.sort()
    return t

def gti_noGRB():
    """ return a Gti object appropriate to mask off the GRB intervals, using Gti.intersection()
    (GRB 080916C: 243216749-243217979, 090510: 263607771-263625987, 090902B: 273582299-273586600, 090926A: 275631598-275632048).
    """
    grbs=[(243216749,243217979), # GRB 080916C
          (263607771,263625987), # GRB 090510
          (273582299,273586600), # GRB 090902B
          (275631598,275632048), # GRB 090926A
          ]
    g = skymaps.Gti()
    tbegin, tend = 0, 999999999
    for grb in grbs:
        g.insertInterval(tbegin, grb[0])
        tbegin = grb[1]
    g.insertInterval(tbegin, tend)
    return g


## basic system configuration defaults
system_config = dict(
    CALDB        = os.environ['CALDB'] if 'CALDB' in os.environ else None, 
    ft1files     = None,
    ft2files     = None,  
    binfile      = None,   
    ltcube       = None,  
    binsperdec   = 4,
    conv_type    = -1,
    emax         = 1000000.0,
    emin         = 100,
    fit_emin     = 175,
    fit_emax     = 600000.,
    event_class  = 3,
    exp_radius   = 180,
    fit_bg_first = False,
    free_radius  = 1.0,
    gtimask      = None,
    irf          = 'P6_v8_diff',
    prune_radius = 0.1,
    quiet        = False,
    roi_dir      = None,
    thetacut     = 66.4,
    tstart       = 0,
    tstop        = 0,
    use_gradient = True,
    verbose      = False,
    zenithcut    = 105,
    pulsar_dict  = None,
    )

class AE(object):
    """ Implement the AnalysisEnvironment interface, just a list of attribures
        kwargs: override only: items must already be present, except for "dataset", 
        which is expanded locally to ft1files, ft2files, binfile, ltcube
        
    """

    def __init__(self,  **kwargs):
        setup() # execute and check system paths, etc.
        self.__dict__.update(system_config)
        self.__dict__.update(self.site_config())
        default_catalog = self.catalog # save 
        dataset = kwargs.pop('dataset', None)
        if dataset is not None: 
            if type(dataset)==types.StringType:
                try:
                    self.__dict__.update(self.data_dict(dataset))
                except KeyError:
                    #raise KeyError, 'dataset "%s" is not in the data dictionary: %s' % (dataset, data_dict.keys())
                    raise KeyError, 'dataset "%s" is not in the data dictionary' % dataset
            elif type(dataset==types.TupleType) and len(dataset)==2:
                # assume a tuple with data, ltcube
                kwargs['binfile'], kwargs['ltcube'] = dataset
            else:
                raise Exception, 'dataset specification %s not recognized: expect string, or tuple of (binfile,ltcube)'%dataset
        for key in kwargs.keys():
            if key not in self.__dict__:
                raise Exception, 'keyword "%s" is not recognized' %key
        self.__dict__.update(**kwargs)
        if self.ft1files is None and self.binfile is None:
            raise Exception('Expect either or both of ft1files and binfile to be set')
        if type( self.ft1files)==types.StringType:
            self.ft1files = glob.glob(self.ft1files)
        if self.ft2files is None and self.ltcube is None:
            raise Exception( 'Expect either or both of ft2files and ltcube to be set')
        if type( self.ft2files)==types.StringType:
            self.ft2files = glob.glob(self.ft2files)
        
        
    def __str__(self):
        s = self.__class__.__name__+'\n'
        for key in sorted(self.__dict__.keys()):
            v = self.__dict__[key]
            if key[:2]=='ft' and v is not None and len(v)>3:
                s += '\t%-20s: %s ... %s\n' %(key, v[0], v[-1])
            else:
                s += '\t%-20s: %s\n' %(key, v)
        return s
        
    def data_dict(self, lookup):
        # basic data files
    
        return  {
        '1FGL': dict( data_name = '11 month 1FGL data set',
            ft1files  = [data_join('flight','%s-ft1.fits')%x for  x in 
                ('aug2008','sep2008','oct2008','nov2008','dec2008',
                 'jan2009','feb2009','mar2009','apr2009','may2009','jun2009','jul2009_1to4')],
            ft2files  = [data_join('flight','%s-ft1.fits')%x for  x in    
                ('aug2008','sep2008','oct2008','nov2008','dec2008',
                 'jan2009','feb2009','mar2009','apr2009','may2009','jun2009','jul2009_1to4')],
            binfile = data_join('catalog_noGRB_4.fits'),
            ltcube =  data_join('catalog_noGRB_livetimes.fits'),
            gtimask = gti_noGRB(),
            ),
        '18M': dict( data_name = '18 month data set for 18M catalog',
            ft1files    = glob.glob(data_join('kerr2', '18M_data','*_ft1.fits')),
            ft2files    = glob.glob(data_join('kerr2', '18M_data','*_ft2.fits')),
            binfile     = data_join('18M', '18months_4bpd.fits'),
            ltcube      = data_join('18M', '18month_livetime.fits'),
            gtimask     = gti_noGRB(),
          ),
        '18M_tev': dict( data_name = '18 month data set for 18M catalog, high energy',
            ft1files    = glob.glob(data_join('kerr2', '18M_data','*_ft1.fits')),
            ft2files    = glob.glob(data_join('kerr2', '18M_data','*_ft2.fits')),
            binfile     = data_join('18M', '18months_tev_4.fits'),
            ltcube      = data_join('18M', '18month_livetime.fits'),
            gtimask     = gti_noGRB(),
          ),
        '20months': dict(data_name = "twenty months, 4 bins/decade to 1 TeV",
            ft1files    = data_glob('monthly','bpd','*_4bpd.fits')[:20],
            ft2files    = data_glob('monthly','lt','*.fits')[:20],
            binfile     = data_join('twenty','20month_4bpd.fits'),
            ltcube      = data_join('twenty','20month_lt.fits'),
            ),
       '2years': dict(data_name = "two years 4 bins/decade to 1 TeV",
            ft1files    = data_glob('monthly','bpd','*_4bpd.fits')[:24],
            ft2files    = data_glob('monthly','lt','*.fits')[:24],
            binfile     = data_join('monthly','2years_4bpd.fits'),
            ltcube      = data_join('monthly','2years_lt.fits'),
            ),
        }[lookup]
        
    def site_config(self):
    # these configuration values are system, or time dependent
        return dict(
            catdir      = catalog_path,      # where to find catalog files
            catalog     = default_catalog,   # the current catalog
            diffuse     = (galprop_path, 
                        'ring_21month_v1.fits','isotrop_21month_v1.txt',),
            aux_cat     = None,              # auxiallary catalog
        )



    