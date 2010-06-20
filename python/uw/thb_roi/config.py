"""
Define an analysis environment for the UW pointlike ROI analysis
"""

import os, glob 
import skymaps 

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
    return os.path.join(data_path, *pars)

def gti_noGRB():
    """ return a Gti object appropriate to mask off the GRB intervals, using Gti.intersection()
    """
    grbs=[(243216749,243217979), # GRB 080916C
          (263607771,263625987), # GRB 090510
          (273582299,273586600), # GRB 090902B:
          ]
    g = skymaps.Gti()
    tbegin, tend = 0, 999999999
    for grb in grbs:
        g.insertInterval(tbegin, grb[0])
        tbegin = grb[1]
    g.insertInterval(tbegin, tend)
    return g

# basic data files
    
data_dict = {
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
    binfile     = data_join('18M', '18months_4.fits'),
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
    ft1files    = glob.glob(data_join('bpd','*.fits')),
    ft2files    = glob.glob(data_join('lt','*.fits')),
    binfile     = data_join('twenty','20month_4bpd.fits'),
    ltcube      = data_join('twenty','20month_lt.fits'),
    ),
}

# these configuration values are system, or time dependent
site_config = dict(
    catdir              = catalog_path,      # where to find catalog files
    catalog             = 'gll_psc_v03.fit', # the current catalog
    diffdir             = galprop_path,      # where to fine diffuse files
    aux_cat             = None,              # auxiallary catalog
)


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
    gti_mask     = None,
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
    )

class AE(object):
    """ Implement the AnalysisEnvironment interface, just a list of attribures
        kwargs: override only: items must already be present, except for "dataset", 
        which is expanded locally to ft1files, ft2files, binfile, ltcube
    """

    def __init__(self,  **kwargs):
        self.__dict__.update(system_config)
        self.__dict__.update(site_config)
        dataset = kwargs.pop('dataset', None)
        if dataset is not None: 
            try:
                self.__dict__.update(data_dict[dataset])
            except KeyError:
                raise KeyError, 'dataset "%s" is not in the data dictionary: %s' % (dataset, data_dict.keys())
        for key in kwargs.keys():
            if key not in self.__dict__:
                raise Exception, 'keyword "%s" is not recognized' %key
        self.__dict__.update(**kwargs)
        #if self.ft1files is None or self.ft2files is None:
        #    raise Exception, 'Expect both ft1files and ft2files to be set to lists of filenames'
        
    def __str__(self):
        s = self.__class__.__name__+'\n'
        for key in sorted(self.__dict__.keys()):
            v = self.__dict__[key]
            if key[:2]=='ft' and len(v)>3:
                s += '\t%-20s: %s ... %s\n' %(key, v[0], v[-1])
            else:
                s += '\t%-20s: %s\n' %(key, v)
        return s
    