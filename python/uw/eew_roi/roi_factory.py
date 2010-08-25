import os
import datetime
from uw.utilities.fermitime import utc_to_met
from uw.like.pointspec import SavedData,SpectralAnalysis

class ROIFactory(SpectralAnalysis):
    """Subclass of SpectralAnalysis as an ROIFactory, using SavedData

    __init__ args:
         tstart: start of the time range of interest (MET)
         tstop:  end of the time range of interest (MET)
    __init__ kwargs:
         Any valid kwargs for SpectralAnalysis (see docstring for SpectralAnalysis.__init__)
         data_dir: Directory in which to look for saved data products (defaults to $DATA_DIR, if defined)
         diffdir:  Directory containing diffuse model files
         catdir:   Directory containing catalogs
         catalogs: Catalogs to use
    """
    defaults = dict( roi_dir     = None,
                     exp_radius  = 20,
                     zenithcut   = 105,
                     thetacut    = 66.4,
                     event_class = 3,
                     conv_type   = -1,
                     binsperdec  = 4,
                     tstart      = utc_to_met(2008,8,4),
                     tstop       = utc_to_met(*datetime.datetime.now().timetuple()[:3]),
                     recalcgti   = False,
                     emin        = 200,
                     emax        = 2e5,
                     use_weighted_livetime = False,
                     lt_phibins  = 0,
                     mc_src_id   = -1,
                     mc_energy   = False,
                     irf         = 'P6_v3_diff',
                     psf_irf     = None,
                     CALDB       = os.environ['CALDB'],
                     background  = '1FGL',
                     maxROI      = 10,
                     minROI      = 5,
                     quiet       = False,
                     verbose     = False,
                     data_dir = (os.environ['DATA_DIR'] if os.environ.has_key('DATA_DIR')
                                 else '/phys/groups/tev/scratch1/users/Fermi/data'),
                     diffdir = (os.environ['DIFF_DIR'] if os.environ.has_key('DIFF_DIR')
                                 else '/phys/groups/tev/scratch1/users/Fermi/data/galprop'),
                     catdir = (os.environ['CAT_DIR'] if os.environ.has_key('CAT_DIR')
                                 else '/phys/groups/tev/scratch1/users/Fermi/catalog/'),
                     catalogs = ['gll_psc_v03.fit'])

    def init(self,**kw):
        d = self.defaults.copy()
        for k,v in d.items():
            d[k] = kw.pop(k,v)
        if kw:
            print("ROIFactory unrecognized kwargs:",
                  '\n'.join(['\t%s'%k for k in kw.keys()]),
                  '\n')
        self.__dict__.update(d)
        self.data = SavedData(self.tstart,self.tstop,binsperdec = self.binsperdec,data_dir = self.data_dir)

    def __init__(self,**kw):
        self.init(**kw)
        SpectralAnalysis.__init__(self,self.data,**kw)

    def __call__(self,skydir,**kwargs):
        catalogs = kwargs.pop('catalogs',self.catalogs)
        diffdir = kwargs.pop('diffdir',self.diffdir)
        return self.roi(roi_dir = skydir,
                        catalogs = [os.path.join(self.catdir,cat) for cat in catalogs],
                        diffdir = diffdir,
                        **kwargs)

