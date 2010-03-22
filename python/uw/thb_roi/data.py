"""  manage Fermi data, set up for spectral analysis


  $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/data.py,v 1.2 2010/03/18 20:36:18 burnett Exp $
"""

import os, sys
import numpy as np
import skymaps


if os.name=='posix':
    fermi_root = '/phys/groups/tev/scratch1/users/Fermi'
elif os.environ['computername']=='GLAST-TS':
    fermi_root = 'f:\\glast'
else: # all others
    fermi_root = 'd:\\fermi'

# expect to find folders for ft1/ft2 data, galprop, and catalog info
data_path    = os.path.join(fermi_root, 'data')
galprop_path = os.path.join(data_path, 'galprop')
catalog_path = os.path.join(fermi_root, 'catalog')
for t in (data_path, galprop_path, catalog_path):
    if not os.path.exists(t):
        raise Exception('path does not exist: "%s"' % t)

# if CALDB is in the environment (set by science tools) use it. Otherwise expect under fermi_root
caldb_version = 'v1r1' # may have to set
if 'CALDB' not in os.environ:
    os.environ['CALDB'] = os.path.join(fermi_root, 'CALDB',caldb_version,'CALDB','data','glast','lat')
if not os.path.exists(os.environ['CALDB']):
    print 'bad CALDB: %s' %os.environ['CALDB']

#default_catalog ='gll_psc11month_v2.fit'   
#default_catalog ='gll_psc11months_combined_v4.fit' # after 11/24
#default_catalog ='gll_psc11month_v4r2_flags.fit' # after 12/04
#default_catalog ='gll_psc11month_v4r4_flags.fit' # after 12/14
#default_catalog ='gll_psc_v01.fit' # after 01/15
#default_catalog ='gll_psc_v02.fit' # after 01/29
default_catalog ='gll_psc_v03.fit' # after 02/10


class BasicData(object):
    def __init__(self,**kwargs):
        super(BasicData,self).__init__(**kwargs)

    def test(self):
        a,b,c,d = self()
        ok = True
        for x in a,b:
            if x is not None:
                for f in x: 
                    if os.path.exists(f): continue
                    print 'file %s not found' %f
                    ok = False
        for f in c,d:
            if f is not None:
                if os.path.exists(f): continue
                print 'file %s not found' %f
                ok = False
                
        return ok


class all_data(BasicData):
    def __init__(self, months=None, bins_per_decade=4, emax=1e5):
        self.bins_per_decade=bins_per_decade
        self.emax=emax  # not implemented
        if months is None: 
            self.months= open(os.path.join(data_path, 'flight','default_months.txt')).read().split()
        else:
            self.months = months
    def __call__(self):
        months = self.months
        event_files  = [os.path.join(data_path,'flight', '%s-ft1.fits'%m) for m in months ]
        history_files= [os.path.join(data_path, 'flight', '%s-ft2.fits'%m) for m in months ]
        livetimefile=os.path.join(data_path,'%s-%s_livetime.fits' %( months[0], months[-1]))
        datafile=os.path.join(data_path,   '%s-%s_%d.fits' %( months[0], months[-1],self. bins_per_decade ) )

        return event_files, history_files, datafile, livetimefile
    def __str__(self):
        return 'months %s to %s' % (self.months[0], self.months[-1])

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

class one_month(all_data):
    def __init__(self, month='aug2008',bins_per_decade=4, noGRB=False):
        self.bins_per_decade=bins_per_decade
        self.months=[month]
        self.noGRB=noGRB

    def __call__(self):
        months = self.months
        event_files  = [os.path.join(data_path,'flight', '%s-ft1.fits'%m) for m in months ]
        history_files= [os.path.join(data_path, 'flight', '%s-ft2.fits'%m) for m in months ]
        if  self.noGRB:
            datafile=os.path.join(data_path,   '%s_noGRB_%d.fits' %( months[0],self. bins_per_decade ) )
            livetimefile=os.path.join(data_path,'%s_noGRB_livetime.fits' %( months[0]))
            return event_files, history_files, datafile, livetimefile, gti_noGRB()
        else:
            datafile=os.path.join(data_path,   '%s_%d.fits' %( months[0],self. bins_per_decade ) )
            livetimefile=os.path.join(data_path,'%s_livetime.fits' %( months[0]))
            return event_files, history_files, datafile, livetimefile
    def __str__(self):
    
        return "one_month: %s (noGRB)" % self.month


class OneYear(BasicData):
    def __init__(self):
        self.bins_per_decade=4
        self.emax=1e6
    def __call__(self):
        months='aug2008 sep2008 oct2008 nov2008 dec2008 jan2009  feb2009 mar2009 apr2009 may2009 jun2009 jul2009'.split()
        event_files  = [os.path.join(data_path,'flight', '%s-ft1.fits'%m) for m in months ]
        history_files= [os.path.join(data_path, 'flight', '%s-ft2.fits'%m) for m in months ]
        livetimefile=os.path.join(data_path,'%s-%s_livetime.fits' %( months[0], months[-1]))
        datafile=os.path.join(data_path,   'one_year_%d.fits' % self.bins_per_decade )
        return event_files, history_files, datafile, livetimefile
    def __str__(self):
        return "one_year, 4 bins/decade to 1 TeV"




class Catalog(BasicData):
    def __init__(self):
        """ special selection, Aug 2008 to July 4 2009 (end of run 268411953)
        """
        self.bins_per_decade=4
        self.emax=1e6
    def __call__(self):
        months='aug2008 sep2008 oct2008 nov2008 dec2008 jan2009  feb2009 mar2009 apr2009 may2009 jun2009 jul2009_1to4'.split()
        event_files  = [os.path.join(data_path,'flight', '%s-ft1.fits'%m) for m in months ]
        history_files= [os.path.join(data_path, 'flight', '%s-ft2.fits'%m) for m in months ]
        livetimefile=os.path.join(data_path, 'catalog_livetimes.fits' )
        datafile=    os.path.join(data_path, 'catalog_%d.fits' %self.bins_per_decade )
        return event_files, history_files, datafile, livetimefile
    def __str__(self):
        return "eleven months, to July 4, 4 bins/decade to 1 TeV"

class Catalog_noGRB(Catalog):
    def __call__(self):
        months='aug2008 sep2008 oct2008 nov2008 dec2008 jan2009  feb2009 mar2009 apr2009 may2009 jun2009 jul2009_1to4'.split()
        event_files  = [os.path.join(data_path,'flight', '%s-ft1.fits'%m) for m in months ]
        history_files= [os.path.join(data_path, 'flight', '%s-ft2.fits'%m) for m in months ]
        livetimefile=os.path.join(data_path, 'catalog_noGRB_livetimes.fits' )
        datafile=    os.path.join(data_path, 'catalog_noGRB_%d.fits' %self.bins_per_decade )
        return event_files, history_files, datafile, livetimefile, gti_noGRB()
    def __str__(self):
        return "eleven months, to July 4, 4 bins/decade to 1 TeV, no GRB"

class FifteenMonths(BasicData):
    def __init__(self):
        self.bins_per_decade=4
        self.emax = 1e6
    def __call__(self):
        event_files =  [os.path.join(data_path, 'kerr', 'merged.fits')]
        history_files= [os.path.join(data_path, 'kerr', 'kerrm-megaFT2-ft2-30s.fits')]
        livetimefile =  os.path.join(data_path, 'kerr', 'fifteen_livetime.fits')
        datafile =      os.path.join(data_path, 'kerr', 'fifteen_4.fits')
        return event_files, history_files,datafile, livetimefile
    def __str__(self):
        return "fifteen months, 4 bins/decade to 1 TeV"

class MyAnalysisEnvironment():
    """
    replacement for uw.like.pointspec.AnalysisEnvironment
    but can be used to initialize a PixelData instance (in uw.like.pixeldata)
 
    """
    def __init__(self,   **kwargs):
        """
    Optional keyword arguments:

      =========    =======================================================
      Keyword      Description
      =========    =======================================================

      =========    KEYWORDS CONTROLLING DATA BINNING AND LIVETIME CALCULATION
      roi_dir      [ None] aperture center; if None, assume all-sky analysis                                                                                  
      exp_radius   [ 20]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)                                                                                
      zenithcut    [ 105]  Maximum spacecraft pointing angle with respect to zenith to allow
      thetacut     [ 66.4] Cut on photon incidence angle
      event_class  [ 3]  select class level (3 - diffuse; 2 - source; 1 - transient; 0 - Monte Carlo)
      conv_type    [ -1] select conversion type (0 - front; 1 - back; -1 = front + back)
      tstart       [0] Default no cut on time; otherwise, cut on MET > tstart
      tstop        [0] Default no cut on time; otherwise, cut on MET < tstop
      binsperdec   [4] energy binning granularity when binning FT1
      emin         [100] Minimum energy                                                                                 
      emax         [3e5] Maximum energy
      gti_mask     [None] If set, will mask the ft1 Gti

      =========    KEYWORDS FOR MONTE CARLO DATA
      mc_src_id    [ -1] set to select on MC_SRC_ID column in FT1
      mc_energy    [False] set True to use MC_ENERGY instead of ENERGY

      =========    KEYWORDS CONTROLLING INSTRUMENT RESPONSE
      irf          ['P6_v8_diff'] Which IRF to use
      psf_irf      [None] specify a different IRF to use for the PSF; must be in same format/location as typical IRF file!
      
    """
        # set defaults
        self.__dict__.update({
            # these allow this to initilaize a PixelData instance, 
            # perhaps generate pixel and exposure cube files
            'dataset'   : 'all_data',
            'roi_dir'   : None,
            'exp_radius': 180,
            'binsperdec':  4,
            'emin'      :  100,
            'emax'      :  3e5,
            'roi_dir'   :  None,
            'zenithcut' :  105,
            'thetacut'  :  66.4,
            'event_class': 3,
            'conv_type' :  -1,
            'tstart'    :   0,
            'tstop'     :   0,
            'mc_energy' :  False,
            'mc_src_id' :  -1,
            'gti_mask'  :  None,  
            
            # specifically affect analysis
            'diffdir'   :  galprop_path,
            'CALDB'     :  os.environ['CALDB'],
            'catdir'    :  catalog_path,
            'catalog'   :  default_catalog,
            'irf'       :  'P6_v8_diff',
            'psf_irf'   :  None,
            
            # others
            'quiet'     :  False,
            'verbose'   :  False,
            })
        self.__dict__.update(kwargs)
        mydata = self.__dict__.pop('dataset')
        if mydata is None: mydata = 'all_data'
       
        if type(mydata)==type(''):
            if mydata=='all_data': 
                args = all_data()()
            else:
                raise Exception('data set id %s not recognized' % mydata)
        else:
            args = mydata()
        event_files, history_files, datafile, livetimefile = args[:4]
        if len(args)==5: self.gti_mask= args[4]
        
        self.ft1files = event_files
        self.ft2files = history_files
        self.ltcube   = livetimefile
        self.binfile  = datafile

    def __str__(self):
        s = self.__class__.__name__+'\n'
        for key in sorted(self.__dict__.keys()):
            v = self.__dict__[key]
            if key[:2]=='ft' and len(v)>3:
                s += '\t%-20s: %s ... %s\n' %(key, v[0], v[-1])
            else:
                s += '\t%-20s: %s\n' %(key, v)
        return s


if __name__=='__main__':
    pass
