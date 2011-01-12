"""
Manage data specification

$Header$

"""
import os, glob
class DataSpec(object):
    """
    This needs to be made local to an installation
    """
    def __init__(self, lookup_key):
        """
        lookup_key: string
            spec
        """
        # basic data files
        def data_join(*pars):
            return os.path.expandvars(os.path.join('$FERMI','data', *pars))

        self.__dict__.update( {
        '1FGL': dict( data_name = '11 month 1FGL data set',
            ft1files  = [data_join('flight','%s-ft1.fits')%x for  x in 
                ('aug2008','sep2008','oct2008','nov2008','dec2008',
                 'jan2009','feb2009','mar2009','apr2009','may2009','jun2009','jul2009_1to4')],
            ft2files  = [data_join('flight','%s-ft1.fits')%x for  x in    
                ('aug2008','sep2008','oct2008','nov2008','dec2008',
                 'jan2009','feb2009','mar2009','apr2009','may2009','jun2009','jul2009_1to4')],
            binfile = data_join('catalog_noGRB_4.fits'),
            ltcube =  data_join('catalog_noGRB_livetimes.fits'),
            #gtimask = gti_noGRB(),
            ),
        '18M': dict( data_name = '18 month data set for 18M catalog',
            ft1files    = glob.glob(data_join('kerr2', '18M_data','*_ft1.fits')),
            ft2files    = glob.glob(data_join('kerr2', '18M_data','*_ft2.fits')),
            binfile     = data_join('18M', '18months_4bpd.fits'),
            ltcube      = data_join('18M', '18month_livetime.fits'),
            #gtimask     = gti_noGRB(),
          ),
        '18M_tev': dict( data_name = '18 month data set for 18M catalog, high energy',
            ft1files    = glob.glob(data_join('kerr2', '18M_data','*_ft1.fits')),
            ft2files    = glob.glob(data_join('kerr2', '18M_data','*_ft2.fits')),
            binfile     = data_join('18M', '18months_tev_4.fits'),
            ltcube      = data_join('18M', '18month_livetime.fits'),
            #gtimask     = gti_noGRB(),
          ),
        '20months': dict(data_name = "twenty months, 4 bins/decade to 1 TeV",
            binfile     = data_join('twenty','20month_4bpd.fits'),
            ltcube      = data_join('twenty','20month_lt.fits'),
            ),
       '2years': dict(data_name = "two years 4 bins/decade to 1 TeV",
            ft1files    = None,
            ft2files    = None,
            binfile     = data_join('monthly','2years_4bpd.fits'),
            ltcube      = data_join('monthly','2years_lt.fits'),
            ),
      '2years_z100': dict(data_name = "two years 4 bins/decade to 1 TeV, z<100",
            binfile     = data_join('monthly_z100','2years_4bpd.fits'),
            ltcube      = data_join('monthly_z100','2years_lt.fits'),
            ),
      '24MP7source': dict(data_name = "two years 4 bins/decade to 1 TeV, p7source (kluged ltcube)",
            binfile     = data_join('pass7.4','monthly','24M7_4bpd.fits'),
            #ltcube      = data_join('pass7.4','monthly','24M7_lt.fits'),
            ltcube      = data_join('monthly_z100','2years_lt.fits'),
            ),
       'skymodel_ISO_GAL_p6v8_18m':  dict(data_name = 'Nicola simulation of 18 months, diffuse only',
            event_class=0,
            ft1files    = glob.glob(data_join('18M', 'obssim_v9r16p1-*_ft1.fits')),
            ft2files    = [data_join('diffuse_simulation', 'ft2_18months.fits')],
            binfile    = data_join('diffuse_simulation','isogal_p6v8_18m_4pbd.fits'),
            ltcube     = data_join('diffuse_simulation', 'isogal_p6v8_18m_lt.fits'),
            ),
       'skymodel_ISO_GAL_p6v8_18m_mcenergy':  dict(data_name = 'Nicola simulation of 18 months, diffuse only (mc energy)', 
            event_class=0,
            ft1files    = '$FERMI/data/diffuse_simulation/ft1/obssim_v9r16p1-*_ft1.fits',
            ft2files    = '$FERMI/data/diffuse_simulation/ft2_18months.fits',
            binfile     = '$FERMI/data/diffuse_simulation/isogal_p6v8_18m_mcenergy_4pbd.fits',
            ltcube      = '$FERMI/data/diffuse_simulation/isogal_p6v8_18m_lt.fits',
            mc_energy=True, 
            ),
        }[lookup_key]
        )

    def __str__(self):
        return self.data_name
