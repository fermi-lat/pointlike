"""
Manage data specification

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/dataspec.py,v 1.9 2011/04/26 16:14:40 burnett Exp $

"""
import os, glob


class DataSpecification(object):
    """ stand-in for current pointspec.DataSpecification 
    """
    def __init__(self, folder, **data):
        """
        folder : string
            the path to the folder where the dictionary was found
        data : dict
            dictionary that contains entries sufficient to describe the data, and how to generate it
        """
        if folder=='.': folder = os.getcwd()
        for key in 'ft1files ft2files binfile ltcube'.split():
            if key in data and data[key] is not None:
                data[key]=os.path.expandvars(data[key]) 
                if not os.path.isabs(key):
                    data[key] = os.path.join(folder,data[key])
                # need a check, but will fail if need to glob
                #assert os.path.exists(data[key]), 'DataSpec: file %s not found' % data[key]
        self.__dict__.update(data)
        print 'data spec:\n', str(self.__dict__)

class DataSpec(object):
    """
    This needs to be made local to an installation
    """

    def data_join(*pars):
        return os.path.expandvars(os.path.join('$FERMI','data', *pars))
    datasets = {
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
        'P7_V4_SOURCE': dict(data_name='pass 7 2FGL data set',
            ft1files   = '$FERMI/data/P7_V4_SOURCE/pass7.3*.fits',
            ft2files   = '$FERMI/data/P7_V4_SOURCE/ft2_2years.fits',
            binfile    = '$FERMI/data/P7_V4_SOURCE/24M7_z100_t90_cl0_4bpd.fits',
            ltcube     = '$FERMI/data/P7_V4_SOURCE/ltcube_24m_pass7.4_source_z100_t90_cl0.fits',
            monthly_bpd= '$FERMI/data/P7_V4_SOURCE/monthly/bpd/*_4bpd.fits',
            monthly_lt = '$FERMI/data/P7_V4_SOURCE/monthly/lt/*_lt.fits'
            ),
       'P7_V4_SOURCE_11M': dict(data_name='pass 7 1FGL data set',
            ft1files   = '$FERMI/data/P7_V4_SOURCE/pass7.3*.fits',
            ft2files   = '$FERMI/data/P7_V4_SOURCE/ft2_2years.fits',
            binfile    = '$FERMI/data/P7_V4_SOURCE/11M7_4bpd.fits',
            ltcube     = '$FERMI/data/P7_V4_SOURCE/11M7_lt.fits',
            ),
        'skymodel_ISO_GAL_p7_2Y':  dict(data_name = 'Nicola simulation of 2 years, diffuse only',
            event_class=-1,
            ft1files    = ['/phys/groups/tev/scratch1/users/Fermi/data/diffuse_simulation/merged_files_ft1.fits'], #glob.glob(data_join('18M', 'obssim_v9r16p1-*_ft1.fits')),
            ft2files    = [], #data_join('diffuse_simulation', 'ft2_18months.fits')],
            binfile    = data_join('diffuse_simulation','2years_4bpd.fits'),
            ltcube     = data_join('monthly','2years_lt.fits'),
            ),
        }

    def __init__(self, lookup_key, month=None):
        """
        lookup_key: string
            specify a key
        """
        # basic data files: will expand here
        t= lookup_key.split('[') #check for format 'name [n]' where n is a division
        if len(t)>1:
            lookup_key = t[0]
            month = int(t[1][:-1])
        data = self.datasets[lookup_key].copy() # copy so changes not kept for subsequent calls
        for key in 'ft1files ft2files binfile ltcube'.split():
            if key in data and data[key] is not None:
                data[key]=os.path.expandvars(data[key]) 
                # need a check, but will fail if need to glob
                #assert os.path.exists(data[key]), 'DataSpec: file %s not found' % data[key]
        if month is not None:
            assert 'monthly_bpd' in data, 'data spec %s does not have subdivisions' %lookup_key
            data['binfile'] = sorted(glob.glob(os.path.expandvars(data['monthly_bpd'])))[month]
            data['ltcube'] = sorted(glob.glob(os.path.expandvars(data['monthly_lt'])))[month]
            data['data_name'] += ' [%d]'%month

        self.__dict__.update(data)

    def __str__(self):
        return self.data_name
