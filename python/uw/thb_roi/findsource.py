"""
$Header$
"""
import os
import numpy as np
from uw.thb_roi.catalog import catalog_root, default_catalog, atnf
from uw.utilities import makerec
from skymaps import SkyDir

alias_list = {'vela'    : 'PSRJ0835-4510',
              'geminga' : 'PSRJ0633+1746',
              'crab'    : 'PSRJ0534+2200',
              'cta1'    : 'PSRJ0007+7303',
              'ic443'   : 'EMS0425',
              }
location_list={ # from catalog somewhere. could also make an alias
            '3C66A':   '035.665048 +43.035500',
            '3C66B':   '035.797547 +42.992051',
            'MRK421':  '166.113808 +38.208833',
            '3C454.3': '343.490616 +16.148211',
            '3C273':   '187.27791  +02.052499', # famuous, first quasar
            'BLLac':   '330.68037  +42.277772', # first of this type, name used for subsequent
            'MRK501':  '253.467569 +39.760169', # BL Lac type
            }


def find_source(spec):
    """ 
    spec can be:
    *   1FGL J...   : lookup from 1FGL list
    *   Jxxxx.y+zzzz: 1FGL pattern at
    *   name ra dec : just return as is
    *   name : look up in alias_list: if found, continue
             : look up in location_list: if found, return ra dec found there
    *   PSR... look up in ATNF
    """
    if spec=='?': 
        print find_source.__doc__
        print 'alias list:' , alias_list.keys()
        print 'source names:' , location_list.keys()
        return 
    t = spec.split()
    if len(t)==2 and t[0]=='1FGL' or len(t)==1 and t[0]=='J' and len(spec)==12:
        sources = makerec.fitsrec(os.path.join(catalog_root, default_catalog), quiet=True)
        q = sources.source_name==spec if len(t)==2 else '1FGL '+spec
        if q.sum()==0: 
            raise Exception('"%s" not found in catalog %s' %(spec, default_catalog))
        s = sources[q][0]
        return t[1], s.ra, s.dec

    if len(t)==3: return t # assume name, ra, dec
    if len(t)==1:
        # a single name: first try some specific names
        name = t[0]
        if name.lower() in alias_list:
            name = alias_list[name.lower()]
        if name in location_list:
            ra,dec = parse_coords(location_list[name])
            return name,ra,dec
        if name[:3]=='PSR':
            z = atnf();
            q = z.jname==name[3:]
            if q.sum()==0:
                raise Exception('name %s not found in ATNF'%name)
            e = z[q][0]
            return (name, e.ra, e.dec)
        if name[:2]=='BZ':
            bzc = bzcat()
            q = bzc.name==name
            if q.sum()==1:
                e = bzc[q][0]
                return (name, e.ra, e.dec)
            else: print '%s starts with "BZ", not found in BZcat' %spec
        sources = makerec.fitsrec(os.path.join(catalog_root, default_catalog), quiet=True)
        q = sources.source_name==name
        if q.sum()==0: 
            raise Exception('"%s" not found in catalog %s' %(spec, default_catalog))
        s = sources[q][0]
        return name, s.ra, s.dec
    else:
        raise Exception('unrecognized: "%s", expect a name, or "name ra dec"' %spec)
