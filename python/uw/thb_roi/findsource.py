"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/findsource.py,v 1.1 2010/04/23 04:13:40 burnett Exp $
"""
import os
import numpy as np
from uw.utilities import makerec
from uw.thb_roi import config
from skymaps import SkyDir

def prune(cat, mindist=0.1):
    """ cat: catalog to prune: must have ra, dec columns
        return a array of bools for entries that are unique
    """
    dirs = map(SkyDir,cat.ra, cat.dec)
    rmin = np.radians(mindist)
    ret = np.ones(len(cat), bool)
    for i in range(1,len(cat)):
        s = dirs[i]
        for j in range(0,i):
            if s.difference(dirs[j])< rmin: 
                ret[i]=False #here if found another source within mindist
                break
    return ret
    
def parse_coords(coords):
    """
    parse various ways to describe ra/dec
    """
    def exception():
        raise Exception('search pattern "%s" not recognized' % coords)

    if type(coords)==type(''):
        toks = coords.split()
        if len(toks)==2:
            ra,dec = (float(toks[0]), float(toks[1])) 
        elif len(toks)==1:
            ip,im = coords.find('+'), coords.find('-')
            if coords[0]=='J' and len(coords)==10 and (ip==5 or im==5):
                ra = float(coords[0:2])*15 + float(coords[2:4])/4
                dec= float(coords[5:7]) + float(coords[7:9])/60.
                if im>0: dec=-dec
            if ip>0:
                ra,dec=(float(coords[:ip]), float(coords[ip+1:]))
            elif im>0:
                ra,dec=(float(coords[:im]),-float(coords[im+1:]))
            else: exception()
        elif len(toks)==6:
            s = toks[3][0]
            if s=='+': sign=1.0
            elif s=='-': sign=-1
            else: exception() #require explicit sign
            p = np.asarray(toks, 'f')
            ra = (p[0]+p[1]/60+p[2]/3600)*15.
            dec= sign*(abs(p[3])+p[4]/60+p[5]/3600)
        else:
            exception()
    else:
        ra,dec=coords
    return ra,dec

def atnf(filename=None):
    if filename is None: filename = os.path.join(config.catalog_path, 'ATNF.txt')
    ret= makerec.textrec(filename, quiet=True)
    assert(len(ret)>0)
    return ret

def bzcat(filename=None):
    """ return a catalog with the BZCAT entries
        initial columns are name, ra, dec, z, mag, type_code
    """
    if filename is None: filename = os.path.join(config.catalog_path, 'BZCAT.txt')
    f = open(filename)
    names = []; ras=[]; decs=[]; zs=[]; mags=[]; types=[]
    for line in f:
        if line[0:2] != 'BZ': continue
        toks = line.split('&')
        names.append(toks[0].strip())
        ra,dec=parse_coords(toks[1]+' '+toks[2])
        ras.append(ra)
        decs.append(dec)
        zs.append(toks[3])
        mags.append(toks[4])
        types.append(toks[5])
    fnames = 'name ra dec z mag type'.split()
    r = np.rec.fromarrays([names, ras,decs,zs,mags,types], names=fnames)
    f.close()
    return r

def crates(filename=None):
    if filename is None: filename = os.path.join(config.catalog_path, 'crates_table5.txt')

    """ return a catalog with the CRATES entries

J000000-002157   116 -0.498 00 00 01.66 -00 22 10.0 V    89.5 00 00 01.66 -00 22 09.8   215.4 -0.490 P
0                                                             63          75 
    """
    f = open(filename)
    names = []; ras=[]; decs=[]
    for line in f:
        if line[0]!='J' or len(line)<86 or line[63]==' ': continue
        names.append(line[0:15])
        ra,dec = parse_coords(line[62:86])
        ras.append(ra)
        decs.append(dec)
    r = np.rec.fromarrays([names, ras, decs], names='name ra dec'.split())
    f.close()
    return r

def cgrabs(filename=None):
    if filename is None: filename = os.path.join(config.catalog_path, 'cgrabs.txt')

    return makerec.textrec(filename)

class Density(object):
    def __init__(self, cat, maxdist=4):
        """ functor to determine the density of sources within the given distance from pos (in deg)
            cat: recarray, with ra,dec
            maxdist: float (deg)
        """
        self.dirs = [SkyDir(float(x),float(y)) for x,y in zip(cat.ra,cat.dec)]
        self.maxdist = np.radians(maxdist)
        self.area = 2*np.pi* (1.-np.cos(self.maxdist))*(180/np.pi)**2

    def __call__(self, pos):
        dist = np.asarray([pos.difference(d) for d in self.dirs])
        count = (dist<self.maxdist).sum()
        return count/self.area

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


def find_source(spec, catalog=None):
    """ 
    spec can be:
    *   1FGL J...   : lookup from 1FGL list
    *   Jxxxx.y+zzzz: 1FGL pattern at
    *   name ra dec : just return as is
    *   name : look up in alias_list: if found, continue
             : look up in location_list: if found, return ra dec found there
    *   PSR... look up in ATNF
    """
    if catalog is None: catalog=config.catalog
    if spec=='?': 
        print find_source.__doc__
        print 'alias list:' , alias_list.keys()
        print 'source names:' , location_list.keys()
        return 
    t = spec.split()
    if len(t)==2 and t[0]=='1FGL' or len(t)==1 and t[0]=='J' and len(spec)==12:
        sources = makerec.fitsrec(os.path.join(config.catalog_path, catalog), quiet=True)
        q = sources.field(sources.dtype.names[0])==spec if len(t)==2 else '1FGL '+spec
        if q.sum()==0: 
            raise Exception('"%s" not found in catalog %s' %(spec, catalog))
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
        sources = makerec.fitsrec(os.path.join(config.catalog_path,catalog), quiet=True)
        q = sources.field(sources.dtype.names[0])==name
        if q.sum()==0: 
            raise Exception('"%s" not found in catalog %s' %(spec, catalog))
        s = sources[q][0]
        return name, s.ra, s.dec
    else:
        raise Exception('unrecognized: "%s", expect a name, or "name ra dec"' %spec)
