"""
 catalog stuff
 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/catalog.py,v 1.3 2010/04/15 05:21:05 burnett Exp $
 """

import sys, os, math
from skymaps import SkyDir, SpectralFunction
import numpy as np
import data
from uw.utilities import makerec

catalog_root = os.path.join(data.fermi_root,'catalog')
default_assoc =  'gll_psc11month_v4r4_flags_v4r3p8.fit' #'gll_psc11month_v1r2.fit'
default_assoc =  'gll_psc11month_v4r4_flags_v4r3p1_v4.fit' # after 14 Jan
default_assoc =  'gll_psc11month_v4r4_flags_v4r4p1.fit'  # final 1FGL

default_catalog = data.default_catalog # use the same one analyhsis

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

def atnf(filename=os.path.join(catalog_root, 'ATNF.txt')):
    ret= makerec.textrec(filename, quiet=True)
    assert(len(ret)>0)
    return ret

def bzcat(filename=os.path.join(catalog_root,'BZCAT.txt')):
    """ return a catalog with the BZCAT entries
        initial columns are name, ra, dec, z, mag, type_code
    """
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

def crates(filename=os.path.join(catalog_root,'crates_table5.txt')):
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

def cgrabs(filename=os.path.join(catalog_root,'cgrabs.txt')):
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

def correlate(cat1, cat2, min_diff=0):
    """ return an array with shape (len(cat1),2) 
        with the first entry the minimum distance from each entry in cat1 to cat2
        the second entry the cat2 index
        expect that each list has ra,dec property
        if mindiff is >0, can exclude duplicate
    """
    mindiff=np.zeros((len(cat1),2))
    cat2d = [SkyDir(s.ra,s.dec) for s in cat2]
    for i,a in enumerate(cat1):
        adir = SkyDir(a.ra, a.dec)
        mind = 99
        bmin = -1
        for j,b in enumerate(cat2d):
            d = math.degrees(adir.difference(b))
            if d< mind and d>min_diff: mind=d; bmin=j
        mindiff[i][0]=mind
        mindiff[i][1]=bmin
    return mindiff

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
            
       
class Blazars(object):
    """ identify sources as blazars
        Usage:
            b = Blazars()
            b.identify(source)
    """
    def __init__(self):
        self.cats = [bzcat(), cgrabs(), crates()]

    def check(self, source):
        """ return stuff 
        """
        ret = []
        for cat in self.cats:
            t= correlate([source], cat)[0]
            ret.append([t[0], cat[t[1]] ])
        return ret

    def identify(self, source, tol=0.25):
        """
            source: object that has ra,dec propery, or a string that can evaluate to ra,dec
        """
        class Source(object):
            def __init__(self,ra,dec):
                self.ra, self.dec= ra,dec
        if type(source)==type(''):
            ra,dec = parse_coords(source)
            source = Source(ra,dec)
        d = self.check(source)
        for r,s in d:
            if r<tol: return s,r
        return None

class Catalog(object):
    def __init__(self, cat_file=default_catalog, assoc=default_assoc, quiet=True):

        self.cat_file = cat_file
        self.sources = makerec.fitsrec(os.path.join(catalog_root, cat_file),quiet=quiet)
        if assoc is None:
            self.assoc=self.sources
        else:
            self.assoc = makerec.fitsrec(os.path.join(catalog_root,assoc ),quiet=quiet) #different version?
        self.asslist = self.assoc[self.assoc.id_name!='']
        self.assdirs = [SkyDir(ass.ra,ass.dec) for ass in self.asslist]

        self.catdirs = [SkyDir(source.ra,source.dec) for source in self.sources]

    def association(self, source):
        adir = SkyDir(source.ra, source.dec)
        diffs = np.array([adir.difference(bdir) for bdir in self.assdirs])
        
        mindif = math.degrees(diffs.min())
        #assert(0)
        if mindif>0.25: return None
        minsource = self.asslist[diffs==diffs.min()]
        name, ra, dec = minsource.id_name, minsource.id_ra[0], minsource.id_dec[0]
        print 'associated with %s: %s, at (%6.3f,%+6.3f) (%6.3f degrees away)' % (source.nickname.strip(), name, ra, dec, mindif)
        return name, SkyDir(ra, dec)

    def neighbor(self, source, maxdiff=2, mindiff=0):
        """
            mindiff>0 can prevent inclusion of a catalog source
        """
        if isinstance(source, SkyDir):
            adir = source
        else:
            adir = SkyDir(source.ra,source.dec)
        diffs = np.array([adir.difference(bdir) for bdir in self.catdirs])
        t = (diffs>math.radians(mindiff)) * (diffs< math.radians(maxdiff))
        if t.sum()==0 : return None
        #print 'found %d neighbors ' % t.sum()
        return [(s.nickname, SkyDir(s.ra,s.dec), math.degrees(adir.difference(SkyDir(s.ra,s.dec)))) for s in self.sources[t]]
        #i = np.arange(len(self.sources))[t][0]
        #s = self.sources[i]
        #return s.nickname, SkyDir(s.ra, s.dec)

    def select(self, name):
        if name.split()[0].lower()=='1fgl':
            ret = self.sources.source_name==name
        elif name[:5]=='1FGL_':
            ret = self.sources.source_name==name[:4]+' '+name[5:]
        else:
            ret = self.sources.nickname==name
        
        if ret.sum() !=1: 
            print '%s not found' % name
            return None
        return ret

    def dump(self, filename=None,  select=np.index_exp[:]):
        """ simple dump of name, position for refit
        """
        outfile = file(filename,'w') if filename is not None else None
        print >> outfile, '#%-19s%10s%10s'%('name','ra','dec')
        for s in self.sources[select]:
            print >> outfile, '%-20s%10.4f%10.4f' % (s.nickname, s.ra, s.dec)

    def mindiff(self):
        """ array of minimum distance to catalog neighbor
        """
        r=[]
        sdir = [SkyDir(s.ra,s.dec) for s in self.sources]
        for i in range(len(sdir)-1):
            z = 99
            a = sdir[i]
            for j in range(i+1, len(sdir)):
                d = a.difference(sdir[j])
                if d<z:
                   z = d
            r.append(math.degrees(z))
        return np.array(r)

    def search(self, coords, maxdiff=0.25, quiet=False):
        ra,dec = parse_coords(coords)
        if not quiet: print 'looking for (%f,%f)' % (ra, dec)
        q = self.neighbor(SkyDir(ra,dec), maxdiff=maxdiff)
        if q is None:
            if not quiet: print 'no source found'
            return None
        for s in q:
            if not quiet: print 'Found %s, at (%5.3f, %+5.3f), within %4.3f deg' % (s[0].strip(),s[1].ra(),s[1].dec(),s[2]) 
        return s[0]



def name_selection(cat, names):
    """ return an index expression into cat for those entries matching one of a list of names
    """
    ret= np.array([s.nickname.strip() in set(names) for s in cat.sources])
    assert(ret.sum()>=1)
    return ret

def file_selection(cat, filename):
    """ return an index expression into cat for those entries matching the first entry(name) in a file
    """
    names = [line.split()[0] for line in open(filename)]
    return name_selection(cat, names)

def search(cat, coords):
    ra,dec = parse_coords(coords)
    print 'looking for (%f,%f)' % (ra, dec)
    q = cat.neighbor(SkyDir(ra,dec), maxdiff=0.25)
    if q is None:
        print 'no source found'
        return
    for s in q:
        print 'Found %s, at (%5.3f, %+5.3f), within %4.3f deg' % (s[0].strip(),s[1].ra(),s[1].dec(),s[2]) 

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
    *   Jxxxx.y+zzzz: 1FGL pattern
    *   name ra dec : just return as is
    *   name : look up in alias_list: if found, continue
             : look up in location_list: if found, return ra dec found there
    """
    t = spec.split()
    if len(t)==2 and t[0]=='1FGL' or len(t)==1 and t[0]=='J' and len(spec)==12:
        sources = makerec.fitsrec(os.path.join(catalog_root, default_catalog), quiet=True)
        q = sources.source_name==spec if len(t)==2 else '1FGL '+spec
        if q.sum()==0: 
            raise Exception('%s not found in catalog %s' %(spec, default_catalog))
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

if __name__=='__main__':
    pass #doit(catalog)
