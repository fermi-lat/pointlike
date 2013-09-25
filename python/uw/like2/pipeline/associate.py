
"""
 Manage the catalog association tables
 $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/associate.py,v 1.1 2011/12/09 16:07:44 burnett Exp $
 author: T. Burnett <tburnett@uw.edu>
"""
import pyfits, os, glob
import cPickle as pickle
import numpy as np
from skymaps import SkyDir
from uw.utilities import makerec
from uw.like import srcid  # gtsrcid work-alike from Eric

def delta(ra1,dec1,ra2,dec2):
    """ return r and theta (assume close)
    """
    cosdec = np.cos(np.radians(dec1))
    ddec = dec2-dec1
    dra  = (ra2-ra1)*cosdec
    if dra<-180: dra+=360
    if dra> 180: dra-=360
    theta = np.degrees(np.arctan2(ddec, dra)) % 180 -90
    return np.sqrt(ddec**2 + dra**2), theta

class Association():
    def __init__(self, acat):
        """ 
        acat: catalog FITS file with associations
                  get default catalog association
        """
        self.acat = acat 
        self.hdu  = pyfits.open(os.path.expandvars(acat))
        self.id_cat= self.hdu[2].data
        self.sources=self.hdu[1].data
        self.n   = len(self.sources)
        self.classes = 'gtsrcid list'
        # the valid keys for associations (todo)
        #keys = 'PSR pwn glc snr mqo xrb bzb bzq bzu rdg agn gal sbg sol'.split() 
        
    def __str__(self):
        return 'class %s: perform associations, generated from %s' % (self.__class__.__name__, self.acat)
    def __getitem__(self, n):
        """
            return a recarray with the associations, or None for source with index n
        """
        ret = []
        s = self.sources
        num = s.ID_Number[n]
        if num==0: return None
        fields = 'name,cat,ang,prob,ra,dec'.split(',')
        nstring= s.ID_Name[n]
        name = np.asarray([nstring[k:k+20].strip() for k in range(0,len(nstring),20)])
        m = len(name)
        ang = s.ID_Angsep[n][:m]
        prob= s.ID_Probability[n][:m]
        cat = s.ID_Catalog[n][:m]
        ra =  s.ID_RA[n][:m]
        dec=  s.ID_DEC[n][:m]
        return np.rec.fromarrays([name,cat,ang,prob,ra,dec], names=fields)


    def __call__(self, name, sdir=None, ellipse=None): # dummy for compatibility
        """ return None or a dict
        """
        if name[0]=='J': name='2FGL '+name
        if name[-1]=='c': name=name[:-1] # drop qualifier character
        for fld in ('Source_Name', 'NickName'):
            check = self.sources.field(fld)==name
            if check.sum()==1: break
        if check.sum()!=1:    
            text = 'failed to find source %s in catalog %s' %( name, self.acat)
            print text
            return None #raise Exception(text)
        i = np.arange(self.n)[check][0]
        d = dict()
        r = self[i]
        if r is None: return None

        for key in r.dtype.names:
            d[key]= r[key]
        d['dir'] = [SkyDir(float(x),float(y)) for x,y in zip(r['ra'],r['dec']) ]
        return d

    def select_cat(self, cats):
        """ return a recarray for a particular catalog number
        """
        fieldnames = 'name nickname ra dec a95 b95 posang ts association id_ra id_dec angsep prob'.split()
        rec = makerec.RecArray(fieldnames)
        sl = self.sources
        for i in range(self.n):
            q = self[i]
            if q is None: continue
            find = q.cat==cat
            if find.sum()==0: continue
            j = np.arange(len(q))[find][0] # if more than one, use first
            a,b,pa,ts= sl.Conf_95_SemiMajor[i], sl.Conf_95_SemiMinor[i], sl.Conf_95_PosAng[i],  sl.Test_Statistic[i],
            rec.append(sl.Source_Name[i][5:], sl.NickName[i], sl.RA[i], sl.DEC[i], a, b,pa ,ts,
                q.name[j],q.ra[j],q.dec[j], q.ang[j], q.prob[j],
                )
        return rec()


class UnPickle():
    """ create dictionary of associations from folder of pickle files
    """
    def __init__(self, path):
        fl = glob.glob(os.path.join(path, '*.pickle'))
        assert( len(fl)>0 )
        self.id={}
        self.info={}
        for f in fl:
            p = pickle.load(open(f))
            name = p['name'] # will use iFGL name as a key
            id = p['id']  # get dictionary
            self.id[name] = id
            self.info[name]=p

        self.fields = self.__dict__.values()[0].keys() # get fields from first file
        self.cat   = catalog.Catalog().sources # get the actual catalog as well
        
    def __getitem__(self, i):
        return self.id.values()[i]

    def __len__(self):
        return len(self.id)

    def key(self,i):
        return self.id.keys()[i]


    def rec(self,i):
        """ return a rec array with association info for source i (or None)

        """
        if self[i] is None: return None

        # note that "dir" is a SkyDir, can't be put into a recarray
        try:
            return np.rec.fromarrays(self[i].values()[:7], names=self[i].keys()[:7])
        except:
            return None

    def fit_tuple(self, name):
        """ return a tuple with a95, b95, phi for source name (with or without 1FGL)
        """
        if name[:5]!='1FGL ': name = '1FGL '+name
        q = self.cat.source_name==name
        assert(q.sum()==1)
        t = self.cat[q]
        return t.conf_95_semimajor.item(), t.conf_95_semiminor.item(), t.conf_95_posang.item()

    def select_cat(self, cats):
        """ return a recarray for a particular catalog number, or set thereoff
        """
        fieldnames = 'name ra dec ts band_ts a b qual a95 b95 phi id_name cat id_ra id_dec ang  prob deltats'.split()
        rec = makerec.RecArray(fieldnames)
        sl = self
        for i in range(len(self)):
            q = self.rec(i)
            if q is None: continue
            if fieldnames is None: fieldnames= q.dtype.names
            for cat in cats:
                find = q.cat==cat
                if find.sum()==0: continue
                # found one.
                name = self.key(i) 
                p = self.info[name]
                j = np.arange(len(q))[find][0] # if more than one, use first
                qf = p['qform_par'] 
                fit = self.fit_tuple(name)
                rec.append(
                    name,  p['ra'], p['dec'], p['ts'], p['band_ts'], qf[3],qf[4], qf[6],
                    fit[0], fit[1], fit[2],
                    q.name[j], q.cat[j], q.ra[j],q.dec[j], q.ang[j],  q.prob[j], q.deltats[j]
                    )
                break
        return rec()

#### add stuff to use Eric's code 

class SrcId(srcid.SourceAssociation):
    """
    adapter to simplify source association
    """
    def __init__(self, catalog_path='$FERMI/catalog/', classes='all_but_gammas', quiet=True):
        """ 
        catalog_path : string
            path to the catalogs, expect to find srcid/classes under it
        clases: ['all' | 'all_but_gammas' | list ]
            list of classes to apply,
        
        """
        self.classes = classes
        catalog_path = os.path.expandvars(catalog_path)
        d = os.path.join(catalog_path, 'srcid', 'classes')
        q = glob.glob(os.path.join(d, '*.py'))
        assert len(q)>0, 'no association classes found in folder %s' % d
        allclasses =[os.path.split(u)[-1].split('.')[0] for u in q if '__init__' not in u]
        if classes=='all': 
            # special tag to really get everything
            q = glob.glob(os.path.join(d, '*.py'))
            self.classes = allclasses
        elif self.classes=='all_but_gammas':
            self.classes = ['agn', 'bllac', 'bzcat', 'cgrabs', 'crates', 'crates_fom', 'dwarfs', 
            'galaxies', 'globular', 'hmxb', 'ibis', 'lbv', 'lmxb',  'ocl', 'ostar', 
             #'pulsar_fom',
            'pulsar_lat', 'pulsar_big', #'msp', 'pulsar_high',  'pulsar_low', 'pulsar_nonATNF', 
            'pwn', 'qso', 'seyfert', 'seyfert_rl', 'snr', 'snr_ext', 'starbursts', 'tev']
        else:
            self.classes=classes
        for c in self.classes:
            if c not in allclasses:
                txt = 'class %s not in set classes: %s' % (c, allclasses)
                raise Exception(txt)
        super(SrcId, self).__init__(os.path.join(catalog_path, 'srcid'),quiet=quiet)
        self.class_list = self.classes # will be used by the id method
     
    def __str__(self):
        return 'SrcId(%s)' %self.classes
    #def id(self, pos, error):
    #    """ the format returned by Srcid:
    #        a dictionary with keys classes and values a dict[sourcename]
    #    """
    #    return super(SrcId,self).id(pos, error, **dict(cpt_class=self.classes))
        
    def __call__(self, name, pos, error):
        """ name: source name, ignored, for reference
            pos: a SkyDir object
            error: [float | tuple]
                a or tuple: (a,b,ang) or, (ra,dec,a,b,ang,)...
                a,b: 1-sigma elipse in deg; ang orientation degrees
            returns None, or a dictionary consistent with Association above. (elements sorted by prob.)
        """
        if not hasattr(error, '__iter__'):
            error = (error,error,0)
        else:
            if len(error)==7: 
                error = error[3:6] 
            assert len(error)==3, 'wrong length for error ellipse specification'
        source_ass = self.id(pos,error)
        # select first association per catalog, rearrange to sort on over-all prob.
        candidates = [(v[0][1],v[0][0],v[0][2],key) for key,v in source_ass.items() if v!={}]
        if len(candidates)==0: return None
        candidates.sort(reverse=True)
        # format as a dict to conform to Association above (TODO: a recarray would be cleaner)
        d = {
            'name':  [a[1] for a in candidates],
            'cat':   [a[3] for a in candidates],
            'dir':   [a[2] for a in candidates],
            'prob':  [a[0] for a in candidates],
            'ra':    [a[2].ra() for a in candidates],
            'dec':   [a[2].dec() for a in candidates],
            'ang':   [np.degrees(a[2].difference(pos)) for a in candidates],
            }
        # stuff depending on the catalog 
        cats = [self.catalogs[a[3]] for a in candidates]
        d['prior'] =   [cat.prior for cat in cats]
        d['density'] = [cat.local_density(pos) for cat in cats]
        #d['fom'] =     [cat.fom for cat in cats]
        return d
        
    def positional_likelihood(self, name, pos, error):
        a,b,ang = error


def run_srcid(r, classes=['agn','bzcat','cgrabs','crates']):
    """ 
    Run the srcid tool on a recarrry, return associations from given set catalogs
        r: recarry with columns name, ra, dec, a, b, ang 
        classes: list of catalog names
        return: dict with key=name, value = sorted list of associations
    
    """
    assoc = SrcId(classes)
    associations = {}

    for s in r:
        pos = SkyDir(s.ra, s.dec)
        error = (s.a, s.b, s.ang)
        associations[s.name]=assoc(pos, error)

    return associations    