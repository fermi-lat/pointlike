
"""
 Manage the catalog association tables
 $Header$
 author: T. Burnett <tburnett@uw.edu>
"""
import catalog, makerec # local
import numpy as np
from skymaps import SkyDir
import pyfits, os, pickle, glob

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
    def __init__(self, acat=None):
        """ acat: [None] catalog FITS file with associations
                  get default catalog association
        """
        self.acat = acat or catalog.default_assoc
        self.hdu  = pyfits.open(os.path.join(catalog.catalog_root, catalog.default_assoc))
        self.id_cat= self.hdu[2].data
        self.sources=self.hdu[1].data
        self.n   = len(self.sources)
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


    def __call__(self, name):
        """ return None or a dict
        """
        if name[0]=='J': name='1FGL '+name
        if name[-1]=='c': name=name[:-1] # drop qualifier character
        check = self.sources.Source_Name==name
        if check.sum()!=1:
            text = 'failed to find source %s in catalog %s' %( name, self.acat)
            print text
            raise Exception(text)
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

