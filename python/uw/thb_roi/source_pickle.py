"""
Manage the UW catalog source pickle storage

$Header$
"""
import os, pickle, glob, math
import numpy as np
from uw.utilities import makerec
from skymaps import SkyDir

def dump_pickle(self, name, outdir, **kwargs):
    """ Write a dictionary with parameters
    
        self: a ROIAnalysis object
        name: name for source, used as filename (unless fname in kwargs)
        outdir: ouput directory
        **kwargs: anything else to add to the dictionanry
    """
    name = name.strip()
    fname = kwargs.pop('fname', name)
    output = dict()
    output['name'] = name
    output['ra']   = self.roi_dir.ra()
    output['dec']  = self.roi_dir.dec()
    
    # get source fit parameters, relative uncertainty
    p,p_unc = self.psm.models[0].statistical()
    output['src_par'] = p 
    output['src_par_unc'] = p*p_unc 
    output['pivot_energy'] = self.psm.models[0].e0 #this is 1000 by default, but catalog will be different
    output['src_model'] = self.psm.models[0] # simpler to save it all
    output['bgm_par'] = np.hstack([m.statistical()[0] for m in self.bgm.models] )
    output['bgm_par_unc'] = np.hstack([m.statistical()[1] for m in self.bgm.models] )
    try:
        output['qform_par'] = self.qform.par if 'qform' in self.__dict__ else None
    except AttributeError:
        output['qform_par'] = None
    output['tsmax'] = None if 'tsmax' not in self.__dict__ else [self.tsmax.ra(),self.tsmax.dec()]
    output.update(kwargs) # add additional entries from kwargs

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    f = file(os.path.join(outdir,fname+'.pickle'),'wb')
    pickle.dump(output,f)
    f.close()


def get_class(adict):
    """Given association dictionary, decide what class to ascribe the source to.  Partly guesswork!
        original version by Eric Wallace
        added tev
    """
    if adict is None: return '   '
    cat_list=adict['cat']
    priority = '''bllac bzcat cgrabs crates crates_fom seyfert seyfert_rl qso agn 
                vcs galaxies pulsar_lat snr snr_ext pulsar_high pulsar_low pulsar_fom
                msp pwn hmxb lmxb globular tev ibis lbv dwarfs
               '''.split()
    others = ['ostar']
    ass_class = ['bzb','bzcat']+['bzq']*3+['agn']*6+['LAT psr']+\
                ['snr']*2 + ['psr']*4 + ['pwn'] + ['hmxb'] + ['lmxb']+ ['glc'] +['tev'] + 3*['None']
    cls = None
    for c,a in zip(priority,ass_class):
        if c in cat_list:
            cls = a
            break
    if cls is None:
        if cls not in others: print 'warning: %s not recognized' % cat_list
        return '   '
    if cls == 'bzcat': #special, use first 3 chars of source name
        cls = adict['name'][cat_list.index(c)][:3].lower()
    return cls

class SourceRecarray(object):
    """ manage creation of a recarray of sources"""
    def __init__(self, procs):
        """get the set of processors, column names from each, and start the recarray"""
        self.procs  = procs
        self.colnames = []
        for proc in self.procs:
            self.colnames += proc.colnames
        self.rec =makerec.RecArray(self.colnames)
    def process(self,pk):
        """pass a dictionary to each of the processors, collect items"""
        data = []
        for proc in self.procs:
            data  += proc(pk)
        assert len(data)== len(self.colnames), 'data length inconsistent with number of columns'
        self.rec.append(*data)
    def __call__(self):
        return self.rec()
    
class Standard(object):
    """ define 'standard' set of items to include in the source recarray
    """
    def __init__(self):
        """ ctor defines column names"""
        self.colnames ="""name ra dec pivot_energy 
            pnorm pnorm_unc pindex pindex_unc cutoff cutoff_unc cc
            ts ts2 band_ts galnorm galindex isonorm
            """.split() 

    def __call__(self, p):
        name, ra, dec   = '%-20s' %p['name'], p['ra'], p['dec']
        # powerlaw fits from pointlike
        pivot_energy = p['pivot_energy'] if 'pivot_energy' in p else 1000.
        src_par, src_par_unc = p['src_par'], p['src_par_unc'] 
        pnorm,pindex  = src_par[:2]
        if 'src_par_unc' in p:
            punc = src_par_unc[:2] 
        else: punc = None
        pnorm_unc, pindex_unc = punc if punc is not None else (0,0)
        cc = 0 # correlation coefficient
        # assume that if there is a third entry in the source parameters, it is exponential cutoff
        cutoff, cutoff_unc = 2*(np.nan,) if len(src_par)==2 else (src_par[2], src_par_unc[2])
        if 'src_model' in p:
            src_model = p['src_model']
            #extract correlation coefficient
            V = src_model.cov_matrix
            cc = V[0,1]/np.sqrt(V[0,0]*V[1,1])
        ts = p['ts']
        ts2 = ts if 'ts2' not in p else p['ts2'] # after fit maybe
        band_ts= p['band_ts']
        galnorm, galindex, isonorm = p['bgm_par']

        return [name, ra,dec, 
            pivot_energy, pnorm, pnorm_unc, pindex, pindex_unc, cutoff, cutoff_unc, cc,
            ts, ts2, band_ts, 
            galnorm, galindex, isonorm, 
            ]

class LocalizationInfo(object):
    """ adds localization info to recarray"""
    def __init__(self):
        self.colnames = 'delta_ts psig dcp qual fit_ra fit_dec a b ang'.split()
    def __call__(self,p):
        delta_ts=p['delta_ts']
        qfp = p['qform_par']
        good = True
        if qfp:
            psig = math.sqrt(qfp[3]*qfp[4])
            ra,dec = p['ra'], p['dec']
            dcp = math.degrees(SkyDir(ra,dec).difference(SkyDir(qfp[0],qfp[1])))
            qual = qfp[6]
            a,b,ang = qfp[3:6]
            fit_ra,fit_dec=qfp[0:2]
        else: 
            psig =dcp= 1    
            a,b,ang= 0,0,0
            fit_ra,fit_dec,qual=999,999,99
            good = False
        if psig>1: 
            psig=1
            good = False
        return (delta_ts, psig, dcp, qual,  fit_ra,fit_dec, a,b,ang,)

class AssociationInfo(object):
    """ adds association info to recarray"""
    def __init__(self):
        self.colnames = 'id_prob aclass'.split()
    def __call__(self, p):
        adict = p.get('adict', None)
        if adict is not None:
            id_prob = adict['prob'][0]
        else:  id_prob =p.get('id_prob', 1e-6) # kluge to make it real
        aclass = '%-7s'%get_class(adict)
        return ( id_prob, aclass,)
    
class FitQualityInfo(object):
    """ manage fit quality items in the rec array"""
    def __init__(self):
        """ ctor appends names to column names"""
        self.colnames='fitq_low fitq_high'.split()
    def __call__(self, pk):
        """ append values to pars tuple from pk dictionary """
        bi = pk['band_info']
        # merge the fit_ts array
        fts = bi['fit_ts']
        fit_ts =fts[::2]+fts[1::2] # fit TS per energy band, assume alternate front, back  
        #dts = bi['band_ts'] -bi['fit_ts']
        dts = bi['ts'] - fit_ts 
        return (dts[:4].sum(), dts[4:].sum()) 
        
class VariabilityInfo(object):
    def __init__(self):
        self.colnames = 'sigma delta prob cprob varindex'.split()
    def __call__(self, p):
        try:
            lc = p['variability']
            pars += (lc['sigma'], lc['delta'], lc['prob'], lc['cprob'], lc['varindex'])
        except KeyError:
            print 'key "variability" not found in file for %s' % p['name']
            pars += 5*(0,)
        return pars    

class HighLowInfo(SourceRecarray):
    def __init__(self):
        colnames+="low_ts high_ts low_counts high_counts low_signal high_signal signal_1GeV".split()
    def __call__(self, p):
        ## todo: restore this if want these guys again
        return (low_ts, high_ts, low_counts, high_counts, low_signal, high_signal, signal_1GeV,)


class OtherKeys(SourceRecarray):
    def __init__(self, keys):
        self.keys = keys
        self.colnames =  keys
    def __call__(self, p):
        return [p[key] for key in self.keys]
            
        
def load_rec_from_pickles(outdir, other_keys=None, **kwargs):
    """
    create recarray from list of pickled sources 
    outdir -- folder for analysis results, expect to find "pickle" below it
    other_keys -- a list of keys to include 
    """
    failed=0
    maxfail = kwargs.pop('maxfail',10)
    ignore_exception = kwargs.pop('ignore_exception',False)
    assert os.path.exists(os.path.join(outdir,'pickle')), 'pickle folder not found under %s' %outdir
    filelist = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
    assert len(filelist)>0, 'no .pickle files found in %s/pickle' % outdir 
    filelist.sort()

    procs = [Standard(), LocalizationInfo(), AssociationInfo(), FitQualityInfo()]
    if kwargs.pop('high_low', False): procs.append(HighLowInfo())
    if kwargs.pop('variability',False) : procs.append(VariabilityInfo())
    if other_keys is not None: procs.append(OtherKeys(other_keys))                        
    rec = SourceRecarray(procs)
    
    for fname in filelist:
        #print 'loading %s...' % fname,
        try:
            p = pickle.load(open(fname))
            rec.process(p)
        except Exception, arg:
            print 'Failed to load file  %s: %s' % (fname, arg)
            if not ignore_exception or failed>maxfail: raise
            failed +=1
    print 'read %d entries from %s (%d failed)' % (len(filelist),outdir,failed)
    return rec()
    
if __name__=='__main__':
    pass
    
