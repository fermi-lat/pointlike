"""
Support for generating output files
$Header$
"""
import os, glob, pickle
import numpy as np
from uw.utilities import makerec
from skymaps import Band, SkyDir

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
    others = ['ostar', 'starbursts', 'ocl']
    ass_class = ['bzb','bzcat']+['bzq']*3+['agn']*6+['LAT psr']+\
                ['snr']*2 + ['psr']*4 + ['pwn'] + ['hmxb'] + ['lmxb']+ ['glc'] +['tev'] + 3*['None']
    cls = None
    for c,a in zip(priority,ass_class):
        if c in cat_list:
            cls = a
            break
    if cls is None:
        if cat_list[0] not in others: print 'warning: ''%s'' not recognized' % cat_list
        return '   '
    if cls == 'bzcat': #special, use first 3 chars of source name
        cls = adict['name'][cat_list.index(c)][:3].lower()
    return cls
        
def create_catalog(outdir, **kwargs):
    """ make a catalog file in the current directory, containing updated values for all the 
        sources found in the pickle folder
            note that order of first 6 parameters is assumed above
        also make a diffuse parameter dictionary
    """
    assert os.path.exists(os.path.join(outdir,'pickle')), 'pickle folder not found under %s' %outdir
    filelist = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
    assert len(filelist)>0, 'no .pickle files found in %s/pickle' % outdir 
    failed,maxfail = 0,kwargs.pop('maxfail',10)
    ignore_exception = kwargs.pop('ignore_exception',False)
    save_local = kwargs.pop('save_local',False) 
    ts_min = kwargs.pop('ts_min', 5)
    assert len(kwargs.keys())==0, 'unrecognized kw %s' %kwargs 
    filelist.sort()
    
    class CatalogRecArray(object):
        def __init__(self, minflux=1e-16, update_position=False, ts_min=ts_min):
            self.minflux=minflux
            self.update_position = update_position
            self.count=self.reject=self.moved=0
            self.rejected=[]
            self.ts_min=ts_min
            self.moved = 0
            self.colnames ="""name ra dec 
                pnorm pindex cutoff 
                pnorm_unc pindex_unc cutoff_unc
                e0 pivot_energy 
                flux flux_unc
                beta beta_unc
                modelname 
                ts band_ts bts10
                fit_ra fit_dec a b ang qual delta_ts
                id_prob aclass hp12
                """.split() 
            self.rec =makerec.RecArray(self.colnames) 

        def process(self,pk, cnan=np.nan):
            sources = pk['sources']
            for name,entry in sources.items():
                self.count+=1
                skydir = entry['skydir']
                data = [name, skydir.ra(), skydir.dec()]
                model  = entry['model']
                p,p_relunc = model.statistical()
                if p[0]< self.minflux or np.any(np.isnan(p[:2])):
                    self.reject+=1
                    self.rejected.append(name)
                    continue
                p_unc = p*p_relunc
                psr_fit =  model.name=='ExpCutoff'
                data += [p[0],     p[1],     p[2] if psr_fit else cnan, ]
                data += [p_unc[0], p_unc[1] ,p_unc[2] if psr_fit else cnan,]
                pivot_energy = entry.get('pivot_energy',model.e0)
                if pivot_energy=='None': pivot_energy=model.e0
                e0 = model.e0 if model.name!='LogParabola' else p[3]
                flux = model(e0)
                flux_unc = flux*p_relunc[0]
                data += [e0, pivot_energy]
                data += [flux, flux_unc]
                if psr_fit:
                    data += [cnan,cnan, 'ExpCutoff']
                else:
                    data += [cnan,cnan, 'PowerLaw'] if p[2]<=0.01 else [p[2], p_unc[2], 'LogParabola']
                ts = entry['ts']
                data += [ts, entry['band_ts']]
                data += [sum(entry['band_info']['ts'][7:])] ### note assumption that 10 GeV starts at 7
                ellipse = entry.get('ellipse', None)
                if ellipse is None or np.any(np.isnan(ellipse)):
                    data += [np.nan]*7
                else:
                    data += ellipse 
                    if self.update_position and ts>self.ts_min:
                        fit_ra, fit_dec, a, b, ang, qual, delta_ts = ellipse
                        #assert False, 'break'
                        if qual<5 and delta_ts<9 and a < 0.2 :
                            data[1:3] = [fit_ra, fit_dec]
                            self.moved +=1
                adict =  entry.get('associations', None)
                data.append( cnan if adict is None else adict['prob'][0])
                data.append('%-7s'%get_class(adict))
                data.append(Band(12).index(skydir))

                assert len(data)==len(self.colnames), 'mismatch between names, data'
                #assert not np.any(np.isinf(data[1:])), 'found inf for source %s' % name 
                self.rec.append(*data)

        def __call__(self): 
            print 'processed %d sources, rejected %d' %(self.count, self.reject)
            if self.reject>0:
                print 'rejected: flux<%.1e ' %self.minflux, self.rejected
            if self.update_position:
                print '\tmoved %d sources according to localization' % self.moved
            return self.rec()
        
    class DiffuseRecArray(object):
    
    
        def __init__(self, ndiff=3, nside=12):
            self.ndiff=ndiff
            if ndiff==3:
                self.colnames = """name galnorm galindex isonorm 
                                galnorm_unc galindex_unc isonorm_unc 
                                 loglike chisq""".split()
            else:
                self.colnames = """name galnorm galindex isonorm  isonorm2
                                galnorm_unc galindex_unc isonorm_unc isonorm2_unc
                                loglike chisq""".split()
            self.rec = makerec.RecArray(self.colnames)
            self.nside=nside
            
        def process(self, pk):
            name = pk['name']
            p, p_relunc = [np.hstack([m.statistical()[i] for m in pk['diffuse']] ) for i in range(2)]
            if len(p)!=self.ndiff:
                msg = 'unexpected number of diffuse parameters, %d vs %d, processing %s' % (len(p),self.ndiff,name)
                #print msg
                p = p[:self.ndiff]; p_relunc = p_relunc[:self.ndiff]
            data = [name] + list(np.hstack((p, p*p_relunc)))
            data += [pk['logl']]
            counts = pk['counts']
            obs,mod = counts['observed'], counts['total']
            data += [((obs-mod)**2/mod).sum()]
            assert len(data)==len(self.colnames), 'mismatch between names, data'
            self.rec.append(*data)
            
        def __call__(self):
            t = self.rec()
            n = 12*self.nside**2
            if len(t)!=n: 
                q = np.array([ 'HP12_%04d' % i not in t.name for i in range(n)])
                msg  = 'pixel file missing entries %s' % np.arange(n)[q]
                print msg
                raise Exception, msg
            return t

        
    crec = CatalogRecArray(**kwargs)
    drec = DiffuseRecArray()
    for fname in filelist:
        try:
            p = pickle.load(open(fname))
            crec.process(p)
            drec.process(p)
        except Exception, arg:
            print 'Failed to load file  %s: %s' % (fname, arg)
            if not ignore_exception or failed>maxfail: raise
            failed +=1
    print 'read %d entries from %s (%d failed)' % (len(filelist),outdir,failed)
    for r,name in ((crec, 'sources'), (drec, 'rois')):
        fname = '%s_%s.rec'%(name, outdir) if not save_local else '%s/%s.rec' % (outdir, name)
        rec = r()
        pickle.dump(rec, open(fname,'wb'))
        print 'saved %d entries to pickle file %s' % (len(rec), fname)
        
