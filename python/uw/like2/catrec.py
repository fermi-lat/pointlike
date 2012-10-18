"""
Support for generating output files
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/catrec.py,v 1.6 2012/09/29 16:05:12 burnett Exp $
"""
import os, glob, types, zipfile
import cPickle as pickle
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
    if type(cat_list[0])==np.typeDict['int16']:
        # apparently from the catalog: don't know how to interpret this number yet, just pass it on: it will be a string
        # suggest  class_index = lambda x: 0 if x=='' else int(x)
        return cat_list[0]
        
    priority = '''bllac bzcat cgrabs crates crates_fom seyfert seyfert_rl qso agn 
                vcs galaxies pulsar_lat snr snr_ext pulsar_high pulsar_low pulsar_fom pulsar_nonATNF
                msp pwn hmxb lmxb globular tev ibis lbv dwarfs
               '''.split()
    others = ['ostar', 'starbursts', 'ocl']
    ass_class = ['bzb','bzcat']+['bzq']*3+['agn']*6+['LAT psr']+\
                ['snr']*2 + ['psr']*5 + ['pwn'] + ['hmxb'] + ['lmxb']+ ['glc'] +['tev'] + 3*['None']
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
    if os.path.exists(os.path.join(outdir,'pickle.zip')):
        pzip = zipfile.ZipFile(os.path.join(outdir,'pickle.zip'))
        filelist = ['pickle/HP12_%04d.pickle' %i for i in range(1728)]
        assert all(f in pzip.namelist() for f in filelist), 'Improper model zip file: expected 1728 pickle files'
        opener = pzip.open
    else:
        assert os.path.exists(os.path.join(outdir,'pickle')), 'pickle folder not found under %s' %outdir
        files = sortec(glob.glob(os.path.join(outdir, 'pickle', '*.pickle')))
        opener = open
 
    failed,maxfail = 0,kwargs.pop('maxfail',10)
    ignore_exception = kwargs.pop('ignore_exception',False)
    save_local = kwargs.pop('save_local',False) 
    ts_min = kwargs.pop('ts_min', 5)
    minflux = kwargs.pop('minflux', 1e-17)
    addtorec = kwargs.pop('addtorec', None)
    assert len(kwargs.keys())==0, 'unrecognized kw %s' %kwargs 
    if 'LATEXTDIR' not in os.environ:
        t = os.path.join(os.environ['FERMI'],'catalog','Extended_archive_v10') 
        assert os.path.exists(os.path.join(t,'Templates')), 'path %s not found' %t
        os.environ['LATEXTDIR']=t
    
    class CatalogRecArray(object):
        def __init__(self, minflux=minflux, update_position=False, ts_min=ts_min):
            self.minflux=minflux
            self.update_position = update_position
            self.count=self.reject=self.moved=0
            self.rejected=[]
            self.ts_min=ts_min
            self.moved = 0
            self.colnames ="""name ra dec 
                ts band_ts
                id_prob aclass id_ts
                hp12
                extended
                pnorm pindex cutoff 
                pnorm_unc pindex_unc cutoff_unc
                e0 pivot_energy 
                flux flux_unc
                eflux eflux_unc
                beta beta_unc
                index2 index2_unc
                modelname 
                fit_ra fit_dec a b ang qual delta_ts
                """.split() 
            if addtorec is not None:
                self.colnames += addtorec.colnames
            self.rec =makerec.RecArray(self.colnames) 

        def process(self,pk, cnan=np.nan):
            sources = pk['sources']
            for name,entry in sources.items():
                self.count+=1
                skydir = entry['skydir']
                data = [name, skydir.ra(), skydir.dec()]
                
                ts = entry['ts']
                data += [ts, entry['band_ts']]
                adict =  entry.get('associations', None)
                data.append( cnan if adict is None else adict['prob'][0])
                data.append('%-7s'%get_class(adict))
                data.append( cnan if adict is None else adict['deltats'][0])
                hp12 = Band(12).index(skydir)
                data.append(hp12)
                data.append(entry['isextended'])
                
                model  = entry['model']
                eflux = list(np.array(model.i_flux(e_weight=1, error=True, emax=1e5, quiet=True))*1e6)\
                                            if ts>10 else [0.,0.]
                #if np.any(np.isnan(eflux)):
                #    print 'bad energy integral, source %s [%d]' % (name, hp12)
                if np.isnan(eflux[0]):
                    import pdb; pdb.set_trace()
                p,p_relunc = model.statistical()
                if p[0]< self.minflux or np.any(np.isnan(p[:2])):
                    self.reject+=1
                    self.rejected.append(name)
                    continue
                p_unc = p*p_relunc
                psr_fit =  model.name.endswith('Cutoff')
                data += [p[0],     p[1],     p[2] if psr_fit else cnan, ]
                data += [p_unc[0], p_unc[1] ,p_unc[2] if psr_fit else cnan,]
                pivot_energy = entry.get('pivot_energy',model.e0)
                if pivot_energy=='None': pivot_energy=model.e0
                e0 = model.e0 if model.name!='LogParabola' else p[3]
                flux = model(e0)
                flux_unc = flux*p_relunc[0]
                data += [e0, pivot_energy]
                data += [flux, flux_unc]
                
                # energy flux from model e < 1e5, 1e-6 MeV units
                data += eflux
                if model.name=='ExpCutoff':
                    data += [cnan,cnan, 1.0, cnan, 'ExpCutoff']
                elif model.name=='PLSuperExpCutoff':
                    data += [cnan,cnan, p[3], p_unc[3], model.name]
                elif p[2]<0.01:
                    data += [0.0, cnan,    cnan, cnan, 'PowerLaw'] 
                else:
                    data += [p[2], p_unc[2], cnan, cnan, 'LogParabola']
                    
                
                #data += [sum(bts[4:]), sum(bts[8:])] ### note assumptions that 1, 10 GeV start at 4,8
                ellipse = entry.get('ellipse', None)
                if ellipse is None:
                    data += [np.nan]*7
                else:
                    # check for temporary mistake
                    if type(ellipse)==types.DictType:
                        data += [ellipse['ra'], ellipse['dec'],ellipse['a'],ellipse['b'],ellipse['ang'],ellipse['qual'], np.nan] #ellipse['delta_ts'] 
                    else: 
                        data += ellipse
                    if self.update_position and ts>self.ts_min:
                        fit_ra, fit_dec, a, b, ang, qual, delta_ts =ellipse
                        #assert False, 'break'
                        if qual<5 and delta_ts<9 and a < 0.2 :
                            data[1:3] = [fit_ra, fit_dec]
                            self.moved +=1
                if addtorec is not None:
                    data += addtorec(entry)
                assert len(data)==len(self.colnames), 'mismatch between names, data, %d, %d' %(len(self.colnames),len(data))
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
    
    
        def __init__(self,  nside=12):
            self.colnames = """name galnorm galindex isonorm  limbnorm
                                galnorm_unc galindex_unc isonorm_unc limbnorm_unc
                                backlimb backlimb_unc
                                loglike chisq glon glat ra dec""".split()
            self.rec = makerec.RecArray(self.colnames)
            self.nside=nside
            
        def process(self, pk):
            name = pk['name']
            diffuse_names = pk['diffuse_names']
            models = pk['diffuse'] 
            norm, norm_relunc = [np.hstack([m.statistical()[i] for m in models] ) for i in range(2)]
            p = norm[:3]
            p_unc = p*norm_relunc[:3] 
            no_limb = len(diffuse_names)==2 or (diffuse_names[2]).split('_')[0]!='limb'
            limbnorm, limbnorm_unc = [0,0] if no_limb else [norm[3], norm[3]*norm_relunc[3]]
            blimbnorm, blimbnorm_unc = [0,0] if no_limb else [norm[4], norm[4]*norm_relunc[4]]
            data = [name] + list(np.hstack((p,[limbnorm], p_unc, [limbnorm_unc])))
            data += [blimbnorm,  blimbnorm_unc]
            data += [pk['logl']]
            counts = pk['counts']
            obs,mod = counts['observed'], counts['total']
            data += [((obs-mod)**2/mod).sum()]
            skydir = pk['skydir']
            data += [skydir.l(), skydir.b(), skydir.ra(), skydir.dec()]
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
            p = pickle.load(opener(fname))
            crec.process(p)
            drec.process(p)
        except Exception, arg:
            print 'Failed to load file  %s: %s' % (fname, arg)
            if not ignore_exception or failed>maxfail: raise
            failed +=1
    print 'read %d entries from %s (%d failed)' % (len(filelist),outdir,failed)
    for r,name in ((crec, 'sources'), (drec, 'rois')):
        fname = '%s_%s.rec'%( outdir,name) if not save_local else '%s/%s.rec' % (outdir, name)
        rec = r()
        pickle.dump(rec, open(fname,'wb'))
        print 'saved %d entries to pickle file %s' % (len(rec), fname)
        
class ZipOpener(object):
    """
        manage a set of files associated with a model
    """
    def __init__(self, outdir, name='pickle'):
        self.name = name
        if os.path.exists(os.path.join(outdir, name+'.zip')):
            self.pzip = zipfile.ZipFile(os.path.join(outdir, name+'.zip'))
            self.opener = self.pzip.open
        else:
            assert os.path.exists(os.path.join(outdir,name)), 'folder %s not found under %s' %(name,outdir)
            #self.files = sorted(glob.glob(os.path.join(outdir, name, '*.pickle')))
            self.opener = open
    def __call__(self, fname):
        """ return the opened file, either from the folder or a zip of the same name"""
        return self.opener(os.path.join(self.name,fname))
        
class FitSource(dict):
    """ get flux information from the pickle """
    def __init__(self, versionfolder):
        self.pickle_opener = ZipOpener(versionfolder,'pickle')
    def __call__(self, src):
        """ for src, a recarray with hp12 and name entries, return the dictionary
        """
        p = pickle.load( self.pickle_opener('HP12_%04d.pickle'%src.hp12) )
        return p['sources'][src.name]
    
