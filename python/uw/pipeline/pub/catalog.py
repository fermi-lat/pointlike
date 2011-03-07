"""
catalog management classes for the UW pipeline
$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/pipeline/catalog.py,v 1.6 2011/01/01 15:50:05 burnett Exp $

"""
import os, pickle, types, copy, glob
import numpy as np
import pyfits
from uw.utilities import makerec, keyword_options, xml_parsers
from uw.like import  Models, pointspec_helpers
from skymaps import SkyDir, Band


def write_reg_file(catrec, filename, color='green'):
    """ generate a 'reg' file from , write to outfile
    """
    have_ellipse = 'Conf_95_SemiMajor' in catrec.dtype.names
    out = open(filename, 'w') if filename is not None else None
    print >>out, "# Region file format: DS9 version 4.0 global color=%s" % color
    for i,s in enumerate(catrec):
        
        if i==len(catrec): break
        if have_ellipse and not np.isnan(s[8]):
            print >>out, "fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) #text={%s}" % \
                            (s[1],s[2],s[8],s[9],s[10],s[0])
        else:
            print >> out, "fk5; point(%.4f, %.4f) # point=cross text={%s}" %\
                            (s[1],s[2],s[0])
    if out is not None: out.close()

class CatRec(object):

    """ manage a recarry used for generating catalogs,
        mainly for updating
    
    """
    def __init__(self, inrec, update_file=None):
        """
        inrec: recarray or string
            if string, load it.
            must have fields needed for catalog:  name, ra, dec, pnorm, pindex, cutoff, band_ts, ts, pivot_energy, e0, 
            the sequence (pnorm, pindex, cutoff) is expected to be in order,  [3:6]
        update_file : None | filename
            if set, must be name of file to use for update
        """
        self.rec = inrec if type(inrec)!=types.StringType else pickle.load(open(inrec))
        print 'CatRec: loaded %d sources' %len(self.rec.name)
        if len(self.rec.name) != len(set(self.rec.name)):
            print 'names are not unique!'
            sn = set()
            for n in self.rec.name:
                if n in sn: print 'duplicate:', n
                sn.add(n)
        self.sdirs = None
        if update_file is not None:
            self.update(update_file)
        
    def update(self, filename):
        """ filename : string
                name of a file containing name, ra, dec fields, or name, glat,glon
        """
        uprecarray = makerec.load(filename)
        cnames = ' '.join(uprecarray.dtype.names[:3])
        gal = False
        if cnames=='name glon glat': gal = True
        elif cnames!='name ra dec':
            raise Exception, 'recarray %s not recognized: column names are %s' %(filename, uprecarray.dtype.names)
        for uprec in uprecarray:
            name = uprec.name.replace('_', ' ') if uprec.name[:4]!='SEED' else uprec.name
            if not gal:
                ra,dec = uprec.ra, uprec.dec 
            else:
                sdir = SkyDir(uprec.glon, uprec.glat, SkyDir.GALACTIC)
                ra,dec = sdir.ra(), sdir.dec()
            if name in self.rec.name:
                if ra<0:
                    self.remove(name)
                else:
                    self.modify_pos(name, ra, dec)
            else:
                self.append(name, ra, dec)
    
    def set_cutoff(self, name, val=1000):
        sel = self.select(name)
        self.rec.cutoff[sel]=val
        
    def __call__(self, name):
        """ return entry (single-element recarray) with name
        """
        sel = self.select(name)
        return self.rec[sel][0]
    
    def select(self, name):
        sel = self.rec.name==name
        assert sum(sel)>0, 'name %s not found' % name
        if sum(sel)>1:
            print 'name %s not unique' %name
        return sel
    
    def modify_pos(self, name, ra, dec):
        sel = self.select(name)
        oldra,olddec = (self.rec.ra[sel], self.rec.dec[sel])
        self.rec.ra[sel]  = ra
        self.rec.dec[sel] = dec
        print 'updated ra,dec for %-18s : %8.3f%8.3f -->  %8.3f%8.3f ' % (name, oldra, olddec ,ra,dec)
        
    def save(self, filename):
        pickle.dump(self.rec, open(filename, 'wb'))
    
    def make_blank(self, name, ra, dec):
        """ create a blank entry, blank name, nans except for cat defauts for rest
        """
        e = copy.copy(self.rec[-1]) #use last entry as template
        for i in range(1,len(e)-1): e[i]=np.nan
        e.name = name
        e.ra = ra
        e.dec = dec
        e.ts=15
        e.band_ts=25
        e.pnorm = 1e-13
        e.pindex = 2.5
        e.e0=1000
        e.pivot_energy=1000
        e.beta = 0.01
        e.hp12 = -1 # not set??
        return e
    
    def remove(self, name):
        sel = self.select(name)
        newrec = self.rec[-sel]
        self.rec = newrec
        print 'removed entry %s from catalog recarray' % name
    
    def append(self, name, ra, dec):
        """ add a new entry, assigning defaults.
        """
        assert name not in self.rec.name, 'name %s not unique' % name
        entry = self.make_blank(name,ra, dec)
        newrec = np.hstack( [self.rec, entry ] ) 
        self.rec= newrec.view(np.recarray)
        print 'added entry with name, ra, dec: %s %8.3f%8.3f' % tuple(entry)[:3]
    
    def rename_seeds(self, seed_names=['SEE'], root='24M'):
        """ rename entries starting with "SEED"
        """
        newrec=np.sort(self.rec, order='ra') # sorting makes a new copy
        last = root+'9999' #bad if need this
        n=1
        for i in range(len(newrec)):
            name = newrec.name[i]
            if name[:len(root)]==root:
                last = name
                n=1
            elif name[:3] in seed_names:
                t = '%s.%d'%(last,n)
                n+=1
                if t in newrec.name:
                    t += '.1' # a little klugy
                print 'renamed %s --> %s' % (name,t)
                newrec.name[i] = t
        self.rec = newrec
        
    def rename_from_list(self, filename):
        fin = open(filename).read().split('\n')
        for line in fin[0:]:
            if line=='' or line[0]=='#': continue
            tokens = line.split()
            if len(tokens)<2: continue
            if len(tokens)<3: 
                old, new = tokens[:2];  tag=''
            else:
                old, tag, new = tokens[:3]
            if tag=='1FGL' or tag=='':
                try:
                    sel = self.select(old)
                    newname = tag+' '+new
                    assert newname not in self.rec.name, '%s already in table' % newname
                    print 'renaming %-15s to %-15s' % (old, newname)
                    self.rec.name[sel]= newname
                except:
                    print 'failed to rename %s' % new
                    raise
            else:
                'unrecognized tag: %s' % tag

    def nearest(self, odir):
        """ look for closest value"""
        if self.sdirs is None:
            self.sdirs = map(SkyDir, self.rec.ra, self.rec.dec)
        dist = np.array(map( odir.difference, self.sdirs))
        mindist = dist.min()
        name = self.rec.name[dist==mindist][0]
        return name, np.degrees(mindist)
        
    def match(self, fgl):
        """ match entries with those from another catalog
        fgl : recarry containing name, ra, dec entries
        if it is a standard catalog, expect 'Source_Name', 'RA', 'DEC'
        """
        cnames = fgl.dtype.names
        snames = fgl.field(cnames[0])
        sdirs = map(SkyDir, np.array(fgl.field(cnames[1]),np.float64), 
                            np.array(fgl.field(cnames[2]),np.float64))
        return [self.nearest(s) for s in sdirs]

class LocalSource(object):

    def __init__(self,skydir,name, model=None,ellipse=None, free_parameters=True, ts=0):
        self.name    = name.strip()
        self.skydir = skydir
        if model is None:
            self.model  = Models.PowerLaw() 
        else:
            self.model = model
        self.model.free[:] = free_parameters
        self.ellipse = None
        self.ts = ts

class CatalogManager(object):
    """Read a  FITS-format catalog and use it to set up a source list for a ROI.
    
    """
    defaults = (
        ('prune_radius',  0.10,  'deg; in a merge, consider sources closer than this duplicates\n'
                                ' or to exclude all sources within a larger region for separate specification'),
        ('free_radius',   1,      'deg; sources within this distance have free spectral parameters'),
        ('min_flux',     2e-9,   'ph/cm2/s; minimum flux for sources beyond a certain proximity'),
        ('max_distance',  5,      'deg; distance inside which sources are returned regardless of flux'),
        ('min_ts',       25,     ''),
        ('quiet',         False,   ''),
        ('verbose',       False,   ''),
        )

    @keyword_options.decorate(defaults)
    def __init__(self,catalog_file, *args, **kwargs):
        """ Create the catalog: 
            catalog_file: a Fermi-LAT FITS format catalog
            """
        keyword_options.process(self, kwargs)
        if not self.quiet: 
            print 'creating a CatalogManager from %s...' %catalog_file
        cdata = pyfits.open(catalog_file)[1].data
        try:
            ts   = np.asarray(cdata.field('Test_Statistic'),float)
        except KeyError:
            ts   = np.asarray(cdata.field('Signif_Avg'))**2
        good = ts>self.min_ts
             
        ras  = np.asarray(cdata.field('RA'),float)[good]
        decs = np.asarray(cdata.field('DEC'),float)[good]
        pens = cdata.field('PIVOT_ENERGY')[good]
        n0s  = cdata.field('FLUX_DENSITY')[good]
        inds = cdata.field('SPECTRAL_INDEX')[good]
        try:
            betas= cdata.field('beta')[good]
        except KeyErro:
            betas = good.sum()*[np.nan]
        try:
            cuts = cdata.field('Cutoff_Energy')[good]
        except KeyError:
            cuts = good.sum()*[np.nan]
        inds = np.where(inds > 0, inds, -inds)
        try:
            self.names  = np.asarray(cdata.field('Source_Name'))[good]
        except KeyError:
            self.names  = np.asarray(cdata.field('NickName'))[good]

        self.dirs   = map(SkyDir,ras,decs)
        
        def load_model(name, n0,ind,beta, pen, cut):
            if not np.isnan(cut) and cut>0:
                return Models.ExpCutoff(p=[n0,ind,cut],e0=pen)
            if beta<=0 or np.isnan(beta) : 
                return Models.PowerLaw(p=[n0,ind], e0=pen)
            return Models.LogParabola(p=[n0, ind, beta , pen ])
            
        #self.models = np.asarray([load_model(name.strip(),n0,ind,beta,pen,cut)\
        #    for name,n0,ind,beta,pen,cut in zip(self.names,n0s,inds,betas,pens,cuts)])
        self.models = map(load_model, self.names,n0s,inds,betas,pens,cuts )
        if not self.quiet: 
            print 'Loaded %d sources  for roi backgrounds' % (len(cdata),)
            if np.sum(-good)>0:
                print '\tNot using %d entries with TS<25' % np.sum(-good)

  
    def append(self, acat):
        """ 
        append sources found in an auxilliary list of sources
        
        acat: a recarray with fields name, ra, dec, pnorm, pindex
        """
        names  = acat.name
        dirs   = map(SkyDir,acat.ra,acat.dec)
        models = [Models.PowerLaw(p=[n0,ind],e0=1000.) for n0,ind in zip(acat.pnorm,acat.pindex)]
        self.names = np.hstack((self.names, names))
        self.dirs= self.dirs+dirs
        self.models = np.hstack((self.models, models))
        if not self.quiet: print 'Added %d sources to roi background catalog, total now: %d'\
            % (len(acat), len(self.names))
 
        
    def __call__(self,skydir, radius, **kwargs):
        """ return sources, as a list of LocalSource objects, within radius of skydir
            set as free those within free_radius
            exclude those within prune_radius
            
        kwargs:
            free_radius : override the free_radius set by the CTOR
            prune_radius : same
        """
        if len(self.names)==0: return [] # special case: empty catalog
        quiet = kwargs.pop('quiet', self.quiet)
        free_radius = kwargs.pop('free_radius', self.free_radius)
        prune_radius = kwargs.pop('prune_radius', self.prune_radius)
        diffs   = np.degrees(np.asarray([skydir.difference(d) for d in self.dirs]))
        #mask    = ((diffs < radius)&(self.ts > self.min_ts)) & \
        #          ((self.fluxes > self.min_flux)|(diffs < self.max_distance))         
        mask    = (diffs < radius) 
        diffs   = diffs[mask]
        sorting = np.argsort(diffs)

        #sort SkyDirs -- pesky vector behavior...
        dirs    = [x for i,x in enumerate(self.dirs) if mask[i]]
        dirs    = [dirs[x] for x in sorting]
        tslist  = self.ts[mask][sorting]

        names   = self.names[mask][sorting]
        models  = [x.copy() for x in self.models[mask][sorting]]

        # now select those in annulus between prune_radius and free_radius for inclusion in the model
        exclude = (diffs[sorting] < prune_radius) 
        numex = exclude.sum()
        makefree = (diffs[sorting] < free_radius) 
        point_sources = map(LocalSource,dirs,names,models,makefree,tslist)
        if not quiet:
            print '...selected %d sources within %.1f deg for roi'   % (len(point_sources), radius)
            print '...selected %d sources within %.1f deg for refit' % (makefree.sum(), free_radius)
            print 'ts values:', tslist
            if numex>0:
                print '...excluded %d sources within %.1f deg : %s ' %\
            (numex, prune_radius, ', '.join([s.name.strip() for s in point_sources[:numex]]))
        
        self.exclude= point_sources[:numex]
        return point_sources[numex:]

    def merge_lists(self,skydir,radius=15,user_list=None):
        """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence.
         """
        cat_list = self(skydir,radius)
        if user_list is None: return cat_list

        from collections import deque
        merged_list = deque(user_list)

        for ncps,cps in enumerate(cat_list):
            merged_list.append(cps)
            for ups in user_list:
                if np.degrees(ups.skydir.difference(cps.skydir)) < self.prune_radius:
                   merged_list.pop(); break

        merged_list = list(merged_list)
        merged_list.sort(key = lambda ps: ps.skydir.difference(skydir))

        return merged_list
        
    def source_recarray(self):
        """ return a recarry of the full list of sources used for models
        """
        return np.rec.fromarrays(
                [   [name.strip() for name in self.names], 
                    [s.ra() for s in self.dirs], 
                    [s.dec() for s in self.dirs], 
                    [10**m.p[0] for m in self.models],
                    [10**m.p[1] for m in self.models],
                    [(10**m.p[2] if m.name=='ExpCutoff' else np.nan) for m in self.models ],
                    [(10**m.p[3] if m.name=='LogParabola' else  m.e0)  for m in self.models],
                    [(10**m.p[2] if m.name=='LogParabola' else np.nan) for m in self.models],
                ],
                names = 'name ra dec pnorm pindex cutoff pivot beta'.split(),
           )

    def write_reg_file(self, filename, color='green'):
        """ generate a 'reg' file from the catalog, write to outfile
        """
        catrec = self.source_recarray()
        have_ellipse = 'Conf_95_SemiMajor' in catrec.dtype.names
        out = open(filename, 'w')
        print >>out, "# Region file format: DS9 version 4.0 global color=%s" % color
        for s in catrec:
            if have_ellipse:
                print >>out, "fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) #text={%s}" % \
                                (s.ra,s,dec,
                                  s.Conf_95_SemiMinor,Conf_95_SemiMajor,Conf_95_PosAng,
                                  s.name)
            else:
                print >> out, "fk5; point(%.4f, %.4f) # point=cross text={%s}" %\
                                (s.ra, s.dec, s.name)
        out.close()

    def write_xml_file(self, filename, title='source_library'):
        
        point_sources = map(LocalSource,self.dirs,self.names,self.models)
        stacks= xml_parsers.unparse_point_sources(point_sources,strict=True)
        xml_parsers.writeXML(stacks, filename, title=title)

class PipelineCatalogManager(CatalogManager):

    """ subclass of CatalogManager specific to pipeline, for now
        overrides the __init__ to set up from a catalog-like rec file, see description below
    
    """
    def __init__(self, filename, **kwargs):
        """ 
        filename : text
            filename of pickled recarray of sources
            assume in order
            name, ra, dec, pnorm, pindex, cutoff, pnorm_unc, pindex_unc, cutoff_unc, 
            pivot_energy, ts, band_ts, 
            fit_ra, fit_dec, a, b, phi, qual, delta_ts' 
                
        additional kwargs:
            band_ts_min : float
                minimum value for selection of a source
            ts_min : float
                minimum value for TS selection
            update_e0 : bool
                set to force use of pivot_energy as reference
            update_positions: float | None
                if not None, a TS threshold to apply to update source positions from the localization info
            cat_update: None | filename
                file is passed to CatRec.update
        """
        # setup parameters
        band_ts_min = kwargs.pop('band_ts_min', 25)
        update_e0   = kwargs.pop('update_e0',  False)
        update_positions=kwargs.pop('update_positions', None)
        ts_min = kwargs.pop('ts_min', None) # used if updating position
        #cat_update = kwargs.pop('cat_update', None)
        extended_catalog_name =kwargs.pop('extended_catalog', None)
        keyword_options.process(self, kwargs)
        
        # load the file
        crec = CatRec(filename)
        inrec = crec.rec
        print 'loaded %d sources from file %s' % (len(inrec), filename)
        names = 'name ra dec pnorm pindex cutoff pnorm_unc pindex_unc cutoff_unc'
        assert ' '.join(inrec.dtype.names[:len(names.split())])==names,\
            'Structure of recarray not as expected'
        if 'e0' not in inrec.dtype.names: update_e0=True # old form had pivot_energy==e0 only
        
        # make selection(s) and generate lists expected by __call__: dirs, names, models
        cut = (inrec.band_ts>=band_ts_min)+(inrec.band_ts==15) # allow exactly 15
        if ts_min is not None:
            cut =cut * (inrec.ts>ts_min)
            print 'Selected %d after applying ts>%f' % (sum(cut), ts_min)
        self.rec = inrec[cut] #*(inrec.ts>band_ts_min-5)] 
        
        # test for bad fit, reset to starting parameters if so
        badfit=  (inrec.band_ts<100)*(inrec.ts<0)
        beta_index = None if 'beta' not in self.rec.dtype.names else list(self.rec.dtype.names).index('beta')
        #if beta_index is not None:
        #    badfit += inrec.pindex[beta_index]>1.5
        if sum(badfit)>0:
            print 'reseting %d bad fit sources %s' %( sum(badfit), inrec.name[badfit])
            inrec.pnorm[badfit] = 1e-15
            inrec.pindex[badfit] = 2.5
            inrec.beta[badfit]   = 0.01
        if sum(-cut)>0:
            print 'reject %d that fail cut band_ts>%.1f' % (len(inrec)-len(self.rec), band_ts_min)
            print '  names  : %s ' % inrec.name[-cut]
            print '  band_ts: %s ' % inrec.band_ts[-cut]
            print '  ts     : %s ' % inrec.ts[-cut]
        self.dirs = map(SkyDir, self.rec.ra, self.rec.dec)
        self.ts   = self.rec.ts
        if update_positions is not None:
            moved = 0
            for i,s in enumerate(self.rec):
                if s.qual<5 and s.a < 0.2 and s.ts>update_positions and s.delta_ts>1:
                    self.dirs[i] = SkyDir(float(s.fit_ra), float(s.fit_dec))
                    moved +=1
            print 'updated positions of %d sources, with s.qual<5 and s.delta_ts>1 and s.a < 0.2 and s.ts>%.0f '\
                % (moved, update_positions)            
                    
        self.names = self.rec.name
        def make_model(rec):
            pars = list(rec)[3:6]
            beta = rec[beta_index] if beta_index else 0.01
            e0 = rec.e0  #always the reference energy
            if np.isnan(pars[2]):
                #return Models.PowerLaw(p=pars[:2], e0=e0) 
                if beta<=0 or np.isnan(beta): beta=0.001
                p = pars[:2]+[beta, e0]
                return Models.LogParabola(p=np.array(p))
            return Models.ExpCutoff(p=pars, e0=e0)
        self.models = np.array([make_model(r) for r in self.rec])
        print 'selected %d for the catalog' % len(self.rec)
        self.extended_catalog = self.__setup_extended( extended_catalog_name)

    def __setup_extended(self, extended_catalog_name):
        if extended_catalog_name is None: return None
        print 'Loaded extended catalog %s' % extended_catalog_name
        return  pointspec_helpers.ExtendedSourceCatalog(extended_catalog_name)
        
    def get_extended_sources(self,skydir, radius):
        """ add any extended sources with center within the outer radius.
            set parameters free if the center is inside the HEALpix region
        """
        if self.extended_catalog is None: return []
        return self.extended_catalog.get_sources(skydir, radius)
        
    def get_point_sources(self,  skydir, radius, **kwargs):
        return self(skydir, radius, **kwargs)
        
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
    assert len(kwargs.keys())==0, 'unrecognized kw %s' %kwargs 
    filelist.sort()
    
    class CatalogRecArray(object):
        def __init__(self, minflux=1e-16, update_position=False, ts_min=25):
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
        
 


