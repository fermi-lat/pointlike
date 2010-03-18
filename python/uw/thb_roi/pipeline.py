"""
basic pipeline setup

Implement processing of a set of sources in a way that is flexible and easy to use with assigntasks

$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/python/analysis/myroi.py,v 1.2 2009/11/20 22:03:54 burnett Exp $

"""
import sys, os, pyfits, glob, pickle, math
import numpy as np
import pylab as plt

from uw.utilities import makerec, fermitime, image, assigntasks
from uw.like import Models
from skymaps import SkyDir
import myroi, data, catalog, associate # for default configuration (local stuff)

def factory(**kwargs):
    """ define the default factory to use
    """
    return myroi.factory(**kwargs)

class Pipeline(object):
    """ base class for pipeline analysis using assigntasks 
    
    """
    def default_data(self):
        return data.all_data() 
    def default_irf(self):
        return 'P6_v8_diff'

    def __init__(self,  **kwargs):
        self.factory = kwargs.pop('factory', None ) or factory(dataset=self.default_data(), irf=self.default_irf())
        self.__dict__.update({'size':0, 
                    'associate': None, 
                    'bgfree':  np.asarray([True,False,True]),
                    'outdir': '.',
                    'tsfits': True,
                    'roifig': False,
                    'sedfig': True,
                    'free_radius': 1.0, # cut back from 2.0
                })
        self.associate = kwargs.pop('associate', None)
        self.__dict__.update(kwargs)
        outdir = self.outdir
        print 'saving files to "%s"' % outdir
        dirs=[outdir]
        self.pickle_dir = os.path.join(outdir,'pickle'); dirs.append(self.pickle_dir)
        self.tsmap_dir  = os.path.join(outdir,'tsmap');dirs.append(self.tsmap_dir)
        if self.tsfits:
            self.tsfits_dir = os.path.join(outdir, 'tsfits'); dirs.append(self.tsfits_dir)
        else: self.tsfits_dir=None
        if self.roifig:
            self.roifig_dir = os.path.join(outdir, 'roifig'); dirs.append(self.roifig_dir)
        else: self.roifig_dir=None
        if self.sedfig:
            self.sedfig_dir = os.path.join(outdir,'sedfig'); dirs.append(self.sedfig_dir)
        else: self.sedfig_dir=None
        for d in (dirs):
            if not os.path.exists(d): os.mkdir(d)
            
    def setup_string(self):
        """return a string that would generate an instance"""
        pass
        
    def source(self, i):
        raise Exception('base class source invoked: must be overridden in derived class')
    def __getitem__(self,i):
        return self.source(i)

    
    def process(self, s, newsource=True):
        """ process a source s
        """
        name,sdir = s.name, SkyDir(s.ra,s.dec)  
        id_class = 'none' if newsource else s.id_class
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        adict = self.associate(name) if self.associate else None 
        r = self.factory([[name, sdir, Models.PowerLaw()]],  bgfree=self.bgfree)


        # if a different direction, we need to disable the original, most likely the nearest
        if r.psm.point_sources[1].name.strip() == '1FGL %s' % r.name.strip():
            r.psm.models[1].p[0]=-20
            print '-----> disabled catalog source with same name found in background <-------'
        r.dump(maxdist=2)

        ret =r.fit()
        r.dump(maxdist=2, title='after fit')
        print 'TS PL, band = %.1f, %.1f' % (r.TS(), r.band_ts())
        #print r
        #r.printSpectrum()
        move=None
        try:
            r.localize()

            tsf = r.tsmap()
            ts_local_max =tsf(r.tsmax) 
            delta_ts=ts_local_max-tsf(sdir) 
            if adict is not None:
                adict['deltats'] = [ts_local_max-tsf(d) for d in adict['dir']]
                print 'associations:'
                print 'cat name                  ra        dec         ang        Delta TS'
                #       15 Mrk 501               253.4897   39.7527    0.0013      0.41
                for i,id_name in enumerate(adict['name']):
                    print '%3d %-20s%10.4f%10.4f%10.4f%10.2f' %\
                        (adict['cat'][i], id_name, adict['ra'][i], 
                            adict['dec'][i], adict['ang'][i], adict['deltats'][i])
            if r.qform is not None:
                par = r.qform.par 
                psig = np.sqrt(par[3]*par[4])
            else:
                psig = 1.
        except:
            print 'localize failed'
            raise
            delta_ts=99; psig=1; ts_local_max=0
            c2=999.

        if self.size >=0:  

            # make a TS map
            if self.size ==0:
                tsize = 20* psig
            else: tsize = self.size
            if tsize>2: tsize=2 # limit to 2 degrees
            plt.ioff()
            tsm=r.plot_tsmap( outdir=None, catsig=0, size=tsize, 
                # todo: fix this
                assoc=adict if adict is not None else None, # either None or a dictionary
                notitle=True, #don't do title
                )
            # override title, add id class, note if flagged
            tname = name.strip()
            if not newsource:
                if id_class=='   ': id_class = 'unid'
                fname = tname if s.flags==0 else tname+'f'
                tsm.zea.axes.set_title('1FGL %s - %s'% (fname, id_class), fontsize=12) 
            else:
                tsm.zea.axes.set_title(tname, fontsize=12) 


            if self.tsmap_dir is not None: 
                fout =os.path.join(self.tsmap_dir,'%s_tsmap.png'%tname)
                plt.savefig(fout)
                print 'saved tsplot to %s' % fout 

            # now check the TS map for another more significant feature
            move = r.rescan(tsm)
            if move is not None:
                print "------> candidate to move to new location ------------"
        else: tsm=None

        if self.tsfits_dir: 
            if tsm:
                # already have the image: save it as fits
                tsm.zea.skyimage.reimage(tsm.zea.center, os.path.join(self.tsfits_dir, 
                '%s_tsmap.fits'%name), tsize/20., tsize)
            else:
                # a bit of magic to save TS maps in FITS format
                z = image.ZEA(sdir, size=tsize, pixelsize=tsize/20, fitsfile=os.path.join(self.tsfits_dir, '%s_tsmap.fits'%name))
                z.fill(r.tsmap())
                del(z)

        if self.roifig_dir:
            t=r.plot_roi( title='%s ROI'%name, lgclim=[0,2], folder=self.roifig_dir, galactic=True,# nside=2**10, pixelsize=0.01,
            )
            
        if newsource:
            r.pickle(name, self.pickle_dir,
                delta_ts=delta_ts,
                ts=r.TS(),
                band_ts = r.band_ts(),
                tsmap_max = None if move is None else [r.tsmap_max - ts_local_max, move.ra(), move.dec()],
                id = None,
                )
        
        else:
            r.pickle( name, self.pickle_dir, ts=s.signif_avg, 
                cat_fit=[s.pivot_energy, s.spectral_index, s.signif_avg, 0], #avoid old "catsig"
                delta_ts=delta_ts,
                ts=r.TS(),
                band_ts = r.band_ts(),
                tsmap_max = None if move is None else [r.tsmap_max - ts_local_max, move.ra(), move.dec()],
                id = adict,
                )
        return

    def make_sed(self, r, tname):
        if self.sedfig_dir is not None:
            r.plot_sed()
            fout = os.path.join(self.sedfig_dir, ('%s_sed.png'%tname).replace('+','p'))
            plt.title(tname)
            plt.savefig(fout)
            print 'saved SED plot to %s' % fout 

        
    def __call__(self, i):
        self.process(self.source(i))
        
    def run(self):
        """ run it using a multiengine client
        """
        for i in range(self.n):
            self(i)
       

class PipeSourceList(Pipeline):
    def __init__(self, infile, **kwargs):
        super(PipeSourceList,self).__init__(**kwargs)
        self.sourcelist = makerec.textrec(infile)
        self.n = len(self.sourcelist)
        assert(self.n>0)
        
    def source(self,i):
        return self.sourcelist[i]

class RefitCatalog(Pipeline):
    """ manage refitting catalog sources, using the 11-month catalog dataset default
    """
    def default_data(self):
        return data.Catalog_noGRB();
    def default_irf(self):
        return 'P6_v3_diff'
        
    def __init__(self, factory=None, subset=None,  **kwargs):
        """ subset: bool array to select sources to refit
        
        """
        factory = factory or myroi.factory(self.default_data(), irf=self.default_irf(), gti_mask=data.gti_noGRB())
        super(RefitCatalog,self).__init__(factory=factory, 
                associate=associate.Association(),
                **kwargs)
        cat =pyfits.open(os.path.join(data.catalog_path,catalog.default_catalog))
        cdata = cat[1].data
        csnames=cdata.Source_Name
#        cjnames=np.asarray([n[5:]  if n[-1]!='c' else n[5:-1]  for n in csnames])
        cjnames=np.asarray([n[5:]  for n in csnames])
        cata = cdata.Conf_95_SemiMajor
        catb = cdata.Conf_95_SemiMinor
        catang=cdata.Conf_95_PosAng
        cindex =cdata.Spectral_Index
        cpivot =cdata.Pivot_Energy
        csignif=cdata.Sqrt_TS300_1000 #place holdef?
        cra    = np.asarray(cdata.RA,float)
        cdec   = np.asarray(cdata.DEC,float)
        cid_class   = cdata.CLASS1
        subset = subset or np.index_exp[:len(cdata)]

        cflags = np.asarray(cdata.Flags,int)
        class Source():
            def __init__(self, *par):
                p = par[0]
                self.name,self.ra,self.dec = p[0],p[1],p[2]
                self.id_class, self.flags = p[3],p[4]
                self.pivot_energy,self.spectral_index, self.signif_avg, self.a95,self.b95,self.ang = p[5:11]

        self.sourcelist = [Source(par) for par  in\
            zip(cjnames[subset], cra[subset], cdec[subset], 
                cid_class[subset],cflags[subset],
                cpivot[subset],cindex[subset],csignif[subset],
                cata[subset],catb[subset],catang[subset])]
        self.n = len(self.sourcelist)
        assert(self.n>0)
        
    def source(self,i):
        return self.sourcelist[i]
        
    def __call__(self, i):
        self.process(self.source(i)) 
        
    def process(self, s):
        name,sdir = s.name, SkyDir(s.ra,s.dec)  
        id_class = s.id_class
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        try:
            adict = self.associate(name) if self.associate else None 
        except:
            adict = None 
        r = self.factory([[name, sdir, Models.PowerLaw()]],  bgfree=self.bgfree)

        # if a different direction, we need to disable the original, most likely the nearest
        if r.psm.point_sources[1].name.strip() == '1FGL %s' % r.name.strip():
            r.psm.models[1].p[0]=-20
            print '-----> disabled catalog source with same name found in background <-------'
        r.dump(maxdist=2)
        ret =r.fit()
        r.dump(maxdist=2, title='after fit')
        r.localize()
        r.find_tsmax()
        tsf = r.tsmap()
        ts_local_max =tsf(r.tsmax) 
        delta_ts=ts_local_max-tsf(sdir) 

        if np.isnan(s.a95):
            tsize = 20*r.qform.par[3]
        else:
            tsize = min(20* s.a95/2.45, 2)
        tsm=r.plot_tsmap( outdir=None, catsig=0, size=tsize, 
            # todo: fix this
            assoc=adict if adict is not None else None, # either None or a dictionary
            notitle=True, #don't do title
            markersize=10,
            primary_markersize=12,
            )
        if not np.isnan(s.a95) and s.a95<1.0:
            print 'overplotting catalog fit: %.3f %.3f %.1f' % (s.a95,s.b95,s.ang)
            tsm.overplot((s.ra,s.dec,s.a95,s.b95,s.ang), contours=[1.0], lw=2, ls='-', color='gray')
        else:
            print 'Source has no catalog error circle'
        # override title, add id class,
        tname = name.strip()
        if id_class=='   ': id_class = 'unid'
        fname = tname #if s.flags==0 else tname+'f'
        tsm.zea.axes.set_title('1FGL %s - %s'% (fname, id_class), fontsize=12) 

        if self.tsmap_dir is not None: 
            fout =os.path.join(self.tsmap_dir, ('%s_tsmap.png'%tname).replace('+','p'))
            plt.savefig(fout)
            print 'saved tsplot to %s' % fout 
        self.make_sed(r, tname)
        r.pickle( name, self.pickle_dir, ts=s.signif_avg, 
                cat_fit=[s.pivot_energy, s.spectral_index, s.signif_avg, 0], #avoid old "catsig"
                delta_ts=delta_ts,
                ts=r.TS(),
                band_ts = r.band_ts(),
                tsmap_max = None ,
                id = adict,
                )
        del r # memory problem???
        return

class FitNewCatalog(RefitCatalog):
    def default_data(self):
        return data.all_data() #needs mask!
    def default_irf(self):
        return 'P6_v8_diff'
    def __init__(self, newcats, **kwargs):
        """ newcats: list of catalog text filenames
        """
        super(FitNewCatalog,self).__init__(**kwargs)
        for cat in newcats:
            self.factory.cb.cm.append(makerec.textrec(cat)) 
        self.sources = self.factory.source_list() #newcat
        self.fit_method = 'simplex' #'minuit'
        self.n = len(self.sources)
    def source(self, i):
        class Source:
            def __init__(self, i, name, ra,dec):
                self.name = name#'UW3\2%05d' % (i+1)
                self.ra,self.dec=ra,dec
            def __str__(self): return '%s  %.3f %+0.3f' % (self.name, self.ra, self.dec)
        return Source(i, self.sources.name[i], self.sources.ra[i], self.sources.dec[i])
    def process(self, s):
        name,sdir = s.name, SkyDir(s.ra,s.dec)  
        adict = None
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        r = self.factory([[name, sdir, Models.PowerLaw()]],  bgfree=self.bgfree)

        # if a different direction, we need to disable the original, most likely the nearest
        if r.psm.point_sources[1].name.strip() == '%s' % r.name.strip():
            r.psm.models[1].p[0]=-20
            print '-----> disabled catalog source with same name found in background <-------'
        r.dump(maxdist=2)
        ret =r.fit(method=self.fit_method)
        r.dump(maxdist=2, title='after fit')
        r.localize()
        r.find_tsmax()
        tsf = r.tsmap()
        ts_local_max =tsf(r.tsmax) 
        delta_ts=ts_local_max-tsf(sdir) 

        tsize = 0.5 if r.qform is None else 20*r.qform.par[3]
        tname = name
        tsm=r.plot_tsmap( outdir=None, catsig=0, size=tsize, 
            # todo: fix this
            assoc=adict if adict is not None else None, # either None or a dictionary
            notitle=True, #don't do title
            markersize=10,
            primary_markersize=12,
            )
        tsm.zea.axes.set_title('%s'% tname, fontsize=12) 

        if self.tsmap_dir is not None: 
            fout =os.path.join(self.tsmap_dir, ('%s_tsmap.png'%tname).replace('+','p'))
            plt.savefig(fout)
            print 'saved tsplot to %s' % fout 
        if self.tsfits_dir: 
            tsm.zea.skyimage.reimage(tsm.zea.center, os.path.join(self.tsfits_dir, 
                '%s_tsmap.fits'%name), tsize/20., tsize)

        r.pickle( name, self.pickle_dir, 
                delta_ts=delta_ts,
                ts=r.TS(),
                band_ts = r.band_ts(),
                tsmap_max = None ,
                id = None,
                )
        if self.sedfig_dir is not None:
            r.plot_sed()
            fout = os.path.join(self.sedfig_dir, ('%s_sed.png'%tname).replace('+','p'))
            plt.title(tname)
            plt.savefig(fout)
            print 'saved SED plot to %s' % fout 
 
        return

        

class UWsourceFits(Pipeline):
    def __init__(self, sourcelist, **kwargs):
        self.sources=makerec.textrec(sourcelist)
        self.sources.sort()
        super(UWsourceFits,self).__init__(**kwargs)
       
    def source(self, i):
        class Source:
            def __init__(self, i, name, ra,dec):
                self.name = name #'1UW%05d' % (i+1)
                self.ra,self.dec=ra,dec
            def __str__(self): return '%s  %.3f %+0.3f' % (self.name, self.ra, self.dec)
        return Source(i, self.sources.name[i], self.sources.ra[i], self.sources.dec[i])
 
    def process(self, s):
        name,sdir = s.name, SkyDir(s.ra,s.dec)  
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        r = self.factory([[name, sdir, Models.PowerLaw()]],  bgfree=self.bgfree)

        # if a different direction, we need to disable the original, most likely the nearest
        if r.psm.point_sources[1].name.strip() == '%s' % r.name.strip():
            r.psm.models[1].p[0]=-20
            print '-----> disabled catalog source with same name found in background <-------'
        r.dump(maxdist=2)
        ret =r.fit()
        r.dump(maxdist=2, title='after fit')
        r.localize()
        r.find_tsmax()
        tsf = r.tsmap()
        ts_local_max =tsf(r.tsmax) 
        delta_ts=ts_local_max-tsf(sdir) 

        tsize = 0.5 # default
        if r.qform is not None: tsize = 20*r.qform.par[3]
        tname=name
        tsm=r.plot_tsmap( outdir=None, catsig=0, size=tsize, 
            # todo: fix this
            assoc=None, #adict if adict is not None else None, # either None or a dictionary
            notitle=True, #don't do title
            markersize=10,
            primary_markersize=12,
            )
        tsm.zea.axes.set_title('%s'% name, fontsize=12) 

        if self.tsmap_dir is not None: 
            fout =os.path.join(self.tsmap_dir, ('%s_tsmap.png'%tname).replace('+','p'))
            plt.savefig(fout)
            print 'saved tsplot to %s' % fout 
        if self.sedfig_dir is not None:
            r.plot_sed()
            fout = os.path.join(self.sedfig_dir, ('%s_sed.png'%tname).replace('+','p'))
            plt.title(tname)
            plt.savefig(fout)
            print 'saved SED plot to %s' % fout 

            
        r.pickle( name, self.pickle_dir, #ts=s.signif_avg, 
                delta_ts=delta_ts,
                ts=r.TS(),
                band_ts = r.band_ts(),
                tsmap_max = None ,
                id = None,
                )
        return

class TrialSourceFits(Pipeline):
    def __init__(self, sources, auxcat=None, **kwargs):
        """
        sourcelist: recarray with list of names, ra dec
        
        
        """
        self.sources=sources
        super(TrialSourceFits,self).__init__(**kwargs)
        if auxcat is not None: self.factory.cb.cm.append(auxcat) 
        self.fit_method = 'minuit'
        self.__dict__.update(kwargs)
        self.n = len(self.sources)

def setup_stuff(classname, **kwargs):
    """  generate strings appropriate for assigntasks
    """
    args = ','.join(['%s=%s' %(key, kwargs[key]) for key in kwargs.keys()])
    setup_string = 'from thb_roi import pipeline; reload(pipeline); g = pipeline.%s(%s)' % (classname, args)
    exec(setup_string)
    n = g.n
    tasklist = ['g(%d)'%i for i in range(n)]
    
    return setup_string, tasklist

def load_rec_from_pickles(outdir):
    """
    create recarray from list of pickled 
    outdir -- folder for analysis results, expect to find "pickle" below it
    """
    failed=0
    assert(os.path.exists(os.path.join(outdir,'pickle')))
    filelist = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
    assert(len(filelist)>0)
    filelist.sort()
    rec = makerec.RecArray("""name ra dec csig psig good dcp qual cindex pnorm pindex
                              delta_ts fit_ra fit_dec a b ang ts band_ts galnorm isonorm
                              tsmap_max_diff tsm_ra tsm_dec
                              id_prob id_dts id_cat""".split())
    for fname in filelist:
        #print 'loading %s...' % fname,
        try:
            p = pickle.load(open(fname))
            name= p['name']
            ra,dec =p['ra'],p['dec']
            qfp = p['qform_par']
            fit_ra,fit_dec=ra,dec  # leave same if bad fit
            good = True
            if qfp:
                ptsig = math.sqrt(qfp[3]*qfp[4])
                dcp = math.degrees(SkyDir(ra,dec).difference(SkyDir(qfp[0],qfp[1])))
                qual = qfp[6]
                a,b,ang = qfp[3:6]
                fit_ra,fit_dec=qfp[0:2]
            else: 
                ptsig = 1
                dcp = 1
                qual=99
                a,b,ang = 0,0,0
                fit_ra,fit_dec,qual=999,999,99
                good = False
            if 'cat_fit' in p:
                csig = p['cat_fit'][3] 
                cindex = p['cat_fit'][1]
            else: csig=cindex=-1
            pnorm,pindex  = p['src_par']
            #if pindex>5 or pindex<1: good=False
            if np.isinf(csig): good=False
            if ptsig>1: 
                ptsig=1
                good = False
            delta_ts=p['delta_ts']
            ts = p['ts']
            band_ts= p['band_ts']
            tsmap_max_diff, tsm_ra, tsm_dec = [0,0,0] if 'tsmap_max' not in p or p['tsmap_max'] is None  else p['tsmap_max']
            galnorm = p['bgm_par'][0]
            isonorm = p['bgm_par'][2]

            id_prob =p.get('id_prob', 0)
            id_dts = p.get('id_dts', 99)
            id_cat = p.get('id_cat', 99) 
            
            rec.append(name, ra,dec, csig, ptsig, \
                good , dcp,  qual, cindex, pnorm, pindex,\
                delta_ts, fit_ra, fit_dec, a, b, ang, ts, band_ts, galnorm, isonorm, tsmap_max_diff, tsm_ra, tsm_dec, 
                id_prob, id_dts, id_cat)
        except Exception:
            #raise
            failed +=1
    print 'read %d entries from %s (%d failed)' % (len(filelist),outdir,failed)
    return rec()

class Logs(object):
    """ access the logs from a run"""
    def __init__(self, fname = 'summary.pickle'):
        self.logs = pickle.load(open(fname))['result']
      
    def __call__(self,name):
        check = cat.source_name=='1FGL %s' % name
        if check.sum()!=1:
            print 'name %s not found in catalog' %name
            return None
        index = np.argmax(check)
        return self.logs[index] 
        
    def saveto(self, folder='logs'):
        if not os.path.exists(folder): os.mkdir(folder)
        for key in self.logs.keys():
            name = cat.source_name[key].split()[1]
            out = open('%s/%s.log' % (folder,name), 'w')
            out.writelines(self.logs[key].__str__())
            out.close()
            


if __name__=='__main__':
    pass
#    f = factory(emin=1000)

