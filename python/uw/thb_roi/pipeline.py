"""
basic pipeline setup

Implement processing of a set of sources in a way that is flexible and easy to use with assigntasks

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/pipeline.py,v 1.17 2010/08/24 22:02:19 burnett Exp $

"""
version='$Revision: 1.17 $'.split()[1]
import sys, os, pyfits, glob, pickle, math, time, types
import numpy as np
import pylab as plt
from scipy import optimize

from uw.utilities import makerec, fermitime, image
from uw.utilities.assigntasks import setup_mec, AssignTasks, get_mec, kill_mec
import uw
from skymaps import SkyDir
from uw.thb_roi import myroi, associate , config, roi_factory, variation

class InvalidArgument(Exception):
    pass
  
 
class Pipeline(object):
    """ base class for pipeline analysis using assigntasks 
    
    """
    def __init__(self,  **kwargs):

        pipeline_opts = dict(
                associate = None, 
                outdir    = '.',
                tsmapfig  = True,
                tsfits    = False,
                roifig    = False,
                sedfig    = True,
                bgfree    =  [True,False,True],
                fit_method= 'simplex',
                seeds     = None,
                seeds_only= False,
                sed_limits= None, 
                size      = 0,
                pivot     = True,
                butterfly = True,
                write_xml = False, # set True to write out model XML description
            )
        for key in pipeline_opts.keys():
            if key in kwargs:
                pipeline_opts[key] = kwargs.pop(key)
        
        outdir = pipeline_opts['outdir']
        if not os.path.exists(outdir): os.mkdir(outdir)
        kwargs['log'] = open(os.path.join(outdir, 'configuration.txt'), 'w')        
                
        self.factory = roi_factory.ROIfactory(**kwargs)
        self.__dict__.update(pipeline_opts)
        
        if type(self.associate)==type(''):
            self.associate = associate.SrcId(self.associate)
            print 'will associate with catalogs %s' % self.associate.classes
        outdir = self.outdir
        print 'saving files to "%s"' % outdir
        dirs=[outdir]
        self.pickle_dir = os.path.join(outdir,'pickle'); dirs.append(self.pickle_dir)
        if self.tsmapfig:
            self.tsmap_dir  = os.path.join(outdir,'tsmap');dirs.append(self.tsmap_dir)
        else: self.tsmap_dir = None
        if self.tsfits:
            self.tsfits_dir = os.path.join(outdir, 'tsfits'); dirs.append(self.tsfits_dir)
        else: self.tsfits_dir=None
        if self.roifig:
            self.roifig_dir = os.path.join(outdir, 'roifig'); dirs.append(self.roifig_dir)
        else: self.roifig_dir=None
        if self.sedfig:
            self.sedfig_dir = os.path.join(outdir,'sedfig'); dirs.append(self.sedfig_dir)
        else: self.sedfig_dir=None
        if self.write_xml:
            self.xml_dir = os.path.join(outdir,'xml'); dirs.append(self.xml_dir)
        else: self.xml_dir=None
        for d in (dirs):
            if not os.path.exists(d): os.mkdir(d)
                
    def source(self, i):
        raise Exception('base class source invoked: must be overridden in derived class')
    def __len__(self):
        return self.n
    def __getitem__(self,i):
        if type(i)==types.StringType:
            i = self.names.index(i)  # will raise value error
        if i==self.n: raise StopIteration
        if i>self.n: raise  IndexError('list index out of range')
        return self.source(i)

        
    def process(self, s):
        name,sdir = str(s.name).strip(), SkyDir(s.ra,s.dec) 
        if name[0]=='1':
            name=name.replace('_',' ').replace('p','+') # convert back if using file version
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        r = self.factory([(name, sdir)], bgfree=self.bgfree) #use default model parameters

        # if a different direction, we need to disable the original, most likely the nearest
        if len(r.psm.point_sources)>1 and r.psm.point_sources[1].name.strip() == '%s' % r.name.strip():
            r.psm.models[1].p[0]=-20
            print '-----> disabled catalog source with same name found in background <-------'
        r.dump(maxdist=2)
        ts = r.fit(method=self.fit_method, pivot=self.pivot)
        r.dump(maxdist=2, title='after fit')
        binfo = r.band_info()
        print '\n     TS: %.1f, %.1f' % ( ts, r.TS())
        print   'band TS: %.1f = ' %r.band_ts(), '+'.join(['%.1f'%ts for ts in binfo['ts']])
        
        print '\n========  LOCALIZATION   =============='
        tsmaxpos, delta_ts = r.localize()
        localized = r.qform is not None
        if localized:
            par = r.qform.par 
            psig = np.sqrt(par[3]*par[4])
            ellipse = par[3:6]
        else:
            # failed: try again with new position
            print 'trying local maximum'
            localmax = image.find_local_maximum(r.tsmap(), sdir)
            tsmaxpos, delta_ts = r.localize(seedpos=localmax)
            localized = r.qform is not None
            if localized:
                print 'succeeded'
                par = r.qform.par 
                psig = np.sqrt(par[3]*par[4])
                ellipse = par[3:6]
            
            else:
                print 'failed'
                psig = 0.5; ellipse = (psig,psig,0) #default 
                tsmaxpos = r.center # just for deltata below
        if psig>0.5 or not localized: #protection
            oldsig = psig
            psig = 0.5; ellipse = (psig,psig,0)
            print 'warning: too large error ellipse, %.1f deg, cut to %.1f' % (oldsig, psig)
            

        print '\n========   ASSOCIATION  =============='
        adict = self.associate(name, sdir, ellipse) if self.associate else None
        if adict is not None:
        
            tsf = r.tsmap()
            ts_local_max=tsf(tsmaxpos)
            adict['deltats'] = [ts_local_max-tsf(d) for d in adict['dir']]
            if localized: print 'associations:'
            else: print 'associations (assuming 0.5 deg error circle at seed position)'
            print '   cat         name                  ra        dec         ang     prob    Delta TS'
            #       15 Mrk 501               253.4897   39.7527    0.0013      0.41
            fmt = '   %-10s %-20s%10.4f%10.4f%10.4f%8.2f%8.1f' 
            for i,id_name in enumerate(adict['name']):
                tup = (adict['cat'][i], id_name, adict['ra'][i], adict['dec'][i], adict['ang'][i], 
                        adict['prob'][i],adict['deltats'][i])
                print fmt % tup
        else:
            print 'No associations found'


        tsize = 1.5 if not localized else 20*psig
        tname = name
        fname = tname.strip().replace('+','p').replace(' ', '_') # filename
        
        if self.tsmap_dir is not None:
            tsm=r.plot_tsmap( outdir=None, catsig=0, size=tsize, 
                # todo: fix this
                assoc=adict if adict is not None else None, # either None or a dictionary
                notitle=True, #don't do title
                markersize=10,
                primary_markersize=12,
                )
            tsm.zea.axes.set_title('%s'% tname, fontsize=18)  # big title
            self.tsmap_append(tsm, s)

            fout =os.path.join(self.tsmap_dir, ('%s_tsmap.png'%fname) )
            plt.savefig(fout)
            print 'saved tsplot to %s' % fout 
            if self.tsfits_dir: 
                tsm.zea.skyimage.reimage(tsm.zea.center, os.path.join(self.tsfits_dir, 
                    '%s_tsmap.fits'%fname), tsize/20., tsize)

        self.pickle( r, name , s,
                fname=fname,
                delta_ts=delta_ts,
                ts=ts, ts2=r.TS(), #inital, current TS
                band_ts = r.band_ts(),
                tsmap_max = None ,
                id = None,
                adict=adict,
                band_info = r.band_info(), 
                parameters = r.get_parameters(),
                )
        
        if self.sedfig_dir is not None:
            self.make_sed(r, tname, fname, limits=self.sed_limits, butterfly = self.butterfly)
  
        if self.write_xml:
            fout = os.path.join(self.xml_dir, fname+'.xml')
            r.toXML( fout)
            print 'saved xml to %s' % fout 
        return

    def pickle(self, roi, name, s, **kwargs):
        """ intercept this to allow subclass to add or subtract stuff """
        roi.pickle( name, self.pickle_dir, **kwargs)
    
    def tsmap_append(self, tsm, s):
        pass

    def make_sed(self, r, tname, fname, limits, **kwargs):
        if self.sedfig_dir is not None:
            r.plot_sed(axis=limits, **kwargs)
            fout = os.path.join(self.sedfig_dir, ('%s_sed.png'%fname) )
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


class LightCurve(Pipeline): 
    """ part of a light curve
    """
    def __init__(self, outdir, monthly, **kwargs):
        """
        monthly: a dict with keys 'number' and 'dataset'
        """
        self.month = monthly['number']
        kwargs['dataset'] = monthly['dataset']
        self.sources = glob.glob(os.path.join(outdir, 'pickle', '*.pickle'))
        self.namelist = [pickle.load(open(source))['name'] for source in self.sources] 
        self.n = len(self.sources)
        assert(self.n>0)
        super(LightCurve,self).__init__(**kwargs)
        
    def source(self,i):
        if type(i)==types.StringType:
            i = self.namelist.index(i.replace('_',' ').replace('p','+'))
        return self.sources[i]
    def names(self):
        return self.namelist


    def roi(self, i): 
        """ return roi for ith source (or source name i)"""
        s = self.source(i)
        pk = pickle.load(open(s))
        parameters = pk['parameters']
        name = pk['name']
        ra, dec = pk['ra'], pk['dec']
        sdir = SkyDir(ra,dec)
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        if name[:4]=='1FGL': name = name[5:]
        roi = self.factory('%-s %.4f %.4f' % (name, ra, dec))
        roi.set_parameters(parameters)
        fr = roi.get_free()
        fr[1:]=False
        roi.set_free(fr)
        model = roi.psm.models[0]
        model.set_e0(pk['pivot_energy'])
        return roi

    def process(self, s):
        """ process the source described by the pickle file s
            save all band info, 
            likelihood function evaluated at 8 points: 
                max, -sigma, +sigma, 95%, 0, avflux, avflux-sigflux, avflux+sigflux
        """
        pk = pickle.load(open(s))
        parameters = pk['parameters']
        name = pk['name']
        ra, dec = pk['ra'], pk['dec']
        sdir = SkyDir(ra,dec)
        print '\n'+60*'=','\nprocessing source %s at %s' % (name, sdir), '\n'+ 60*'='
        if name[:4]=='1FGL': name = name[5:]
        
        roi = self.factory('%-s %.4f %.4f' % (name.replace(' ','_'), ra, dec))
        roi.set_parameters(parameters)
        model = roi.psm.models[0]
        model.set_e0(pk['pivot_energy'])
        
        # freeze all but first parameter, then refit
        fr = roi.get_free()
        fr[1:]=False
        roi.set_free(fr)
        roi.fit(pivot=False, ignore_exception=True, use_gradient=False)
        roi.dump(maxdist=2, title='after refit to flux')
        
        avflux =pk['src_par'][0]
        sigflux=pk['src_par_unc'][0]
        if avflux<1e-14: avflux, siglfux = 1e-14,0 #avoid original failure, slog on
        w = variation.LogLikelihood(roi) #, maxflux=10*avflux
        print 'Likelihood: \n', w
        t = (w(avflux), w(avflux-sigflux),w(avflux+sigflux))
        print (' likelihood values at avflux=%.2e, +/-sigma=%.2e: '+3*'%.2e, ') %\
            (avflux, sigflux, t[0],t[1],t[2])

        if 'months' not in pk: pk['months'] = dict()
        pk['months'][self.month] = dict( 
            band_info =roi.band_info(),
            likelihood = (w.maxl,) +  w.errors() + (w.upper_limit(),  w.TS())+ t
            )
        print 'updating file %s' % s
        pickle.dump(pk, open(s, 'wb'))
        
        #w.plot(fignum = self.month+100) # see what is going on
        

class PipeSourceList(Pipeline):
    def __init__(self, infile, **kwargs):
        super(PipeSourceList,self).__init__(**kwargs)
        self.sourcelist = makerec.textrec(infile)
        self.n = len(self.sourcelist)
        assert(self.n>0)
        
    def source(self,i):
        return self.sourcelist[i]

#class RefitCatalog(Pipeline):
#    """ manage refitting catalog sources, using the 11-month catalog dataset default
#    """
#    def default_data(self):
#        return data.Catalog_noGRB();
#    def default_irf(self):
#        return 'P6_v3_diff'
#        
#    def __init__(self, factory=None, subset=None,  **kwargs):
#        """ subset: bool array to select sources to refit
#        
#        """
#        factory = factory or uw.factory(dataset=self.default_data(), irf=self.default_irf(), gti_mask=data.gti_noGRB())
#        self.fit_method = 'simplex' #'minuit'
#        super(RefitCatalog,self).__init__(factory=factory, 
#                associate=associate.Association(),
#                **kwargs)
#        cat =pyfits.open(os.path.join(data.catalog_path,catalog.default_catalog))
#        cdata = cat[1].data
#        csnames=cdata.Source_Name
##        cjnames=np.asarray([n[5:]  if n[-1]!='c' else n[5:-1]  for n in csnames])
#        cjnames=np.asarray([n[5:]  for n in csnames])
#        cata = cdata.Conf_95_SemiMajor
#        catb = cdata.Conf_95_SemiMinor
#        catang=cdata.Conf_95_PosAng
#        cindex =cdata.Spectral_Index
#        cpivot =cdata.Pivot_Energy
#        csignif=cdata.Sqrt_TS300_1000 #place holdef?
#        cra    = np.asarray(cdata.RA,float)
#        cdec   = np.asarray(cdata.DEC,float)
#        cid_class   = cdata.CLASS1
#        subset = subset or np.index_exp[:len(cdata)]
#
#        cflags = np.asarray(cdata.Flags,int)
#        class Source():
#            def __init__(self, *par):
#                p = par[0]
#                self.name,self.ra,self.dec = p[0],p[1],p[2]
#                self.id_class, self.flags = p[3],p[4]
#                self.pivot_energy,self.spectral_index, self.signif_avg, self.a95,self.b95,self.ang = p[5:11]
#
#        self.sourcelist = [Source(par) for par  in\
#            zip(cjnames[subset], cra[subset], cdec[subset], 
#                cid_class[subset],cflags[subset],
#                cpivot[subset],cindex[subset],csignif[subset],
#                cata[subset],catb[subset],catang[subset])]
#        self.sourcenames = cjnames[subset]
#        self.n = len(self.sourcelist)
#        assert(self.n>0)
#        
#    def source(self, i):
#        if type(i)==types.StringType:
#            try:
#                i = list(self.sourcenames).index(i)
#            except:
#                raise InvalidArgument('source %s not found in list of sources' % i)
#        return self.sourcelist[i]
#    def names(self):
#        return self.sourcenames
#        
#    def __call__(self, i):
#        self.process(self.source(i)) 
#
#    def tsmap_append(self, tsm, s):
#        if not np.isnan(s.a95) and s.a95<1.0:
#            print 'overplotting catalog fit: %.3f %.3f %.1f' % (s.a95,s.b95,s.ang)
#            tsm.overplot((s.ra,s.dec,s.a95,s.b95,s.ang), contours=[1.0], lw=2, ls='-', color='gray')
#        else:
#            print 'Source has no catalog error circle'
#
#    def pickle(self, roi, name, s, **kwargs):
#        """ intercept this to add fit info from catalog """
#        
#        roi.pickle( name, self.pickle_dir,
#            cat_fit=[s.pivot_energy, s.spectral_index, s.signif_avg, 0], **kwargs)
#
#
class FitNewCatalog(Pipeline):
    """
    pipeline with extended catalog
    """
    def default_data(self):
        return data.all_data() #needs mask!
    def default_irf(self):
        return 'P6_v8_diff'
    def __init__(self, newcats, **kwargs):
        """ newcats: list of catalog text filenames
        """
        super(FitNewCatalog,self).__init__(**kwargs)
        for cat in newcats:
            self.factory.catman.append(makerec.textrec(cat)) 
        self.sources = self.factory.source_list() #newcat
        self.ncat= len(self.sources)
        #self.fit_method = 'simplex' #'minuit'
        seedfile = kwargs.pop('seeds', None)
        self.n = len(self.sources) 
        if seedfile is not None:
            self.seeds = makerec.load(seedfile)
            self.n += len(self.seeds)
            print 'Loaded %d seeds from file %s' % (len(self.seeds), seedfile)
        
    def source(self, i):
        """ i is the index, or the name of a source"""
        class Source:
            def __init__(self, i, name, ra,dec):
                self.name = str(name)#'UW3\2%05d' % (i+1)
                self.ra,self.dec=ra,dec
            def __str__(self): return '%s  %.3f %+0.3f' % (self.name, self.ra, self.dec)
        if type(i)==types.StringType:
            try:
                i = self.names().index(i)
            except ValueError:
                raise InvalidArgument('Source "%s" not found in list of names' %i)
        if i< len(self.sources):
            return Source(i, self.sources.name[i], self.sources.ra[i], self.sources.dec[i])
        else:
            return self.seeds[i-len(self.sources)]
    def names(self):
        return [self.source(i).name for i in range(self.n)]

class UWsourceFits(Pipeline):
    def __init__(self, sourcelist, **kwargs):
        if type(sourcelist)==types.StringType:
            # assume filename
            ft = os.path.splitext(sourcelist)[-1]
            if ft=='.txt':
                t=makerec.textrec(sourcelist)
                if 'name' not in t.dtype.names:
                    self.sources = plt.mlab.rec_append_fields(t, 'name', t.source_name)
                else: self.sources=t    
            elif ft=='.fit' or ft=='.fits':
                t=makerec.fitsrec(sourcelist)
                self.sources = plt.mlab.rec_append_fields(t, 'name', t.source_name)
            else:
                raise InvalidArgument('Source filetype "%s" not recognized' % ft)
        else: self.sources=sourcelist    
        self.sources.sort()
        self.n = len(self.sources)
        super(UWsourceFits,self).__init__(**kwargs)

    def names(self):
        return self.sources.name
        
    def source(self, i):
        if type(i)==types.StringType:
            try:
                i = list(self.sources.name).index(i)
            except:
                raise InvalidArgument('source %s not found in list of sources' % i)
        class Source:
            def __init__(self, i, name, ra,dec):
                self.name = name #'1UW%05d' % (i+1)
                self.ra,self.dec=ra,dec
            def __str__(self): return '%s  %.3f %+0.3f' % (self.name, self.ra, self.dec)
        return Source(i, self.sources.name[i], self.sources.ra[i], self.sources.dec[i])
 

 
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

def getg(setup_string):
    exec(setup_string, globals()) # should set g
    return g
            

def main( setup_string, outdir, mec=None, startat=0, n=0, local=False,
        machines='tev1 tev2 tev3 tev4'.split(), engines=None,
        seeds_only=False,
        ignore_exception=True, 
        logpath='log'):
        
    if not os.path.exists(outdir): 
        os.mkdir(outdir)
        #raise Exception('outdir folder, "%s", not found' % outdir)
    if setup_string[0]=='@':
        setup_string = open(setup_string).read()
    g = getg(setup_string)
    names = g.names(); 
    if seeds_only:
        startat = g.ncat+n
        print 'Processing seeds only, will start at source %d' % startat
    else:
        print 'Start at source %d' % startat
    if len(names)==0:
        raise InvalidArgument('no tasks defined by Pipeline object')
    else:
        print 'found %d sources to process' % len(names)

    if n==0: endat = len(names)
    else: endat = min(startat+n, len(names))
         
    tasks = ['g(%d)'% i for i in range(len(names))]
    del g #do not need this instance
    
    def callback(id, result):
        try:
            name = names[id+startat]
            logdir = os.path.join(outdir, logpath)
            if not os.path.exists(logdir): os.mkdir(logdir)
            out = open(os.path.join(logdir, '%s.txt' % name.strip().replace(' ','_').replace('+','p')), 'w')
            print >>out, result
            out.close()
        except:
            print 'got exception writing %s' % name
            raise
            
    if not local: setup_mec(machines=machines, engines=engines)
    time.sleep(10)
    lc= AssignTasks(setup_string, tasks[startat:endat], mec=mec, timelimit=1000, local=local, callback=callback, 
    ignore_exception=ignore_exception)

    lc(5)
    if not local: get_mec().kill(True)
    return lc
   
if __name__=='__main__':
    from optparse import OptionParser
    usage = """\n
usage: %prog [options] setup_string outdir \n
Run the UW pipeline
    setup_string: python executable string that sets a variable g; if start with "@", name of file
    outdir: where save the files, under folders log, tsmap, sed"""
    parser = OptionParser(usage, version=version)
    parser.add_option('-l', '--local', help='run locally', action='store_true', dest='local',default=False)
    parser.add_option('-x', '--stop_on_exception', help='do not ignore exceptions', 
                    action='store_false', dest='ignore_exception',default=True)
    parser.add_option('s', '--start_at', help='initial task to run (default %default)', dest='startat', default=0, type='int')
    parser.add_option('n', '--number', help='number of tasks to run', dest='n', default=0, type = 'int')
    options, args = parser.parse_args()
    if len(args)!=2: 
        parser.print_usage()
        sys.exit(-1)
    main(args[0],args[1], local=options.local, ignore_exception=options.ignore_exception, startat=option.startat, n=option.n)