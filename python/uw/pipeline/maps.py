"""
Code to generate a set of maps for each ROI
$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/pipeline/ts_map.py,v 1.6 2011/01/01 15:50:05 burnett Exp $

"""
import os, glob, pickle, types
import numpy as np
from uw.like import roi_tsmap, roi_plotting # for TSCalc
from uw.utilities import image
from uw.pipeline import skymodel
from skymaps import Band, SkyDir, PySkyFunction, Hep3Vector 
import pylab as plt

#from pipeline import roi_maps
def table_array(name, outdir):
    files =glob.glob(os.path.join(outdir, '%s_table'%name,'*.pickle'))
    nf = len(files)
    assert nf>0, 'no pickle files found in %s' % os.path.join(outdir, '%s_table'%name)
    if nf<1728: print 'warning: missing %d files' % (1728-nf)
    files.sort()

    pklist = np.array([pickle.load(open(f)) for f in files])
    return pklist
def hpname(index):
    return 'HP12_%04d'%index


def skyplot(crec, title='', axes=None, fignum=30, ait_kw={}, **kwargs):
    """ make an AIT skyplot of a HEALpix array
    crec : array
        must be sorted according to the HEALpix index
    title : string
        set the figure title
    ait_kw : dict
        to set kwargs for image.AIT, perhaps pixelsize
    
    Other args passed to imshow
    """
    n = len(crec)
    nside = int(np.sqrt(n/12))
    assert n==12*nside**2, 'wrong length to be healpix array'
    band = Band(nside)
    def skyplotfun(v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        index = band.index(skydir)
        return crec[index]
    if axes is None:
        plt.close(fignum)
        fig = plt.figure(fignum, figsize=(12,6))
    ait=image.AIT(PySkyFunction(skyplotfun) ,axes=axes, **ait_kw)
    ait.imshow(title=title, **kwargs)
    return ait

  
    
#class DataMap(dict):
#    """ Implement a SkyFunction that returns binned data"""
#    
#    def __init__(self, roi, nside=12*32, emin=1000):
#        """ roi: an ROIAnalysis object
#            nside: int 
#                HEALPix nside
#        """
#        self.nside = nside
#        self.emin = emin
#        self.data = dict()
#        self.bindex = Band(nside).index
#        for band in roi.bands:
#            if band.e <emin: continue
#            # cut on class? t = band.b.event_class() & 3 # note: mask off high bits!
#            for wsd in band.wsdl:
#                index = self.bindex(wsd)
#                self[index] = self.get(index,0) + int(wsd.weight())
#                
#    def __call__(self, v, isskydir=False):
#        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2])) if not isskydir else v
#        index = self.bindex(skydir)
#        return self.get(index,0) 
#
class KdeMap(dict):
    """ Implement a SkyFunction that returns KDE data for a given ROI"""
    
    def __init__(self, roi, nside=12, factor=32, **kwargs):
        """ roi: an ROIAnalysis object
            nside: int 
                HEALPix nside
            factor: int
                factor defining nside of subpixel
        """
        self.nside = nside
        # invoke code in like/roi_tsmap that will create the full map
        hr = roi_tsmap.ROIWrapper(roi,nside)
        self.kde = roi_tsmap.HealpixKDEMap(hr,factor=factor, **kwargs)
                
    def __call__(self, v, isskydir=False):
        #skydir = SkyDir(Hep3Vector(v[0],v[1],v[2])) if not isskydir else v
        return self.kde.call2(v)
        #return self.kde.multipix_call(skydir)

class ResidualTS(object):
    """ manage a residual TS plot, using roi_tsmap.TSCalc 
        and image.ZEA to make a plot of the "residual" TS
    """

    def __init__(self, roi, **kwargs):
        """
        roi : a ROIanalysis object
        kwargs:
            index : float
                default 2.0 the photon index to use
            nside : int
                default 12, for the HEALpix nside parameter
        """
        self.roi = roi
        self.center = roi.roi_dir
        self.tsfun = roi_tsmap.TSCalc(roi, photon_index=kwargs.pop('index', 2.0))
        self.band = Band(kwargs.pop('nside',12))
        self.hp_index = self.band.index(self.center)
        self.setup_zea(**kwargs)
        
    def get_pyskyfun(self):
        return PySkyFunction(self)
    
    def inside(self, skydir):
        return self.band.index(skydir)==self.hpindex
        
    def __call__(self, v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        if not self.zea.inside(skydir): return np.nan
        if self.band.index(skydir) != self.hp_index: return np.nan
        return self.tsfun(skydir)
        
    def setup_zea(self, **kwargs):
        if 'size' not in kwargs.keys(): kwargs['size']= 8
        if 'pixelsize' not in kwargs.keys(): kwargs['pixelsize']=0.10
        self.zea = image.ZEA(self.center, **kwargs)
        
    def fill_zea(self):
        print 'filling %d pixels...' % (self.zea.size/self.zea.pixelsize)**2 
        self.zea.fill(self.get_pyskyfun())
        
    def show_all(self, vmin=0, vmax=50):
        self.zea.imshow(vmin=vmin,vmax=vmax)
        self.zea.grid(color='grey')
        self.zea.colorbar()
        self.show_sources()
        self.zea.axes.set_title(self.roi.name)
        plt.draw_if_interactive()
        
    def show_sources(self):
        for s in self.roi.psm.point_sources:
            if not self.band.index(s.skydir)==self.hp_index: break # assume inside guys first
            self.zea.plot_source(s.name, s.skydir, color='grey')
        
    def savefig(self, outdir):
        if not os.path.exists(outdir): os.mkdir(outdir)
        outfile = os.path.join(outdir, self.roi.name.replace(' ', '_').replace('+','p')+'_tsmap.png')
        self.zea.axes.figure.savefig(outfile)
 
       
def process(roi, outdir='tsmap', **kwargs):
    plt.figure(10); clf();
    rts = ResidualTS(roi, axes=plt.gca(), **kwargs)
    rts.fill_zea()
    rts.show_all()
    rts.savefig(outdir)

def process_table(skyfun,name, pos_list, outfile=None, **kwargs):
    print 'filling table %-5s with %d entries...' % (name, len(pos_list)),
    skytable = np.array([skyfun(p) for p in pos_list])
    print ' min=%6.1f, max=%6.1f, mean=%6.1f ' \
        % (skytable.min(), skytable.max(),skytable.mean()) ,
    if outfile is not None:
        print '--> %s' % outfile
        pickle.dump(skytable, open(outfile,'wb'))
    else: print
            
def make_index_table(nside, subnside):
    band, subband = Band(nside), Band(subnside)
    npix, nsubpix = 12*nside**2, 12*subnside**2
    t=np.array([band.index(subband.dir(i)) for i in xrange(nsubpix)])
    a = np.arange(nsubpix)
    index_table = [a[t==i] for i in xrange(npix)]
    return index_table
 
   
    
def load_tables(name, outdir, nside=12, subnside=12*32):
    files =glob.glob(os.path.join(outdir, '%s_table'%name,'*.pickle'))
    nf = len(files)
    assert nf>0, 'no pickle files found in %s' % os.path.join(outdir, '%s_table'%name)
    if nf<1728: print 'warning: missing %d files' % (1728-nf)
    files.sort()

    r = np.zeros(12*subnside**2)
    r.fill(np.nan)
    pklist = [pickle.load(open(f)) for f in files]
    i12 = [int(f[-11:-7]) for f in files]
    index_table = make_index_table(nside, subnside)
    for index, pk in zip(i12,pklist):
        indeces = index_table[index]
        for i,v in enumerate(pk):
            r[indeces[i]]=v
            
    return r
def dump_tstable(tstable):
    pickle.dump(tstable, open('tstable%02d.pickle'%version, 'wb'))
    

class ROItables(object):
    """ manage one or more tables of values subdividing a HEALpix roi
    
        skyfuns : list of 3-tuples
            the 3-tuples must have: (skyfunction, tablename, dict)
            Implemented now, and default here:
                (roi_tsmap.TSCalc, 'ts', dict(photon_index=2.0),) , 
                (KdeMap,    'kde', dict()),
    """
    
    def __init__(self, outdir, **kwargs):
        nside = kwargs.pop('nside', 12)
        subnside = kwargs.pop('subnside', 384)
        self.make_index_table(nside,subnside)
        self.subdirfun = Band(subnside).dir
        self.skyfuns = kwargs.pop('skyfuns', 
              ( (roi_tsmap.TSCalc, 'ts', dict(photon_index=2.0),) , 
                (KdeMap,    'kde', dict()),
              ),
            )
        self.subdirs = [os.path.join(outdir, name+'_table') for s, name, kw in self.skyfuns]
        for subdir in self.subdirs: 
            if not os.path.exists(subdir):  os.makedirs(subdir)
            
    def make_index_table(self, nside, subnside):
        band, subband = Band(nside), Band(subnside)
        npix, nsubpix = 12*nside**2, 12*subnside**2
        t=np.array([band.index(subband.dir(i)) for i in xrange(nsubpix)])
        a = np.arange(nsubpix)
        self.index_table = [a[t==i] for i in xrange(npix)]
        
    def process_table(self, skyfun,name, pos_list, outfile=None, **kwargs):
        print 'filling table %-5s with %d entries...' % (name, len(pos_list)),
        skytable = np.array([skyfun(p) for p in pos_list])
        print ' min=%6.1f, max=%6.1f, mean=%6.1f ' \
            % (skytable.min(), skytable.max(),skytable.mean()) ,
        if outfile is not None:
            print '--> %s' % outfile
            pickle.dump(skytable, open(outfile,'wb'))
        else: print    
    
    def __call__(self, roi):
        index = int(roi.name[5:])
        pos_list = [self.subdirfun(int(i)) for i in self.index_table[index]]
                
        for i,fun in enumerate(self.skyfuns):
            self.process_table(fun[0](roi, **fun[2]), fun[1], pos_list, 
                os.path.join(self.subdirs[i], roi.name+'.pickle'))
 
class Display_map(object):

    def __init__(self, table, outdir,
                map_dir = None,
                nside=12, subnside=12*32,
                imshow_kw=dict(interpolation='bilinear', vmin=0, vmax=100, ),
                **kwargs):
        self.v = pickle.load(open(table))
        self.subband = Band(subnside)
        self.band = Band(nside)
        self.n = 12*nside**2
        self.imshow_kw=imshow_kw
        self.scale = kwargs.pop('scale', lambda x: x)
        if type(self.scale) == types.StringTypes:
            if self.scale=='sqrt': self.scale= lambda x: np.sqrt(max(x,0))
            elif self.scale=='log': self.scale=lambda x: np.log10(max(x,0.1))
            else:
                raise Exception, 'unrecognized scale function, %s' %self.scale
        print 'Can generate %d map figures' %(self.n)
        self.outdir = outdir
        self.ZEA_kw = kwargs.pop('ZEA_kw', dict())
        if map_dir is not None:
            self.map_path = os.path.join(outdir,map_dir) 
            if not os.path.exists(self.map_path):
                os.makedirs(self.map_path)
            print 'will save figures in folder %s' % self.map_path
        else: self.map_path = None
        self.sources = skymodel.SkyModel(outdir).sources
        print 'loaded %d sources from skymodel %s' % (len(self.sources),outdir)
         
    def get_pyskyfun(self):
        return PySkyFunction(self)

    def skyfun(self, v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        return self.v[self.subband.index(skydir)]
        
    def __call__(self,v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        t =self.v[self.subband.index(skydir)]
        #if   self.scale=='sqrt': return np.sqrt(max(t,0))
        #elif self.scale=='log':  return np.log10(max(t, 1e-1))
        return self.scale(t) 

    def fill_ait(self, fignum=11, axes=None, show_kw={}, **kwargs):
        if axes is None:
            plt.close(fignum)
            fig=plt.figure(fignum, figsize=(12,8));
            axes=fig.gca()
        pixelsize = kwargs.pop('pixelsize', 0.25)
        ait = image.AIT(self.get_pyskyfun(),axes=axes, pixelsize=pixelsize, **kwargs)
        self.imgplot=ait.imshow(**show_kw)
        return ait

    def fill_zea(self, index, fignum=12, axes=None, **kwargs):
        """ sources is a recarray with name, ra, dec
        """
        if axes is None:
            plt.close(fignum)
            fig = plt.figure(fignum,figsize=(6,6)); 
            axes = fig.gca()
        size = kwargs.pop('size', 10)
        pixelsize = kwargs.pop('pixelsize', 0.1)
        title = kwargs.pop('title', hpname(index))
        label = kwargs.pop('label', '')
        zea = image.ZEA(self.band.dir(index), size=size, pixelsize=pixelsize,**self.ZEA_kw)
        zea.grid()
        zea.fill(self.get_pyskyfun())
        imshow_kw = self.imshow_kw #
        imshow_kw.update(kwargs)
        zea.imshow(**imshow_kw)
        zea.colorbar(label=label)
        axes.set_title(title)
        if self.sources is not None:
            count = 0
            for s in self.sources:
                sdir = s.skydir
                if not zea.inside(sdir):continue
                count += 1
                inside =self.band.index(sdir)==index
                zea.plot_source(s.name, sdir, symbol='*' if inside else 'd', 
                    markersize=14 if inside else 8,
                    color='w')
            print 'found %d sources to plot' %count        
        
        if self.map_path is not None:
            fout = os.path.join(self.map_path,hpname(index)+'.png')
            plt.savefig(fout)
            print 'saved figure to %s' % fout
        plt.draw_if_interactive()

class Wrapper(object):
    """ adapt the class that generates a ZEA display for multiprocessing"""
    def __init__(self, title, *pars, **kwargs):
        self.label = kwargs.pop('label')
        self.display = Display_map(*pars, **kwargs)
        self.n = self.display.n
        self.title = title  
    def names(self):
        return map(hpname,range(self.n))
    def __call__(self, index):
        ait=self.display.fill_zea(index, title=self.title+' ('+hpname(index)+')', label=self.label)
  