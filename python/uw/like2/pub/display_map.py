"""
Manage generation of maps from HEALpix tables
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pub/display_map.py,v 1.3 2011/06/24 04:53:06 burnett Exp $
"""
import os,sys, types, pickle
import numpy as np
import pylab as plt
from uw.utilities import image
from skymaps import Band, SkyDir, PySkyFunction, Hep3Vector, SkyImage

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

class DisplayMap(object):
    """ display the contents of a HEALpix table as ait or zea
    
    """
    def __init__(self, table,
                sources=None, 
                imshow_kw=dict(interpolation='bilinear',  ),
                **kwargs):
        """table : string or iterable
                If a string, the name of a pickled file
           sources : None or a string 
                if a string, the name of a pickled rec with name, ra, dec fields
        """
        if type(table)==types.StringType:        
            self.v = pickle.load(open(table))
            print ('Loaded HEALpix table from file %s' %table)
        else: self.v=table
        self.nside = int(np.sqrt(len(self.v)/12))
        assert len(self.v)==12*self.nside**2, 'size of map not consistent with expected nside %d' % nside 
        self.band = Band(self.nside)
        self.imshow_kw=imshow_kw
        self.scale = kwargs.pop('scale', lambda x: x)
        if type(self.scale) == types.StringTypes:
            if self.scale=='sqrt': self.scale= lambda x: np.sqrt(max(x,0))
            elif self.scale=='log': self.scale=lambda x: np.log10(max(x,0.1))
            else:
                raise Exception('unrecognized scale function, %s' %self.scale)
        self.ZEA_kw = kwargs.pop('ZEA_kw', dict(galactic=True, size=10, pixelsize=0.1))
        if sources is not None:
            self.sources = pickle.load(open(sources))
            print ('loaded %d sources from %s' % (len(self.sources),sources))
        else:self.sources=None
        
        self.map_path = kwargs.pop('map_path',None)
        
    def get_pyskyfun(self):
        return PySkyFunction(self)

    def skyfun(self, v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        return self.v[self.band.index(skydir)]
        
    def __call__(self,v):
        skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
        t =self.v[self.band.index(skydir)]
        return self.scale(t) 

    def fill_ait(self, fignum=11, axes=None, show_kw={}, source_kw={}, figwidth=12, margin=0.15, **kwargs):
        if axes is None:
            # set up a figure for 2x1 image with equal margins
            plt.close(fignum)
            figheight = figwidth*(1.+2*margin)/(1+margin)/2.
            fig=plt.figure(fignum, figsize=(figwidth, figheight));
            axes=plt.gca()
            plt.subplots_adjust(left=0.05, right=0.95) #gives reasonable equal margins
        pixelsize = kwargs.pop('pixelsize', 0.25)
        ait = image.AIT(self.get_pyskyfun(),axes=axes, pixelsize=pixelsize, **kwargs)
        self.imgplot=ait.imshow(**show_kw)
        ait.axes.set_autoscale_on(False)

        if self.sources is not None:
            sdirs = map(SkyDir, self.sources.ra, self.sources.dec)
            ait.plot(sdirs, **source_kw)
            print ('found %d sources to plot' % len(sdirs) )
        plt.draw_if_interactive()
        return ait

    def fill_zea(self, index, fignum=12, axes=None, show_kw=None, **kwargs):
        """ index: integer, or a SkyDir
                the HP12 index if integer
            figmun: integer
                used if axes is None
            show_kw : dict
                override imshow keywords
            kwargs 
                size
                pixelsize
                galactic
        """
        if axes is None:
            plt.close(fignum)
            fig = plt.figure(fignum,figsize=(6,6)); 
            axes = fig.gca()
        if type(index) == types.IntType:
            sdir = Band(12).dir(index)
            title = 'HP12_%4d'%index
        else:
            sdir = index
            title = 'l = %.1f, b=%.1f' % (sdir.l(), sdir.b())
        title = kwargs.pop('title',title)
        kw = self.ZEA_kw
        kw.update(kwargs)
        zea = image.ZEA(sdir, **kw)
        zea.grid()
        zea.fill(self.get_pyskyfun())
        zea.imshow( **(show_kw if show_kw is not None else self.imshow_kw))
        zea.colorbar()
        if title is not None: axes.set_title(title)
        if self.sources is not None:
            count = 0
            for s in self.sources:
                sdir = SkyDir(s.ra,s.dec)
                if not zea.inside(sdir):continue
                count += 1
                inside =self.band.index(sdir)==index
                zea.plot_source(s.name, sdir, symbol='*' if inside else 'd', 
                    markersize=14 if inside else 8,
                    color='w')
            print ('found %d sources to plot' %count        )
        
        if self.map_path is not None:
            fout = os.path.join(self.map_path,hpname(index)+'.png')
            plt.savefig(fout, bbox_inches='tight')
            print ('saved figure to %s' % fout)
        plt.draw_if_interactive()
        return zea

class SourceDensity(object):
    """ create source density HEALpix array from a list of locations
    """
    def __init__(self, nside=12):
        """
        nside: integer
            the HEALpix nside parameter
        """
        self.v = np.zeros(12*nside**2, float)
        self.index = Band(nside).index
        
    def fill(self, sdirs):
        """ sdirs: a list of SkyDir objects
        """
        for s in sdirs:
            self.v[self.index(s)]+=1
            
    def fill_rec(self, rec, cut=None):
        """ rec: a recarry with ra, dec columns
            cut : None or a mask arrray
            
        """
        if cut is None:
            sdirs = map(SkyDir, rec.ra, rec.dec)
        else:
            sdirs = map(SkyDir, rec.ra[cut], rec.dec[cut])
        self.fill(sdirs)
        
    def save(self, fn):
        pickle.dump(self.v, open(fn, 'wb'))
        print ('saved file %s' % fn)

class SourceMap(DisplayMap):
    """ subclass of DisplayMap to display point source positions on a photon density map
    """
    def __init__(self,  kde,
                sources ,
                show_kw=dict(fun = lambda x:np.sqrt(x/1e6), vmax=4, cmap='hot'),
                plot_kw=dict(nocolorbar=False,),
                pos=None, size=180,
                ):
        
        super(SourceMap,self).__init__(kde)
        if type(sources) == types.StringType:
            self.s = pickle.load(sources)
            print ('loaded %5d sources from %s' %(len(self.s), fn))
        else: self.s = sources

        self.show_kw = show_kw
        
    def fill_ait(self, fignum=20, axes=None, **kwargs):
        ait = super(SourceMap, self).fill_ait( fignum=fignum, axes=axes, show_kw= self.show_kw, **kwargs)
        ait.axes.set_autoscale_on(False) # prevent rescaling when adding points
        self.ait=ait
        return ait
    
    def fill_zea(self, pos, fignum=21, axes=None, which=-1, savefn=None, **kwargs):
        sfactor = kwargs.pop('sfactor', 1)
        zea = super(DMap, self).fill_zea(pos, fignum=fignum, axes=axes, show_kw= self.show_kw, **kwargs)
        s = self.s
        for subset, marker, color, size, label in self.subsets(s, which):
            zea.plot(map(SkyDir, s.ra[subset], s.dec[subset]), edgecolor='grey',
                marker=marker, c=color, s=size*sfactor, label=label)
            print ('plotted %4d sources, subset "%s"' %(sum(subset), label))
        plt.legend(scatterpoints=1, loc=2)
        if savefn is not None:
            self.savefig(savefn % i); i+=1

        return zea

    def legend(self):
        plt.legend(frameon=False,scatterpoints=1, loc=(-0.05,-0.05))
    
    def savefig(self, fn):
        plt.savefig(fn, bbox_inches='tight', pad_inches=0, dpi=160)
    
    def subsets(self, s, which):
        assoc = s.id_prob>0.8
        ts25=s.ts>=25
        lt25=(s.ts<25) 
        t =(((-assoc)*(lt25),'+', 'grey',   8, 'no id, TS<25'),
            ((-assoc)*(ts25), 's', 'red',   10, 'no id, TS>25'),
            (assoc,           'o', 'green', 12, 'associated'  ),
            )
        return t if which <0 else (t[which],)

    def add_sources(self, which=-1, sfactor=1):
        s = self.s
        print ('loaded %5d sources' %(len(s),))
        i=0 if which<0 else which+10
        plt.rcParams['legend.fontsize']= 8.0
        for subset, marker, color, size, label in self.subsets(s, which):
            self.ait.plot(map(SkyDir, s.ra[subset], s.dec[subset]), edgecolor='grey',
                marker=marker, c=color, s=size*sfactor, label=label)
            print ('plotted %4d sources, subset "%s"' %(sum(subset), label))
            self.legend()


def load_skyspect(fn = r'T:\data\galprop\ring_21month_P6v11.fits', 
# r'D:\fermi\data\galprop\gll_iem_v02.fit', 
        nside=192, 
        show_kw = dict(fun=np.log10, cmap='hot'),
        ):
    """
    load a galactic diffuse distribution.
    Save the HEALpix respresentation at an energy (1 GeV default)
    
    fn : string
        filename for the FITS representaion of a  SKySpectrum  
    nside: int
        HEALpix nside to use for represenation -- note that 192 is 12*16, about 0.25 deg
    show_kw : dict
        fun: weighting function, cmap, vmin, vmax
    """
    t = SkyImage(fn)
    galname = os.path.split(fn)[-1]
    print ('%s: nx, ny, layers: %d %d %d' %(galname, t.naxis1(), t.naxis2(), t.layers()))
    hpdir = Band(nside).dir
    dmap = map(lambda i:t(hpdir(i)), xrange(12*nside**2))
    tdm=DisplayMap(dmap)
    tdm.fill_ait(fignum=12, source_kw=dict(edgecolor='w',), show_kw=show_kw )
    plt.title(galname+' (1 GeV)')
    sfn = galname.split('.')[0]+'.png'
    plt.savefig(galname.split('.')[0]+'.png', bbox_inches='tight', pad_inches=0)
    print ('saved figure to %s' % sfn)
    return tdm