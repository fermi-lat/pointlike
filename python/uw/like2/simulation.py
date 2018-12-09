"""

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/simulation.py,v 1.2 2018/01/27 15:37:17 burnett Exp $

"""
import os, sys,  types, glob
import cPickle as pickle
import numpy as np
import pandas as pd
from skymaps import Band
from ..data import binned_data
from . import configuration
from astropy.io import fits

def make_index_table(nside=12, subnside=512, usefile=True):
    """create, and/or use a table to convert between different nside pixelizations
    """
    filename = os.path.expandvars('$FERMI/misc/index_table_%02d_%03d.pickle' % (nside, subnside) )
    if os.path.exists(filename) and usefile:
        return pickle.load(open(filename))
    print 'generating index table for nside, subnside= %d %d' % (nside, subnside)
    band, subband = Band(nside), Band(subnside)
    npix, nsubpix = 12*nside**2, 12*subnside**2
    t=np.array([band.index(subband.dir(i)) for i in xrange(nsubpix)])
    a = np.arange(nsubpix)
    index_table = [a[t==i] for i in xrange(npix)]
    if usefile:
        pickle.dump(index_table, open(filename,'w'))
    return index_table
    
def default_geom():
    """return a DataFrame with indexed by the band, containing emin, emax, event_type, nside
    """
    energies = np.logspace(2, 6, 17) # 100 MeV to 1Tev, 4/decade
    elow = energies[:-1]
    ehigh=energies[1:]
    nside_list = ModelCountMaps.nside_array
    geom = dict()
    for ib in range(32):
        ie = ib/2; it=ib%2
        geom[ib]=dict(elow=elow[ie],ehigh=ehigh[ie], 
                      event_type=it, nside=nside_list[ib])
    g = pd.DataFrame(geom).T['elow ehigh event_type nside'.split()]
    g.event_type = g.event_type.astype(int)
    g.nside = g.nside.astype(int)
    return g

class ModelCountMaps(object):
    """ This makes, and saves, a list of predicted counts for all the pixels within the central
    HEALPix pixel.
    """
    nside_array = np.ones(32,int) * 1024
    nside_array[:10] = 64,64, 128,64 ,256,128, 256,256, 512,512

    def __init__(self, roi, nbands=None, bandlist=None,  subdir='model_counts'):
        """
        roi : ROI object
            Expect it to behave like an array, indexed by the band index,
            and each element is a SkyFunction, returning the count density.
        bandlist : list of int | None
            bands to process. If None, do them all

        subdir : string | None
            folder name to write results to
        """
        roi_index = Band(12).index(roi.roi_dir)  

        if bandlist is None:
            if nbands is not None: bandlist=range(nbands)
            else: bandlist = range(len(roi))
        for ebi in bandlist:
            eb = roi[ebi] #EnegyBand object

            # choose nside to be power of 2, override the data
            nside =  ModelCountMaps.nside_array[ebi]
            dirfun = Band(nside).dir
            pixel_area = Band(nside).pixelArea()
            index_table = make_index_table(12, nside)
            pix_ids = index_table[roi_index]
            dirs = map(dirfun, pix_ids)
            cnts = np.array(map(eb, dirs),np.float32) * pixel_area
            
            print '{:4d} {:4d} {:6d} {:8.2e} {:8.2e} {:8.2e}'.format(
                ebi,nside,len(cnts), cnts.mean(), cnts.min(), cnts.max())
            if subdir is not None:
                subsubdir = subdir+'/{}'.format(ebi)
                if not os.path.exists(subsubdir): os.makedirs(subsubdir)
                outfile = subsubdir+'/HP12_{:04d}.pickle'.format(roi_index)
                pickle.dump(cnts, open(outfile, 'w'))
                print '\t\t--> {}'.format(outfile)

class BandCounts(object):
    """Manage results of the Model counts for a Band
    """
    def __init__(self, band_index, path='model_counts', reload=False):
        """Assume in a skymodel folder, containg a model_counts subfolder
        """
        self.bi = band_index
        self.path=path+'/{:02d}'.format(band_index)
        nside_list = ModelCountMaps.nside_array
        self.nside = nside_list[band_index]
        self.filename = self.path+'/combined.pickle'
        if os.path.exists(self.filename) and not reload:
            self.load()
        else:
            self.combine()
            self.dump()

    def combine(self):
        index_table = make_index_table(12, self.nside)
        d = dict()
        ff = sorted(glob.glob(self.path+'/HP12*.pickle'.format(self.bi)));
        assert len(ff)==1728, 'only found {} pickle files in {}'.format(len(ff), self.path)
        for i,f in enumerate(ff):
            ids = index_table[i]
            try:
                values = pickle.load(open(f))
            except Exception, msg:
                print 'Failed to load file {}: {}'.format(f, msg)
                raise
            assert len(ids)==len(values), 'oops: {} ids, but {} values'.format(len(ids), len(values))
            d.update(zip(ids, values))
        self.counts= np.array(d.values(),np.float32)

    def plot(self, **kwargs):
        from uw.like2.pub import healpix_map as hpm
        name = 'band{:02d}'.format(self.bi)
        hpm.HParray(name, self.counts).plot(log=True, title=name, **kwargs)

    def dump(self):
        pickle.dump(self.counts, open(self.filename, 'w'))
        print 'Saved file {}'.format(self.filename)

    def load(self):
        self.counts = pickle.load(open(self.filename))
        
    def simulate(self):
        """Return sparsified Poisson simulation
        """
        sim =np.random.poisson(self.counts)
        nonzero = sim>0
        ids=np.arange(len(sim))[nonzero]
        counts = sim[nonzero]
        return np.array([ids, counts], np.int32)

class SimulatedPixels(binned_data.Pixels):
    """Generate simulated pixel data, using pixel-based count predictions
    Inherits from the Pixels class in binned_data to export the simuulation to a sparse FITS representation
    """
    def __init__(self, model_path='.', subfolder='model_counts', numchan=None):
        """
        model_path : str
            expect to find a sky model folder, containing itself a folder with 
        """
        self.countsfolder=os.path.join(model_path, subfolder)
        assert os.path.exists(self.countsfolder), 'did not find folder {}'.format(countsfolder)
        band_folders = sorted(glob.glob(self.countsfolder+'/*'))
        if numchan is None:
            numchan=len(band_folders)
        sim = dict()
        print 'Simulating from model predictions in\n  {}\n  chan   pixels    counts'.format(
            os.path.abspath(self.countsfolder))
        for f in  band_folders[:numchan]:
            channel = int(os.path.split(f)[-1])
            print '{:6d}'.format(channel), 
            sim[channel]= s= BandCounts(channel).simulate()
            print '{:8d} {:9d}'.format(s.shape[1], sum(s[1,:]))
            
        # set these to be consistent with base class    
        keys = sorted(sim.keys())
        self.cnt = np.hstack([sim[k][1,:] for k in keys]); 
        self.pix = np.hstack([sim[k][0,:] for k in keys]);
        self.chn = np.hstack([np.ones(sim[k].shape[1], dtype=np.int16)*k for k in keys])
        self.counter=None
        self._sorted=True
        channels = sorted(list(set(self.chn)))
        indexchan = list(np.searchsorted(self.chn, channels))+[len(self.chn)]
        self.lookup = dict(zip(channels,zip(indexchan[:-1], indexchan[1:])))

class DefaultBands(binned_data.BandList):
    """Define a BANDS table corresponding to traditional pointlike setup, but with power-of-2 nside values

    It inherits from Bandlist in binned data, allowing creation of a FITS HDU
    """
    def __init__(self, channels=None, emin=(100,100)):
        """
        channels : list of int | None
            if a list, the channel IDs to include
        emin : 2-tuple of float
            minimum energy for front or back 
        """
        if channels is None:
            channels = range(32)
        energies = np.logspace(2, 6, 17) # 100 MeV to 1Tev, 4/decade
        elow = energies[:-1]
        ehigh=energies[1:]
        nside_array = np.ones(32,int) * 1024
        ### This was designed for F,B,F,B,...
        nside_array[:10] = 64,64, 128,64 ,256,128, 256,256, 512,512
        geom = dict()
        ic=0
        for ib in range(32): 
            ie = ib/2; it=ib%2
            if elow[ie]<emin[it] or ic not in channels: continue # skip according to front or back emin
            geom[ic]=dict(e_min=elow[ie],e_max=ehigh[ie], 
                          event_type=it, nside=nside_array[ic])
            ic +=1
        g = pd.DataFrame(geom).T['e_min e_max event_type nside'.split()]
        g.event_type = g.event_type.astype(int)
        g.nside = g.nside.astype(int)
        self.df=g

    def dataframe(self):
        return self.df

class Simulate(binned_data.BinFile):
    """
    """
    def __init__(self, model_path='.', subfolder='model_counts', numchan=None):
        # get the GTI from the original file? Do we need it?
        config = configuration.Configuration('.', postpone=True, quiet=True)
        bf = config.dataset.binfile
        self.hdus = hdus = fits.open(bf)
        self.hdu0 = hdus[0]
        self.gti = binned_data.GTI(hdus['GTI'])

        self.pixels=SimulatedPixels(model_path, subfolder, numchan=numchan)
        self.bands = DefaultBands(self.pixels.lookup.keys(), 
                emin=config['input_model'].get('emin',(100,100)))

    def writeto(self, filename, clobber=False):
        """ Override base class to avoid clobber existing files
        """
        super(Simulate, self).writeto(filename, clobber)

