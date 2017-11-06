"""
Manage the sparse pointlike data. 
Duplicates the functionality of the C++ class BinnedPhotonData
Implements the new standard data format
http://gamma-astro-data-formats.readthedocs.io/en/latest/skymaps/healpix/index.html#hpx-bands-table

$Header$

"""
import os, glob, StringIO
from collections import Counter 
import numpy as np
from astropy.io import fits
import pandas as pd

class GTI(object):
    def __init__(self, gti_hdu):
        self.hdu=gti_hdu
        data = gti_hdu.data 
        self.start=data.START
        self.stop =data.STOP

    def add(self, other):
        self.start=np.hstack([self.start, other.start])
        self.stop =np.hstack([self.stop, other.stop])
        g = self.hdu
        g.columns['START'].array=self.start
        g.columns['STOP'].array = self.stop

    def make_hdu(self):
        cols = [
            fits.Column('START', format='D', unit='s',array=self.start),
            fits.Column('STOP', format='D', unit='s', array=self.stop)
            ]
        return fits.BinTableHDU.from_columns(cols,header=self.hdu.header,name='GTI')
 
    def __repr__(self):
        return '{} intervals from {:.0f} to {:.0f}, {:,.0f} s'.format(len(self.start), 
            self.start[0], self.stop[-1], sum(self.stop-self.start))


class BandList(object):
    def __init__(self, bands_hdu):
        self.hdu=bands_hdu
        self.bands=np.asarray(bands_hdu.data)
        self.version = bands_hdu.header['VERSION']
        if self.version<3:
            self.pixelcnt = np.array(map( lambda n: 0 if n<0 else n,
                 bands_hdu.data.field('COUNT')),int)
            if 'NBRBANDS' not in self.hdu.header:
                print 'Warning: header in BANDS HDU missing cards'
        else:
            self.pixelcnt=None
     
    def __repr__(self):
        df = self.dataframe()
        return '{} bands from {:.0f} to {:.0f} MeV'.format(
            len(df), df.e_min[0],df.e_max.max())

    def __getitem__(self, index):
        df = self.dataframe()
        return df.iloc[index]
        
    def make_hdu(self, photons=None, version=3):
        """
        """
        # new format
        df = self.dataframe()
        band_cols = [
            fits.Column(name='NSIDE', format='J', array=df.nside),
            fits.Column(name='E_MIN', format='D', array=df.e_min*1e-3, unit='keV'),
            fits.Column(name='E_MAX', format='D', array=df.e_max*1e-3, unit='keV'),
            fits.Column(name='EVENT_TYPE', format='J', array=df.event_type),
        ]
        bands_hdu=fits.BinTableHDU.from_columns(band_cols, name='BANDS')
        bands_hdu.header.update(VERSION=version)
        return bands_hdu

    def dataframe(self):
        data = self.hdu.data
        if self.version>2:
            cdata = [data.E_MIN*1e3, data.E_MAX*1e3, data.EVENT_TYPE, data.NSIDE]
        else:
            cdata = [self.hdu.data.field(cname) for cname in 'emin emax event_class nside'.split()]
        df = pd.DataFrame(cdata, index='e_min e_max event_type nside'.split()).T
        df['event_type']= df.event_type.astype(int)
        df.nside = df.nside.astype(int)
        return df
    
class Pixels(object):
    """The list of pixels
    Each line is a pixel, with a HEALPix index, a channel index, and the number of photons in the pixel 
    """
    def __init__(self, pixel_hdu, pixel_count=None):
        """pixel_hdu : HDU 
           pixel_count : array of int | None
                number of pixels per band, from the band list in old format
        """
        self.hdu = pixel_hdu
        pixeldata = pixel_hdu.data
        if pixel_count is not None:
            # old format: make a list of channel numbers from the bands HDU
            # the old list of pixels was sorted by channel
            chn= []
            for i,c in enumerate(pixel_count):
                if c<0: continue
                chn = chn + [i]*c 
            self.pix = pixeldata.field('INDEX')
            self.cnt = pixeldata.field('COUNT')
            self.chn = chn
        else:
            # read new format
            self.chn = pixeldata.field('CHANNEL')  # band index
            self.pix = pixeldata.field('PIX')     # pixel index (depends on channel)
            self.cnt = pixeldata.field('VALUE')   # number of photons in bin
        
        self.counter = None
        self._sorted = False

    def _make_counter(self):
        # Make a Counter, with keys combined from channel and pixel id
        if self.counter is None:
            chn = np.array(self.chn,int) # convert to int for shift
            keys = list(np.left_shift(32) + self.pix)
            self.counter = Counter(dict(zip(keys,cnt)))
        return self.counter
 
    def add(self, other):
        """combind the current list of pixels with another
        other : Pixels object
        """
        # combine the pixels by updating the Counter dictionary-like objects
        self._make_counter() # create local Counter only if adding another
        self.counter.update(other._make_counter())
    
    def _decode_counter(self):
        # Update the three arrays following adding another set of data
        assert self.counter is not None, 'logic failure'
        items = np.array(self.counter.items(), int)
        keys = items[:,0]
        self.chn = np.right_shift(keys,32)
        self.pix = np.bitwise_and(keys, 2**32-1)
        self.cnt = items[:,1]

    def dataframe(self):
        """return a DataFrame with number of pixels and photons per band
        """
        if self.counter is not None:
            self._decode_counter()
        channels = sorted(list(set(self.chn))); 
        d = dict()
        for channel in channels:
            c = self.cnt[self.chn==channel]
            d[channel] = {'pixels': len(c), 'photons': sum(c)}
        df = pd.DataFrame(d).T[['pixels', 'photons']]
        return df

    def __getitem__(self, channel):
        """return a list of (pixel, count) pairs for the band 
        """
        if not self._sorted:
            if self.counter is not None:
                self._decode_counter()
            # sort the list of pixels according to channel number (band)
            # create a lookup dictionary with limits for the pixel and count lists
            csort = self.chn.argsort()
            self.chn = self.chn[csort]
            self.pix = self.pix[csort]
            self.cnt = self.cnt[csort]
            channels = sorted(list(set(self.chn)))
            indexchan = list(np.searchsorted(self.chn, channels))+[len(self.chn)]
            self.lookup = dict(zip(channels,zip(indexchan[:-1], indexchan[1:])))
            self._sorted = True
        try:
            a,b = self.lookup[channel]
            return zip(self.pix[a:b], self.cnt[a:b])
        except KeyError:
            return [] # empty list of no entry

    def make_hdu(self):
        """ create a new HDU in new format
            
        """
        if self.counter is not None:
            # needed it result of combining
            self._decode_counter()
        skymap_cols = [
            fits.Column(name='PIX', format='J',    array=self.pix),
            fits.Column(name='CHANNEL', format='I',array=self.chn),
            fits.Column(name='VALUE', format='J',  array=self.cnt),
        ]
        skymap_hdu=fits.BinTableHDU.from_columns(skymap_cols, name='SKYMAP')
        skymap_hdu.header.update(
            PIXTYPE='HEALPIX',
            INDXSCHM='SPARSE',
            ORDERING='RING',
            COORDSYS='GAL',
            BANDSHDU='BANDS',
            )
        return skymap_hdu
    
    def __repr__(self):
        if self.counter is not None:
            npix, nphot = len(self.counter), sum(self.counter.values())
        else:
            npix, nphot = len(self.cnt), sum(self.cnt)
        return '{:,} pixels, {:,} photons'.format(npix, nphot) 

    
class BinFile(object):
    """ A Binned photon data file
    Manages Two FITS tables: 
        * Pixels
        * BandList
    Implements an indexing interface. Now returns a special Band object
    """
    def __init__(self, filenames, outfile=None, adding=False):
        """
        filenames : a FITS file name, or a list
            if a list, combine them
        outfile : filename | Null [default]
            write the corresponding fits file
        """
        if not hasattr(filenames, '__iter__'):
            filenames = [filenames]
        for i,filename in enumerate(filenames):
            if i==0:
                print '\n"{}" '.format(filename),
                self.hdus=fits.open(filename)
                self.gti=GTI(self.hdus['GTI'])
                if self.hdus[1].name=='BANDS':
                    self.bands=BandList(self.hdus['BANDS'])
                    self.pixels=Pixels(self.hdus['PIXELS'], self.bands.pixelcnt)
                else:
                    self.bands=BandList(self.hdus['BANDS'])
                    self.pixels=Pixels(self.hdus['SKYMAP'], self.bands.pixelcnt)
                print self.pixels.__repr__(),
            else:
                self.add(BinFile(filename, adding=True))
                print self.pixels.__repr__(),
        if not adding: print

        if outfile is not None:
            self.writeto(outfile)

    def fits_info(self):
        output= StringIO.StringIO()
        self.hdus.info(output)
        return output.getvalue()

    def __repr__(self):
        out = self.fits_info()
        out += 'GTI: '+ self.gti.__repr__()  
        out += '\nBands: '+self.bands.__repr__()
        out += '\nPixels: '+self.pixels.__repr__()
        return out

    def __getitem__(self, index):
        """ return a skymaps.Band C++ object corresponding to the band index
        This object implements query_disk for extracting the pixels within a given ROI.
        """
        from skymaps import Band
        b = self.bands[index]
        bb = Band(int(b.nside), int(b.event_type), b.e_min, b.e_max, 0,0)
        # this is an unfortunate loop: need to consider an interface for adding a set of pixels
        for pix, cnt in self.pixels[index]:
            bb.add(int(pix), int(cnt))
        return bb
    def __len__(self): return len(self.bands.bands)

    def add(self, other):
        """ other : BinFile object
        """
        # combine the pixel and GTI arrays
        self.pixels.add(other.pixels)
        self.gti.add(other.gti)

    def dataframe(self):
        """return a DataFrame with Band info and Pixel summary
        """
        dfb = self.bands.dataframe()
        dfp = self.pixels.dataframe()

        df = pd.concat([dfb,dfp], axis=1) 
        df.index.name='band' 
        # this just to replace NaN's for bands missing in pixel list with zeros
        pixcnt = df.pixels.copy()
        gcnt = df.photons.copy()
        missing = pd.isnull(pixcnt)
        pixcnt[missing]=0
        gcnt[missing]=0
        df.pixels=pixcnt.astype(int)
        df.photons=gcnt.astype(int)

        return df

    def writeto(self, filename, clobber=True):
        """write to a file

        """
        gti_hdu = self.gti.make_hdu()
        bands_hdu=self.bands.make_hdu() 
        pixels_hdu = self.pixels.make_hdu()
        hdus=[self.hdus[0], pixels_hdu, bands_hdu, gti_hdu]
        fits.HDUList(hdus).writeto(filename, clobber=clobber)
        print 'wrote file {}'.format(filename)
        
       