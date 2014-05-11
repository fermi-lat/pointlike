"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/phasedata.py,v 1.7 2013/04/10 01:06:24 kerrm Exp $

Handle loading of FT1 file and phase folding with polycos.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""

import numpy as np
import pyfits
import sys
from uw.utilities import keyword_options
from timeman import METConverter

def rad_mask(ras,decs,ra,dec,radius,mask_only=False):
    # make a cut on radius
    ra0,dec0 = np.radians(ra),np.radians(dec)
    ras,decs = np.radians(ras),np.radians(decs)
    cos_diffs = np.sin(decs)*np.sin(dec0)+np.cos(decs)*np.cos(dec0)*np.cos(ras-ra0)
    mask = cos_diffs > np.cos(np.radians(radius))
    return mask

class PhaseData(object):
    """ Object to manage loading data. """

    defaults = (
        ('emin',100,'(MeV) minimum energy to accept'),
        ('emax',30000,'(MeV) maximum energy to accept'),
        ('tmin',0,'Minimum time (MET) to use; 0 is no cut'),
        ('tmax',0,'Maximum time (MET) to use; 0 is no cut'),
        ('rmax',0,'Maximum radius (deg) to use; 0 is no cut.'),
        ('we_col_name','WEIGHT','name for column with weight data in FT1 file'),
        ('use_weights',False,'load weights if available'),
        ('wmin',0,'minimum weight to use'),
        ('ft2',None,'an FT2 file that can be used for on-the-fly time correction'),
        ('bary',False,'use barycenter time instead of geocenter'),
        ('timecol','TIME','name of the time column in the FT1 file'),
        ('noprocess',False,'do not geo- or barycenter times in FT1 file'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,ft1file,polyco,**kwargs):
        keyword_options.process(self,kwargs)
        self.ft1file = ft1file
        self.polyco = polyco
        self.process_ft1()

    def process_ft1(self):
        mc = self.mc = METConverter(self.ft1file,ft2=self.ft2,
            ra=self.polyco.ra,dec=self.polyco.dec,bary=self.bary,
            noprocess=self.noprocess)
        self.mjd_start = mc.MJDSTART; self.mjd_stop = mc.MJDSTOP
        f    = pyfits.open(self.ft1file)
        ens  = np.asarray(f['EVENTS'].data.field('ENERGY'))
        mask = (ens >= self.emin) & (ens < self.emax)
        tcol = self.timecol
        if self.tmin > 0:
            mask &= np.asarray(f['EVENTS'].data.field(tcol)) >= self.tmin
        if self.tmax > 0:
            mask &= np.asarray(f['EVENTS'].data.field(tcol)) < self.tmax
        if self.rmax > 0:
            ra = self.polyco.ra; dec = self.polyco.dec
            if (ra is None) or (dec is None):
                raise ValueError('Cannot make cut without position!')    
            ras = np.asarray(f['EVENTS'].data.field('RA'))
            decs = np.asarray(f['EVENTS'].data.field('Dec'))
            mask &= rad_mask(ras,decs,ra,dec,self.rmax)
        if self.use_weights:
            weights = np.asarray(f['EVENTS'].data.field(self.we_col_name))
            mask = mask & (weights > self.wmin)
            self.weights = weights[mask]
        else:
            self.weights = None
        mets = np.asarray(f['EVENTS'].data.field(tcol))[mask]
        self.mjds = mc(mets)
        self.ph = self.polyco.vec_evalphase(self.mjds)
        try:
            self.orig_ph = f['EVENTS'].data.field('PULSE_PHASE')[mask]
        except KeyError:
            pass
        f.close()

        print >>sys.stderr, "PhaseData: Cuts left %d out of %d events." % (mask.sum(), len(mask))

    def write_phase(self,col_name='PULSE_PHASE'):
        f = pyfits.open(self.ft1file)
        mjds = self.mc(np.asarray(f['EVENTS'].data.field(self.timecol)))
        ph = self.polyco.vec_evalphase(mjds)
        
        try:
            # clobber old values if there
            f['EVENTS'].data.field(col_name)[:] = ph
            print 'Clobbered old %s column.'%(col_name)
        except KeyError:
            c    = pyfits.Column(name=col_name,format='D',array=ph)
            cols = f['EVENTS'].columns
            newcols = cols.add_col(c)
            t = pyfits.new_table(newcols,header=f['EVENTS'].header)
            t.name = 'EVENTS'
            f[1] = t
        f.writeto(self.ft1file,clobber=True)

    def toa_data(self,mjd_start,mjd_stop):
        mask = (self.mjds >= mjd_start)&(self.mjds < mjd_stop)
        if self.weights is None:
            return self.ph[mask],None
        return self.ph[mask],self.weights[mask]
