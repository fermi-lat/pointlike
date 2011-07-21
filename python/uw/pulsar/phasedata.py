"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/phasedata.py,v 1.1 2011/04/27 18:32:03 kerrm Exp $

Handle loading of FT1 file and phase folding with polycos.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""

import numpy as np
import pyfits
import sys
from uw.utilities import keyword_options
from timeman import METConverter

class PhaseData(object):
    """ Object to manage loading data. """

    defaults = (
        ('emin',100,'(MeV) minimum energy to accept'),
        ('emax',30000,'(MeV) maximum energy to accept'),
        ('we_col_name','WEIGHT','name for column with weight data in FT1 file'),
        ('use_weights',False,'load weights if available'),
        ('wmin',0,'minimum weight to use'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,ft1file,polyco,**kwargs):
        keyword_options.process(self,kwargs)
        self.ft1file = ft1file
        self.polyco = polyco
        self.process_ft1()

    def process_ft1(self):
        mc = METConverter(self.ft1file)
        self.mjd_start = mc.MJDSTART; self.mjd_stop = mc.MJDSTOP
        f    = pyfits.open(self.ft1file)
        ens  = np.asarray(f['EVENTS'].data.field('ENERGY'))
        mask = (ens >= self.emin) & (ens < self.emax)
        if self.use_weights:
            weights = np.asarray(f['EVENTS'].data.field(self.we_col_name))
            mask = mask & (weights > self.wmin)
            self.weights = weights[mask]
        else: self.weights = None
        mets = np.asarray(f['EVENTS'].data.field('TIME'))[mask]
        self.mjds = mc(mets)
        self.ph = self.polyco.vec_evalphase(self.mjds)
        f.close()

        print >>sys.stderr, "Cuts left %d out of %d events." % (mask.sum(), len(mask))

    def write_phase(self,col_name='PULSE_PHASE'):
        f = pyfits.open(self.ft1file)
        mc = METConverter(self.ft1file)
        mjds = mc(np.asarray(f['EVENTS'].data.field('TIME')))
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
        if self.weights is None: return self.ph[mask],None
        return self.ph[mask],self.weights[mask]
