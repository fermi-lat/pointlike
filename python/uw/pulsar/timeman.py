"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/timeman.py,v 1.1 2011/04/27 18:32:03 kerrm Exp $

Handle MET(TT) to MJD(UTC) conversions.

Author: Paul S. Ray <paul.ray@nrl.navy.mil>
"""

import os
import numpy as np
import pyfits

SECSPERDAY = 86400.0

class ClockCorr(object):
    '''Converts MJD among UTC, TT, TAI using tempo2 clock correction file.'''
    def __init__(self,fname=None):
        if fname is None:
            fname = os.environ['TEMPO2'] + '/clock/utc2tai.clk'
        mjd,dt = np.loadtxt(fname,unpack=True)
        self.mjds = mjd
        self.dt = dt
    def getcorr(self,utc):
        '''Return value of TAI-UTC at a given MJD(UTC)'''
        idx = np.where(self.mjds<utc)[0][-1]
        #print utc, idx, self.mjds[idx], self.dt[idx]
        corr = self.dt[idx]
        return(corr)
    def tt2tai(self,tt):
        return(tt - 32.184/SECSPERDAY)
    def tai2tt(self,tai):
        return(tai + 32.184/SECSPERDAY)
    def tt2utc(self,tt):
        # BUG!! getcorr MUST be called with UTC, not tai!
        # This will cause a problem if you are within ~1 min of a leapsec
        tai = self.tt2tai(tt)
        corr = self.getcorr(tai)
        utc = tai-corr/SECSPERDAY
        return(utc)
    def utc2tai(self,utc):
        corr = self.getcorr(utc)
        tai = utc + corr/SECSPERDAY
        return(tai)


class METConverter(object):
    """Convert LAT arrival times (in Mission Elapsed Time) to MJD(UTC)."""

    def __init__(self,ft1file):

        clockcorr = ClockCorr()

        # Read FT1 file
        hdulist = pyfits.open(ft1file)
        ft1hdr  = hdulist['EVENTS'].header
        ft1dat  = hdulist['EVENTS'].data

        if ft1hdr['TIMEREF'] != 'GEOCENTRIC':
              print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              print "# !!!!!!!!! WARNING !!!!!!!!!! TIMEREF is not GEOCENTRIC! This code is intended for GEOCENTERED times!"
              print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"

        if ft1hdr['TIMESYS'] != 'TT':
              print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              print "# !!!!!!!!! WARNING !!!!!!!!!! TIMESYS is not TT.  We are expecting TT times!"
              print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"

       # Collect TIMEZERO and MJDREF
        try:
            TIMEZERO = ft1hdr['TIMEZERO']
        except KeyError:
            TIMEZERO = ft1hdr['TIMEZERI'] + ft1hdr['TIMEZERF']
        #print"# TIMEZERO = ",TIMEZERO
        try:
            MJDREF = ft1hdr['MJDREF']
        except KeyError:
            # Here I have to work around an issue where the MJDREFF key is stored
            # as a string in the header and uses the "1.234D-5" syntax for floats, which
            # is not supported by Python
            MJDREF = ft1hdr['MJDREFI'] + \
            float(ft1hdr['MJDREFF'].replace('D','E'))

        self.clockcorr = clockcorr
        self.MJDREF    = MJDREF
        self.TIMEZERO  = TIMEZERO

        TSTART = float(ft1hdr['TSTART'])
        TSTOP = float(ft1hdr['TSTOP'])

        # Compute MJDSTART and MJDSTOP in MJD(UTC)
        self.MJDSTART = clockcorr.tt2utc(TSTART/86400.0 + MJDREF + TIMEZERO)
        self.MJDSTOP  = clockcorr.tt2utc(TSTOP/86400.0 + MJDREF + TIMEZERO)

        hdulist.close()

    def __call__(self,times):
        times = np.asarray([times] if not hasattr(times,'__iter__') else times)
        times = times/SECSPERDAY + self.MJDREF + self.TIMEZERO
        for i in xrange(len(times)):
            times[i] = self.clockcorr.tt2utc(times[i])
        if len(times) == 1: return times[0]
        return times

# standalone function for converting non-FT1 based METs
def met2mjd_utc(times,mjdref=51910+7.428703703703703e-4,tzero=0):
    clockcorr = ClockCorr()
    times = np.asarray([times] if not hasattr(times,'__iter__') else times)
    times = times/SECSPERDAY + mjdref + tzero # copy
    for i in xrange(len(times)):
        times[i] = clockcorr.tt2utc(times[i])
    if len(times) == 1: return times[0]
    return times
    
