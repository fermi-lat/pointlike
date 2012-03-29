"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/timeman.py,v 1.2 2011/07/08 23:05:48 kerrm Exp $

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

class GeoConverter(object):
    def __init__(self,ft2,ra,dec):
        self.ft2 = ft2
        self.ra = ra
        self.dec = dec
    def __call__(self,times):
        if not self.can_geo():
            raise Exception('Cannot geocenter!  Must provide FT2 and position.')
        else:
            print 'Attempting to geocenter on-the-fly.'
        from skymaps import PythonUtilities
        times = times.astype(np.float64)
        PythonUtilities.met2geo(times,self.ra,self.dec,self.ft2)
        return times
    def can_geo(self):
        return (self.ft2 is not None) and \
               (self.ra is not None) and \
               (self.dec is not None)

class METConverter(object):
    """Convert LAT arrival times (in Mission Elapsed Time) to MJD(UTC)."""

    def __init__(self,ft1file,ft2=None,ra=None,dec=None):

        self.clockcorr = ClockCorr()

        # Read FT1 file
        hdulist = pyfits.open(ft1file)
        ft1hdr  = hdulist['EVENTS'].header
        ft1dat  = hdulist['EVENTS'].data

        if ft1hdr['TIMEREF'] != 'GEOCENTRIC':
            self.geocon = GeoConverter(ft2,ra,dec)
            if not self.geocon.can_geo():
                print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print "# !!!!!!!!! WARNING !!!!!!!!!! TIMEREF is not GEOCENTRIC! This code is intended for GEOCENTERED times!"
                print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            else:
                print 'TIMEREF is not GEOCENTRIC! But will attempt to correct times on-the-fly.'
        else:
            self.geocon = lambda x: x

        if ft1hdr['TIMESYS'] != 'TT':
              print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              print "# !!!!!!!!! WARNING !!!!!!!!!! TIMESYS is not TT.  We are expecting TT times!"
              print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"

       # Collect TIMEZERO and MJDREF
        try:
            self.TIMEZERO = ft1hdr['TIMEZERO']
        except KeyError:
            self.TIMEZERO = ft1hdr['TIMEZERI'] + ft1hdr['TIMEZERF']
        #print"# TIMEZERO = ",TIMEZERO
        try:
            self.MJDREF = ft1hdr['MJDREF']
        except KeyError:
            # Here I have to work around an issue where the MJDREFF key is stored
            # as a string in the header and uses the "1.234D-5" syntax for floats, which
            # is not supported by Python
            self.MJDREF = ft1hdr['MJDREFI'] + float(ft1hdr['MJDREFF'].replace('D','E'))

        TSTART = float(ft1hdr['TSTART'])
        TSTOP = float(ft1hdr['TSTOP'])

        # Compute MJDSTART and MJDSTOP in MJD(UTC)
        self.MJDSTART = self.clockcorr.tt2utc(TSTART/86400.0 + self.MJDREF + self.TIMEZERO)
        self.MJDSTOP  = self.clockcorr.tt2utc(TSTOP/86400.0 + self.MJDREF + self.TIMEZERO)

        hdulist.close()

    def __call__(self,times):
        times = np.asarray([times] if not hasattr(times,'__iter__') else times)
        times = self.geocon(times)
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
    
