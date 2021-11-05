"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/timeman.py,v 1.9 2013/06/26 19:40:08 kerrm Exp $

Handle MET(TT) to MJD(UTC) conversions.

Author: Paul S. Ray <paul.ray@nrl.navy.mil>
"""

import os
import numpy as np
import astropy.io.fits as pyfits
import warnings

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
        #print (utc, idx, self.mjds[idx], self.dt[idx])
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
        self._setup()

    def _setup(self):
        from skymaps import PythonUtilities
        self.label = 'geocenter'
        self.func = PythonUtilities.met2geo

    def __call__(self,times):
        if not self.can_process():
            raise Exception('Cannot %s!  Must provide FT2 and position.'%self.label)
        else:
            #print ('Attempting to %s on-the-fly.'%self.label)
            pass
        if not hasattr(times,'__len__'):
            times = np.array([times])
        elif not hasattr(times,'dtype'):
            times = np.asarray(times)
        from skymaps import PythonUtilities
        times = times.astype(np.float64)
        self.func(times,self.ra,self.dec,self.ft2)
        return times

    def can_process(self):
        return (self.ft2 is not None) and \
               (self.ra is not None) and \
               (self.dec is not None)

class BaryConverter(GeoConverter):

    def _setup(self):
        from skymaps import PythonUtilities
        self.label = 'barycenter'
        self.func = PythonUtilities.met2tdb

class IdentityConverter(GeoConverter):
    def __init__(self,*args):
        pass
    def __call__(self,times):
        return times
    def can_process(self):
        return True

class METConverter(object):
    """ Convert LAT arrival times (in Mission Elapsed Time) to MJD(UTC)
        in a geocentric reference frame.  The user may alternatively
        specify to use a barycentric frame."""

    def __init__(self,ft1file,ft2=None,ra=None,dec=None,bary=False,noprocess=False):

        self.bary = bary
        self.noprocess = noprocess
        self.clockcorr = ClockCorr()

        # Read FT1 file
        hdulist = pyfits.open(ft1file)
        ft1hdr  = hdulist['EVENTS'].header
        ft1dat  = hdulist['EVENTS'].data

        if bary:
            frame = 'SOLARSYSTEM'
            timesys = 'TDB'
            timecon = BaryConverter
        else:
            frame = 'GEOCENTRIC'
            timesys = 'TT'
            timecon = GeoConverter

        if (not noprocess) and (ft1hdr['TIMEREF'] != frame):
            self.timecon = timecon(ft2,ra,dec)
            if not self.timecon.can_process():
                s = 'BARYCENTERED' if bary else 'GEOCENTERED'
                print ("# !!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print ("# !!!!!!!!! WARNING !!!!!!!!!! TIMEREF is not %s! This code is intended for %s times!"%(frame,s))
                print ("# !!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else:
                print ('TIMEREF is not %s! But will attempt to correct times on-the-fly.'%frame)
        else:
            self.timecon = IdentityConverter()

        if (ft1hdr['TIMESYS'] != timesys) and (not noprocess):
              print ("# !!!!!!!!!!!!!!!!!!!!!!!!!!!!")
              print ("# !!!!!!!!! WARNING !!!!!!!!!! TIMESYS is not %s.  We are expecting %s times!"%(timesys,timesys))
              print ("# !!!!!!!!!!!!!!!!!!!!!!!!!!!!")

       # Collect TIMEZERO and MJDREF
        try:
            self.TIMEZERO = ft1hdr['TIMEZERO']
        except KeyError:
            self.TIMEZERO = ft1hdr['TIMEZERI'] + ft1hdr['TIMEZERF']
        #print("# TIMEZERO = ",TIMEZERO)
        try:
            self.MJDREF = ft1hdr['MJDREF']
        except KeyError:
            # some FT1 files store these values as floats; others as
            # strings with the 'D' exponential syntax
            try:
                self.MJDREF = ft1hdr['MJDREFI'] + float(ft1hdr['MJDREFF'])
            except ValueError:
                self.MJDREF = ft1hdr['MJDREFI'] + float(ft1hdr['MJDREFF'].replace('D','E'))

        # timecon added Mar 28 2013; correction needs to be applied in
        # on the fly case; minor for GEO but important for BARY
        TSTART = self.timecon(float(ft1hdr['TSTART']))
        TSTOP = self.timecon(float(ft1hdr['TSTOP']))

        # Compute MJDSTART and MJDSTOP in MJD(UTC); these are used by the
        # binning to compute uniform intervals
        if bary:
            #raise NotImplementedError('Cannot yet convert this to bary.')
            warnings.warn('Cannot yet convert this to bary.')
        self.MJDSTART = self.clockcorr.tt2utc(TSTART/SECSPERDAY + self.MJDREF + self.TIMEZERO)
        self.MJDSTOP  = self.clockcorr.tt2utc(TSTOP/SECSPERDAY + self.MJDREF + self.TIMEZERO)

        hdulist.close()

    def __call__(self,times):
        """ Convert MET to MJD(UTC)/GEO."""
        times = np.asarray([times] if not hasattr(times,'__iter__') else times)
        times = self.timecon(times)
        times = times/SECSPERDAY + self.MJDREF + self.TIMEZERO
        # disable clock corrections for SSB times
        if (not self.bary) and (not self.noprocess):
            for i in range(len(times)):
                times[i] = self.clockcorr.tt2utc(times[i])
        if len(times) == 1: return times[0]
        return times

# standalone function for converting non-FT1 based METs
def met2mjd_utc(times,mjdref=51910+7.428703703703703e-4,tzero=0):
    clockcorr = ClockCorr()
    times = np.asarray([times] if not hasattr(times,'__iter__') else times)
    times = times/SECSPERDAY + mjdref + tzero # copy
    for i in range(len(times)):
        times[i] = clockcorr.tt2utc(times[i])
    if len(times) == 1: return times[0]
    return times
    
