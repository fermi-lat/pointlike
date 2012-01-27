"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/polyco.py,v 1.2 2011/07/07 19:10:38 kerrm Exp $

Mange polycos from tempo2.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""
from __future__ import division
import math
import numpy as np
import datetime
import os

class PolycoEntry:
    STATIC_COUNTER = 0
    def __init__(self,tmid,mjdspan,rphase,f0,ncoeff,coeffs,obs):
        self.tmid = tmid
        self.mjdspan = mjdspan
        self.rphase = rphase
        self.f0 = f0
        self.ncoeff = ncoeff
        self.coeffs = coeffs
        self.obs = obs
        self.uid = PolycoEntry.STATIC_COUNTER
        PolycoEntry.STATIC_COUNTER += 1
    def __str__(self):
        return("PE: "+repr(self.tmid)+" "+repr(self.mjdspan)+" "+repr(self.rphase)+" "+repr(self.ncoeff)+" "+repr(self.coeffs))
        
    def valid(self,t):
        '''Return True if this polyco entry is valid for the time given (MJD)'''
        return t>=(self.tmid-self.mjdspan/2.0) and t<(self.tmid+self.mjdspan/2.0)

    def evalphase(self,t):
        '''Return the phase at time t, computed with this polyco entry'''
        dt = (t-self.tmid)*1440.0
        # Compute polynomial by factoring out the dt's
        phase = self.coeffs[self.ncoeff-1]
        for i in range(self.ncoeff-2,-1,-1):
            phase = self.coeffs[i] + dt*phase
        # Add DC term
        phase += self.rphase + dt*60.0*self.f0
        phase -= math.floor(phase)
        if phase < 0.0 or phase >= 1.0:
            print "BAD PHASE ",phase
        return(phase)

    def evalfreq(self,t):
        '''Return the freq at time t, computed with this polyco entry'''
        dt = (t-self.tmid)*1440.0
        sum = 0.0
        for i in range(1,self.ncoeff):
            sum += float(i) * self.coeffs[i] * dt**(i-1)
        freq = self.f0 + sum/60.0
        return(freq)

class Polyco:
    def __init__(self, fname, psrname=None, recalc_polycos=True,mjd0=51544):

        if fname.endswith( ".par" ) or recalc_polycos:
            from uw.pulsar.parfiles import ParFile
            pf = ParFile(fname)
            self.ra = pf.get_ra()
            self.dec = pf.get_dec()
            fname = self.gen_polycos(fname,recalc_polycos=recalc_polycos,mjd0=mjd0)
        else:
            self.ra = self.dec = None
        
        VERBOSE= False
        self.entries = []
        f = open(fname,"r")
        set = 0
        while True:
            line1 = f.readline()
            if len(line1) == 0:
                break
            sp = line1.split()
            psrname = sp[0].strip()
            date = sp[1].strip()
            utc = sp[2]
            tmid = float(sp[3])
            dm = float(sp[4])
            #doppler = float(sp[5])
            logrms = float(sp[6])
            if VERBOSE:
                print "- - - - - - -"
                print "psrname %s date %s utc %s tmid %s dm %f doppler %f logrms %f" % (psrname,date,utc,tmid,dm,doppler,logrms)
            line2 = f.readline()
            rphase = float(line2[0:20])
            f0 = float(line2[20:38])
            obs = line2[38:43].strip()
            nspan = int(line2[43:49])
            mjdspan = float(nspan)/(60*24)
            ncoeff = int(line2[49:54])
            obsfreq = float(line2[54:64])
            if len(line2[75:80].strip()) > 0:
                binphase = float(line2[75:80])
            else:
                binphase = 0.0
            if VERBOSE:
                print "rphase %s f0 %s obs %s ncoeff %d nspan %d obsfreq %f binphase %f" % (repr(rphase),repr(f0),obs,ncoeff,nspan,obsfreq,binphase)
            nlines = ncoeff//3
            nlast = ncoeff%3
            if nlast > 0:
                nlines += 1
            coeffs = []
            for i in range(nlines):
                line = f.readline()
                for c in line.split():
                    coeffs.append(float(c))
            coeffs = np.array(coeffs)
            if VERBOSE:
                print "COEFFS: ",coeffs
            pe = PolycoEntry(tmid,mjdspan,rphase,f0,ncoeff,coeffs,obs)
            self.entries.append(pe)
        self.make_keys()

    def gen_polycos(self,polyconame,recalc_polycos=True,mjd0=51544):
        """If par file passed in, generate polyco file on the fly."""

        # get MJDs
        nDays=(datetime.date.today()-datetime.date(2000,1,1)).days+(51544-mjd0)
        endMJD=mjd0+nDays+2
        print "MJD limits: %s %s"%(str(mjd0),str(endMJD))
        if recalc_polycos:
            os.system( "rm polyco_new.dat newpolyco.dat polyco.tim" )
            os.system( "tempo2 -f " + polyconame + " -polyco \"54628 "+ str(endMJD) + " 360 12 12 coe 0 0\"" )
        polyconame="polyco_new.dat"
        return polyconame

    def make_keys(self):
        """Keys for a binary search."""
        keys = [e.tmid - e.mjdspan/2. for e in self.entries]
        sorting = np.argsort(keys)
        self.entries = np.asarray(self.entries)[sorting]
        keys = np.asarray(keys)[sorting]
        self.keys = np.append(keys,[self.entries[-1].tmid + self.entries[-1].mjdspan/2.])

    def getentry(self,t,use_keys=True):
        '''Returns the polyco entry corresponding to time t (in MJD)'''
        if use_keys:
            idx = np.searchsorted(self.keys,t)
            if np.any(idx == len(self.keys)):
                print 'Could not find a valid entry for MJD(s)...'
                print t[idx == len(self.keys)] if type(t) is type(np.array([1])) else t
                raise ValueError
            return self.entries[idx-1]
        for pe in self.entries:
            if pe.valid(t):
                return pe
        # This should probably throw an exception instead.
        print "No valid polyco found for MJD ",t
        sys.exit(9)
        return None

    def vec_evalphase(self,times):
        """ Return the phases for a vector of times; NB times should be in
            MJD @ GEO."""
        pces = self.getentry(times,use_keys=True)
        self.ids = np.asarray([pc.uid for pc in pces])
        return np.asarray([pce.evalphase(time) for pce,time in zip(pces,times)])
