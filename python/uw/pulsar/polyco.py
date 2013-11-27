"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/polyco.py,v 1.14 2013/06/30 01:24:53 kerrm Exp $

Mange polycos from tempo2.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""
from __future__ import division
import math
import numpy as np
import datetime
import os
import subprocess

class PolycoEntry:
    STATIC_COUNTER = 0
    def __init__(self,tmid,mjdspan,rphase,f0,ncoeff,coeffs,obs):
        self.tmid = tmid
        self.mjdspan = mjdspan
        self.tstart = tmid - float(mjdspan)/2
        self.tstop = tmid + float(mjdspan)/2
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

    def evalabsphase(self,t):
        '''Return the phase at time t, computed with this polyco entry'''
        dt = (t-self.tmid)*1440.0
        # Compute polynomial by factoring out the dt's
        phase = self.coeffs[self.ncoeff-1]
        for i in range(self.ncoeff-2,-1,-1):
            phase = self.coeffs[i] + dt*phase
        # Add DC term
        phase += self.rphase + dt*60.0*self.f0
        return(phase)

    def evalfreq(self,t):
        '''Return the freq at time t, computed with this polyco entry'''
        dt = (t-self.tmid)*1440.0
        s = 0.0
        for i in range(1,self.ncoeff):
            s += float(i) * self.coeffs[i] * dt**(i-1)
        freq = self.f0 + s/60.0
        return(freq)

    def evalfreqderiv(self,t):
        """ Return the frequency derivative at time t."""
        dt = (t-self.tmid)*1440.0
        s = 0.0
        for i in range(2,self.ncoeff):
            s += float(i) * float(i-1) * self.coeffs[i] * dt**(i-2)
        freqd = s/(60.0*60.0)
        return(freqd)

class Polyco:
    def __init__(self, fname, psrname=None, recalc_polycos=True,
        mjd0=51544,bary=False, working_dir=None, output=None, ndays=None,
        verbose=False):
        """ Create an object encapsulating a set of polynomial coefficients             for evaluating phase.

            fname -- either an existing polyco .dat file or an ephemeris
                     with which to generate the polycos
            recalc_polycos -- force generation of polycos; fname must be an
                     ephemeris in this case
            mjd0 -- start of polyco validity; default Jan 1, 2000
            bary -- generate polycos at barycenter if set (default geo)
            working_dir -- change to this directory to generate polycos
            output -- use this stem to prepend to polyco .dat files
            ndays -- number of days to include; default spans 2000 to 
                     present but can be lengthy to compute
            verbose -- spit out more diagnostic output
        """

        self.bary = bary
        self.working_dir = working_dir
        self.output = output
        self.verbose = verbose
        if fname.endswith( ".par" ) or recalc_polycos:
            from uw.pulsar.parfiles import ParFile
            pf = ParFile(fname)
            self.ra = pf.get_ra()
            self.dec = pf.get_dec()
            fname = self.gen_polycos(fname,recalc_polycos=recalc_polycos,mjd0=mjd0,ndays=ndays)
        else:
            self.ra = self.dec = None
        
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
            if verbose:
                print "- - - - - - -"
                #print "psrname %s date %s utc %s tmid %s dm %f doppler %f logrms %f" % (psrname,date,utc,tmid,dm,doppler,logrms)
                print "psrname %s date %s utc %s tmid %s dm %f logrms %f" % (psrname,date,utc,tmid,dm,logrms)
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
            if verbose:
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
            if verbose:
                print "COEFFS: ",coeffs
            pe = PolycoEntry(tmid,mjdspan,rphase,f0,ncoeff,coeffs,obs)
            self.entries.append(pe)
        if len(self.entries)==0:
            raise ValueError('No polycos generated!')
        self.make_keys()

    def gen_polycos(self,polyconame,recalc_polycos=True,mjd0=51544,ndays=None):
        """If par file passed in, generate polyco file on the fly."""

        # get MJDs
        if ndays is None:
            nDays=(datetime.date.today()-datetime.date(2000,1,1)).days+(51544-mjd0)
        else:
            nDays = ndays
        endMJD=mjd0+nDays+2
        if (endMJD-mjd0) < 2:
            raise ValueError('Unacceptable MJD bounds.')
        if self.verbose:
            print "MJD limits: %s %s"%(str(mjd0),str(endMJD))
        curdir = os.getcwd()
        if self.working_dir is not None:
            os.chdir(self.working_dir)
        prefix = self.output or ''
        if recalc_polycos:
            fnames = ['%s%s'%(prefix,x) for x in 
                ['polyco.tim','polyco_new.dat','newpolyco.dat']]
            map(os.remove,filter(os.path.isfile,fnames))
            obs_string = '@' if self.bary else 'coe'
            out_string = '' if self.output is None else ' -polyco_file %s'%self.output
            t2cmd = 'tempo2 -f %s%s -polyco "%s %s 360 12 12 %s 0 0\"'%(
                polyconame,out_string,mjd0,endMJD,obs_string)
            o = subprocess.check_output(t2cmd,shell=True)
            if self.verbose:
                print 'Creating polycos with command:\n',t2cmd
                print o
        fname = '%spolyco_new.dat'%(prefix)
        polyconame=os.path.abspath(fname)
        DEVNULL = open(os.devnull,'wb')
        subprocess.call('rm %snewpolyco.dat polyco.tim'%(prefix),
            shell=True,stderr=DEVNULL)
        os.chdir(curdir)
        return polyconame

    def make_keys(self):
        """Keys for a binary search.  Use the edges."""
        keys = np.asarray([e.tstop for e in self.entries])
        sorting = np.argsort(keys)
        self.entries = np.asarray(self.entries)[sorting]
        self.keys = np.append(self.entries[0].tstart,keys[sorting])

    def getentry(self,t,use_keys=True):
        '''Returns the polyco entry corresponding to time t (in MJD)'''
        if use_keys:
            idx = np.searchsorted(self.keys,t)
            if np.any(idx == len(self.keys)) or np.any(idx==0):
                print 'Could not find a valid entry for MJD(s)...'
                print t[idx == len(self.keys)] if type(t) is type(np.array([1])) else t
                print t[idx == 0] if type(t) is type(np.array([1])) else t
                raise IndexError
            return self.entries[idx-1]
        for pe in self.entries:
            if pe.valid(t):
                return pe
        # This should probably throw an exception instead.
        print "No valid polyco found for MJD ",t
        sys.exit(9)
        return None

    def _vec_eval(self,times,func):
        if not hasattr(times,'__len__'):
            times = [times]
        pces = self.getentry(times,use_keys=True)
        self.ids = np.asarray([pc.uid for pc in pces])
        return np.asarray(map(func,pces,times))

    def vec_evalphase(self,times):
        """ Return the phases for a vector of times; NB times should be in
            MJD @ GEO."""
        return self._vec_eval(times,PolycoEntry.evalphase)

    def vec_evalabsphase(self,times):
        """ Return the phases for a vector of times; NB times should be in
            MJD @ GEO."""
        return self._vec_eval(times,PolycoEntry.evalabsphase)

    def vec_evalfreq(self,times):
        """ Return the phases for a vector of times; NB times should be in
            MJD @ GEO."""
        return self._vec_eval(times,PolycoEntry.evalfreq)

    def vec_evalfreqderiv(self,times):
        """ Return the phases for a vector of times; NB times should be in
            MJD @ GEO."""
        return self._vec_eval(times,PolycoEntry.evalfreqderiv)

    def invert_phase_shift(self,t0,phi):
        """ Compute the time lapse (in s) corresponding to phi at t0."""
        pe = self.getentry(t0)
        f = pe.evalfreq(t0)
        return phi/f

    def get_obs(self):
        return self.entries[0].obs
            
