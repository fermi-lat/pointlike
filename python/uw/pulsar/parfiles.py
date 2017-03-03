"""
Module reads and manipulates tempo2 parameter files.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/parfiles.py,v 1.76 2017/01/17 15:23:22 kerrm Exp $

author: Matthew Kerr
"""

import numpy as np
import os
import subprocess
from uw.utilities.coords import ec2eq,eq2ec
from collections import deque,defaultdict
import tempfile
import time

C = 29979245800.

def sex2dec(s,mode='ra'):
    """ Convert a colon- or space-delimited sexagesimal string to decimal 
        degrees.
    """
    if ':' in s:
        toks = s.split(':')
    else:
        toks = s.split(' ')
    multi = (-1 if '-' in s else 1) * (15 if mode == 'ra' else 1)
    return multi*sum( (60**(-i) * abs(float(toks[i])) for i in xrange(len(toks)) ) )

def ra2dec(s): return sex2dec(s,mode='ra')
def decl2dec(s): return sex2dec(s,mode='decl')
def dec2sex(d,mode='ra',places=6):
    """ Convert decimal degree to column separated format.
    
    Strip out each entry (deg/min/sec) separately, then take fractional
    seconds as a fixed precision string (specified by places kwarg).
    """
    d = np.float128(d)
    orig = np.float128(d)
    scale = 15. if mode=='ra' else 1. # scale to deg to hour (ra) or deg
    sign = '' if mode=='ra' else ('-' if orig < 0 else '+')
    d = np.float128(abs(d))
    dg = int(d/scale)
    mi = int((d/scale-dg)*60.)
    ss = ((d/scale-dg)*60-mi)*60
    ssi = int(ss)
    ssf = int(round(10**places*(ss-ssi)))
    if ssf==10**places:
        ssf = 0 # prevent overflow at this precision
        ssi += 1
    ssfs = '0'*(places-len(str(ssf)))+str(ssf)
    s = '%s%02d:%02d:%02d.%s'%(sign,dg,mi,ssi,ssfs)
    # check the routine is inverse to sex2dec to indicated precision
    precision = 2.0*scale/3600.*10**-places # prec. in deg. with 2x fudge
    diff = abs(sex2dec(s,mode=mode)-orig)
    assert(diff<precision)
    return s

def pad(s,n,c=' '):
    if len(s) >= n:
        return s + c # always have a trailing space!
    return s + c*(n-len(s))
def pad26(s): 
    """ Specialized version to handle single bit (free/fixed) flags."""
    if len(s)==1:
        return s+'  '
    return pad(s,26)

class StringFloat(object):
    """ Use strings and python longs to handle float strings with arbitrary
        precision."""

    def __init__(self,s):
        s = self._s = self._parse_enotation(str(s))
        self._i = int(s.replace('.',''))
        # "places" keeps track of places right of decimal pt
        if '.' in s:
            self.places = len(s.split('.')[-1])
        else:
            self.places = 0

    def _parse_enotation(self,s):
        s = s.replace('E','e').replace('d','e').replace('D','e')
        if 'e' in s:
            l,r = s.split('e')
            if '.' not in l:
                l = l + '.'
            sl,sr = l.split('.')
            if r[0] == '-':
                places = int(r[1:])
                s = '0.' + '0'*(places-len(sl)) + sl + sr
            else:
                if r[0] == '+': r = r[1:]
                places = int(r)
                s = sl + sr + '0'*(places-len(sr))
        return s

    def copy(self):
        return StringFloat(self._s)

    def __str__(self):
        return self._s

    def __float__(self):
        return float(self._s)

    def __add__(self,other):
        if not hasattr(other,'places'):
            other = StringFloat(str(other))
        if self.places==0 and other.places==0:
            return StringFloat(str(self._i+other._i))
        t = max(self.places,other.places)
        i1 = int(str(self._i) + '0'*(t-self.places))
        i2 = int(str(other._i) + '0'*(t-other.places))
        s = i1+i2
        sign = '-' if s < 0 else ''
        s = str(abs(s))
        return StringFloat(sign+s[:-t]+'.'+s[-t:])

    def __mul__(self,other):
        if not hasattr(other,'places'):
            other = StringFloat(str(other))
        if self.places==0 and other.places==0:
            return StringFloat(str(self._i*other._i))
        t = (self.places+other.places)
        i = self._i*other._i
        sign = '-' if i < 0 else ''
        s = str(abs(i))
        if t > len(s):
            return StringFloat(sign+'0.'+'0'*(t-len(s))+s)
        return StringFloat(sign+s[:-t]+'.'+s[-t:])

    def __sub__(self,other):
        if not hasattr(other,'places'):
            other = StringFloat(str(other))
        return self + other*StringFloat(-1)

    def sigfig_string(self,sigfigs):
        """ Return a string version with specified number of sig figs."""
        sign = self._i < 0
        icop = abs(self._i)
        mysigfigs = len(str(icop))
        if sigfigs >= mysigfigs:
            return self._s
        i = str(icop)[:(sigfigs+1)]
        if int(i[-1]) < 5:
            i = int(i[:-1])
        else:
            i = int(i[:-1])+1
        places = self.places
        # now figure out where decimal goes
        i = str(i)
        if places >= mysigfigs:
            i = '0.' + '0'*(places-mysigfigs) + i
        elif places < mysigfigs:
            dp = mysigfigs-places # decimal point goes in this position
            if dp < sigfigs:
                # dp goes in middle of #
                i = i[:dp] + '.' + i[dp:]
            else:
                i = i + '0'*(dp-sigfigs) + '.'
        i =  i.rstrip('.')
        if sign:
            return '-'+i
        return i


class ParFile(dict):

    def __init__(self,parfile):
        if not os.path.exists(parfile):
            raise IOError('Indicated file %s does not exist.'%parfile)
        self.parfile = parfile
        self.init()

    def init(self,buf=None):
        # make this a little more bulletproof on the cluster
        self.ordered_keys = []
        self.duplicates = defaultdict(list)
        comment_counter = 0
        for i in xrange(3):
            if buf is None:
                f = open(self.parfile,'r')
                lines = f.readlines()
                f.close()
            else:
                buf.seek(0)
                lines = buf.readlines()
            if len(lines) > 0:
                break
            else:
                time.sleep(1)
        if len(lines) == 0:
            raise IOError('Could not read data from %s.'%self.parfile)
        for line in lines:
            tok = line.strip().split()
            if len(tok)==0: continue
            if line.strip()[0] == '#': # handle comments
                tok = ['#COMMENT%d'%(comment_counter),line.strip()[1:]]
                comment_counter += 1
            key = tok[0]
            val = tok[1:] if (len(tok[1:]) > 1) else tok[1:][0]
            if key not in self.ordered_keys:
                self.ordered_keys.append(tok[0])               
                self[key] = val
            else:
                self.duplicates[key].append(val)
        self.degree = self._calc_degree()
    
    def _calc_degree(self):
        degree = -1
        for i in xrange(12):
            if 'F%d'%(i) in self.keys(): degree += 1
            else: break
        return degree

    def get(self,key,first_elem=True,type=str):
        if type==float:
            type = lambda x: float(x.replace('D','E'))
        t = self[key]
        if hasattr(t,'__iter__'):
            if first_elem: return type(t[0])
            return [type(x) for x in t]
        return type(t)

    def set(self,key,val):
        """ Replace the FIRST ELEMENT of a field with val."""
        t = self[key];
        if not type(val)==type(''):
            val = '%.18g'%(val)
        if hasattr(t,'__iter__'):
            self[key][0] = val
        else: self[key] = val

    def get_psrname(self,add_j=True):
        try:
            name = self.get('PSR')
        except KeyError:
            name = self.get('PSRJ')
        if (name[0] != 'J') and add_j:
            name = 'J' + name
        return name

    def get_ra(self):
        try:
            return ra2dec(self.get('RAJ'))
        except KeyError:
            elong = self.get("ELONG",type=float)
            elat = self.get("ELAT",type=float)
            ra, dec = ec2eq(elong,elat)
            return ra[0]
        
    def get_dec(self):
        try:
            return decl2dec(self.get('DECJ'))
        except KeyError:
            elong, elat = self.get("ELONG",type=float), self.get("ELAT",type=float)
            ra, dec = ec2eq(elong,elat)
            return dec[0]                                            
        
    def get_skydir(self):
        from skymaps import SkyDir
        ra,dec = self.get_ra(), self.get_dec()
        return SkyDir(ra,dec)

    def get_binary_period(self):
        """ Return binary period in days."""
        try:
            return self.get('PB',type=float)
        except KeyError:
            return None

    def is_ecliptic(self):
        try:
            elon = self['ELONG']
            return True
        except Exception:
            pass
        return False

    def get_astrometry(self,epoch=None):
        """ Return the astrometry as a 4-tuple with RA, error, Dec, error
            in degrees.
            
            If the epoch is provided, and the source has a proper motion,
            then correct the position to given epoch and add the errors
            in quadrature.
            
            NB that for tempo2, both PMRA and PMDEC are physical scales,
            so PMRA shoud be divided by cos(DEC) to get the "coordinate
            speed".
        """
        ra = ra2dec(self.get('RAJ'))
        dec = decl2dec(self.get('DECJ'))
        if hasattr(self['RAJ'],'__iter__'):
            idx = 1 if (len(self['RAJ'])==2) else 2
            rae = float(self['RAJ'][idx]) * (15./3600) # deg
        else:
            rae = 0.
        if hasattr(self['DECJ'],'__iter__'):
            idx = 1 if (len(self['DECJ'])==2) else 2
            dece = float(self['DECJ'][idx]) * (1./3600) # deg
        else:
            dece = 0.
        
        if epoch is not None:
            try:
                dt = (epoch - self.get('POSEPOCH',type=float))/365.24
            except KeyError:
                dt = 0
            if 'PMRA' in self.keys():
                #in TEMPO2,  proper motion is in physical units; convert
                # to coordinate units by correcting for cos(dec)
                cdec = np.cos(np.radians(dec))
                pmra = self.get('PMRA',type=float)/1000/3600./cdec # deg
                if hasattr(self['PMRA'],'__iter__'):
                    idx = 1 if (len(self['PMRA'])==2) else 2
                    pmrae = float(self['PMRA'][idx])/1000/3600./cdec # deg
                else:
                    pmrae = 0
                ra += pmra*dt
                rae = (rae**2 + (pmrae*dt)**2)**0.5
            if 'PMDEC' in self.keys():
                pmdec = self.get('PMDEC',type=float)/1000/3600 # deg
                if hasattr(self['PMDEC'],'__iter__'):
                    idx = 1 if (len(self['PMDEC'])==2) else 2
                    pmdece = float(self['PMDEC'][idx])/1000/3600 # deg
                else:
                    pmdece = 0
                dec += pmdec * dt
                dece = (dece**2 + (pmdece*dt)**2)**0.5
        return [ra,rae,dec,dece]

    def get_astrometry_ecliptic(self,epoch=None):
        """ Return the astrometry as a 4-tuple with ELONG, error, ELAT, 
            error, in degrees.

            NB in tempo2 convention ELONG and ELAT are given in degr.
            
            If the epoch is provided, and the source has a proper motion,
            then correct the position to given epoch and add the errors
            in quadrature.
            
            NB that for tempo2, both PMELONG and PMELAT are physical
            scales, so PMELONG shoud be divided by cos(ELAT) to get the 
            "coordinate speed".
        """
        elon = float(self.get('ELONG'))
        elat = float(self.get('ELAT'))
        if hasattr(self['ELONG'],'__iter__'):
            idx = 1 if (len(self['ELONG'])==2) else 2
            elone = float(self['ELONG'][idx])
        else:
            elone = 0.
        if hasattr(self['ELAT'],'__iter__'):
            idx = 1 if (len(self['ELAT'])==2) else 2
            elate = float(self['ELAT'][idx])
        else:
            elate = 0.

        if epoch is None:
            return [elon,elone,elat,elate]
        
        try:
            dt = (epoch - self.get('POSEPOCH',type=float))/365.24
        except KeyError:
            dt = 0
        if 'PMELONG' in self.keys():
            #in TEMPO2,  proper motion is in physical units; convert
            # to coordinate units by correcting for cos(ELAT)
            celat = np.cos(np.radians(elat))
            pmelon = self.get('PMELONG',type=float)/1000/3600./celat
            if hasattr(self['PMELONG'],'__iter__'):
                idx = 1 if (len(self['PMELONG'])==2) else 2
                pmelone = float(self['PMELONG'][idx])/1000/3600./celat
            else:
                pmelone = 0
            elon += pmelon*dt
            elone = (elone**2 + (pmelone*dt)**2)**0.5
        if 'PMELAT' in self.keys():
            pmelat = self.get('PMELAT',type=float)/1000/3600 # deg
            if hasattr(self['PMELAT'],'__iter__'):
                idx = 1 if (len(self['PMELAT'])==2) else 2
                pmelate = float(self['PMELAT'][idx])/1000/3600 # deg
            else:
                pmelate = 0
            elat += pmelat * dt
            elate = (elate**2 + (pmelate*dt)**2)**0.5
        return [elon,elone,elat,elate]

    def convert_ecliptic(self):
        """ Convert ecliptic coordinates to celestial coordinates."""
        if not self.is_ecliptic():
            return
        if ('PMELONG' in self.keys()) or ('PMELAT' in self.keys()):
            elon1,a,elat1,b = self.get_astrometry_ecliptic(epoch=55000)
            elon2,a,elat2,b = self.get_astrometry_ecliptic(epoch=55365.24)
            raj1,decj1 = ec2eq(elon1,elat1)
            raj2,decj2 = ec2eq(elon2,elat2)
            pmra = (raj2-raj1)*3600*1000*np.cos(
                np.radians(0.5*(decj1+decj2)))
            pmdec = (decj2-decj1)*3600*1000
            self.delete_key('PMELONG')
            self.delete_key('PMELAT')
            self.add_key('PMRA',pmra[0])
            self.add_key('PMDEC',pmdec[0])
        raj = self.get_ra()
        decj = self.get_dec()
        self.add_key('RAJ',dec2sex(raj,mode='ra'))
        self.add_key('DECJ',dec2sex(decj,mode='dec'))
        self.delete_key('ELONG')
        self.delete_key('ELAT')

    def set_posepoch(self,epoch):
        """ Update POSEPOCH, evolving RAJ/DECJ if necessary."""
        try:
            posepoch = self['POSEPOCH']
        except KeyError:
            self.add_key('POSEPOCH',str(epoch))
            return
        if self.is_ecliptic():
            elon,elone,elat,elate = self.get_astrometry_ecliptic(
                epoch=epoch)
            self.set('ELONG',elon)
            self.set('ELAT',elat)
        else:
            raj,raje,decj,decje = self.get_astrometry(
                epoch=epoch)
            self.set('RAJ',dec2sex(raj,mode='ra'))
            self.set('DECJ',dec2sex(decj,mode='dec'))
        self.set('POSEPOCH',str(epoch))

    def get_icrs_coord(self):
        """ Return an astropy ICRSCoordinate object."""
        from astropy.coordinates import ICRSCoordinates
        ra = self.get('RAJ')
        dec = self.get('DECJ')
        return ICRSCoordinates(ra=ra,dec=dec,unit=('hour','deg'))

    def p(self,error=False):
        f0 = self.get('F0',first_elem=not error,type=float)
        if error and hasattr(f0,'__iter__'):
            f0,err = f0[0],f0[-1]
            return 1./f0,err/f0**2
        return 1./f0
 
    def pdot(self,error=False):
        f0 = self.get('F0',first_elem=not error,type=float)
        f1 = self.get('F1',first_elem=not error,type=float)
        if error and hasattr(f0,'__iter__'):
            f0,f0_err = f0[0],f0[-1]
            f1,f1_err = f1[0],f1[-1]
            pdot = -f1/f0**2
            pdot_err = 1./f0**2*(f1_err**2+4*pdot*f0_err**2)**0.5
            return pdot,pdot_err
        return -f1/f0**2

    def edot(self,mom=1e45,distance=None):
        om = self.get('F0',type=float)*(np.pi*2)
        omdot = self.get('F1',type=float)*(np.pi*2)
        edot = -mom*om*omdot # say this without laughing
        if distance is not None:
            corr = 1-self.get_shklovskii_pdot(distance)/self.pdot()
        else: corr = 1
        return edot*corr

    def is_msp(self,maxp=0.03,maxpdot=1e-18):
        # hard cut on period to catch energetic msps
        if self.p() < 0.018:
            return True
        # otherwise, restrict to lower rectangle of P/Pdot space
        return (self.p() < maxp) and (self.pdot() < maxpdot)

    def get_time_cuts(self):
        import datetime
        from uw.utilities.fermitime import utc_to_met
        tomorrow = datetime.date.fromtimestamp(time.time()+86400)
        tmin = 239557418 # first photon timestamp
        tmax = utc_to_met(tomorrow.year,tomorrow.month,tomorrow.day)
        if self.is_msp():
            return tmin,tmax
        t0 = self.get('START',type=float)
        t1 = self.get('FINISH',type=float)
        mjd_0 = 51910 + 7.428703703703703e-4
        met_0 = (t0 - mjd_0) * 86400
        met_1 = (t1 - mjd_0) * 86400
        return max(met_0,0), max(met_1,0)

    def get_binary_dict(self):
        try:
            binary_dict = {}
            keys = ['A1','PB','T0','OM'] + ['E' if 'E' in self.keys() else 'ECC']
            for key in keys:
                binary_dict[key] = self.get(key,first_elem=True,type=float)
            binary_dict['PB'] *= 86400
            return binary_dict
        except: return None

    def add_binary_model(self,idx=None,clobber=False):
        """ Add default binary parameters to model.
        
        NB -- for the T2 model, the accepted way of having more than
        one binary companion seems to be
        (1) T0/PB etc. for first companion
        (2) T0_1/PB_1 etc. for second companion
        (2) T0_2/PB_2 etc. for third companion
        etc.
        """
        if (not clobber) and ('BINARY' in self.ordered_keys):
            return
        self.add_key('BINARY','T2')
        if (idx is None):
            label = ''
        else:
            label = '_%d'%idx
        self.add_key('PB'+label,'1e6')
        self.add_key('A1'+label,'0')
        epoch = self.get('PEPOCH')
        self.add_key('T0'+label,epoch)
        self.add_key('ECC'+label,'0')
        return self

    def zero_binary_models(self):
        try:
            self.set('A1',0)
        except:
            return self
        for i in xrange(2,10):
            try:
                self.set('A1_%d'%i,0)
            except:
                return self

    def comp_mass(self,m1=1.4,sini=0.866):
        from scipy.optimize import fsolve
        d = self.get_binary_dict()
        if d is None:
            print 'No binary companion in ephemeris.'
            return
        a1 = d['A1']*29979245800 # semi-major in cm
        pb = d['PB'] # period in s
        kappa = 1.32712440018e26 # heliocentric grav. const in cm^3/s^2
        f = 4*np.pi**2/kappa*(a1)**3/pb**2 # mass function in Msun
        m2_0 = (f*m1**2)**(1./3)/sini # initial guess; good for light comps
        return fsolve( lambda m2: (m2*sini)**3/(m1+m2)**2-f,m2_0)

    def t2_mass_range(self):
        # emulate tempo2's calculation
        print self.comp_mass(m1=1.35,sini=1.0)
        print self.comp_mass(m1=1.35,sini=0.86603)
        print self.comp_mass(m1=1.35,sini=0.43589)

    def eval_freq(self,times,degree=None,t0=None):
        if degree is None:
            degree = self.degree
        if t0 is None:
            t0 = self.get('PEPOCH',type=np.longdouble)
        else: t0 = np.asarray(t0,dtype=np.longdouble)
        times = np.asarray(times,dtype=np.longdouble)
        dts = (times-t0)*86400
        print dts[0]
        multi = np.ones_like(times)
        freq = np.zeros_like(times)
        for i in xrange(0,degree+1):
            tterm = self.get('F%d'%i,type=np.longdouble)
            print tterm,multi[0]
            freq += tterm*multi
            multi *= (dts/(i+1))
        return freq

    def shift_epoch(self,new_epoch,degree=None,t0=None):
        if degree is None:
            degree = self.degree
        if t0 is None:
            t0 = self.get('PEPOCH',type=np.longdouble)
        else: t0 = np.asarray(t0,dtype=np.longdouble)
        dt = (new_epoch-t0)*86400
        multis = 1
        fns = np.asarray([self.get('F%d'%i,type=np.longdouble) for i in xrange(0,degree+1)])
        multis = np.empty_like(fns)
        results = np.empty_like(fns)
        multis[0] = 1
        for i in xrange(1,degree+1):
            multis[i] = multis[i-1] * dt/(i+1)
        for i in xrange(0,degree+1):
            print fns[i:]
            print multis[:degree+1-i]
            print fns[i:]*multis[:degree+1-i]
            results[i] = (fns[i:]*multis[:degree+1-i]).sum()
        return results

    def shift_tzr(self,phase_shift):
        from uw.pulsar import polyco
        t0 = self.get('TZRMJD',type=float)
        pc = polyco.Polyco(self.parfile)
        dt = pc.invert_phase_shift(t0,phase_shift)
        print dt
        sf1 = StringFloat(self.get('TZRMJD'))
        print sf1
        sf2 = StringFloat(dt)
        print sf2
        self['TZRMJD'] = str(sf1+sf2)

    def eval_phase(self,times,degree=None,t0=None):
        degree = degree or self.degree
        if t0 is None: t0 = self.get('PEPOCH',type=np.longdouble)
        else: t0 = np.asarray(t0,dtype=np.longdouble)
        times = np.asarray(times,dtype=np.longdouble)
        dts = (times-t0)*86400
        multi = dts.copy()
        phase = np.zeros_like(times)
        for i in xrange(0,degree+1):
            phase += self.get('F%d'%i,type=np.longdouble)*multi
            multi *= (dts/(i+2))

        # glitch part
        # NB -- no exponential recovery supported (yet)
        for i in xrange(10):
            if 'GLEP_%d'%i not in self.keys(): continue
            g0 = np.asarray(self.get('GLEP_%d'%i,type=np.longdouble))
            if 'GLPH_%d'%i in self.keys():
                phase += self.get('GLPH_%d'%i,type=np.longdouble)
            dts = (times - g0)*86400
            mask = dts > 0
            multi = dts.copy()*mask
            for j in xrange(3):
                key = 'GLF%d_%d'%(j,i)
                if key in self.keys():
                    phase += self.get(key,type=np.longdouble)*multi
                    multi *= (dts/(j+2))
        return phase

    def write(self,output):
        # write to buffer until all keys successfully parsed
        f = deque()
        for key in self.ordered_keys:
            vals = [self[key]]
            if key in self.duplicates.keys():
                vals += self.duplicates[key]
            if key[0] == '#': # handle comments
                key = key.split('COMMENT')[0]
            else:
                key = pad(key,20)
            # ensure a space for long keys; ignore comments
            if (key[-1] != ' ') and (key[0]!='#'):
                key += ' '
            for val in vals:
                if not hasattr(val,'__iter__'):
                    val = [val]
                try:
                    val = ''.join(map(pad26,val))
                except TypeError as e:
                    print 'Failed writing key/val pair %s/%s.'%(key,val)
                    raise e
                f.append('%s%s\n'%(key,val))
        try:
            file(output,'w').write(''.join(f))
        except TypeError:
            # perhaps output is a buffer
            output.write(''.join(f))

    def get_pm(self):
        """ Get the proper motion, return value is of form
            PMRA, PMRA_ERR, PMDEC, PMDEC_ERR.  If no proper
            motion is in ephemeris, all 0s are returned."""
        vals = [0.]*4
        for ikey,key in enumerate(['PMRA','PMDEC']):
            if key in self.keys():
                pmval = self.get(key,first_elem=False)
                if hasattr(pmval,'__iter__'):
                    vals[ikey*2] = float(pmval[0])
                    vals[ikey*2+1] = float(pmval[-1])
                else: vals[ikey*2] = float(pmval)
        return vals

    def get_comp_pm(self):
        """ Return the composite proper motion."""
        vals = self.get_pm()
        v = (vals[0]**2+vals[2]**2)**0.5
        e = ((vals[0]*vals[1])**2+(vals[2]*vals[3])**2)**0.5/v
        return v,e

    def get_shklovskii_pdot(self,dist,velocity=None,
        get_corrected_pdot=False):
        """ Given dist (kpc), try to compute Shklovskii from proper
            motion.  Alternatively, if velocity (km/s) is set, use
            this for computation.  Return the shift in Pdot expected
            from Shklovskii."""
        dist *= 3.08568025e21 # in cm
        if velocity is not None:
            velocty *= 1e5 # cm/2
            return self.p()*velocity**2/dist/C
        v,e = self.get_comp_pm()
        if v == 0.: return 0.
        # convert to rad/s
        s1 = 365.*86400.*1000.*3600.*180/np.pi
        v /= s1
        e /= s1
        print v,e
        s2 = self.p()*dist/C
        if get_corrected_pdot:
            pdot,epdot = self.pdot(error=True)
            cpdot = pdot - s2*v**2
            cepdot = (epdot**2 + (s2*e**2)**2)**0.5
            return cpdot,cepdot
        return s2*v**2

    def get_transverse_velocity(self,dist):
        """ Given dist (kpc), return transverse velocity in km/s. """
        dist *= 3.08568025e21 # in cm
        v,e = self.get_comp_pm()
        scale = (np.pi/180)/3600/(365.24*86400)*dist/1e5
        return v*scale,e*scale

    def get_bfield(self,distance=None):
        """ Return characteristic surface field (G).  If distance (kpc)
            is provided, correct for Shklovskii."""
        pdot = self.pdot()
        if distance is not None:
            pdot -= self.get_shklovskii_pdot(distance)
        return 3.2e19*(pdot*self.p())**0.5

    def get_age(self,distance=None):
        """ Return characteristic age in Myr.  If distance (kpc)
            is provided, correct for Shklovskii."""
        pdot = self.pdot()
        if distance is not None:
            pdot -= self.get_shklovskii_pdot(distance)
        return self.p()/(2.*pdot)/(365*86400)/1e6

    def get_dm_alignment_uncertainty(self,use_last_sig_fig=True,dme=None):
        """ Attempt to determine an uncertainty on the gamma-ray phase based on
            the specified uncertainty in DM.  If no uncertainty is specified,
            return 0.

            If dme is not None, use this value for alignment uncertainty.
            
            Return value is in rotational phase."""
        def get_error(x,factor=3):
            # This function is a somewhat fragile -- assumes Python does the
            # rounding correctly during the cast to float
            l,r = str(float(x)).split('.')
            if len(r)==1 and r=='0':
                return factor * 1
            else:
                return float(factor) / 10**len(r)
        try:
            if dme is None:
                dm = self.get('DM',first_elem=False,type=float)
                if not hasattr(dm,'__len__'):
                    if not use_last_sig_fig: return 0
                    dme = get_error(dm)
                    print 'Using %s as DM error for %s'%(dme,dm)
                else:
                    dme = dm[-1] # should account for .par files with "1/0" in the line
            f0 = self.get('F0',type=float)
            freq = self.get('TZRFRQ',type=float) # in MHz
            if freq==0:
                print 'TZRFRQ == 0!  Assuming 1400 MHz...'
                freq = 1400.
        except KeyError:
            return 0
        return dme/2.41e-4/freq**2*f0

    def get_tres_uncertainty(self):
        """ Attempt to determine the contribution of the timing residuals to
            peak width uncertainties.

            Returns 0 if not able to determine such a thing.

            Return value is in rotational phase."""
        try:
            tres = self.get('TRES',type=float,first_elem=True)
        except KeyError:
            return 0.
        return tres/self.p()*1e-6 # residual in mus

    def is_binary(self):
        return 'BINARY' in self.keys()

    def add_key(self,key,val,allow_duplicates=False,stringify=True):
        """ Insert a key, placed at bottom of file.  If key already
            present, update its entry.  NB will NOT make a duplicate."""
        if key not in self.ordered_keys:
            self.ordered_keys.append(key)
            if stringify:
                self[key] = str(val)
            else:
                self[key] = val
        elif allow_duplicates:
            self.duplicates[key].append(val)
        else:
            if stringify:
                self[key] = str(val)
            else:
                self[key] = val

    def delete_key(self,key):
        try:
            idx = self.ordered_keys.index(key)
            self.ordered_keys.pop(idx)
            self.pop(key)
            if key in self.duplicates.keys():
                self.duplicates.pop(key)
        except ValueError:
            pass

    def change_key(self,key,newkey,swap=True):
        """ Replace one key label with a new one.  If the new key is
            already populated, replace the old key values with those that
            originally belonged to the new key.  (Swap.)"""
        try:
            idx1 = self.ordered_keys.index(key)
            if swap and newkey in self.ordered_keys:
                idx2 = self.ordered_keys.index(newkey)
                self.ordered_keys[idx1],self.ordered_keys[idx2] = \
                    newkey,key
                self[newkey],self[key] = self.pop(key),self.pop(newkey)
            else:
                self.ordered_keys[idx1] = newkey
                self[newkey] = self.pop(key)
            if key in self.duplicates.keys():
                if swap and (newkey in self.duplicate_keys):
                    self.duplicates[newkey],self.duplicates[key] = \
                        self.duplicates.pop(key),self.duplicates.pop(newkey)
                else:
                    self.duplicates[newkey] = self.duplicates.pop(key)
        except ValueError:
            print 'Could not find key %s.'%key
            pass

    def delete_val(self,val):
        """ Attempt to remove an entry by value."""
        while True:
            try:
                vals = [self[k] for k in self.ordered_keys]
                idx = vals.index(val) 
                key = self.ordered_keys.pop(idx)
                self.pop(key)
                if key in self.duplicates.keys():
                    self.duplicates.pop(key)
            except ValueError:
                return

    def has_glitches(self):
        for key in self.keys():
            if 'GLEP' in key:
                return True
        return False

    def get_glitch_index(self,glepoch=None,tol=1e-2):
        """ Get the integer(s) corresponding to the provided glitch 
            epoch.  If no epoch is provided, return all indices"""
        indices = []
        for key in self.keys():
            if key[:4]=='GLEP':
                k_glepoch = self.get(key,type=float)
                if (glepoch is None) or (abs(glepoch-k_glepoch) < tol):
                    indices.append(int(key.split('_')[-1]))
        return sorted(indices)

    def zero_glitches(self,glepoch=None,transients_only=False):
        """ If glepoch is provided, only zero glitches near (within 1d) of
            that epoch.  Otherwise, remove all glitches."""
        no_glitches = True
        indices = self.get_glitch_index(glepoch=glepoch)
        for key in self.keys():
            if (key[:2]=='GL') and (int(key.split('_')[-1]) in indices):
                no_glitches = False
                if key[2:4] != 'EP':
                    if ((not transients_only) or 
                        ((key[2:5]=='F0D') or (key[2:4]=='TD'))):
                        self.set(key,0)
        return no_glitches

    def zero_glitch_transients(self,glepoch=None):
        self.zero_glitches(glepoch=glepoch,transients_only=True)

    def get_glitch_parameters(self,index,strip_index=True):
        indices = self.get_glitch_index(glepoch=None)
        if index not in indices:
            raise IndexError('Do not have glitch %d!'%index)
        pars = dict()
        for key in self.keys():
            if (key[:2]=='GL') and (int(key.split('_')[-1])==index):
                val = self[key]
                if hasattr(val,'__iter__'):
                    if len(val) > 1:
                        pars[key] = [float(val[0]),float(val[-1])]
                    else:
                        pars[key] = [float(val[0]),None]
                else:
                    pars[key] = [float(val),None]
        if strip_index:
            d = dict()
            for key in pars.keys():
                d[key.split('_')[0]] = pars[key]
            pars = d
        return pars

    def sort_glitches(self):
        """ Put glitch indices in time order."""
        indices = []
        epochs = []
        for key in self.ordered_keys:
            if key.startswith('GLEP'):
                print key
                indices.append(int(key.split('_')[-1]))
                epochs.append(self.get(key,type=float))

        # perform a stable sort on epoch
        new_indices = np.argsort(epochs,kind='mergesort')+1

        for key in self.keys():
            if key.startswith('GL'):
                base,idx = key.split('_')
                new_idx = new_indices[indices.index(int(idx))]
                self.change_key(key,base+'_%d'%new_idx)
                print key,base+'_%d'%new_idx

    def glitch_phase_jump(self,glepoch):
        """ Compute the phase jump associated with a particular glitch."""
        indices = self.get_glitch_index(glepoch=glepoch)
        phase = 0
        for idx in indices:
            if 'GLPH_%d'%idx in self.keys():
                phase += self.get('GLPH_%d'%idx,type=float)
            if 'GLTD_%d'%idx in self.keys():
                phase += self.get('GLTD_%d'%idx,type=float)* \
                    86400*self.get('GLF0D_%d'%idx,type=float)
        return phase

    def delete_glitch(self,glepoch):
        indices = self.get_glitch_index(glepoch=glepoch)
        glphase = self.glitch_phase_jump(glepoch)
        if len(indices)==0:
            return glphase
        # remove overlapping keys
        for key in self.keys():
            if (key[:2]=='GL'):
                idx = int(key.split('_')[-1])
                if idx in indices:
                    self.delete_key(key)
        # finally, relabel & sort remaining glitches
        self.sort_glitches()
        return glphase

    def delete_glitches(self):
        for key in self.ordered_keys[:]:
            if key[:2] == 'GL':
                self.pop(key)
                self.ordered_keys.remove(key)

    def get_glepochs(self,unique=True):
        glepochs = []
        for key in self.keys():
            if key[:4] == 'GLEP':
                glepochs.append(self.get(key,type=np.float128))
        if unique:
            return sorted(list(set(glepochs)))
        return glepochs

    def get_nglitch(self,thresh=0):
        """ Return number of glitches above a df/f threshold."""
        count = 0
        f0 = self.get('F0',type=float)
        for i in xrange(100):
            try:
                pars = self.get_glitch_parameters(i+1)
                df_on_f = abs(pars['GLF0'][0] / f0)
                if df_on_f > thresh:
                    count += 1
            except IndexError:
                break
        return count
            

    def has_waves(self):
        return self.num_waves() > 0

    def num_waves(self):
        f = lambda s: s.startswith('WAVE') and s[-1].isdigit()
        return len(filter(f,self.keys()))

    def delete_waves(self):
        for key in self.ordered_keys[:]:
            if ('WAVE' in key):
                self.pop(key)
                self.ordered_keys.remove(key)

    def get_wave_string(self,epoch=None):
        """ Return a sting suitable for appending to an ephemeris giving
            a zeroed wave contribution with the same number of terms as
            the current ephemeris."""
        keys = []
        has_waves = False
        for key in self.ordered_keys:
            if 'WAVE' in key:
                has_waves = True
                if 'EPOCH' not in key:
                    keys.append(key)
        if not has_waves:
            return ''
        epoch = epoch or self['WAVEEPOCH']
        return 'WAVEEPOCH %s\n'%(epoch) + \
            '\n'.join(map(lambda s: '%s 0 0'%s,keys))

    def replace_astrometry(self,block,output=None,clobber=True):
        """ Replace the astrometry in the ephemeris with the string
            indicated in block.  Note that we wish to preserve it verbatim
            since it may include comments and errors.  Therefore, delete
            any existing conflicting entries, then prepend this to the
            file and reload it.  Example of an astrometry block:

            # from Dodson et al. VLBI
            RAJ    08:35:20.61149  0.00002
            DECJ   -45:10:34.8751  0.00030

            NB this WILL clobber the ephemeris file unless otherwise
            specified.
        """
        if block is None:
            return
        if len(block.strip()) == 0:
            return
        lines = filter(lambda l: len(l) > 0, block.split('\n'))
        # remove all duplicate keys
        keys = [l.split()[0] for l in lines if l[0] != '#']
        map(self.delete_key,keys)
        # remove duplicate comments (e.g. if this has been done before)
        comm = [l.strip()[1:] for l in lines if l[0] == '#']
        map(self.delete_val,comm)

        import StringIO
        buf = StringIO.StringIO()
        self.write(buf)
        buf.seek(0)
        lines = map(str.strip,map(str,buf.readlines()))
        idx = 0
        for line in lines:
            if 'PSRJ' in line:
                break
            idx += 1
        buf = StringIO.StringIO()
        buf.write('\n'.join(
            lines[:idx+1]+block.split('\n')+lines[idx+1:]))

        self.init(buf=buf)

        if output is not None:
            self.write(output)
        elif clobber:
            self.write(self.parfile)

    def freeze_params(self):
        for key in self.keys():
            if hasattr(self[key],'__iter__') and (len(self[key])>1) and self[key][1]=='1':
                self[key][1] = '0'
        return self # for chaining

    def delete_smooth(self):
        keys = filter(lambda x: 'IFUNC' in x, self.ordered_keys)
        for key in keys:
            self.delete_key(key)
        return self # for chaining

    def set_smooth(self,on=False):
        if 'SIFUNC' in self.ordered_keys:
            self['SIFUNC'][0] = '2' if on else '0'
        return self

class TimFile(object):

    def __init__(self,timfile):
        efac = 1
        self.toas = toas = []
        self.toa_lines = toa_lines = []
        lines = filter(lambda x:len(x) > 0,
            map(str.strip,file(timfile).readlines()))
        for line in lines:
            toks = line.strip().split()
            if toks[0] == 'EFAC':
                efac = float(toks[1])
                continue
            if id_toa_line(line):
                # is a TOA
                toa =   {
                        'TEL':str(toks[0]),
                        'FREQ':float(toks[1]),
                        'MJD':float(toks[2]),
                        'ERR':float(toks[3])*efac,
                        }
                toas += [toa] 
                toa_lines += [line]

    def get_mjds(self):
        return np.asarray([toa['MJD'] for toa in self.toas],dtype=np.float128)

    def get_errs(self,days=False):
        """ days or microseconds"""
        if not days:
            return np.asarray([toa['ERR'] for toa in self.toas],dtype=float)
        return np.asarray([toa['ERR'] for toa in self.toas],dtype=np.float128)*(1./(86400*1e6))

    def time_order(self,output=None):
        a = np.argsort(self.get_mjds())
        toa_lines = np.asarray(self.toa_lines)[a]
        if output is not None:
            file(output,'w').write('FORMAT 1\n'+''.join(toa_lines))

    def get_frqs(self):
        return np.asarray([toa['FREQ'] for toa in self.toas],dtype=float)

def parfunc(par,func):
    """ Instantiate a temporary ParFile and execute the function,
        then write it out.
    """
    pf = ParFile(par)
    eval('pf.%s()'%func)
    pf.write(par)

def add_key(par,key,val):
    """ Shortcut method for adjusting an ephemeris."""
    pf = ParFile(par)
    pf.add_key(key,val)
    pf.write(par)

def copy_tzr(par1,par2):
    """ Update the TZR parameters of par2 with those of par1."""
    keys = ['TZRMJD','TZRFRQ','TZRSITE']
    pf1 = ParFile(par1)
    pf2 = ParFile(par2)
    for key in keys:
        if key in pf1.keys():
            pf2.add_key(key,pf1.get(key))
    pf2.write(par2)
        
def get_bats_etc(par,tim,output=None,full_output=False,binary=False,
    nofit=False):
    """ Use the tempo2 general plugin to compute the bats and absolute
        phases of a set of TOAs."""
    if not os.path.isfile(par):
        raise IOError('Ephemeris %s is not a valid file!'%par)
    if not os.path.isfile(tim):
        raise IOError('TOA collection %s is not a valid file!'%tim)
    binary = binary and ParFile(par).is_binary()
    cmd = """tempo2 -output general2 -s "onerous\t{bat}\t{err}\t{npulse}\t{pre_phase}\t{sat}\n" -f %s %s"""%(par,tim)
    if binary:
        cmd = cmd.replace('bat','bbat')
    if nofit:
        cmd += ' -nofit'
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    # TODO -- should check output...
    toks = map(lambda l:l.split('\t')[1:],
        filter(lambda l:l[:7]=='onerous',proc.stdout))
    #toks = [line.split('\t')[1:] for line in proc.stdout if line[:7]=='onerous']
    bats = np.array([x[0] for x in toks],dtype=np.float128)
    errs = np.array([x[1] for x in toks],dtype=np.float128)*(1e-6/86400)
    phas = np.array([x[2] for x in toks],dtype=np.float128)
    offs = np.array([x[3] for x in toks],dtype=np.float128)
    sats = np.array([x[4] for x in toks],dtype=np.float128)
    #phas += np.round(offs)
    if output is not None:
        outstring = '\n'.join(['%.20f %.20f %d'%(b,e,p) for b,e,p in zip(bats,errs,phas)])
        file(output,'w').write('# %s %s\n'%(par,tim)+outstring)
    #phas += offs # total phase
    if full_output:
        return bats,errs,phas,offs,sats
    else:
        return bats,errs,phas

def get_resids(par,tim,emax=None,phase=False,get_mjds=False,latonly=False,
    jitter=None,select=None,nofit=True,echo=False):
    """ Use the tempo2 general2 plugin to obtain residuals.  By default, 
        the residuals and the errors are given in microseconds.

        phase -- convert residuals and errors to phase units
        get_mjds -- return the site arrival times in MJD too
        latonly -- filter TOAs that aren't at 0 frequency
    """
    if not os.path.isfile(par):
        raise IOError('Ephemeris %s is not a valid file!'%par)
    if not os.path.isfile(tim):
        raise IOError('TOA collection %s is not a valid file!'%tim)
    cmd = """tempo2 -output general2 -s "onerous\t{err}\t{post}\t{bat}\t{freq}\n" -f %s %s"""%(par,tim)
    if nofit:
        cmd += ' -nofit'
    if select is not None:
        cmd += ' -select %s'%select
    if echo:
        print cmd
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    toks = [line.split('\t')[1:] for line in proc.stdout if line[:7]=='onerous']
    if len(toks)==0:
        # tempo2 failed us
        if get_mjds:
            return [None]*6
        return [None]*4
        
    # NB -- both residuals and errors in microseconds
    frqs = np.array([x[3] for x in toks],dtype=np.float128)
    if latonly:
        m = frqs==0
    else:
        m = np.asarray([True]*len(frqs))
    frqs = frqs[m]
    errs = np.array([x[0] for x in toks],dtype=np.float128)[m]
    resi = np.array([x[1] for x in toks],dtype=np.float128)[m]*1e6
    mjds = np.array([x[2] for x in toks],dtype=np.float128)[m]
    # add any jitter noise if specified
    if jitter is not None:
        # get observation lengths
        tobs = np.asarray(map(lambda x: float(x) if len(x)>0 else 1e6,
            flag_values(tim,'length')))
        print errs.shape,jitter,tobs.shape
        x = (errs**2 + (1e6*jitter)**2/tobs)**0.5
        errs = x
    # if we are restricting large error bars, remove their contribution
    # to the RMS
    if phase:
        # convert to phase
        p = ParFile(par).p()*1e6 # period in mus
        errs /= p
        resi /= p
    # apply emax *after* phase correction to be consistent
    if emax is not None:
        mask = errs < emax
        dof = mask.sum()
    else:
        dof = len(resi)
        mask = np.asarray([True]*dof)
    if not np.any(mask):
        chi2 = 0
    else:
        # I think we always want to subtract the weighted mean
        resi -= np.average(resi[mask],weights=1./errs[mask]**2)
        chi2 = ((resi[mask]/errs[mask])**2).sum()
    if get_mjds:
        return resi,errs,chi2,dof,mjds,frqs
    return resi,errs,chi2,dof

def set_sifuncs(par,tim):
    """ Use interpolation to remove any residual.

    The idea with this at the moment is simply to whiten a file with an
    large red noise component.
    """

    # (1) get site arrival times and residuals
    cmd = """tempo2 -output general2 -nofit -s "onerous\t{pre}\t{sat}\n" -f %s %s"""%(par,tim)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    toks = [l.split('\t')[1:] for l in proc.stdout if l[:7]=='onerous']
    # NB -- both residuals and errors in microseconds
    resi = np.array([x[0] for x in toks],dtype=np.float128)
    sats = np.array([x[1] for x in toks],dtype=np.float128)

    # clone end points to handle edge effects
    sats = np.concatenate([[sats[0]-1],sats,[sats[-1]+1]])
    resi = np.concatenate([[resi[0]],resi,[resi[-1]]])

    pf = ParFile(par)
    pf.add_key('SIFUNC','2 0') # tells tempo2 to do linear interp
    for i in xrange(len(resi)):
        if i==0:
            sat = sats[i] - 1e-8/86400 # make 10ns early
        elif i== len(resi)-1:
            sat = sats[i] + 1e-8/86400 # make 10ns late
        else:
            sat = sats[i]
        pf.add_key('IFUNC%d'%(i+1),'%.20f %.10f 0.0'%(sat,-resi[i]))
    pf.write(par)

def get_toa_strings(tim):
    """ Return a list of strings for all TOA entries."""
    lines = filter(lambda l: len(l) > 0,
        map(str.strip,file(tim).readlines()))
    if not lines[0].startswith('FORMAT 1'):
        raise ValueError('Cannot parse non-tempo2 style TOA files.')
    lines = filter(id_toa_line,lines)
    if (len(lines)>0) and (not lines[-1].endswith('\n')):
        lines[-1] = lines[-1] + '\n'
    return lines

def tim_filter(tim,thresh=5,output=None):
    """ Parse a likelihood file and attempt to comment out TOAs that do
        not pass a likelihood (detection) threshold.
    """
    output = output or tim
    lines = file(tim).readlines()
    new_lines = deque()
    for iline,line in enumerate(lines):
        new_lines.append(line)
        if line[0] != ' ':
            continue
        toks = line.split()
        try:
            idx = toks.index('-logl')
            logl = -float(toks[idx+1])
            # if doesn't meet threshold, comment out old value and replace
            # with one with a very large error; this allows us to keep
            # "uniform" TOAs if necessary but deweight them in fits
            if logl < thresh:
                new_lines.pop()
                new_lines.append('C ' + line[1:])
                toks[3] = '1e7'
                nl = ' '+' '.join(toks)+'\n'
                new_lines.append(nl)
        except ValueError:
            continue
    file(output,'w').write(''.join(new_lines))

def tim_faker(tim,tim_shifts,thresh=5,output=None):
    """ Parse a likelihood file and attempt to comment out TOAs that do
        not pass a likelihood (detection) threshold.
    """
    output = output or tim
    lines = file(tim).readlines()
    new_lines = deque()
    toa_counter = -1
    for iline,line in enumerate(lines):
        new_lines.append(line)
        if line[0] != ' ':
            continue
        toa_counter += 1
        toks = line.split()
        try:
            idx = toks.index('-logl')
            logl = -float(toks[idx+1])
            # if doesn't meet threshold, comment out old value and replace
            # with one with a very large error; this allows us to keep
            # "uniform" TOAs if necessary but deweight them in fits
            if logl < thresh:
                new_lines.pop()
                new_lines.append('C ' + line[1:])
                toks[2] = '%.15f'%(np.array(toks[2],dtype=np.float128) + tim_shifts[toa_counter])
                toks[3] = '1e7'
                nl = ' '+' '.join(toks)+' -fake \n'
                new_lines.append(nl)
        except ValueError:
            continue
    file(output,'w').write(''.join(new_lines))

def strip_radio(tim,output=None):
    """ Remove non-LAT TOAs from a file based on frequency."""
    toas = filter(lambda toa: int(float(toa.split()[1]))==0,
        get_toa_strings(tim))
    output = output or tim
    file(output,'w').write(''.join(toas))

def id_toa_line(line):
    """ Determine if the line from a .tim file corresponds to an 
        (uncommented) TOA.
        
        Primarily, check to see if the third column is consistent with an MJD.
    """
    if line[0]=='#' or line[0]=='C':
        return False
    if 'FORMAT' in line:
        return False
    if 'MODE' in line:
        return False
    toks = line.split()
    if len(toks) < 3:
        return False
    mjd = toks[2].split('.')
    if len(mjd)==2 and len(mjd[0])==5:
        return True
    return False

def id_toa_comment(line):
    return line[0]=='#' or line[0]=='C'

def add_flag(tim,flag,vals,output=None):
    """ Add a flag to the TOAs in an output file with the specified value.
        If the flag is already present, the value is replaced.
        
        The user may specify either a single value for all TOAs (e.g. a 
        flag to specify a jump) or a value for each and every TOA.
    """
    lines = file(tim).readlines()
    new_lines = deque()
    counter = 0
    flagstr = str(flag)
    if not flag.startswith('-'):
        flagstr = '-%s'%flagstr
    if not hasattr(vals,'__iter__'):
        vals = [vals]*len(lines)
        one2one = False
    else:
        one2one = True
    for iline,line in enumerate(lines):
        new_lines.append(line)
        if not id_toa_line(line):
            continue
        new_lines.pop()
        ph = str(vals[counter])
        toks = line.split()
        try:
            idx = toks.index(flagstr)
            toks[idx+1] = ph
            new_lines.append(' '+' '.join(toks)+'\n')
        except ValueError:
            new_lines.append(line.strip('\n')+' %s %s\n'%(flagstr,ph))
        counter += 1
    if one2one and (counter != len(vals)):
        raise ValueError('Found %d vals for %d TOAs!'%(len(vals),counter))
    output = output or tim
    file(output,'w').write(''.join(new_lines))

def del_flag(tim,flag,output=None):
    """ Remove a flag from a set of observations.  It is assumed that the
        flags follow the standard flag/val pattern."""
    lines = file(tim).readlines()
    new_lines = deque()
    if flag[0] != '-':
        flag = '-%s'%flag
    for iline,line in enumerate(lines):
        new_lines.append(line)
        if not id_toa_line(line):
            continue
        new_lines.pop()
        toks = line.split()
        try:
            print toks
            idx = toks.index(flag)
            print idx
            if idx==len(toks)-2:
                line = ' '.join(toks[:idx])
            else:
                line = ' '.join(toks[:idx] + toks[idx+1:])
            new_lines.append(' %s\n'%line)
        except ValueError:
            new_lines.append(line)
    output = output or tim
    file(output,'w').write(''.join(new_lines))

def add_phase(tim,abs_phase,output=None):
    """ Add an absolute phase flag to a TOA file.  abs_phase is a vector
        containing the (integer) pulse number for each TOA."""
    if type(abs_phase)==type(''):
        # phase is an ephemeris
        bats,errs,abs_phase,offs,sats = get_bats_etc(
            abs_phase,tim,full_output=True)
    abs_phase = [int(x) for x in abs_phase]
    add_flag(tim,'pn',abs_phase,output=output)

def mask_tim(tim,mask,output=None,verbose=False,comment=False):
    """ Make TOAs in a .tim file.  Two methods available:

    (1) [default] set error to a large value
    (2) comment out TOA
    """
    output = output or tim
    lines = file(tim).readlines()
    new_lines = deque()
    masked_lines = deque()
    counter = 0
    for iline,line in enumerate(lines):
        new_lines.append(line)
        if not id_toa_line(line):
            continue
        if mask[counter]:
            if verbose:
                print 'Masking line:'
                print line
            masked_lines.append(new_lines.pop())
            if not comment:
                toks = line.split()
                toks[3] = '1e7'
                nl = ' '+' '.join(toks)+'\n'
            else:
                nl = '#' + line
            new_lines.append(nl)
        counter += 1
    file(output,'w').write(''.join(new_lines))
    return masked_lines

def replace_nondetections(tim,output=None):
    """ Replace nondetections with errors drawn from rest of population."""
    tf = TimFile(tim)
    mask = tf.get_errs() < 1e6
    errs = tf.get_errs()[tf.get_errs() < 1e6]
    mask = ~mask

    output = output or tim
    lines = file(tim).readlines()
    new_lines = deque()
    counter = 0
    for iline,line in enumerate(lines):
        new_lines.append(line)
        if not id_toa_line(line):
            continue
        if mask[counter]:
            new_lines.pop()
            toks = line.split()
            toks[3] = '%.5f'%(errs[int(np.random.rand(1)*len(errs))])
            nl = ' '+' '.join(toks)+'\n'
            new_lines.append(nl)
        counter += 1
    file(output,'w').write(''.join(new_lines))

def cut_tim(tim,tmin=None,tmax=None,output=None):
    """ Trim out TOAs before tmin / after tmax (MJD).  This differs slighly
        from get_toa_strings + time cut in that commented lines and other
        instructions are preserved."""
    tmin = tmin or -np.inf
    tmax = tmax or np.inf
    tim_strings = deque()
    lines = file(tim).readlines()
    if 'FORMAT 1' not in lines[0]:
        raise ValueError('Cannot parse non-tempo2 style TOA files.')
    for line in lines:
        tim_strings.append(line)
        if id_toa_line(line):
            mjd = float(line.split()[2])
            if not((mjd >= tmin) and (mjd < tmax)):
                tim_strings.pop()
    if output is not None:
        file(output,'w').write(''.join(tim_strings))
    return tim_strings

def flag_filter(tim,flag,valfunc=lambda x: True,return_lines=False,
    converse=False):
    """ Sort through a TOA file and return a mask containing those lines
        that have a flag with the correct value."""
    if flag[0] != '-':
        flag = '-'+flag
    toa_strings = get_toa_strings(tim)
    def f(s):
        s = s.split()
        try:
            idx = s.index(flag)
            return idx==(len(s)-1) or valfunc(s[idx+1])
        except ValueError:
            return False
    mask = np.argwhere(map(f,toa_strings)).flatten()
    if converse:
        mask = ~mask
    if return_lines:
        return [toa_strings[m] for m in mask]
    return mask

def flag_values(tim,flag):
    """ Return the "argument" of a flag if present in a TOA, else an
        empty string."""
    if flag[0] != '-':
        flag = '-'+flag
    toa_strings = get_toa_strings(tim)
    def f(s):
        s = s.split()
        try:
            idx = s.index(flag)
            if idx >= (len(s)-1):
                return ''
            return s[idx+1]
        except ValueError:
            return ''
    return map(f,toa_strings)

def flag_filter_tim(tim,flag,output,val=None,converse=False):
    """ Write out a new .tim file with TOAs with the specified flag.

    Optionally, require the flag to be set to a particular value.
    """
    if val is None:
        flagfunc = lambda x: True
    else:
        flagfunc = lambda x: x == val
    tim_strings = flag_filter(tim,flag,flagfunc,return_lines=True,
        converse=converse)
    file(output,'w').write('FORMAT 1\n%s'%('\n'.join(tim_strings)))

def merge_tim(timfiles,output,tmin=None,tmax=None):
    """ Merge TOA files, applying time cuts if desired."""
    f = lambda x: 'FORMAT' not in x
    g = lambda x: ''.join(filter(f,cut_tim(x,tmin=tmin,tmax=tmax)))
    tim_strings = ''.join(map(g,timfiles))
    file(output,'w').write('FORMAT 1\n%s'%(tim_strings))

def delete_jumps(par,flags,vals):
    """ Delete one or more JUMP parameter from an ephemeris matching the
        arguments.  E.g., flags=['-i'],vals=['parkes_20cm'] would remove
        JUMP -i parkes_20cm XXXXXX
    """
    pf = ParFile(par)
    if 'JUMP' not in pf.keys():
        return
    flags = map(lambda x: x if x[0]=='-' else '-'+x,flags)
    # read all of the jumps into a list
    jumps = [pf['JUMP']]
    if 'JUMP' in pf.duplicates.keys():
        for x in pf.duplicates['JUMP']:
            jumps.append(x)
    pf.delete_key('JUMP')
    new_jumps = []
    for jump in jumps:
        ok = True
        for f,v in zip(flags,vals):
            if jump[0]==f and jump[1]==v:
                ok = False
                break
        if ok:
            new_jumps.append(jump)
    for jump in new_jumps:
        pf.add_key('JUMP',' '.join(jump),allow_duplicates=True)
    pf.write(par)

def compute_jump(par,tim,flags,vals,tmin=None,tmax=None,add_jumps=False,
    clear_jumps=False):
    """ Compute the jump parameter (microseconds) by computing the mean
        (error-weighted) offset of the residuals.
        
        How this will work -- we have a .par file and a merged TOA file,
        and each TOA to be included must have a flag/value pair to allow
        us to select on them.  By default, all jumps are computed relative
        to the first flag/value pair.
        
        add_jumps -- add the resulting jumps to the .par file
        clear_jumps -- clear any conflicting jumps from the .par file
    """
    if not os.path.isfile(par):
        raise IOError('Ephemeris %s is not a valid file!'%par)
    if not os.path.isfile(tim):
        raise IOError('TOA collection %s is not a valid file!'%tim)
    if len(flags)<2:
        print 'Fewer than two flags.'
        return 0
    
    if clear_jumps:
        delete_jumps(par,flags,vals)
    # construct a temporary TOA file to handle case of MWL files with TOAs
    # extending before/beyond the domain of validity
    tf = tempfile.mkstemp()
    handle,fname = tf
    cut_tim(tim,tmin=tmin,tmax=tmax,output=fname)
    cmd = """tempo2 -output general2 -s "onerous\t{sat}\t{err}\t{pre}\n" -f %s %s"""%(par,fname)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    toks = [line.split('\t')[1:] for line in proc.stdout if line[:7]=='onerous']
    sats = np.array([x[0] for x in toks],dtype=np.float128)
    errs = np.array([x[1] for x in toks],dtype=np.float128) # in mus
    resi = np.array([x[2] for x in toks],dtype=np.float128) # in s

    # now, make the masks for each set of flags
    means = np.empty(len(flags))
    for i in xrange(len(flags)):
        if vals[i] is not None:
            valfunc = lambda s: s==vals[i]
        else:
            valfunc = lambda s: True
        mask = flag_filter(fname,flags[i],valfunc)
        if not np.any(mask):
            means[i] = np.nan
        else:
            means[i] = np.average(resi[mask],weights=(errs**-2)[mask])
    jumps = [means[0]-means[i] for i in xrange(1,len(means))]
    os.remove(fname)
    if add_jumps:
        # Add the parameters to the .par file.
        pf = ParFile(par)
        for i in xrange(1,len(flags)):
            # don't write bad values or re-write 0 JUMPs
            if np.isnan(jumps[i-1]) or (abs(jumps[i-1])<1e-9):
                continue
            flag = str(flags[i])
            if flag[0] != '-':
                flag = '-'+flag
            val = str(vals[i] or 1)
            t = [flag,val,'%.8f'%jumps[i-1]]
            print 'Adding %s.'%(' '.join(t))
            pf.add_key('JUMP',' '.join(t),allow_duplicates=True)
        pf.write(par)
    return jumps

def istempo2(tim):
    return file(tim).next().startswith('FORMAT 1')

def parkes2tempo2(tim,output):
    """ VERY crude conversion of a Parkes-style TOA file to tempo2 format.
        Currently ignores everything but crucial data and is only tested
        on Ryan's P574 TOAs."""
    toa_strings = filter(id_toa_line,file(tim).readlines())
    new_strings = deque()
    for toa in toa_strings:
        toks = toa.split()
        new_strings.append(' %s %s %s %s %s\n'%(
            toks[0],toks[3],toks[4],toks[6],toks[7]))
    file(output,'w').write('FORMAT 1\n'+''.join(new_strings)) 

def modify_pn(tim,offset,output=None):
    """ Add a constant to the pulse numbers in a TOA file.

    Useful for joining timing solutions.
    """
    strings = deque()
    init = None
    for line in file(tim).readlines():
        toks = line.split()
        if '-pn' in toks:
            idx = toks.index('-pn')
            if init is None:
                init = int(toks[idx+1])
            toks[idx+1] = str(int(toks[idx+1])-init+offset)
            strings.append(' '.join(toks))
        else:
            strings.append(line)
    if output is None:
        output = tim
    file(output,'w').write('\n'.join(strings))

def add_instrument_flag(tim,output=None):
    """ Do a simple flagging."""
    lines = deque()
    for line in map(str.strip,file(tim).readlines()):
        if not id_toa_line(line):
            lines.append(line)
            continue
        flag = (line.split()[0][0]).upper()
        lines.append(line + ' -i %s'%flag)
    if output is None:
        output = tim
    file(output,'w').write('\n'.join(lines))

def process_ft_data(tim,output=None):
    """ Flag older timing data to allow for JUMPs."""
    lines = deque()
    insts = ['C','P','D','F']
    for line in map(str.strip,file(tim).readlines()):
        if not id_toa_line(line):
            lines.append(line)
            continue
        toks = line.split()
        inst = (toks[0][0]).upper()
        if inst not in insts:
            lines.append(line)
            continue
        if '-i' not in toks:
            toks.extend(('-i',inst))

        if '-ft' not in toks:
            mjd = int(float(toks[2]))
            if inst=='C':
                if (mjd < 49888):
                    # baseline jump
                    pass
                elif  (mjd < 50462):
                    toks.extend(('-ft','C95'))
                else:
                    toks.extend(('-ft','C97'))
            elif inst=='D':
                if mjd < 50462:
                    # baseline jump
                    pass
                else:
                    toks.extend(('-ft','D97'))
        lines.append(' '.join(toks))
    if output is None:
        output = tim
    file(output,'w').write('\n'.join(lines))

def edit_start_finish(par,output=None):
    """ Delete the "free parameter" flag after START/FINISH."""
    lines = file(par).readlines()
    for iline,line in enumerate(lines):
        if ('START' in line) or ('FINISH' in line):
            toks = line.strip().split()
            if not ((toks[0] == 'START') or (toks[0] == 'FINISH')):
                continue
            if len(toks) > 2:
                otoks = line.split(toks[1])
                lines[iline] = otoks[0] + toks[1] + '\n'
    if output is None:
        output = par
    file(output,'w').write(''.join(lines))
