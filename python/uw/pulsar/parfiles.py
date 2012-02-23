import numpy as np
import os
from uw.utilities.coords import ec2eq

C = 29979245800.

def sex2dec(s,mode='ra'):
    toks = s.split(':')
    multi = (-1 if '-' in s else 1) * (15 if mode == 'ra' else 1)
    return multi*sum( (60**(-i) * abs(float(toks[i])) for i in xrange(len(toks)) ) )

def ra2dec(s): return sex2dec(s,mode='ra')
def decl2dec(s): return sex2dec(s,mode='decl')

class ParFile(dict):

    def __init__(self,parfile):
        if not os.path.exists(parfile):
            raise IOError('Indicated file %s does not exist.'%parfile)
        self.parfile = parfile
        self.ordered_keys = []
        for line in file(parfile):
            tok = line.strip().split()
            if len(tok)==0: continue
            known_key = False
            for key in self.ordered_keys:
                if tok[0] == key: known_key = True
            if known_key:
                if len(self[tok[0]][-1]) == 1: self[tok[0]] = [self[tok[0]]]
                self[tok[0]] += [tok[1:]]
            else:
                self.ordered_keys.append(tok[0])               
                self[tok[0]] = tok[1:] if (len(tok[1:]) > 1) else tok[1:][0]
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
        t = self[key]; v = '%.18g'%(val)
        if hasattr(t,'__iter__'):
            self[key][0] = v
        else: self[key] = v

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
            elong, elat = self.get("ELONG",type=float), self.get("ELAT",type=float)
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

    def has_waves(self):
        return np.any(['WAVE' in key for key in self.keys()])

    def is_msp(self,maxp=0.03,maxpdot=1e-18):
        return (self.p() < maxp) and (self.pdot() < maxpdot) and (not self.has_waves())

    def get_time_cuts(self):
        import datetime
        import time
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
        f = file(output,'w')
        for key in self.ordered_keys:
            val = self[key]
            if hasattr(val,'__iter__'):
                val = '  '.join(val)
            key = key + ' '*(15-len(key))
            f.write('%s%s\n'%(key,val))
        f.close()

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

    def get_shklovskii_pdot(self,dist,velocity=None):
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
        v /= 365.*86400.*1000.*3600.*180/np.pi # rad/s
        return self.p()*v**2*dist/C

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

class TimFile(object):

    def __init__(self,timfile):
        efac = 1
        toas = []
        for line in file(timfile):
            toks = line.strip().split()
            if toks[0] == 'EFAC':
                efac = float(toks[1])
                continue
            if line[0] == ' ':
                # is a TOA
                toa =   {
                        'TEL':str(toks[0]),
                        'FREQ':float(toks[1]),
                        'MJD':float(toks[2]),
                        'ERR':float(toks[3])*efac,
                        }
                toas += [toa] 
        self.toas = toas

    def get_mjds(self):
        return np.asarray([toa['MJD'] for toa in self.toas],dtype=float)

    def get_errs(self):
        return np.asarray([toa['ERR'] for toa in self.toas],dtype=float)
