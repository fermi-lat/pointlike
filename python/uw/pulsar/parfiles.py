import numpy as np

def sex2dec(s,mode='ra'):
    toks = s.split(':')
    multi = (-1 if '-' in s else 1) * (15 if mode == 'ra' else 1)
    return multi*sum( (60**(-i) * abs(float(toks[i])) for i in xrange(len(toks)) ) )

def ra2dec(s): return sex2dec(s,mode='ra')
def decl2dec(s): return sex2dec(s,mode='decl')

class ParFile(dict):

    def __init__(self,parfile):
        self.parfile = parfile
        self.ordered_keys = []
        for line in file(parfile):
            tok = line.strip().split()
            if len(tok)==0: continue
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

    def get_ra(self): return ra2dec(self.get('RAJ'))
    def get_dec(self): return decl2dec(self.get('DECJ'))
    def get_skydir(self):
        from skymaps import SkyDir
        ra = ra2dec(self.get('RAJ'))
        dec = decl2dec(self.get('DECJ'))
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

    def edot(self,mom=1e45):
        om = self.get('F0',type=float)*(np.pi*2)
        omdot = self.get('F1',type=float)*(np.pi*2)
        return -mom*om*omdot # say this without laughing

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
        if degree is None:
            degree = self.degree
        if t0 is None:
            t0 = self.get('PEPOCH',type=np.longdouble)
        else: t0 = np.asarray(t0,dtype=np.longdouble)
        times = np.asarray(times,dtype=np.longdouble)
        dts = (times-t0)*86400
        multi = dts.copy()
        phase = np.zeros_like(times)
        for i in xrange(0,degree+1):
            tterm = self.get('F%d'%i,type=np.longdouble)
            print tterm,multi[0]
            phase += tterm*multi
            multi *= (dts/(i+2))
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

def shift_pepoch(self,newepoch):
    import pylab as pl
    t0 = self.get('START',type=np.float128)
    t1 = self.get('FINISH',type=np.float128)
    T0 = self.get('PEPOCH',type=np.float128)
    dom = np.linspace(t0,t1,100).astype(np.longdouble)
    cod = self.eval_freq(dom)
    #dom2 = dom-newepoch
    cod2 = self.eval_freq(dom,t0=newepoch)
    return cod,cod2
    pl.plot(dom,cod2-cod)
    """
    p = np.polyfit(dom2.astype(float),cod.astype(float),self.degree)
    print p
    pl.plot(dom,(cod-np.polyval(p,dom2))*86400)
    for i in xrange(self.degree+1):
        self.set('F%d'%i,p[-(i+1)])     
    self.set('PEPOCH',newepoch)
    """

def shift_pepoch2(self,newepoch):
    N = self.degree + 1
    SECSPERDAY = 86400.
    if N == 1:
        self.set('PEPOCH',newepoch)
        return
        
    import pylab as pl
    t0 = self.get('START',type=np.longdouble)
    t1 = self.get('FINISH',type=np.longdouble)
    T0 = self.get('PEPOCH',type=np.longdouble)

    # sample phase on num points equal to num freqs
    dom = np.linspace(t0,t1,N).astype(np.longdouble)
    dt1 = SECSPERDAY*(dom - T0)
    dt2 = SECSPERDAY*(dom - newepoch)
    phi = np.zeros_like(dom)

    # build matrix to invert
    mat = np.empty([N,N],dtype=np.longdouble)
    m1 = dt1.copy(); m2 = dt2.copy()
    for i in xrange(N):
        mat[:,i] = m2
        tterm = self.get('F%d'%i,type=np.longdouble)
        phi += tterm * m1
        m1 *= dt1/(i+2)
        m2 *= dt2/(i+2)

    return mat,phi


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
