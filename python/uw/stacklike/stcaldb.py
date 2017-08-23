import astropy.io.fits as pf
import numpy as np
import os as os
import scipy.integrate as si
from uw.stacklike.angularmodels import PSF

class IrfLoader(object):
    
    def __init__(self,name):
        self.cdir = os.environ['CALDB']+'/data/glast/lat/bcf/'
        self.fpsffile = self.cdir+'psf/psf_%s_front.fits'%name
        self.bpsffile = self.cdir+'psf/psf_%s_back.fits'%name
        self.feafile = self.cdir+'ea/aeff_%s_front.fits'%name
        self.beafile = self.cdir+'ea/aeff_%s_back.fits'%name
        self.pars = [self.load(self.fpsffile,self.feafile,0),self.load(self.bpsffile,self.beafile,1)]
        self.dims = [len(self.pars[0][0]),len(self.pars[0][2])]
        self.effdims = [len(self.pars[0][-6]),len(self.pars[0][-4])]

    def load(self,psfname,eaname,eclass):
        ff = pf.open(psfname)
        tb = ff[1].data
        emin = tb.field('ENERG_LO')[0]
        emax = tb.field('ENERG_HI')[0]
        ctlo = tb.field('CTHETA_LO')[0]
        cthi = tb.field('CTHETA_HI')[0]
        try:
            s1 = tb.field('SCORE')[0]
            s2 = tb.field('STAIL')[0]
            self.old=False
        except:
            s1 = tb.field('SIGMA')[0]
            s2 = tb.field('SIGMA')[0]
            self.old=True
        g1 = tb.field('GCORE')[0]
        g2 = tb.field('GTAIL')[0]
        try:
            n1 = tb.field('NCORE')[0]
            n2 = tb.field('NTAIL')[0]
        except:
            n1 = tb.field('NCORE')[0]
            n2 = 1-tb.field('NCORE')[0]
        p0 = ff[2].data[0][0]
        ff.close()
        ff = pf.open(eaname)
        tb = ff[1].data
        effemi = tb.field('ENERG_LO')[0]
        effema = tb.field('ENERG_HI')[0]
        effctl = tb.field('CTHETA_LO')[0]
        effcth = tb.field('CTHETA_HI')[0]
        eff = tb.field('EFFAREA')[0]
        ff.close()
        #scales = np.array([uf.scale2(np.sqrt(emin[it]*emax[it]),p0[2*eclass],p0[2*eclass+1],p0[4]) for it2 in range(len(ctlo)) for it in range(len(emin))])
        return [emin,emax,ctlo,cthi,s1,g1,n1,s2,g2,n2,effemi,effema,effctl,effcth,eff,p0]

    def effaverage(self,emin,emax,ctmin,ctmax,eclass):
        emiidx = min(len(self.pars[eclass][-5][self.pars[eclass][-5]<emin]),self.effdims[0]-1)
        emaidx = min(len(self.pars[eclass][-5][self.pars[eclass][-5]<emax]),self.effdims[0]-1)
        ctmidx = len(self.pars[eclass][-3][self.pars[eclass][-3]<ctmin])
        ctmadx = len(self.pars[eclass][-3][self.pars[eclass][-3]<ctmax])
        erange = np.arange(emiidx,emaidx+1,1)
        crange = np.arange(ctmidx,ctmadx+1,1)
        effarea = np.array([self.pars[eclass][-2][it][it2] for it in crange for it2 in erange])
        return sum(effarea)/len(effarea)

    def params(self,energy,eclass,ct=-1):
        eidx = min(len(self.pars[eclass][1][self.pars[eclass][1]<energy]),self.dims[0]-1)
        effidx = min(len(self.pars[eclass][-5][self.pars[eclass][-5]<energy]),self.effdims[0]-1)
        emi = self.pars[eclass][0][eidx]
        ema = self.pars[eclass][1][eidx]
        #average pars with effarea
        if ct<0:
            eweights = np.array([self.effaverage(emi,ema,self.pars[eclass][2][it],self.pars[eclass][3][it],eclass) for it in range(self.dims[1])])
            effave = sum(eweights)/len(eweights)
            eweights/=sum(eweights)
            rpars = [sum(eweights*np.array([self.pars[eclass][it2+4][it][eidx] for it in range(self.dims[1])])) for it2 in range(6)]
            rpars.append(effave)
            scale = scale2(energy,self.pars[eclass][-1][2*eclass],self.pars[eclass][-1][2*eclass+1],self.pars[eclass][-1][4])
            rpars[0]*=scale
            rpars[3]*=scale
            if self.old:
                f0 = 2*np.sqrt(5)
                rpars[2]=1./(1+psfbase(f0*rpars[0],rpars[0],rpars[1])/psfbase(f0*rpars[0],rpars[3],rpars[4]))
            else:
                rpars[2]=1./(1+rpars[5]*(rpars[3]/rpars[0])**2)
            rpars[5]=1-rpars[2]
            return np.array(rpars)
        #find cos bin
        else:
            cidx = len(self.pars[eclass][3][self.pars[eclass][3]<ct])
            ceffidx = len(self.pars[eclass][-3][self.pars[eclass][-3]<ct])
            rpars = [self.pars[eclass][it2+4][self.dims[0]*cidx+eidx] for it2 in range(6)]
            rpars.append(self.pars[eclass][-2][self.effdims[0]*ceffidx+effidx])
            rpars = np.array(rpars)
            scale = scale2(energy,self.pars[eclass][-1][2*eclass],self.pars[eclass][-1][2*eclass+1],self.pars[eclass][-1][4])
            rpars[0]*=scale
            rpars[3]*=scale
            rpars[2]=1./(1+rpars[5]*(rpars[3]/rpars[0])**2)
            rpars[5]=1-rpars[2]
            return rpars

    def average_psf(self,emin,emax,ctmi,ctma,eclass,ltc=None,sd=None):
        eiidx = min(len(self.pars[eclass][1][self.pars[eclass][1]<emin]),self.dims[0]-1)
        eaidx = min(len(self.pars[eclass][1][self.pars[eclass][1]<emax]),self.dims[0]-1)
        ciidx = len(self.pars[eclass][3][self.pars[eclass][3]<ctmi])
        caidx = len(self.pars[eclass][3][self.pars[eclass][3]<ctma])
        era = np.arange(eiidx,eaidx+1,1)
        cra = np.arange(ciidx,caidx+1,1)
        eweights = np.array([self.pars[eclass][-2][it][it2]/(self.pars[eclass][1][it2]*self.pars[eclass][0][it2]) for it in cra for it2 in era])
        if ltc is not None:
            lweights = np.array([ltc.value(sd,float(min(self.pars[eclass][3][it],1.-(1e-9)))) for it in cra for it2 in era])
            eweights*=lweights
        effave = sum(eweights)/len(eweights)
        eweights/=sum(eweights)
        rpars = [sum(eweights*np.array([self.pars[eclass][it3+4][self.dims[0]*it+it2] for it in cra for it2 in era])) for it3 in range(6)]
        scale = scale2(np.sqrt(emin*emax),self.pars[eclass][-1][2*eclass],self.pars[eclass][-1][2*eclass+1],self.pars[eclass][-1][4])
        rpars[0]*=scale
        rpars[3]*=scale
        if self.old:
            f0 = 2*np.sqrt(5)
            rpars[2]=1./(1+psfbase(f0*rpars[0],rpars[0],rpars[1])/psfbase(f0*rpars[0],rpars[3],rpars[4]))
        else:
            rpars[2]=1./(1+rpars[5]*(rpars[3]/rpars[0])**2)
        rpars[5]=1-rpars[2]
        return np.array(rpars)

    def rcontain(self,pars,frac):
        import scipy.optimize as so
        guess = (pars[0]*pars[2]+pars[3]*pars[5])/(2.*(1.-frac))
        fint = doublepsfint(np.pi,pars)
        bestp = so.fmin_powell(lambda x: ((doublepsfint(x,pars)/fint-frac))**2 if x>0 else np.Infinity,guess,disp=0,full_output=1)
        #print bestp[1]
        return bestp[0].item()

## 3 parameter scale function
#  @param e energy in MeV
#  @param c0 sigma at 100 MeV
#  @param c1 sigma from instrument pitch
#  @param b power law of multiple scattering
def scale2(e,c0,c1,b):
    return np.sqrt((c0*((e/100.)**b))**2+c1**2)

def doublepsfint(dm,pars):
    s1,g1,n1,s2,g2,n2 = pars
    return n1*intg(dm,s1,g1)+n2*intg(dm,s2,g2)

def intg(dm,s1,g1):
    if dm>0.5:
        fint = intgbase(0.1,s1,g1)
        fint += si.quad(lambda x: psfbase(x,s1,g1)*np.sin(x),0.1,dm)[0]
        #print fint/intgbase(dm,s1,g1)
        return fint
    else:
        return intgbase(dm,s1,g1)

def intgbase(dm,s1,g1):
    return (1.-(1+((dm/s1)**2)/(2*g1))**(-g1+1))*s1**2

def psfbase(dm,s1,g1):
    return (1-1./g1)*(1+((dm/s1)**2)/(2*g1))**(-g1)