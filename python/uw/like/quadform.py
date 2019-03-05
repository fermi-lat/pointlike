""" 
Manage quadratic fitting, elliptical representation of the TS surface

see Localize

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/quadform.py,v 1.2 2010/08/03 03:48:03 lande Exp $

"""
from numpy import asarray, array, arange, zeros, matrix, sign, rec, linspace, isnan
from math import sqrt, atan2, sin, cos, atan, pi, degrees, radians
import os
from skymaps import SkyDir
import numpy as np


def quadfun(r, p):
    x,y=r
    return p[0]*x*x+p[1]*x+p[2]*y*y+p[3]*y+p[4]*x*y+p[5]

flip= True  # set to measure from N
#major = True # set to make semi-major axis first

class QuadForm(object):
    d = 1/sqrt(2)
    # points at the center, then ccw from x axis on ring of unit radius
    points = [(0,0), (1,0), (d,d), (0,1), (-d,d), (-1,0), (-d,-d),  (0,-1),   (d,-d)]
    F = matrix([[x*x, x, y*y, y, x*y, 1] for x,y in points])
    Q = (F.T*F) **-1
    Q[ abs(Q)<1e-9] = 0

    def __init__(self, u):
        """ u: array of values evaluated at the ring of points """
        if len(u)==8: # assume only the ring, differences from center, so add a zero
            self.u = asarray([0]+list(u))
        else:
            self.u = asarray(u)
        self.p = array((self.u*QuadForm.F)*QuadForm.Q)[0] #magically does the least squares fit!
        self.v = array([self(r) for r in QuadForm.points]) #evaluate the fit function at the points
        self.chisq = ((self.u -self.v)**2).sum() # can be used to check

    def __call__(self, r):
        return quadfun(r, self.p)

    def fun(self, f):
       return array([f((x,y),self.p) for x,y in QuadForm.points])


class Ellipse(object):
    def __init__(self, q, convert=False, raw=True):
        """ q: ellipical parameters
            convert [False] if True, convert from quadratic form
            raw [False]   if True, the array is data to be fit

        """
        if raw:
            self.qf = QuadForm(q)
            self.chisq = self.qf.chisq
            self.q = self.convert(self.qf.p)
        else:
            self.q = q if not convert else self.convert(q)
            self.chisq=-1 #NA

    def convert(self, p):
        """ convert to elliptical parameters
            input: p defined as,
                f((x,y),p)= p[0]*x**2+p[1]*x+p[2]*y**2+p[3]*y+p[4]*x*y+p[5]
            output:
                a,b semimajor,-minor axis
                phi angle to semimajor axis (rad)
                x0,y0 origin
                k  constant offset
        """
        p = asarray(p,dtype='float64') # int can mess up angle
        if 4*p[0]*p[2] <= p[4]**2:
                raise Exception('Poorly formed quadratic form')
        if p[0] <0 and p[2] < 0:
                p = -p  # following assumes positive
        X = p[2]-p[0]
        Y = p[4]
        if abs(X)<1e-9: X=0
        Z = p[0]+p[2]
        if abs(Y)<1e-9: Y=0
        D = sqrt(X*X+Y*Y)
        if X<0: D=-D
        phi = pi/4 if X==0 else atan(Y/X)/2.
        phi = phi-pi/2
        alpha, beta = abs(Z-D)/2, abs(Z+D)/2
        a = 1/sqrt(alpha)
        b = 1/sqrt(beta)
        # convert to major axis, map to range (-90 to 90 deg)
        if  b>a:
            a,b = b,a
            phi += pi/2
        if phi>= pi/2: phi-=pi
        if phi< -pi/2: phi+=pi

        det = 4*p[0]*p[2] - p[4]**2
        x0 = (p[4]*p[3]-2*p[2]*p[1])/det
        y0 = (p[4]*p[1]-2*p[0]*p[3])/det
        # constant term is messy
        s=-sin(phi); c = cos(phi)
        s,c = c,s
        K = x0**2*((c/a)**2+(s/b)**2)\
           +y0**2*((c/b)**2+(s/a)**2)\
           -2*c*s*x0*y0*(1/a**2-1/b**2)

        return [a,b,phi, x0,y0, p[5]-K ]

    def __call__(self,r):
        """ quadratic function with angular parameters 
            r is x,y, angular offset from current orign in ra,dec direction
            q: a, b, phi, x0, y0, k (elliptical parameters)
        """
        x,y = r
        phi = self.q[2]
        s,c = sin(-self.q[2]), cos(self.q[2])
        s,c = c,s
        dx,dy = x-self.q[3], y-self.q[4]
        return ((c*dx-s*dy)/self.q[0])**2\
             + ((c*dy+s*dx)/self.q[1])**2\
             + self.q[5]

    def contour(self, r=1, count=50):
        """ return set of points in around closed figure, offset from x0,yo"""
        s,c = sin(-self.q[2]), cos(self.q[2])
        a,b = self.q[0],self.q[1]
        x0,y0 = self.q[3],self.q[4]
        s,c = c,s
        x = []
        y = []
        for t in linspace(0, 2*pi, count):
            ct,st = cos(t), sin(t)
            x.append( r*(a*ct*c - b*st*s))
            y.append( r*(a*ct*s + b*st*c))
        return x,y      

    def draw(self, data=None, scale=2):
        import pylab 
        x,y = self.circuit()
        pylab.plot(x,y, '-')
        pylab.plot([x0],[y0], '+') 


        pylab.axis((-scale,scale,-scale,scale))
        pylab.axvline(0, color='k')
        pylab.axhline(0, color='k')
        pylab.grid()

            
def testit(p=[ 1, 0, 2., 0, 0, 0 ]):
    print 'testing with quad pars=' ,p
    points =  QuadForm.points
    u1 = asarray([quadfun(r,p) for r in points]) # generate data
    qf = QuadForm(u1)
    pfit = qf.p   # fit quadratic form coefficients
    check = ((p-pfit)**2).sum()
    print 'fit chisq, check: %10.1g %10.1g' %  (qf.chisq, check)
    if check>1e-9: 
        print 'failed to fit: output parameters ', pfit
    ell = Ellipse(pfit,True,raw=False)
    print ('elliptical pars: '+6*'%10.3f') % tuple(ell.q)
    u2= asarray([ell(r) for r in points]) 
    check = ((u1-u2)**2).sum()
    print 'compare two functions: %10.1g' %check
    if check>1e-10:
        print 'Failed comparison!'
        print (5*'%+10s') % ('x   ','y   ', 'quad ', 'eliptical','diff ')
        for i, (x,y) in enumerate(points):
            print (4*'%10.3f'+'%10.1g') % ( x,y , u1[i], u2[i], u1[i]-u2[i] )

class Localize(object):
    fit_radius=2.5 #### modified from 2, works better for many weak sources
    def __init__(self, psl, verbose=True):
        self.verbose = verbose
        self.psl = psl
        self.dir = psl.dir()
        self.ra,self.dec = self.dir.ra(),self.dir.dec()
        self.sigma = psl.errorCircle()
        self.qual_cache=-1
        if verbose: print ('initial: ra,dec, sigma:' +3*'%10.4f') % (self.ra,self.dec,self.sigma)

        #self.fit(update=True)
        
        try:
            self.fit(update=True) # needed?
        except:
            if self.verbose: print 'update failed: center on highest TS and try again'
            self.recenter()
            self.fit(update=True)
            
    def recenter(self):
        ts = array(self.ts)
        tsmax = ts.max()
        if isnan(tsmax):
            print 'really lost'
            raise Exception('Localize: reallylost')
        idir = arange(9)[tsmax==ts][0]
        mdir = self.rcirc[idir]
        if self.verbose: print 'try ra,dec,ts =', mdir.ra(), mdir.dec(), tsmax 
        self.ra= mdir.ra()
        self.dec = mdir.dec()

        
    def TS(self,sdir):
        return self.psl.TSmap(sdir)
    
    def fit(self, update=True):
        verbose = self.verbose
        self.rcirc = self.circle()
        self.qual_cache = -1
        self.ts = [self.TS(r) for r in self.rcirc]
        if verbose: print ('ts:   ' + ' '.join(9*['%9.2f'])) % tuple(self.ts)
        self.ellipse = Ellipse(self.ts)
        self.chisq = self.ellipse.chisq
        if verbose: print ('resid:' + ' '.join(9*['%9.2f']))% tuple(self.ts-self.ellipse.qf.v)
        if verbose: print ('fit:  ' +len(self.ellipse.q)*'%9.2f') % tuple(self.ellipse.q)
        if verbose: print 'chisq: %9.2f' % self.ellipse.chisq
        radius = Localize.fit_radius
        if update:
            self.ra += self.ellipse.q[3]*self.sigma*radius
            self.dec+= self.ellipse.q[4]*self.sigma*radius

            self.dir = SkyDir(self.ra,self.dec)
            self.sigma= sqrt(self.ellipse.q[0]*self.ellipse.q[1])*self.sigma*radius
            if verbose: print ('update:  ra,dec, sigma:' +3*'%10.4f') % (self.ra,self.dec,self.sigma)

        self.par=[self.ra, self.dec, self.ts[0], radius*self.sigma*self.ellipse.q[0],
                radius*self.sigma*self.ellipse.q[1], degrees(self.ellipse.q[2]),
                self.quality(),  # insert quality, but leave after
                self.ellipse.chisq, # this was interpreted as a quality
                ]

    def circle(self):
        """ make a circle at the radius (in sigma units) """
        d = 1/sqrt(2.)
        points = [(0,0), (1,0), (d,d), (0,1), (-d,d), (-1,0), (-d,-d),  (0,-1),   (d,-d)]
        ddec = Localize.fit_radius*self.sigma
        dra  = ddec/cos(radians(self.dec))
        return [SkyDir(self.ra+x*dra,  self.dec+y*ddec) for x,y in points]

    def quality(self, radius=2.5):
        """ return a quality factor for the fit,
            the sqrt of the sum of squares of the residuals at the given radius
        """
        if self.qual_cache>0: return self.qual_cache
        qf = self
        xp,yp = qf.ellipse.contour(qf.fit_radius, 8)  # get points at standard radius
        ddec = radius*qf.sigma
        dra  = ddec/np.cos(np.radians(qf.dec))
        points = [SkyDir(qf.ra-x*dra,  qf.dec+y*ddec) for x,y in zip(xp,yp)]

        tszero = qf.TS(SkyDir(qf.ra,qf.dec))-radius**2;
        ts = np.asarray([qf.TS(p) for p in points]) #evaluate TS at the points
        qual = np.sqrt( ((ts-tszero)**2).sum())
        self.qual_cache=qual
        return qual





if __name__=='__main__':
    testit([ 1, 0, 2., 0, 1, 0 ])
    testit([ 2, 0, 1., 0, 1, 0 ])
    testit([1,1,1,1,1,1])
    testit([9.,0,9.,0,0,0]) # should give sigmas of 1/3
