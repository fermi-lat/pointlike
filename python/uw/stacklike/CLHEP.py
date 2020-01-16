"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/stacklike/CLHEP.py,v 1.6 2012/03/12 20:31:59 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

import numpy as np

############################################ START HEPROTATION CLASS ##########################################################

## Realization of the rotation matrix in CLHEP (takes either three angles, three vectors, or another matrix)
class HepRotation(object):

    ## Constructor
    #
    # @param c HepRotation, [Hep3Vectors] , [rotation angles in radians]
    # @param axes if a list, Hep3Vectors?
    def __init__(self,c,axes=True):
        if len(c)==1:
            self.matrix = c[0]
        if len(c)==3 and axes:
            vecx,vecy,vecz=c[0],c[1],c[2]
            self.matrix = np.matrix([[vecx.x(),vecy.x(),vecz.x()],[vecx.y(),vecy.y(),vecz.y()],[vecx.z(),vecy.z(),vecz.z()]])
        if len(c)==3 and not axes:
            x,y,z=c[0],c[1],c[2]
            #non-Euler
            self.matrix = np.matrix([[1,0,0],[0,np.cos(x),np.sin(x)],[0,-np.sin(x),np.cos(x)]])
            self.matrix = self.matrix * np.matrix([[np.cos(y),0,-np.sin(y)],[0,1,0],[np.sin(y),0,np.cos(y)]])
            self.matrix = self.matrix * np.matrix([[np.cos(z),np.sin(z),0],[-np.sin(z),np.cos(z),0],[0,0,1]])
            #Euler
            """self.matrix = np.matrix([[np.cos(x),np.sin(x),0],[-np.sin(x),np.cos(x),0],[0,0,1]])
            self.matrix = self.matrix * np.matrix([[1,0,0],[0,np.cos(y),np.sin(y)],[0,-np.sin(y),np.cos(y)]])
            self.matrix = self.matrix * np.matrix([[np.cos(z),np.sin(z),0],[-np.sin(z),np.cos(z),0],[0,0,1]])"""

    ## displays matrix
    def echo(self):
        for i in range(3):
                print (self.matrix[i,0],self.matrix[i,1],self.matrix[i,2])

    ## multiplication with either matrix or vector
    #  @param other matrix to be multiplied
    def m(self,other):
        if other.matrix.shape[1]==3:
            return HepRotation([self.matrix * other.matrix])
        else:
            return Hep3Vector([self.matrix * other.vec])

    ## make inverse
    def inverse(self):
        return HepRotation([np.linalg.inv(self.matrix)])

################################################## END HEPROTATION CLASS ################################################


################################################## START HEP3VECTOR CLASS ###############################################

## Realization of the Vector class in CLHEP with SkyDir support

class Hep3Vector(object):

    def __init__(self,c):
        if len(c)==0:
            self.vec=c
        if len(c)==1:
            self.vec = c[0]
        if len(c)==2:
            phi,theta = c[0],c[1]
            self.vec = np.matrix([[np.cos(phi)*np.cos(theta)],[np.sin(phi)*np.cos(theta)],[np.sin(theta)]])
        if len(c)==3:
            x,y,z=c[0],c[1],c[2]
            self.vec = np.matrix([[x],[y],[z]])/np.sqrt(x*x+y*y+z*z)
        self.matrix =self.vec

    ## convert from a skydir
    #  @param sd skymaps SkyDir object
    def __call__(self,sd):
        phi,theta = sd.ra()*np.pi/180,sd.dec()*np.pi/180
        self.vec = np.matrix([[np.cos(phi)*np.cos(theta)],[np.sin(phi)*np.cos(theta)],[np.sin(theta)]])

    def __str__(self):
        return 'vec: (%1.4f, %1.4f, %1.4f)'%(self.x(),self.y(),self.z())

    def add(self,other):
        return Hep3Vector([self.x()+other.x(),self.y()+other.y(),self.z()+other.z()])

    def subt(self,other):
        return Hep3Vector([self.x()-other.x(),self.y()-other.y(),self.z()-other.z()])

    def phi(self):
        x,y=self.x(),self.y()
        if x>0 and y>0:
            return np.arctan(y/x)
        if x<0 and y>0:
            return np.pi+np.arctan(y/x)
        if x<0 and y<0:
            return np.pi+np.arctan(y/x)
        if x>0 and y<0:
            return 2*np.pi+np.arctan(y/x)
        if x==0 and y>0:
            return np.pi/2
        if x==0 and y<0:
            return 3.*np.pi/2
        if y==0 and x>0:
            return 0
        if y==0:
            return np.pi
    
    def ra(self):
        return self.phi()*180/np.pi
    
    def dec(self):
        return self.theta()*180/np.pi
    
    def theta(self):
        return np.arcsin(self.z())

    def x(self):
        return self.vec[0][0].item()
    
    def y(self):
        return self.vec[1][0].item()
    
    def z(self):
        return self.vec[2][0].item()

    def cross(self,other):
        return Hep3Vector([np.matrix([[self.y()*other.z()-self.z()*other.y()],[self.z()*other.x()-self.x()*other.z()],[self.x()*other.y()-self.y()*other.x()]])])

    def dot(self,other):
        return self.x()*other.x()+self.y()*other.y()+self.z()*other.z()
 
    def diff(self,other):
        return np.arccos(self.dot(other))

    def norm(self):
        x,y,z=self.x(),self.y(),self.z()
        self.vec = np.matrix([[x],[y],[z]])/np.sqrt(x*x+y*y+z*z)
    
    def scale(self,scl):
        x,y,z=self.x(),self.y(),self.z()
        return Hep3Vector([np.matrix([[x*scl],[y*scl],[z*scl]])])
 
    def echo(self):
        for i in range(3):
            print (self.vec[i,0])

################################################# END HEP3VECTOR CLASS #######################################################


################################################# START PHOTON CLASS #########################################################

## "Smart photon": associated with a source, and can be rotated in the GLAST frame

class Photon(object):

    def __init__(self,ra,dec,en,time,ec,hepx,hepz,srcdir,weight=1,ct=1.0):
        self.vec=Hep3Vector([ra*np.pi/180,dec*np.pi/180])
        self.energy=en
        self.time=time
        self.event_class=ec
        if hepx!=[]:
            self.hepx=hepx
            self.hepz=hepz
            self.hepy=hepz.cross(hepx)
            self.rot=HepRotation([hepx,self.hepy,hepz])
        self.srcdir=srcdir
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        self.sdiff = np.arccos(h1.dot(self.vec))
        self.weight=weight
        self.ct=ct

    ## return rotated vector
    #  @param rot CLHEP HepRotation matrix
    def rotate(self,rot):
        return self.rot.m(rot.m(self.rot.inverse().m(self.vec)))

    ## return angular separation from source (after rotation in GLAST frame)
    #  @param rot CLHEP HepRotation matrix
    def diff(self,rot):
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        return np.arccos(h1.dot(self.rotate(rot)))
    
    def difftheta(self,dt,verb=False):
        from time import sleep
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        sc = self.rot.inverse().m(self.vec)
        ct = sc.z()
        st = np.sqrt(1-ct*ct)
        cp = sc.x()/st
        sp = sc.y()/st
        cdt = np.cos(dt)
        sdt = np.sin(dt)
        nw = Hep3Vector([cp*(st*cdt+ct*sdt),sp*(st*cdt+ct*sdt),ct*cdt-st*sdt])
        diff = np.arccos(h1.dot(self.rot.m(nw)))
        if verb:
            print (sc.x(),sc.y(),sc.z())
            print (nw.x(),nw.y(),nw.z())
            print (ct,st,cp,sp,cdt,sdt,diff)
        return diff

    ## return angular separation from source
    def srcdiff(self):
        return self.sdiff

###################################################  END PHOTON CLASS #######################################################