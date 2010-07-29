"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/CLHEP.py,v 1.0 2010/07/29 13:53:17 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

import skymaps as s
import numpy as N

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
            self.matrix = N.matrix([[vecx.x(),vecy.x(),vecz.x()],[vecx.y(),vecy.y(),vecz.y()],[vecx.z(),vecy.z(),vecz.z()]])
        if len(c)==3 and not axes:
            x,y,z=c[0],c[1],c[2]
            #non-Euler
            self.matrix = N.matrix([[1,0,0],[0,N.cos(x),N.sin(x)],[0,-N.sin(x),N.cos(x)]])
            self.matrix = self.matrix * N.matrix([[N.cos(y),0,-N.sin(y)],[0,1,0],[N.sin(y),0,N.cos(y)]])
            self.matrix = self.matrix * N.matrix([[N.cos(z),N.sin(z),0],[-N.sin(z),N.cos(z),0],[0,0,1]])
            #Euler
            """self.matrix = N.matrix([[N.cos(x),N.sin(x),0],[-N.sin(x),N.cos(x),0],[0,0,1]])
            self.matrix = self.matrix * N.matrix([[1,0,0],[0,N.cos(y),N.sin(y)],[0,-N.sin(y),N.cos(y)]])
            self.matrix = self.matrix * N.matrix([[N.cos(z),N.sin(z),0],[-N.sin(z),N.cos(z),0],[0,0,1]])"""

    ## displays matrix
    def echo(self):
        for i in range(3):
                print self.matrix[i,0],self.matrix[i,1],self.matrix[i,2]

    ## multiplication with either matrix or vector
    #  @param other matrix to be multiplied
    def m(self,other):
        if other.matrix.shape[1]==3:
            return HepRotation([self.matrix * other.matrix])
        else:
            return Hep3Vector([self.matrix * other.vec])

    ## make inverse
    def inverse(self):
        return HepRotation([N.linalg.inv(self.matrix)])

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
            self.vec = N.matrix([[N.cos(phi)*N.cos(theta)],[N.sin(phi)*N.cos(theta)],[N.sin(theta)]])
        if len(c)==3:
            x,y,z=c[0],c[1],c[2]
            self.vec = N.matrix([[x],[y],[z]])/N.sqrt(x*x+y*y+z*z)
        self.matrix =self.vec

    ## convert from a skydir
    #  @param sd skymaps SkyDir object
    def __call__(self,sd):
        phi,theta = sd.ra()*N.pi/180,sd.dec()*N.pi/180
        self.vec = N.matrix([[N.cos(phi)*N.cos(theta)],[N.sin(phi)*N.cos(theta)],[N.sin(theta)]])

    def phi(self):
        x,y=self.x(),self.y()
        if x>0 and y>0:
            return N.arctan(y/x)
        if x<0 and y>0:
            return N.pi-N.arctan(y/x)
        if x<0 and y<0:
            return N.pi+N.arctan(y/x)
        if x>0 and y<0:
            return 2*N.pi-N.arctan(y/x)
        if x==0 and y>0:
            return N.pi/2
        if x==0 and y<0:
            return 3.*N.pi/2
        if y==0:
            return 0

    def dir(self):
        return s.SkyDir(self.ra(),self.dec())
    
    def ra(self):
        return self.phi()*180/N.pi
    
    def dec(self):
        return self.theta()*180/N.pi
    
    def theta(self):
        return N.arcsin(self.z())

    def x(self):
        return self.vec[0][0].item()
    
    def y(self):
        return self.vec[1][0].item()
    
    def z(self):
        return self.vec[2][0].item()

    def cross(self,other):
        return Hep3Vector([N.matrix([[self.y()*other.z()-self.z()*other.y()],[self.z()*other.x()-self.x()*other.z()],[self.x()*other.y()-self.y()*other.x()]])])

    def dot(self,other):
        return self.x()*other.x()+self.y()*other.y()+self.z()*other.z()

    def echo(self):
        for i in range(3):
            print self.vec[i,0]

################################################# END HEP3VECTOR CLASS #######################################################


################################################# START PHOTON CLASS #########################################################

## "Smart photon": associated with a source, and can be rotated in the GLAST frame

class Photon(object):

    def __init__(self,ra,dec,en,time,ec,hepx,hepz,srcdir):
        self.vec=Hep3Vector([ra*N.pi/180,dec*N.pi/180])
        self.energy=en
        self.time=time
        self.event_class=ec
        self.hepx=hepx
        self.hepz=hepz
        self.hepy=hepz.cross(hepx)
        self.rot=HepRotation([hepx,self.hepy,hepz])
        self.srcdir=srcdir
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        self.sdiff = N.arccos(h1.dot(self.vec))

    ## return rotated vector
    #  @param rot CLHEP HepRotation matrix
    def rotate(self,rot):
        return self.rot.m(rot.m(self.rot.inverse().m(self.vec)))

    ## return angular separation from source (after rotation in GLAST frame)
    #  @param rot CLHEP HepRotation matrix
    def diff(self,rot):
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        return N.arccos(h1.dot(self.rotate(rot)))
    
    ## return angular separation from source
    def srcdiff(self):
        return self.sdiff

###################################################  END PHOTON CLASS #######################################################