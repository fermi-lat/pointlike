import pyfits
import sys
from numpy import *
import pointlike as pl

class PeakFinder:
 
   def __init__(self,mints=25,mindist=0.2,highE_startbin=4):
      self.mindist=mindist
      self.mints=mints
      self.hE_start=highE_startbin
      self.projector=0
      self.nx=0
      self.ny=0
      self.peaks=[]
      self.pixels=[]
      
   def findPeaks(self,tsmap,verbose=0):   
       self.projector=pl.SkyProj(tsmap)
       f=pyfits.open(tsmap)
       pic=f["PRIMARY"].data
       head=f["PRIMARY"].header
       csys,proj=head["CTYPE1"].split('-')
       if csys=='GLON': regstring='galactic'
       elif csys=='RA': 
          regstring='fk5'
          raise RuntimeError,"Currently only galactic coordinates. Sorry."
       else: raise RuntimeError,"Coordinate system %s not understood."%(csys)

       self.le_ts=pic[1:,:,:].sum(0).transpose()
       self.he_ts=pic[self.hE_start:,:,:].sum(0).transpose()
       self.ts=self.le_ts
       self.nx,self.ny=self.le_ts.shape
       self.cl,self.cb=self.projector.pix2sph(self.nx/2+0.5,self.ny/2+0.5)
       nepix=[]
       for dx in range(-2,3):
	  for dy in range(-2,3):
	    if(dx!=0 or dy!=0): nepix.append((dx,dy))
       nepix=array(nepix,'Int32')     
       circ=[]

       peaks_hxy,peaks_hlb=[],[]
       peaks_lxy,peaks_llb=[],[]
       self.pixels=[]
       for x in range(2,self.nx-2):
	 for y in range(2,self.ny-2): 
	   xy=tile((x,y),24).reshape(-1,2)
	   xy=xy-nepix
	   l,b=self.projector.pix2sph(x+1.,y+1.)
	   
	   ts_h=self.he_ts[x,y]
	   ts_l=self.le_ts[x,y]
	   self.pixels.append((x,y,l,b,ts_l,ts_h,ts_l))
	   neighb_h=[self.he_ts[p[0],p[1]] for p in xy]
	   neighb_l=[self.le_ts[p[0],p[1]] for p in xy]
	   if(ts_h>self.mints and max(neighb_h)<ts_h and min(neighb_h)<0.9*ts_h ): 
	      if(verbose): print "Found peak (high): x=",x,"y=",y," l=",l," b=",b
	      peaks_hxy.append((x,y))
	      peaks_hlb.append((l,b))
	   if(ts_l>self.mints and max(neighb_l)<ts_l and min(neighb_l)<0.9*ts_l): 
	      if(verbose): print "Found peak (low): x=",x,"y=",y," l=",l," b=",b
	      peaks_lxy.append((x,y))
	      peaks_llb.append((l,b))

       self.peaks=[]
       for hl,hb in peaks_hlb:
           remove=False
	   p1dir=pl.SkyDir(hl,hb,pl.SkyDir.GALACTIC)
           for l,b in self.peaks:
   	      p2dir=pl.SkyDir(l,b,pl.SkyDir.GALACTIC)
	      if p1dir.difference(p2dir)<self.mindist*pi/180.: remove=True
           if(not remove): self.peaks.append((hl,hb))

       for ll,lb in peaks_llb:
	 remove=False
	 p1dir=pl.SkyDir(ll,lb,pl.SkyDir.GALACTIC)
	 for hl,hb in peaks_hlb:
 	     p2dir=pl.SkyDir(hl,hb,pl.SkyDir.GALACTIC)
	     if p1dir.difference(p2dir)<self.mindist*pi/180.: remove=True
	 if(not remove): self.peaks.append((ll,lb))

       self.peakts=[]
       for l,b in self.peaks:
          px,py=self.projector.sph2pix(l,b)
	  px,py=int(px-0.5),int(py-0.5)
          self.peakts.append(self.ts[px,py])
     
       return self.peaks

   def writeReg(self,outfile):
       out=open(outfile,'w')   
       out.write('global color=white point=diamond font="helvetica 11 normal" select=1 highlite=1 edit=1 move=0 delete=1 include=1 fixed=0 source\ngalactic\n\n') 
       for l,b in self.peaks:
          px,py=self.projector.sph2pix(l,b)
	  px,py=int(px-0.5),int(py-0.5)
	  out.write("circle(%.4f,%.4f,%.4f) # text={TS=%.1f}\n"%(float(l),float(b),0.1,float(self.ts[px,py])))
       out.close()


class ResultsFileParser:
   def __init__(self,resultfile,deltats=8.,skip=4):
      try:
        pfile=pyfits.open(resultfile)
        fname=pfile["SOURCES"].data.field('NAME')
      except AttributeError: return
      ftype=pfile["SOURCES"].data.field('TYPE')
      fl=pfile["SOURCES"].data.field('L')
      fb=pfile["SOURCES"].data.field('B')
      fts=pfile["SOURCES"].data.field('TS')
      fr=pfile["SOURCES"].data.field('R0') 
      
      self.srclist={}
      for name,stype,l,b,ts,r in zip(fname,ftype,fl,fb,fts,fr):
          basename=name.split(':')[0]
	  ts=ts[skip:].sum()
	  if basename not in self.srclist.keys() and stype in ['point','pseudopiont']:
	      self.srclist[basename]={'type':stype,'l':l,'b':b,'ts':ts,'r':0.,'ext':0,
	                              'dts':0.,'dts_disk':0.,'r_disk':0.,'ts_point':ts,
				      'l_point':l,'b_point':b}
	  
      for name,stype,l,b,ts,r in zip(fname,ftype,fl,fb,fts,fr):
         basename=name.split(':')[0]
         ts=ts[skip:].sum()
         if basename in self.srclist.keys() and stype not in ['point','pseudopoint']: 
	     if ts > self.srclist[basename]['ts']+deltats: self.srclist[basename]['ext']+=1

      for name,stype,l,b,ts,r in zip(fname,ftype,fl,fb,fts,fr):
         basename=name.split(':')[0]
         ts=ts[skip:].sum()
         if basename in self.srclist.keys() and self.srclist[basename]['ext']==2:
	     if stype=='disk':
	         self.srclist[basename]['dts_disk']=ts- self.srclist[basename]['ts_point']
	         self.srclist[basename]['r_disk']=r
	     if stype=='gauss':
	         self.srclist[basename]['dts']=ts- self.srclist[basename]['ts_point']
	         self.srclist[basename]['r']=r
	         self.srclist[basename]['l']=l
	         self.srclist[basename]['b']=b
	         self.srclist[basename]['ts']=ts
	         self.srclist[basename]['type']=stype
	     
 		 
   def sourceList(self):
       slist=[]
       for k in self.srclist.keys():
         elem=self.srclist[k] 
         slist.append((k,elem["l"],elem["b"],elem["type"],elem["ts"],elem["r"],elem['dts'],elem['dts_disk'],elem['l_point'],elem['b_point']))
       return slist		    
   
   def writeReg(self,outfile):
       out=open(outfile,'w')   
       out.write('global color=red point=diamond font="helvetica 11 normal" select=1 highlite=1 edit=1 move=0 delete=1 include=1 fixed=0 source\ngalactic\n\n') 
       for k in self.srclist.keys(): 
          elem=self.srclist[k] 
	  out.write("circle(%.4f,%.4f,%.4f) # text={%s(TS=%.1f)}\n"%(elem['l'],elem['b'],max(0.05,elem['r']*180./pi),k,elem['ts']))
       out.close()

if __name__ == "__main__":
    import ROOT as R

    p=ResultsFileParser(sys.argv[1])
    for name,l,b,stype,ts,rad,dts,dtsd,lp,bp in p.sourceList():
        if(ts>20): print "%20s %8s %8.3f %8.3f %8.1f %5.1f %6.1f %6.1f %8.3f %8.3f"%(name,stype,l,b,ts,rad*60.*180./pi,dts,dtsd,lp,bp)    
    p.writeReg(sys.argv[2])
    sys.exit(1)
    
    pf=PeakFinder()
    peaks= pf.findPeaks(sys.argv[1],1)
   
    h1=R.TH2F("heTS","TS (E>6.4GeV)",pf.nx,-0.5,pf.nx-0.5,pf.ny,-0.5,pf.ny-0.5)
    h2=R.TH2F("leTS","TS (E>400MeV)",pf.nx,-0.5,pf.nx-0.5,pf.ny,-0.5,pf.ny-0.5)

    for x,y,l,b,lts,hts,ts in pf.pixels:
       h1.Fill(x,y,hts)
       h2.Fill(x,y,lts)
       
    R.gStyle.SetPalette(1)

    cv1=R.TCanvas()
    h1.Draw("CONTZ")    

    circ=[]
    for l,b in peaks:
       px,py=pf.projector.sph2pix(l,b)
       circ.append(R.TEllipse(px-1,py-1,1.,1.))
       circ[-1].SetLineColor(1)
       circ[-1].SetLineWidth(2)
       circ[-1].Draw("same")


    cv2=R.TCanvas()
    h2.Draw("CONTZ")    

    for l,b in peaks:
       px,py=pf.projector.sph2pix(l,b)
       circ.append(R.TEllipse(px-1,py-1,1.,1.))
       circ[-1].SetLineColor(1)
       circ[-1].SetLineWidth(2)
       circ[-1].Draw("same")

    if(len(sys.argv)>2): pf.writeReg(sys.argv[2])   
    raw_input()

