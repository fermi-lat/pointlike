import numpy as N
import pyfits as PF
from types import ListType,FunctionType,MethodType

try:
   import uw.pointlike
except: pass


import pointlike as pl
from pointlike import SkyDir

def phase_circ(eventfiles,center=None,radius=6,phaseranges=[[0,1]],\
   erange=None,event_class=-1,ctbclasslevel=None,drift_check=False,phase_col_name='PULSE_PHASE'):
    """Construct a histogram of events interior to a circular region."""

    if not (type(radius) is FunctionType or type(radius) is MethodType):
        rval=radius
        def radius(e,event_class=0): return e/e*rval

    if not type(eventfiles) is ListType: eventfiles=[eventfiles]

    center=center or SkyDir(128.83646,-45.17658) #Vela

    ef = [PF.open(e,memmap=1) for e in eventfiles]
    columns = {}
    columns['ra']  = N.concatenate([N.asarray(e['EVENTS'].data.field('RA')).astype(float) for e in ef])
    columns['dec'] = N.concatenate([N.asarray(e['EVENTS'].data.field('DEC')).astype(float) for e in ef])
    columns['ph']  = N.concatenate([N.asarray(e['EVENTS'].data.field(phase_col_name)).astype(float) for e in ef])
    columns['en']  = N.concatenate([N.asarray(e['EVENTS'].data.field('ENERGY')).astype(float) for e in ef])
    columns['ec']  = N.concatenate([N.asarray(e['EVENTS'].data.field('EVENT_CLASS')).astype(int) for e in ef])
    if ctbclasslevel is not None:
      columns['ctb'] = N.concatenate([N.asarray(e['EVENTS'].data.field('CTBCLASSLEVEL')).astype(int) for e in ef])
    if drift_check:
      columns['dc'] = N.concatenate([N.asarray(e['EVENTS'].data.field('TIME')).astype(float) for e in ef])
    for e in ef: e.close()

    columns['radii'] = radius(columns['en'],event_class=columns['ec'])
    #Fast, quick and dirty cut
    mask = N.abs(columns['dec']-center.dec())<=columns['radii']
    for x in ['ra','dec','ph','en','ec','radii','ctb','dc']:
      if x in columns: columns[x] = columns[x][mask]


    #Slow, accurate cut
    scale=180/N.pi
    ra,dec,radii = columns['ra'],columns['dec'],columns['radii']
    mask = N.array([center.difference(SkyDir(ra[i],dec[i]))*scale<radii[i] for i in xrange(len(ra))])
    for x in ['ra','dec','ph','en','ec','radii','ctb','dc']:
      if x in columns: columns[x] = columns[x][mask]

    ph = columns['ph']
    mask=N.array([False]*len(ph))
    for r in phaseranges:
        for i in xrange(len(ph)):
            if r[0]<= ph[i] and ph[i] < r[1]: mask[i]=True
    if erange is not None:
        mask = N.logical_and( N.logical_and(columns['en']>=erange[0],columns['en']<erange[1]), mask)
    if event_class >= 0:
        mask = N.logical_and(columns['ec']==event_class,mask)
    if ctbclasslevel is not None:
        mask = N.logical_and(columns['ctb']>=ctbclasslevel,mask)
    
    if drift_check:
        return [ph[mask],columns['dc'][mask]] #Return phase and time
    return ph[mask] #Just return phase

def phase_cut(eventfile,outputfile=None,phaseranges=[[0,1]],phase_col_name='PULSE_PHASE'):
    """Select phases within a set of intervals.
    
        outputfile  - set to change the default output (eventfile_PHASECUT.fits)
        phaseranges - a set of ranges on which to make inclusive cuts"""


    from numarray import array as narray #important -- on SLAC, this should be changed to numpy.array
    from numpy import array
    import pyfits as PF
    ef=PF.open(eventfile)
    ph=array(ef['EVENTS'].data.field(phase_col_name)).astype(float)
    mask=narray([False]*len(ph))
    for r in phaseranges:
        for i in xrange(len(ph)):
            if r[0]<= ph[i] and ph[i] <= r[1]: mask[i]=True                    
    hdu=PF.new_table(ef['EVENTS'].columns,nrows=len(mask[mask]))
    for i in xrange(len(ef['EVENTS'].columns)):
        hdu.data.field(i)[:]=ef['EVENTS'].data.field(i)[mask]
    ef['EVENTS'].data=hdu.data
    if outputfile:
        ef.writeto(outputfile,clobber=True)
    else:
        ef.writeto(eventfile.replace('.fits','_PHASECUT.fits'),clobber=True)
    ef.close()

def constant_count_histogram(phases,photons_per_bin=100):
   """Return a set of bins and rates for a 'constant count' lightcurve.
      phases -- a list of photon phases"""
   from collections import deque
   counter = 0
   ev = N.asarray(phases)
   ev.sort()

   edges = deque()
   edges.append(0)
   for i,j in enumerate(ev):
      if counter == photons_per_bin-1:
         edges.append(j)
         counter = 0
      else: counter += 1
   edges.append(1)
   edges = N.asarray(edges)
   delta_phi = edges[1:]-edges[:-1]
   rates = photons_per_bin/delta_phi
   rates[-1] = counter/delta_phi[-1]
   return [edges[:-1],rates,delta_phi] #Edges, rates, and widths, appropriate for a bar chart

class PulsarLightCurve(object):

   def init(self):
      self.cookie_cutter_radius = None
      self.energy_range         = [100,2e4]
      self.phase_col_name       = 'PULSE_PHASE'
      self.use_mc_psf           = True
      self.max_radius           = 20.

   def __init__(self,**kwargs):
      self.init()
      self.__dict__.update(kwargs)

   def get_phases(self,event_files,skydir):

      if self.cookie_cutter_radius is not None:
         self.phases,self.times = phase_circ(event_files,center=skydir,erange=self.energy_range,\
                                  phase_col_name=self.phase_col_name, drift_check = True, radius = self.cookie_cutter_radius)
      else:
         from psf import PSF         
         psf_obj = PSF(use_mc=self.use_mc_psf)

         b = N.logspace(N.log10(self.energy_range[0]),N.log10(self.energy_range[1]))
         bin_cents = (b[:-1]*b[1:])**0.5
      
         radii = N.append([psf_obj.conf_region(e,event_class=0) for e in bin_cents],\
                 [psf_obj.conf_region(e,event_class=1) for e in bin_cents])

         def r(e,event_class=0):
            rads = N.zeros_like(e)
            for i in xrange(len(e)):
               for j in xrange(len(b)-1):
                  if b[j]<= e[i] and e[i] < b[j+1]:
                     rads[i] = radii[j*(1+event_class[i])]
                     break
            rads = N.minimum(self.max_radius,rads)
            return rads
         
         self.phases,self.times = phase_circ(event_files,center=skydir,erange=self.energy_range,\
                                  phase_col_name=self.phase_col_name, drift_check = True, radius = r)

   def plot(self,fignum=1,show_drift=True,show_errorbars=False,outfile=None,show=True):
      import pylab as P
      if show_drift: P.figure(fignum,(6,10))
      else: P.figure(fignum)

      npb = max(int(round(len(self.phases)/100.)),20)

      ev_cop = self.phases.copy()
      times = (self.times-self.times.min())/(24*3600.)
      phases,rates,delta_phi = constant_count_histogram(self.phases,npb)
      rates_pre = rates.copy()
      rates = rates/rates[rates>0].min()

      ax1 = P.axes([0.15,0.1,0.75,0.4]) if show_drift else P.axes()
      ax1.step(N.append(phases,1),N.append(rates[0],rates),where='pre',color='black',label='%d cts/bin'%npb)
      ax1.set_ylabel('Relative Count Rate')
      ax1.grid(b=True)
      if show_errorbars:
         tph = N.append(phases,1)
         delta_phi = tph[1:]-tph[:-1]
         errs = ((rates_pre*delta_phi)**0.5/rates_pre[rates_pre>0].min())/delta_phi
         tph = (tph[:-1]+tph[1:])/2.
         P.errorbar(tph,rates,yerr=errs,ls=' ',capsize=0)

      if show_drift:
         ax2 = P.axes([0.15,0.55,0.75,0.4])
         stride = max(1,int(len(ev_cop)/500.))
         ax2.scatter(ev_cop[::stride],times[::stride],marker='+',color='red')
         a = ax2.axis()
         ax2.axis([0,1,0,a[3]])
         ax2.grid(b=True)
         ax2.set_ylabel('Days')

      P.xlabel('Phase')
      if outfile is not None: P.savefig(outfile)
      if show: P.show()


#--------------------------------------------------------
    
if __name__=='__main__':
   a=Phase('ft1_first_diff_v2.fits')
   a.write()
   from pylab import hist
   hist(a.phaseHist(), bins = 20)
    

        
"""Same, but with periods.  This is harder, don't know why I bothered.
if pd==0:
    phases = [t/p for t in self.times] #Exact to 0th
elif pdd==0:
    phases = [1/pd*N.log(1.+pd/p*t) for t in self.times] #Exact to 1st
else:
    d=(4*pdd*p-pd**2)**0.5 #Discriminant
    phases = [2/d*(N.arctan( (pd+2*pdd*t)/d ) - N.arctan( pd/d )) for t in self.times] #Exact to 2nd
"""