"""
Plotting routines to display results of an ROI analysis.
Given an ROIAnalysis object roi:
    plot_spectra(roi)
    ROIDisplay(roi).show()
    plot_counts(roi)


$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_plotting.py,v 1.6 2010/02/23 19:13:57 burnett Exp $

author: Matthew Kerr
"""
import exceptions
import numpy as N
from skymaps import PySkyFunction,Background,Band,SkyDir,Hep3Vector,SkyIntegrator
from roi_bands import ROIEnergyBand
from pypsf import OldPsf
from uw.utilities.fitstools import get_fields
from uw.utilities import colormaps
from uw.utilities.image import ZEA

from collections import deque

from scipy.stats import poisson,norm
from scipy.optimize import fmin,fsolve

import pylab as P
from matplotlib import rcParams,mpl,pyplot,ticker
from matplotlib.patches import FancyArrow

def band_spectra(r,source=0):
   
   en = N.asarray([band.e for band in r.bands])

   groupings = [deque() for x in xrange(len(r.bin_centers))]

   #group slw by energy
   for i,ei in enumerate(r.bin_centers):
      for band in r.bands:
         if band.e == ei:
            groupings[i].append(band)
      groupings[i] = list(groupings[i])

   scales    = N.zeros(len(groupings))
   scales_hi = N.zeros(len(groupings))
   scales_lo = N.zeros(len(groupings))

   m      = r.psm.point_sources[source].model
   ps     = N.asarray([sum([band.ps_counts[source] for band in g]) for g in groupings])
   exp    = N.asarray([sum([band.expected(m)/m.i_flux(band.emin,band.emax) for band in g]) for g in groupings])

   for i,gi in enumerate(groupings):
      #print pslw.bin_centers[i]
      obs = sum([band.photons for band in gi])
      """ #way to handle 0 count bins
      
      if obs == 0:
         bg = sum([slw.gal_exp*slw.gal_counts + slw.iso_exp*slw.iso_counts for slw in gi])
         so = sum([slw.ps_counts[source]*slw.overlaps[source] for slw in gi])
         eta = (1.92 - bg)*so
         print 'multi = %.2f'%(eta)
      """
      #define a joint likelihood for the energy band
      f = lambda scale: sum( [band.bandLikelihood(scale,source) for band in gi] )
      result = fmin(f,[1.],disp=0,full_output=1)
      if result[4] != 1 and result[4] != 2:
         scales[i] = result[0]
         f68 = lambda scale: -result[1] + f(scale) - 0.5 #68% confidence
         f95 = lambda scale: -result[1] + f(scale) - 1.92 #95% confidence
         
         if scales[i] < 0.05 or scales[i]*ps[i] < 0.1: #upper limit
            scales[i] = 0
            scales_lo[i] = 0
            scales_hi[i] = fsolve(f95,[max(obs/ps[i],scales[i])])

         else:
            scales_lo[i] = fsolve(f68,[0])
            scales_hi[i] = fsolve(f68,[scales[i]*1.1])

   cts    = scales*ps
   cts_hi = scales_hi*ps
   cts_lo = scales_lo*ps
   

   return r.bin_edges,cts,cts_lo,cts_hi,exp,r.psm.point_sources[source]

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def band_fluxes(r,which=0,axes=None,axis=None,outfile=None, **kwargs):

   plot_kwargs = {'color':'red', 'linewidth':1,}
   plot_kwargs.update(kwargs)
   if axes is None:
      ax = P.gca()
   else: ax = axes

   ax.set_xscale('log')
   ax.set_yscale('log')
   ax.grid(True)

   r.setup_energy_bands()

   for eb in r.energy_bands:
      eb.bandFit(which=which)
      eb.merged = False

   # count in from high end to find how many high energy bins have only upper limits
   for n,neb in enumerate(r.energy_bands[::-1]):
      if neb.flux is not None: break
   # count in from low end to find how many low energy bins have only upper limits
   for nlo,neb in enumerate(r.energy_bands):
      if neb.flux is not None: break
   #reg_slice    = slice(0,len(r.bin_centers) - n,1)
   reg_slice    = slice(nlo,len(r.bin_centers) - n,1)
   merged_slice = slice(len(r.bin_centers)-n, len(r.bin_centers), 1)
   lo_merged_slice = slice(0,nlo,1)

   if n > 0:
      bands = []
      for eb in r.energy_bands[merged_slice]:
         for band in eb.bands:
            bands.append(band)
      ebhi = ROIEnergyBand(bands,r.energy_bands[merged_slice][0].emin,r.energy_bands[merged_slice][-1].emax)
      ebhi.bandFit(which=which)
      ebhi.merged = True
      ebhi.num = n
      ebhi = [ebhi]
   else:
      ebhi = []

   if nlo > 0:
      bands = []
      for eb in r.energy_bands[lo_merged_slice]:
         for band in eb.bands:
            bands.append(band)
      eblo = ROIEnergyBand(bands,r.energy_bands[lo_merged_slice][0].emin,r.energy_bands[lo_merged_slice][-1].emax)
      eblo.bandFit(which=which)
      eblo.merged = True
      eblo.num = nlo
      eblo = [eblo]
   else:
      eblo = []

   r.energy_bands = eblo + r.energy_bands[reg_slice] + ebhi

   for eb in r.energy_bands:

      xlo,xhi = eb.emin,eb.emax
      bc  = (xlo*xhi)**0.5
      fac = xlo*xhi 
      hw  = (xhi-xlo)/2.
      if eb.merged: hw /= eb.num
      
      if eb.flux is not None: 
         ax.plot([xlo,xhi],[fac*eb.flux,fac*eb.flux], **plot_kwargs)
         ax.plot([bc,bc]  ,[fac*eb.lflux,fac*eb.uflux],**plot_kwargs)
      else:
         ax.plot([xlo,xhi],[fac*eb.uflux,fac*eb.uflux],color='k')
         dy = -(fac*eb.uflux - 10**(N.log10(fac*eb.uflux) - 0.6))
         hl = 10**(N.log10(fac*eb.uflux) - 0.4) - 10**(N.log10(fac*eb.uflux) - 0.6)
         a = FancyArrow(bc,fac*eb.uflux,0,dy,width=(xhi-xlo)/100.,head_width=hw,
                                         head_length=hl, length_includes_head=True,
                                         facecolor='k',edgecolor='k',fill=False)
         ax.add_patch(a)

   if axis is not None: ax.axis(axis)
   if outfile is not None: P.savefig(outfile)


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def make_sed(r,which=0,axes=None,axis=None,plot_model=True, 
        data_kwargs=None, fit_kwargs=None):
    if data_kwargs is None: data_kwargs={}
    if fit_kwargs is None: fit_kwargs={'color':'blue', 'linewidth':1}
    if axes is None: axes = P.gca()
    band_fluxes(r,which=which,axes=axes, **data_kwargs)
    axes.set_xlabel('Energy (MeV)')
    axes.set_ylabel('Energy Flux (MeV/cm2/s)')
    if plot_model:
        try:
            dom = N.logspace(N.log10(r.fit_emin[0]),N.log10(r.fit_emax[0]),51)
            cod = r.psm.models[which](dom)*dom**2
            axes.plot(dom,cod, **fit_kwargs)
        except exceptions.OverflowError:
            pass # failed, at least plot axes
        except:
            raise # anything else
    if axis is None:
        axes.axis([1e2,1e5,1e-10,1e-2])
    else: axes.axis(axis)
    axes.grid(True)
    

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
   

def counts(r,integral=False):

   groupings = [deque() for x in xrange(len(r.bin_centers))]
   p = r.phase_factor

   #group slw by energy
   for i,ei in enumerate(r.bin_centers):
      for band in r.bands:
         if band.e == ei:
            groupings[i].append(band)
      groupings[i] = list(groupings[i])

   #iso = N.asarray([ sum((band.bg_counts[1] for band in g)) for g in groupings]) * p
   #gal = N.asarray([ sum((band.bg_counts[0] for band in g)) for g in groupings]) * p
   dif = N.asarray([ N.asarray([band.bg_counts for band in g]).sum(axis=0) for g in groupings])*p
   obs = N.asarray([ sum((band.photons for band in g)) for g in groupings])
   src = N.asarray([ N.asarray([band.ps_counts*band.overlaps for band in g]).sum(axis=0) for g in groupings])*p
   
   if integral:
      for i in xrange(len(iso)):
         #iso[i] = iso[i:].sum()
         #gal[i] = gal[i:].sum()
         dif[i] = dif[i:].sum(axis=0)
         obs[i] = obs[i:].sum()
         src[i] = src[i:].sum(axis=0)

   return r.bin_edges,dif,src,obs,[b.name for b in r.bgm.bgmodels],[p.name for p in r.psm.point_sources]


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


def plot_counts(r,fignum=1,outfile=None,integral=False,max_label=10,merge_non_free=True):

   colors = ['blue','green','red','orange']

   en,dif,src,obs,bg_names,ps_names = counts(r,integral=integral)
   en = (en[1:]*en[:-1])**0.5

   # optionally merge all of the "frozen" sources for a cleaner plot
   if merge_non_free:
      free_mask = N.asarray([N.any(m.free) for m in r.psm.models])
      new_src   = N.zeros([len(en),free_mask.sum()+1])
      counter = 0
      for i in xrange(len(free_mask)):
         if free_mask[i]:
            new_src[:,counter] = src[:,i]
            counter += 1
         else:
            new_src[:,-1] += src[:,i]
      #return en,dif,src,obs,bg_names,ps_names,new_src
      src      = new_src
      ps_names = [ps_names[i] for i in xrange(len(ps_names)) if free_mask[i]]
      ps_names += ['Other Point Sources']
      

   P.clf()
   P.figure(fignum,(14,8))
   P.subplot(121)
   P.gca().set_xscale('log')
   P.gca().set_yscale('log')
   for i,name in enumerate(ps_names):
      #if i < max_label:
      P.loglog(en,src[:,i],linestyle='-',marker='',label=name)
      #else:
      #   P.loglog(en,src[:,i],linestyle='-',marker='')
   for i,name in enumerate(bg_names):
      if N.any(dif[:,i]==0): continue
      P.loglog(en,dif[:,i],linestyle='-',marker='',label=name)
   #if not N.any(gal==0.):
   #   P.loglog(en,gal,linestyle='-',marker='',label='Galactic')

   #if not N.any(iso==0.):
   #   P.loglog(en,iso,linestyle='-',marker='',label='Isotropic')

   #tot = src.sum(axis=1)+iso+gal
   tot = src.sum(axis=1) + dif.sum(axis=1)
   P.loglog(en,tot,linestyle='steps-mid',label='Total Model',color='black')
   err = obs**0.5
   low_err = N.where(obs-err <= 0, 0.99*obs, err)
   P.errorbar(en,obs,yerr=[low_err,err],linestyle=' ',marker='x',label='Counts',color='black')
   ax = P.axis()
   P.axis([ax[0],ax[1],max(0.1,ax[2]),ax[3]])
   if integral:
      P.ylabel('Counts > E')
   else:
      P.ylabel('Counts per Bin')
   P.xlabel('Energy (MeV)')
   P.legend(loc=0)
   P.grid(b=True)

   P.subplot(122)
   P.gca().set_xscale('log')
   P.errorbar(en,(tot-obs)/(tot),yerr=tot**-0.5,linestyle=' ',marker='o',color='black')
   P.plot(N.linspace(P.axis()[0],P.axis()[1],100),[0]*100,color='black')
   ax = P.axis()
   ybound = min( 0.5, max( abs(ax[2]), abs(ax[3]) ))
   P.axis([ax[0],ax[1],-ybound,ybound])
   P.ylabel('Fractional Deviation (model-obs)/model')
   P.xlabel('Energy (MeV)')
   P.grid(b=True)

   if outfile is not None: P.savefig(outfile)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#



def plot_spectra(r, which=0, eweight=2,fignum=1,outfile=None,merge_bins=False,
                 return_vals=False,axis=None,use_legend=True,axes=None):

   colors = ['blue','green','red','orange']

   en,cts,ctslo,ctshi,exp,ps = band_spectra(r,which)

   
   if merge_bins:
      new_en    = []
      new_cts   = []
      new_ctslo = []
      new_ctshi = []
      new_exp   = []
      for i in range(len(cts))[::2]:
         new_en    += [en[i]]         
         new_cts   += [cts[i] + (cts[i+1] if i < (len(cts) - 1) else 0)]
         new_ctslo += [((cts[i]-ctslo[i])**2 + ( (cts[i+1]-ctslo[i+1])**2 if i < (len(ctslo) - 1) else 0))**0.5]
         new_ctshi += [((ctshi[i]-cts[i])**2 + ( (ctshi[i+1]-cts[i+1])**2 if i < (len(ctshi) - 1) else 0))**0.5]
         new_exp   += [(exp[i]*exp[i+(1 if i < (len(cts) - 1) else 0)])**0.5]
      new_en += [en[-1]]

      en    = N.asarray(new_en)
      cts   = N.asarray(new_cts)
      ctslo = cts - N.asarray(new_ctslo)
      ctshi = N.asarray(new_ctshi) + cts
      exp   = N.asarray(new_exp)
   
   bc = (en[1:]*en[:-1])**0.5

   if axes is None:
      P.figure(fignum)
      P.clf()
      axes = P.gca()

   fluxes    = bc**eweight*cts/exp/(en[1:]-en[:-1])
   fluxes_lo = bc**eweight*ctslo/exp/(en[1:]-en[:-1])
   fluxes_hi = bc**eweight*ctshi/exp/(en[1:]-en[:-1])
   P.gca().set_xscale('log')
   P.gca().set_yscale('log')
   if N.all(fluxes==0):
      axes.plot(bc[fluxes==0],fluxes_hi[fluxes==0],marker='v',ls=' ',mfc='black',mec='black')
   else:
      f   = fluxes[fluxes>0]
      flo = fluxes_lo[fluxes>0]
      flo[flo==0] = 1e-20#0.999*f
      fhi = fluxes_hi[fluxes>0]
      axes.errorbar(x=bc[fluxes>0],y=f,yerr=[f-flo,fhi-f],linestyle=' ',marker='o',mfc = 'white', mec = 'black',\
                       color='black',label='Band Fits',ms=10)
      axes.plot(bc[fluxes==0],fluxes_hi[fluxes==0],marker='v',ls=' ',mfc='black',mec='black')
      domain = N.logspace(N.log10(en[0]),N.log10(en[-1]),50)
      axes.plot(domain,domain**eweight*ps.model(domain),color='red',lw=2,label='Model')
      if use_legend: P.legend(loc=0,numpoints=1)
      axes.grid(b=True)
      axes.set_xlabel('Energy (MeV)')
      if axis is None:
         axes.axis([0.7*min(bc),1.3*max(bc),max(min(fluxes_lo[fluxes>0])*0.7,max(fluxes)/100.),max(fluxes_hi)*1.3])
      else:
         axes.axis(axis)

   if outfile is not None: P.savefig(outfile)
   if return_vals: return en,fluxes,fluxes_lo,fluxes_hi
      
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

class ROIModelSkyFunction(object):

   def init(self):
      self.emin  = 1e2
      self.emax  = 1e5
      self.nside = 2**7

      #1 for pixelixed display, #2 for density
      self.mode  = 1
      self.std_weight = True #weight model-obs by 1./root(model)

      self.event_class = -1

   def __init__(self,roi_manager,**kwargs):
      self.init()
      #self.__dict__.update(kwargs)
      for k,v in kwargs.iteritems():
         if k in self.__dict__: self.__dict__[k] = v
   
      self.roi = roi_manager
      sa = self.roi.sa
      self.psf = OldPsf(CALDB = sa.CALDB,irf = sa.irf, )

      self.simpsn = simpsn = (int(round((N.log10(self.emax)-N.log10(self.emin))/0.2)) >> 1) << 1 #5 per decade
      self.points = N.logspace(N.log10(self.emin),N.log10(self.emax),simpsn+1)
      self.simpsf = self.points*N.log(self.points[-1]/self.points[0])*\
                    N.asarray([1.] + ([4.,2.]*(simpsn/2))[:-1] + [1.])/(3.*simpsn)

      self.psf.set_cache(N.append(self.points,self.points),N.append([0]*len(self.points),[1]*len(self.points)))

      self.exp = exp = sa.exposure.exposure
      ps  = self.roi.psm.point_sources
      self.exp_f  = N.asarray([[exp[0].value(p.skydir,point) for point in self.points] for p in ps])
      self.exp_b  = N.asarray([[exp[1].value(p.skydir,point) for point in self.points] for p in ps])
      self.model  = N.asarray([p.model(self.points) for p in ps])

      
      if self.event_class < 0: self.event_classes = [0,1]
      else:
         #setting exposure to 0 a succinct but inefficient way to select on event class
         [self.exp_b,self.exp_f][self.event_class] = N.zeros_like([self.exp_f])
         self.event_classes = [self.event_class]

      bg = self.roi.bgm

      self.e_models  = [m.smodel(self.points) for m in bg.bgmodels]
      self.b_models  = [ [Background(m.get_dmodel(0),exp[0]),Background(m.get_dmodel(1),exp[1])] for m in bg.bgmodels]

      self.cache_pix = dict()
      self.band      = Band(self.nside)
      self.pskyf     = PySkyFunction(self)

   def __call__(self,v,v_is_skydir=False):
   
      skydir = v if v_is_skydir else SkyDir(Hep3Vector(v[0],v[1],v[2]))

      if self.mode == 1:

         hindex = self.band.index(skydir)

         if hindex in self.cache_pix: 
            hval = self.cache_pix[hindex]
         else:
            self.mode = 2
            hval = SkyIntegrator.pix_int(self.pskyf,skydir,self.nside)
            self.cache_pix[hindex]=hval
            self.mode = 1

         return hval

      simpsn,points,simpsf = self.simpsn,self.points,self.simpsf
      ps = self.roi.psm.point_sources

      integrand = N.zeros_like(points)

      #point source contributions
      for i,p in enumerate(ps):
         delta  = skydir.difference(p.skydir)
         psf_vals = self.psf.get_cache(delta,density=True) #gets cached vals
         np = len(points)
         integrand += self.model[i,:] * (self.exp_f[i,:] * psf_vals[:np] + self.exp_b[i,:] * psf_vals[np:])
      
      #diffuse contributions
      for em,bm in zip(self.e_models,self.b_models):
         for ec in self.event_classes:
            integrand += em * N.asarray([bm[ec].value(skydir,point) for point in self.points])

      return (simpsf * integrand).sum()

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

class DataSkyFunction(object):
   """Build a SkyFunction from FT1 files.

   Arguments:
      sa_or_data  --- provide either a list of FT1 files or a SpectralAnalysis object
      
   Keyword Arguments:
      nside       --- [2^7] Healpix resolution parameter; resolution ~ 60deg / nside
      mode        --- [1] == 1 for pixelized representation -- mode 2 currently broken!
      cuts        --- [None] a list of cuts to apply to the data, e.g., ['ENERGY > 1000','CTBCLASSLEVEL == 3']
      emin        --- [1e2] minimum energy
      emax        --- [1e5] maximum energy

   Note: implements the PySkyFunction interface, so. e.g, to use with SkyImage, just cast it
   """

   def init(self):

      self.nside   = 2**7
      self.mode    = 1
      self.cuts    = None
      self.emin    = 1e2
      self.emax    = 1e5
      self.event_class = -1
      self.mc_src_id = None

   def __init__(self,sa_or_data,**kwargs):

      self.init()
      for k,v in kwargs.iteritems():
         if k in self.__dict__: self.__dict__[k] = v

      self.sa_or_data = sa_or_data
      self.mode1setup = self.mode2setup = False

      #if self.mode == 1: self.mode_1_setup()
      #else: self.mode_2_setup()
      self.mode_1_setup()

      self.cache_pix = dict()

   def process_filedata(self,fields):

      base_cuts = ['ENERGY > %s'%(self.emin),'ENERGY < %s'%(self.emax),'ZENITH_ANGLE < 105']
      if self.event_class >= 0:      base_cuts += ['EVENT_CLASS == %d'%(self.event_class)]
      if self.mc_src_id is not None: base_cuts += ['MC_SRC_ID == %d'%(self.mc_src_id)]
      cuts = base_cuts if self.cuts is None else self.cuts + base_cuts
      event_files = self.sa_or_data if type(self.sa_or_data) == type([]) else self.sa_or_data.pixeldata.ft1files
      
      data = get_fields(event_files,fields,cuts)

      return data

   def mode_1_setup(self):

      data = self.process_filedata(['RA','DEC'])

      self.band = Band(self.nside)

      for i in xrange(len(data['RA'])):
         self.band.add(SkyDir(float(data['RA'][i]),float(data['DEC'][i])))

      #for skydir in map(SkyDir,data['RA'],data['DEC']): self.band.add(skydir)

      self.mode1setup = True

   """
   def mode_2_setup(self):

      data = self.process_filedata(['RA','DEC','ENERGY','EVENT_CLASS'])

      self.skydirs = map(SkyDir,data['RA'],data['DEC'])
      
      from psf import PSF
      if type(self.sa_or_data) != type([]):
         sa = self.sa_or_data         
         psf = PSF(use_mc_psf = sa.use_mc_psf, use_psf_init = sa.use_psf_init, irf = sa.irf)
      else:
         psf = PSF(use_mc_psf = True, use_psf_init = True, irf = 'P6_v3_diff')

      self.sigmas = psf.sigma(data['ENERGY'],data['EVENT_CLASS'])*N.pi/180.
      self.gammas = psf.gamma(data['ENERGY'],data['EVENT_CLASS'])

      self.mode2setup = True
   """
   
   def __call__(self,v,v_is_skydir=False):

      v = v if v_is_skydir else SkyDir(Hep3Vector(v[0],v[1],v[2]))
      
      #if self.mode == 1:
      if not self.mode1setup: self.mode_1_setup()
      ind = self.band.index(v)
      
      #need to store these values, so this check is necessary
      if not ind in self.cache_pix:
         val = self.band(v)
         self.cache_pix[ind] = val
         return val
      
      else: return self.cache_pix[ind] #but need to check whether the O(1) hash quicker than the Band call -- suspect it is!

      #else:
      #   if not selt.mode2setup: self.mode_2_setup()
      #   diffs = N.asarray([v.difference(sd) for sd in self.skydirs])
      #   us    = 0.5 * (diffs/self.sigmas)**2
      #   jac   = (2*N.pi)*(self.sigmas)**2
      #   return ( (1.-1./self.gammas)*(1+u/self.gammas)**-self.gammas )/jac


class ROIDisplay(object):
   """Manage the plotting of ROI info."""

   def init(self): 
      
      self.figsize = (12,8)
      self.fignum  = 3
      
      self.nside = 2**7
      self.pixelsize = 0.25
      self.size = 10
      self.nticks = 5

      self.log_counts_min = 1.0
      self.log_counts_max = 2.5

      self.label_sources = False

      self.galactic = True

      self.std_weight = True

      self.emin = 1e2
      self.emax = 1e5

      self.event_class = -1

      self.mm = self.cm = None

   def mathrm(self,st):
      return  r'$\mathrm{'+st.replace(' ','\/')+'}$'

   def __init__(self, roi_manager, **kwargs):
      
      self.init()
      self.__dict__.update(kwargs)
      self.rm = roi_manager

      rcParams['xtick.major.size']=10 #size in points
      rcParams['xtick.minor.size']=6
      rcParams['ytick.major.size']=10 #size in points
      rcParams['ytick.minor.size']=6
      rcParams['xtick.labelsize']=12
      rcParams['ytick.labelsize']=12
      rcParams['font.family']='serif'
      #rcParams['text.usetex']=True

      try:
         self.cmap_sls = colormaps.sls
         self.cmap_b   = colormaps.b
      except:
         self.cmap_sls = self.cmap_b = mpl.cm.jet

      P.ioff()
      fig = P.figure(self.fignum,self.figsize)
      P.clf()

      voff = -0.05
      imw = ( 0.55, 0.55, 0.55,  0.55,     0.25, 0.25) #0.38, 0.38)
      imh = ( 0.35, 0.35, 0.35,  0.35,     0.17, 0.17) # picture sizes
      imx = (-0.09, 0.23, 0.55, -0.09,     0.40, 0.40) #0.55, 0.55)
      imy = ( 0.60+voff, 0.60+voff, 0.60+voff,  0.15+voff,     0.32, 0.08)


      titles  = ['Modeled Counts', 'Observed Counts', 'Weighted Residuals', 'P-values','', '']
      xlabels = ['RA']*4 + ['p-values','weighted residuals'] 
      ylabels = ['Dec']*4 + ['Frequency','Frequency']

      for y in [titles,xlabels,ylabels]:
         for i in xrange(len(y)):
            y[i] = self.mathrm(y[i])

      for i,(t,xl,yl) in enumerate(zip(titles,xlabels,ylabels)):
         ax = self.__dict__['axes%d'%(i+1)] = fig.add_axes([imx[i],imy[i],imw[i],imh[i]])
         if i < 4: ax.set_aspect(1)
         ax.set_title(t)
         ax.set_xlabel(xl) 
         ax.set_ylabel(yl)

      self.norm1 = mpl.colors.Normalize(vmin=0,vmax=5)
      self.norm2 = mpl.colors.Normalize(vmin=self.log_counts_min,vmax=self.log_counts_max)
      self.norm3 = mpl.colors.Normalize(vmin=-5,vmax=5)

      #resid colorbar
      rcb_axes = fig.add_axes([imx[3]+imw[3]*0.72,imy[3],0.01,imh[3]])
      #rcb_norm = mpl.colors.Normalize(vmin=0,vmax=5)
      rcb      = mpl.colorbar.ColorbarBase(rcb_axes,
                                             norm=self.norm1,
                                             ticks = ticker.MaxNLocator(4),
                                             cmap  = colormaps.b,
                                             orientation='vertical')
      #rcb.set_label('')

      #counts colorbar
      ccb_axes1 = fig.add_axes([imx[0]+imw[0]*0.72,imy[0],0.01,imh[0]])
      ccb1      = mpl.colorbar.ColorbarBase(ccb_axes1,
                                             norm=self.norm2,
                                             ticks = ticker.MaxNLocator(4),
                                             cmap  = colormaps.b,
                                             orientation='vertical')
      ccb_axes2 = fig.add_axes([imx[1]+imw[1]*0.72,imy[1],0.01,imh[1]])
      ccb2      = mpl.colorbar.ColorbarBase(ccb_axes2,
                                             norm=self.norm2,
                                             ticks = ticker.MaxNLocator(4),
                                             cmap  = colormaps.b,
                                             orientation='vertical')

      ccb_axes3 = fig.add_axes([imx[2]+imw[2]*0.72,imy[2],0.01,imh[2]])
      ccb3      = mpl.colorbar.ColorbarBase(ccb_axes3,
                                             norm=self.norm3,
                                             ticks = ticker.MaxNLocator(4),
                                             #cmap  = colormaps.b,
                                             orientation='vertical')
      #rcb.set_label('')

      self.mm = ROIModelSkyFunction(roi_manager,mode=1,**self.__dict__) if self.mm is None else self.mm
      self.cm = DataSkyFunction(roi_manager.sa,mode=1,**self.__dict__) if self.cm is None else self.cm

   def model_plot(self):

      self.mm_zea = ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes1)
      self.mm_zea.fill(PySkyFunction(self.mm))
      self.axes1.imshow(N.log10(self.mm_zea.image),origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm2)
      self.mm_zea.grid()
      self.mm_zea.scale_bar(color='white')
      self.plot_sources(self.mm_zea,mc='k')
      
   def counts_plot(self):

      self.cm_zea = ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes2)
      self.cm_zea.fill(PySkyFunction(self.cm))
      self.axes2.imshow(N.log10(self.cm_zea.image),origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm2)
      self.cm_zea.grid()
      self.cm_zea.scale_bar(color='white')
      self.plot_sources(self.cm_zea,mc='k')

   def resids_plot(self):

      self.resids_zea = ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes3)
      self.resids_zea.image = (self.mm_zea.image - self.cm_zea.image)/self.mm_zea.image**0.5
      self.axes3.imshow(self.resids_zea.image,origin='lower',interpolation='nearest',norm=self.norm3)
      self.resids_zea.grid()
      self.plot_sources(self.resids_zea,mc='k')
      
      pvals = 1 - poisson.cdf(self.cm_zea.image,self.mm_zea.image ) #0 problem?
      pvals = N.abs(N.where(pvals < 0.5, N.log10(pvals), N.log10(1-pvals)))
      #pvals = N.abs( N.tan ( pvals*N.pi - N.pi/2 ) )

      self.pvals_zea = ZEA(self.rm.psm.ROI_dir(),self.size,self.pixelsize,galactic=self.galactic,axes=self.axes4)
      self.pvals_zea.image = pvals
      self.axes4.imshow(self.pvals_zea.image,origin='lower',interpolation='nearest',cmap=self.cmap_b,norm=self.norm1)
      self.pvals_zea.grid()

      self.plot_sources(self.pvals_zea,mc='k')

   def hist_plot(self):

      mc = N.asarray(self.mm.cache_pix.values())
      cc = N.asarray(self.cm.cache_pix.values())

      pvals = N.asarray([1-poisson.cdf(cc[i],mc[i]) for i in xrange(len(mc))])

      nb = 20
      av = float(len(pvals)) / nb
      self.axes5.hist(pvals,bins=N.linspace(0,1,20),histtype='step')
      self.axes5.axhline( av, color='red')
      lo,hi = ppf( (50.-95/2)/100., av), ppf( (50. + 95/2)/100.,av)
      self.axes5.axhspan( lo , hi , facecolor='red', alpha=0.3,label='95% Conf.')
      self.axes5.legend(loc='upper right')

      self.axes6.hist( (mc - cc)/mc**0.5, bins=N.linspace(-5,5,20), histtype='step')
      self.axes6.axvline(0,color='red')
      self.axes6.grid(True)
      self.axes6.set_xbound(lower=-5,upper=5)


   def show(self,to_screen=True,out_file=None):

      t =self.label_sources
      self.model_plot()
      self.label_sources=False #only label the model plot
      self.counts_plot()
      self.resids_plot()
      self.hist_plot()
      self.label_sources = t

      if out_file is not None: P.savefig(out_file)
      if to_screen: P.show()

   def plot_sources(self, image, symbol='+',  fontsize=8, markersize=10, fontcolor='w', mc= 'green',**kwargs):
      nx = image.nx
      ny = image.ny
      ps = self.rm.psm.point_sources

      def allow(nx, ny, px, py, padx = 0.15, pady = 0.15):
         padx = padx * nx
         pady = pady * ny

         return (px > padx and px < nx - padx) and (py > pady and py < ny - pady)

      for p in ps:
         x,y = image.pixel(p.skydir)
         if self.label_sources:
            padx = pady = 0.15
         else:
            padx = pady = 0.02
         if not allow(nx,ny,x,y,padx,pady): continue
         image.axes.plot([x], [y], symbol, markersize=markersize, mec = mc, mfc = mc,**kwargs)
         image.axes.plot([x], [y], 'x', markersize=markersize, mec = 'white', mfc = 'white',**kwargs)
         if self.label_sources:
            image.axes.text( x+nx/100., y+nx/100., p.name, fontsize=fontsize, color=fontcolor)




#===============================================================================================#

def ppf(prob,mean):
   """Return the (approximate) Poisson percentage point function for given distribution.  Klugey."""
   if mean > 100: #normal approximation
      n = norm(mean,mean**0.5)
      return n.ppf(prob)      
   d = poisson(mean)
   prev = 0
   for i in xrange(1,200):      
      new = d.cdf(i)
      if new >= prob: break
      prev = new
   return (i-1)+(prob-prev)/(new-prev) #linear interpolation

#===============================================================================================#

def int2bin(n, count=24):
   """returns the binary of integer n, using count number of digits"""
   return "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])



