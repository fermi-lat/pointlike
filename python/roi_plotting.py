"""
Plotting routines to display results of an ROI analysis.
Given an ROIAnalysis object roi:
    plot_spectra(roi)
    ROIDisplay(roi).show()
    plot_counts(roi)


$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/roi_plotting.py,v 1.6 2009/07/22 02:56:54 burnett Exp $

author: Matthew Kerr
"""

import numpy as N
from skymaps import PySkyFunction
import colormaps
from image import ZEA

def band_spectra(r,source=0):

   from scipy.optimize import fmin,fsolve
   from collections import deque
   
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
      f = lambda scale: sum( [r.bandLikelihood(scale,band,source) for band in gi] )
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

   

def counts(r,integral=False):

   from collections import deque
   groupings = [deque() for x in xrange(len(r.bin_centers))]

   #group slw by energy
   for i,ei in enumerate(r.bin_centers):
      for band in r.bands:
         if band.e == ei:
            groupings[i].append(band)
      groupings[i] = list(groupings[i])

   #kluge for now, assumes standard ordering of backgrounds...

   iso = N.asarray([ sum((band.bg_counts[1] for band in g)) for g in groupings])
   gal = N.asarray([ sum((band.bg_counts[0] for band in g)) for g in groupings])
   obs = N.asarray([ sum((band.photons for band in g)) for g in groupings])
   src = N.asarray([ N.asarray([band.ps_counts*band.overlaps for band in g]).sum(axis=0) for g in groupings])
   
   if integral:
      for i in xrange(len(iso)):
         iso[i] = iso[i:].sum()
         gal[i] = gal[i:].sum()
         obs[i] = obs[i:].sum()
         src[i] = src[i:].sum(axis=0)
   return r.bin_edges,iso,gal,src,obs,[p.name for p in r.psm.point_sources]


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


def plot_counts(r,fignum=1,outfile=None,integral=False,):

   colors = ['blue','green','red','orange']

   en,iso,gal,src,obs,ps_names = counts(r,integral=integral)
   en = (en[1:]*en[:-1])**0.5

   import pylab as P
   P.clf()
   P.figure(fignum,(14,8))
   P.subplot(121)
   P.gca().set_xscale('log')
   P.gca().set_yscale('log')
   for i,name in enumerate(ps_names):
      P.loglog(en,src[:,i],linestyle='-',marker='',label=name)

   if not N.any(gal==0.):
      P.loglog(en,gal,linestyle='-',marker='',label='Galactic')

   if not N.any(iso==0.):
      P.loglog(en,iso,linestyle='-',marker='',label='Isotropic')

   tot = src.sum(axis=1)+iso+gal
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
            return_vals=False,axis=None,use_legend=True):

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


   import pylab as P
   
   P.figure(fignum)
   P.clf()

   fluxes    = bc**eweight*cts/exp/(en[1:]-en[:-1])
   fluxes_lo = bc**eweight*ctslo/exp/(en[1:]-en[:-1])
   fluxes_hi = bc**eweight*ctshi/exp/(en[1:]-en[:-1])
   P.gca().set_xscale('log')
   P.gca().set_yscale('log')
   if N.all(fluxes==0):
      P.plot(bc[fluxes==0],fluxes_hi[fluxes==0],marker='v',ls=' ',mfc='black',mec='black')
   else:
      f   = fluxes[fluxes>0]
      flo = fluxes_lo[fluxes>0]
      flo[flo==0] = 1e-20#0.999*f
      fhi = fluxes_hi[fluxes>0]
      P.errorbar(x=bc[fluxes>0],y=f,yerr=[f-flo,fhi-f],linestyle=' ',marker='o',mfc = 'white', mec = 'black',\
                       color='black',label='Band Fits',ms=10)
      P.plot(bc[fluxes==0],fluxes_hi[fluxes==0],marker='v',ls=' ',mfc='black',mec='black')
      domain = N.logspace(N.log10(en[0]),N.log10(en[-1]),50)
      P.plot(domain,domain**eweight*ps.model(domain),color='red',lw=2,label='Model')
      if use_legend: P.legend(loc=0,numpoints=1)
      P.grid(b=True)
      P.xlabel('Energy (MeV)')
      if axis is None:
         P.axis([0.7*min(bc),1.3*max(bc),max(min(fluxes_lo[fluxes>0])*0.7,max(fluxes)/100.),max(fluxes_hi)*1.3])
      else:
         P.axis(axis)

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
   
      from psf import PSF
      self.roi = roi_manager
      sa = self.roi.sa
      self.psf = PSF(psf_irf = sa.psf_irf, CALDB = sa.CALDB)

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
      from skymaps import Background

      self.e_models  = [m.smodel(self.points) for m in bg.bgmodels]
      self.b_models  = [ [Background(m.get_dmodel(0),exp[0]),Background(m.get_dmodel(1),exp[1])] for m in bg.bgmodels]

      self.cache_pix = dict()
      from skymaps import Band
      self.band      = Band(self.nside)
      from skymaps import PySkyFunction
      self.pskyf     = PySkyFunction(self)

   def __call__(self,v,v_is_skydir=False):
   
      from skymaps import SkyDir,Hep3Vector,Band,SkyIntegrator
      
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
         psf_vals = self.psf(1,delta,1,radians=True,density=True) #gets cached vals
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
      mode        --- [1] ==1 for a pixelized representation, ==2 for a density estimator (slow!)
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

      if self.mode == 1: self.mode_1_setup()
      else: self.mode_2_setup()

      self.cache_pix = dict()

   def process_filedata(self,fields):

      from fitstools import get_fields
      base_cuts = ['ENERGY > %s'%(self.emin),'ENERGY < %s'%(self.emax),'ZENITH_ANGLE < 105']
      if self.event_class >= 0:      base_cuts += ['EVENT_CLASS == %d'%(self.event_class)]
      if self.mc_src_id is not None: base_cuts += ['MC_SRC_ID == %d'%(self.mc_src_id)]
      cuts = base_cuts if self.cuts is None else self.cuts + base_cuts
      event_files = self.sa_or_data if type(self.sa_or_data) == type([]) else self.sa_or_data.pixeldata.event_files
      
      data = get_fields(event_files,fields,cuts)

      return data

   def mode_1_setup(self):

      data = self.process_filedata(['RA','DEC'])
      
      from skymaps import Band,SkyDir
      self.band = Band(self.nside)

      for i in xrange(len(data['RA'])):
         self.band.add(SkyDir(float(data['RA'][i]),float(data['DEC'][i])))

      #for skydir in map(SkyDir,data['RA'],data['DEC']): self.band.add(skydir)

      self.mode1setup = True

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

   def __call__(self,v,v_is_skydir=False):
      
      from skymaps import Hep3Vector,SkyDir
      v = v if v_is_skydir else SkyDir(Hep3Vector(v[0],v[1],v[2]))
      
      if self.mode == 1:
         if not self.mode1setup: self.mode_1_setup()
         ind = self.band.index(v)
         
         #need to store these values, so this check is necessary
         if not ind in self.cache_pix:
            val = self.band(v)
            self.cache_pix[ind] = val
            return val
         
         else: return self.cache_pix[ind] #but need to check whether the O(1) hash quicker than the Band call -- suspect it is!

      else:
         if not selt.mode2setup: self.mode_2_setup()
         diffs = N.asarray([v.difference(sd) for sd in self.skydirs])
         us    = 0.5 * (diffs/self.sigmas)**2
         jac   = (2*N.pi)*(self.sigmas)**2
         return ( (1.-1./self.gammas)*(1+u/self.gammas)**-self.gammas )/jac


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

      from matplotlib import rcParams
      rcParams['xtick.major.size']=10 #size in points
      rcParams['xtick.minor.size']=6
      rcParams['ytick.major.size']=10 #size in points
      rcParams['ytick.minor.size']=6
      rcParams['xtick.labelsize']=12
      rcParams['ytick.labelsize']=12
      rcParams['font.family']='serif'
      #rcParams['text.usetex']=True
   
      import pylab as P
      from matplotlib import mpl,pyplot,ticker
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

      from scipy.stats import poisson

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

      import pylab as P
      from scipy.stats import poisson

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

      import pylab as P
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




#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

class PSFChecker(object):

   def init(self):

      self.emin = 100
      self.emax = 1e5
      self.bins_per_decade = 3
      self.umax = 100
      self.psf_irf = 'P6_v3_diff'
      self.CALDB   = None
      self.nsimps = 16
   
   def __init__(self,roi,**kwargs):

      self.init()
      self.__dict__.update(**kwargs)
      self.roi = roi
      self.sa  = roi.sa

      self.CALDB = self.CALDB or self.sa.CALDB

      lmin,lmax    = N.log10([self.emin,self.emax])
      self.bands   = N.logspace(lmin,lmax,round( (lmax-lmin)*self.bins_per_decade ) + 1)
      self.bands_c = (self.bands[:-1]*self.bands[1:])**0.5
      
      from collections import defaultdict
      self.keys     = self.bands_c.astype(int)
      self.bin_dict = defaultdict(dict)

      from psf import PSF
      self.psf = PSF(psf_irf = self.psf_irf, CALDB = self.CALDB)

      n = len(self.bands_c)
      self.sigmas = self.psf.sigma(N.append(self.bands_c,self.bands_c), N.append([0] * n, [1] * n) ) * (N.pi/180.)
      self.gammas = self.psf.gamma(N.append(self.bands_c,self.bands_c), N.append([0] * n, [1] * n) )

      for i,key in enumerate(self.keys):
         
         #front events, key with positive energy
         self.bin_dict[key]['sigma'] = self.sigmas[i]
         self.bin_dict[key]['gamma'] = self.gammas[i]

         #back events, key with negative energy
         self.bin_dict[-key]['sigma'] = self.sigmas[i+n]
         self.bin_dict[-key]['gamma'] = self.gammas[i+n]

      self.u_bins = N.logspace(-1.5,2,15)

      self.bindata()
      self.ps_models()
      self.bg_models()
        
   
   def bindata(self):

      from fitstools import rad_extract
      max_rad   = self.psf.sigma(self.bands_c[0],1)*(2*self.umax)**0.5 #rad_extract expects degrees
      data_dict = rad_extract(self.sa.pixeldata.event_files, self.sa.roi_dir,
                              max_rad, return_cols=['ENERGY','EVENT_CLASS'])
      
      diffs = data_dict['DIFFERENCES']
      ens   = data_dict['ENERGY']
      ecs   = data_dict['EVENT_CLASS']

      for i,key in enumerate(self.keys):

         emin,emax = self.bands[i],self.bands[i+1]
         mask      = N.logical_and(ens >= emin, ens < emax)
         
         for k,ec in zip([key,-key],[0,1]):

            us    = 0.5 * (diffs[N.logical_and(mask,ecs==ec)] / self.bin_dict[k]['sigma']) ** 2
            self.bin_dict[k]['data'] = N.histogram(us,bins=self.u_bins,new=True)[0]        

   
   def ps_models(self):
      
      psm = self.roi.psm
      exp = self.sa.exposure.exposure
      #need to calculate the overlap for all bin radii and all sources

      from roi_modules import ROIOverlap

      f  = self.psf.overlap
      dirs = [ps.skydir for ps in psm.point_sources]
      d1   = dirs[0]

      for i,key in enumerate(self.keys):

         emin,emax = self.bands[i:i+2]

         for k,ec in zip([key,-key],[0,1]):

            radii = ((2*self.u_bins)**0.5 * self.bin_dict[k]['sigma'])
            self.bin_dict[k]['ne_counts'] = N.zeros_like(radii)

            for n,ps in enumerate(psm.point_sources):
               fracs = [float(f(d1,ps.skydir,self.bands_c[i],ec,r)) for r in radii]
               expec = ps.model.expected(emin,emax,exp,ps.skydir,event_class=ec)
               if n == 0:
                  self.bin_dict[k]['ps_counts'] = N.asarray(fracs) * expec
               else:
                  self.bin_dict[k]['ne_counts'] += N.asarray(fracs) * expec

            #adjust to annular counts
            self.bin_dict[k]['ps_counts'] = self.bin_dict[k]['ps_counts'][1:] - self.bin_dict[k]['ps_counts'][:-1]
            self.bin_dict[k]['ne_counts'] = self.bin_dict[k]['ne_counts'][1:] - self.bin_dict[k]['ne_counts'][:-1]

  
   def bg_models(self):

      bg  = self.roi.bgm
      sa  = self.sa
      rd  = self.sa.roi_dir
      exp = sa.exposure.exposure
      from skymaps import Background,SkyIntegrator

      simps_vec = (N.asarray([1,4,2,4,2,4,2,4,1]))/(3.*8.)

      e_models = [m.smodel for m in bg.bgmodels]
      b_models = [Background(m.dmodel,exp[0],exp[1]) for m in bg.bgmodels]

      
      for i,key in enumerate(self.keys):

         emin, emax  = self.bands[i:i+2]
         lemin,lemax = N.log10([emin,emax])
         simps_pts   = N.logspace(lemin,lemax,9)

         for k,ec in zip([key,-key],[0,1]):

            radii = ((2*self.u_bins)**0.5 * self.bin_dict[k]['sigma'])
            self.bin_dict[k]['bg_counts'] = N.zeros_like(radii)

            for nrad,rad in enumerate(radii):

               solid = N.pi * rad**2

               for n in xrange(len(e_models)):
                  
                  em = e_models[n]
                  bm = b_models[n]

                  bm.set_event_class(ec)

                  bm_pts = N.empty_like(simps_pts)
                  for ne,e in enumerate(simps_pts):
                     bm.setEnergy(e)
                     bm_pts[ne] = SkyIntegrator.ss_average(bm,rd,rad)

                  self.bin_dict[k]['bg_counts'][nrad] += \
                     (bm_pts * em(simps_pts) *simps_pts * simps_vec).sum() * (N.log(emax/emin) * solid )
                     
            self.bin_dict[k]['bg_counts'] = self.bin_dict[k]['bg_counts'][1:] - self.bin_dict[k]['bg_counts'][:-1]


         
   def show(self):

      for i,key in enumerate(self.keys):
         import pylab as P
         from math import floor

         emin,emax = self.bands[i:i+2]

         for k,fignum in zip([key,-key],[5,6]):
            if i > 5: continue
            P.figure(fignum,figsize=(12,10))
            ax = P.subplot( 2, 3, i + 1)
            ax.set_xscale('log')
            ax.set_yscale('log')
            x = N.zeros(len(self.u_bins) * 2 )
            
            #data
            y = N.ones(len(self.u_bins) * 2 ) * 1e-10
            x[0::2],x[1::2] = self.u_bins,self.u_bins
            y[1:-1:2],y[2::2] = self.bin_dict[k]['data'],self.bin_dict[k]['data']
            #ax.fill(x,y,closed=False,fill=False,edgecolor='k',label='%d-%d MeV'%(int(emin),int(emax)))
            ax.errorbar(x = (self.u_bins[1:]*self.u_bins[:-1])**0.5, y = self.bin_dict[k]['data'],
                        yerr = (self.bin_dict[k]['data'])**0.5, ls= ' ', capsize=4,
                        label='%d-%d MeV'%(int(emin),int(emax)),color='k')
            
            #models
            ps_diffs = self.bin_dict[k]['ps_counts']
            ne_diffs = self.bin_dict[k]['ne_counts']
            bg_diffs = self.bin_dict[k]['bg_counts']
            to_diffs = ne_diffs + bg_diffs + ps_diffs
            
            y = N.zeros(len(self.u_bins) * 2 )
            y[1:-1:2],y[2::2] = ps_diffs,ps_diffs
            ax.fill(x,y,closed=False,fill=False,edgecolor='blue')

            y = N.zeros(len(self.u_bins) * 2 )
            y[1:-1:2],y[2::2] = to_diffs,to_diffs
            ax.fill(x,y,closed=False,fill=False,edgecolor='red')

            a = ax.axis()
            ymin = 10**(floor(N.log10(to_diffs.min())))
            ax.axis([self.u_bins[0],self.u_bins[-1],ymin,a[3]])
            ax.legend(loc='upper left')
            ax.set_xlabel('u (bot), deg (top)')
            P.suptitle('Radial Counts for Event Class: %d, Black: Data, Red: Total Model, Blue: Point Source'%(0 if k > 0 else 1),size='large')

            ax2 = ax.twiny()
            rad = self.bin_dict[k]['sigma']*(2*self.u_bins)**0.5*180/N.pi
            ax2.axis([rad[0],rad[-1],ymin,a[3]])


#===============================================================================================#

def ppf(prob,mean):
   """Return the (approximate) Poisson percentage point function for given distribution.  Klugey."""
   if mean > 100: #normal approximation
      from scipy.stats import norm
      n = norm(mean,mean**0.5)
      return n.ppf(prob)      
   
   from scipy.stats import poisson
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

"""       
                  
def loglikelihood(v,psfc):

   self = psfc

   ll = 0;
   nk = len(self.keys)
   bd = self.bin_dict

   for i,key in enumerate(self.keys):
      for k,n in zip([key,-key],[0,nk]):
         scale = v[i+n]
         rate = (scale *bd[k]['ps_counts']) + bd[k]['bg_counts'] + bd[k]['ne_counts']
         ll += (rate - bd[k]['data']*N.log(rate)).sum()
   return ll

def fit_scale(psfc):

   self = psfc
   psfc.loglikelihood = loglikelihood

   scale = N.asarray([1]*(2*len(self.keys)))
   from scipy.optimize import fmin
   f = fmin(self.loglikelihood,scale,args=(psfc,))
   return f

"""

def uhists(roi,bins=N.arange(0,11),show=True,integral=False,residuals=False):
   """Assume only source is central."""

   bands = roi.bands
   rd    = roi.sa.roi_dir

   from skymaps import WeightedSkyDirList,PsfSkyFunction

   for band in bands:
      if not band.has_pixels: continue

      if band.ps_pix_counts[:,0].sum() > 0.8:
         print 'using pixel values for energy %d and event class %d'%(band.e,band.ec)
         wsdl          = WeightedSkyDirList(band.b,rd,band.radius_in_rad,True)
         psf           = PsfSkyFunction(rd,band.g,band.s)
         ps_pix_counts = N.asarray(psf.wsdl_vector_value(wsdl)) * (band.b.pixelArea()/(2*N.pi*band.s**2))
         diffs         = N.asarray([rd.difference(di) for di in wsdl])
         counts        = N.asarray([d.weight() for d in wsdl])
         us            = 0.5 * (diffs/band.s)**2
         rates         = N.histogram(us,bins=bins,new=True,weights=band.ps_counts[0]*ps_pix_counts)[0]
      
      else:
         diffs         = N.asarray([rd.difference(di) for di in band.wsdl])
         us            = 0.5 * (diffs/band.s)**2
         counts        = band.pix_counts
         fs            = 1 - (1 + bins/band.g)**(1-band.g)
         rates         = band.ps_counts[0] * (fs[1:] - fs[:-1])
      
      band.urates = rates 
      band.curates = N.asarray([rates[x:].sum() for x in xrange(len(rates))])
      band.uhist  = N.histogram(us,bins=bins,new=True,weights=counts)
      band.cucounts = N.asarray([band.uhist[0][x:].sum() for x in xrange(len(rates))])

   if show:

      import pylab as P

      x = N.empty(2*len(bins))
      x[::2],x[1::2] = bins,bins
      xc = (bins[1:] + bins[:-1]).astype(float)/2

      P.ioff()

      for ec in [0,1]:
         P.figure(10+ec,(14,12))
         fbands = [band for band in bands if band.ec == ec]
         n = len(fbands)
         nside  = min(int ( n**0.5 ) + 1 , 4)
         for i in xrange(nside):
            for j in xrange(nside):
               ind = nside*i + j
               #print i,j,ind
               if (ind >= n or ind > nside**2): break
               b = fbands[ind]
               if not b.has_pixels: continue
               P.subplot(nside,nside,ind+1)
               y = N.zeros(2*len(bins))

               if residuals:
                  
                  if not integral:
                     P.errorbar(xc,(b.uhist[0]-b.urates)/b.urates,yerr=b.urates**-0.5,ls=' ',marker='o',color='red')
                  else:
                     P.errorbar(xc,(b.cucounts-b.curates)/b.curates,yerr=b.curates**-0.5,ls=' ',marker='o',color='red')
                  P.axhline(0,color='k')
                  ax = P.axis()
                  P.axis([bins[0],bins[-1],-1,1])
               else:
                  if not integral:
                     y[1:-1:2],y[2::2] = b.urates,b.urates #b.uhist[0],b.uhist[0]
                  else:
                     y[1:-1:2],y[2::2] = b.curates,b.curates
                  P.fill(x,y,closed=False,fill=False,edgecolor='blue')
                  if not integral:
                     P.errorbar(xc,b.uhist[0],yerr=b.urates**0.5,ls=' ',marker='o',color='red')
                  else:
                     P.errorbar(xc,b.cucounts,yerr=b.curates**0.5,ls=' ',marker='o',color='red')
               
                  ax = P.axis()
                  P.axis([bins[0],bins[-1],0,ax[3]])
               P.title('%d-%d MeV'%(b.emin,b.emax))
              


