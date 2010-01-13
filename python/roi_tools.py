
import numpy as N

def profile(roi,param,points):
   """Find likelihood profile for given profile.  Assume we are talking about the source at the center of the ROI."""
   
   likelihood_vals = N.zeros_like(points)
   m = roi.ps_manager.models[0]
   save = (m.p[param],m.free[param])

   for i,point in enumerate(points):
      m.p[param] = N.log10(point)
      m.free[param] = True
      roi.logLikelihood(roi.get_parameters())
      m.free[param] = False
      likelihood_vals[i] = roi.fit(save_values = False)

   m.p[param] = save[0]
   m.free[param] = save[1]
   roi.logLikelihood(roi.get_parameters())

   return likelihood_vals

###====================================================================================================###


from types import FunctionType,MethodType
from pypsf import PsfOverlap
from skymaps import Background,SkyIntegrator

def snr_setup(roi,radius_function):
   """Calculate the SNR (or TS) as a function of aperture/energy.
   
      radius_func = a function of energy, conversion_type giving the aperture in deg
      can be a fixed aperture (a scalar)"""

   #for key in ['bands','psm','bgm']:
   #   locals()[key] = roi.__dict__[key]
   bands         = roi.bands
   psm           = roi.psm
   bgm           = roi.bgm
   point_sources = psm.point_sources
   roi_dir = rd  = psm.roi_dir
   nm            = len(bgm.models)
   ns            = bgm.nsimps
   exp           = bgm.sa.exposure.exposure

   front_bgs = [Background(model.get_dmodel(0),exp[0]) for model in bgm.bgmodels]
   back_bgs  = [Background(model.get_dmodel(1),exp[1]) for model in bgm.bgmodels]


   if not (type(radius_function) is FunctionType or type(radius_function) is MethodType):
      simple_scalar = True
      rval = radius_function
      radius_function = lambda e,event_class: rval
   else: simple_scalar = False

   
   overlap = PsfOverlap()
   for nband,band in enumerate(bands):

      # establish new overlaps for aperture
      radius_in_rad          = N.radians(radius_function(band.e,band.ct))
      denom                  = band.exp.value(roi_dir,band.e)
      exposure_ratios        = N.asarray([band.exp.value(ps.skydir,band.e)/denom for ps in point_sources])
      band.snr_overlaps      = N.asarray([overlap(band,roi_dir,ps.skydir,radius_in_rad) for ps in point_sources])
      #band.snr_overlaps     *= (band.solid_angle/band.solid_angle_p) #small correction for ragged edge
      band.snr_overlaps     *= exposure_ratios

      # make a mask for pixels
      if band.has_pixels:
         if 'pixel_diffs' not in band.__dict__.keys():
            band.pixel_diffs = N.asarray([wsd.difference(roi_dir) for wsd in band.wsdl])
         band.snr_mask = band.pixel_diffs <= radius_in_rad

      # do diffuse counts
      bgs       = back_bgs if band.ct else front_bgs
      ap_evals  = N.empty([nm,ns + 1])      
               
      for ne,e in enumerate(band.bg_points):
         for nbg,bg in enumerate(bgs):
            bg.setEnergy(e)
            ap_evals[nbg,ne]      = SkyIntegrator.ss_average(bg,rd,radius_in_rad)
      ap_evals *= (N.pi*radius_in_rad**2) * band.bg_vector
      
      mo_evals = N.empty([nm,ns + 1])
      for n,m in enumerate(bgm.models):
         mo_evals[n,:] = m(band.bg_points)

      band.snr_bg_counts     = (ap_evals * mo_evals).sum(axis = 1)
      band.snr_bg_all_counts = band.snr_bg_counts.sum()
 
###====================================================================================================###

def snr_profile(roi,radii):

   signals = N.zeros([len(radii),len(roi.energy_bands)])
   backs   = N.zeros_like(signals)
   altlogs  = N.zeros_like(signals)
   nulllogs = N.zeros_like(signals)

   if 'energy_bands' not in roi.__dict__.keys(): roi.setup_energy_bands()

   for nrad,rad in enumerate(radii):
      snr_setup(roi,rad)
      for neb,eb in enumerate(roi.energy_bands):
         for nb,b in enumerate(eb.bands):
            spscounts = b.ps_counts[0] * b.snr_overlaps[0]
            bpscounts = (b.ps_counts[1:] * b.snr_overlaps[1:]).sum()
            dcounts   = b.snr_bg_all_counts
            m         = b.snr_mask

            signals[nrad,neb] += spscounts
            backs[nrad,neb]   += bpscounts + dcounts

            tot_pix = b.ps_all_pix_counts[m] + b.bg_all_pix_counts[m]
            null_tot_pix = b.null_ps_all_pix_counts[m] + b.null_bg_all_pix_counts[m]
            altlogs[nrad,neb]  += (b.pix_counts[m] * N.log(tot_pix)).sum()
            #nulllogs[nrad,neb] += (b.pix_counts[m] * N.log(tot_pix - b.ps_counts[0]*b.ps_pix_counts[:,0][m])).sum()
            nulllogs[nrad,neb] += (b.pix_counts[m] * N.log(null_tot_pix)).sum()
            
   ts = 2* ( altlogs - nulllogs - signals)
   return signals,backs,ts,altlogs,nulllogs

###====================================================================================================###

def TS(self,quick=True,which=0,save_null=False):
   """Calculate the significance of the central point source.
      
      quick -- if set True, just calculate likelihood with source flux set to 0
               if set False, do a full refit of all other free sources
               
      which -- the index of source to calculate -- default to central."""

   if quick:
      self.zero_ps(which)
      ll_0 = self.logLikelihood(self.get_parameters())
      self.unzero_ps(which)
      ll_1 = self.logLikelihood(self.get_parameters())
      return 2*(ll_0 - ll_1)

   save_params = self.parameters().copy() # save free parameters
   self.zero_ps(which)
   ll_0 = self.fit(save_values = False)
   print self
   
   if save_null:
      base_keys = ['ps_all_counts','bg_all_counts']
      for band in self.bands:
         keys = base_keys + (['ps_all_pix_counts','bg_all_pix_counts'] if band.has_pixels else [])
         for key in keys: band.__dict__['null_'+key] = band.__dict__[key]

   self.unzero_ps(which)
   self.set_parameters(save_params) # reset free parameters
   self.__pre_fit__() # restore caching
   ll = -self.logLikelihood(save_params)
   return -2*(ll_0 - ll)

###====================================================================================================###

def read_fit(outfile):
   from cPickle import load
   from psmanager import PointSource
   from skymaps import SkyDir
   ps = []
   d = load(file(outfile))
   for m in d['point_sources']:
      sd = SkyDir(m.ra,m.dec)
      p  = PointSource(sd,m.source_name,m)
      ps.append(p)
   return ps,d['backgrounds']

###====================================================================================================###

class PulsarFitter(object):
   """Provide a master class to manage two ROIAnalysis object to roughly divide
      a light curve into "pulsed" and "unpulsed" emission and do a joint fit.
      
      The first ROI object should contain *two* models for the pulsar, one for
      the pulsed emission and one for the (full-period) DC background.  The
      second ROI object should contain only the DC background.
      
      REDO whole thing to use ROIEnergyBand and possibly to fit the log fluxes in
      the model-independent plot."""

   def __init__(self,roi1,roi2):
      self.roi1,self.roi2 = roi1,roi2
      self.offset1,self.offset2 = 0,0

      self.ll1 = self.roi1.logLikelihood
      self.ll2 = self.roi2.logLikelihood

      self.roi1.setup_energy_bands()
      self.roi2.setup_energy_bands()

      self.on =  N.asarray([ sum( (b.photons for b in eb.bands) ) for eb in roi1.energy_bands])
      self.off = N.asarray([ sum( (b.photons for b in eb.bands) ) for eb in roi2.energy_bands])

      self.phi = self.roi1.phase_factor

      self.on_mask  = self.on  > 0
      self.off_mask = self.off > 0

   def logLikelihood(self,p,*args):
      
      p1 = p
      p2 = N.append(p[:self.offset1],p[self.offset1 + self.offset2:])
      
      ll1 = self.ll1(p1)
      ll2 = self.ll2(p2)

      return ll1 + ll2

   def fit(self,*args,**kwargs):

      # make sure each ROI has the same parameter set
      bm1,bm2 = self.roi1.bgm.models,self.roi2.bgm.models
      for i in xrange(len(bm1)):
         bm2[i].free[:] = bm1[i].free[:]
         bm2[i].p[:]    = bm1[i].p[:]

      pm1,pm2 = self.roi1.psm.models,self.roi2.psm.models
      for n in xrange(len(pm2)):
         pm2[n].free[:] = pm1[n+1].free[:]
         pm2[n].p[:]    = pm1[n+1].p[:]

      self.roi2.__pre_fit__()

      self.offset1 = int(sum([m.free.sum() for m in bm1]))
      self.offset2 = int(pm1[0].free.sum())
      self.roi1.logLikelihood = self.logLikelihood
      self.roi1.fit(*args,**kwargs)
      self.roi1.logLikelihood = self.ll1

   def proflogl(self,p,*args):
      """Calculate the profile likelihood by maximizing the likelihood with
         respect to the background; that is, the background fit is nonparametric."""

      r  = self.roi1
      eb = r.energy_bands
      m  = self.prof_model
      m.set_parameters(p)
      o1  = self.on_mask
      o2  = self.off_mask

      s    = N.asarray( [sum( b.expected(m) * b.overlaps[0] for b in my_eb.bands) for my_eb in eb] )
      s   /= self.phi # gives phase-averaged flux since we are using the full livetime
      b    = s - self.on - self.off
      b    = 0.5 * ((b**2 + 4*self.off*s)**0.5 - b) #positive root

      return (b + s*self.phi).sum() - \
             (self.on[o1]*N.log(b[o1] + s[o1])).sum() - (self.off[o2]*N.log(b[o2])).sum()


   def proffit(self):

      from scipy.optimize import fmin
      self.prof_model = self.roi1.psm.models[0].copy()
      fit = fmin(self.proflogl,self.prof_model.get_parameters(),disp=0,full_output=1,ftol=1e-11,xtol=1e-11)

      from specfitter import SpectralModelFitter as sfm
      h = sfm.hessian(self.prof_model,self.proflogl)
      from numpy.linalg import inv
      self.prof_model.set_cov_matrix(inv(h[0]))
      
   def profplot(self,axes = None, models = [], erange = [1e2,1e4]):

      self.proffit()

      from collections import defaultdict
      from pointlike import DoubleVector

      d = defaultdict(float)
      for b in self.roi1.bands:
         d[b.e] += b.ps_counts[0] * b.overlaps[0]
      modeled = N.asarray([d[e] for e in self.roi1.bin_centers])

      d = defaultdict(float)
      m = self.prof_model
      for b in self.roi1.bands:
         exp = b.expected(m)/m.i_flux(b.emin,b.emax)*b.overlaps[0]
         d[b.e] += exp
      
      exp = N.asarray([d[e] for e in self.roi1.bin_centers])

      w       = self.phi/(1 - self.phi)
      s       = self.on - w * self.off
      b       = self.off/(1-self.phi)
      sigma_s = (self.on + w**2 * self.off)**0.5
      sigma_b = (self.off/(1 - self.phi)**2)**0.5
      
      mask    = s > 0
      dom     = self.roi1.bin_centers[mask]
      delta_e = (self.roi1.bin_edges[1:] - self.roi1.bin_edges[:-1])[mask]
      s       = s[mask]
      sigma_s = sigma_s[mask]
      exp     = exp[mask]

      dom2    = N.logspace(N.log10(erange[0]),N.log10(erange[1]),100)
      
      if axes is None:
         import pylab as P
         axes = P.gca()

      axes.set_xscale('log')
      axes.set_yscale('log')

      y    = dom**2*s/exp/delta_e
      yerr = dom**2*sigma_s/exp/delta_e
      xel  = (self.roi1.bin_centers - self.roi1.bin_edges[:-1])[mask]
      xeh  = (self.roi1.bin_edges[1:] - self.roi1.bin_centers)[mask]
      axes.errorbar(x=dom,xerr=[xel,xeh],y=y,yerr=[N.minimum(0.99*y,yerr),yerr],ls=' ',marker='o',label='Band Fits')
      axes.plot(dom2,dom2**2*m(dom2),label='ON-OFF prof. like. fit')
      
      for n,model in enumerate(models):
         if 'label' not in model.__dict__.keys():
            model.label = 'Additional model %d'%(n+1)
         axes.plot(dom2,dom2**2*model(dom2),label=model.label)

      axes.legend(numpoints=1,loc=0)
      axes.axis([erange[0],erange[1],1e-7,1e-2])
      axes.grid()
      axes.set_xlabel('$\mathrm{Energy\ (MeV)}$')
      axes.set_ylabel('$\mathrm{E^2\ dN/dE\ (MeV\/cm^{-2}\/s^{-1}}$')


"""


def band_fits(self,plot=False,verbose=False):

   from scipy.optimize import newton
   pf  = self.roi.phase_factor
   bands = self.bands

   if plot:
      import pylab as P
      P.figure(4)
      P.gca().set_xscale('log')
      #P.gca().set_yscale('log')

   for i,key in enumerate(self.keys):
      for k,ec in zip([key,-key],[0,1]):

         band = bands[k]

         mask  = ((2*self.u_bins)**0.5 * band.sigma)[1:] <= self.max_radius*N.pi/180
         d  = band.data[mask]
         s  = band.ps_counts[mask] * pf
         b  = (band.bg_counts + band.ne_counts)[mask] * pf
         ns = float(s.sum())
         nb = float(b.sum())

         a0 = ns / (ns + nb)
         
         # fit the signal using extended likelihood

         q  = b / (s/ns)
         f  = lambda a: ( d/(a + q) ).sum() - 1
         fp = lambda a:-( d/(a + q)**2 ).sum()
         sb  = newton(f,ns,fprime=fp)
         dsb = 1./(-fp(sb))**0.5
         asb = sb / (nb + sb)
         dasb = dsb/sb * asb

         # fit the signal fraction using spatial likelihood
         
         q  = b/(s*nb/ns - b) # note both signal and background normalized
         f  = lambda a: ( d/(a + q) ).sum()
         ab  = newton(f,a0,fprime=fp) # estimate from spatial only
         abe = 1./(-fp(ab))**0.5

         if verbose:
            print '\nBand with e: %d and ec: %d'%(key,ec)
            print 'Signal: %.1f (broadband), %.1f +/- %.1f (band only), %.1f +/- %.1f (spatial * counts)'%(ns,sb,dsb,d.sum()*ab,((d.sum()*abe)**2 + d.sum())**0.5)
            print 'Signal Fraction: %.3f (broadband), %.3f +/- %.3f (band only), %.3f +/- %.3f (spatial)'%(a0,asb,dasb,ab,abe)
         
         bd[k]['a0'] = [a0,a0e]
         bd[k]['ab'] = [ab,abe]

         if plot:
            color = 'red' if ec else 'green'
            if i == 0:
               label0 = 'EC %d, Broadband'%(ec)
               label1 = 'EC %d, Band-by-band'%(ec)
            else:
               label0 = None; label1 = None;
            P.errorbar([key],[a0],yerr=[a0e],color=color,marker='x',label=label0)
            P.errorbar([key],[ab],yerr=[abe],color=color,marker='o',label=label1)

         if self.roi_off is not None:
            od = (bd[k]['off_data'])[mask]
            ap = (d.sum() - od.sum()*pf/self.roi_off.phase_factor)/d.sum() # phase estimate
            ape = (1-ap)*(1./d.sum() + 1./od.sum())**0.5
            if verbose: print '     ON/OFF fit      : alpha = %.3f +/- %.3f'%(ap,ape)
            bd[k]['ap'] = [ap,ape]
            if plot:
               color = 'red' if ec else 'green'
               label0 = 'EC %d, ON/OFF'%(ec) if i == 0 else None
               P.errorbar([key],[ap],yerr=[ape],color=color,marker='v',label=label0)
         

def r68(self,plot=False,alpha_mode = 1):

   bbb_fit(self)

   bg_fac = 2
   ae_sys = 0.00

   #ub_r68(self)
   
   bands = self.bands
   pf  = self.roi.phase_factor
   import pylab as P
   if plot:
      P.gca().set_xscale('log')
      P.gca().set_yscale('log')

   for i,key in enumerate(self.keys):
      #if i > 0: continue
      for k,ec in zip([key,-key],[0,1]):
         band = bands[k]
         rad  = ((2*self.u_bins)**0.5 * band.sigma)*180./N.pi
         radm = (rad[1:]*rad[:-1])**0.5
         mask = rad[1:] < self.max_radius

         #akey = ['a0','ab','ap'][alpha_mode]
         #a,ae = bd[k][akey]
         #ae   = (ae**2 + ae_sys**2)**0.5

         d  = band.data[mask]
         n  = d.sum()
         bg = (band.bg_counts + band.ne_counts)[mask] * pf
         bg_hi = bg + bg**0.5
         bg_lo = bg - bg**0.5
         #bg = bg / bg.sum() * (1 - a) * n

         #bg_lo = bg * (1 - a - ae) / (1 - a)
         #bg_hi = bg * (1 - a + ae) / (1 - a)
         
         #s  = a * n
         #s_hi = (a + ae) * n
         #s_lo = (a - ae) * n
         #ds = s * ( (ae/a)**2 + 1./n )**0.5

         #s = (bd[k]['ps_counts'][mask] * pf).sum()
         s = band.ps_total * pf
         s_hi = s + s**0.5
         s_lo = s - s**0.5
         
         integral = N.cumsum(d - bg)/s
         for nv,v in enumerate(integral):
            if v > 0.68: break
         r68 = radm[nv-1] + (0.68 - integral[nv-1])*(radm[nv]-radm[nv-1])/(integral[nv]-integral[nv-1])

         integral = N.cumsum(d - bg_hi)/(s_lo)
         for nv,v in enumerate(integral):
            if v > 0.68: break
         r68_lo = radm[nv-1] + (0.68 - integral[nv-1])*(radm[nv]-radm[nv-1])/(integral[nv]-integral[nv-1])

         integral = N.cumsum(d - bg_lo)/(s_hi)
         for nv,v in enumerate(integral):
            if v > 0.68: break
         r68_hi = radm[nv-1] + (0.68 - integral[nv-1])*(radm[nv]-radm[nv-1])/(integral[nv]-integral[nv-1])

         band.r68 = N.asarray([r68,r68_lo,r68_hi])

         #r68,r68_lo,r68_hi = bd[k]['ub_r68']
         
         if plot:
            l1 = ('Data  Back' if ec else 'Data  Front') if i == 0 else None
            l2 = ('P6_v3 Back' if ec else 'P6_v3 Front') if i == 0 else None
            P.errorbar([key],[r68],yerr=[[r68-r68_lo],[r68_hi-r68]],color='green' if ec else 'red',marker='o',label=l1,capsize=0)
            g = band.gamma
            s = band.sigma
            mc_r68 = (s*180/N.pi) * (2* g*((1-0.68)**(1./(1-g)) - 1))**0.5
            P.plot([key],[mc_r68],color='green' if ec else 'red',marker='^',label=l2)

   if plot:
      P.legend(loc='upper right',numpoints=1)
      P.axis([self.emin,self.emax,0.05,5])
      P.grid()
      P.xlabel('Energy (MeV)')
      P.ylabel('r68 (deg)')



def smearing(self):

   delta1 = 0.067 / 2
   delta2 = 0.300 / 2
   p1     = 0.95
   p2     = 0.05

   def smearprob1(delta):
      return (delta < delta1)/delta1
   def smearprob2(delta):
      return (delta < delta2)/delta2

   bd  = self.bin_dict

   for i,key in enumerate(self.keys):
      for k,ec in zip([key,-key],[0,1]):
         
         rad = ((2*self.u_bins)**0.5 * bd[k]['sigma'])*180./N.pi
         rad = (rad[:-1]*rad[1:])**0.5 # bin centers
         pr1 = smearprob1(rad)
         pr2 = smearprob2(rad)
         ps  = bd[k]['ps_counts']
         nps = N.zeros_like(ps)
         for j,b in enumerate(ps):
            
            if N.all(pr1[j:]<1e-5):
               nps[j]  += p1*ps[j]
            else:
               smear    = pr1[j:] / pr1[j:].sum()
               nps[j:] += p1 * ps[j] * smear
            
            if N.all(pr2[j:]<1e-5):
               nps[j]  += p2*ps[j]
            else:
               smear    = pr2[j:] / pr2[j:].sum()
               nps[j:] += p2 * ps[j] * smear

         
         bd[k]['sps_counts'] = nps

         

"""       
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

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

class PSFBand(object):
   """Container class for data and methods for checking PSF."""

   def init(self):
      self.use_off = True

   def __init__(self,emin,emax,ec,**kwargs):

      self.init()
      self.__dict__.update(kwargs)

      self.emin = emin; self.emax = emax; self.ec = ec

      self.free = N.asarray([True,True,True,True])
      self.p    = N.empty(4)

   def fit(self,psfc,fit_psf=True):

      self.pf = pf = psfc.roi.phase_factor
      
      bins = (2*psfc.u_bins)**0.5 * self.sigma
      mask = bins[1:] < (psfc.max_radius * N.pi / 180)

      self.lefts  = bins[:-1][mask]
      self.rights = bins[1:][mask]
      
      self.bis  = (self.bg_counts[mask] + self.ne_counts[mask])
      self.bis /= self.bis.sum()
      self.counts = self.data[mask]
      if self.use_off:
         self.ocounts= self.off_data[mask]
         self.on2off  = psfc.roi_off.phase_factor/pf

      nb = self.bg_counts[mask].sum() + self.ne_counts[mask].sum()
      ns = self.ps_counts[mask].sum()
      self.p = N.asarray([nb*pf,ns*pf,self.sigma,self.gamma])

      if fit_psf:
         self.free = N.asarray([True]*4)
      else:
         self.free = N.asarray([True]*2 + [False]*2)

      from scipy.optimize import fmin
      self.my_fit = fmin(self.logLikelihood,self.p[self.free],disp=0,full_output=1)

      from specfitter import SpectralModelFitter
      h = SpectralModelFitter.hessian(self,self.logLikelihood)
      from numpy.linalg import inv
      self.cov_matrix = inv(h)

   def scalar_fit(self,psfc):

      mask  = ((2*psfc.u_bins)**0.5 * band.sigma)[1:] <= psfc.max_radius*N.pi/180
      pf    = psfc.roi.phase_factor

      d  = self.data[mask]
      s  = self.ps_counts[mask] * pf
      b  = (self.bg_counts + self.ne_counts)[mask] * pf
      ns = float(s.sum())
      nb = float(b.sum())
      
      # fit the signal using extended likelihood
      q  = b / (s/ns)
      f  = lambda a: ( d/(a + q) ).sum() - 1
      fp = lambda a:-( d/(a + q)**2 ).sum()
      sb  = newton(f,ns,fprime=fp)
      dsb = 1./(-fp(sb))**0.5

      # fit the signal fraction using spatial likelihood      
      q  = b/(s*nb/ns - b) # note both signal and background normalized
      f  = lambda a: ( d/(a + q) ).sum()
      ab  = newton(f,a0,fprime=fp) # estimate from spatial only
      abe = 1./(-fp(ab))**0.5

   def get_parameters(self):
      return self.p[self.free]

   def set_parameters(self,p):
      self.p[self.free ] = p

   def logLikelihood(self,p,*args):
      self.set_parameters(p)
      b,s,sig,gam = self.p

      rus  = 0.5 * (self.rights/sig) ** 2
      lus  = 0.5 * (self.lefts/sig)**2
      fis  = (1 + lus/gam)**(1-gam) - (1 + rus/gam)**(1-gam)
      fmax = fis.sum()

      self.rates = b*self.bis + s*fis

      logl = b + s*fmax - (self.counts*N.log(b*self.bis + s*fis)).sum()
      if self.use_off:
         logl += b*self.on2off- (self.ocounts*N.log(b*self.bis)).sum()

      return logl

   def display(self,ax,psfc,residuals=False,first_plot=False,fit_mode=1):

      #fit_mode: 0 -> use broadband normalization
      #          1 -> use band fits for signal/background normalization
      #          2 -> fit PSF parameters and s/b normalization


      sig = self.sigma
      gam = self.gamma
      ub  = psfc.u_bins
      jac = (2*N.pi*sig**2)*(ub[1:]*ub[:-1])**0.5*N.log(ub[1:]/ub[:-1]) #dn/domega
      pf  = psfc.roi.phase_factor
      mr  = psfc.max_radius

      mask  = ((2*ub)**0.5 *sig)[1:] <= psfc.max_radius*N.pi/180


      def data_plot(dkey,color='black',label=None,pf=1.):
         x = ((ub[1:]*ub[:-1])**0.5)
         y = self.__dict__[dkey]
         hiyerr = y**0.5
         loyerr = N.where(hiyerr >= y,0.99*y,hiyerr)
         ax.errorbar(x = x, y = (y/jac/pf),
                     yerr = [loyerr/jac/pf,hiyerr/jac/pf], ls= ' ', capsize=4,
                     label=label,color=color,marker='o',ms=4)

      def line_plot(ubins,codomain,ec='blue',ls='solid',label=None):

         x = N.zeros(len(ubins) * 2 )
         y = N.zeros(len(ubins) * 2 )
         x[0::2],x[1::2] = ubins,ubins
         y[1:-1:2],y[2::2] = codomain,codomain
         ax.fill(x,y,closed=False,fill=False,edgecolor=ec,ls=ls,label=label)

      def resid_plot(data,model,color='black',label=None):
         x = (ub[1:]*ub[:-1])**0.5
         y = (model-data)/model
         yerr = 1/model**0.5
         ax.errorbar(x = x, y = y, yerr = yerr, ls= ' ', capsize=4,
                     label=label,color=color,marker='o',ms=4)
         ax.axhline(0,color='red')
      
      s_scale,b_scale = 1.,1.
      if fit_mode == 1:
         self.fit(psfc,fit_psf=False) # use fits to renormalize signal and background
         s_scale = self.p[1]/self.ps_counts[mask].sum()/pf
         b_scale = self.p[0]/((self.ne_counts + self.bg_counts)[mask]).sum()/pf

      ps_diffs = s_scale * self.ps_counts / jac
      bg_diffs = b_scale * (self.ne_counts + self.bg_counts) / jac
      to_diffs = bg_diffs + ps_diffs

      if fit_mode == 2:
         self.fit(psfc)
         fsig,fgam = self.p[2:]
         rus  = 0.5 * (self.rights/fsig) ** 2
         lus  = 0.5 * (self.lefts/fsig)**2
         fis  = (1 + lus/fgam)**(1-fgam) - (1 + rus/fgam)**(1-fgam)
         rep_ps = fis * self.p[1] / jac[mask] / pf
         ps_diffs[:len(rep_ps)] = rep_ps
         rep_bg = self.bis * self.p[0] / jac[mask] / pf
         bg_diffs[:len(rep_bg)] = rep_bg
         to_diffs[:len(rep_bg)] = rep_ps + rep_bg


      if not residuals:
         ax.set_yscale('log')
         data_plot('data','black','%d-%d MeV'%(int(self.emin),int(self.emax)),pf)
         if psfc.roi_off is not None:
            data_plot('off_data','red','Offpulse' if first_plot else None,psfc.roi_off.phase_factor)

         line_plot(ub,ps_diffs,'blue','solid',label=None if not first_plot else 'Point Source')
         line_plot(ub,bg_diffs,'red' ,'solid',label=None if not first_plot else 'Background')
         line_plot(ub,to_diffs,'black','dotted',label=None if not first_plot else 'Total Model')

      else:
         resid_plot(self.data/pf,to_diffs*jac,label='%d-%d MeV'%(int(psfc.emin),int(psfc.emax)))


      umin = ub[0]
      rmin = sig*(2*umin)**0.5*180/N.pi
      umax = min(0.5 * (mr*N.pi/180/sig)**2,psfc.umax)
      rmax = min(mr, sig*(2*umax)**0.5*180/N.pi )
      
      mc_r68 = (2* gam*((1-0.68)**(1./(1-gam)) - 1))**0.5
      mc_r95 = (2* gam*((1-0.95)**(1./(1-gam)) - 1))**0.5
      ax.axvline(mc_r68,color='green',ls='-')
      ax.axvline(mc_r95,color='green',ls='-.')

      return [umin,umax,rmin,rmax]



#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


class PSFChecker(object):

   def init(self):

      self.emin = 200
      self.emax = 2e4
      self.ebins_per_decade = 4
      self.umin = 1e-2
      self.umax = 1e3
      self.ubins_per_decade = 4
      self.psf_irf = 'P6_v3_diff'
      self.CALDB   = None
      self.use_lt = True
      self.weight_psf = True
      self.nsimps = 16
      self.roi_off = None
      self.max_radius = 10 # degrees
   
   def __init__(self,roi,**kwargs):

      self.init()
      self.__dict__.update(**kwargs)
      self.roi = roi
      self.sa  = roi.sa

      self.CALDB = self.CALDB or self.sa.CALDB

      lmin,lmax     = N.log10([self.emin,self.emax])
      self.ebands   = N.logspace(lmin,lmax,round( (lmax-lmin)*self.ebins_per_decade ) + 1)
      self.ebands_c = (self.ebands[:-1]*self.ebands[1:])**0.5
      
      self.keys     = self.ebands_c.astype(int)
      self.bands    = dict()

      from psf import PSF
      self.psf = PSF(psf_irf = self.psf_irf, CALDB = self.CALDB, lt = None, sd = None)
      n = len(self.ebands_c)
      self.sigmas = self.psf.sigma(N.append(self.ebands_c,self.ebands_c), N.append([0] * n, [1] * n) ) * (N.pi/180.)
      self.gammas = self.psf.gamma(N.append(self.ebands_c,self.ebands_c), N.append([0] * n, [1] * n) )
      
      if self.use_lt:
         # use exposure weighting for PSF -- this is optimal for a single source analysis
         self.psf = PSF(psf_irf = self.psf_irf, CALDB = self.CALDB, lt = self.sa.livetimefile, sd = self.sa.roi_dir)
         sigmas = self.psf.sigma(N.append(self.ebands_c,self.ebands_c), N.append([0] * n, [1] * n) ) * (N.pi/180.)
         gammas = self.psf.gamma(N.append(self.ebands_c,self.ebands_c), N.append([0] * n, [1] * n) )

         print 'Difference from using exposure instead of Aeff:'
         print 'Sigmas:'
         print (sigmas - self.sigmas)/self.sigmas
         print 'Gammas:'
         print (gammas - self.gammas)/self.gammas

         self.sigmas = sigmas
         self.gammas = gammas

      if self.weight_psf:
         m   = self.roi.psm.models[0]
         exp = self.roi.sa.exposure.exposure
         sd  = self.roi.sa.roi_dir

         sigmas = N.empty( 2 * n )
         gammas = N.empty( 2 * n )
         for ne,(emin,emax)in enumerate(zip(self.ebands[:-1],self.ebands[1:])):
            sigmas[ne]     = m.expected(emin,emax,exp,sd,0,lambda e: self.psf.sigma(e,[0]*len(e)))
            sigmas[ne + n] = m.expected(emin,emax,exp,sd,1,lambda e: self.psf.sigma(e,[1]*len(e)))
            gammas[ne]     = m.expected(emin,emax,exp,sd,0,lambda e: self.psf.gamma(e,[0]*len(e)))
            gammas[ne + n] = m.expected(emin,emax,exp,sd,1,lambda e: self.psf.gamma(e,[1]*len(e)))
         
         sigmas *= (N.pi/180)

         print 'Comparing weighted versus standard parameters'
         print 'Sigmas:'
         print (sigmas - self.sigmas)/self.sigmas
         print 'Gammas:'
         print (gammas - self.gammas)/self.gammas

         self.sigmas = sigmas
         self.gammas = gammas

      for i,(emin,emax,key) in enumerate(zip(self.ebands[:-1],self.ebands[1:],self.keys)):
         
         fband = PSFBand(emin,emax,0)
         fband.sigma = self.sigmas[i]
         fband.gamma = self.gammas[i]

         bband = PSFBand(emin,emax,1)
         bband.sigma = self.sigmas[i+n]
         bband.gamma = self.gammas[i+n]

         self.bands[key]  = fband
         self.bands[-key] = bband
      
      lumin,lumax  = N.log10([self.umin,self.umax])
      self.u_bins  = N.logspace(lumin,lumax,round( (lumax-lumin)*self.ubins_per_decade ) + 1)

      self.bindata(self.sa)
      self.ps_models()
      self.bg_models()
      if self.roi_off is not None:
         self.bindata(self.roi_off.sa,dkey='off_data')
        
   
   def bindata(self,sa,dkey='data'):

      from fitstools import rad_extract
      max_rad   = self.psf.sigma(self.ebands_c[0],1)*(2*self.umax)**0.5 #rad_extract expects degrees
      data_dict = rad_extract(sa.pixeldata.event_files, sa.roi_dir,
                              max_rad, return_cols=['ENERGY','EVENT_CLASS'])
      
      diffs = data_dict['DIFFERENCES']
      ens   = data_dict['ENERGY']
      ecs   = data_dict['EVENT_CLASS']

      bands = self.bands

      for i,key in enumerate(self.keys):

         emin,emax = self.ebands[i],self.ebands[i+1]
         mask      = N.logical_and(ens >= emin, ens < emax)
         
         for k,ec in zip([key,-key],[0,1]):

            us    = 0.5 * (diffs[N.logical_and(mask,ecs==ec)] / bands[k].sigma) ** 2
            bands[k].__dict__[dkey+'_us'] = us #save differences for possible unbinned use
            bands[k].__dict__[dkey]       = N.histogram(us,bins=self.u_bins,new=True)[0]

   
   def ps_models(self):
      
      psm = self.roi.psm
      exp = self.sa.exposure.exposure
      #need to calculate the overlap for all bin radii and all sources

      from roi_modules import ROIOverlap

      f  = self.psf.overlap
      dirs = [ps.skydir for ps in psm.point_sources]
      d1   = dirs[0]
      bands = self.bands

      for i,key in enumerate(self.keys):

         emin,emax = self.ebands[i:i+2]

         for k,ec in zip([key,-key],[0,1]):

            band = bands[k]

            radii = ((2*self.u_bins)**0.5 * bands[k].sigma)
            band.ne_counts = N.zeros_like(radii)

            for n,ps in enumerate(psm.point_sources):
               fracs = [float(f(d1,ps.skydir,self.ebands_c[i],ec,r)) for r in radii]
               expec = ps.model.expected(emin,emax,exp,ps.skydir,event_class=ec)
               if n == 0:
                  band.ps_counts = N.asarray(fracs) * expec
                  band.ps_total  = expec
               else:
                  band.ne_counts += N.asarray(fracs) * expec

            #adjust to annular counts
            band.ps_counts = band.ps_counts[1:] - band.ps_counts[:-1]
            band.ne_counts = band.ne_counts[1:] - band.ne_counts[:-1]

  
   def bg_models(self):

      bg  = self.roi.bgm
      sa  = self.sa
      rd  = self.sa.roi_dir
      exp = sa.exposure.exposure
      from skymaps import Background,SkyIntegrator

      simps_vec = (N.asarray([1,4,2,4,2,4,2,4,1]))/(3.*8.)

      e_models = [m.smodel for m in bg.bgmodels]
      b_models = [ [Background(m.dmodel[i],exp[i]) for i in xrange(len(m.dmodel))] for m in bg.bgmodels]

      
      for i,key in enumerate(self.keys):

         emin, emax  = self.ebands[i:i+2]
         lemin,lemax = N.log10([emin,emax])
         simps_pts   = N.logspace(lemin,lemax,9)

         for k,ec in zip([key,-key],[0,1]):

            band = self.bands[k]

            radii = ((2*self.u_bins)**0.5 * band.sigma)
            band.bg_counts = N.zeros_like(radii)

            for nrad,rad in enumerate(radii):

               solid = N.pi * rad**2

               for n in xrange(len(e_models)):
                  
                  em = e_models[n]
                  bm = b_models[n]
                  bm = bm[ec if len(bm) > 1 else 0]
                  
                  bm_pts = N.empty_like(simps_pts)
                  for ne,e in enumerate(simps_pts):
                     bm.setEnergy(e)
                     bm_pts[ne] = SkyIntegrator.ss_average(bm,rd,rad)

                  band.bg_counts[nrad] += \
                     (bm_pts * em(simps_pts) *simps_pts * simps_vec).sum() * (N.log(emax/emin) * solid )
                     
            band.bg_counts = band.bg_counts[1:] - band.bg_counts[:-1]

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def ub_r68(self):

   bands = self.bands
   w     = self.roi.phase_factor/self.roi_off.phase_factor
   from scipy.optimize import newton

   for i,key in enumerate(self.keys):
      for k,ec in zip([key,-key],[0,1]):
         
         band = bands[k]
         umax = self.u_bins[-1]

         d1 = band.data_us
         d2 = band.off_data_us

         d1 = d1[d1 < umax]; d2 = d2[d2 < umax]
         d1.sort();d2.sort()
         nd1 = len(d1); nd2 = len(d2);
         nd2p = w * nd2
         ns  = nd1 - nd2p

         f = lambda sig: nd1 * (N.log( nd1/(sig + nd2p) ) - 1) + (sig + nd2p) - 2.71/2
         ns_hi = newton(f,ns + ns**0.5)
         ns_lo = newton(f,max(0,ns - ns**0.5))

         n68    = ns    * 0.68
         n68_lo = ns_lo * 0.68
         n68_hi = ns_hi * 0.68

         u68 = 0; u68_lo = 0; u68_hi = 0;
         for nu,u in enumerate(d1):
            nback = N.searchsorted(d2,u)
            if (not u68) and (nu - w*nback >= n68):
               u68 = u
            if (not u68_lo) and (nu - w*nback >= n68_lo):
               u68_lo = u
            if (not u68_hi) and (nu - w*nback >= n68_hi):
               u68_hi = u
            if (u68 and u68_lo and u68_hi): break

         band.ub_r68 = (2*N.asarray([u68,u68_lo,u68_hi]))**0.5*band.sigma*(180/N.pi)

def show(self,residuals=False,outfilestem=None,fit_mode=1):

   import pylab as P
   P.rcParams['legend.fontsize']='small'
   from math import floor

   ub = self.u_bins
   pf = self.roi.phase_factor

   first_plot = True

   for i,key in enumerate(self.keys):

      emin,emax = self.ebands[i:i+2]

      #if i < 4 or i > 7: continue # show subset of panels
      if i > 7: continue

      for k,fignum in zip([key,-key],[5,6]):

         band = self.bands[k]

         P.figure(fignum,figsize=(18,12))
         ax = P.subplot( 2, 4, i + 1)
         left_plot = (i + 1) % 4 == 1
         top_plot  = (i + 1) <= 4
         bot_plot  = (i + 1) > 1 * 4
         #ax = P.subplot( 2, 2, i - 3)
         ax.set_xscale('log')

         if not residuals:
            ax.set_yscale('log')
            ymin = 2e3
            ymax = 2e8

         else:
            ymin = -0.5
            ymax = 0.5

         umin,umax,rmin,rmax = band.display(ax,self,
                                            residuals  = residuals,
                                            first_plot = first_plot,
                                            fit_mode   = fit_mode,)

         if bot_plot:
            ax.set_xlabel('u (bottom axis)',size='large')
         if left_plot:
            ax.set_ylabel('$\mathrm{dN/d\Omega\ (ph\/sr^{-1})}$',size='large')
         ax.legend(loc='lower left' if first_plot else 'upper right',numpoints=1)

         ax2 = ax.twiny()
         ax2.set_xscale('log')
         if top_plot:
            ax2.set_xlabel('degrees (top axis)',size='large')

         ax.axis([umin,umax,ymin,ymax])
         ax2.axis([rmin,rmax,ymin,ymax])

         if first_plot:
            if not residuals:
               P.suptitle('Radial Counts for %s'%('Front' if k > 0 else 'Back'),size='large')
            else:
               P.suptitle('Fractional Residuals for %s'%('Front' if k > 0 else 'Back'),size='large')
         

      first_plot = False

   if outfilestem is not None:
      P.figure(5)
      P.savefig(outfilestem+'_ec0.png')
      P.figure(6)
      P.savefig(outfilestem+'_ec1.png')


def all_panels(self,outfilestem,fit_mode=1):
   import pylab as P
   P.close('all')
   P.ioff()
   show(self,residuals = False, outfilestem = outfilestem+'_density_', fit_mode = fit_mode)
   P.close('all')
   P.ioff()
   show(self,residuals = True, outfilestem = outfilestem+'_residuals_', fit_mode = fit_mode)

def show_PSF_params(self):

   ub_r68(self)

   import pylab as P
   P.ioff()
   P.close('all')
   P.figure(7)
   ax1 = P.gca()
   P.figure(8)
   ax2 = P.gca()
   P.figure(9)
   ax3 = P.gca()

   for ax in [ax1,ax2,ax3]: ax.set_xscale('log')
   for ax in [ax1,ax3]: ax.set_yscale('log')

   for i,key in enumerate(self.keys):
      for k,ec in zip([key,-key],[0,1]):

         band = self.bands[k]
         band.fit(self,fit_psf = True)

         sig,gam = band.p[2:]
         sige,game = (N.diag(band.cov_matrix)**0.5)[2:]

         sig  *= (180/N.pi)
         sige *= (180/N.pi)

         mcsig = band.sigma * (180/N.pi)
         mcgam = band.gamma

         r68   = sig * (2*gam* ( (1-0.68)**(1./(1-gam)) - 1) )**0.5
         mcr68 = mcsig * (2*mcgam* ( (1-0.68)**(1./(1-mcgam)) - 1) )**0.5

         l1 = '%s: Vela Fit'%('Back' if ec else 'Front') if i == 0 else None
         l2 = '%s: P6_v3'%('Back' if ec else 'Front') if i == 0 else None
         l4 = '%s: Unbinned'%('Back' if ec else 'Front') if i == 0 else None
         
         ax1.errorbar([key],[sig],yerr = [sige],color='red' if ec else 'green',label = l1,marker='^')
         ax1.plot([key],[mcsig],color='red' if ec else 'green',label = l2,marker='o',ls=' ')

         ax2.errorbar([key],[1./gam],yerr = [game/gam**2],color='red' if ec else 'green',label = l1,marker='^')
         ax2.plot([key],[1./mcgam],color='red' if ec else 'green',label = l2,marker='o',ls=' ')

         ax3.plot([key],[r68],color='red' if ec else 'green',label = l1,marker='^',ls=' ')
         ax3.plot([key],[mcr68],color='red' if ec else 'green',label = l2,marker='o',ls=' ')
         ax3.plot([key],[band.ub_r68[0]],color='red' if ec else 'green',label = l4,marker='x',ls=' ')

   
   emin,emax,r68f,r68ferr,r68b,r68berr,sf,sferr,gf,gferr,sb,sberr,gb,gberr = parse_toby_psf()
   l3 = '%s: Toby Fit'%('Back' if ec else 'Front') if i == 0 else None
   for i,(e1,e2) in enumerate(zip(emin,emax)):
      for ec in [0,1]:
         tkey = (e1*e2)**0.5*1000
         l3 = '%s: Toby Fit'%('Back' if ec else 'Front') if i == 0 else None
         tsig,tsige = (sb[i],sberr[i]) if ec else (sf[i],sferr[i])
         tgam,tgame = (gb[i],gberr[i]) if ec else (gf[i],gferr[i])
         tr68,tr68e = (r68b[i],r68berr[i]) if ec else (r68f[i],r68ferr[i])
         ax1.errorbar([tkey],[tsig],yerr = [tsige],color='red' if ec else 'green',label = l3,marker='s')
         ax2.errorbar([tkey],[1./tgam],yerr = [tgame/tgam**2],color='red' if ec else 'green',label = l3,marker='s')
         ax3.plot([tkey],[tr68],color='red' if ec else 'green',label = l3,marker='s',ls=' ')

   ax1.set_xlabel('Energy (MeV)')
   ax1.set_ylabel('sigma (deg)')
   ax1.legend(numpoints=1,loc='lower left')
   ax1.axis([1e2,1e5,1e-2,2])
   ax1.grid()

   ax2.set_xlabel('Energy (MeV)')
   ax2.set_ylabel('1/gamma')
   ax2.legend(numpoints=1,loc='lower left')
   ax2.axis([1e2,1e5,0,1])
   ax2.grid()

   ax3.set_xlabel('Energy (MeV)')
   ax3.set_ylabel('r68 (deg)')
   ax3.legend(numpoints=1,loc='lower left')
   ax3.axis([1e2,1e5,5e-2,5])
   ax3.grid()

def parse_toby_psf(filename = r'd:\users\kerrm\python\analyses\vela2\toby_psf.txt'):
   toks = N.asarray([line.strip().split() for line in file(filename)][1:]).astype(float)
   emin,emax = toks[:,0],toks[:,1]
   r68f,r68ferr = toks[:,2],toks[:,3]
   r68b,r68berr = toks[:,4],toks[:,5]
   sf,sferr = toks[:,6],toks[:,7]
   gf,gferr = toks[:,8],toks[:,9]
   sb,sberr = toks[:,10],toks[:,11]
   gb,gberr = toks[:,12],toks[:,13]

   return emin,emax,r68f,r68ferr,r68b,r68berr,sf,sferr,gf,gferr,sb,sberr,gb,gberr

def comp_spec_fits(self):

   nband = len(self.bands)
   bands = self.bands

   from collections import deque
   p6_signal = deque()
   fi_signal = deque()

   for i,key in enumerate(self.keys):
      for k,ec in zip([key,-key],[0,1]):

         band = bands[k]
         band.fit(self,fit_psf=False)
         p6_signal.append(band.p[1])
         band.fit(self,fit_psf=True)
         fi_signal.append(band.p[1])

   p6_signal = N.asarray(p6_signal)
   fi_signal = N.asarray(fi_signal)

   from Models import ExpCutoff
   mp6 = ExpCutoff()
   mfi = ExpCutoff()

   exp = self.roi.sa.exposure.exposure
   sd  = self.roi.psm.point_sources[0].skydir
   
   def logLikelihood(p,m,*args):
      signal = args[0]
      m.set_parameters(p)
      expected = deque()
      for i,key in enumerate(self.keys):
         for k,ec in zip([key,-key],[0,1]):
            band = self.bands[k]
            expected.append(m.expected(band.emin,band.emax,exp,sd,event_class=ec))
      
      expected = N.asarray(expected)
      ll = (expected - signal * N.log(expected)).sum()
      return ll

   from scipy.optimize import fmin
   fit1 = fmin(logLikelihood,mp6.get_parameters(),args=(mp6,p6_signal))
   fit2 = fmin(logLikelihood,mfi.get_parameters(),args=(mfi,fi_signal))

   from specfitter import SpectralModelFitter
   h1 = SpectralModelFitter.hessian(mp6,logLikelihood,p6_signal)
   h2 = SpectralModelFitter.hessian(mfi,logLikelihood,fi_signal)

   from numpy.linalg import inv
   mp6.set_cov_matrix(inv(h1))
   mfi.set_cov_matrix(inv(h2))

   return mp6,mfi

###====================================================================================================###

class NSMap(object):
   """Attempt to add and fit a new source to an existing model and display the
      results as a TS map."""

   def __init__(self,roi):
      self.roi = roi
      self.phase_factor = roi.phase_factor
      self.ro = ROIOverlap()
      self.mo = PowerLaw(free=[True,False])
      self.ll = roi.logLikelihood(roi.parameters())

   def __call__(self,v):
      from skymaps import SkyDir,Hep3Vector,PsfSkyFunction
      skydir = SkyDir(Hep3Vector(v[0],v[1],v[2]))
      #skydir = v
      ro = self.ro
      rd = self.roi.sa.roi_dir

      self.mo.p[0] = -11.5

      # calculate new spatial structure for source
      for i,band in enumerate(self.roi.bands):

         sigma,gamma,en,exp,pa = band.s,band.g,band.e,band.exp,band.b.pixelArea()
         exposure_ratio        = exp.value(skydir,en)/exp.value(rd,en) #needs fix? -- does not obey which?
         psf                   = PsfSkyFunction(skydir,gamma,sigma)
         
         band.ns_overlap       = ro(band,rd,skydir) * exposure_ratio * band.solid_angle / band.solid_angle_p

         if band.has_pixels:
            band.ns_ps_pix_counts = N.asarray(psf.wsdl_vector_value(band.wsdl))*((pa/(2*N.pi*sigma**2))*exposure_ratio)

         band.base_counts = band.expected(self.mo)
      from scipy.optimize import fmin
      f = fmin(self.logLikelihood,N.asarray([-11.5]),full_output=1,disp=0)
      return 2*(self.ll - f[1])

   def logLikelihood(self,parameters,*args):

      bands = self.roi.bands
      ll    = 0
      mo    = self.mo
      mo.set_parameters(parameters)

      for b in bands:

         new_ps_counts = b.base_counts * 10**(mo.p[0] + 11.5)


         ll +=  ( 
                   #integral terms for ROI (go in positive)
                   (b.bg_all_counts + b.ps_all_counts + new_ps_counts*b.ns_overlap)*self.phase_factor

                   -

                   #pixelized terms (go in negative)
                   (b.pix_counts *
                       N.log(  (b.bg_all_pix_counts + b.ps_all_pix_counts + new_ps_counts*b.ns_ps_pix_counts)*self.phase_factor )
                   ).sum() if b.has_pixels else 0.
                )

      return 1e6 if N.isnan(ll) else ll

