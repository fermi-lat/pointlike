
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

def snr(roi,aperture=1):
   #could ultimately use an arbitrary function of energy to calculate SNR; for now, fixed

   aperture *= (N.pi/180.)

   ps_counts = 0
   bg_counts = 0
   
   psm = roi.ps_manager
   sd  = psm.ROI_dir()
   from roi_modules import ROIOverlap
   ro = ROIOverlap()

   for slw in roi.fitter.pslw:
      for i,ps in enumerate(psm.point_sources):
         contrib = ro(slw,sd,ps.skydir,fixed_radius=aperture) * slw.ps_counts[i]
         ps_counts += (contrib if i == 0 else 0)
         bg_counts += (contrib if i >  0 else 0)

   from skymaps import Background,BandBackground

   bgm = roi.bg_manager
   for bgt in ['galactic','isotropic']:

      if bgm.__dict__['use_%s'%bgt]:

         sa = bgm.spectral_analysis
         bg = Background(sa.background.__dict__['%s_diffuse'%bgt],sa.exposure.exposure[0],sa.exposure.exposure[1])
         sh = bgt[:3]

         for slw in roi.fitter.pslw:
            band            = slw.sl.band()
            band_background = BandBackground(bg,band)
            counts          = band_background.average(bgm.roi_dir,aperture,0.01) * N.pi*(aperture)**2
            bg_counts       += counts/slw.__dict__[sh+'_ic']*slw.__dict__[sh+'_exp']

   return ps_counts,bg_counts,ps_counts/(ps_counts + bg_counts)**0.5

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
      self.prof_model.set_cov_matrix(inv(h))
      
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
      b_models = [ [Background(m.dmodel[i],exp[i]) for i in xrange(len(m.dmodel))] for m in bg.bgmodels]

      
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

                  #bm.set_event_class(ec)
                  bm = bm[ec if len(bm) > 1 else 0]
                  

                  bm_pts = N.empty_like(simps_pts)
                  for ne,e in enumerate(simps_pts):
                     bm.setEnergy(e)
                     bm_pts[ne] = SkyIntegrator.ss_average(bm,rd,rad)

                  self.bin_dict[k]['bg_counts'][nrad] += \
                     (bm_pts * em(simps_pts) *simps_pts * simps_vec).sum() * (N.log(emax/emin) * solid )
                     
            self.bin_dict[k]['bg_counts'] = self.bin_dict[k]['bg_counts'][1:] - self.bin_dict[k]['bg_counts'][:-1]


         
   def show(self):

      pf = self.roi.phase_factor

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
            ps_diffs = self.bin_dict[k]['ps_counts'] * pf
            ne_diffs = self.bin_dict[k]['ne_counts'] * pf
            bg_diffs = self.bin_dict[k]['bg_counts'] * pf
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
            if i == 0:
               P.suptitle('Radial Counts for Event Class: %d, Black: Data, Red: Total Model, Blue: Point Source'%(0 if k > 0 else 1),size='large')

            ax2 = ax.twiny()
            rad = self.bin_dict[k]['sigma']*(2*self.u_bins)**0.5*180/N.pi
            ax2.axis([rad[0],rad[-1],ymin,a[3]])

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
              


