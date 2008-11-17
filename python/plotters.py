"""
A module to generate plots for analysis of point sources.  The input for all classes below is an
instance of PointSourceLikelihoodWrapper from which information on the photon distribution and the
spectral fits is obtained.  Spectral models can be furnished separately to provide an independent
estimate of the source counts and for display purposes.

$Header$
"""

import numpy as N

class GraphicAnalysis(object):
   """Parent class for studies of the angular distribution of photons."""

   def init(self):
      pass

   def __init__(self,pslw,**kwargs):
      """pslw - an instance of PointSourceLikelihoodWrapper."""
      self.init()
      self.__dict__.update(**kwargs)
      self.pslw = pslw

      import pylab
      self.pylab = pylab

class TwoDAnalysis(GraphicAnalysis):
   """Analyze the photon distribution in two dimensions: (RA,DEC) or (L,B).
   
Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  scale       [8] image extent in degrees                                                                    
  resolution  [100] pixels per axis                                                                            
  equatorial  [True] boolean; True --> equatorial, False --> galactic                                                                  
  model       [None] an optional spectral model to use to predict source counts                                                                
  relative    [False] mode for residuals; absolute or relative (scaled by model counts)                                                                    
  units       [1e6] scale factor for photon density                                                
  =========   =======================================================
  """

   def init(self):
      """Set default values for keyword arguments and internal variables."""
      key = ['scale','resolution','equatorial','model','relative','units']
      val = [8,100,True,None,False,1e6]      
      for i in xrange(len(key)): self.__dict__[key[i]] = val[i]
      self._images = [None]*5
      
   def __make_image__(self,mode):
      """Generate image with appropriate call."""
      if mode == 3 or mode == 4:
         for i in [0,1,2]: self.__make_image__(i)
         if mode == 3: self._images[3] = self._images[0] + self._images[2]
         if mode == 4:
            self.__make_image__(mode=3)
            self._images[4] = self._images[1] - self._images[3]
            if self.relative:
               self._images[4] /= self._images[3]
               self._images[4][N.abs(self._images[4]) > 1] = 1
      if type(self._images[mode]) is not type(None): return

      from skymaps import SkyDir
      grid=N.linspace(-self.scale, self.scale, self.resolution)
      d = self.pslw.psl.dir()
      coords = d.ra(),d.dec() if self.equatorial else d.l(),d.b()
      flag = SkyDir.EQUATORIAL if self.equatorial else SkyDir.GALACTIC
      coslat = N.cos(coords[1]*N.pi/180.)
      self.corners = [coords[0] - self.scale/coslat,coords[0] + self.scale/coslat, coords[1] - self.scale, coords[1] + self.scale]
      self._images[mode] = N.array([self.pslw.display(SkyDir(coords[0]-dx/coslat, coords[1]-dy,flag),\
                          mode=mode,model=self.model) for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
      self._images[mode]/=self.units

   def __label__(self,mode):
      """Label the plot."""
      if self.equatorial: self.pylab.xlabel('Right Ascension');self.pylab.ylabel('Declination')
      else: self.pylab.xlabel('Galactic Longitude'); self.pylab.ylabel('Galactic Latitude')
      clabel = 'Photon Density / $10^6$'
      if mode == 4 and self.relative: clabel = '1-obs/model'
      c=self.pylab.colorbar(orientation='horizontal')
      c.set_label(clabel)
      self.pylab.title(['Source Prediction','Counts','Background Prediction', 'Total Prediction', 'Residuals'][mode])

   def show(self, mode = 4, **kwargs):
      """Display various spatial plots with mode defined as follows:
         mode == 0: point source expectation
         mode == 1: observed values
         mode == 2: background prediction
         mode == 3: total prediction
         mode == 4: residuals.
         """
      self.__make_image__(mode)
      self.pylab.imshow(self._images[mode], extent=self.corners,interpolation='nearest', **kwargs)
      self.__label__(mode)     

   def grand_slam(self, **kwargs):
      """Convenience method for a four-panel plot with source, background, total, and residuals."""
      for i,mode in enumerate([0,2,1,4]):
         self.pylab.subplot(2,2,i+1)
         self.show(mode=mode)

   def histogram(self, mode = 4, bins=20):
      """Distribution of image pixels in given mode."""
      pixels = N.ravel(self.image)
      self.pylab.hist(pixels,bins=bins)
      self.pylab.axvline(0, color='white')
      self.pylab.axvline(pixels.mean(), color='red',label='Mean')
      self.pylab.legend()

class OneDAnalysis(GraphicAnalysis):
   """Analyze the photon distribution integrated over azimuth.

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  bins        [20] bins in angular deviation                                                                  
  density     [True] if True, divide counts by annulus solid angle; makes isotropic background a constant                                                                           
  normalized  [True] if True, use absolute normalization of background; otherwise, use observed counts                                                                  
  model       [None] an optional spectral model to use to predict source counts                                                                                                        
  =========   =======================================================
  """
   def init(self):
      self.bins = 20
      self.density = True
      self.normalized = True
      self.model = None
      self._edges = None

   def __calc_counts__(self):
      from skymaps import WeightedSkyDirList
      src_dir = self.pslw.psl.dir()
      bands = [slw.sl.band() for slw in self.pslw]
      wsdl  = [WeightedSkyDirList(b,src_dir,self.pslw.psl.maxROI()/180.*N.pi) for b in bands]
      diffs,weights = N.asarray([[src_dir.difference(x),x.weight()] for y in wsdl for x in y]).transpose()
      diffs*= (180./N.pi)
      bins = N.array([x**2 for x in xrange(1,self.bins+1)])/float(self.bins**2)*self.pslw.psl.maxROI()
      bins = N.append([0],bins[bins>=0.1])
      self.counts,self._edges = N.histogram(diffs,weights=weights,bins=bins,new=True)

   def __calc_expected__(self):
      src_dir = self.pslw.psl.dir()
      radii   = self._edges/180.*N.pi
      disc_solid_angles = 2*N.pi*N.array([1-N.cos(r) for r in radii])
      self.annuli_solid_angles = (disc_solid_angles[1:]-disc_solid_angles[:-1])
      self.bg_predicted  = N.zeros_like(self.counts)
      self.src_predicted = N.zeros_like(self.bg_predicted)

      for i,slw in enumerate(self.pslw):
         sl = slw.sl
         background = sl.background_function()
         sigma,gamma = sl.band().sigma(),sl.band().gamma()
         alpha = min(sl.alpha(),1-1e-5)
         u = 0.5*(radii/sigma)**2

         #Average background over annuli -- painfully slow
         prev = background.average(src_dir,radii[0],0.05) if self._edges[0] > 0 else 0
         for j in xrange(len(radii)-1):
            disc = background.average(src_dir,radii[j+1],0.05)
            sa1,sa2 = disc_solid_angles[j],disc_solid_angles[j+1]
            self.bg_predicted[j] += (sa2*disc - sa1*prev)
            prev = disc

         #Calculate source counts in each bin by integrating PSF and multiplying by normalization
         if alpha < 0: continue
         
         integ_psf = (1+u[:-1]/gamma)**(1-gamma) - (1+u[1:]/gamma)**(1-gamma)
         fmax = 1- (1+u[-1]/gamma)**(1-gamma)
         norm = slw.expected(model) if self.model is not None else slw.signal()[0]
         self.src_predicted += norm/fmax*integ_psf

   def __set_axis__(self,mi,ma):
      """Set axis appropriately, with a minimum of 1 decade."""
      diff = N.log10(ma) - N.log10(mi)
      if diff < 1:
         ma = 1.3*10**(N.log10(ma)+(1-diff)/2.)
         mi = 0.7*10**(N.log10(mi)-(1-diff)/2.)
      self.pylab.axis([self._edges[0],self._edges[-1],mi,ma])

   def __label__(self):
      """Label plot."""
      self.pylab.legend(loc='upper right')
      self.pylab.grid(b=True)
      self.pylab.xlabel('Angular Deviation (deg)')
      if self.density: self.pylab.ylabel('Photon Density')
      else: self.pylab.ylabel('Counts')

   def show(self,**kwargs):
      """Display azimuthally-integrated observed and predicted counts."""
      
      if self._edges is None: self.__calc_counts__();self.__calc_expected__()
      x,obs,bg,src,sa = (self._edges[1:]+self._edges[:-1])/2.,self.counts,self.bg_predicted,self.src_predicted,self.annuli_solid_angles
      tot = bg+src
      if self.density: self.pylab.gca().set_yscale('log')
      else: sa = 1 #solid angle to transform counts -> dn/domega

      self.pylab.errorbar(x,y=obs/sa,yerr=(tot**0.5)/sa,ls=' ',marker='o',color='k',mfc='white',label='Observed')
      self.pylab.plot(x,tot/sa,color='k',label='Model')
      
      self.__set_axis__(min(obs/sa)*0.7,max(obs/sa)*1.3)
      self.__label__()

class ScatterAnalysis(GraphicAnalysis):
   """Look at invidual photon locations.  Only useful for a small energy slice or plot is too confused.

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  emin        [2000] minimum energy to include                                                                 
  emax        [1e6] maximum energy to include                                                                          
  equatorial  [True] if True, use equatorial, otherwise, galactic                                                             
  roi_factor  [1] multiplier for ROI; get photons from roi_factor*ROI
  conf_radius [None] if set to a value between 0 and 1, set ROI to confidence radius; otherwise, pointlike default
  draw_rois   [True] if True, draw ROI boundary on plot
  event_class [-1] select events of event_class only
  =========   =======================================================
  """

   def init(self):
      key = ['emin','emax','equatorial','roi_factor','conf_radius','draw_rois','event_class']
      val = [2000,1e6,True,1,None,True,-1]
      for i in xrange(len(key)): self.__dict__[key[i]] = val[i]
      self._wsdl = None
      
   def __calc_counts__(self):
      from skymaps import WeightedSkyDirList
      if self.conf_radius is not None and (self.conf_radius > 1 or self.conf_radius <= 0): self.conf_radius = None
      src_dir = self.pslw.psl.dir()
      self.rois,self._wsdl = [],[]
      for slw in self.pslw:
         sl,b = slw.sl,slw.sl.band()
         if not (b.emin() >= self.emin and b.emax() <= self.emax and \
                (self.event_class < 0 or b.event_class() == self.event_class)): continue
         u = ((1-self.conf_radius)**(1./(1-sl.gamma()))-1)*sl.gamma() if self.conf_radius is not None else sl.umax()
         self.rois += [self.roi_factor*(2*u)**0.5*sl.sigma()]
         self._wsdl += [WeightedSkyDirList(b,src_dir,self.rois[-1])]
         self._wsdl[-1].energy = (b.emin()*b.emax())**0.5

   def __process_dirs__(self):
      lon_func,lat_func = ('wsd.ra','wsd.dec') if self.equatorial else ('wsd.l','wsd.b')
      exec('lon  = [%s() for wsdl in self._wsdl for wsd in wsdl ]'%lon_func)
      exec('lat  = [%s() for wsdl in self._wsdl for wsd in wsdl ]'%lat_func)
      ens = [wsdl.energy for wsdl in self._wsdl for wsd in wsdl]
      self.single_energy = min(ens)==max(ens)
      return lon,lat,ens

   def __draw_rois__(self,s):
      d = self.pslw.psl.dir
      center = (d().ra(),d().dec()) if self.equatorial else (d().l(),d().b())
      coslat = N.cos(center[1]/180.*N.pi)
      self.pylab.axvline(center[0],color='k');self.pylab.axhline(center[1],color='k')
      
      scale = max(self.rois)*180/N.pi
      self.pylab.axis([center[0]-scale/coslat,center[0]+scale/coslat,center[1]-scale,center[1]+scale])

      from matplotlib.patches import Ellipse
      for i,roi in enumerate(self.rois):
         radius = 2*roi*180/N.pi/self.roi_factor #Ellipse convention is the major axis, not semi-major
         edgecolor = 'k' if self.single_energy else s.cmap( s.norm( N.log10(self._wsdl[i].energy))) 
         e = Ellipse(center,width=radius/coslat,height=radius,ec=edgecolor,fill=False)
         self.pylab.gca().add_artist(e)

   def __label__(self,s,c):
      self.pylab.xlabel('Right Ascension' if self.equatorial else 'Galactic Longitude')
      self.pylab.ylabel('Declination' if self.equatorial else 'Galactic Latitude')
      if c is not None: c.set_label('$\mathrm{Log_{10} (Energy/MeV)}$')

   def show(self,**kwargs):
      """Display the photon locations with optional ROI boundaries."""
      self.__calc_counts__()
      lon,lat,ens = self.__process_dirs__()
      s = self.pylab.scatter(lon,lat,c=N.log10(ens),s=25,alpha=1.0)
      c = self.pylab.colorbar() if not self.single_energy else None
      if self.draw_rois: self.__draw_rois__(s)
      self.__label__(s,c)
      return s
      
class EnergyAnalysis(GraphicAnalysis):
   """Analyze data in energy space, either in flux or in counts.

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  residuals   [True] display a residuals panel (not yet implemented)                                                                 
  models      [] list of spectral models to display                                                                           
  grid        [True] if True, display grid                                                              
  legend      [True] if True, display legend                                                               
  mode        ['sed'] 'sed' --> spectral energy density; 'counts' --> counts                                                                    
  cgs         [True] if True, use ergs for abscissa rather than MeV         
  e_weight    [2] energy weighting for spectral energy density plot
  normalized [True] if True, use background normalization; if false, use observed counts (not yet implemented)
  =========   =======================================================
  """

   def init(self):
      key = ['residuals','models','grid','legend','mode','cgs','e_weight','normalized']
      val = [True,[],True,True,'sed',False,2,True]
      for i in xrange(len(key)): self.__dict__[key[i]] = val[i]
      self._bounds = None

   def __setup_axes__(self):
      from pylab import axes
      self.ax = axes()
      self.ax.clear()
      self.ax.set_xscale('log')
      self.ax.set_yscale('log')

   def __drawData__(self):
      if self.mode == 'sed':
         x,y,yerr = self.pslw.spectrum(e_weight=self.e_weight,cgs=self.cgs)
         low_yerr = N.where(y-yerr <= 0,0.99*y,yerr)
         self.ax.errorbar(x=x,y=y,yerr=[low_yerr,yerr],linestyle=' ',marker='o',mfc = 'white', mec = 'black',\
                          color='black',ms=8,capsize=0,label='Observed')
         self._bounds = [0.7*min(x),1.3*max(x[y>0]),max(min(y[y>0])*0.7,max(y)/1e3),max(y)*1.3]
      else:
         x,obs = self.pslw.counts()
         yerr = N.asarray([ [y-ppf(.16,y) for y in obs], [ppf(1-.16,y)-y for y in obs] ])
         self.ax.errorbar(x=x,y=obs,yerr=yerr,linestyle=' ',marker='o',mfc = 'white', mec = 'black',\
                          color='black',ms=8,capsize=0,label='Observed')
         x,self.bg  = self.pslw.counts(background=True)
         self.ax.plot(x,self.bg,color='black',label='Background',ls=':')  
         self._bounds = [0.7*min(x),1.3*max(x),min(obs)*0.7,max(obs)*1.3]

   def __drawModel__(self,model):
      if self.mode == 'sed':
         domain = N.logspace(N.log10(self.pslw.emin),N.log10(self.pslw.emax),200)
         units = (1.60218e-6)**(self.e_weight-1) if self.cgs else 1.
         return self.ax.plot(domain,domain**self.e_weight*units*model(domain),label=model.name)
      else:
         x,y = self.pslw.counts(model=model)
         return self.ax.plot(x,y+self.bg,label=model.name,color='black')

   def __label__(self):
      ax = self.ax
      if self.legend: ax.legend(loc='lower left')
      ax.grid(b=self.grid)
      ax.set_title(self.pslw.psl.name())
      if self.mode == 'sed':
         ax.set_xlabel('$\mathrm{E\ (MeV)}$')
         en_tag = 'ergs' if self.cgs else 'MeV'
         if self.e_weight==1: ax.set_ylabel(r'$\rm{E\ dN/dE\ (ph\ cm^{-2}\ s^{-1})}$')
         elif self.e_weight==2: ax.set_ylabel(r'$\rm{E^2\ dN/dE\ (%s\ ph\ cm^{-2}\ s^{-1})}$'%en_tag)
         else: ax.set_ylabel(r'$\rm{E^%d\ dN/dE\ (%s^%d ph\ cm^{-2}\ s^{-1})}$'%(self.e_weight,en_tag,self.e_weight-1))
      else:
         ax.set_xlabel('$\mathrm{E\ (MeV)}$')
         ax.set_ylabel('$\mathrm{Counts}$')
      if self._bounds is not None: ax.axis(self._bounds)

   def show(self,**kwargs):
      """Display the spectrum in counts or flux."""
      self.__setup_axes__()
      self.__drawData__()
      for i,m in enumerate(self.models):
         self.__drawModel__(m)
         if self.mode == 'counts': self.ax.lines[-1].set_linestyle(['-','--','-.',':'][i%4])
      self.__label__()
      from pylab import show; show()

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
