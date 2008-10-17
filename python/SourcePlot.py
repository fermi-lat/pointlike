"""A variety of routines for plotting pointlike fits.

"""

#import uw.pointlike
import pointlike as pl
from SourceLib import *
import numpy as N
import pylab as P
import math as M
#from Response import *
#import image as I
from types import *
from matplotlib.font_manager import FontProperties

#---------------------------------------------------------------------------
class Image(object):
    def __init__(self, fun, s, emin=100, emax=200, scale = 0.5, resolution = 100, anchor = 0):
        """fun is a SkySpectrum object, s either has a dir() or ra(),dec() method"""

        import math
        grid=N.linspace(-scale, scale, resolution)
        #emin, emax = energy_range(level)
        if 'ra' in dir(s): ra, dec = s.ra(), s.dec()
        else: ra,dec = s.dir().ra(), s.dir().dec()
        cosdec = math.cos(math.radians(dec))
        self.image = N.array([fun.integral(pl.SkyDir(ra-dx/cosdec, dec-dy), emin, emax)\
                   for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
        delta = self.image.max() - self.image.min()
        self.offset=(anchor-delta/2)
        self.scale, self.ra, self.dec, self.emin, self.emax, self.resolution =\
            scale, ra, dec, emin, emax, resolution
    
    def show(self, **kwargs):

        scale=self.scale
        if 'radial' in kwargs and kwargs['radial']==True:
            d = pl.SkyDir(self.ra,self.dec)
            grid = N.linspace(-scale, scale, self.resolution)
            import math
            ra,dec = self.ra,self.dec
            cosdec = math.cos(math.radians(dec))
            grid = 180/N.pi*N.array([d.difference(pl.SkyDir(ra-dx/cosdec, dec-dy))\
                   for dy in grid for dx in grid]).reshape((len(grid),len(grid)))
            temp_xs = grid.ravel()
            sorting = N.argsort(temp_xs)
            temp_xs = temp_xs[sorting]
            xs = N.linspace(temp_xs.min(),temp_xs.max(),8)
            temp_ys = self.image.flatten()
            temp_ys = temp_ys[sorting]
            from collections import deque
            ys = [deque() for x in xrange(len(xs)-1)]
            marker=0
            for i in xrange(len(xs)-1):
               while temp_xs[marker]>=xs[i] and temp_xs[marker]<xs[i+1]:
                  ys[i].append(temp_ys[marker])
                  marker+=1
               ys[i] = N.mean(ys[i])    
               
            xs = (xs[1:]+xs[:-1])/2.
            P.subplot(122)
            P.scatter(xs,ys)
            P.subplot(121)

        if 'radial' in kwargs: kwargs.pop('radial')
        
        if 'cmap' not in kwargs:
           #Build a sweet, ad hoc color map.
           if 'anchor' in kwargs:
               anchor = kwargs['anchor']
           else: anchor=0
           vmin,vmax = self.image.min(),self.image.max()
           if vmin<=anchor and anchor<=vmax:
               breakpoint=(anchor-vmin)/(vmax-vmin)
               kwargs['cmap'] = self.properColormap(breakpoint)

        P.imshow(self.image, extent=[-scale, scale, -scale, scale],interpolation='nearest', **kwargs)
        P.axvline(0, color='white')
        P.axhline(0, color='white')
        P.colorbar()
        P.xlabel('RA Offset (deg)')
        P.ylabel('DEC Offset (deg)')
        P.title('Energy Range %d-%d'%(self.emin, self.emax))#, size=10)
    
    def properColormap(self,breakpoint=0.5):
        brp = breakpoint #make the formatting pretty
        cdict = {'red':   ((0.0, 0.0, 0.0),
                         (brp, 1.0, 1.0),
                         (1.0, 1.0, 1.0)),
               'green': ((0.0, 0.0, 0.0),
                         (brp, 1.0, 1.0),
                         (1.0, 0.0, 0.0)),
               'blue':  ((0.0, 1.0, 1.0),
                         (brp, 1.0, 1.0),
                         (1.0, 0.0, 0.0))
              }
        return P.matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#


class NormedPSF():
   """Used for plotting the PSF matched to the estimated signal fraction."""
   def __init__(self,psl,level):
      self.solid_angle=pl.Healpix(2**level).pixelArea()
      self.psl=psl

   def integral(self,in_dir,e1,e2):
      return self.psl.integral(in_dir,e1,e2)*self.solid_angle

   def __call__(self):
      return self.psl

   def dir(self):
      return self.psl.dir()

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PySkyFun(object):

   def __init__(self,psl,sls):
      self.psl = psl
      self.sls = sls if (type(sls) is ListType or type(sls) is N.ndarray) else [sls]
      self.__myinit__()
   
   def __call__(self):
      return self.psl

   def dir(self):
      return self.psl.dir()

   def __myinit__(self):
      pass
   
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#


class PValue(PySkyFun):
   """Calculate a "p-value" map for "testing source hypothesis"."""
 
   def __myinit__(self):
      
      from scipy.stats import poisson
      self.poisson_cdf = poisson.cdf
      self.previous_val=[-2,-2]
      self.p_val=0

   def integral(self,in_dir,e1,e2):
      cv = [0,0] #Current value

      for sl in self.sls:
         cv[0] += sl.display(in_dir,3) #predicted mean counts
         cv[1] += sl.display(in_dir,1) #actual counts
      if cv[0] == 0.: return 0. #Outside of u_max
      if cv!=self.previous_val:
         self.previous_val[:]=cv[:]
         self.p_val=1-poisson.cdf(cv[1],cv[0])
      #if self.p_val<0.01: print cv
      return self.p_val

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#


class Residuals(PySkyFun):
   """Calculate a difference weighted by root(N)."""

   def integral(self,in_dir,e1=100,e2=200):
      counts,obs_counts=[0,0]
      for sl in self.sls:
         counts += sl.display(in_dir,3) #predicted mean counts
         obs_counts += sl.display(in_dir,1) #actual counts
      if counts == 0: return 0.
      #std_diff = (obs_counts-counts)/counts**0.5
      #return std_diff if abs(std_diff)<5 else 5
      return 1 - counts/obs_counts


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def resid_root(source,mode='PValue',emin=1000,emax=10000,event_class=0):
   
   bands = source.response.bands(infinity=True)
   args = []
   for i in xrange(len(bands)-1):
      if (  (bands[i]>=emin and bands[i+1]<emax) or 
            (bands[i]<emin and bands[i+1]>emin)  or
            (bands[i]<emax and bands[i+1]>emax)   ) :  args += [i]
   args=N.array(args)
   #Test whether both events are present
   n = len(source.photons)
   if n > len(bands)-1:
      print 'True'
      if event_class == 0:
         index = args
      elif event_class == 1:
         index = n/2+args
      else:
         index = N.concatenate(args,args+n/2)
   else: index = args
   slikes = source.slikes[index]
   if args.max()<len(bands): args = N.append(args, args[-1]+1)
   energies = bands[args]
   print 'Energy range: %d-%d'%(energies.min(),energies.max())
   umax = pl.SimpleLikelihood.defaultUmax()
   scale = 0.5*(2*umax)**0.5*pl.IParams.sigma(energies.min(),event_class)*180/N.pi
   print 'Scale: %.2f'%scale

   exec('p=%s(source.psl,slikes)'%mode)
   r=Image(p, p, scale=scale, resolution=100, emin=energies.min(), emax=energies.max())
   if mode!='Residuals':
      r.show(cmap=P.get_cmap('PuOr'))
   else: r.show()


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

#Convenience functions
def pvalue(source,emin=1000,emax=10000,event_class=0):
   resid_root(source,mode='PValue',emin=emin,emax=emax,event_class=event_class)

def residuals(source,emin=1000,emax=10000,event_class=0):
   resid_root(source,mode='Residuals',emin=emin,emax=emax,event_class=event_class)


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def display_map(source,level,mode=1):
   """Just a wrapper for making an image with the various display modes in SimpleLikelihood.
      
         mode: 0: display psf density
               1: counts map by HealPixel
               2: background counts map by HealPixel
               3: predicted counts map by HealPixel
               4: difference counts map... you guessed it, by HealPixel"""   
   p=source()
   s=pl.PSLdisplay(p,mode)
   r=I.Image(s,p,level=level,scale=2**(9-level),step=2**(9-level)/50.)
   r.show()



#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
 

def log_like(source,level,exp=False):
   """Plot the log likelihood for alpha for a level of the source.  If 'exp' is
      true, plot the likelihood (times an arbitrary constant)."""
   domain=N.linspace(0,5,100)
   if not exp:
      codomain=[-source()[level](x) for x in domain]
   else:
      import math as M
      codomain=[M.exp(-source()[level](x)+source()[level]()) for x in domain]
   P.plot(domain,codomain)
   P.xlabel(r'$\alpha$')
   if not exp: P.ylabel('Log Likelihood')
   else: P.ylabel('Likelihood')
   P.title('Level %i'%(level))

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
 

def visi(source,level,nmodel,log = False):
   """Plot the integrand for the marginalization over alpha."""
   n=500
   domain=N.linspace(1e-3,1,n+1)
   photons=source.photons[level]
   alpha=source.alphas[level,0]
   slike=source.slikes[level]
   codomain1=N.asarray([M.exp(-slike(x)+slike()) for x in domain])
   norm=(N.append( ([1.]+[4.,2.]*(n/2))[:-1],1.)*codomain1).sum()/(3*n)
   codomain1/=norm
   P.subplot(221)
   P.plot(domain,codomain1)
   from scipy.stats import poisson
   codomain2=poisson.pmf(photons,nmodel/domain)
   P.subplot(222)
   P.plot(domain,codomain2)
   P.subplot(223)
   P.plot(domain,codomain1*codomain2)
   P.axvline([alpha])

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def rad_resid(source,energy=1000,event_class=0,resids_only=False):
   """The coolest resid in town."""

   bands = source.response.bands(infinity=True)
   for i in xrange(len(bands)-1):
      if bands[i] <= energy and bands[i+1] > energy: break
   emin,emax = bands[i],bands[i+1]
   sl = source.slikes[i+len(source.slikes)/2*event_class]
   sigma,gamma,alpha = sl.band().sigma(),sl.band().gamma(),sl.alpha()
   max_ROI = sigma*(2*sl.umax())**0.5 #max_ROI in radians
   wsdl = pl.WeightedSkyDirList(sl.band(), source().dir(), max_ROI)
   diff = source().dir().difference
   dists = N.concatenate([ [diff(x)]*int(x.weight()) for x in wsdl ])
   y,x = N.histogram(dists,bins=N.linspace(0,max_ROI,10),new=True)#int(max(10,len(dists)/14.)))
   print x
   us = (0.5/sigma**2)*( x**2 )
   integ_psf = (1+us[:-1]/gamma)**(1-gamma) - (1+us[1:]/gamma)**(1-gamma)
   fmax = 1- (1+us[-1]/gamma)**(1-gamma)
   prediction = alpha/fmax*integ_psf + (1-alpha)/us[-1]*(us[1:]-us[:-1])
   centers = ((us[1:]+us[:-1])/2)**0.5
   if not resids_only:
      a=P.axes([.14,.40,.8,.55])
      P.cla()
      a.errorbar(centers,y,yerr = (y**0.5),linestyle=' ')
      print x
      print y
      #a.bar(us[:-1]**0.5,y,yerr = y**0.5,width=centers[1]-centers[0],facecolor='white')
      a.plot(centers,prediction*len(dists), linestyle=' ',marker='x')
      locs=a.get_xticks()
      a.set_xticklabels(['']*len(locs))
      a.set_ylabel('Counts')

   b=P.axes([0.14,0.1,0.8,0.30]) if not resids_only else P.gca()
   #P.cla()
   resids = -(prediction*len(dists)-y)/(prediction*len(dists))
   b.errorbar(centers,resids,yerr=y**0.5/(prediction*len(dists)),linestyle=' ')
   b.plot(N.linspace(0,us[-1]**0.5,100),[0]*100,color='black')
   b.set_xlabel('root(u)')
   b.set_ylabel('Relative Residuals')
   lims = b.axis() if resids_only else a.axis()
   b.axis([lims[0],lims[1],-0.5,+0.5])
   the_ax = b if resids_only else a
   P.axes(the_ax)
   P.title('%d-%d MeV'%(emin,emax))
   

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def spectrum(s,models=None,flag='sed',fignum=30,sedweight=1,cgs=False):
   """Plot estimated spectral density, along with any model fits."""

   colors = ['green','red','orange','purple'] #Go Tigers!
   if models is None: models = []
   from types import ListType
   if type(models) is not ListType: models = [models]
   model = (lambda e: 1) if flag=='sed' else (lambda e: 1e-9*(100/e)**2)
   domain,domain_err,codomain,up_codomain_err,down_codomain_err =\
      s.spectrum(energy_weight=sedweight,model=model,cgs=cgs)[:5]
   low_eng = N.round(domain-domain_err[0])[codomain>0][0]
   #Adjust codomains for plotting, two lines below
   if flag=='sed': codomain[codomain==0]=uplim=1e-20
   else: codomain[codomain<1e-5]=uplim=1e-5
   down_codomain_err=N.where((codomain-down_codomain_err)>0,down_codomain_err,codomain*0.99)
   f7=s.i_flux()
   fit_flag=N.any([m.good_fit for m in models])
   
   P.figure(fignum,(6,8))
   a=P.axes([.14,.40,.8,.55]) if fit_flag else P.gca()
   a.set_xscale('log')
   a.set_yscale('log')
   #codomain_err = N.asarray([[up_codomain_err[i],down_codomain_err[i]] for i in xrange(len(up_codomain_err))])
   codomain_err = N.array([down_codomain_err,up_codomain_err])
   a.errorbar(domain/1000.,codomain,yerr=codomain_err,xerr=domain_err,\
               linestyle=' ',marker='d',markersize=4,capsize=0,color='blue') 
   if fit_flag:
      b=P.axes([0.14,0.1,0.8,0.30])
      b.set_xscale('log')
      b.set_yscale('linear')
      domain_fit=N.logspace(M.log(domain[0]/10.,10),M.log(domain[-1]*10,10),200) if flag=='sed' else domain
   
   for i,model in enumerate([x for x in models if x.good_fit]):
      if flag=='sed':
         cgs_scale = sedweight*1.60218e-6 if cgs else 1
         codomain_fit=domain_fit**sedweight*model(domain_fit)*cgs_scale
         codomain_res=domain**sedweight*model(domain)*cgs_scale
         ls,ms='-',0
      else:
         codomain_fit=codomain_res=s.spectrum(energy_weight=1,model=model)[-1]
         ls,ms=' ',4
      residuals=N.nan_to_num((codomain-codomain_res)/codomain_res)
      residuals=N.where(N.abs(residuals)>10,10,residuals)
      a.plot(domain_fit/1000.,codomain_fit,label=model.pretty_name,linestyle=ls,marker='s',markersize=ms,color=colors[i%len(colors)],lw=2)
      bbox={'lw':1.0,'pad':10,'ec':'black','fc':'white'}
      a.text(0.025,0.025+0.15*i,str(model),size='small',transform=a.transAxes,bbox=bbox,family='monospace')
      b.errorbar(domain/1000.,residuals,yerr=(codomain_err/codomain_res),\
                  xerr=domain_err,linestyle=' ',marker='d',markersize=4, \
                  capsize=0, label='Model Residuals', color=colors[i%len(colors)])         
   P.axes(a)
   
   #low_eng=round(domain[0]*(domain[0]/domain[1])**0.5) #Assumes logarithmic bin spacing -- fix sometime
   P.title(r'%s [$\mathrm{\sigma=%.2f,\, f7>%i=%.2f (1 \pm %.2f)}$]'\
               %(s().name(),(s().TS())**0.5,low_eng,f7[0],f7[1]/f7[0]))
   #P.title('%s'%s().name())
   if flag=='sed':
      en_tag = 'ergs' if cgs else 'MeV'
      if sedweight==1: P.ylabel(r'$\rm{E\ dN/dE\ (ph\ cm^{-2}\ s^{-1})}$')
      elif sedweight==2: P.ylabel(r'$\rm{E^2\ dN/dE\ (%s\ ph\ cm^{-2}\ s^{-1})}$'%en_tag)
      else: P.ylabel(r'$\rm{E^%d\ dN/dE\ (%s^%d ph\ cm^{-2}\ s^{-1})}$'%(sedweight,en_tag,sedweight-1))
   else: P.ylabel(r'Signal and Model Counts')
   if fit_flag:
	ticks,locs=P.xticks()
	P.xticks(ticks,tuple(['' for x in locs]))
   ma=codomain.max()*1.8
   mi=max(codomain[codomain>uplim].min()*0.3,ma/(1.8)/1e4)
   #dma,dmi=0.64*domain[N.argmin(N.abs(codomain-ma))],1.32*domain[N.argmin(N.abs(codomain-mi))]
   #dma,dmi=0.64*domain[0],1.32*domain[N.argmin(N.abs(codomain-mi))]
   dma,dmi=domain[0]*0.6/1000,domain[-1]*1.6/1000
   P.axis([dma,dmi,mi,ma])
   P.grid()
   P.legend()
   ax=P.axis()
   if fit_flag:
      P.axes(b)
      low,high=P.axis()[2:]
      P.plot(N.linspace(ax[0],ax[1],50),[0]*50,color='black')
      P.ylabel('$\mathrm{(obs-model)/model}$')
      P.axis([ax[0],ax[1],max(low,-.55),min(high,0.55)])
      P.grid()
   P.xlabel('$\mathrm{E\ (GeV)}$')   
   
   

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def logn_logts(sl,cut_ts=10,bins=30):
   cut_ts=cut_ts**0.5
   ts=N.array([x().TS()**0.5 for x in sl])
   ts=ts[ts>cut_ts] #Cut on ts
   bin_width=(ts.max()/cut_ts)**(1./bins)
   bins=N.array([cut_ts*bin_width**x for x in xrange(bins+1)]) #Make logarithmic bins
   widths=bins[1:]-bins[:-1]
   counts=N.array([N.where(ts>=x,1,0).sum() for x in bins])
   #P.loglog(bins,counts,linestyle='-',marker='s')
   a=P.gca()
   a.set_xscale('log')
   a.set_yscale('log')
   P.bar(bins[:-1],counts[:-1],width=widths,facecolor='w')
   P.ylabel('N>=root(TS)')
   P.xlabel('root(TS)')
   P.title('LogN-LogTS (TS>%i)'%(int(cut_ts**2)))

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def logn_logs(sl,cut_ts=10,nbins=30):
   flux=N.fromiter((x.flux() for x in sl),float)
   ts=N.fromiter((x().TS() for x in sl),float)
   mask=ts>=cut_ts
   flux=flux[mask]
   bin_width=(flux.max()/flux.min())**(1./nbins)
   bins=N.array([flux.min()*bin_width**x for x in xrange(nbins+1)]) #Make logarithmic bins
   widths=bins[1:]-bins[:-1]
   counts=N.array([N.where(flux>=x,1,0).sum() for x in bins])
   #P.loglog(bins,counts,linestyle='-',marker='s')
   a=P.gca()
   a.set_xscale('log')
   a.set_yscale('log')
   P.bar(bins[:-1],counts[:-1],width=widths,facecolor='w')
   P.ylabel('N>=S')
   P.xlabel('S > 100MeV (ph cm^-2 s^-1)')
   P.title('LogN-LogS (TS>%i)'%(int(cut_ts)))


#Convenience functions
def sed(s,models=[],fignum=30,cgs=False):
   spectrum(s,models=models,flag='sed',fignum=fignum,cgs=cgs,sedweight=2)

def counts(s,models=[],fignum=30,cgs=False):
   spectrum(s,models=models,flag='counts',fignum=fignum,cgs=cgs)