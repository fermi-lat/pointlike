#import uw.pointlike
import pointlike as pl
import numpy as N
import pylab as P
import math as M
#from Response import *
import image as I
from types import *
from matplotlib.font_manager import FontProperties

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


class PValue():
   """Calculate a "p-value" map for "testing source hypothesis"."""
   def __init__(self,psl,level):
      
      self.sl=psl[level]
      self.psl=psl
      self.level=level
      from scipy.stats import poisson
      self.poisson=poisson
      self.previous_val=[-2,-2]
      self.current_val=[-1,-1]
      self.p_val=0

   def integral(self,in_dir,e1,e2):
      cv=self.current_val
      sl=self.sl
      cv[0]=sl.display(in_dir,3) #predicted mean counts
      if cv[0]==0.: #Outside of u_max
         return 0.
      cv[1]=sl.display(in_dir,1) #actual counts
      if cv!=self.previous_val:
         self.previous_val[:]=cv[:]
         self.p_val=1-poisson.cdf(cv[1],cv[0])
      return self.p_val

   def __call__(self):
      return self.psl

   def dir(self):
      return self.psl.dir()

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#


class Residuals():
   """Calculate a difference weighted by root(N)."""
   def __init__(self,psl,level):
      self.psl=psl
      self.level=level
   
   def integral(self,in_dir,e1,e2):
      sl=self.psl[self.level]
      sl.setDisplayMode(1)
      counts=sl(in_dir)
      if counts==0: return 0 #Not sure what to do about this
      sl.setDisplayMode(4)
      return sl(in_dir)/counts**0.5

   def __call__(self):
      return self.psl

   def dir(self):
      return self.psl.dir()

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def pvalue(source,level): 
   p=PValue(source(),level)
   r=I.Image(p,p,level=level,scale=2**(9-level),step=2**(9-level)/50.)
   r.show()

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

def residuals(source,level): 
   p=Residuals(source(),level)
   r=I.Image(p,p,level=level,scale=2**(9-level),step=2**(9-level)/50.)
   r.show()


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
 

def visi(source,level,nmodel):
   """Plot the integrand for the marginalization over alpha."""
   n=500
   domain=N.linspace(1e-3,1,n+1)
   photons=source.get('photons')[level]
   alpha=source.get('alphas')[level][0]
   slike=source.get('slikes')[level]
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

def spectrum(s,models=[],flag='sed',fignum=30):
   """Plot estimated spectral density, along with any model fits."""

   model = (lambda e: 1) if flag=='sed' else (lambda e: 1e-9*(100/e)**2)
   domain,domain_err,codomain,codomain_err=s.spectrum(energy_weight=1,model=model)[:4]
   #Adjust codomains for plotting, two lines below
   if flag=='sed': codomain[codomain==0]=uplim=1e-20
   else: codomain[codomain<1e-5]=uplim=1e-5
   codomain_err=N.array([N.where((codomain-codomain_err)>0,codomain_err,codomain*0.99),codomain_err])
   f9=s.f9()
   fit_flag=N.any([m.good_fit for m in models])
   
   P.figure(fignum,(6,8))
   a=P.axes([.14,.40,.8,.55]) if fit_flag else P.gca()
   a.set_xscale('log')
   a.set_yscale('log')
   a.errorbar(domain,codomain,yerr=codomain_err,xerr=domain_err,\
               linestyle=' ',marker='d',markersize=4,capsize=0) 
   if fit_flag:
      b=P.axes([0.14,0.1,0.8,0.30])
      b.set_xscale('log')
      b.set_yscale('linear')
      domain_fit=N.logspace(M.log(domain[0]/10.,10),M.log(domain[-1]*10,10),200) if flag=='sed' else domain
   
   for i,model in enumerate([x for x in models if x.good_fit]):
      if flag=='sed':
         codomain_fit=domain_fit*model(domain_fit)
         codomain_res=domain*model(domain)
         ls,ms='-',0
      else:
         codomain_fit=codomain_res=s.spectrum(energy_weight=1,model=model)[-1]
         ls,ms=' ',4
      residuals=N.nan_to_num((codomain-codomain_res)/codomain_res)
      residuals=N.where(N.abs(residuals)>10,10,residuals)
      a.plot(domain_fit,codomain_fit,label=model.name,linestyle=ls,marker='s',markersize=ms)
      bbox={'lw':1.0,'pad':10,'ec':'black','fc':'white'}
      a.text(0.025,0.025+0.15*i,model.error_string(),size='small',transform=a.transAxes,bbox=bbox,family='monospace')
      b.errorbar(domain,residuals,yerr=(codomain_err/codomain_res),\
                  xerr=domain_err,linestyle=' ',marker='d',markersize=4, capsize=0, label='Model Residuals')         
   P.axes(a)
   low_eng=domain[0]*(domain[0]/domain[1])**0.5 #Assumes logarithmic bin spacing -- fix sometime
   P.title('%s [$\mathrm{\sigma=%.2f,\, f9>%i=%.2f (1 \pm %.2f)}$]'\
               %(s().name(),(s().TS())**0.5,low_eng,f9[0],f9[1]/f9[0]))
   #P.title('%s'%s().name())
   if flag=='sed': P.ylabel(r'$\rm{E\ dN/dE\ (ph\ cm^{-2}\ s^{-1})}$')
   else: P.ylabel(r'Signal and Model Counts')
   if fit_flag:
	ticks,locs=P.xticks()
	P.xticks(ticks,tuple(['' for x in locs]))
   ma=codomain.max()*1.8
   mi=max(codomain[codomain>uplim].min()*0.3,ma/(1.8)/1e4)
   #dma,dmi=0.64*domain[N.argmin(N.abs(codomain-ma))],1.32*domain[N.argmin(N.abs(codomain-mi))]
   #dma,dmi=0.64*domain[0],1.32*domain[N.argmin(N.abs(codomain-mi))]
   dma,dmi=domain[0]*0.6,domain[-1]*1.6
   P.axis([dma,dmi,mi,ma])
   P.grid()
   P.legend()
   ax=P.axis()
   if fit_flag:
      P.axes(b)
      low,high=P.axis()[2:]
      P.plot(N.linspace(ax[0],ax[1],50),[0]*50)
      P.ylabel('$\mathrm{(obs-model)/model}$')
      P.axis([ax[0],ax[1],max(low,-.55),min(high,0.55)])
      P.grid()
   P.xlabel('$\mathrm{E\ (MeV)}$')   
   
   

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
def sed(s,models=[],fignum=30):
   spectrum(s,models=models,flag='sed',fignum=fignum)

def counts(s,models=[],fignum=30):
   spectrum(s,models=models,flag='counts',fignum=fignum)