"""Docstring for CVS here?."""

import numpy as N

class SimpleLikelihoodWrapper(object):
   """Wrap up SimpleLikelihood as the basic object for spectral fitting."""

   def init(self):
      """Establish defaults as class members."""
      self.n_simps = 8 #sampling points in log space for Simpson's rule
      self.emulate_unbinned = True

   def __init__(self,sl,exposure,sky_dir,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.sl,self.sky_dir,self.exposure = sl,sky_dir,exposure
      emin,emax,ec = sl.band().emin(),sl.band().emax(),sl.band().event_class()
      self.emin,self.emax=emin,emax
      
      #Pieces for calculating the expected counts under a given model
      self.sampling_points  = N.logspace(N.log10(emin),N.log10(emax),self.n_simps+1)
      self.exposure_points  = N.fromiter( (exposure.value(sky_dir,e,ec) for e in self.sampling_points) , float)
      self.simpsons_weights = N.log(self.sampling_points[-1] / self.sampling_points[0]) * \
                              N.asarray([1.] + ([4.,2.]*(self.n_simps/2))[:-1] + [1.])/(3.*self.n_simps)
      self.psf_correction   = 1 - (1+sl.umax()/sl.gamma())**(1-sl.gamma()) #fraction of signal in ROI

      #Pieces for calculating marginalization integral; needed for speed!
      if not self.sl.extended_likelihood():
         n_points = max((int(100*(self.sl.photons()/10.)**0.5) >> 1) << 1,100) #force n_points to be even
         print n_points
         self.points = N.linspace(0.001,1,n_points+1)
         self.weights = N.asarray([1.] + ([4.,2.]*(n_points/2))[:-1] + [1.])/(3.*n_points)
         self.point_likes = N.exp( N.nan_to_num(N.fromiter( (-self.sl(a)+self.sl() for a in self.points) , float)) )

   def expected(self,model):
      """Return expected counts in the ROI for a source model.
      
         Integration is by Simpson's rule in log space."""

      return self.psf_correction*\
            (self.simpsons_weights*self.sampling_points*model(self.sampling_points)*self.exposure_points).sum()

   def avg_exposure(self):
      """Return the logarithmic average of the exposure."""
      return self.expected(model=lambda e: 1./e)/N.log(self.emax/self.emin)/self.psf_correction

   def logLikelihood(self,model):
      """Return the (negative) log likelihood.  If extended likelihood is enabled, use that.
         Otherwise, return the alpha marginalization."""

      expected = self.expected(model)
      if self.sl.photons() == 0 and not self.emulate_unbinned: return expected

      if self.sl.extended_likelihood(): return self.sl.logLikelihood(expected)

      from scipy.stats import poisson
      poiss_likes = N.nan_to_num(poisson.pmf(self.sl.photons(),expected/self.points))
      integral = (self.weights*self.point_likes*poiss_likes).sum()
      if integral == 0.: return 1e6
      return -N.log(integral)

   def flux(self,e_weight=0,cgs=True):
      """Return the differential flux, multiplied by energy**e_weight.

         The estimator here is chosen to be optimal with constant exposure, i.e., if one expands
         an arbitrary spectral model about some energy e0 and assumes constant exposure, then the
         first-order integral vanishes if e0 is the arithmetic mean.  Thus, the flux is to be
         evaluated at the arithmetic mean, properly."""
      units = (1.60218e-6)**(e_weight-1) if cgs else 1 #convert MeV to ergs
      multi = units*self.energy()**e_weight/self.avg_exposure()/(self.emax-self.emin)/self.psf_correction
      a,sa,signal,photons = self.sl.alpha(),self.sl.sigma_alpha(),self.sl.signal(),self.sl.photons()
      error = (a**2*photons + (sa*photons)**2)**0.5
      return (multi*signal,multi*error)

   def energy(self):
      """Return energy 'associated' with this band, taken to mean where the flux should optimally be evaluated"""
      #return (self.sl.band().emax() + self.sl.band().emin())/2.
      return (self.emin*self.emax)**0.5



class PointSourceLikelihoodWrapper(object):

   def init(self):

      self.emin = 100
      self.emax = 1e6
      self.quiet = False

   def __init__(self,psl,exposure,**kwargs):

      self.init()
      self.__dict__.update(kwargs)
      self.psl = psl
      self.sl_wrappers = [SimpleLikelihoodWrapper(sl,exposure,psl.dir()) for sl in psl \
                          if sl.band().emin() >= self.emin and sl.band().emax() < self.emax]
      self.bin_centers = N.sort(list(set([x.energy() for x in self.sl_wrappers])))
      ens = N.asarray([ [x.emin,x.emax] for x in self.sl_wrappers]).flatten()
      self.emin,self.emax = ens.min(),ens.max() #replace defaults with observed
      

   def least_squares(self,model,quiet=False):
      from scipy.optimize import leastsq

      #Products for fit
      photons = N.fromiter( (sl.sl.photons() for sl in self.sl_wrappers ), int )
      alphas  = N.array([ [sl.sl.alpha(),sl.sl.sigma_alpha()] for sl in self.sl_wrappers ]).transpose()
      errors  = (alphas[0,:]**2*photons+alphas[1,:]**2*photons.astype(float)**2)**0.5
      signals = N.column_stack([alphas[0,:]*photons,errors]).transpose()

      def chi(parameters,*args):
         """Routine for use in least squares routine."""
         model,photons,signals,sl_wrappers = args
         model.p = parameters
         expected = N.fromiter( (sl.expected(model) for sl in sl_wrappers ), float)
         return ((expected - signals[0,:])/signals[1,:])[photons > 0]

      try:
         fit = leastsq(chi,model.p,args=(model,photons,signals,self.sl_wrappers),full_output=1)
      except:
         fit = [0]*5

      if fit[4] == 1:

         model.good_fit=True #Save parameters to model
         model.cov_matrix=fit[1] #Save covariance matrix
         model.p=fit[0]
         vals = chi(fit[0],model,photons,signals,self.sl_wrappers)
         model.dof=len(vals)
         model.chi_sq=(vals**2).sum() #Save chi-sq information

         if not self.quiet and not quiet:
            print '\nFit converged!  Chi-squared at best fit: %.4f'%model.chi_sq
            print str(model)+'\n'
      return model

   
   def poisson(self,model,prefit=True):
      from scipy.optimize import fmin
      if prefit: self.least_squares(model,quiet=True)

      sl_wrappers = self.sl_wrappers

      def logLikelihood(parameters,*args):
         """Routine for use in minimum Poisson likelihood."""
         model = args[0]
         model.p = parameters
         return sum( (sl.logLikelihood(model) for sl in sl_wrappers) )
      
      fit=fmin(logLikelihood,model.p,args=(model,),full_output=1,disp=0,maxiter=1000,maxfun=2000)
      warnflag=(fit[4]==1 or fit[4]==2)
      
      if not warnflag: #Good fit (claimed, anyway!)

         def hessian(m,mf,*args):
            """Calculate the Hessian; f is the minimizing function, m is the model,args additional arguments for mf."""
            delt=0.01
            p = m.p.copy()
            hessian=N.zeros([len(p),len(p)])
            for i in xrange(len(p)):
               for j in xrange(i,len(p)): #Second partials by finite difference
                  
                  xhyh,xhyl,xlyh,xlyl=p.copy(),p.copy(),p.copy(),p.copy()

                  xhyh[i]*=(1+delt)
                  xhyh[j]*=(1+delt)

                  xhyl[i]*=(1+delt)
                  xhyl[j]*=(1-delt)

                  xlyh[i]*=(1-delt)
                  xlyh[j]*=(1+delt)

                  xlyl[i]*=(1-delt)
                  xlyl[j]*=(1-delt)

                  hessian[i][j]=hessian[j][i]=(mf(xhyh,m,*args)-mf(xhyl,m,*args)-mf(xlyh,m,*args)+mf(xlyl,m,*args))/\
                                                (p[i]*p[j]*4*delt**2)

            m.p = p #Restore parameters
            return hessian
            
         try:
            from numpy.linalg import inv
            model.cov_matrix = inv(hessian(model,logLikelihood))
            model.good_fit   = True
            model.p          = fit[0]
            model.logl       = -fit[1] #Store the log likelihood at best fit
            if not self.quiet:
               print '\nFit converged!  Function value at minimum: %.4f'%fit[1]
               print str(model)+'\n'
         except:
            print 'Hessian inversion failed!'
         return model

   def spectrum(self,e_weight=0,cgs=True):
      keys,d = [str(e) for e in self.bin_centers.astype(int)],dict()
      for key in keys: d[key] = N.zeros(3).astype(float) #flux, error, num
      for x in self.sl_wrappers:
         t = d[str(int(x.energy()))]
         f,e = x.flux(e_weight=e_weight,cgs=cgs)
         t += N.asarray([f,e**2,1.])
      results = N.asarray([ [d[key][0]/d[key][2],d[key][1]**0.5/d[key][2]] for key in keys])
      return (self.bin_centers,results[:,0],results[:,1])
      #TODO -- guard against 0

   def counts(self,model,background=True):
      keys,d = [str(e) for e in self.bin_centers.astype(int)],dict()
      for key in keys: d[key] = N.zeros(4).astype(float)
      for x in self.sl_wrappers:
         t = d[str(int(x.energy()))]
         model_source,model_background,counts,signal_frac,signal_frac_err = \
            x.expected(model),x.sl.background(),x.sl.photons(),x.sl.alpha(),x.sl.sigma_alpha()
         if background:
            t += N.asarray([model_source+model_background,counts,(model_source+model_background),1.])
         else:
            t += N.asarray([model_source,signal_frac*counts,signal_frac_err*signal_frac*counts,1.]) #Not right error
      results = N.asarray([ [d[key][0]/d[key][3],d[key][1]/d[key][3],d[key][2]**0.5/d[key][3] ] for key in keys])
      return (self.bin_centers,results[:,0],results[:,1],results[:,2])

class PointspecPlotter(object):

   def init(self):
      self.residuals = True
      self.ax = None
      self.models = []
      self.grid = True
      self.legend = True
      self.mode = 'sed' #alternative is 'counts'
      self.cgs = False
      self.e_weight = 2
      self.background_counts = False #use background model to calculate counts
   
   def __init__(self,pslw,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.pslw = pslw
      self.models = [x for x in self.models] #copy
      self.replot()

   def replot(self):
      exec('self.%s()'%self.mode)

   def addModel(self,model):
      if model not in self.models:self.models += [model]
      ax = self.ax
      if N.any([line.model == model for line in ax.lines]): return
      t_axis = ax.axis()
      if self.mode == 'sed':
         domain = N.logspace(N.log10(self.pslw.emin),N.log10(self.pslw.emax),200)
         units = (1.60218e-6)**(self.e_weight-1) if self.cgs else 1.
         lines = ax.plot(domain,domain**self.e_weight*units*model(domain),label=model.name)
         lines[0].model = model

      else:
         x,y,obs,yerr = self.pslw.counts(model,background=self.background_counts)
         low_yerr = N.where(y-yerr <=0,0.99*y,yerr)
         lines = ax.errorbar(x=x,y=y,yerr=[low_yerr,yerr],linestyle=' ',label=model.name)
         lines[0].model,lines[1][0].model,lines[1][1].model = model,model,model
      if self.legend: ax.legend(loc='lower left')
      ax.axis(t_axis)


   def delModel(self,model):
      for i in xrange(len(self.models)):
         if self.models[i] == model:
            self.models.pop(i)
            break
      for i in xrange(len(self.sed_ax.lines)):
         if self.sed_ax.lines[i].model == model:
            self.sed_ax.lines.pop(i)
            break

   def common(self):
      from pylab import axes
      self.ax = axes()
      self.ax.clear()
      self.ax.set_xscale('log')
      self.ax.set_yscale('log')

   def sed(self):
      self.common()
      ax = self.ax
      e_weight,cgs = self.e_weight,self.cgs
      x,y,yerr = self.pslw.spectrum(e_weight=self.e_weight,cgs=self.cgs)
      low_yerr = N.where(y-yerr <= 0,0.99*y,yerr)

      lines = ax.errorbar(x=x,y=y,yerr=[low_yerr,yerr],linestyle=' ',marker='o',capsize=0,label='Observed')
      lines[0].model = None,None,None
      for m in self.models:
         self.addModel(m)
      ax.grid(b=self.grid)
      if self.legend: ax.legend(loc='lower left')
      ax.axis([0.7*min(x),1.3*max(x),min(y[y>0])*0.7,max(y)*1.3])
      ax.set_xlabel('$\mathrm{E\ (MeV)}$')
      en_tag = 'ergs' if cgs else 'MeV'
      if e_weight==1: ax.set_ylabel(r'$\rm{E\ dN/dE\ (ph\ cm^{-2}\ s^{-1})}$')
      elif e_weight==2: ax.set_ylabel(r'$\rm{E^2\ dN/dE\ (%s\ ph\ cm^{-2}\ s^{-1})}$'%en_tag)
      else: ax.set_ylabel(r'$\rm{E^%d\ dN/dE\ (%s^%d ph\ cm^{-2}\ s^{-1})}$'%(e_weight,en_tag,e_weight-1))
      ax.set_title(self.pslw.psl.name())
      ax.e_weight,ax.cgs,ax.counts = e_weight,cgs,False
      return ax

   def counts(self):
      self.common()
      ax = self.ax
      x,y,obs,yerr = self.pslw.counts(lambda e: 1,background=self.background_counts)
      lines = ax.plot(x,obs,linestyle=' ',marker='o',label='Observed')
      lines[0].model = None
      markers=['s','d','h','x','^','+','>']
      for i,m in enumerate(self.models):
         self.addModel(m)
         ax.lines[-1].set_marker(markers[i])

      ax.grid(b=self.grid)
      if self.legend: ax.legend(loc='lower left')

      ax.axis([0.7*min(x),1.3*max(x),min(obs)*0.7,max(obs)*1.3])
      ax.set_xlabel('$\mathrm{E\ (MeV)}$')
      ax.set_ylabel('Counts')
      ax.set_title(self.pslw.psl.name())
      return ax

