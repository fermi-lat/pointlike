"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/toagen.py,v 1.16 2013/07/16 03:00:22 kerrm Exp $

Calculate TOAs with a variety of methods.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""

from __future__ import division
import sys
import numpy as np
import pylab as pl
import scipy.stats
from scipy.optimize import fmin,golden,brentq
from stats import hm,hmw,sf_hm
from edf import EDF,find_alignment
from collections import deque

try:
    import fftfit
    import psr_utils
except:
    pass

# other module to deal with weights -- find when the significance stops increasing and use
# that set of photons for fits (via min weight);

# also, might want to use a posterior probability?
# i.e., factoring in the timing template; iterate?  will it converge?
# do not want to do this -- taken care of in likelihood itself

SECSPERDAY = 86400.

class TOAGenerator(object):
    """Manage a data set and set of options to produce the required TOAs from LAT events."""

    def init(self):
        pass

    def __init__(self,data,polyco,template,**kwargs):

        self.data = data; self.polyco = polyco; self.template = template
        self.init()
        self.__dict__.update(**kwargs)
        self.mean_err = 1
        # keep track of errors in phase shift for later analysis
        self.phases = deque()
        self.phase_errs = deque()

    def get_toas(self,binner,use_midpoint=True):
        """ Calculate the TOAs specified by the binner."""

        self.phases.clear(); self.phase_errs.clear()

        # Compute observation duration for each TOA
        toas = np.empty(binner.ntoa)
        err_toas = np.empty(binner.ntoa)
        tim_strings = ['FORMAT 1']
        self.counter = 0

        for ii,(mjdstart,mjdstop) in enumerate(binner):

            # Compute freq and period at middle of observation
            tmid = (mjdstop + mjdstart)/2.
            pe   = self.polyco.getentry(tmid)
            freq = pe.evalfreq(tmid)
            period = 1.0/freq

            # Compute phase at start of observation or at midpoint
            phase_time = tmid if use_midpoint else mjdstart
            pe = self.polyco.getentry(phase_time)
            polyco_phase0 = pe.evalphase(phase_time)
            
            # Select phases
            phases,weights = self.data.toa_data(mjdstart,mjdstop)
            if len(phases) == 0: continue
            
            tau,tau_err,prob,logl = \
                self.get_phase_shift(phases,weights,polyco_phase0)
            self.mean_err = (self.mean_err*ii + tau_err)/(ii+1)

            # compute RMS between linear assumption and actual polyco
            # over the integration period
            dom = np.linspace(mjdstart,mjdstop,100)
            ph = self.polyco.vec_evalabsphase(dom)
            ph -= ph.min()
            p = np.polyfit(dom-tmid,ph,1)
            rms = (ph-np.polyval(p,dom-tmid)).std()*(1e6/freq) # mus
            
            # Prepare a string to write to a .tim file or to send to STDOUT
            toa = phase_time + (tau*period)/SECSPERDAY
            toa_err = tau_err*period*1.0e6
            frac_err = tau_err
            frame_label = 'BAT' if self.data.bary else 'GEO'
            weight_string = '' if (weights is None) else '-nwp %.2f'%(weights.sum())
            duration_string = '-tstart %s -tstop %s'%(mjdstart,mjdstop)
            rms_string = '-drms %.2f'%(rms)
            logl_string = '-logl %.2f'%(logl)
            s = " %s 0.0 %.12f %.2f %s -i LAT %s -np %d %s -chanceprob %.2e -pherr %.3f %s %s" % (frame_label,toa,toa_err,pe.obs,duration_string,len(phases),weight_string,prob,frac_err,rms_string,logl_string)
            toas[ii] = toa
            err_toas[ii] = toa_err
            tim_strings.append(s)
            self.counter += 1

        # Note TOAS in MJD, err_toas in microseconds, tim_strings a line for a FORMAT 1 .tim file
        return toas,err_toas,tim_strings

class UnbinnedTOAGenerator(TOAGenerator):

    def init(self):
        #self.seeds = np.arange(0.01,0.991,0.01)
        self.seeds = np.arange(0.00,0.991,0.01)
        self.good_ephemeris = True
        self.phi0 = self.template.get_location()
        self.prev_peak = self.phi0
        self.phase_drift = False
        self.weights = None
        self.plot_stem = None
        self.display = True

    def __toa_error__(self,val,*args):
        f      = self.__toa_loglikelihood__
        delta  = 0.01
        d2 = (f( [val + delta], *args) - 2*f([val], *args) + f( [val - delta], *args))/(delta)**2
        if d2 < 0:
            print 'WARNING! Could not estimate error properly.  Setting to a large value...'
            return 0.2
        return d2**-0.5

    def __toa_loglikelihood__(self,p,*args):
        self.template.set_overall_phase(p[0])
        if args[1] is None:
            return -np.log(self.template(args[0])).sum()
        #return -np.log(1+args[1]*(self.template(args[0],suppress_bg=True)-1)).sum()
        return -np.log(1+args[1]*(self.template(args[0])-1)).sum()

    def get_phase_shift_old(self,phases,weights,polyco_phase0):

        f   = self.__toa_loglikelihood__
        jump = False
        # if the ephemeris is deemed good, we use it to track the 
        # solution such that we don't adopt a fluctuation for the TOA 
        if (self.good_ephemeris):
            seed = [self.prev_peak]
            fit  = fmin(f,seed,args=(phases,weights),disp=0,ftol=1e-9,full_output=True)
            best_ll = fit[1]
            jump = abs(self.prev_peak - fit[0][0])/self.mean_err
            if jump > 10:
                print 'Found a jump, doing a blind search now.'
                print self.prev_peak,fit[0][0],self.mean_err
            jump = jump > 10
            peak_shift = (fit[0][0] - self.phi0)
            if abs(peak_shift) > 0.25:
                self.phase_drift = True
            elif (not jump):
                self.prev_peak = fit[0][0]
                tau     = (peak_shift - polyco_phase0)
                tau_err = self.__toa_error__(fit[0][0],phases,weights)
                if self.display: print 'Peak Shift: %.5f +/- %.5f'%(peak_shift,tau_err)
                self.phases.append(peak_shift)
            
        if (not self.good_ephemeris) or self.phase_drift or jump:
            # look for the absolute best fit without any prior information

            best_ll,tau = np.inf,0
            seed_vals = [f([x],phases,weights) for x in self.seeds]
            top10     = self.seeds[np.argsort(seed_vals)][:2] #NB change
            for seed in top10:
                seedfit   = fmin(f,[seed],args=(phases,weights),disp=0,ftol=1e-9,full_output=True)
                if seedfit[1] < best_ll:
                    fit = seedfit
                    best_ll = fit[1]
                    tau = fit[0][0]
                    if jump: self.prev_peak = tau

            tau_err = self.__toa_error__(tau,phases,weights)         
            tau -= (self.phi0 + polyco_phase0)
            if self.display: 
                print '(Blind) Peak Shift: %.5f +/- %.5f'%(
                    tau+polyco_phase0,tau_err)
            self.phases.append(tau+polyco_phase0)        

        if (self.plot_stem is not None):
            dom1 = np.linspace(0,1,100)
            cod1 = np.asarray([f([x],phases,weights) for x in dom1])
            dom2 = np.linspace(fit[0][0]-0.04,fit[0][0]+0.04,30)
            cod2 = np.asarray([f([x],phases,weights) for x in dom2])
            pl.figure(10); pl.clf();
            ax1 = pl.gca()
            ax1.axhline(0,color='k')
            # calculate coordinates for inset; should use transAxes?
            ax2 = pl.axes([0.2+0.4*(fit[0][0]<0.5),0.60,0.25,0.25])
            for i,(dom,cod,ax) in enumerate(zip([dom1,dom2],[cod1,cod2],[ax1,ax2])):
                ax.plot(dom,cod)
                ax.axvline(fit[0][0],color='red')
                ax.axvline(self.phi0,color='k',ls='-')
                if i==1:
                    ax.xaxis.set_major_locator(pl.matplotlib.ticker.MaxNLocator(4))
                    ax.axis([dom[0],dom[-1],cod.min(),cod.min()+5])
                    ax.axhline(cod.min()+0.5,color='blue',ls='--')
                    ax.axhline(cod.min()+2.0,color='blue',ls='-.')
                    ax.axvline(fit[0][0]-tau_err,color='red',ls='--')
                    ax.axvline(fit[0][0]+tau_err,color='red',ls='--')
                    ax.axvline(fit[0][0]-2*tau_err,color='red',ls='-.')
                    ax.axvline(fit[0][0]+2*tau_err,color='red',ls='-.')
            name = ('%3d'%(self.counter+1)).replace(' ','0')
            ax1.set_xlabel('(Relative) Phase')
            ax1.set_ylabel('Negative Log Likelihood')
            pl.savefig('%s%s.png'%(self.plot_stem,name))

        self.phase_errs.append(tau_err)
        h = hm(phases) if (weights is None) else hmw(phases,weights)
        return tau,tau_err,sf_hm(h),best_ll

    def get_phase_shift(self,phases,weights,polyco_phase0):

        f   = self.__toa_loglikelihood__
        jump = False
        if self.plot_stem is not None:
            name = ('%3d'%(self.counter+1)).replace(' ','0')
            plot_output = '%s%s.png'%(self.plot_stem,name)
        else:
            plot_output = None
        # if the ephemeris is deemed good, we use it to track the 
        # solution such that we don't adopt a fluctuation for the TOA 
        if (self.good_ephemeris):
            seed = self.prev_peak
        else:
            seed = None
        x0,x0_err,best_ll = profile_analysis(
            f,(phases,weights),pred_phase=seed,plot_output=plot_output)
        if x0_err < 1e2:
            self.prev_peak = x0
        else:
            x0 = self.prev_peak # track solution with "fake" TOA
        peak_shift = (x0 - self.phi0)
        tau     = (peak_shift - polyco_phase0)
        tau_err = x0_err
        if self.display: 
            print 'Peak Shift: %.5f +/- %.5f'%(peak_shift,tau_err)
        self.phases.append(peak_shift)
        self.phase_errs.append(tau_err)
        h = hm(phases) if (weights is None) else hmw(phases,weights)
        return tau,tau_err,sf_hm(h),best_ll


class BinnedTOAGenerator(TOAGenerator):

    def init(self):
        self.profile_only = False
        self.nbins        = 32
        self.bins         = np.linspace(0,1,self.nbins + 1)

    def _prof_chisq(self,p):
        """
        Compute the Reduced Chi^2 of a binned profile p.  Each bin should hold
        the number of counts in this bin.
        Returns a tuple containing both the reduce Chi^2 and the probability
        that this value occurred by change.
        """
        # Compute mean rate
        mean = p.mean()
        # Compute chi squared
        chisq = ((p-mean)**2/mean).sum()
        # Return reduced chi squared.  DOF is nbins-1 because we are fitting
        # for the mean rate.
        dof = len(p)-1
        redchi = chisq/dof
        return (redchi,scipy.stats.chisqprob(chisq,dof))

    def _measure_phase(self,profile, template, rotate_prof=True):
        """
        measure_phase(profile, template):
            Call FFTFIT on the profile and template to determine the
                following parameters: shift,eshift,snr,esnr,b,errb,ngood
                (returned as a tuple).  These are defined as in Taylor's
                talk at the Royal Society.
        """
        c,amp,pha = fftfit.cprof(template)
        print "template", template
        pha1 = pha[0]
        if (rotate_prof):
            pha = np.fmod(pha-np.arange(1,len(pha)+1)*pha1,2.0*np.pi)
        print "profile",len(profile),profile
        print "amp", amp
        print "pha ",pha
        shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
        if 0:
            c2,amp2,pha2=fftfit.cprof(profile)
            pl.ion()
            pl.figure(1)
            pl.title('Amplitudes')
            pl.plot(amp/amp.sum())
            pl.plot(amp2/amp2.sum())
            pl.figure(2)
            pl.title('Phases')
            pl.plot(pha)
            pl.plot(pha2)
            pl.figure(3)
            #pl.title('Profile')
            pl.plot(profile/profile.sum())
            pl.plot(template/template.sum())
            pl.draw()
            tmp = raw_input()
            pl.close('all')
            print "shift, eshift, snr = ",shift,eshift,snr
        return shift,eshift,snr,esnr,b,errb,ngood

    def get_phase_shift(self,phases,weights,polyco_phase0):
            
        # Compute phase at start of this block
        if options.profile:
            polyco_phase0 = 0.0

        else:
            phases -= polyco_phase0
            phases[phases < 0.] += 1

        profile = np.histogram(phases,bins=self.bins)[0]
        # This is older version from Matthew Kerr
        #profile = np.histogram(phases,bins=self.bins[:-1])[0]


        if not self.profile_only:
            
            # Compute shift from profile    

            # Try using FFTFIT first
            shift,eshift,snr,esnr,b,errb,ngood = self._measure_phase(profile, self.template, True)
            # tau and tau_err are the predicted phase of the pulse arrival
            tau, tau_err = shift/len(profile), eshift/len(profile)
            # Note: "error" flags are shift = 0.0 and eshift = 999.0

            # If that failed, use a time-domain correlation
            if (np.fabs(shift) < 1e-7 and np.fabs(eshift-999.0) < 1e-7):
                   print >>sys.stderr, "Warning!  Bad return from FFTFIT.  Using PRESTO correlation..."
                   # Not enough structure in the template profile for FFTFIT
                   # so use time-domain correlations instead
                   tau = psr_utils.measure_phase_corr(profile, self.template)
                   # This needs to be changed
                   tau_err = 0.1/len(profile)
            
            #tau = np.random.rand(1)[0];tau_err = np.random.rand(1)[0]*0.01 # testing
            redchi,prob = self._prof_chisq(profile)
            return tau,tau_err,prob,0


        else:
            # This doesn't work very well since it doesn't pause after the any
            # but the first plot.
            pl.ioff() #does this help the pausing problem?
            pl.figure(1)
            bbins = np.arange(2.0*len(profile))/len(profile)
            pl.bar(bbins,np.concatenate((profile,profile)),width=bbins[1])
            #bins = np.arange(2.0*len(profile)+1)/len(profile)
            #pl.plot(bins,list(profile)+list(profile)+[0.0],linestyle='steps-post')
            pl.grid(1)
            pl.xlabel('Pulse Phase')
            pl.ylabel('Number of Photons')
            #pl.title('Profile')
            pl.show()
            of = file('profile.bestprof',"w")
            print >>of, "# Profile generated by polyfold.py"
            for i,np in zip(range(len(profile)),profile):
                print >>of,"%d %d" % (i,np)
            of.close()

            return 0,0,0,0

class EDFTOAGenerator(TOAGenerator):

    def init(self):
        self.edf = EDF(self.phases)

    def get_phase_shift(self,phases,polyco_phase0):

        peak_shift,sups = find_alignment(self.edf,phases)
        pl.figure(3);pl.clf();pl.plot(sups,marker='.');pl.axvline(50);pl.show()

        tau_err = float(raw_input('Estimate error width in delta_phi:'))
        #tau_err = 0.02
        return peak_shift-polyco_phase0,tau_err,sf_hm(hm(phases)),0

def profile_analysis(logl,logl_args,pred_phase=None,nsamp=100,thresh=5,
    plot_output=None,max_jump=0.25):

    # (0) establish profile
    f = lambda x: logl([x],*logl_args)
    dom = np.linspace(0,1,nsamp+1)[:-1]
    cod = np.asarray(map(f,dom))

    # (1) find all local minima
    mask = (cod < np.roll(cod,1)) & (cod < np.roll(cod,-1))

    # (2) require that all local minima surpass a likelihood threshold
    m2 = mask & (cod < -abs(thresh))

    # (3) if no peaks satisfy conditions, allow global minimum
    if m2.sum() > 0:
        mask = m2

    # (4) if given a predicted phase, choose the peak closest to it as
    # TOA; otherwise, the global minimum
    if (mask.sum() > 1) and (pred_phase is not None):
        d1 = np.abs(dom-pred_phase)
        d2 = np.abs(pred_phase+(1-dom))
        d3 = np.abs(dom-(pred_phase-1))
        diffs = np.minimum(np.minimum(d1,d2),d3)
        idx = np.argmin(diffs[mask])
    else:
        idx = np.argmin(cod[mask])
    idx = np.arange(nsamp)[mask][idx] # index into main array

    # TODO -- something more sophisticated for 0.5 aliasing -- perhaps
    # allow less significant peaks... or possibly "search" at half the
    # frequency...

    # (5) find the minimum
    # define a shifted likelihood function to avoid phase wraps
    phi0 = dom[idx]
    g = lambda x: f(x+phi0)
    xmin = golden(g,brack=[-1./nsamp,0,1./nsamp])
    fmin = g(xmin)
    phi0 += xmin

    # (6) find the error bounds, slowly but surely
    ldiffs = cod - cod[idx] - 2 # cod[idx] vs. fmin is conservative
    h = lambda x: f(x+phi0) - fmin - 2

    rt_mask = np.roll(ldiffs,-idx) > 0
    if not np.any(rt_mask):
        rt = float(nsamp-1)/2
    else:
        rt_diff = np.arange(nsamp)[rt_mask][0]
        rt_brack = float(rt_diff)/nsamp-xmin
        if h(rt_brack) < 0:
            rt_brack *= 1.1 # hopefully deal with numerical slop
        rt = brentq(h,0,rt_brack)

    lt_mask = np.roll(ldiffs,nsamp-idx-1)[::-1] > 0
    if not np.any(lt_mask):
        lt = float(nsamp-1)/2
    else:
        lt_diff = np.arange(nsamp)[lt_mask][0]
        lt_brack = -float(lt_diff)/nsamp-xmin
        if h(lt_brack) < 0:
            lt_brack *= 1.1 # hopefully deal with numerical slop
        lt = brentq(h,lt_brack,0)

    # (7) construct TOA as the mean of the error positions
    tau = (rt+lt)/2 + phi0
    tau_err = (rt-lt)/4 # by 2 for average, by 2 again for 2->1 sigma

    if plot_output is not None:
        dom1 = dom
        cod1 = cod
        inset_delta = max(abs(rt),abs(lt))+0.01
        dom2 = np.linspace(tau-inset_delta,tau+inset_delta,30)
        cod2 = np.asarray(map(f,dom2))
        pl.figure(10); pl.clf();
        ax1 = pl.gca()
        ax1.axhline(0,color='k')
        # calculate coordinates for inset; should use transAxes?
        ax2 = pl.axes([0.2+0.4*(tau<0.5),0.60,0.25,0.25])
        for i,(dom,cod,ax) in enumerate(zip([dom1,dom2],[cod1,cod2],[ax1,ax2])):
            ax.plot(dom,cod)
            ax.axvline(phi0,color='red',ls='--')
            ax.axvline(tau,color='red')
            ax.axvline(tau-tau_err,color='green',ls='-')
            ax.axvline(tau+tau_err,color='green',ls='-')
            if pred_phase is not None:
                ax.axvline(pred_phase,color='k',ls='-')
            if i==1:
                ax.xaxis.set_major_locator(pl.matplotlib.ticker.MaxNLocator(4))
                ax.axis([dom[0],dom[-1],cod.min(),cod.min()+5])
                ax.axhline(fmin+0.5,color='blue',ls='--')
                ax.axhline(fmin+2.0,color='blue',ls='-.')
                ax.axvline(tau-tau_err,color='green',ls='--')
                ax.axvline(tau+tau_err,color='green',ls='--')
                #ax.axvline(xmin-2*tau_err,color='green',ls='-.')
                #ax.axvline(xmin+2*tau_err,color='green',ls='-.')
                pl.axvline(np.mod(rt+phi0,1),color='purple',ls='-')
                pl.axvline(np.mod(lt+phi0,1),color='purple',ls='-')
        ax1.set_xlabel('(Relative) Phase')
        ax1.set_ylabel('Negative Log Likelihood')
        pl.savefig(plot_output)

    if (m2.sum() == 0):
        tau_err = 100
    # this is to catch aliases and prevent following TOAs from having
    # incorrect seed phase
    if (abs(tau-pred_phase)>max_jump):
        tau_err = 100
        tau = pred_phase
    return tau,tau_err,fmin

