"""A module providing tools for building light curves.

   Provides a single class, LightCurve.

   Author: Eric Wallace
"""

__version__ = "$Revision: 1.0 $".split()[0]
#$Header: $

import datetime
import cPickle
from calendar import monthrange
import numpy as np
import pylab as pl
from scipy import integrate
from skymaps import SkyDir
from uw.like.pointspec_helpers import FermiCatalog
from uw.utilities.fermitime import utc_to_met,MET
from uw.eew_roi.roi_factory import ROIFactory

class LightCurve(object):
    """A class to handle building light curves.

    The standard procedure is to first perform a fit for the background
    parameters, and then use refit a small subset of those parameters
    (generally only normalization parameters for one or a few sources)
    for each subinterval.

    The background fit is performed using all the
    data for the full period covered by the light curve, including all
    all catalog point sources within a given radius, and leaving all
    spectral parameters for all point sources and diffuse backgrounds
    free. The same set of sources is then used as background for the
    subinterval fits.

    For each source in each time bin, save the flux with upper and lower
    error bars, log likelihood at four points (the background fit values,
    the maximum for the interval, and plus and minus 1 sigma), and the TS.
    For cases where the lower flux error bar overlaps 0, or where the TS
    is less than the specified threshold (5 by default), the lower flux
    value is set to zero, and the upper to an upper limit (95% confidence
    by default).
    """

    defaults = dict(bg_interval =
                        (utc_to_met(2008,8,4),
                         utc_to_met(*datetime.datetime.now().timetuple()[:3]))
                   ,refit_radius = 2
                   ,refit_flux = 0
                   ,index_free = False
                   ,free_mask = None
                   ,timebin = 'month'
                   ,ts_threshold = 5
                   ,upper_limit_confidence = .95
                   )
    def init(self):
        self.factory_kwargs = ROIFactory.defaults.copy()
        self.factory_kwargs.update(dict(emin=150
                                       ,tstart = self.defaults['bg_interval'][0]
                                       ,tstop = self.defaults['bg_interval'][1]
                                       ,use_weighted_livetime = True
                                       ,irf = 'P6_v7_diff'
                                       ,catalogs = ['gll_psc18month_uw11b.fits']
                                       ,minROI = 12
                                       ,maxROI = 12
                                       ,exp_radius = 30
                                       ))
        self.roi_kwargs = dict(fit_emin = [150,150]
                              ,fit_emax = [2e5,2e5]
                              )
        self.catalog_kwargs = dict(free_radius = 8
                                  ,prune_radius = .1
                                  )
        self.fit_kwargs = dict(method = 'minuit'
                              ,tolerance = .0001
                              ,fit_bg_first = True
                              ,use_gradient = True
                              ,error_for_steps = False
                              )

    def __init__(self,skydir,**kwargs):
        """Create and set up a LightCurve.

        Parameters:
            bg_interval: tuple [(239500801, today)]
                Start and stop times in MET for the background fit. Also
                determines the full range of the light curve. The default
                ending value is 0:00 of the current day, in UTC.
            refit_radius: float [2.]
            refit_flux: float [0.]
            free_mask: array(Bool) [None]
                Sources within |refit radius| of the central source and with
                fluxes >= |refit_flux| will have their flux parameters left
                free in the fits for individual time bins. If provided,
                |free_mask| overrides the mask determined by |refit_flux| and
                |refit_radius|.
            index_free: Bool [False]
                If true, all spectral parameters for the central source will
                be left free in the fits for the individual time bins
            timebin: string or int ['month']
                The size of the individual time bins for the light curve. An
                integer specifies the length of a bin in days. A string can
                be one of 'day', 'week', or 'month', indicating bins of one
                day, one week, or one *calendar* month.
            ts_threshold: float [5]
                TS below which to calculate an upper limit for a source in a
                given time bin.
            upper_limit_confidence: float (0-1.) [.95]
                Confidence level for upper limit calculations

        """
        self.init()
        self.skydir = skydir
        d = self.defaults.copy()
        for kw in [self.factory_kwargs,self.catalog_kwargs,
                   self.roi_kwargs,self.fit_kwargs,d]:
            for k,v in kw.items():
                kw[k] = kwargs.pop(k,v)
        self.__dict__.update(d)
        if kwargs:
            print('\n'.join(['Unrecognized kwargs:']+
                            ["\t%s"%k for k in kwargs.keys()]+['']))
        self.bins = self._get_bins()
        self.bg_roi = self.background_fit()
        if self.free_mask is None:
            self.free_mask = self._get_free_mask()
        self.data = np.array((self.bins.shape[0],     #axis 0 = time bins
                              self.free_mask[0],  #axis 1 = free sources
                              6))                      #axis 2 = values


    def _get_bins(self):
        """Set up the time bins for the light curve."""
        start,stop = self.factory_kwargs['tstart'],self.factory_kwargs['tstop']
        if type(self.timebin)==type(1):
            bin_edges = np.arange(start,stop+1,86400*self.timebin,dtype='float')
        elif type(self.timebin)==type(''):
            if self.timebin == 'day':
                self.timebin = 1
                return self._get_bins()
            elif self.timebin == 'week':
                self.timebin = 7
                return self._get_bins()
            elif self.timebin == 'month':
                bin_edges = np.array([start],dtype='float')
                date = MET(start).time
                date += datetime.timedelta(monthrange(date.year,
                                                      date.month)[1]
                                           ,0,0)
                while utc_to_met(*date.timetuple()[:6])<stop:
                    bin_edges = np.append(bin_edges,
                                          utc_to_met(*date.timetuple()[:6]))
                    date += datetime.timedelta(monthrange(date.year,
                                                          date.month)[1]
                                               ,0,0)
                if bin_edges[-1]<stop:
                    bin_edges = np.append(bin_edges,stop)
            else:
                raise ValueError('''%s is not a valid time binning:
                                    must be "day", "week", "month",
                                    or an integer number of days.''')
        else:
            raise TypeError('Keyword argument "timebin" must be\
                             an int or a string.')
        return np.array(zip(bin_edges[:-1],bin_edges[1:]))

    def _get_free_mask(self):
        """Determine which sources to leave free in the individual fits"""
        diffs = np.degrees([ps.skydir.difference(self.skydir)
                            for ps in self.bg_roi.psm.point_sources])
        fluxes = np.asarray([m.i_flux() for m in self.bg_roi.psm.models])
        return np.logical_and(diffs<self.refit_radius,fluxes>self.refit_flux)

    def background_fit(self):
        """Perform a fit over the full time range."""
        kw = self.factory_kwargs.copy()
        kw['tstart'],kw['tstop'] = self.bg_interval
        f = ROIFactory(**kw)
        mapper = lambda x: FermiCatalog(x,**self.catalog_kwargs)
        roi = f(self.skydir,catalog_mapper = mapper,**self.roi_kwargs)
        roi.fit(**self.fit_kwargs)
        return roi

    def fit(self,start,stop,return_likelihood = False):
        """Perform a fit over the specified interval."""
        fact_kw = self.factory_kwargs.copy()
        fact_kw['tstart'],fact_kw['tstop'] = start,stop
        fit_kwargs = self.fit_kwargs.copy()
        fit_kwargs['fit_bg_first'] = False
        f = ROIFactory(**fact_kw)
        mapper = lambda x: FermiCatalog(x,**self.catalog_kwargs)
        roi = f(self.skydir,catalog_mapper = mapper,**self.roi_kwargs)
        #roi.set_parameters(self.bg_roi.parameters().copy())
        ll_0 = roi.logLikelihood(self.bg_roi.parameters().copy())
        for m in roi.bgm.models:
            m.free[:] = False
        for m in roi.psm.models:
            m.free[:] = False
        for ps in roi.psm.point_sources[self.free_mask]:
            ps.model.free[0] = True
        roi.psm.models[0].free[1] = self.index_free
        roi.fit(**fit_kwargs)
        if return_likelihood:
            return roi,ll_0,roi.logLikelihood(roi.parameters())
        else:
            return roi

    def build_light_curve(self):
        """Construct the light curve."""
        nfree = self.free_mask.sum()
        self.data = np.empty((self.bins.shape[0],nfree,11), dtype='object')
        for i,bin in enumerate(self.bins):
            roi,ll_0,ll_1 = self.fit(*bin,return_likelihood=True)
            #Not saving individual model parameters for now, but might want to
            pars = roi.parameters()
            hi_errs = (10**(pars + roi.get_free_errors()) - 10**pars)
            lo_errs = (10**(pars - roi.get_free_errors()) - 10**pars)
            ts = np.array([roi.TS(which = int(i))
                           for i in np.where(self.free_mask)[0]])
            if self.index_free:
                index = 10**pars.pop(1)
                index_hi, index_lo = hi_errs.pop(1), lo_errs.pop(1)
            dfluxes = 10**pars
            dflux_hi_errs = hi_errs
            dflux_lo_errs = lo_errs
            ifluxes = np.array([m.i_flux(error=True, two_sided=True)
                                for m in roi.psm.models[self.free_mask]])
            iflux_hi_errs = ifluxes[:,1]
            iflux_lo_errs = ifluxes[:,2]
            ifluxes = ifluxes[:,0]
            ulimit_mask = (iflux_lo_errs<=0) | (ts < self.ts_threshold)
            for j in np.where(ulimit_mask)[0]:
                ulimit = self.upper_limit(roi = roi
                                          ,which = j
                                          ,confidence = .95)
                iflux_hi_errs[j] = ulimit
                iflux_lo_errs[j] = 0
            ll_plus1sigma = roi.logLikelihood(np.log10(hi_errs))
            ll_minus1sigma = roi.logLikelihood(np.log10(lo_errs))
            self.data[i] = np.vstack([[bin[0]]*nfree
                             ,[bin[1]]*nfree
                             ,[ps.name.replace(' ','_') for ps in
                               roi.psm.point_sources[self.free_mask]]
                             ,ifluxes
                             ,iflux_hi_errs
                             ,iflux_lo_errs
                             ,[ll_0]*nfree
                             ,[ll_1]*nfree
                             ,[ll_plus1sigma]*nfree
                             ,[ll_minus1sigma]*nfree
                             ,ts
                             ]).transpose()


    def upper_limit(self,**kwargs):
        """Compute an upper limit on the source flux, by the "PDG Method"

        This method computes an upper limit on the flux of a specified source
        for a given time interval. The limit is computed by integrating the
        likelihood over the flux of the source, via Simpson's Rule, up to the
        desired percentile (confidence level). As such, it is essentially a
        Bayesian credible interval, using a uniform prior on the flux
        parameter.

        Note that the default integral limits are determined assuming that
        the relevant parameter is the normalization parameter of a PowerLaw
        model. For other models, especially PowerLawFlux, the limits should
        be specified appropriately.

        Arguments:
            roi: ROIAnalysis [None]
                ROIAnalysis object for which to compute the limit. If None,
                limit will be calculated for the full time interval covered by
                the light curve (i.e., using the background fit).
            which: integer [0]
                Index of the point source for which to compute the limit.
            confidence: float [.95]
                Desired confidence level of the upper limit.
            integral_min: float [-15]
                Lower limit of the likelihood integral *in log space*.
            integral_max: float [-8]
                Upper limit of the likelihood integral *in log space*.
            simps_points: int [10]
                Number of evaluation points *per decade* for Simpson's rule.
        """
        kw = dict(roi=None,
                  which=0,
                  confidence=0.95,
                  integral_min=-15,
                  integral_max =-8,
                  simps_points = 100)
        for k,v in kw.items():
            kw[k] = kwargs.pop(k,v)
        if kwargs:
            for k in kwargs.keys():
                print("Invalid keyword argument for LightCurve.upper_limit: %s"%k)
        if kw['roi'] is None:
            roi = self.bg_roi
        else:
            roi = kw['roi']
        params = roi.parameters().copy()
        ll_0 = roi.logLikelihood(roi.parameters())
        def like(norm):
            roi.psm.models[kw['which']].p[0] = norm
            return np.exp(ll_0-roi.logLikelihood(roi.parameters()))
        npoints = kw['simps_points'] * (kw['integral_max'] - kw['integral_min'])
        points = np.log10(np.logspace(kw['integral_min'],kw['integral_max'],npoints*2+1))
        y = np.array([like(x)*10**x for x in points])
        trapz1 = integrate.cumtrapz(y[::2])
        trapz2 = integrate.cumtrapz(y)[::2]
        cumsimps = (4*trapz2 - trapz1)/3.
        cumsimps /= cumsimps[-1]
        i1 = np.where(cumsimps<.95)[0][-1]
        i2 = np.where(cumsimps>.95)[0][0]
        x1, x2 = points[::2][i1], points[::2][i2]
        y1, y2 = cumsimps[i1], cumsimps[i2]
        #Linear interpolation should be good enough at this point
        limit = x1 + ((x2-x1)/(y2-y1))*(kw['confidence']-y1)
        roi.psm.models[0].p[0] = limit
        uflux = roi.psm.models[0].i_flux()
        roi.logLikelihood(params)
        return uflux

    def save(self,outfile = ''):
        """Save results to a fits file."""
        if outfile == '':
            if type(self.timebin)==type(int):
                timebin = '%iday'%self.timebin
            else:
                timebin = self.timebin
            name = self.bg_roi.psm.point_sources[self.free_mask][0].name
            outfile = ('%s_light_curve_%i_%i_%s.pickle'%
                       (name,self.factory_kwargs['tstart'],
                             self.factory_kwargs['tstop'],timebin))
        file_ = open(outfile,'w')
        bg_data = dict(parameters = self.bg_roi.parameters()
                      ,logl_0 = self.bg_roi.logLikelihood(self.bg_roi.parameters())
                      ,bg_interval = self.bg_interval
                      ,bg_roi_radius = (self.bg_roi.sa.minROI,self.bg_roi.sa.maxROI)
                      ,bg_free_radius = self.catalog_kwargs['free_radius']
                      ,roi_dir = (self.bg_roi.roi_dir.ra(),self.bg_roi.roi_dir.dec()))
        dict_ = dict(data=self.data,bg_data=bg_data)
        cPickle.dump(dict_,file_)
        file_.close()

    def plot(self):
        """Make a nice plot of a light curve.
        If outfile == '', make up a name; if outfile is None, don't save."""
        edges = np.append(self.bins[:,0],self.bins[-1,-1])
        x = self.bins.mean(axis=1)
        xlo = self.bins[:,0]
        xhi = self.bins[:,1]
        for i in xrange(self.free_mask.sum()):
            name = self.bg_roi.psm.point_sources[self.free_mask][i].name
            outfile = '%s_lightcurve_%i-%i.png'%(name,
                                                 self.factory_kwargs['tstart'],
                                                 self.factory_kwargs['tstop'])
            y = self.ifluxes[i]
            ylo,yhi = y-self.iflux_lo_errs[i],y+self.iflux_hi_errs[i]
            mask = (y>1e-11)&(ylo<y)&~np.isnan(ylo)
            ax = pl.gca()
            ax.plot([xlo[mask],xhi[mask]],[y[mask],y[mask]],color='r',linewidth=1)
            ax.plot([x[mask],x[mask]],[ylo[mask],yhi[mask]],color='r',linewidth=1)
            nmask = ~mask
            ax.plot([xlo[nmask],xhi[nmask]],
                    [np.array([self.upper_limit(bin = j,
                                                which = i,
                                                confidence = .68)
                     for j in np.where(nmask)[0]])]*2,
                    color='k',linewidth=1)
            ax.set_xbound(xlo[0]-86400,xhi[-1]+86400)
            ticks = ax.get_xticks()
            labels = []
            labels = [MET(t).time.isoformat()[:10] for t in ticks]
            ax.set_xticklabels(labels)
            pl.title('Light Curve for %s'%name)
            ax.set_ylabel(r'Flux (10^-6 ph cm^-2 s^-1)')
            if outfile is not None:
                pl.savefig(outfile)
            pl.delaxes(ax)



if __name__=='__main__':
    start,stop = (utc_to_met(2008,8,4),utc_to_met(2010,8,4))
    roi_dir =  SkyDir(250.745,39.810)
    lc = LightCurve(roi_dir,
                    tstart = start,
                    tstop = stop,
                    bg_interval=(start,stop),
                    refit_radius = 1,
                    timebin = 7)
    lc.build_light_curve()
    self.save()
    #lc.plot()

