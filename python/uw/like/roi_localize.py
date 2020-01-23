"""
Module implements localization based on both broadband spectral models and band-by-band fits.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_localize.py,v 1.38 2012/03/06 04:52:45 lande Exp $

author: Matthew Kerr
"""

import quadform
import numpy as np
from skymaps import SkyDir,Hep3Vector
from . import pypsf, roi_extended
from uw.utilities import rotations
from uw.utilities import keyword_options

def localizer(roi, which, **kwargs):
    """
    roi : a ROIanalysis object
    which: int or string
            if int, the index of the point source; if string, the name of a point or extended source
    --> a Localizer object
    """
    manager,index = roi.mapper(which)

    source=roi.get_source(which)

    if not hasattr(source,'skydir'):
        raise Exception("Can only localize Point and Extended Sources")

    return ROILocalizer(roi, index, **kwargs) if manager == roi.psm\
      else ROILocalizerExtended(roi, index, **kwargs)

###====================================================================================================###
class ROILocalizer(object):

    defaults = (
        ('tolerance',1e-3),
        ('verbose',False),
        ('update',False,"Update the source position after localization"),
        ('max_iteration',10,"Number of iterations"),
        ('bandfits',True,"Default use bandfits"),
        ('maxdist',1,"fail if try to move further than this")
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,which=0,**kwargs):
        keyword_options.process(self, kwargs)

        self.roi,self.which = roi, which
        self.quiet = roi.quiet

        # we need a reference both to the center of the ROI, for use in calculating overlap,
        # and to the original position of the point source, which may or may not be the same
        self.rd = roi.roi_dir

        self.set_source_info()
        if self.bandfits: self.do_bandfits()
        
        self.tsref=0
        self.tsref = self.TSmap(self.sd) # source position not necessarily ROI center

    def set_source_info(self):
        """ Set information about the point source that will be localized. """
        self.sd   = self.roi.psm.point_sources[self.which].skydir
        self.name = self.roi.psm.point_sources[self.which].name

    def do_bandfits(self):
        roi=self.roi
        if 'energy_bands' not in roi.__dict__.keys(): roi.setup_energy_bands()
        for eb in roi.energy_bands: eb.bandFit(which=self.which,saveto='bandfits')
        if np.all(np.asarray([b.bandfits for b in roi.bands]) < 0):
            if not self.quiet: print ('Warning! No good band fits.  Reverting to broadband fit...')
            self.bandfits = False

    def TSmap(self,skydir):
        return (-2)*self.spatialLikelihood(skydir)-self.tsref

    def dir(self):
        #return self.rd
        return self.sd # changed 6 Feb to account for non-central point sources

    def errorCircle(self):
        return 0.05 #initial guess

    def __call__(self,v):
        return self.TSmap(SkyDir(Hep3Vector(v[0],v[1],v[2])))


    def localize(self):
        """Localize a source using an elliptic approximation to the likelihood surface.

            return fit position, number of iterations, distance moved, delta TS
        """
        roi    = self.roi
        which = self.which
        bandfits = self.bandfits
        verbose  = self.verbose
        update    = self.update
        tolerance= self.tolerance
        l    = quadform.Localize(self,verbose = verbose)
        ld  = SkyDir(l.dir.ra(),l.dir.dec())

        ll0 = self.spatialLikelihood(self.sd,update=False)

        if not self.quiet:
            fmt ='Localizing source %s, tolerance=%.1e...\n\t'+7*'%10s'
            tup = (self.name, tolerance,)+tuple('moved delta ra     dec    a     b  qual'.split())
            print (fmt % tup)
            print (('\t'+4*'%10.4f')% (0,0,self.sd.ra(), self.sd.dec()))
            diff = l.dir.difference(self.sd)*180/np.pi
            print (('\t'+7*'%10.4f')% (diff,diff, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6]))
        
        old_sigma=1.0
        for i in xrange(self.max_iteration):
            try:

                l.fit(update=True)
            except:
                #raise
                l.recenter()
                if not self.quiet: print ('trying a recenter...')
                continue
            diff = l.dir.difference(ld)*180/np.pi
            delt = l.dir.difference(self.sd)*180/np.pi
            sigma = l.par[3]
            if not self.quiet: (print ('\t'+7*'%10.4f')% (diff, delt, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6]))
            if delt>self.maxdist:
                if not self.quiet: print ('\t -attempt to move beyond maxdist=%.1f' % self.maxdist)
                raise Exception('roi_localize failure: -attempt to move beyond maxdist=%.1f' % self.maxdist)
            if (diff < tolerance) and (abs(sigma-old_sigma) < tolerance):
                break
            ld = SkyDir(l.dir.ra(),l.dir.dec())
            old_sigma=sigma

        roi.qform    = l
        roi.ldir     = l.dir
        roi.lsigma   = l.sigma

        ll1 = self.spatialLikelihood(l.dir,update=update)
        if not self.quiet: print ('TS change: %.2f'%(2*(ll0 - ll1)))

        roi.delta_loc_logl = (ll0 - ll1)

        # this is necessary in case the fit always fails.
        delt = l.dir.difference(self.sd)*180/np.pi

        return l.dir, i, delt, 2*(ll0-ll1)

    def spatialLikelihood(self,skydir,update=False):
        """Calculate log likelihood as a function of position a point source.
        
            which    -- index of point source; default to central
                          ***if localizing non-central, ensure ROI is large enough!***
        """
        
        ro  = pypsf.PsfOverlap()
        rd  = self.rd
        roi = self.roi
        ll  = 0
        wh  = self.which

        for i,band in enumerate(roi.bands):

            # N.B. -- needs to be the ratio with ROI_DIR to be consistent with ROIPointSourceManager
            exposure_ratio  = band.exp.value(skydir,band.e)/band.exp.value(rd,band.e)
            
            nover           = ro(band,rd,skydir)
            oover           = band.overlaps[wh]
            psnc            = (band.bandfits if self.bandfits else band.ps_counts[wh])*exposure_ratio/band.er[wh]                                          
            psoc            = band.ps_counts[wh] # N.B. -- ps_counts includes exposure ratio

            if psnc < 0: continue # skip potentially bad band fits, or bands without appreciable flux

            tot_term = (band.bg_all_counts + band.ps_all_counts + psnc*nover - psoc*oover ) * band.phase_factor

            if band.has_pixels:
                
                rvals = np.empty(len(band.wsdl))
                band.psf.cpsf.wsdl_val(rvals,skydir,band.wsdl)
                ps_pix_counts = rvals*band.b.pixelArea()

                pix_term = (band.pix_counts * np.log(
                                band.bg_all_pix_counts + band.ps_all_pix_counts -
                                psoc*band.ps_pix_counts[:,wh] + psnc*ps_pix_counts
                            ) ).sum()

            else: pix_term = 0

            ll += tot_term - pix_term
            if np.isnan(ll):
                raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f, band %d' %(skydir.ra(),skydir.dec(),i))

            if update:
                # update the cached counts with new location -- note that this used the _broadband_ spectral
                # model rather than the band-by-band fit; otherwise, subsequent fits for broadband parameters
                # would fail

                # here, only update stuff set in setup_initial_counts which changes w/ localization
                band.overlaps[wh] = nover
                band.er[wh] = exposure_ratio
                if band.has_pixels: band.ps_pix_counts[:,wh]   = ps_pix_counts

        if update:
            # need to update frozen_pix_counts/unfrozen_pix_counts
            roi.psm.cache(roi.bands)
            # need to update ps_counts, ps_all_counts, ps_all_pix_counts
            roi.psm.update_counts(roi.bands)
            
            # update source position
            roi.psm.point_sources[wh].skydir = skydir

        if np.isnan(ll):
            raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f' %(skydir.ra(),skydir.dec()))
        return ll

def print_ellipse(roi, label=True, line=True):
    """ print the ellipical parameters (all deg units):
            ra, dec
            a, b  : major and minor 1-sigma axes
            ang   : ellipse orientation, E of N
            qual  : measure of fit quality.
    Optional parameters:
        label [True] print a label line
        line  [True] print the line corresponding to the current fit
              (only knows about one at a time)
    """
    self=roi
    if not self.qform: return
    labels = 'ra dec a b ang qual'.split()
    if label: print ((len(labels)*'%10s') % tuple(labels))
    if not line: return
    p = self.qform.par[0:2]+self.qform.par[3:]
    print (len(p)*'%10.4f' % tuple(p))

def get_ellipse(roi):
    """ Returns a dictionary specifying the elliptical 
        localiztion parameters. """
    d={}
    if hasattr(roi,'qform'):
        q=roi.qform.par
        d.update(
            ra=float(q[0]), dec=float(q[1]),
            a=float(q[3]), b=float(q[4]),
            ang=float(q[5]), qual=float(q[6])
            )
    if hasattr(roi,'lsigma'): d['lsigma']=roi.lsigma
    return d

class ROILocalizerExtended(ROILocalizer):

    def set_source_info(self):
        self.sd = self.roi.dsm.diffuse_sources[self.which].spatial_model.center
        self.name = self.roi.dsm.diffuse_sources[self.which].name

    def do_bandfits(self):
        if 'energy_bands' not in self.roi.__dict__.keys(): self.roi.setup_energy_bands()
        for eb in self.roi.energy_bands: 
            bfe=roi_extended.BandFitExtended(self.which,eb,self.roi)
            bfe.fit(saveto='bandfits')
        if np.all(np.asarray([b.bandfits for b in self.roi.bands]) < 0):
            if not self.quiet: print ('Warning! No good band fits.  Reverting to broadband fit...')
            self.bandfits = False

    def spatialLikelihood(self,skydir,update=False):
        """Calculate log likelihood as a function of position of an extended source.
        
            which    -- index of the extended source; default to central
                          ***if localizing non-central, ensure ROI is large enough!***
        """
        
        rd  = self.rd
        roi = self.roi
        ll  = 0
        wh  = self.which

        es = roi.dsm.bgmodels[wh]
        sm  = es.extended_source.spatial_model

        sm.modify_loc(skydir)

        for i,(band,myband) in enumerate(zip(roi.bands,es.bands)):

            es.set_state(band)

            # N.B. -- needs to be the ratio with ROI_DIR to be consistent with ROIPointSourceManager
            exposure_ratio  = band.exp.value(skydir,es.current_energy)/band.exp.value(rd,es.current_energy)
            
            nover           = es._overlaps(skydir,band)
            oover           = myband.overlaps
            esnc            = (band.bandfits if self.bandfits else myband.es_counts)*exposure_ratio/myband.er
            esoc            = myband.es_counts # N.B. -- es_counts includes exposure ratio

            if esnc < 0: continue # skip potentially bad band fits, or bands without appreciable flux

            tot_term = (band.bg_all_counts + band.ps_all_counts + esnc*nover - esoc*oover ) * band.phase_factor

            if band.has_pixels:
                
                es_pix_counts = es._pix_value(band.wsdl)*band.b.pixelArea()

                pix_term = (band.pix_counts * np.log(
                                band.bg_all_pix_counts + band.ps_all_pix_counts -
                                esoc*myband.es_pix_counts + esnc*es_pix_counts
                            ) ).sum()

            else: pix_term = 0

            ll += tot_term - pix_term
            if np.isnan(ll):
                raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f, band %d' %(skydir.ra(),skydir.dec(),i))


        if update:
            # update the cached counts with new location -- note that this used the _broadband_ spectral
            # model rather than the band-by-band fit; otherwise, subsequent fits for broadband parameters
            # would fail
            es.initialize_counts(roi.bands)
            roi.update_counts()

            # Update covariance matrix of the extended source with the best fit error.
            # Note that this update doesn't deal with correlation in the fit parameters,
            # but should be good enough for the purpose of displaying sources.
            # Since we can't fit lon/lat in log space, filling to covaraince matrix is easy.
            sm.cov_matrix[0:2,:] = 0
            sm.cov_matrix[:,0:2] = 0
            sm.cov_matrix[0][0] = sm.cov_matrix[1][1] = roi.lsigma**2

        else:
            sm.modify_loc(self.sd)
            # Note that since initialize_counts wasn't previous called,
            # all of the previous values were still set at the correct
            # values, so there is no reason to reinitialize back to the
            # original values.

            # You have to set_state to redo the convolution. Note that
            # since caching has been enabled, all this really does is
            # reset the center of the extended source for the convolution.
            for band in roi.bands: es.set_state(band)

        return ll


class DualLocalizer():
    """ Fit two point sources at the same time. Fit the center of position
        and relative difference since they are more robust parameters.

        This method is only suitable for sources relativly nearby (~<1 degree).

        Note, bandfits is not allowed for fitting because the bandfits
        algorithm does not really work when fitting two really nearby
        sources. """

    defaults = (
            ('use_gradient',     True, "Analytic gradient of spectral parameters when fitting."),
            ('tolerance',        0.01, "Fit tolerance to use when fitting"),
            ('verbose',          True, "Print more stuff during fit.")
    )

    @staticmethod
    def symmetric_mod(a,b):
        """ Just like mod, but is symmetric for + and - values
            
            >>> DualLocalizer.symmetric_mod(10, 360)
            10
            >>> DualLocalizer.symmetric_mod(-10, 360)
            -10
            >>> DualLocalizer.symmetric_mod(370, 360)
            10
            >>> DualLocalizer.symmetric_mod(-370, 360)
            -10
        """
        temp=a%b
        if temp > b/2.: temp-=b
        return temp

    @staticmethod
    def approx_mid_point(skydir1,skydir2):
        """ This method is only valid in the small distance limit, in which
            case space is approximatly flat and the rhomb line is equal
            to the great circle line. """
        if abs(skydir1.b()) < abs(skydir1.dec()):
            return SkyDir(skydir1.l() + DualLocalizer.symmetric_mod(skydir2.l()-skydir1.l(),360)/2,
                          skydir1.b() + DualLocalizer.symmetric_mod(skydir2.b()-skydir1.b(),360)/2,
                          SkyDir.GALACTIC)
        else:
            return SkyDir(skydir1.ra() + DualLocalizer.symmetric_mod(skydir2.ra()-skydir1.ra(),360)/2,
                          skydir1.dec() + DualLocalizer.symmetric_mod(skydir2.dec()-skydir1.dec(),360)/2,
                          SkyDir.EQUATORIAL)


    @keyword_options.decorate(defaults)
    def __init__(self, roi, which1, which2, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi

        self.p1 =roi.get_source(which1)
        self.p2 =roi.get_source(which2)

    @staticmethod
    def print_flux(source,roi):
        return source.model.i_flux(emin=min(roi.fit_emin),emax=max(roi.fit_emax),e_weight=0)

    def fit(self,p):
        m_x,m_y,d_x,d_y=p
        roi=self.roi

        s2 = SkyDir(m_x+d_x,m_y+d_y)
        s1 = SkyDir(m_x-d_x,m_y-d_y)

        rot_back_1=rotations.anti_rotate_equator(s1,self.middle)
        rot_back_2=rotations.anti_rotate_equator(s2,self.middle)

        roi.modify(which=self.p1,skydir=rot_back_1)
        roi.modify(which=self.p2,skydir=rot_back_2)

        ll=roi.fit(use_gradient=self.use_gradient,estimate_errors=False)

        if ll < self.ll_best:
            prev = roi.parameters().copy()

            roi.set_parameters(self.best_spectral.copy())
            roi.__update_state__()

            ll_alt=roi.fit(use_gradient=self.use_gradient,estimate_errors=False)

            if ll_alt > ll: 
                ll=ll_alt
            else: 
                roi.set_parameters(prev.copy())
                roi.__update_state__()

        if ll < self.ll_0:

            prev = roi.parameters().copy()

            roi.set_parameters(self.init_spectral.copy())
            roi.__update_state__()

            ll_alt=roi.fit(use_gradient=self.use_gradient,estimate_errors=False)

            if ll_alt > ll: 
                ll=ll_alt
            else: 
                roi.set_parameters(prev.copy())
                roi.__update_state__()

        if ll > self.ll_best: 
            self.ll_best = ll
            self.best_spectral = roi.parameters().copy()

        if self.verbose: print ('d=%s f=%.1e, d2=%s, f=%.1e, dist=%.3f logL=%.3f dlogL=%.3f' % \)
                (rot_back_1, rotations.print_flux(self.p1,roi), 
                 rot_back_2, rotations.print_flux(self.p2,roi), 
                 np.degrees(rot_back_1.difference(rot_back_2)),
                 ll,ll-self.ll_0)

        return -ll # minimize negative log likelihood

    def localize(self):
        roi=self.roi

        p1,p2=self.p1, self.p2

        if not roi.quiet: print ('Dual localizing source %s and %s' % (self.p1.name,self.p2.name))

        d1=self.p1.skydir
        d2=self.p2.skydir

        self.middle=rotations.approx_mid_point(d1,d2)

        # Points rotated so that the middle is at the equator
        rot1=rotations.rotate_equator(d1,self.middle)
        rot2=rotations.rotate_equator(d2,self.middle)

        # Fit average point and distance between them
        # Wrap coordiantes to vary between -180 and 180
        x1 = rot1.ra(); x1 = x1 - (x1>180)*360
        x2 = rot2.ra(); x2 = x2 - (x2>180)*360
        y1 = rot1.dec(); y2 = rot2.dec()

        m_x,m_y = (x2+x1)/2, (y2+y1)/2
        d_x,d_y = (x2-x1)/2, (y2-y1)/2

        p0 = [m_x,m_y,d_x,d_y]

        self.ll_0=self.ll_best=-1*roi.logLikelihood(roi.parameters())

        self.init_spectral = self.best_spectral = roi.parameters().copy()

        old_quiet= roi.quiet
        roi.quiet=True

        from uw.utilities.minuit import Minuit
        steps=[0.1,0.1,0.1,0.1] # expect to fit sources ~ 0.1 degrees away.
        m = Minuit(self.fit,p0,
                   tolerance = self.tolerance,
                   maxcalls  = 500,
                   printMode = True, 
                   steps     = steps)

        best_spatial,fval = m.minimize(method="SIMPLEX")

        roi.quiet = old_quiet

        return

if __name__ == "__main__":
    import doctest
    doctest.testmod()
