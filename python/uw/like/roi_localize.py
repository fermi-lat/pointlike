"""
Module implements localization based on both broadband spectral models and band-by-band fits.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_localize.py,v 1.6 2010/04/09 17:51:49 kerrm Exp $

author: Matthew Kerr
"""

import quadform
import numpy as N
from skymaps import SkyDir,Hep3Vector
from pointlike import DoubleVector
from pypsf import PsfOverlap

###====================================================================================================###
class ROILocalizer(object):

   def init(self):
      self.tolerance = 1e-3
      self.verbose   = False
      self.update    = False
      self.max_iteration=10
      self.bandfits  = True  #default use bandfits
      self.maxdist   = 1    #fail if try to move further than this

   def __init__(self,roi,which=0,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      self.roi,self.which = roi, which
      self.quiet = roi.quiet

      # we need a reference both to the center of the ROI, for use in calculating overlap,
      # and to the original position of the point source, which may or may not be the same
      self.rd = roi.psm.roi_dir
      self.sd = roi.psm.point_sources[which].skydir
      
      if self.bandfits:
         if 'energy_bands' not in roi.__dict__.keys(): roi.setup_energy_bands()
         for eb in roi.energy_bands: eb.bandFit(which=which,saveto='bandfits')
         if N.all(N.asarray([b.bandfits for b in roi.bands]) < 0):
            if not self.quiet: print 'Warning! No good band fits.  Reverting to broadband fit...'
            self.bandfits = False
      self.tsref=0
      self.tsref = self.TSmap(self.rd)

   def TSmap(self,skydir):
      return (-2)*self.spatialLikelihood(skydir,which=self.which,bandfits=self.bandfits)-self.tsref

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
      roi   = self.roi
      which = self.which
      bandfits = self.bandfits
      verbose  = self.verbose
      update   = self.update
      tolerance= self.tolerance
      l   = quadform.Localize(self,verbose = verbose)
      ld  = SkyDir(l.dir.ra(),l.dir.dec())
      ps  = roi.psm.point_sources[which]


      ll0 = self.spatialLikelihood(ps.skydir,update=False,which=which,bandfits=bandfits)

      if not self.quiet:
         fmt ='Localizing source %s, tolerance=%.1e...\n\t'+7*'%10s'
         tup = (ps.name, tolerance,)+tuple('moved delta ra    dec   a    b  qual'.split())
         print fmt % tup
         print ('\t'+4*'%10.4f')% (0,0,ps.skydir.ra(), ps.skydir.dec())
         diff = l.dir.difference(ps.skydir)*180/N.pi
         print ('\t'+7*'%10.4f')% (diff,diff, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
      
      old_sigma=1.0
      for i in xrange(self.max_iteration):
         try:

            l.fit(update=True)
         except:
            #raise
            l.recenter()
            if not self.quiet: print 'trying a recenter...'
            continue
         diff = l.dir.difference(ld)*180/N.pi
         delt = l.dir.difference(ps.skydir)*180/N.pi
         sigma = l.par[3]
         if not self.quiet: print ('\t'+7*'%10.4f')% (diff, delt, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
         if delt>self.maxdist:
            if not self.quiet: print '\t -attempt to move beyond maxdist=%.1f' % self.maxdist
            raise Exception('roi_localize failure: -attempt to move beyond maxdist=%.1f' % self.maxdist)
         if (diff < tolerance) and (abs(sigma-old_sigma) < tolerance):
            break
         ld = SkyDir(l.dir.ra(),l.dir.dec())
         old_sigma=sigma

      if update:
         ps.skydir = l.dir
         
      ll1 = self.spatialLikelihood(l.dir,update=update,which=which,bandfits=bandfits)
      if not self.quiet: print 'TS change: %.2f'%(2*(ll0 - ll1))

      roi.qform   = l
      roi.ldir    = l.dir
      roi.lsigma  = l.sigma
      roi.delta_loc_logl = (ll0 - ll1)
      return l.dir, i, delt, 2*(ll0-ll1)

   def spatialLikelihood(self,skydir,update=False,which=0,bandfits=False):
      """Calculate log likelihood as a function of position a point source.
      
         which   -- index of point source; default to central
                    ***if localizing non-central, ensure ROI is large enough!***
      """
      
      ro  = PsfOverlap()
      roi = self.roi
      rd  = self.rd
      sd  = self.sd 
      ll  = 0
      pf  = roi.phase_factor
      bf  = self.bandfits
      wh  = self.which
      dv  = DoubleVector()

      for i,band in enumerate(roi.bands):

         en,exp,pa             = band.e,band.exp,band.b.pixelArea()
         exposure_ratio        = exp.value(skydir,en)/exp.value(sd,en) # changed to "sd" Feb. 6 2010
         
         nover                 = ro(band,rd,skydir) * exposure_ratio * band.solid_angle / band.solid_angle_p
         oover                 = band.overlaps[wh]
         psnc                  = band.bandfits if bf else band.ps_counts[wh]
         psoc                  = band.ps_counts[wh]

         if psnc < 0: continue # skip potentially bad band fits, or bands without appreciable flux

         tot_term              = (band.bg_all_counts + band.ps_all_counts + psnc * nover - psoc * oover ) * pf

         if band.has_pixels:

            #ps_pix_counts = band.psf(N.asarray([skydir.difference(x) for x in band.wsdl]),density=True)*pa*exposure_ratio
            band.wsdl.arclength(skydir,dv)
            ps_pix_counts = band.psf(N.fromiter(dv,dtype=float),density=True)*pa*exposure_ratio

            pix_term = (band.pix_counts * N.log
                           (
                              band.bg_all_pix_counts + 
                              band.ps_all_pix_counts -
                              psoc * band.ps_pix_counts[:,wh] +
                              psnc * ps_pix_counts
                           )
                       ).sum()

         else:
            pix_term = 0

         ll += tot_term - pix_term
         if N.isnan(ll):
            raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f, band %d' %(skydir.ra(),skydir.dec(),i))

         if update:
            # update the cached counts with new location -- note that this used the _broadband_ spectral
            # model rather than the band-by-band fit; otherwise, subsequent fits for broadband parameters
            # would fail
            band.overlaps[which] = nover
            band.ps_all_counts  += psoc * (nover - oover)
            if band.has_pixels:
               band.ps_all_pix_counts   += psoc * (ps_pix_counts - band.ps_pix_counts[:,wh])                              
               band.ps_pix_counts[:,wh]  = ps_pix_counts

      if N.isnan(ll):
         raise Exception('ROIAnalysis.spatialLikelihood failure at %.3f,%.3f' %(skydir.ra(),skydir.dec()))
      return ll

