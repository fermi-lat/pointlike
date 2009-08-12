
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
