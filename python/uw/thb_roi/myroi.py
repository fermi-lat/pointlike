"""
User interface to SpectralAnalysis
----------------------------------
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/myroi.py,v 1.2 2010/03/18 20:36:18 burnett Exp $

"""

import numpy as np
import numpy as N  # for Kerr compatibility
from scipy import optimize
import pylab as plt
import os, pickle, math 

from uw.like import roi_analysis, roi_localize, roi_plotting,  Models
from uw.utilities import makerec, fermitime, image
from skymaps import SkyDir,  PySkyFunction


def spectralString(band,which=None):
  """Return a string suitable for printSpectrum.

     which -- if not None, an array of indices for point sources to fit
  """
  self=band
  which = which or [0]
  r     = []
  for w in which:
     try:
        self.bandFit(which=w)
     except np.linalg.LinAlgError:
        return 'singular matrix'
     self.m.p[0] = N.log10(self.uflux)
     ul = sum( (b.expected(self.m) for b in self.bands) )
     if self.flux is None:
        r += [0,ul,0]
     else:
        n = ul*self.flux/self.uflux
        r += [n,ul - n, n - ul*self.lflux/self.uflux]
     r += [self.ts]

  rois = self.__rois__()
  ph   = int(sum( (b.photons for b in self.bands) ))
  gal  = sum( (b.bg_counts[0] for b in self.bands) )
  iso  = sum( (b.bg_counts[1] for b in self.bands) )
  values = tuple([int(round(self.emin)),rois[0],rois[1],ph,gal,iso, ph-gal-iso] + r)
  format = '  '.join(['%6i','%6.1f','%6.1f','%7i','%8.1f','%9.1f','%6.1f']+['%7.1f +%5.1f -%5.1f %6.1f' ]*len(which))
  return format%values


class MyROI(roi_analysis.ROIAnalysis):
    """ create an ROIAnalysis subclass 

        Parameters

        ps_manager :  PointSourceManager
        bg_manager :  BackgroundManager
        roifactory :  ROIfactory

    """

    def __init__(self,ps_manager,bg_manager,roifactory,**kwargs):
        """
        Parameters

        ps_manager :  
        bg_manager : 
        roifactory : 
       
        change default free_radius to turn this off

        optional:
            bgfree [True, False, True]
            prune_radius [0.1]
            free_radius  [0]

        """
        bgfree = [True,False,True]
        if 'bgfree' in kwargs:
            bgfree = np.array(kwargs['bgfree'])
            kwargs.pop('bgfree')
        if 'fit_bg_first' in kwargs:
            self.fit_bg_first = kwargs['fit_bg_first']
        super(MyROI, self).__init__(ps_manager, bg_manager, roifactory, **kwargs)
        self.bgm.models[0].free = np.array(bgfree[:2])
        self.bgm.models[1].free = np.array([bgfree[2]])
        self.name = self.psm.point_sources[0].name # default name
        self.center= self.sa.roi_dir

    def fit(self, **kwargs):
        """ invoke base class fitter, but insert defaults first 
        """
        fit_bg_first = self.fit_bg_first
        if 'fit_bg_first' in kwargs: 
            fit_bg_first = kwargs.pop('fit_bg_first')
        if 'use_gradient' not in kwargs: kwargs['use_gradient']=self.use_gradient
        ret = super(MyROI, self).fit(fit_bg_first=fit_bg_first, **kwargs)
        if not self.quiet: print self

    def band_ts(self, which=0):
        """ return the sum of the individual band ts values
        """
        self.setup_energy_bands()
        ts = 0
        for eb in self.energy_bands:
            eb.bandFit(which)
            ts += eb.ts
        return ts

    def localize(self,which=0, tolerance=1e-3,update=False, verbose=False, bandfits=True):
        """Localize a source using an elliptic approximation to the likelihood surface.

          which     -- index of point source; default to central 
                      **if localizing non-central, ensure ROI is large enough!**
          tolerance -- maximum difference in degrees between two successive best fit positions
          update    -- if True, update localization internally, i.e., recalculate point source contribution
          bandfits  -- if True, use a band-by-band (model independent) spectral fit; otherwise, use broabband fit

         return fit position
        """
        try:
            loc = super(MyROI,self).localize(which=which,bandfits=bandfits,tolerance=tolerance,update=update,verbose=verbose)
            if not self.quiet: self.print_ellipse()
        except:
            #raise
            self.qform=None
            loc = None 
        self.find_tsmax()
        return loc

    def print_ellipse(self, label=True, line=True):
        if not self.qform: return
        labels = 'ra dec a b phi qual'.split()
        if label: print (len(labels)*'%10s') % tuple(labels)
        if not line: return
        p = self.qform.par[0:2]+self.qform.par[3:]
        print len(p)*'%10.3f' % tuple(p)

    def find_tsmax(self, bandfits=True):
        """ 
            very simple function that looks for maximum position
        """
        class TSfun(object):
            """ helper class """
            def __init__(self, roi, which=0, bandfits=bandfits):
                self.tsf=roi.tsmap(which,bandfits=bandfits)
                self.sdir = roi.center
                self.ra,self.dec = self.sdir.ra(), self.sdir.dec()
                self.cdec= math.cos(math.degrees(self.dec))
            def __call__(self,par):
                ra = self.ra+par[0]/self.cdec
                dec= self.dec+par[1]
                return -self.tsf(SkyDir(ra,dec))
            def maximize(self):
                dx,dy = optimize.fmin(self, (0,0),disp=0)
                return SkyDir(self.ra+dx/self.cdec, self.dec+dy)

        tsf = TSfun(self)
        
        self.tsmax= tsf.maximize()
        return self.tsmax
        

    def dump(self, sdir=None, galactic=False, maxdist=5, title=''):
        """ formatted table point sources positions and parameter in the ROI
        Parameters
        ----------
            sdir : SkyDir, optional, default None for center
                for center: default will be the first source
            galactic : bool, optional, default False
               set true for l,b
            maxdist : float, optional, default 5
               radius in degrees

        """
        self.print_summary(sdir, galactic, maxdist, title)

    

    def tsmap(self, which=0, bandfits=True):
        """ return function of likelihood in neighborhood of given source
            tsm = roi.tsmap(which)
            size=0.25
            tsp = image.TSplot(tsm, center, size, pixelsize =size/20, axes=plt.gca())
            tsp.plot(center, label=name)
            tsp.show()

        """
        self.localizer = roi_localize.ROILocalizer(self, which, bandfits=bandfits)
        return PySkyFunction(self.localizer)


    def plot_spectra(self, which=0, axes=None,axis=None, outfile=None, fignum=1, **kwargs):
        """generate a spectral plot
        sets up a figure, if fignum is specified.
        """

        if axes is None:
            fig=plt.figure(fignum);
            plt.clf()
            axes = plt.gca()
        return roi_plotting.band_fluxes(self, which,axes, axis,outfile, **kwargs)

    def pickle(self, name, outdir, **kwargs):
        """ name: name for source, used as filename
            outdir: ouput directory
        """
        name = name.strip()
        output = dict()
        output['name'] = name
        output['ra']   = self.center.ra()
        output['dec']  = self.center.dec()
        output['src_par'] = 10**self.psm.models[0].p
        output['bgm_par'] = np.hstack((10**self.bgm.models[0].p, 10**self.bgm.models[1].p))
        output['qform_par'] = self.qform.par if self.qform is not None else None
        output['tsmax'] = None if 'tsmax' not in self.__dict__ else [self.tsmax.ra(),self.tsmax.dec()]
        output.update(kwargs) # add additional entries from kwargs

        if not os.path.exists(outdir):
            os.mkdir(outdir)
        f = file(os.path.join(outdir,name+'.pickle'),'wb')
        pickle.dump(output,f)
        f.close()

    def printSpectrum(self,sources=None):
        """Print total counts and estimated signal in each band for a list of sources.

        Sources can be specified as PointSource objects, source names, or integers
        to be interpreted as indices for the list of point sources in the roi. If
        only one source is desired, it needn't be specified as a list. If no sources
        are specified, all sources with free fit parameters will be used."""
        if sources is None:
         sources = [s for s in self.psm.point_sources if N.any(s.model.free)]
        elif type(sources) != type([]): 
         sources = [sources]
        bad_sources = []
        for i,s in enumerate(sources):
         if type(s) == roi_analysis.PointSource:
            if not s in self.psm.point_sources:
               print 'Source not found in source list:\n%s\n'%s
               bad_sources += [s]
         elif type(s) == int:
            try:
               sources[i] = self.psm.point_sources[s]
            except IndexError:
               print 'No source #%i. Only %i source(s) specified.'\
                     %(s,len(self.psm.point_sources))
               bad_sources += [s]
         elif type(s) == type(''):
            names = [ps.name for ps in self.psm.point_sources]
            try:
               sources[i] = self.psm.point_sources[names.index(s)]
            except ValueError:
               print 'No source named %s'%s
               bad_sources += [s]            
         else:
            print 'Unrecognized source specification:', s
            bad_sources += [s]
        sources = set([s for s in sources if not s in bad_sources])
        indices = [self.psm.point_sources.index(s) for s in sources]
        self.setup_energy_bands()

        fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic','Excess']\
                +[' '*15+'Signal']*len(sources)
        outstring = 'Spectra of sources in ROI about %s at ra = %.2f, dec = %.2f\n'\
                    %(self.psm.point_sources[0].name, self.center.ra(), self.center.dec())
        outstring += ' '*54+'  '.join(['%21s'%s.name for s in sources])+'\n'
        outstring += '  '.join(fields)+'\n'
        print outstring
        for eb in self.energy_bands:
        #         print eb.spectralString(which=indices)
         # local
            print spectralString(eb, which=indices)

    def get_spectrum(self):
        """
        return a recarry with spectral details
        """
        
        fields = 'emin  f_ROI b_ROI events galactic isotropic excess'.split();
        rec = makerec.RecArray(fields)
        def get_band_info(band):
            """ copied from spectralString """
            self = band
            rois = self.__rois__()
            ph   = int(sum( (b.photons for b in self.bands) ))
            gal  = sum( (b.bg_counts[0] for b in self.bands) )
            iso  = sum( (b.bg_counts[1] for b in self.bands) )
            values = [int(round(self.emin)),rois[0],rois[1],ph,gal,iso, ph-gal-iso]
            return values
        self.setup_energy_bands()
        for eb in self.energy_bands:
            rec.append( *get_band_info(eb) )
        return rec()

    def printMySpectrum(self,sources=None, out=None):
        """Print total counts and estimated signal in each band for a list of sources.

        Sources can be specified as PointSource objects, source names, or integers
        to be interpreted as indices for the list of point sources in the roi. If
        only one source is desired, it needn't be specified as a list. If no sources
        are specified, all sources with free fit parameters will be used.
        
        """
        #super(MyROI,self).printSpectrum(sources)
        if sources is None:
         sources = [s for s in self.psm.point_sources if N.any(s.model.free)]
        elif type(sources) != type([]): 
         sources = [sources]
        bad_sources = []
        for i,s in enumerate(sources):
            if type(s) == roi_analysis.PointSource:
               if not s in self.psm.point_sources:
                  print 'Source not found in source list:\n%s\n'%s
                  bad_sources += [s]
            elif type(s) == int:
               try:
                  sources[i] = self.psm.point_sources[s]
               except IndexError:
                  print 'No source #%i. Only %i source(s) specified.'\
                        %(s,len(self.psm.point_sources))
                  bad_sources += [s]
            elif type(s) == type(''):
               names = [ps.name for ps in self.psm.point_sources]
               try:
                  sources[i] = self.psm.point_sources[names.index(s)]
               except ValueError:
                  print 'No source named %s'%s
                  bad_sources += [s]            
            else:
               print 'Unrecognized source specification:', s
               bad_sources += [s]
        sources = set([s for s in sources if not s in bad_sources])
        spectra = dict.fromkeys(sources)
        for s in sources:
         spectra[s] = roi_plotting.band_spectra(self,source=self.psm.point_sources.index(s))
        iso,gal,src,obs = roi_plotting.counts(self)[1:5]
        fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic','   excess']\
                +[' '*10+'Signal']*len(sources)
        outstring = 'Spectra of sources in ROI about %s at ra = %.3f, dec = %.3f\n'\
                    %(self.psm.point_sources[0].name, self.center.ra(), self.center.dec())\

        outstring += ' '*68+'  '.join(['%-18s'%s.name for s in sources])+'\n'
        outstring += '  '.join(fields)+'\n'
        for i,band in enumerate(zip(self.bands[::2],self.bands[1::2])):
            values = (band[0].emin, band[0].radius_in_rad*180/N.pi,
                      band[1].radius_in_rad*180/N.pi,obs[i],gal[i],iso[i], obs[i]-gal[i]-iso[i])
            for s in sources:
                values+=(spectra[s][1][i],.5*(spectra[s][3][i]-spectra[s][2][i]))          
            string = '  '.join(['%6i','%6.2f','%6.2f','%7i','%8.1f','%9.1f','%9.1f']+
                               ['%8.1f +/-%6.1f']*len(sources))%values
            outstring += string+'\n'
        print >>out, outstring

    def get_photons(self, emin=1000): 
        """ access binned photon data

        """
        roi = self
        ra,dec = roi.center.ra(), roi.center.dec() 
        energy=[]; etype=[]; dra=[]; ddec=[]

        cosdec=math.cos(math.radians(dec))
        for band in roi.bands:
            e = band.e
            if e<emin: continue
            t = band.b.event_class() & 3 # note: mask off high bits!
            for wsd in band.wsdl:
                for i in range(wsd.weight()):
                    energy.append(e)
                    etype.append(t)
                    x =wsd.ra()-ra
                    if x<-180: x+=360
                    dra.append(x*cosdec)
                    ddec.append(wsd.dec()-dec)
                    #print '%6.0f %8.3f %8.3f' % (e, (wsd.ra()-ra)*cosdec,wsd.dec()-dec)

        return np.rec.fromarrays([energy,etype,dra,ddec],
                           names='energy etype dra ddec'.split())         

    def plot_tsmap(self, name=None, center=None, size=0.5, pixelsize=None, outdir=None, 
            which=0, catsig=99, axes=None, fignum=99, 
            bandfits=True,
            galmap=True, galactic=False,
            assoc = None,
            notitle = False,
            nolegend = False,
            markersize=12,
            primary_markersize=14,
            ): #, **kwargs):
        """ create a TS map for the source

        Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  name        [None]  -- provide name for title, and to save figure if outdir set
  center      [None] -- center, default the roi center
  outdir      [None] -- folder name to save the image in, outdir/<name>_tsmaps.png
  catsig      [99]  -- if set and less than 1.0, draw cross with this size (degrees)
  size        [0.5]  -- half width=height (deg)
  pixelsize   [None] -- if not set, will be 20 x20 pixels
  galmap      [True] -- if set, draw a galactic coordinate image with the source position shown
  galactic    [False] -- plot using galactic coordinates
  which       [0]    -- chose a different source in the ROI to plot
  assoc       [None] -- if set, a list of tuple of associated sources 
  notitle     [False] -- set to turn off (allows setting the current Axes object title)
  nolegend    [False]
  markersize  [12]
  =========   =======================================================

        returns the image.TSplot object for plotting positions, for example
        """
        roi = self
        kwargs={} #fix later
        name = self.name if name is None else name
        tsm = roi.tsmap(which=which, bandfits=bandfits)
        sdir = center if center is not None else self.center
        if axes is None: 
            plt.figure(fignum,figsize=(5,5)); plt.clf()
        
        tsp = image.TSplot(tsm, sdir, size, pixelsize =pixelsize if pixelsize is not None else size/20. , 
                    axes=axes, galactic=galactic, galmap=galmap)
        if 'qform' in roi.__dict__ and roi.qform is not None:
            sigma = math.sqrt(roi.qform.par[3]*roi.qform.par[4]) # why do I need this?
            qual = roi.qform.par[6]
            if sigma<1 and qual <50:
                tsp.overplot(roi.qform, sigma)
            else:
                print 'bad fit sigma %g, >1 or qual %.1f >50' % (sigma, qual)
        tsp.show(colorbar=False)
        if catsig<1:
            tsp.cross(sdir, catsig, lw=2, color='grey')
            
        # plot the primary source, any nearby from the fit
        x,y = tsp.zea.pixel(sdir)
        tsp.zea.axes.plot([x],[y], '*', color='grey', label=name, markersize=primary_markersize)
        marker = 'ov^<>1234sphH'; i=k=0
        for ps in self.psm.point_sources: # skip 
            x,y = tsp.zea.pixel(ps.skydir)
            if ps.name==name or x<0 or x>tsp.zea.nx or y<0 or y>tsp.zea.ny: continue
            tsp.zea.axes.plot([x],[y], marker[k%12], color='blue', label=ps.name, markersize=markersize)
            k+=1
        
        if 'tsmax' in self.__dict__ and self.tsmax is not None:
            tsp.plot(self.tsmax, symbol='x')
        if not notitle: plt.title( name, fontsize=12)

        if assoc is not None:
            # eventually move this to image.TSplot
            last_loc=SkyDir(0,90)
            for aname, loc, prob, catid in zip(assoc['name'],assoc['dir'],assoc['prob'],assoc['cat']):
                print 'associate with %s, prob=%.2f' % (aname.strip(),prob)
                if catid in (1,6,9,10, 11): 
                    print '---skip gamma cat #%d' % catid
                    continue
                if i>8:
                    print '---skip because too many for display'
                    continue
                x,y = tsp.zea.pixel(loc)
                diff = np.degrees(loc.difference(last_loc)); last_loc=loc
                if diff>1e-3: k+=1 # new marker only if changed place
                tsp.zea.axes.plot([x], [y], marker=marker[k%12], color='green', linestyle='None',
                    label='%s[%d] %.2f'%(aname.strip(), catid, prob ), markersize=markersize)
                i+=1
        
        fs = plt.rcParams['font.size']
        plt.rcParams.update({'legend.fontsize':8, 'font.size':8})
        # put legend on left.
        if not nolegend: tsp.zea.axes.legend(loc=2, numpoints=1, bbox_to_anchor=(-0.15,1.0))
        plt.rcParams['font.size'] = fs

        #if galmap: # moved
        #    axi = plt.gcf().add_axes((0.75, 0.80, 0.20, 0.10))
        #    ait_insert=image.AIT_grid(axes=axi, labels=False, color='w')
        #    ait_insert.plot([self.center], 'sr')
        
        if outdir is not None: plt.savefig(os.path.join(outdir,'%s_tsmap.png'%name.strip()))
        return tsp

    def plot_FITS(self, fitsfile,  size=0.5, bandfits=True,):
        """ kluge to make a fits file
        """
        z = image.ZEA(self.center, size=size, pixelsize=size/20, fitsfile=fitsfile)
        z.fill(self.tsmap(bandfits=bandfits))
        del(z) 

    def rescan(self, tsm, threshold=1):
        """ scan the tsmap for a secondary peak
        tsm is a map generated by tsm = roi.plot_tsmap...
        return position of peak position, if delta TS >threshold
        """
        M = tsm.zea.image
        self.tsmap_max = M.max()
        maxpixel = None
        if self.tsmap_max >threshold:
            # there is a secondary peak
            t = np.where(M==M.max())
            x,y=t[1][0],t[0][0]
            maxpixel = tsm.zea.skydir(x,y)
            if not self.quiet: print 'found maximum %.1f at (%.3f,%.3f)'\
                % (self.tsmap_max, maxpixel.ra(),maxpixel.dec())
        return maxpixel
        
    def plot_sed(self, fignum=5, axes=None,
            axis=(1e2,1e5,1e-8,1e-2),
            data_kwargs={'linewidth':2, 'color':'k',},
            fit_kwargs={'lw':2, 'color':'r',},
            ):
        """Plot a SED, a thin interface to uw.like.roi_plotting.make_sed
        ========     ===================================================
        keyword      description
        ========     ===================================================
        fignum       [5] if set, use (and clear) this figure. If None, use current Axes object
        axes         [None] If set use this Axes object
        axis         (1e2, 1e5, 1e-8, 1e-2)
        data_kwargs  a dict to pass to the data part of the display
        fit_kwargs   a dict to pass to the fit part of the display
        ========     ===================================================
        
        """
        oldlw = plt.rcParams['axes.linewidth']
        plt.rcParams['axes.linewidth'] = 2
        if fignum is not None: 
            fig=plt.figure(fignum, figsize=(4,4)); plt.clf()
            fig.add_axes((0.2,0.15,0.75,0.72))
        roi_plotting.make_sed(self,axes=axes, axis=axis, data_kwargs=data_kwargs, fit_kwargs=fit_kwargs)
        plt.rcParams['axes.linewidth'] = oldlw
        plt.title(self.name)

if __name__=='__main__':
    pass
