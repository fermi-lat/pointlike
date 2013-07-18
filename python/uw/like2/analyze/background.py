"""
background analysis
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/background.py,v 1.1 2013/07/09 02:05:43 burnett Exp $

"""
import os
import pandas as pd
import numpy as np
import pylab as plt

from . import roi_info

class Background(roi_info.ROIinfo):
    r"""This is an analysis of the sensitivity of point source analysis to the presence of 
    diffuse background. The ingredients are the PSF and the galactic and isotropic background fluxes.
    <br>
    The PSF is represented by the integral over its square, see the 
    <a href='../environment/index.html?skipDecoration'>environment summary</a>. 
    <br>We examine two ROIS: %(rois)s, representing low latitude (the galactic center) and high latitude, the north Galactic pole.
    In this analysis, the ROI size is %(roi_size)d degrees.
    """
    
    def introduction(self):
        r"""Introduction: summary of likelihood definition
        Consider $K$ energy bands, each with $N_k$ angular bins in an ROI, with $M$ component sources. 
        There are separate bands for front and back, the angular bin size is small compared with the PSF, 
        and there are four energy bands per decade, 14 total from 100 MeV to 316 GeV. 
        The bin sizes are small enough that the resolution is equivalent to unbinned likelihood, 
        the limit for zero-sized bins.

        Let:  
        <ul>
        <li>$ m_{kj}$: The number of measured events in energy band $k$ and angular bin $j$  </li> 
        <li>$ f_{ikj}$: the unnomalized exposure-weighted flux for source $i$ in band $k$, bin $j$ </li>    
        <li>$ F_{ik} = \sum\limits_j f_{ikj}$ </li>   
        <li>$ \alpha_i$: The prefactor, or multiplicative parameter for source $i$, to be determined. 
        (There are also of course shape parameters, not specified here.) Thus $\alpha_i f_{ijk}$ 
        is the predicted counts from source $i$ for band  $k$, bin $j$. </li>
        </ul>
        <p>
        Using Poisson statistics, the maximum likelihood solution for the $\alpha$ parameters is the 
        maximum with respect to $\alpha_i$ of the log likelihood,   

        $$ w(\alpha) = \sum\limits_k \left[ \sum\limits_j m_{kj} \ln(\sum\limits_i \alpha_i f_{ikj})  -\sum\limits_i\alpha_i F_{ik} \right] $$  

        or the set of $M$ equations  

        $$\  \frac{\partial w}{\partial \alpha_i} = 0 =\sum\limits_k \left[ \sum\limits_j m_{kj}  \frac{f_{ikj}}{\sum\limits_l\alpha_l f_{lkj}}  - F_{ik} \right]$$
        <h4>Caveats</h4><p>
        This procedure assumes perfect knowledge of the model coefficients $f_{ikj}$. 
        I want to deal with the special case of the IEM, or Galactic diffuse, for which there is an 
        additional uncertainty in the spectral values; that is, a systematic uncertainty for each energy band. 
        Since this involves an integral of the diffential flux over the range of energies, a factor of 1.8 in this case, 
        the dispersion may also be a factor.
        
        <h4>A simple example</h4>
        <p>A simple, but relevant, "toy" example is a single point source with a constant background, 
        for one energy band and conversion type. Center the ROI, with solid angle $\Delta \Omega$, 
        on the source and assume azimuthal symmetry and that the PSF is contained in the ROI. 
        Thus ignore the $j$ and $k$ indeces. For simplicity, fold the exposure, a constant, into 
        the $\alpha$ definitions and let $\alpha_0\equiv\alpha, \alpha_1\equiv\beta$. 
        Then, the log likelihood for a single photon at an angle $\theta$ with 
        respect to the point source is 
        
        $$ w(\alpha; \theta) = \ln(\alpha P(\theta) + \beta ) -(\alpha + \beta \Delta \Omega) $$

        Now $\alpha$ and $\beta$ are respectively the number of photons from the source and the background photon density.

        This form can be used to determine the expected resolution for the parameters. In the following,
        let $\bar{P}=\int \mathrm{d}\Omega P(\theta)^2$, the average value of the PSF, corresponding to the inverse of the
        solid angle extent. The inverse of this is a solid angle that can be interpreted as the effective size of the PSF, or
        "footprint".
        Assuming that the background dominates, the case of interest for weak sources,  $ \beta >> \alpha  \bar{P} $, 
        the expected variance matrix for $N=\beta \Delta \Omega$ photons is:

        $$ V = \frac{1}{\beta(\bar{P} \Delta \Omega -1)} \begin{pmatrix}
        \Delta\Omega & -1  \\
        -1 & \bar{P}  
        \end{pmatrix} $$
        
        <p> This simplifies further if the point source extent is small compared with the ROI, 
        $\bar{P} \Delta\Omega  >> 1$,  which is the case for all but the lowest energy. In particular, 
        $\sigma_\alpha^2=\beta / \bar{P}$, which is the number of background photons in the PSF footprint.
        """
        return None
    
    def setup(self, **kw):
        super(Background, self).setup(**kw)
        self.plotfolder='background'
        
        s = [[ x[1]['counts']['models'][modelnumber][1][:16] for x in self.df.iterrows()] for modelnumber in range(2)]
        self.bdf = [pd.DataFrame(y, index=self.df.index) for y in s]
        self.roi_size=5 # wired-in for now
        self.rois=[888,0]

    def get_background(self, roi):
        roiname='HP12_%04d' % roi
        return [t.ix[roiname] for t in self.bdf]
    
    def diffuse_flux(self, rois=None):
        """Diffuse flux
        Predicted counts for the low latitude and high latitude ROIs.
        """
        fig, ax = plt.subplots(1,1, figsize=(4,4), dpi=150, sharey=True)
        egev = np.array(self.energy)/1e3
        if rois is None: rois = self.rois

        for r in rois:
            gal, iso = self.get_background(r)
            ax.plot(egev, gal, '-d', label='gal %d'%r)
            ax.plot(egev, iso, '--o', label='iso %d'%r)
        plt.setp(ax, xscale='log', xlim=(0.1,300), xlabel='Energy (GeV)',
            yscale='log',  ylabel='Diffuse counts/ROI')
        ax.legend(prop=dict(size=10)); ax.grid()
        return fig

    def psf_background(self, rois=None, outfile='psf_bkg.csv',):
        """Background counts in PSF
        For galactic and isotropic backgrounds, and front and back, determine the number of counts 
        in the effective "footprint" of the PSF.
        
        """
        # get effective PSF size, and area from environment analysis
        epsf = pd.read_csv('plots/environment/psf.csv', index_col=0)
        psa = [np.pi*np.radians(epsf[ct].values[:14])**2 for ct in ['front', 'back']]
        
        # assume 5 degrees for solid angle (should get the value actually used)
        solid_angle = 2*np.pi*(1-np.cos(np.radians(self.roi_size)))
        
        fig, axx = plt.subplots(2,2, figsize=(7,7), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.,top=0.9)
        energy = self.energy 
        
        if rois is None: rois = self.rois
        flist = []
        blist = []
        
        for k,roi in enumerate(rois):
            gal, iso = self.get_background(roi)
            for ax, diffuse, what in zip(axx[k,:], (gal, iso), 'galactic isotropic'.split()):
                
                front, back = [ diffuse*psa[i]/solid_angle for i in range(2)]
                flist.append(front.values)
                blist.append(back.values)
                
                ax.plot(energy, front, '-o', label='front')
                ax.plot(energy, back, '-Dr', label='back')
                ax.legend(prop=dict(size=10)); ax.grid()
                ax.text(300,2e5, what, fontsize=12)
                plt.setp(ax, xscale='log', xlim=(100, 12000),
                    yscale='log', ylim=(1.0, 1e6), ylabel='counts in PSF' if what=='galactic' else '',)
            #plt.suptitle('ROI %04d' % roi)
            
        axx[0,1].set_xticklabels(['0.1', '1', '10'])
        fig.text(0.4, 0.05, 'Energy (GeV)')
        #make a DataFrame with the 8 plots
        self.psf_back_df = pd.DataFrame(dict(
            hgf=flist[0],  hgb=blist[0], 
            hif=flist[1],  hib=blist[1], 
            lgf=flist[2],  lgb=blist[2], 
            lif=flist[3],  lib=blist[3], 
             ), 
        index=(np.array(energy)*1e-3).round(3))
        self.psf_back_df.index.name='energy'
        self.psf_back_df.to_csv(os.path.join(self.plotfolder, outfile))
        print 'wrote file %s' % os.path.join(self.plotfolder, outfile)
    
        return fig
            
    def all_plots(self):
        self.runfigures([self.introduction, self.diffuse_flux, self.psf_background,])
