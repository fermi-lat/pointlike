"""
Manage data and livetime information for an analysis


    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pixeldata.py,v 1.21 2009/12/08 01:55:19 wallacee Exp $

"""
version='$Revision: 1.21 $'.split()[1]
import os
import math
import skymaps
import pointlike


class PixelData(object):
    """
Create a new PixelData instance, managing data and livetime.  

    analysis_environment: an instance of AnalysisEnvironment correctly configured with
                          the location of files needed for spectral analysis (see its
                          docstring for more information.)

Optional keyword arguments:

   see docstring for SpectralAnalysis
"""

    def __init__(self, spectral_analysis ):
        """
Create a new PixelData instance, managing data and livetime.  

    analysis_environment: an instance of AnalysisEnvironment correctly configured with
                          the location of files needed for spectral analysis (see its
                          docstring for more information.)

Optional keyword arguments:

   see docstring for SpectralAnalysis
"""

        self.__dict__.update(spectral_analysis.__dict__)

        from numpy import arange,log10
        self.my_bins = 10**arange(log10(self.emin),log10(self.emax*1.01),1./self.binsperdec)
        self.binner_set = False
        
        # check explicit files
        for filelist in [self.ft1files, self.ft2files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise Exception('PixelData setup: file name or path "%s" not found'%filelist[0])

        self.lt   =  self.get_livetime()
        self.PSF_setup()
        self.data = self.get_data()
        self.dmap = self.data.map()
        self.dmap.updateIrfs()


    def reload_data(self,ft1files):
        self.ft1files = ft1files if type(ft1files)==type([]) or ft1files is None else [ft1files]
        for filelist in  [self.ft1files, self.ft2files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise Exception('PixelData setup: file name or path "%s" not found'%filelist[0])
        self.data = self.get_data()
        self.dmap = self.data.map()
        self.dmap.updateIrfs()

    def fill_empty_bands(self,bpd):
      
        from pointlike import Photon
        from skymaps import SkyDir
        dummy = SkyDir(0,0)
        bands = self.my_bins
        
        for bin_center in (bands[:-1]*bands[1:])**0.5:
             ph_f = Photon(dummy,bin_center,2.5e8,0)
             ph_b = Photon(dummy,bin_center,2.5e8,1)
             bpd.addPhoton(ph_f,0)
             bpd.addPhoton(ph_b,0)
             
    def Data_setup(self):

        from numpy import arccos,pi
        pointlike.Data.set_class_level(self.event_class)
        pointlike.Data.set_zenith_angle_cut(self.zenithcut)
        pointlike.Data.set_theta_cut(self.thetacut)
        pointlike.Data.set_use_mc_energy(self.mc_energy)
        if not self.quiet: print '.....set Data theta cut at %.1f deg'%(self.thetacut)

        if not self.binner_set:
            from pointlike import DoubleVector
            self.binner = skymaps.PhotonBinner(DoubleVector(self.my_bins))
            pointlike.Data.setPhotonBinner(self.binner)
            self.binner_set = True

    def PSF_setup(self):
   
        # modify the psf parameters in the band objects
        # note there is a small discrepancy since the PSF average does not
        # currently know about the theta cut
        ip = skymaps.IParams
        #if self.ltcube is not None and self.roi_dir is not None:
           #ip.set_livetimefile(self.ltcube)
           #ip.set_skydir(self.roi_dir)
        ip.set_CALDB(self.CALDB)
        ip.init('_'.join(self.irf.split('_')[:-1]),self.irf.split('_')[-1])

    def get_data(self):
        
        #if no binned object present, create; apply cuts
        if self.binfile is None or not os.path.exists(self.binfile):

            self.Data_setup()

            if self.verbose: print 'loading file(s) %s' % self.ft1files
            data = pointlike.Data(self.ft1files,self.conv_type,self.tstart,self.tstop,self.mc_src_id,'')

            self.fill_empty_bands(data.map())     # fill any empty bins           

            if self.verbose: print 'done'
            if self.binfile is not None:
                if not self.quiet: print '.....saving binfile %s for subsequent use' % self.binfile
                data.write(self.binfile)
        
        else:
            data = pointlike.Data(self.binfile)
            if not self.quiet: print '.....loaded binfile %s ' % self.binfile
        
        if self.verbose:
            data.map().info()
            print '---------------------'

        return data

    def get_livetime(self,   pixelsize=1.0):
    
        gti = skymaps.Gti(self.ft1files[0])
        for ef in self.ft1files[1:]:
            gti.combine(skymaps.Gti(ef))
        tmax = self.tstop if self.tstop > 0 else gti.maxValue()

        gti = self.gti = gti.applyTimeRangeCut(self.tstart,tmax) #save gti for later use

        if self.ltcube is None or not os.path.exists(self.ltcube):
            if self.roi_dir is None:
                # no roi specified: use full sky
                self.roi_dir = skymaps.SkyDir(0,0)
                self.exp_radius = 180

            lt = skymaps.LivetimeCube(
                cone_angle =self.exp_radius,
                dir        =self.roi_dir,
                zcut       =math.cos(math.radians(self.zenithcut)),
                pixelsize  =pixelsize,
                quiet      =self.quiet )

            for hf in self.ft2files:
                lt_gti = skymaps.Gti(hf,'SC_data')
                if not ((lt_gti.maxValue() < self.gti.minValue()) or 
                        (lt_gti.minValue() > self.gti.maxValue())):
                   lt.load(hf,gti)

            # write out ltcube if requested
            if self.ltcube is not None: lt.write(self.ltcube)
        else:
            # ltcube exists: just use it! (need to bullet-proof this)
            lt = skymaps.LivetimeCube(self.ltcube)
            print 'loaded LivetimeCube %s ' % self.ltcube
        return lt


#--------------------------------------------------------------------
# test for this application
#--------------------------------------------------------------------

def setup(bins=8):
    data_path = r'f:\glast\data\flight'
    months = ['aug', 'sep', 'oct']
    event_files  = [os.path.join(data_path, '%s2008-ft1.fits'%m) for m in months ]
    history_files= [os.path.join(data_path, '%s2008-ft2.fits'%m) for m in months ]
    return PixelData( event_files, history_files,  
        livetimefile='../data/aug-oct_livetime.fits', 
        datafile='../data/aug-oct2008_%d.fits'%bins, 
        binsperdec=bins)
if __name__=='__main__':
    pass
