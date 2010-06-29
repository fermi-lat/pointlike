"""
Manage data and livetime information for an analysis


    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pixeldata.py,v 1.12 2010/06/28 20:32:58 kerrm Exp $

"""
version='$Revision: 1.12 $'.split()[1]
import os
import math
import skymaps
import pointlike
import pyfits
import numpy as N

class NsideMapper(object):
    """
Manage the mapping between energy bands and pixel size, as parameterized
by nside, the number of subdivisions of the base pixels in HEALPix scheme.

Roughly, the side of a (quadrilateral) pixel is 1/Nside radians, or
60/Nside degrees.

The default scheme is hardwired based on the pre-scaling of the PSF
for Pass 6.  This is a power law with slope -0.8.  This prescaling
gives, approximately, r68.

The default pixelization is chosen so that 5 pixels fit within r68.

This pixelization continues until about 1 GeV, when nside goes rapidly
to 8192, the maximum value for 32-bit architecture, with a pixel size
of 0.007 deg.  We do this because with pixels appropriate to the PSF
the mean pixel occupation becomes 1 at about 1 GeV, so there is nothing
to be gained by binning.  Using such small pixels means the data are
essentially unbinned in position for E > a few GeV.
    """

    norms    = [0.0116,0.0192] # pix size at 100 MeV in radians
    slopes   = [-0.8,-0.8]     # slope for pix size with energy
    cuts     = [20.,20.]   # "cutoff" energy, 2 GeV = 20 in E_100 units
    maxnside = [8192,8192]

    @staticmethod
    def nside(en,ct=0):
        """Return nside for provide energy and conversion type.
           en -- energy in MeV
           ct -- conversion type (0/1)
        """
        en  = N.asarray(en)/100.
        nsm = NsideMapper
        mns = nsm.maxnside[ct]
        t = nsm.norms[ct]*(en)**nsm.slopes[ct]*N.exp(-(en/nsm.cuts[ct])**2)
        return N.round(float(mns)/(1+mns*t)).astype(int)


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

        self.__dict__.update(spectral_analysis)

        from numpy import arange,log10
        self.my_bins = 10**arange(log10(self.emin),log10(self.emax*1.01),1./self.binsperdec)
        self.binner_set = False

        # check explicit files
        for filelist in [self.ft1files, self.ft2files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise Exception('PixelData setup: file name or path "%s" not found'%filelist[0])

        # order of operations: GTI is needed for livetime; livetime is needed for PSF
        self.gti  =  self._GTI_setup()
        self.lt   =  self.get_livetime()
        self.weighted_lt = self.get_livetime(weighted=True) if self.use_weighted_livetime else None
        self.data = self.get_data()
        self.dmap = self.data.map()
        self.dmap.updateIrfs()

    def __str__(self):
        return 'Pixeldata: %.3f M events, %.3f Ms' % (self.dmap.photonCount()/1e6, self.gti.computeOntime()/1e6)

    def __repr__(self):
        return str(self)

    def reload_data(self,ft1files):
        self.ft1files = ft1files if type(ft1files)==type([]) or ft1files is None else [ft1files]
        for filelist in  [self.ft1files, self.ft2files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise Exception('PixelData setup: file name or path "%s" not found'%filelist[0])
        self.data = self.get_data()
        self.dmap = self.data.map()
        self.dmap.updateIrfs()

    def fill_empty_bands(self,bpd):

        dummy = skymaps.SkyDir(0,0)
        bands = self.my_bins

        for bin_center in (bands[:-1]*bands[1:])**0.5:
             ph_f = pointlike.Photon(dummy,bin_center,2.5e8,0)
             ph_b = pointlike.Photon(dummy,bin_center,2.5e8,1)
             bpd.addBand(ph_f)
             bpd.addBand(ph_b)

    def _Data_setup(self):

        # check emin and bpd for consistency with CALDB
        c1 = N.abs(self.my_bins - 100).min() > 1
        c2 = (self.binsperdec % 4) > 0
        if c1 or c2:
            print """
            ###################WARNING!!!##########################
            It is STRONGLY recommended that you use a binning comm-
            ensurate with CALDB binning.  This can be achieved by 
            making sure 100 MeV appears in the bin edges (an easy 
            way is to set emin=100) and that binsperdec is a multi-
            ple of 4.  You can choose a subset of these bins in the
            likelihood fitting if you want to narrow the energy
            bounds.
            #######################################################
            """

        from numpy import arccos,pi
        pointlike.Data.set_class_level(self.event_class)
        pointlike.Data.set_zenith_angle_cut(self.zenithcut)
        pointlike.Data.set_theta_cut(self.thetacut)
        pointlike.Data.set_use_mc_energy(self.mc_energy)
        pointlike.Data.set_Gti_mask(self.gti)
        if not self.quiet: print '.....set Data theta cut at %.1f deg'%(self.thetacut)

        if not self.binner_set:
            from pointlike import DoubleVector,IntVector
            f_nside = IntVector(NsideMapper.nside(self.my_bins,0))
            b_nside = IntVector(NsideMapper.nside(self.my_bins,1))
            self.binner = skymaps.PhotonBinner(DoubleVector(self.my_bins),f_nside,b_nside)
            pointlike.Data.setPhotonBinner(self.binner)
            self.binner_set = True


    def _GTI_setup(self):
        """Create the GTI object that will be used to filter photon events and
           to calculate the livetime.

           If the user has provided ltcube or binfile arguments, the saved GTI
           will be used.  The ltcube GTI takes precedence over the binfile GTI;
           if both are provided, it is the users responsibility to ensure that
           the two sets of GTIs are consistent.

           Generally, one will have run gtmktime on the FT1 file(s) provided
           to filter out large excursions in rocking angle, and/or make a
           a cut for when a (small) ROI is too close to the horizon.

           gti_mask is an optional kwarg to provide an additional GTI object that
           can, e.g., be used to filter out GRBs.  An intersection with the GTIs
           from the FT1 files and this gti_mask is made.
        """
        if not self.quiet: print('applying GTI')

        if self.ltcube is not None and os.path.exists(self.ltcube):
            if self.verbose: print('Using gti from %s'%(self.ltcube))
            return skymaps.Gti(self.ltcube)
        elif self.binfile is not None and os.path.exists(self.binfile):
            if self.verbose: print('Using gti from %s'%(self.binfile))
            return skymaps.Gti(self.binfile)

        gti = skymaps.Gti(self.ft1files[0])

        # take the union of the GTI in each FT1 file
        for ef in self.ft1files[1:]:
            gti.combine(skymaps.Gti(ef))
        tmax = self.tstop if self.tstop > 0 else gti.maxValue()

        gti = gti.applyTimeRangeCut(self.tstart,tmax) #save gti for later use

        if 'gti_mask' in self.__dict__ and self.gti_mask is not None:
            before = gti.computeOntime()
            gti.intersection(self.gti_mask)
            if not self.quiet: print 'applied gti mask, before, after times: %.1f, %.1f' % (before, gti.computeOntime())

        return gti

    def get_data(self):

        #if no binned object present, create; apply cuts
        if self.binfile is None or not os.path.exists(self.binfile):

            self._Data_setup()

            if not self.quiet: print 'loading file(s) %s' % self.ft1files
            data = pointlike.Data(self.ft1files,self.conv_type,self.tstart,self.tstop,self.mc_src_id,'')

            self.fill_empty_bands(data.map())     # fill any empty bins

            if self.verbose: print 'done'
            if self.binfile is not None:
                if not self.quiet: print '.....saving binfile %s for subsequent use' % self.binfile
                data.write(self.binfile)

                """
                # a start on adding keywords -- not yet finished, but need to merge in CVS
                # now, add the appropriate entries to the FITS header
                f = pyfits.open(self.binfile)
                h = f[1].header
                pos = len(h.keys)
                entries = []
                entries += ['EMIN', self.emin, 'Minimum energy in MeV']
                entries += ['EMAX', self.emax, 'Maximum energy in MeV']
                entries += ['BINS_PER_DECADE', self.binsperdec, 'Number of (logarithmic) energy bins per decade']
                entries += ['TMIN', self.tmin, 'Exclude all data before this MET']
                entries += ['TMAX', self.tmax, 'Exclude all data after this MET']
                entries += ['EVENT_CLASS', self.event_class, 'Exclude all data with class level < this value']
                entries += ['CONVERSION_TYPE', self.event_class, '0 = Front only, 1 = Back only, -1 = Front + Back']

                for entry in entries:
                    k,v,c = entry
                    h.update(k,str(v),c)

                f.writeto(self.binfile,clobber=True)
                f.close()
                """
        else:
            data = pointlike.Data(self.binfile)
            if not self.quiet: print '.....loaded binfile %s ' % self.binfile

        if self.verbose:
            data.map().info()
            print '---------------------'

        return data

    def test_pixelfile_consistency(self):
        """Check the keywords in a binned photon data header for consistency with the analysis environment."""
        pass

    def get_livetime(self,   pixelsize=1.0,weighted = False):

        gti = self.gti
        ltcube = self.weighted_ltcube if weighted else self.ltcube
        if ltcube is None or not os.path.exists(ltcube):
            if self.roi_dir is None:
                # no roi specified: use full sky
                self.roi_dir = skymaps.SkyDir(0,0)
                self.exp_radius = 180
            if not self.quiet:
                print 'LivetimeCube file %s does not exist: will generate it from the ft2 files' % self.ltcube
            lt = skymaps.LivetimeCube(
                cone_angle =self.exp_radius,
                dir        =self.roi_dir,
                zcut       =math.cos(math.radians(self.zenithcut)),
                pixelsize  =pixelsize,
                quiet      =self.quiet,
                weighted   =weighted)

            for hf in self.ft2files:
                if not self.quiet: print 'loading FT2 file %s' %hf ,
                lt_gti = skymaps.Gti(hf,'SC_DATA')
                if not ((lt_gti.maxValue() < self.gti.minValue()) or
                        (lt_gti.minValue() > self.gti.maxValue())):
                   lt.load(hf,gti)

            # write out ltcube if requested
            if self.ltcube is not None: lt.write(ltcube)
        else:
            # ltcube exists: just use it! (need to bullet-proof this)
            lt = skymaps.LivetimeCube(ltcube)
            if not self.quiet: print 'loaded LivetimeCube %s ' % ltcube
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
