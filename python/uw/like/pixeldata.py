"""
Manage data and livetime information for an analysis


    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pixeldata.py,v 1.25 2012/07/16 16:44:08 lande Exp $

"""
version='$Revision: 1.25 $'.split()[1]
import os, math, pyfits, types, glob
import numpy as N; np = N
import pointlike, skymaps

from uw.utilities import path

class NsideMapper(object):
    """
    Manage the mapping between energy bands and pixel size, as parameterized
    by nside, the number of subdivisions of the base pixels in HEALPix 
    scheme.

    Roughly, the side of a (quadrilateral) pixel is 1/Nside radians, or
    60/Nside degrees.

    The default scheme is hardwired based on the pre-scaling of the PSF
    for Pass 6.  This is a power law with slope -0.8.  This prescaling
    gives, approximately, r68.  The default pixelization is so that, at
    100 MeV, 5 pixels fit within the sigma pre-scale.

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
    minnside = [0,0]

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
        nside = N.round(float(mns)/(1+mns*t)).astype(int)
        return N.maximum(nside,nsm.minnside[ct]).tolist()

class PixelDataException(Exception): pass

class PixelData(object):
    """
Manage data and livetime for spectral analysis

"""

    def __init__(self, aedict=dict(), **kwargs ):
        """
Create a new PixelData instance, managing data and livetime.
    aedict: a dictionary from the analysis environment
    kwargs: set one of the following parameters.
            binfile   = None,
            binsperdec= 4,
            conv_type = -1,
            emin      = 100,
            emax      = 1e6,
            event_class = 3,
            ft1files  = [],
            ft2files  = [],
            gti_mask  = None,
            ltcube    = None,
            recalcgti = False,
            thetacut  = 66.4,
            tstart    = 0,
            tstop     = 0,
            quiet     = False,
            roi_dir   = None,
            verbose   = False,
            use_weighted_livetime = False,
            zenithcut = 105,
            mc_src_id = -1,
            mc_energy = False


    """
        defaults =  dict(
            binfile   = None,   # if specified and exists, will be used for binned photon data
                                # if specified and does not exist, will be generated
            binsperdec= 4,      # bins per decade if generating data
            conv_type = -1,
            emin      = 100,
            emax      = 1e6,
            event_class = 3,
            ft1files  = [],
            ft2files  = [],
            gti_mask  = None,
            ltcube    = None,
            recalcgti = False,
            thetacut  = 66.4,
            tstart    = 0,
            tstop     = 0,
            quiet     = False,
            roi_dir   = None, exp_radius=180, # for selecting a cone when generating a LT cube
            verbose   = False,
            use_weighted_livetime = False,
            zenithcut = 105,
            mc_src_id = -1,
            mc_energy = False
        )
        self.__dict__.update(defaults)
        for key, value in aedict.items():
            if key in self.__dict__:
                self.__dict__[key] = value
        for key, value in kwargs.items():
            if key in self.__dict__:
                self.__dict__[key] = value
            else:
                raise PixelDataException, 'keyword %s not recognized' % key
        self._binner_set = False

        self._setup_files()
        # order of operations: GTI is needed for livetime; livetime is needed for PSF
        self.gti  =  self._GTI_setup()
        self.lt   =  self.get_livetime()
        self.weighted_lt = self.get_livetime(weighted=True,clobber=False) if self.use_weighted_livetime else None
        self.get_data()
        self.dmap.updateIrfs()

    def _setup_files(self):
        """
        use glob and expand to expand wildcards and convert to absolute paths (in backward-compatible manner)
        """
        if type(self.ft1files)==types.StringType:
            self.ft1files = glob.glob(path.expand(self.ft1files))
            self.ft1files.sort()
        if type(self.ft2files)== types.StringType:
            self.ft2files = glob.glob(path.expand(self.ft2files))
            self.ft2files.sort()
            
        # check explicit files
        for filelist in [self.ft1files, self.ft2files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise PixelDataException('PixelData setup: file name or path "%s" not found'%filelist[0])
        if type(self.binfile)==types.StringType:
            self.binfile = path.expand(self.binfile)
        if type(self.ltcube)==types.StringType:
            self.ltcube=path.expand(self.ltcube)


    def __str__(self):
        s ='Pixeldata: %.3f M events, %.3f Ms' % (self.dmap.photonCount()/1e6, self.gti.computeOntime()/1e6)
        for key in sorted(self.__dict__.keys()):
            if key[0]=='_': continue
            v = self.__dict__[key]
            if key[:2]=='ft' and len(v)>3:
                s += '\t%-20s: %s ... %s\n' %(key, v[0], v[-1])
            else:
                s += '\t%-20s: %s\n' %(key, v)
        return s

    def __repr__(self):
        return str(self)

    def reload_data(self,ft1files):
        self.ft1files = ft1files if type(ft1files)==type([]) or ft1files is None else [ft1files]
        for filelist in  [self.ft1files, self.ft2files] :
            if filelist is not None and len(filelist)>0 and not os.path.exists(filelist[0]):
                raise PixelDataException('PixelData setup: file name or path "%s" not found'%filelist[0])
        #self.data = self.get_data()
        #self.dmap = self.data.map()
        self.dmap = self.get_data()
        self.dmap.updateIrfs()

    def fill_empty_bands(self,bpd,bands):

        dummy = skymaps.SkyDir(0,0)

        for bin_center in (bands[:-1]*bands[1:])**0.5:
             ph_f = skymaps.Photon(dummy,bin_center,2.5e8,0)
             ph_b = skymaps.Photon(dummy,bin_center,2.5e8,1)
             bpd.addBand(ph_f)
             bpd.addBand(ph_b)

    def _Data_setup(self,bins):

        # check emin and bpd for consistency with CALDB
        c1 = N.abs(bins - 100).min() > 1
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

        pointlike.Data.set_class_level(self.event_class)
        pointlike.Data.set_zenith_angle_cut(self.zenithcut)
        pointlike.Data.set_theta_cut(self.thetacut)
        if 'mc_energy' in self.__dict__:
            pointlike.Data.set_use_mc_energy(self.mc_energy)
        pointlike.Data.set_Gti_mask(self.gti)
        if not self.quiet: print '.....set Data theta cut at %.1f deg'%(self.thetacut)

        if not self._binner_set:
            from pointlike import DoubleVector,IntVector
            f_nside = IntVector(NsideMapper.nside(bins,0))
            b_nside = IntVector(NsideMapper.nside(bins,1))
            self.binner = skymaps.PhotonBinner(DoubleVector(bins),f_nside,b_nside)
            pointlike.Data.setPhotonBinner(self.binner)
            self._binner_set = True

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

        if not self.recalcgti:
            if self.ltcube is not None and os.path.exists(self.ltcube):
                if not self.quiet: print('Using gti from %s'%(self.ltcube))
                return skymaps.Gti(self.ltcube)
            elif self.binfile is not None and os.path.exists(self.binfile):
                if not self.quiet: print('Using gti from %s'%(self.binfile))
                return skymaps.Gti(self.binfile)
        if len(self.ft1files)==0 or not os.path.exists(self.ft1files[0]):
            raise PixelDataException('Cannot process GTI: no ft1 files found, nor %s nor %s'%(self.ltcube,self.binfile))
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
            my_bins = 10**N.arange(N.log10(self.emin),N.log10(self.emax*1.01),1./self.binsperdec)

            self._Data_setup(my_bins)

            if not self.quiet: 
                print 'loading %d FT1 file(s) %s...%s' % (len(self.ft1files), self.ft1files[0], self.ft1files[-1])
                if self.event_class>-1: print 'selecting event_class %d' %self.event_class
            src_id = -1 if 'mc_src_id' not in self.__dict__ else self.mc_src_id
            data = pointlike.Data(self.ft1files,self.conv_type,self.tstart,self.tstop, src_id,'')

            self.dmap =data.map()
            self.fill_empty_bands(self.dmap, my_bins)     # fill any empty bins

            if self.verbose: print 'done'
            if self.binfile is not None:
                if not self.quiet: print '.....saving binfile %s for subsequent use' % self.binfile
                data.write(self.binfile)
        if self.binfile is not None:
            if not self.quiet: print '.....loading binfile %s ...' % self.binfile ,
            self.dmap = skymaps.BinnedPhotonData(self.binfile)
        if not self.quiet: print 'found %d bands, energies %.0f-%.0f MeV'\
                % (len(self.dmap), self.dmap[1].emin(), self.dmap[len(self.dmap)-1].emax())

        if self.verbose:
            self.dmap.info()
            print '---------------------'

    def get_livetime(self,pixelsize=1.0,weighted = False,clobber = True):

        gti = self.gti
        if self.ltcube is not None and os.path.exists(self.ltcube):
            try:
                lt = skymaps.LivetimeCube(self.ltcube,weighted=weighted)
                if not self.quiet: print 'loaded LivetimeCube %s ' % self.ltcube
                return lt
            except RuntimeError:
                if not self.quiet:
                    ext = 'WEIGHTED_EXPOSURE' if weighted else 'EXPOSURE'
                    print('no extension %s in file %s: will generate it from the ft2files'%(ext,self.ltcube))
        elif not self.quiet:
            print 'LivetimeCube file %s does not exist: will generate it from the ft2 files' % self.ltcube
        if self.roi_dir is None:
            # no roi specified: use full sky
            self.roi_dir = skymaps.SkyDir(0,0)
            self.exp_radius = 180
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
        if self.ltcube is not None:
            extension = 'WEIGHTED_EXPOSURE' if weighted else 'EXPOSURE'
            lt.write(self.ltcube,extension,clobber)
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
