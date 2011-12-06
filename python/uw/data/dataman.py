"""
Module implements classes and functions to specify data for use in pointlike analysis.

author(s): Matthew Kerr, Eric Wallace
"""

__version__ = '$Revision: 1.11 $'
#$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/data/dataman.py,v 1.11 2011/12/02 20:12:52 wallacee Exp $

import os
import collections
import glob
import warnings
from cPickle import dump,load

import numpy as np
import pyfits

import pointlike
import skymaps
from uw.utilities import keyword_options,fitstools
from uw.data import dssman

dataman_version=__version__.split()[1]
# TODO energy sanity check on DSS keywords
# -- scheme would be to cut all data that isn't commensurate
# -- or flag bad bins
# TODO event class check against Pass6 or Pass7
# TODO check to make sure FT2 and FT1 agree (compare Gti to FT2 times?)
# -- think this is too onerous if we allow FT1/FT2 lists
# -- eew: might also cause problems for some cases where the FT1 Gti is 
# --      subset of the FT2 Gti, which I think can be valid.
# TODO idea for binning -- just have it save both 4 and 8?
# TODO -- check exposure radius against ROI kwargs...
#      -- eew: perhaps better done in the ROI constructor?

def get_pass(ft1):
    """ Try to determine if the provided FT1 file is Pass6 or Pass7.
        If the algo. can't figure it out, Pass7 is assumed."""
    f = pyfits.open(ft1)
    h = f[1]._header
    v = h.get('PASS_VER')
    if (v is not None):
        if str(v).strip()[:2]=='P7': return 7
    if 'CTBCLASSLEVEL' in h.keys(): return 6
    return 7

def SimpleCut(vmin,vmax,vuni,colname):
    """ Specify a simple (inclusive or exclusive) cut on a single
        FT1 column, e.g. ZENITH_ANGLE < 100 or 100 < ENERGY < 100000.
        
        To specify an open-ended cut, use None for the upper or lower
        bound.
        
        This is a wrapper/factory for DSSSimpleRange."""

    return dssman.make_simple_dss(colname,vuni,vmin,vmax)

def get_default(colname,pass7=True):
    if colname == 'ZENITH_ANGLE':
        return SimpleCut(None,100,'deg','ZENITH_ANGLE')
    if colname == 'THETA':
        return SimpleCut(None,66.4,'deg','THETA')
    if colname == 'EVENT_CLASS':
        if pass7:
            d = dict(TYP='BIT_MASK(EVENT_CLASS,2)',UNI='DIMENSIONLESS',
                     VAL='1:1')
            return DSSBitMask(d)
        else:
            return SimpleCut(3,None,'dimensionless','EVENT_CLASS')

class NsideMapper(object):
    """
Manage the mapping between energy bands and pixel size, as parameterized
by nside, the number of subdivisions of the base pixels in HEALPix scheme.

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
        en  = np.asarray(en)/100.
        nsm = NsideMapper
        mns = nsm.maxnside[ct]
        t = nsm.norms[ct]*(en)**nsm.slopes[ct]*np.exp(-(en/nsm.cuts[ct])**2)
        nside = np.round(float(mns)/(1+mns*t)).astype(int)
        return np.maximum(nside,nsm.minnside[ct]).tolist()

class DataSpec(object):
    """ Class to specify the data for use in a spectral analysis.
        This includes both RAW data (e.g. FT1) and CUTS.
       
        Required cuts are on ZENITH, THETA, and EVENT_CLASS
        In order of precedence, these are taken from
        (1) the FT1 DSS keywords
        (2) user specification
        (3) defaults        
        An exception is raised if a set of FT1 files have differing
        DSS keywords. This is done to ensure consistency with gtmktime.

        Some other cuts are recognized by the pointlike machinery.
        (1) PULSE_PHASE (treated as a table, a la GTI) / TODO

        If a binned photon file (binfile) and livetime cube (ltcube)
        do not yet exist (see keywords), they will be created.  If the
        destination for these files isn't specified, files with sensible
        names will be created in the same directory as the base FT1 file.

        DSS keywords are saved both to the binfile and the ltcube,
        conventionally in the primary header (0).  On loading these
        data, they keywords will be compared with those saved at the
        time of their creation as a sanity check.
        """

    defaults = (
        ('ft1',None,'a file, list of files, or wildcard expression'),
        ('ft2','$FERMI/ft2.fits','a file, list of files, or wildcard expression'),
        ('binfile',None,'(a) destination for new binfile or (b) location of existing one'),
        ('ltcube',None,'(a) destination for new ltcube or (b) location of existing one'),
        ('binsperdec',8,'energy bins per decade; must be 8 or 4'),
        ('zenith_cut',None,'a SimpleCut wrapper giving zenith cuts'),
        ('theta_cut',None,'a SimpleCut wrapper giving theta cuts'),
        ('event_class_cut',None,'a SimpleCut wrapper giving event cuts'),
        ('gti_mask',None,'a GTI mask to apply to the data (intersection); note this can be used to set tstart/tstop for the data'),
        ('mc_src_id',-1,'select only photons from MC source ID; default is no selection'),
        ('mc_energy',False,'bin on MC_ENERGY instead of ENERGY'),
        ('clobber',False,'if True, will attempt to produce new binfile and ltcube and replace any existing ones'),
        ('quiet',True,'control verbosity, ever so coarsely'),
        ('use_weighted_livetime',True,'if True, calculate the weighted livetime for use in livetime-dependent corrections to the effective area'),
        ('livetime_buffer',10,'radius in degrees by which livetime cube cone is larger than ROI cone'),
        ('livetime_pixelsize',1,'pixel size to use for livetime calculation'),
        ('data_name', '', 'descriptive name for the data set'),
        ('legacy', False,  'relax DSS requirements for legacy files'),
        ('data_pass',7,'the generation (Pass6, Pass7,...) of the data'),
        # keyword controlling livetimecube pixel size? and buffer cone?
    )
    binner = None # static variable for PhotonBinner

    @keyword_options.decorate(defaults)
    def __init__(self,output=None,**kwargs):
        """ **NB -- if ft1 is None, binfile MUST be set to a real file
            **NB -- if ft2 is None, either ltcube must be set to a real file or $FERMI/ft2.fits must exist
        """
        keyword_options.process(self,kwargs)
        self._set_bins()
        self.dss = None # initialize

        self.ft1files = self._parse_filename(self.ft1)
        if self.ft1files is not None:
            # Register FT1 DSS keywords
            self._get_ft1_dss()
            self._make_cuts()
            self.gti = self._get_GTI()
            if not self._check_binfile():
                self._make_binfile()
        elif not self._check_binfile():
            raise ValueError('No FT1 files or valid binned data found.')

        self.ft2files = self._parse_filename(self.ft2)
        if self.ft2files is not None:
            if not self._check_ltcube():
                self._make_ltcube()
        elif not self._check_ltcube():
            raise ValueError('No FT2 files or valid livetime found.')

        # save version to allow custom processing for backwards compat.
        self.version = dataman_version
        if output is not None: self.dump(output)

    def __str__(self):
        """ Pretty print of cuts/data."""
        s = collections.deque()
        s.append('Bins per decade: {0}'.format(self.binsperdec))
        s.append('DSS keywords:\n{0}'.format(self.dss))
        def process_ft(files):
            if files is None: 
                s.append('\t\tNone')
                return
            if len(files) < 10:
                s.append('\n\t'.join(files))
            else:
                s.append('\n\t'.join(files[:5]))
                s.append('...')
                s.append('\n\t'.join(files[-5:]))
        s.append('FT1 files: ')
        process_ft(self.ft1files)
        s.append('FT2 files: ')
        process_ft(self.ft2files)
        s.append('Binned data:   {0}'.format(self.binfile))
        s.append('Livetime cube: {0}'.format(self.ltcube))
        return '\n'.join(s)
        
    def __getstate__(self):
        #If you're pickling, you shouldn't be clobbering
        self.clobber = False
        self.binner = None
        return self.__dict__
    def __setstate__(self,dict):
        """ Override default unpickle to perform a few sanity checks."""
        for t in self.defaults:
            self.__dict__[t[0]] = self.__dict__.pop(t[0],t[1])
        self.__dict__.update(dict)
        if not self._check_binfile():
            raise ValueError('Binned photon data not compatible!')
        if not self._check_ltcube():
            raise ValueError('Livetime cube not compatible!')
    def dump(self,filename):
        """ Output as a pickle."""
        dump(self,file(filename,'w'))
    def _set_bins(self):
        """ Check binning and set in Data class."""
        # TODO -- check if need Singleton
        if (self.binsperdec!=4) and (self.binsperdec!=8):
            raise ValueError('Only support 4 or 8 energy bins per decade.')
        self.bins = bins = np.logspace(1,6,5*self.binsperdec+1)
        f_nside = pointlike.IntVector(NsideMapper.nside(bins,0))
        b_nside = pointlike.IntVector(NsideMapper.nside(bins,1))
        DataSpec.binner = skymaps.PhotonBinner(pointlike.DoubleVector(bins),f_nside,b_nside)
        pointlike.Data.setPhotonBinner(DataSpec.binner)
        
    def _parse_filename(self,ft):
        """ Convert the argument into a list of filenames."""
        if ft is None: return None
        if isinstance(ft,str):
            # parse as a wildcard; should work for any string
            files = glob.glob(os.path.expandvars(ft))
        elif (isinstance(ft,list) or isinstance(ft,tuple)):
            files = list(ft)
        else:
            raise ValueError('Unable to parse FT input')
        # now check for validity of files
        if len(files) < 1:
            print 'FT file(s) "%s" not found: assume None to test for valid binfile' % ft
            return None
            raise DataError('FT file(s) "%s" not found' % ft)
        files = [os.path.expandvars(f) for f in files]
        files.sort()
        for f in files:
            if not os.path.exists(f):
                raise DataError('Unable to find file {0}'.format(f))
        return files
        
    def _make_cuts(self):
        """ Check existing FT1 cuts and apply default cuts.
            Order of precedence: user specified cut, existing FT1 cut, default cuts."""
        basic_cuts = [('ZENITH_ANGLE','zenith_cut'),
                      ('THETA','theta_cut'),
                      ('EVENT_CLASS','event_class_cut')]
        for col,mycut in basic_cuts:
            ft1_cut,index = self.dss.get_simple_dss(col)
            if self.__dict__[mycut] is not None:
                print('Working on {0}'.format(col))
                if ft1_cut is not None:
                    # apply more restrictive of two cuts
                    self.__dict__[mycut].intersection(ft1_cut)
                    self.__dict__[mycut]['index'] = index
                    self.dss[index] = self.__dict__[mycut] # replace in DSS
                else:
                    self.__dict__[mycut]['index'] = len(self.dss)+1
                    self.dss.append(self.__dict__[mycut])
            else:
                if ft1_cut is not None: self.__dict__[mycut] = ft1_cut
                else:
                    self.__dict__[mycut] = get_default(col)
                    self.__dict__[mycut]['index'] = len(self.dss)+1
                    self.dss.append(self.__dict__[mycut])
    def _get_ft1_dss(self):
        """ Get existing DSS keywords from FT1 files.
            If there is more than one FT1 file present, the protocol
            is that all DSS entries must agree."""
        self.data_pass = get_pass(self.ft1files[0])
        if len(self.ft1files) == 1:
            self.dss = dssman.DSSEntries(self.ft1files[0])
        else:
            all_dss = [dssman.DSSEntries(ft1) for ft1 in self.ft1files]
            if not np.all([all_dss[0] == x for x in all_dss]):
                raise ValueError('DSS keywords are inconsistent for FT1 files.')
            self.dss = all_dss[0]
            for ft1 in self.ft1files[1:]:
                # sanity check that all files have same pass
                if get_pass(ft1) != self.data_pass:
                    raise DataError('%s is not Pass %d'%(ft1,self.data_pass))
    def _make_default_name(self,prefix='bpd'):
        """ Come up with a default name for binfile or ltcube."""
        left,right = os.path.split(self.ft1files[0])
        return os.path.join(left,'{0}_{1}'.format(prefix,right))
    def _check_binfile(self):
        """ Verify binfile exists and is consistent with any existing data cuts."""
        if (self.clobber or self.binfile is None or
            (not os.path.exists(self.binfile))): return False
        dss = dssman.DSSEntries(self.binfile,header_key=0)
        gti = skymaps.Gti(self.binfile)
        if (dss is None):
            print("dss is None")
            return False
        if self.dss is None: # legacy case -- bpd with no FT1
            print('Legacy case -- accepting binfile on its own merits.')
            self.dss = dss
            self.gti = gti
            return True
        # check against FT1 info
        if (dss == self.dss) and (gti == self.gti):
            if (not self.quiet): print('Verified binfile {0}'.format(self.binfile))
            return True

        return False
    def _make_binfile(self):
        """ Generate the binned photon file."""
        self.binfile = self.binfile or self._make_default_name(prefix='bpd')
        self._Data_setup() # set up Data to use cuts
        def fill_empty_bands(bpd,bands):
            dummy = skymaps.SkyDir(0,0)
            for bin_center in (bands[:-1]*bands[1:])**0.5:
                 ph_f = pointlike.Photon(dummy,bin_center,2.5e8,0)
                 ph_b = pointlike.Photon(dummy,bin_center,2.5e8,1)
                 bpd.addBand(ph_f); bpd.addBand(ph_b)
        data = pointlike.Data(self.ft1files,-1,0,0, self.mc_src_id,'')
        dmap = data.map() # local reference to avoid segfaults
        fill_empty_bands(dmap, self.bins)
        dmap.write(self.binfile)
        self.dss.write(self.binfile,header_key=0)
    def _check_ltcube(self):
        """ Verify ltcube exists and is consistent with any existing data cuts."""
        if os.path.exists(self.ltcube) and self.legacy : 
            print('Accepting ltcube without dss checks since legacy specified')
            return True
        if self.clobber or (not os.path.exists(self.ltcube or '')): return False
        # check for presence of important history
        # eew -- primary header doesn't have info about extensions, so just
        #   open the file and use that handle to check header keys and
        #   extensions.
        # h = pyfits.getheader(self.ltcube,0) 
        lt = pyfits.open(self.ltcube)
        try:
            lt[0].header['RADIUS']; lt[0].header['PIXSIZE']
        except KeyError: pass #return False
        # check for weighted extension if we are using it
        if self.use_weighted_livetime:
            #try: h['WEIGHTED_EXPOSURE']
            #except KeyError: return False
            try: assert(len(lt)>3)
            except AssertionError: return False
        dss = dssman.DSSEntries(self.ltcube,header_key=0)
        gti = skymaps.Gti(self.ltcube)
        if (dss is None): return False
        # DSS *must* be set by this point, in binfile process
        # check against FT1 info
        if (dss == self.dss) and (gti == self.gti):
            if (not self.quiet): print('Verified ltcube {0}'.format(self.ltcube))
            return True
        return False
    def _make_ltcube(self):
        """ Generate the livetime cube."""
        #TODO -- unify weighted livetime
        import sys
        self.ltcube = self.ltcube or self._make_default_name(prefix='lt')
        roi_info = self.dss.roi_info()
        if (roi_info is None) or (roi_info[2] > 90):
            # full sky ltcube
            roi_dir = skymaps.SkyDir(0,0)
            exp_radius = 180
        else:
            roi_dir = skymaps.SkyDir(roi_info[0],roi_info[1])
            exp_radius = roi_info[2] + self.livetime_buffer
        zenithcut = self.zenith_cut.get_bounds()[1]
        if not self.quiet:
            print('Constructing livetime cube about RA,decl = ({0:0.3f},{1:0.3f}) with a radius of {2:0.3f} deg.'.format(roi_dir.ra(),roi_dir.dec(),exp_radius))
        for i in xrange(1+self.use_weighted_livetime):
            print('on iteration {0}'.format(i))
            sys.stdout.flush()
            lt = skymaps.LivetimeCube(
                cone_angle =exp_radius,
                dir        =roi_dir,
                zcut       =np.cos(np.radians(zenithcut)),
                pixelsize  =self.livetime_pixelsize,
                quiet      =self.quiet,
                weighted   =i)

            for hf in self.ft2files:
                if not self.quiet: print('loading FT2 file {0}'.format(hf))
                lt_gti = skymaps.Gti(hf,'SC_DATA')
                if not ((lt_gti.maxValue() < self.gti.minValue()) or
                        (lt_gti.minValue() > self.gti.maxValue())):
                   lt.load(hf,self.gti)

            # write out ltcube
            extension = 'WEIGHTED_EXPOSURE' if i else 'EXPOSURE'
            lt.write(self.ltcube,extension,not bool(i))
        self.dss.write(self.ltcube,header_key=0)
        # write some info to livetime file
        f = pyfits.open(self.ltcube)
        f[0]._header.update('RADIUS',exp_radius)
        f[0]._header.update('PIXSIZE',self.livetime_pixelsize)
        f.writeto(self.ltcube,clobber=True)
        f.close()
    def _get_GTI(self):
        """ Apply GTI cuts and get resulting merged GTI."""
        # get GTI for a set of FT1 files
        gti = skymaps.Gti(self.ft1files[0])
        for ef in self.ft1files[1:]:
            gti.combine(skymaps.Gti(ef))
        if self.gti_mask is not None:
            before = gti.computeOntime()
            gti.intersection(self.gti_mask)
            if not self.quiet:
                print('applied gti mask, before, after times: {0:.1f}, {1:.1f}'
                      .format(before, gti.computeOntime()))
        return gti
    def _Data_setup(self):
        """ Set static variables in Data, GTI, etc. """
        zenithcut = self.zenith_cut.get_bounds()[1]
        thetacut = self.theta_cut.get_bounds()[1]
        event_class = self.event_class_cut.get_bounds()[0]
        pointlike.Data.set_class_level(event_class)
        pointlike.Data.set_zenith_angle_cut(zenithcut)
        pointlike.Data.set_theta_cut(thetacut)
        pointlike.Data.set_use_mc_energy(self.mc_energy)
        pointlike.Data.set_Gti_mask(self.gti)

    def check_consistency(self,other):
        """Check compatibility to combine with another DataSpec
        
        In order to be considered compatible for combination the two DataSpec
        instances must have:
            consistent DSS entries
            consistent binning (bins/decade, and MC specs)
            consistent livetime specifications (pixelsize, radius, weighting)
        If incompatible, return an instance of DataError which can then be
        raised if desired, else return None.
        """
        if self.dss!=other.dss:
            return DataError("DataSpec instances have inconsistent DSS keywords")
        if self.binsperdec!=other.binsperdec:
            return DataError("DataSpec instances have inconsistent binning")
        if self.livetime_pixelsize!=other.livetime_pixelsize:
            return DataError("DataSpec instances have inconsistent livetime pixelsize")
        if self.livetime_buffer!=other.livetime_buffer:
            return DataError("DataSpec instances have inconsistent livetime buffer")
        if self.use_weighted_livetime!=other.use_weighted_livetime:
            return DataError("DataSpec instances have inconsistent use_weighted_livetime")
        if self.mc_energy!=other.mc_energy:
            return DataError("DataSpec instances have inconsistent mc_energy")
        if self.mc_src_id!=other.mc_src_id:
            return DataError("DataSpec instances have inconsistent mc_src_id")
        return

    def add(self,others,output,binfile,ltcube):
        """Combine this DataSpec instance with another and return the result

        The two instances must have consistent definitions (DSS, binning,
        and livetime), and must have non-overlapping Gtis. The binfiles and
        ltcubes will be combined and written to the provided destinations;
        perhaps in the future, I will come up with sensible defaults.
        """
        if not hasattr(others,'__iter__'):
            others = [others]
        for other in others:
            exc = self.check_consistency(other)
            if exc is not None:
                raise(exc)
        gti = skymaps.Gti(self.gti)
        ft1 = self.ft1files
        ft2 = self.ft2files
        bpd = skymaps.BinnedPhotonData(self.binfile)
        for other in others:
            gti.intersection(other.gti)
            ft1 += other.ft1files
            ft2 += other.ft2files
            bpd.add(skymaps.BinnedPhotonData(other.binfile))
        if gti.computeOntime()>0:
            raise DataError("DataSpec instances have overlapping GTIs")
        ft2 = sorted(list(set(ft2)))
        bpd.write(binfile)
        fitstools.merge_lt([self.ltcube]+[other.ltcube for other in others],
                            outfile=ltcube,weighted=self.use_weighted_livetime)
        dssman.DSSEntries(self.binfile,header_key=0).write(binfile,header_key=0)
        dssman.DSSEntries(self.ltcube,header_key=0).write(ltcube,header_key=0)
        #TODO: move the copying of DSS entries into the merge_bpd and merge_lt functions
        gti_mask = skymaps.Gti(self.gti_mask)
        for other in others:
            gti_mask.combine(other.gti_mask)
        kwargs = dict(ft1=ft1,
                      ft2=ft2,
                      binfile=binfile,
                      ltcube=ltcube,
                      gti_mask = gti_mask,
                      binsperdec=self.binsperdec,
                      mc_src_id=self.mc_src_id,
                      mc_energy=self.mc_energy,
                      use_weighted_livetime = self.use_weighted_livetime,
                      livetime_buffer = self.livetime_buffer,
                      livetime_pixelsize = self.livetime_pixelsize,
                      clobber = False)
        return DataSpec(output,**kwargs) 
    #TODO -- replace/modify ExposureManager
    
    def __call__(self):
        """Return a DataManager created from this DataSpec"""
        return DataManager(self)


class DataManager(object):
    """A wrapper class holding references to data objects

    From the filenames specified in a DataSpec, loads and holds references
    to a BinnedPhotonData (stored in the bpd member) and a LivetimeCube (or
    two, if there's a weighted one; stored as lt and weighted_lt).
    """
    def __init__(self,ds):
        """Initialize a DataManager instance.
        
            ds: a DataSpec instance
        """
        self.dataspec = ds
        self.bpd = skymaps.BinnedPhotonData(ds.binfile)
        self.lt = skymaps.LivetimeCube(ds.ltcube,weighted=False)
        if ds.use_weighted_livetime:
            self.weighted_lt = skymaps.LivetimeCube(ds.ltcube,weighted=True)
        else:
            self.weighted_lt = None
        self.gti = self.lt.gti() #Just to provide a reference.

    @property
    def dmap(self):
        """Alias for backward compatibility"""
        warnings.warn(DeprecationWarning('DataManager.bpd is the preferred name for the BinnedPhotonData object'))
        return self.bpd

class DataSet(object):
    """A helper class to manage DataSpecs

    The basic function of this class is to retrieve pickled DataSpecs from a 
    standard directory structure. If a dataset contains multiple DataSpecs
    (e.g. for individual months), this class will store and provide access to
    the list. It can also be used to manage one or more pre-loaded DataSpecs.
    """
    defaults = (('pickle',None,
                 '''A filename, list of filenames, or wildcard expression
                    indicating a set of pickled DataSpecs''')
               ,('dataspec',None,'One or more DataSpec instances')
               ,('data_dir',os.path.join('$FERMI','data'),
                 'Path to main data directory')
               ,('clobber',False,'If True, overwrite existing DataSet pickle')
               )
    @keyword_options.decorate(defaults)
    def __init__(self,name=None,**kwargs):
        """Instantiate a DataSet

        If name is not None, it will be converted to a filename as
        os.path.expandvars(os.path.join('$FERMI','data',name.replace(' ','_')+'.pickle')).
        If this file exists, it will be taken to be a pickled DataSpec or
        list of DataSpec pickles. If it does not exist, it will be created with
        the appropriate format for reuse with this class. If the $FERMI
        environment variable is not set, files will be looked up and saved in
        the current working directory.
        Passing a value for name is the preferred use. If name is None, the
        user must provide a value for either pickle or dataspec but NOT both.
        In this case, the instance will manage the DataSpecs specified by
        whichever of those kwargs is used, but will not be saved for reuse.
        A single file passed via the pickle kwarg should point to either a
        pickled DataSpec, or a pickled list of DataSpec pickle files. A list of
        filenames, or a wildcard that matches a list, should each point to
        a pickled DataSpec.
        """

        keyword_options.process(self,kwargs)
        self.name = name
        self.filename = self._parse_name(name)
        if name is None or not os.path.exists(self.filename):
            if self.pickle is None and self.dataspec is None:
                raise DataError("Must specify either pickle files or DataSpecs")
        if os.path.exists(self.filename or '') and not self.clobber:
            self.dataspec = self._load_files(self.filename)
        else:
            if self.pickle is not None:
                if self.dataspec is not None:
                    raise DataError("Must specify dataspec OR pickle, not both.")
                self.dataspec = self._load_files(self.pickle)
            if self.filename is not None:
                dump(self.dataspec,open(self.filename,'w'))

    def _parse_name(self,name):
        """Return the filename corresponding to a DataSet name"""
        if name is None: return name
        data_dir = os.path.expandvars(self.data_dir)
        if data_dir=="$FERMI/data":
            #$FERMI undefined
            print("$FERMI environment variable not set - using current directory")
            data_dir = os.getcwd()
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        basename = name.replace(' ','_')
        #in case we're parsing the name of a pickle file
        if not basename.endswith('.pickle'):
            basename = '.'.join([basename,'pickle'])
        return os.path.join(data_dir,basename)

    def _load_files(self,pfile):
        """Load a pickle file and return the resulting DataSpec(s)"""
        if hasattr(pfile,'__iter__'):
            return [self._load_files(pf) for pf in pfile]
        pfile = os.path.expandvars(pfile)
        if not os.path.isabs(pfile):
            pfile = self._parse_name(pfile)
        gfiles = sorted(glob.glob(pfile))
        if len(gfiles)==0:
            raise DataError("No matches found for wildcard {0}".format(pfile))
        elif len(gfiles)>1:
            return [self._load_files(gf) for gf in gfiles]
        else:
            pfile = gfiles[0]
        pdat = load(open(pfile))
        if isinstance(pdat,DataSpec):
            return pdat
        elif hasattr(pdat,'__iter__'):
            if isinstance(pdat[0],DataSpec):
                return pdat
            return self._load_files(pdat)
        else:
            #temporary kludge for my local version
            try:
                from users.kerrm import tools
                if isinstance(pdat,tools.dataman.DataSpec):
                    return pdat
            except ImportError:
                pass
            raise DataError("Invalid DataSet pickle: {0}".format(pfile))

    def __getitem__(self,i):
        if hasattr(self.dataspec,'__iter__'):
            return self.dataspec[i]
        else:
            raise DataError('DataSet only contains one DataSpec')
    def __getattr__(self,x):
        return getattr(self.dataspec,x)
    def __call__(self):
        if not hasattr(self.dataspec,'__iter__'):
            return self.dataspec()
        else:
            raise TypeError('"DataSpec" object is not callable.')

class DataError(Exception):
    """Exception class for data management errors"""
    pass
