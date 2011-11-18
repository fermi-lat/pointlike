"""
Module implements classes and functions to specify data for use in pointlike analysis.

author(s): Matthew Kerr, Eric Wallace
"""

__version__ = '$Revision: 1.1 $'
#$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/data/dataman.py,v 1.1 2011/11/10 23:54:34 wallacee Exp $

import os
from collections import deque
from cPickle import dump,load
from glob import glob

import numpy as np
import pyfits

import pointlike
import skymaps
from uw.utilities import keyword_options,fitstools
from uw.like.pixeldata import NsideMapper
import dssman

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

def SimpleCut(vmin,vmax,vuni,colname):
    """ Specify a simple (inclusive or exclusive) cut on a single
        FT1 column, e.g. ZENITH_ANGLE < 100 or 100 < ENERGY < 100000.
        
        To specify an open-ended cut, use None for the upper or lower
        bound.
        
        This is a wrapper/factory for DSSSimpleRange."""

    return dssman.make_simple_dss(colname,vuni,vmin,vmax)

def get_default(colname):
    if colname == 'ZENITH_ANGLE':
        return SimpleCut(None,100,'deg','ZENITH_ANGLE')
    if colname == 'THETA':
        return SimpleCut(None,66.4,'deg','THETA')
    if colname == 'EVENT_CLASS':
        return SimpleCut(3,None,'dimensionless','EVENT_CLASS')

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
        ('weighted_livetime',True,'if True, calculate the weighted livetime for use in livetime-dependent corrections to the effective area'),
        ('livetime_buffer',10,'radius in degrees by which livetime cube cone is larger than ROI cone'),
        ('livetime_pixelsize',1,'pixel size to use for livetime calculation')
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
        s = deque()
        s.append('Bins per decade: {0}'.format(self.binsperdec))
        s.append('DSS keywords:\n{0}'.format(self.dss))
        def process_ft(files):
            if len(files) < 10:
                s.append('\n'.join(files))
            else:
                s.append('\n'.join(files[:5]))
                s.append('...')
                s.append('\n'.join(files[-5:]))
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
            files = glob(os.path.expandvars(ft))
        elif (isinstance(ft,list) or isinstance(ft,tuple)):
            files = list(ft)
        else:
            raise ValueError('Unable to parse FT input')
        # now check for validity of files
        if len(files) < 1:
            raise IOError('Unable to find indicated FT input')
        files = [os.path.expandvars(f) for f in files]
        files.sort()
        for f in files:
            if not os.path.exists(f):
                raise IOError('Unable to find file {0}'.format(f))
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
        if len(self.ft1files) == 1:
            self.dss = dssman.DSSEntries(self.ft1files[0])
        else:
            all_dss = [dssman.DSSEntries(ft1) for ft1 in self.ft1files]
            if not np.all([all_dss[0] == x for x in all_dss]):
                raise ValueError('DSS keywords are inconsistent for FT1 files.')
            self.dss = all_dss[0]
    def _make_default_name(self,prefix='bpd'):
        """ Come up with a default name for binfile or ltcube."""
        left,right = os.path.split(self.ft1files[0])
        return os.path.join(left,'{0}_{1}'.format(prefix,right))
    def _check_binfile(self):
        """ Verify binfile exists and is consistent with any existing data cuts."""
        #NOTE: this if statement causes a DataSpec created with clobber==True
        #to be unloadable from the pickle file! For now, fixed by setting
        #clobber=False in __set_state__.
        if self.clobber or (not os.path.exists(self.binfile)): return False
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
        if self.weighted_livetime:
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
        weighted = self.weighted_livetime
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
        for i in xrange(1+weighted):
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
        if self.weighted_livetime!=other.weighted_livetime:
            return DataError("DataSpec instances have inconsistent weighted_livetime")
        if self.mc_energy!=other.mc_energy:
            return DataError("DataSpec instances have inconsistent mc_energy")
        if self.mc_src_id!=other.mc_src_id:
            return DataError("DataSpec instances have inconsistent mc_src_id")
        return

    def add(self,other,output,binfile,ltcube):
        """Combine this DataSpec instance with another and return the result

        The two instances must have consistent definitions (DSS, binning,
        and livetime), and must have non-overlapping Gtis. The binfiles and
        ltcubes will be combined and written to the provided destinations;
        perhaps in the future, I will come up with sensible defaults.
        """
        exc = self.check_consistency(other)
        if exc is not None:
            raise(exc)
        gti = skymaps.Gti(self.gti)
        gti.intersection(other.gti)
        if gti.computeOntime()>0:
            raise DataError("DataSpec instances have overlapping GTIs")
        ft1 = self.ft1files+other.ft1files
        if self.ft2files == other.ft2files:
            ft2 = self.ft2files
        else:
            ft2 = self.ft2files+other.ft2files
        bpd = skymaps.BinnedPhotonData(self.binfile)
        bpd.add(skymaps.BinnedPhotonData(other.binfile))
        bpd.write(binfile)
        dssman.DSSEntries(self.binfile,header_key=0).write(binfile,header_key=0)
        fitstools.merge_lt([self.ltcube,other.ltcube],outfile=ltcube,weighted=self.weighted_livetime)
        dssman.DSSEntries(self.ltcube,header_key=0).write(ltcube,header_key=0)
        #TODO: move the copying of DSS entries into the merge_bpd and merge_lt functions
        gti_mask = skymaps.Gti(self.gti_mask)
        gti_mask.combine(other.gti_mask)
        kwargs = dict(ft1=ft1,
                      ft2=ft2,
                      binfile=binfile,
                      ltcube=ltcube,
                      gti_mask = gti_mask,
                      binsperdec=self.binsperdec,
                      mc_src_id=self.mc_src_id,
                      mc_energy=self.mc_energy,
                      weighted_livetime = self.weighted_livetime,
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
        self.bpd = skymaps.BinnedPhotonData(ds.binfile)
        self.lt = skymaps.LivetimeCube(ds.ltcube,weighted=False)
        if ds.weighted_livetime:
            self.weighted_lt = skymaps.LivetimeCube(ds.ltcube,weighted=True)
        self.gti = self.lt.gti() #Just to provide a reference.

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
               ,('dataspec',None,'One or more DataSpec instances'))
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
        if name is not None:
            self.filename = self._parse_name(name)
        if name is None or not os.path.exists(self.filename):
            if self.pickle is None and self.dataspec is None:
                raise DataError("Must specify either pickle files or DataSpecs")
        if os.path.exists(self.filename):
            self.dataspec = self._load_files(self.filename)
        else:
            if self.pickle is not None:
                if self.dataspec is not None:
                    raise DataError("Must specify dataspec OR pickle, not both.")
                self.dataspec = self._load_files(self.pickle)
            dump(self.dataspec,open(self.filename,'w'))

    def _parse_name(self,name):
        """Return the filename corresponding to a DataSet name"""
        if os.environ.has_key("FERMI"):
            data_dir = os.path.expandvars(os.path.join("$FERMI","data"))
            if not os.path.exists(data_dir):
                os.mkdir(data_dir)
        else:
            data_dir = os.getcwd()
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
        gfiles = sorted(glob(pfile))
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

class DataError(Exception):
    """Exception class for data management errors"""
    pass
