"""
Module to provide access to CALDB information

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/irfs/caldb.py,v 1.2 2016/06/27 23:06:35 wallacee Exp $
Author: Eric Wallace

"""
__version__ = '$Revision: 1.2 $'

import os

import numpy as np
from astropy.io import fits

from uw.utilities import keyword_options
from . import IrfError

class CALDBError(IrfError):
    pass

class CALDB(object):

    """An object to keep track of the LAT CALDB information.
    
    Parameters
    ----------
        CALDB_dir
            Path to CALDB directory.  Environment variables will be expanded.
            Defaults to $CALDB. The specified directory must have a
            subdirectory 'data/glast/lat', containing a CALDB index file
            named 'caldb.indx'. 

    Attributes
    ----------
        CALDB_dir
            Absolute path to the directory specified by the `CALDB_dir` init
            parameter.
        index
            FITS HDU containing CALDB index data

    Methods
    -------
        list_irfs
            Prints a list of available irfs from the CALDB index, optionally
            filtered by selections specified with the same keyword arguments
            as __call__.
    """

    event_type_names = ('front','back', 'psf0','psf1','psf2','psf3','edisp0',
                         'edisp1','edisp2','edisp3')
    event_type_partitions = dict(fb = ('front','back'),
                                 psf = ('psf0','psf1','psf2','psf3'),
                                 edisp = ('edisp0','edisp1','edisp2','edisp3'))

    event_type_partitions = dict(fb = (0,1),
                                 psf = (2,3,4,5),
                                 edisp = (6,7,8,9))
    def __init__(self,CALDB_dir="$CALDB"):
        self.CALDB_dir = os.path.abspath(os.path.expandvars(CALDB_dir))
        self.index = self._load_caldb_index()
    
    def __repr__(self):
        return "{self.__class__}(CALDB_dir={self.CALDB_dir})".format(self=self)

    def _load_caldb_index(self):
        if not os.path.exists(self.CALDB_dir):
            raise CALDBError('CALDB directory {} not found.'.format(self.CALDB_dir))
        index_file = os.path.join(self.CALDB_dir,'data','glast','lat','caldb.indx')
        if os.path.exists(index_file):
            return fits.getdata(index_file,'CIF')
        elif os.path.exists(os.path.join(self.CALDB_dir,'caldb.indx')):
            return fits.getdata(os.path.join(self.CALDB_dir,'caldb.indx'),'CIF')
        raise CALDBError('No CALDB index found in {}'.format(self.CALDB_dir))
    
    def __call__(self,irf,version="P8R2_V6",event_class="source",event_type='fb'):
        """Return the filenames and FITS extensions for a given irf.
        
        Parameters
        ----------

        irf : {'psf','aeff','edisp'}
            The IRF to find files for.
        version : str
            IRF version to look up files for (e.g., 'P6_V1','P7_V6',
            'P7REP_V6','P8R2_V15'). Default is "P8R2_V6".
        event_class : str
            Event class selection to find IRFs for (e.g., "source", "clean",
            "diffuse"). Default is "source". Bitmask values are not currently
            recognized.
        event_type : str or int
            Event type ('front','back','psf0', etc) to find files for. Integer
            indices are also accepted (0-9 for ('front', 'back', 'psf0', ...,
            'edisp3'), respectively). May also be one of ('fb','psf','edisp'),
            in which case information will be returned for all event types in
            the specified set. Default is 'fb'.

        Returns
        -------

        irf_info : dict
            If a single event type is requested, the return value is a dict
            with keys "filename" and "extensions". "filename" is the path to
            the FITS file containing the description of the requested IRF.
            "extensions" is another dict in which the keys and values are the
            names and extensions numbers for each of the FITS extensions for
            the specified irf (e.g., 'RPSF', 'PSF_SCALING', and 'FISHEYE_CORR'
            for the PSF file).  

            If multiple event types are requested, the return value is a dict
            in which the keys are the event type names and the values are dicts
            of the form described above.
        """

        irf_entries = self._filter(irf=irf, version=version,
                                   event_class=event_class, event_type=event_type)
        def irf_dict(entries):
            #Assume all entries passed to this function are for the same file.
            filename = os.path.join(self.CALDB_dir,entries[0]['CAL_DIR'],
                                    entries[0]['CAL_FILE'])
            extensions = {cname:extension for cname,extension in 
                          zip(entries['CAL_CNAM'],entries['CAL_XNO'])}
            return dict(filename=filename,extensions=extensions)

        ets = np.asarray(self._parse_event_type(np.char.strip(irf_entries['DETNAM'])))
        ets_unique = np.unique(ets)
        if ets_unique.shape[0]==1:
            return irf_dict(irf_entries)
        else:
            return {et:irf_dict(irf_entries[ets==et]) for et in ets_unique}
        
        if event_type is not None:
            event_type = self._parse_event_type(event_type)
            if hasattr(event_type,'__iter__'):
                return {et:irf_dict(et) for et in event_type}
            else:
                return irf_dict(event_type)
    
    def _filter(self,irf=None,version=None,event_class=None,event_type=None):
        """Return CALDB index entries filtered according to the provided criteria."""
        def _parse_filenames(filenames):
            filenames = np.char.rpartition(np.char.lower(filenames),'.')[:,0]
            #Add a missing underscore for the P7 irfs for easier processing
            for ec in ('transient','source','clean','ultraclean'):
                filenames = np.char.replace(filenames,'p7'+ec,'p7_'+ec)
            tokens = np.char.strip(np.vstack(np.char.split(filenames,'_')))
            p6mask = (tokens[:,1]=='p6')
            #Handle the reversed version/event class order for P6 irfs
            tokens[p6mask,2],tokens[p6mask,3] = tokens[p6mask,3],tokens[p6mask,2]
            irf,pass_ver,ec,ver,et = tokens.T
            ver = np.array(['_'.join(x) for x in zip(pass_ver,ver)])
            return irf,ver,ec

        mask = np.ones_like(self.index,dtype='bool') 
        

        for val,selection in zip(_parse_filenames(self.index.CAL_FILE),(irf,version,event_class)):
            if selection is None:
                selection = val
            mask &= (val==np.char.lower(selection))

        if event_type is None:
            et = list(range(10))
        else:
            et = self._parse_event_type(event_type)
        if hasattr(et,"__iter__"):
            etmask = np.zeros_like(mask)
            for e in et:
                e = self.event_type_names[e]
                etmask |= (e==np.char.lower(np.char.strip(self.index.DETNAM)))
        else:
            et = self.event_type_names[et]
            etmask = (et==np.char.lower(np.char.strip(self.index.DETNAM)))

        mask &= etmask

        return self.index[mask]

    def _event_type_lookup(self,event_type):
        """Find the event type or types for a given event_type selection""" 
        if hasattr(event_type,'__iter__'):
            return [self._event_type_lookup(et) for et in event_type]
        try:
            et = event_type.lower()
        except AttributeError: #Not a string, assume index
            try:
                et = self.event_type_names[event_type]
            except IndexError:
                msg = "Invalid event type index {}. Should be 0-9.".format(
                        event_type, self.event_type_names)
                raise IrfError(msg)

        if et in self.event_type_names:
            return et

        try:
            return self.event_type_partitions[et]
        except KeyError:
            raise IrfError("Invalid event type {}".format(et))

    def _parse_event_type(self,event_type):
        """Find the event type or types for a given event_type selection""" 
        if hasattr(event_type,'__iter__'):
            return [self._parse_event_type(et) for et in event_type]

        try:
            event_type = event_type.lower()
        except AttributeError: #Not a string, assume index
            if event_type in range(10):
                return event_type
            else:
                raise IrfError("Invalid event type index {}. Should be 0-9.".format(event_type))

        try:
            return self.event_type_names.index(event_type)
        except ValueError as exc:
            try: return self.event_type_partitions[event_type]
            except KeyError:
                raise exc

    def list_irfs(self,**selections):
        """List available IRFs matching the given selections.

            Available selections are:

                version
                    IRF version (e.g. 'P8R2_V6','P6_V3','P7REP_V6')
                event_class
                    Event class (e.g. 'source','diffuse','clean')
                event_type_partition
                    Event type partition (one of 'fb', 'psf', 'edisp').
                    Overridden by event_type selection, if provided. 
                    Otherwise, all types of the selected partition are displayed.
                event_type
                    Event type ('front','back','psf0', etc). Integer indices
                    are also accepted (0-9 for ('front','back','PSF0',...,'EDISP3'))
        """

        def _format_names(names,indent=2,width=100):
            """Pretty print irf names in columns.

            irfnames
                Names to print.
            indent [2]
                Spaces to indent the list
            width [100]
                Total width of the display
            """
            wordwidth = names.dtype.itemsize
            names = np.unique(['{:<{w}s}'.format(n,w=wordwidth) for n in names])
            names.sort()
            sep = ' '*4
            ncols = (width-indent)//(wordwidth+4)
            names = np.append(names,[' '*wordwidth]*(ncols-names.shape[0]%ncols))
            names = names.reshape(ncols,names.shape[0]//ncols).T
            return [' '*indent+sep.join(ns) for ns in names]

        entries = self._filter(**selections)
        if entries.shape[0]==0:
            print("No available IRFs matching the provided selections.")
            return
        names = np.char.partition(entries.CAL_FILE,'_')[:,2]
        names = np.char.partition(names,'.')[:,0]
        if np.all([s is None for s in selections.values()]):
            st = 'Displaying all available IRFs:\n'
        else:
            st = 'Displaying available IRFs with ' 
            st += ', '.join(['{}={}'.format(k,v) for k,v in selections.items()
                              if v is not None])
            st += ":\n"
        st += "\n".join(_format_names(names))
        print(st)





