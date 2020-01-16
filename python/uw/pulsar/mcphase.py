"""A module for adding phase to photon events via Monte Carlo.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/mcphase.py,v 1.1 2010/11/01 19:54:18 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>

"""

from uw.utilities.fitstools import get_fields
from lcfitters import TaylorMapper
import numpy as N
import astropy.io.fits as pyfits

class MCPhase(object):
    """Support adding phase to photon events.
       
       The user *must* provide a light curve template that provides amplitude as
       a function of rotational phase.

       Optionally, the user may supply an instance of PhaseMapper.  This object
       implements a mapping between MET and rotational phase, e.g., by a .par
       file or by specification of F, Fdot, etc.  If this object is provided, the
       MET of the photon arrival will be "tweaked" to emulate its arrival from the
       pulsar.

       Can read/write from FT1 files, adding PULSE_PHASE if desired, modifying TIME
       if desired.
    """
    
    def __init__(self,lc_template,phase_mapper=None):
        """lc_template  -- an instance of LCTemplate
           phase_mapper -- an (optional) instance of PhaseMapper
        """
        self.lct = lc_template
        self.pm  = phase_mapper

    def process_edict(self,ed,new_lct=None,mc_src_id=0):
        """Replace (or add) the PULSE_PHASE entry of an edict using a
           (possibly new) light curve.  Useful for studying dependence
           on light curve morphology.
           
           Usage note: an "edict" is a dictionary with keys typical of
           an FT1 file.  In particular, the key MC_SRC_ID must exist.
        """
        lct = new_lct or self.lct
        m = ed['MC_SRC_ID']==mc_src_id
        pp = N.empty(len(m))
        pp[m] = lct.random(m.sum())
        m = ~m
        pp[m] = N.random.rand(m.sum())
        ed['PULSE_PHASE'] = pp

    def loadFT1(self,ft1files,mc_src_id):
        """Access data from one or more FT1 files.
           NB -- pyfits may choke on files bigger than a few hundred MB.

           ft1files  -- a single file name for an FT1 file, or a list of file names;
                        NOTA BENE -- writeback to multiple files is not supported

           mc_src_id -- the user should specify the MC_SRC_ID of the pulsar;
                        only its events will receive a "real" phase.
                        The user may pass "None", in which case all events
                        will have phase drawn from the light curve template.
        """

        keys = ['TIME'] + (['MC_SRC_ID'] if mc_src_id is not None else [])
        data = get_fields(ft1files,keys)
        for key in keys:
            data[key] = data[key].copy()

        self.ft1files = ft1files
        self.data = data
        self.mc_src_id = mc_src_id

        if self.pm is not None:
            # set phase0 to the middle of the observation
            self.pm.epoch = (self.data['TIME'].max() + self.data['TIME'].min())/2

    def processFT1(self):
        """Generate phases from the light curve template and adjust the photon
           arrival times if so specified."""

        # first, calculate the phase for the observations from the template
        times = self.data['TIME']
        self.phases = phases = N.empty(len(times))
        self.timeshifts = N.zeros(len(times))

        if self.mc_src_id is not None:
            mcids = self.data['MC_SRC_ID']
            self.mask = mask = mcids == self.mc_src_id

        phases[mask]  = self.lct.random(mask.sum())
        phases[~mask] = N.random.rand((~mask).sum())
        
        # next, if the user has provided a phase mapping, calculate the
        # adjustments to the time that need to be made
        # this calculation is done by calculating dphi/dt numerically

        # additionally, change the random pulse phases to the appropriate
        # ones from the ephemeris for the background photons
        if self.pm is not None:
            self.timeshifts[mask]  = (phases[mask]-self.pm(times[mask]))/self.pm.frequency(times[mask])
            phases[~mask] = self.pm(times[~mask])

    
    def writeFT1(self,phase_col_name   ='PULSE_PHASE',
                     adj_time_col_name ='TIME',
                     new_file_name     = None):
        """Write out data to FT1 file.

           phase_col_name  -- if not None, will add a pulse phase column of given
                              name.  The background  photons (as defined by 
                              MC_SRC_ID) will be folded at their default arrival time, 
                              whie those belonging to the pulsed source will be 
                              adjusted to come from the provided light curve template.

           adj_time_col_name  -- if not None, and if the user has provided a PhaseMapper
                              object, the algorithm will adjust the time of arrival
                              of the pulsed source photons to correspond to the
                              folded phase.  After this treatment, the FT1 file will
                              be suitable, e.g., for blind searches.  On the other
                              hand, if just checking for the efficacy of pulse-
                              detection algorithms for pulsars with timing solutions,
                              this step is not necessary.  The default is 'TIME', but
                              the user may choose a different name.

           new_file_name   -- if not None, will attempt to clobber the existing
                              FT1 file; note, if multiple FT1 files were provided,
                              the user *must* provided a new file name for output.
        """
        
        if hasattr(self.ft1files,'len'):
            if len(self.ft1files) == 1:
                self.ft1files = self.ft1files[0]
            else:
                print ('ERROR!')
                print ('Multiple FT1 files were passed in.')
                print ('Writeback to multiple files is not supported.')
                print ('Please process multiple FT1 files individually.')
                return

        new_file_name = new_file_name or self.ft1files

        try:
            f = pyfits.open(self.ft1files,memmap=True)
        except:
            f = pyfits.open(self.ft1files,memmap=False)

        newcols = []

        if adj_time_col_name is not None:
            if self.pm is not None:
                print ('Writing out adjusted event times to %s'%(adj_time_col_name))
                newtimes = self.data['TIME'] + self.timeshifts
                try:
                    f[1].data.field(adj_time_col_name)
                    f[1].data.field(adj_time_col_name)[:] = newtimes
                except KeyError:
                    col  = pyfits.Column(name=adj_time_col_name,format='D',array=newtimes)
                    newcols += [col]
            else:
                print ('ERROR!')
                print ('Found a column to write out adjusted time, but no')
                print ('PhaseMapper was provided to calculate the time shifts.')
                print ('No time entry will be written.')

        if phase_col_name is not None:
            print ('Writing out pulse phase column to %s'%(phase_col_name))
            try:
                f[1].data.field(phase_col_name)
                f[1].data.field(phase_col_name)[:] = self.phases
            except KeyError:
                col  = pyfits.Column(name=phase_col_name,format='D',array=self.phases)
                newcols += [col]
            
        if len(newcols) > 0:
            cols = f['EVENTS'].columns
            for newcol in newcols:
                cols = cols.add_col(newcol)
            t = pyfits.new_table(cols,header=f['EVENTS'].header)
            t.name = 'EVENTS'
            f[1] = t
        
        f.writeto(new_file_name,clobber=True)
        #tname = 'temp%d'%(int(N.random.rand(1)[0]*2**31))
        #f.writeto(tname)
        f.close()
        """
        # Windows' file protection is exceptionally stupid
        del(f)
        del(newcols)
        import gc
        gc.collect()
        gc.collect()
        import os
        try:
            os.remove(new_file_name)
        except:
            os.remove(new_file_name)
        os.rename(tname,new_file_name)
        """
