
import numpy as np
import os

class Profile(object):
    """ Encapsulate the profile with its identifying information."""

    """
    defaults = (
        ('ncol',2,'The column containing the Stokes I parameter.'),
        ('freq',1.4,'The observing frequency in GHz.')
    )
    """

    #@keyword_options.decorate(defaults)
    def __init__(self,profile,jname=None,**kwargs):
        """ profile -- the (ASCII) radio profile file """
        #keyword_options.process(self,kwargs)
        self.pfile = profile
        self.jname = jname
        self.obs = None

        ### default values -- will be updated automatically for profile type
        self.ncol = 2 # the ascii column with the Stokes I parameter (from 1)
        self.freq = 1.4 # observing freq in GHz
        self.fidpt = 0 # fiducial point of TOAs relative to profile
        ###

        self._process()

    def _process(self):
        """ Determine properties of profiles. """
        pfile = os.path.split(self.pfile)[-1]
        #0 -- special cases from collected profiles
        if self._process_special(pfile):
            pass
        #1 -- if file name ends with '.asc', it's a Nancay profile
        elif pfile.endswith('.asc'):
            self._process_nan()
        #2 -- JBO/PKS in the file name are giveaways
        elif 'JBO' in pfile:
            self._process_jbo()
        elif 'PKS' in pfile:
            self._process_pks()
        #3 -- one instance of DFB2 (parkes)
        elif 'PDFB2' in pfile:
            self.ncol = 4
            self._process_pks()
        #4 -- a Camilo/GBT style bestprof file
        elif pfile.endswith('bestprof') or ('gbt' in pfile):
            self._process_bestprof()
        else:
            raise ValueError('Could not discern type of %s'%pfile)
        self.fidpt = self.fidpt % 1 # just in case

    def _process_special(self,pfile):
        if 'J1446-4701' in pfile:
            # oddball parfile -- must be PKS?
            self.ncol = 4
            self.obs = 'PKS'
        elif '2043+1711' in pfile:
            # L-band AO profile from Lucas
            self.ncol = 2
            self.obs = 'AO'
            self.fidpt = 1-0.793468
        elif 'J0102+4839' in pfile:
            # profile from Megan, unknown convention
            self.ncol = 3
            self.obs = 'GBT'
        elif '2051-0827' in pfile:
            # WSRT profile from Cristobal, fiducial point at center
            self.ncol = 4
            self.fidpt = 0.5
        elif '1741+1351' in pfile:
            # AO profile from Cristobal, first harmonic convention
            self.ncol = 2
            self._first_harmonic()
        elif 'J2047+1053' in pfile:
            # 820 MHz profile from Scott, fiducial at 0
            self.ncol = 1
        elif 'J1124-5916' in pfile:
            # 1PC profile from PKS
            self.ncol = 4
        elif 'J1907+0602' in pfile:
            # profile from paper; AO L-band
            self.ncol = 1
        elif 'J1124-3653' in pfile:
            # profile from Paul, fiducial point at 0
            self.ncol = 2
        elif (('J0023+0923' in pfile) or
              ('J1810+1744' in pfile) or
              ('J2215+5135' in pfile)):
            # profiles by M. Kerr from PSRCHIVE files from J. Hessels
            self.ncol = 4
            self.obs = 'GBT'
            self.freq = 2.0
            self.fidpt = 0.5 # I *think* this is right
        else:
            return False
        return True

    def _process_nan(self):
        """ Nancay profiles properties:
            oo 1st harmonic convention (approximate)
            oo 1.4 GHz observing frequency
            oo fidicual_point specified in comments 
            oo amplitude given by 2nd column
        """
        self.obs = 'NAN'
        if 'J2241-5236' in self.pfile: # special case -- why?
            self.ncol = 4; self.fidpt = 0; return
        if 'J1125-5825' in self.pfile: # pcm
            self.ncol = 4; self.fidpt = 0; return
        for line in file(self.pfile):
            if 'fiducial_point' in line: break
        line = line.strip().split()[-2]
        self.fidpt = float(line)

    def _process_jbo(self):
        """ Jodrell Bank profiles properties:
            oo 1.4 GHz observing frequency
            oo fidicual_point specified in comments (first line; typically 0)
            oo amplitude given by 2nd column
        """
        self.obs = 'JBO'
        self.fidpt = float((file(self.pfile).next()).split()[-1])
        #if self.fidpt != 0:
            #print 'Found a JBO profile (%s) without fidpt=0 (at %.4f)'%(self.pfile,self.fidpt)

    def _process_pks(self):
        """ Parkes profiles properties:
            oo 1.4 GHz observing frequency
            oo fidicual_point is always 0
            oo amplitude given by 2nd column
        """
        # nothing to do for PKS
        self.obs = 'PKS'

    def _process_bestprof(self):
        """ 'bestprof' profiles properties:
            oo follow first harmonic conventions
            oo arbitrary observing frequency, typically 1.4
            oo arbitrary observatory, typically PKS
            oo amplitude given by 2nd column
        """
        # look for information encoded in .profile
        comments = [line for line in file(self.pfile) if line[0] == '#']
        for comment in comments:
            if 'obs' in comment:
                self.obs = comment.split()[-1]
            elif 'freq' in comment:
                self.freq = comment.split()[-1]
        # otherwise, assume defaults
        self._first_harmonic()

    def _first_harmonic(self):
        """ Compute the zero of phase of a radio profile by determining the 
            position of the fundamental peak."""
        TWOPI = 2*np.pi
        self.fidpt = 0
        vals = self.get_amplitudes(align_to_peak=False)[0]
        ph = np.linspace(0,TWOPI,len(vals)+1)[:-1] # LEFT bin edges
        a1 = (np.sin(ph)*vals).sum()
        a2 = (np.cos(ph)*vals).sum()
        self.fidpt = np.arctan2(a1,a2)/TWOPI

    def get_amplitudes(self,align_to_peak=True,bin_goal=None,phase_shift=0,phases=None):
        """ Produce an ordered list of amplitudes as a function of phase,
            rotated such that the fiducial point is at 0.  The profile 
            is re-binned with linear interpolation.

            align_to_peak [True] -- will rotate the profile to place the
            peak bin at zero, if not already there
            
            bin_goal [None] -- will attempt to average the light curve to
                between between bin_goal and 2xbin_goal bins"""
        if phases is None:
            try:
                phases = np.loadtxt(self.pfile,comments='#')[:,self.ncol-1]
            except IndexError:
                phases = np.loadtxt(self.pfile,comments='#')
            if len(phases.shape) > 1:
                raise ValueError('Could not read profile values.')

        if align_to_peak:
            a = np.argmax(phases)
            align_shift = float(a) / len(phases) # NB -- really peak position
            phases = np.roll(phases,-a)
        else:
            align_shift = 0

        rvals = phases
        if phase_shift != 0:
            # interpolate the profile so that the fiducial point is a zero
            x0 = np.linspace(0,1,len(phases)+1)[:-1]
            x = np.concatenate((x0-1,x0,x0+1,[2]))
            y = np.concatenate((phases,phases,phases,[phases[0]]))
            rvals = np.interp(x0-phase_shift,x,y)

        if bin_goal > 0:
            if len(rvals) % 2 > 0: 
                rvals = np.append(rvals,rvals[0])
            while len(rvals) > bin_goal:
                rvals = (rvals[:-1:2]+rvals[1::2])/2
        #return rvals,self.fidpt+align_shift,self.fidpt-align_shift
        return rvals,self.fidpt-align_shift

class FullStokesProfile(Profile):
    """ Handle a pulsar profile with all Stokes parameters."""

    def init(self):
        # the columns in the text file for I, Q, U, V
        self.fscol = [4,5,6,7,8,9]

    def __init__(self,profile,jname=None,**kwargs):
        """ profile -- the (ASCII) radio profile file """
        # TO DO -- switch this so that it calls pdv
        #keyword_options.process(self,kwargs)
        self.init()
        self.__dict__.update(kwargs)
        self.pfile = profile
        self.jname = jname
        self.obs = None

        ### default values -- will be updated automatically for profile type
        self.ncol = 2 # the ascii column with the Stokes I parameter (from 1)
        self.freq = 1.4 # observing freq in GHz
        self.fidpt = 0 # fiducial point of TOAs relative to profile
        ###

    def _get_rms(self):
        # way fragile
        return float(file('/tmp/tmp').readlines()[0].split('RMS: ')[-1])
    
    def _get_data(self):
        #try:
        data = np.loadtxt(self.pfile,comments='#')
        return np.asarray([data[:,col-1] for col in self.fscol])
        #except:
        #    print 'Could not read data!'
        #    raise ValueError

    def get_intensity(self,**kwargs):
        i,q,u,v,pa,pae = self._get_data()
        return self.get_amplitudes(phases=i,**kwargs)
        
    def get_position_angles(self,**kwargs):
        i,q,u,v,pa,pae = self._get_data()
        #return self.get_amplitudes(phases=0.5*np.degrees(np.arctan2(u,q)),**kwargs)
        return pa,pae

    def get_linear_intensity(self,**kwargs):
        i,q,u,v,pa,pae = self._get_data()
        return self.get_amplitudes(phases=(q**2+u**2)**0.5,**kwargs)

    def get_circular_intensity(self,**kwargs):
        i,q,u,v,pa,pae = self._get_data()
        return self.get_amplitudes(phases=v,**kwargs)
        
        
