'''
Python module to plot phaseograms in different energy bands
'''
__author__  = 'Damien Parent'
__version__ = 1.5

import os, sys
import numpy, pyfits
from numpy import array, vectorize, sqrt, exp, log10, power
from itertools import izip
from ROOT import gStyle, TH1F, TH2F, TGraph, TText, TCanvas, TPad, TLegend, TLine

'''
from coords import sepangle_deg
import fits_utils as utilfits
import root_utils as root
from phasetools import h_statistic, weighted_h_statistic, z2m, sig2sigma
'''

from uw.utilities.coords import sepangle_deg
import uw.utilities.fits_utils as utilfits
import uw.utilities.root_utils as root
from uw.utilities.phasetools import weighted_h_statistic, h_statistic, z2m, sig2sigma

# ===================================================
# print in colors
# ===================================================
HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
OKRED = '\033[31m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'


#function to calculate the weight of a given event based on Philippe's tool
def get_weight(energy,angsep,theta,logeref,logesig):
    return exp(-(log10(energy)-logeref)**2/(2*logesig**2))*(1+angsep**2/(4./9.*theta**2))**-2.

def get_theta(par,energy):
    '''Get theta cut = (ROI x (e/100)**(alpha))**2 + beta**2'''
    return sqrt(par[0]**2 * power(100./energy,par[1]*2.) + par[2]**2)

def h_sig(h):
       """Convert the H-test statistic to a chance probability."""
       prob = exp(-0.4*h)
       if prob < 1e-130: prob = 1e-130 # to avoid saturation
       return prob

def progressbar( nentries, jentry ):
    # step --> 10%                                                                                                                                                                                 
    step1 = int(0.101 * nentries)
    step2 = int(0.02 * nentries)

    progress = int(float(jentry) / float(nentries) * 100.)

    if( jentry%step1 == 0 ):
        sys.stdout.write("|%i" %(progress)),
        sys.stdout.write('%'),
        sys.stdout.flush()

    if( jentry%step2 == 0 and jentry%step1 != 0 ):
        sys.stdout.write("."),; sys.stdout.flush()

    if( jentry == (nentries-1) ):
        print "|100% \n"
        sys.stdout.flush()

class Event(object):
    '''Declaration of an object called Event.'''
    def __init__(self, ra, dec, time, energy, evtclass, phase, zenith_angle, angsep):
        self.time         = time
        self.ra           = ra
        self.dec          = dec
        self.energy       = energy
        self.phase        = phase
        self.evtclass     = evtclass
        self.zenith_angle = zenith_angle
        self.angsep       = angsep

class PulsarLightCurve:
    '''Constructor of PulsarLightCurve'''
    def init(self):
        self.psrname       = 'PSR'
        self.ra_src        = None
        self.dec_src       = None
        self.radius        = None
        self.emin          = None
        self.emax          = None
        self.tmin          = None
        self.tmax          = None
        self.eclsmin       = 0
        self.eclsmax       = 1e9
        self.zenithcut     = 100.
        self.phase_colname = 'PULSE_PHASE'
        self.psfcut        = False
        self.convtype      = -1

    def __init__(self, ft1name, **kwargs):
        '''
        ============   =============================================
        keyword        description
        ============   =============================================
        ft1name        FT1 filename (necessary)
        ra_src         Right Ascension (deg)           [FT1 header]
        dec_src        Declination (deg)               [FT1 header]
        psrname        source name                     ["PSR"]
        radius         Maximum radius (deg)            [FT1 header]
        emin           Minimum energy (MeV)            [FT1 header]
        emax           Maximum energy (MeV)            [FT1 header]
        tmin           Start time (MET)                [FT1 header]
        tmax           End time (MET)                  [FT1 header]
        eclsmin        Minimum event class             [0]
        eclsmax        Maximum event class             [1e9]
        zenithcut      Zenith Cut (deg)                [100.]
        phase_colname  Phase column name               ["PULSE_PHASE"]
        psfcut         PSF-cut selection [True/False]  [False]
        convtype       -1=bothP6, 0=FrontP6, 1=BackP6  [-1]
                       -2=bothP7Source
        '''
        self.init()
        self.__dict__.update(kwargs)

        self.ft1name = ft1name

        print "================="
        print "  ", OKBLUE + self.psrname + ENDC
        print "================="

        # ====== time range =====
        ft1range = utilfits.ft1_time_range(ft1name)
        self.tmin = ft1range['START'] if self.tmin is None else self.tmin
        self.tmax = ft1range['STOP'] if self.tmax is None else self.tmax

        # ===== energy range =====
        self.__energy_range  = [ [100.,3e5] , [3000.,3e5], [1000.,3e5], [300.,1000.], [100.,300.], [3e3,3e5] ]
        self.emin = utilfits.get_header_erange(ft1name)[0] if self.emin is None else self.emin
        self.emax = utilfits.get_header_erange(ft1name)[1] if self.emax is None else self.emax

        # ====== position of the source (ra,dec) ======
        self.ra_src = utilfits.get_header_position(ft1name)[0] if self.ra_src is None else self.ra_src
        self.dec_src = utilfits.get_header_position(ft1name)[1] if self.dec_src is None else self.dec_src

        # ===== radius selection =====
        self.radius = utilfits.get_header_position(ft1name)[2] if self.radius is None else self.radius
        self.__radius_range = []
        
        # ======= psf selection =======
        self.__psf_selection = [5.12,0.8,0.07]
        if self.convtype == -1: self.__psf_selection = [5.12,0.8,0.07]
        elif self.convtype ==  0: self.__psf_selection = [3.3,0.8,0.05]
        elif self.convtype ==  1: self.__psf_selection = [5.5,0.8,0.1]
        elif self.convtype == -2: self.__psf_selection = [5.3,0.745,0.09]
        else: print "The selected convtype does not exist. Use the default convtype -1."

        # ======= print the selection =======
        self.print_selection()
        print ""

        # ===== declaration of histograms (ROOT ...) =====
        self.nbins = 32.
        self.phase_shift = 0.
        self.phaseogram = []
        self.pulse_phase = []
        self.phaseogram2D = []
        for i in range(6):
            histname = "phaseogram_" + str(i)
            self.phaseogram.append(TH1F())
            self.pulse_phase.append([])
            self.phaseogram2D.append(TH2F())
            self.__radius_range += [self.radius]

        self.weight = False
        self.__logeref = 0.
        self.__logesig = 0.
        
        # ===== Get the events from the FT1 file =====
        print "___pending___: loading photons from FT1 file ..."
        ft1file = pyfits.open(ft1name)
        hduname = 'EVENTS'

        ra = ft1file[hduname].data.field('RA')
        dec = ft1file[hduname].data.field('DEC')
        time = ft1file[hduname].data.field('TIME')
        energy = ft1file[hduname].data.field('ENERGY')
        evtclass = ft1file[hduname].data.field('EVENT_CLASS')
        zenith_angle = ft1file[hduname].data.field('ZENITH_ANGLE')
        angsep = sepangle_deg(self.ra_src,self.dec_src,ra,dec)
        
        try:
            phase = ft1file[hduname].data.field(self.phase_colname)
        except KeyError:
            print "\t You need to assign a phase for each event in the fits file. Exiting ..."
            exit()

        event_iter = izip(ra,dec,time,energy,evtclass,phase,zenith_angle,angsep)
        self.eventlist = [Event(*item) for item in event_iter]
        # ===========================================================================================
        # Before creating the event list, strip the FT1 file by selecting only events in a box/circle
        # centered at the pulsars location. This will greatly increase the speed of the code if
        # the ft1 file covers a large region of the sky compared to the maximum roi seeked.
        # ===========================================================================================
        self.eventlist = [item for item in self.eventlist if (
            (item.angsep <= self.radius) and (self.emin <= item.energy <= self.emax) and
            (self.tmin <= item.time <= self.tmax) and (self.eclsmin <= item.evtclass <= self.eclsmax) and
            (item.zenith_angle <= self.zenithcut) ) ]

    def set_psf_param(self,psfpar0=5.12,psfpar1=0.8,psfpar2=0.07):
        '''Set PSF parameterization: sqrt( (psfpar0x(100/E)^psfpar1)^2 + psfpar2^2 )'''
        self.__psf_selection = [psfpar0,psfpar1,psfpar2]
        print "Parameterization of the PSF has been modified:"
        print "PSF-cut: %.2f x (100/E)^%.2f ++ %.2f (++: in quadrature)\n" %(psfpar0,psfpar1,psfpar2)

    def set_energy_range(self,which=0,emin=100.,emax=300000.):
        if emin<self.emin:
            print OKRED + "warning: emin (for histo=%i) cannot be lower than PulsarLightCurve EMIN" %which + ENDC
            emin = self.emin
        if emax>self.emax:
            print OKRED + "warning: emax (for histo=%i) cannot be higher than PulsarLightCurve EMAX" %which + ENDC
            emax = self.emax
        self.__energy_range[which] = [emin,emax]
        print "Energy range for histo=%i has been changed: (%.0f,%.0f) MeV" %(which,emin,emax)

    def set_radius_range(self,which=0, radius=100.):
        '''Set the radius (deg) for the histo [which]'''
        if radius>self.radius:
            print OKRED + "warning: radius (for histo=%i) cannot be larger than PulsarLightCurve RADIUS" %which + ENDC 
            radius = self.radius
        self.__radius_range[which] = radius
        print "Radius selection for histo=%i has been changed: %.2f deg" %(which,radius)
        
    def set_weight_param(self,logeref=2.,logesig=0.5):
        self.__logeref = logeref
        self.__logesig = logesig

    def modify_loc(self,ra,dec):
        '''Set the source coordinates (ra,dec)'''
        self.ra_src  = ra
        self.dec_src = dec

    def fill_phaseogram(self, radius=None, nbins=32, bin_range=[0.,2.], phase_shift=0., weight=False):
        '''Fill both phaseograms (TH1F, TH2F) and pulse_phase objects.
        ===========    ===============================================
        keyword        description
        ===========    ===============================================
        nbins          number of bins in your histogram       [20]
        radius         Maximum radius (deg)                   [FT1 radius]
        bin_range      min and max phase interval             [0,2]
        phase_shift    add a shift to the pulse phase values  [0]
        weight         Use a weight for each photon           [False]
                       Must use the set_weight_param function
        ===========    ===============================================
        '''
        #print "\t(radmin,radmax): (", self.radmin, ",", self.radius, ") deg"
        
        print "\n___pending___: loading histograms ..."
        self.nbins = nbins
        self.phase_shift = phase_shift
        self.binmin = bin_range[0]
        self.binmax = bin_range[1]
        self.weight = weight

        # === radius ===
        if radius is None: radius = self.radius
        else:
            print "(radius): (", radius, ") deg"
            for i, item in enumerate(self.__radius_range):
                self.set_radius_range(i,radius)
        
        if radius>self.radius:
            print OKRED + "warning: radius cannot be larger than PulsarLightCurve RADIUS" + ENDC
            radius = self.radius
        
        # === initialization of histograms ===
        for item in self.phaseogram:
            item.Reset()
            item.SetBins(nbins*int(self.binmax-self.binmin),self.binmin,self.binmax)

        for item in self.phaseogram2D:
            item.Reset()
            item.SetBins(nbins,self.binmin,self.binmax,int((self.tmax-self.tmin)/(60*86400)),self.tmin,self.tmax)

        # ===========================================================================================
        # Before creating the event list, strip the FT1 file by selecting only events in a box/circle
        # centered at the pulsars location. This will greatly increase the speed of the code if
        # the ft1 file covers a large region of the sky compared to the maximum roi seeked.
        # ===========================================================================================  
        evtlist = [item for item in self.eventlist if item.angsep <= radius]
        
        energy = array([item.energy for item in evtlist])
        angsep = array([item.angsep for item in evtlist])
        theta  = array([get_theta(self.__psf_selection,item.energy) for item in evtlist])
        psfcut = (angsep<theta)
                
        if weight:
            print "--> using the weighting method with (LOGEREF,LOGESIG) = (%.1f,%.1f)" %(self.__logeref,self.__logesig)
            self.print_psf()
            vget_weight = vectorize(get_weight)
            myweight = vget_weight(energy,angsep,theta,self.__logeref,self.__logesig)
            if self.psfcut: print OKRED + "You are using the weighting methog and applying an energy-dependent cut" + ENDC

        erange = self.__energy_range
        rrange = self.__radius_range
        # ====================== LOOP OVER THE EVENTS  ======================
        #
        # ===================================================================
        for i, itevt in enumerate(evtlist):
            progressbar(len(evtlist),i)

            phase = itevt.phase + phase_shift
            if phase > 1.: phase -= 1.
            if phase < 0.: phase += 1.

            w = 1 if not weight else myweight[i]
            psf_filter = psfcut[i] if self.psfcut else True

            # --- light curve > energy threshold ---
            for j in range(len(self.phaseogram)):
                energy_filter = erange[j][0] <= itevt.energy <= erange[j][1]
                radius_filter = angsep[i] <= rrange[j]
                if energy_filter and radius_filter and psf_filter:
                    self.phaseogram[j].Fill(phase,w); self.phaseogram[j].Fill(phase+1.,w)
                    self.phaseogram2D[j].Fill(phase,itevt.time,w); self.phaseogram2D[j].Fill(phase+1,itevt.time,w)
                    self.pulse_phase[j] += [[phase,w]]

        for i, item in enumerate(self.pulse_phase):
            self.pulse_phase[i] = array(item)
        
	# to make a better display due to errorbars
        for histo in self.phaseogram:
            histo.SetMinimum(0.01)
            max_histo = histo.GetBinContent(histo.GetMaximumBin())
            if max_histo != 0:
                if weight:
                    factor = 1. + sqrt(max_histo)/max_histo
                    histo.SetMaximum(int(max_histo*factor))
                else:
                    factor = 1.1 + sqrt(max_histo)/max_histo
                    histo.SetMaximum(int(max_histo*factor))
                    if histo.GetMaximum()%5 == 0: histo.SetMaximum(histo.GetMaximum()+3)
            else:
                histo.SetMaximum(1.1)
        
    def get_basic_filter(self, event):
        '''Return True if the event goes through the cuts.'''
        return ( event.evtclass >= self.eclsmin and event.evtclass <= self.eclsmax and
                 event.energy >= self.emin and event.energy <= self.emax and
                 event.time >= self.tmin and event.time <= self.tmax and
                 event.zenith_angle <= self.zenithcut )

    def get_background(self, which=0, phase_range=[0.,1.], ring_range=[1.,2.], method='ring'):
        '''Return background/bin for a phaseogram.
        ===========   ======================================
        keyword       description
        ===========   ======================================
        which         phaseogram number    def:0
        phase_range   phase selection      def:[0,1]
        ring_range    ring dimension       def:[1,2]
        method        from a ring/weight   def:"ring"
        '''
        print "___pending___: evaluating background using %s ..." %method

        basic_events, psf_events, ring_events = 0., 0., 0.
        pmin, pmax = phase_range
        radius = self.__radius_range[which]
        theta_radmax = get_theta(self.__psf_selection,self.__energy_range[which][0])
        pulse_phase = []       
        
        if self.__psf_selection is None and self.psfcut:
            print OKRED + "Error, you must fill histograms before. Exiting..." + ENDC; sys.exit()
            
        if self.radius<ring_range[1]:
            print OKRED + "Error. Ring radius is larger than PulsarLightCurve RADIUS. Exiting ..." + ENDC; sys.exit()
        
        if self.psfcut and radius>theta_radmax:
            print "Warning your maximum PSF-cut radius is larger than your maximum radius"
            print "Setting the radius to the maximum PSF-cut radius:", theta_radmax, "deg"
            radius = theta_radmax

        # =============================================================
        # initialization of the list of events
        # =============================================================
        evtlist = self.eventlist
        angsep = array([item.angsep for item in evtlist])
        theta  = array([get_theta(self.__psf_selection,item.energy) for item in evtlist])
        psfcut = (angsep < theta)

        if method == 'weight':
            print "--> (LOGEREF,LOGESIG) = (%.1f,%.1f)" %(self.__logeref,self.__logesig)
            self.print_psf()
            vget_weight = vectorize(get_weight)
            energy = array([item.energy for item in evtlist])
            myweight = vget_weight(energy,angsep,theta,self.__logeref,self.__logesig)
        
        erange = self.__energy_range
        rrange = self.__radius_range
        # ====================== LOOP OVER THE EVENTS  ======================
        #
        # ===================================================================
        for i, itevt in enumerate(evtlist):
            progressbar(len(evtlist),i)

            phase = itevt.phase + self.phase_shift
            if phase > 1.: phase -= 1.
            if phase < 0.: phase += 1.
            
            energy_filter = erange[which][0] <= itevt.energy <= erange[which][1]
            phase_filter  = (phase>=pmin and phase<=pmax) if pmin<pmax else (phase>=pmin or phase<=pmax)
            radius_filter = (angsep[i] <= radius)
            psf_filter    = psfcut[i] if self.psfcut else True
            ring_filter   = ring_range[0] <= angsep[i] <= ring_range[1]

            if method is 'ring':
                if energy_filter and phase_filter:
                    if radius_filter: basic_events += 1
                    if psf_filter and radius_filter: psf_events += 1
                    if ring_filter: ring_events += 1
            if method == 'weight':
                if energy_filter and radius_filter and psf_filter:
                    pulse_phase += [[phase,myweight[i]]]

        # ======> results
        if method is 'ring':
            phase_norm = (pmax-pmin) if pmin<pmax else (1.+pmax-pmin)
            surf_norm = (power(ring_range[1],2)-power(ring_range[0],2)) / power(radius,2)
            return ring_events * ((psf_events/basic_events) / phase_norm / float(self.nbins) / surf_norm)
        elif method is 'weight':
            item = array(pulse_phase)
            return (len(item[:,0])-sum(item[:,1])) / float(self.nbins)
        else:
            print "The selected method does not exist."
            return 0.

    def add_weight_column(self,colname='WEIGHT',new_filename=None):
        '''add a column with the weight information in your FT1 file'''

        print "___pending___: adding the weight information in the FT! file ..."
        print "--> using (LOGEREF,LOGESIG) = (%.1f,%.1f)" %(self.__logeref,self.__logesig)
        
        # ============> calculate the weight for each event
        evtlist = [item for item in self.eventlist]
        energy = array([item.energy for item in evtlist])
        angsep = array([item.angsep for item in evtlist])
        theta  = array([get_theta(self.__psf_selection,item.energy) for item in evtlist])
        psfcut = (angsep<theta)
        vget_weight = vectorize(get_weight)
        myweight = vget_weight(energy,angsep,theta,self.__logeref,self.__logesig)

        ft1name = new_filename or self.__ft1name
        file = pyfits.open(ft1name)

        # clobber old values
        try:
            file['EVENTS'].data.field(colname)[:] = myweight
            print 'Clobbered old WEIGHT column.'
        except KeyError: 
            cols = file['EVENTS'].columns
            newcols = cols.add_col(pyfits.Column(name=colname,format='D',array=myweight))
            
            table = pyfits.new_table(newcols,header=file['EVENTS'].header)
            table.name = 'EVENTS'
            file[1] = table
            
        file.writeto(self.__ft1name,clobber=True)
        file.close()

    def print_psf(self):
        p0,p1,p2 = self.__psf_selection
        print "PSF-cut: %.2f x (100/E)^%.2f ++ %.2f (++: in quadrature)" %(p0,p1,p2)

    def print_selection(self):
        print "___SELECTION___:"
        print "\tFT1: %s" %self.ft1name
        print "\t(ra,dec): (", self.ra_src, ",", self.dec_src, ") deg"
        print "\t(tmin,tmax): (", self.tmin, ",", self.tmax, ") MET"
        print "\t(emin,emax): (", self.emin, ",", self.emax, ") MeV"
        print "\t(radius): (", self.radius, ") deg"
        print "\t(eclsmin,eclsmax): (", self.eclsmin, ",", self.eclsmax, ")"
        print "\t(zenith angle): ( <", self.zenithcut, ") deg"
        if self.psfcut:
            print "\t",
            self.print_psf()

    def print_summary(self):
        self.print_selection()
        print ""
        print "__parameters__| ",
        for i in range(len(self.phaseogram)): print "PHASO:" +  str(i), "\t",
        print ""
        for x in range(37): print "- ",
        print ""
        print "(emin,emax)MeV| ",
        for it, histo in enumerate(self.phaseogram):
            print '%(emin)d %(emax)d \t' %{'emin':self.__energy_range[it][0], 'emax':self.__energy_range[it][1]},
        print ""
        print "(radius)deg   | ",
        for it, histo in enumerate(self.phaseogram):
            print '%(radius).1f \t\t' %{'radius':self.__radius_range[it]},
        print ""
        print "number of evts| ",
        for histo in self.phaseogram:
            print "%.0f \t \t" %(histo.GetEntries()/2),
        print ""

        # ===========> htest
        htest = []
        for pulse in self.pulse_phase:
            if len(pulse) > 0.:
                H = weighted_h_statistic(pulse[:,0],pulse[:,1])
                Hsig = sig2sigma(h_sig(H))
                htest += [[H,Hsig]]
            else:
                htest += [[0,0]]
        print "H-test        | ",
        for item in htest:
            print "%.0f \t \t" %item[0],
        print ""
        print "H-test(sig)   | ",
        for item in htest:
            print "%.1f \t \t" %item[1],
        print ""
        print ""

    # ================ LOAD RADIO or XRAY PROFILE ================
    def load_profile(self, profile, fidpt=0., norm=True, mean_sub=False, ncol=1,
                     ylabel="Radio Flux (au)", comment="", histo=False):
        '''Load radio/xray profiles
        ___arguments___:
        profile     : radio or x-ray template file (ascii)
        fidpt       : fiducial point
        ncol        : column number for the radio template file
        mean_sub    :
        norm [True] : normalization
        ylabel      : ylabel title for the bottom panel
        comment     : comment on radio observation (e.g. 1.4 GHz)
        histo[False]: Change the display
        '''
        try:
            phase_ext = numpy.loadtxt(profile,comments="#")[:,ncol-1]
        except IndexError:
            phase_ext = numpy.loadtxt(profile,comments="#")
            
        phase_ext = numpy.roll(phase_ext,-int(len(phase_ext)*fidpt))  # fiducial point
        if mean_sub: phase_ext -= numpy.mean(phase_ext)
        if norm: phase_ext /= numpy.max(phase_ext)   # normalization

        if histo:
            template = TH1F()
            template.SetBins(len(phase_ext)*int(self.binmax-self.binmin),self.binmin,self.binmax)
            for i, item in enumerate(phase_ext):
                template.Fill(float(i)/len(phase_ext),item)
                template.Fill(1+float(i)/len(phase_ext),item)
            template.SetMinimum(0.01)
            max_histo = template.GetBinContent(template.GetMaximumBin())
            factor = 1.12 + sqrt(max_histo)/max_histo
            template.SetMaximum(int(max_histo*factor))
        else:
            template = TGraph(2*len(phase_ext))
            for i, item in enumerate(phase_ext):
                template.SetPoint(i,float(i)/len(phase_ext),item)
                template.SetPoint(i+len(phase_ext),1.+float(i)/len(phase_ext),item)
                template.GetXaxis().SetRangeUser(self.binmin,self.binmax)

        return template, ylabel, comment

    def plot_lightcurve( self, nbands=4, xdim=1100, ydim=1500, background=None,
                         zero_sup=False, inset=False, profile=None, profile2=None,
                         outdir="", outfile=None):
        '''ROOT function to plot gamma-ray phaseograms + (radio,x-ray)
        ===========   ================================================
        keyword       description
        ===========   ================================================
        nbands        number of gamma-ray bands in the figure [2,..,5]
        profile       radio or x-ray template (use load_profile)
        profile2      radio or x-ray template (use load_profile)
        backgroung    add a background level on the first panel [None]
        zero_sup      active the zero-suppress on the top panel [false]
        inset         add an inset on the top panel [false]
        xdim, ydim    plot dimensions
        outdir        output directory
        outfile       Name of the output file [None]
        '''
        print "___pending___: generating figure ..."
        root.initialization(batch=True)
        phaseogram = self.phaseogram
        erange = self.__energy_range

        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        pad = []
        text = TText()

        if nbands is 3:
            ylow, yup, ystep = 0.68, 1., 0.31
            if profile is not None: ylow, yup, ystep = 0.64, 1., 0.34
            
        if nbands is 4:
            ylow, yup, ystep = 0.67, 1., 0.21
            if profile is not None: ylow, yup, ystep = 0.67, 1., 0.22

        if nbands is 5:
            text.SetTextSize(0.09)
            ylow, yup, ystep = 0.7, 1., 0.16

        for i in range(nbands):
            padname = "pad" + str(i)
            if i+1 is nbands: ylow = 0
            pad.append(TPad(padname,padname,0,ylow,1.,yup)); pad[i].Draw()
            if i == 0: yup = ylow; ylow -= ystep
            else: yup -= ystep; ylow -= ystep
            pad[i].SetBottomMargin(0.); pad[i].SetTopMargin(0.);

        # ==========> labels
        xlabel_title = "Pulsar Phase"
        ylabel_title = "W. Counts" if self.weight else "Counts"
            
        # ==========> (0) TOP PANEL
	ecol = 0
        pad[ecol].cd(); pad[ecol].SetTopMargin(0.05);

        # title size, title offset, label size
        if nbands == 3: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.08,0.62,0.06); text.SetTextSize(0.06)
        if nbands == 4: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.072,0.69,0.055); text.SetTextSize(0.07)
        if nbands == 5: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.073,0.69,0.06); text.SetTextSize(0.08)

        if zero_sup:
            min_histo = phaseogram[ecol].GetBinContent(phaseogram[ecol].GetMinimumBin())
            factor = .8 - sqrt(min_histo)/min_histo
            min_level = int(min_histo*factor)
            if min_level%5 == 0 or min_level%5 == 4: min_level += 3
            phaseogram[ecol].SetMinimum(min_level)
            
        phaseogram[ecol].Draw()
        if not self.weight: phaseogram[ecol].Draw("Esame")
        tt = "> %.2f GeV" %(erange[ecol][0]/1e3)
        factor_level = 0.9 if zero_sup is False else 0.92
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*factor_level,tt)
	text.DrawText(1.5,phaseogram[ecol].GetMaximum()*factor_level,self.psrname)

        if background is not None:
            line = TLine(0,background,2,background)
            line.SetLineStyle(2); line.Draw()

        if inset:
            ecol = -1
            color = 14
            phaseogram[ecol].SetFillColor(color)
            phaseogram[ecol].Draw("same")
            text.SetTextColor(color)
            tt = "> " + str(erange[ecol][0]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[0].GetMaximum()*0.82,tt)
            text.SetTextColor(1)

        # =========> (1) PANEL
        ecol = 1
        if nbands > ecol:
            pad[ecol].cd();
            if nbands is 3: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.083,0.6,0.07); text.SetTextSize(0.09)
            if nbands is 4: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.11,0.45,0.09); text.SetTextSize(0.11)
            if nbands is 5: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.13,0.39,0.11); text.SetTextSize(0.13)
            phaseogram[ecol].Draw()
            if not self.weight: phaseogram[ecol].Draw("Esame")
            if erange[ecol][1] > 50000: tt = "> " + str(erange[ecol][0]/1e3) + " GeV"
            else: tt = str(erange[ecol][0]/1e3) + " - " + str(erange[ecol][1]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)

        # =========> (2) PANEL
        ecol = 2
        if nbands > ecol:
            pad[ecol].cd()
            if ecol == (nbands-1): pad[ecol].SetBottomMargin(0.2)
            if nbands is 3: root.set_axis_thx(phaseogram[ecol],xlabel_title,0.09,0.95,0.06,ylabel_title,0.068,0.72,0.05); text.SetTextSize(0.07)
            if nbands is 4: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.11,0.45,0.09); text.SetTextSize(0.11)
            if nbands is 5: root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,ylabel_title,0.13,0.39,0.11); text.SetTextSize(0.13)
            phaseogram[ecol].Draw()
            if not self.weight: phaseogram[ecol].Draw("Esame")
            if erange[ecol][1] > 50000: tt = "> " + str(erange[ecol][0]/1e3) + " GeV"
            else: tt = str(erange[ecol][0]/1e3) + " - " + str(erange[ecol][1]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)

        # ============> (3) PANEL
        ecol = 3
        if nbands > ecol:
            pad[ecol].cd()
            if ecol == (nbands-1): pad[ecol].SetBottomMargin(0.2) 
            if nbands is 4: root.set_axis_thx(phaseogram[ecol],xlabel_title,0.1,0.95,0.09,ylabel_title,0.092,0.55,0.08); text.SetTextSize(0.09)
            if nbands is 5: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,ylabel_title,0.13,0.39,0.11); text.SetTextSize(0.13)
            phaseogram[ecol].Draw()
            if not self.weight: phaseogram[ecol].Draw("Esame")
            tt = str(erange[ecol][0]/1e3) + " - " + str(erange[ecol][1]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.86,tt)

        # ==========> (4) PANEL
        ecol = 4
        if nbands > ecol:
            pad[ecol].cd()
            if ecol == (nbands-1): pad[ecol].SetBottomMargin(0.3)
            if nbands is 5: root.set_axis_thx(phaseogram[ecol],xlabel_title,0.11,1.1,0.09,ylabel_title,0.098,0.52,0.08); text.SetTextSize(0.095)
            phaseogram[ecol].Draw()
            if not self.weight: phaseogram[ecol].Draw("Esame")
            tt = str(erange[ecol][0]/1e3) + " - " + str(erange[ecol][1]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.85,tt)

        # ===========> TEMPLATE 1 (e.g. RADIO)
        if profile is not None:
            pad[-1].Clear()
            pad[-1].cd(); pad[-1].SetTopMargin(0.03); pad[-1].SetBottomMargin(0.)

            histo, ylabel, comment = profile
            # title size, title offset, label size
            if nbands is 3:
                pad[-1].SetTopMargin(0.02); pad[-1].SetBottomMargin(0.22)
                root.set_axis_thx(histo,xlabel_title,0.095,1.,0.08,ylabel,0.09,0.53,0.075)
                textsize = 0.095
            if nbands is 4:
                pad[-1].SetBottomMargin(0.22)
                root.set_axis_thx(histo,xlabel_title,0.1,0.95,0.09,ylabel,0.095,0.5,0.09)
                textsize = 0.11
            if nbands is 5:
                pad[-1].SetBottomMargin(0.28)
                root.set_axis_thx(histo,xlabel_title,0.11,1.1,0.09,ylabel,0.09,0.55,0.08)
                textsize = 0.12

            text.SetTextSize(textsize)
            
            if histo.__class__.__name__ == 'TGraph':
                histo.Draw("al")
                text.DrawText(0.1,0.9,comment)
            elif histo.__class__.__name__ == 'TH1F':
                histo.Draw(); histo.Draw("Esame")
                text.DrawText(0.1,histo.GetMaximum()*0.86,comment)
            else:
                print "Error in format. Exiting ..."; print histo.__class__.__name__
                sys.exit()
            
        # ===========> TEMPLATE 2 (e.g. X-RAY PROFILE)
        if profile2 is not None:
            pad[-2].Clear()
            pad[-2].cd(); pad[-2].SetTopMargin(0.02); pad[-2].SetBottomMargin(0.)

            histo, ylabel, comment = profile2
            # title size, title offset, label size
            if nbands is 3: root.set_axis_thx(histo,xlabel_title,0.11,1.0,0.11,ylabel,0.1,0.44,0.09)
            if nbands is 4:
                root.set_axis_thx(histo,xlabel_title,0.1,0.95,0.09,ylabel,0.11,0.43,0.09)
            if nbands is 5:
                root.set_axis_thx(histo,"",0.,0.,0.,"Counts",0.13,0.39,0.11)
                textsize = 0.13

            text.SetTextSize(textsize)
                
            if histo.__class__.__name__ == 'TGraph':
                histo.Draw("al")
                text.DrawText(0.1,0.8,comment)
            elif histo.__class__.__name__ == 'TH1F':
                histo.Draw(); histo.Draw("Esame")
                text.DrawText(0.1,histo.GetMaximum()*0.86,comment)
            else:
                print "Error in format. Exiting ..."; print histo.__class__.__name__
                sys.exit()

        # ===========> OUTPUT
        if outfile is None: outfilename = self.psrname + "_phaseogram_gamma_" + str(nbands) + "bands"
        else: outfilename = outfile
        canvas.Print(os.path.join(outdir,outfilename+".eps"))
        canvas.Print(os.path.join(outdir,outfilename+".root"))
        canvas.Print(os.path.join(outdir,outfilename+".png"))
        print ""
        
    def plot_phase_time(self, which=0, background=None, zero_sup=False, xdim=500, ydim=900,
                        outdir="", outfile=None):
        '''Plot phases as a function of the time.
        ___arguments___:
        which          : histogram number     Default: 0
        backgroung     : add a background level on the phaseogram [None]
        zero_sup       : active the zero-suppress on the top panel [False]
        xdim, ydim     : plot dimensions      Default: [500,800]
        outdir         : output directory
        outfile        : Name of the output file [None]   
        '''
        root.initialization()
        phaseogram = self.phaseogram
        phaseogram2D = self.phaseogram2D
        erange = self.__energy_range
        
        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        text = TText()
        pad1 = TPad("pad1","pad1",0,0.7,1,1); pad1.Draw()        
        pad2 = TPad("pad2","pad2",0,0,1,0.7); pad2.Draw()

        # ==========> labels
        xlabel_title = "Pulsar Phase"
        ylabel_title = "W. Counts" if self.weight else "Counts"
        
        # ==========> pad 1: phaseogram
        pad1.cd()
        pad1.SetTopMargin(0.05); pad1.SetLeftMargin(0.12); pad1.SetRightMargin(0.11)
        
        # title size, title offset, label size
        root.set_axis_thx(phaseogram[which],"",0.,0.,0.,ylabel_title,0.08,0.75,0.055)

        if zero_sup:
            min_histo = phaseogram[which].GetBinContent(phaseogram[which].GetMinimumBin())
            factor = .8 - sqrt(min_histo)/min_histo
            min_level = int(min_histo*factor)
            if min_level%5 == 0 or min_level%5 == 4: min_level += 3
            phaseogram[which].SetMinimum(min_level)

        phaseogram[which].Draw()
        if not self.weight: phaseogram[which].Draw("Esame")
        
        text.SetTextSize(0.07)
        tt = "> %.2f GeV" %(erange[which][0]/1e3)
        factor_level = 0.9 if zero_sup is False else 0.92
        text.DrawText(0.1,phaseogram[which].GetMaximum()*factor_level,tt)
        text.DrawText(1.4,phaseogram[which].GetMaximum()*factor_level,self.psrname)

        if background is not None:
            line = TLine(0,background,2,background)
            line.SetLineStyle(2)
            line.Draw()
            
        # ==========> pad2: phase vs time 
        pad2.cd()
        pad2.SetTopMargin(0.04); pad2.SetBottomMargin(0.1)
        pad2.SetLeftMargin(0.12); pad2.SetRightMargin(0.11)
        root.set_axis_thx(phaseogram2D[which],"Pulsar Phase",0.04,0.94,0.03,"Time [MET]",0.044,1.25,0.03)
        phaseogram2D[which].GetZaxis().SetLabelSize(0.03)
        phaseogram2D[which].Draw("CONT4 Z")
        
        # =======> OUTPUT
        if outfile is None: outfilename = self.psrname + "_phaseogram2D"
        else: outfilename = outfile
        canvas.Print(os.path.join(outdir,outfilename+".eps"))
        canvas.Print(os.path.join(outdir,outfilename+".root"))
        canvas.Print(os.path.join(outdir,outfilename+".png"))
        print ""

    def plot_ptest_time(self, phases, xdim=1000, ydim=600, outdir="."):
        '''Plot periodic test results as a funtion of the time.'''

        from scipy.stats import chi2
        root.initialization()
        
        nphotons = len(phases)*0.05
        nbins    = int(len(phases)/nphotons)
        htest    = TGraph(nbins)
        z2m_test = TGraph(nbins)

        for i in range(nbins):
            step = nphotons*(i+1)
            z2m_prob = chi2.sf(z2m(phases[0:step])[1],4)
            z2m_test.SetPoint(i,step,sig2sigma(z2m_prob))
            htest.SetPoint(i,step,sig2sigma(h_sig(h_statistic(phases[0:step]))))

        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        text = TText()
        pad1 = TPad("pad1","pad1",0,0,1,1); pad1.Draw()
        
        pad1.cd()
        pad1.SetTopMargin(0.05); pad1.SetBottomMargin(0.1)
        pad1.SetLeftMargin(0.1); pad1.SetRightMargin(0.05)
        root.set_axis_thx(htest,"Number of Photons",0.04,0.94,0.03,"Significance [\sigma]",0.05,0.8,0.03)
        htest.GetHistogram().SetMaximum(htest.GetHistogram().GetMaximum()*1.1)
        htest.SetLineWidth(2); htest.SetLineStyle(5); htest.Draw("apl")
        z2m_test.SetLineWidth(2); z2m_test.Draw("same")

        # significance
        if htest.GetHistogram().GetMaximum()>5: 
            line = TLine(0,5,step,5)
            line.SetLineStyle(4)
            line.SetLineColor(4)
            line.Draw()

        # legend
        xl1=.14; yl1=0.75; xl2=xl1+.15; yl2=yl1+.15
        leg = TLegend(xl1,yl1,xl2,yl2)
        leg.SetHeader(self.psrname)
        leg.AddEntry(htest,"H-test","l")
        leg.AddEntry(z2m_test,"z2m","l")
        leg.Draw()

        # =======> OUTPUT
        outfilename = self.psrname + "_ptest_time.eps"
        canvas.Print(os.path.join(outdir,outfilename))
        
    def toASCII(self, outdir="."):
        '''Copy the pulse phase of each phaseogram to an ascii file.

        ___arguments___:
        outdir : output directory
        '''
        erange = __energy_range
        
        filename = os.path.join(outdir,self.psrname+"_phaseogram.asc")
        outfile = open(filename,"w")
        outfile.write("# ------------------------------------------------------------------ \n")
        outfile.write("# Private results: LAT Collaboration and PTC ONLY \n")
        outfile.write("# Gamma-Ray Phase Histogram for %s \n" %self.psrname)
        outfile.write("# column: number of counts per energy band [MeV] \n")
        outfile.write("#%.0f-%.0f   %.0f-%.0f   %.0f-%.0f   %.0f-%.0f   %.0f-%.0f\n"
                      %(erange[0][0],erange[0][1],
                        erange[1][0],erange[1][1],
                        erange[2][0],erange[2][1],
                        erange[3][0],erange[3][1],
                        erange[4][0],erange[4][1])
                      )
        outfile.write("# ------------------------------------------------------------------ \n")

        for i in range(self.nbins):
            for histo in self.phaseogram:
                outfile.write("%.0f\t" %histo.GetBinContent(i+1))
            outfile.write("\n")

        print "INFO:", filename, "has been created."            

