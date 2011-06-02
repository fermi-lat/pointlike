'''
Python module for plotting phaseograms in different energy bands
'''
__author__  = 'Damien Parent'
__version__ = 1.3

import os, sys
import numpy, pyfits
import itertools
from ROOT import gStyle, TH1F, TH2F, TGraph, TText, TCanvas, TPad, TLegend, TLine

from uw.utilities.coords import sepangle_deg
import uw.utilities.fits_utils as utilfits
import uw.utilities.root_utils as root                                                                                                                                                            
#import root_utils as root
from uw.utilities.phasetools import h_statistic, z2m, sig2sigma


def set_psfcut(par,energy):
    '''Apply a psf cut of the form (ROI x (e/100)**(alpha))**2 + beta**2'''
    theta = numpy.sqrt( numpy.power(par[0],2.) * numpy.power(energy/100.,par[1]*2.) + numpy.power(par[2],2.) )
    return theta

def h_sig(h):
       """Convert the H-test statistic to a chance probability."""
       prob = numpy.exp(-0.4*h)
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
    def __init__(self, ra, dec, time, energy, evtclass, phase, zenith_angle ):
        self.time         = time
        self.ra           = ra
        self.dec          = dec
        self.energy       = energy
        self.phase        = phase
        self.evtclass     = evtclass
        self.zenith_angle = zenith_angle

class PulsarLightCurve:
    '''Constructor of PulsarLightCurve
    Put more information ...
    '''

    def init(self):
        '''
        Initialization:
        ---------------
             psrname, radus, emin, emax, tmin, tmax,
             eclsmin, eclsmax, zenithcut, phase_colname
        '''
        self.psrname       = 'PSR'
        self.radius        = None
        self.radmin        = 0.
        self.emin          = None
        self.emax          = None
        self.tmin          = None
        self.tmax          = None
        self.eclsmin       = 0
        self.eclsmax       = 1e5
        self.zenithcut     = 105.
        self.phase_colname = 'PULSE_PHASE'
        self.psfcut        = True
        self.convtype      = -1

    def __init__(self, ft1name, ra_src=None, dec_src=None, **kwargs):
        '''
        ___arguments___:
        ft1name   : FT1 filename (necessary)
        ra_src    : Right Ascension (deg). Default: FT1 header.
        dec_src   : Declination (deg).     Default: FT1 header.

        options:
        --------
        psrname       : source name                      Default: PSR
        radius        : Maximum radius (deg).            Default: FT1 header.
        radmin        : Keep all the photons in x deg    Default: 0 (deg)
        emin          : Minimum energy (MeV).            Default: FT1 header.
        emax          : Maximum energy (MeV).            Default: FT1 header.
        tmin          : Start time (MET).                Default: FT1 header.
        tmax          : End time (MET).                  Default: FT1 header.
        eclsmin       : Minimum event class              Default: 0
        eclsmax       : Maximum event class              Default: 1e5
        zenithcut     : Zenith Cut (deg)                 Default: 105.
        phase_colname : Phase column name.               Default: PULSE_PHASE
        psfcut        : PSF-cut selection [True/False]   Default: True
        convtype      : -1=both, 0=Front, 1=Back         Defalut: -1
        '''

        self.init()
        self.__dict__.update(kwargs)
        self.__ra_src  = ra_src
        self.__dec_src = dec_src

        self.psf_selection = None
        if self.convtype == -1: self.psf_selection = [ 5.12 , -0.8 , 0.07 ]
        elif self.convtype ==  0: self.psf_selection = [ 3.3 , -0.8 , 0.05 ]
        elif self.convtype ==  1: self.psf_selection = [ 5.5 , -0.8 , 0.1 ]
        else: print "convtype that you selected does not exist. Exiting ..."
                        
        self.energy_range  = [ [100.,3e5] , [3000.,3e5], [1000.,3e5], [300.,1000.], [100.,300.], [3e3,3e5] ]

        print "================="
        print "  ", self.psrname
        print "================="
        print "___PARAMETERS___:"

        # ====== time range =====
        ft1range = utilfits.ft1_time_range(ft1name)
        if self.tmin is None: self.tmin = ft1range['START']
        if self.tmax is None: self.tmax = ft1range['STOP']
        print "\t(tmin,tmax): (", self.tmin, ",", self.tmax, ") MET"

        # ===== energy range =====
        if self.emin is None: self.emin = utilfits.get_header_erange(ft1name)[0]
        if self.emax is None: self.emax = utilfits.get_header_erange(ft1name)[1]
        print "\t(emin,emax): (", self.emin, ",", self.emax, ") MeV"

        # ====== position of the source (ra,dec) ======
        if ra_src is None: self.__ra_src = utilfits.get_header_position(ft1name)[0]
        if dec_src is None: self.__dec_src = utilfits.get_header_position(ft1name)[1]
        print "\t(ra,dec): (", self.__ra_src, ",", self.__dec_src, ") deg"

        # ===== radius selection =====
        if self.radius is None: self.radius = utilfits.get_header_position(ft1name)[2]
        print "\t(radmin,radmax): (", self.radmin, ",", self.radius, ") deg"

        # ===== PSF ====
        if self.psfcut:
            self.print_psf()

        print ""

        # ===== declaration of histograms (ROOT ...) =====
        self.phaseogram = []
        self.pulse_phase = []
        self.phaseogram2D = []
        for i in range(6):
            histname = "phaseogram_" + str(i)
            self.phaseogram.append(TH1F())
            self.pulse_phase.append(numpy.array([]))
            self.phaseogram2D.append(TH2F())

        self.nbins = 0
        self.binmin, self.binmax = 0, 0

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

        try:
            phase = ft1file[hduname].data.field(self.phase_colname)
        except KeyError:
            print "\t You need to assign a phase for each event in the fits file. Exiting ..."
            exit()

        event_iter = itertools.izip(ra, dec, time, energy, evtclass, phase, zenith_angle)
        self.eventlist = [Event(*row) for row in event_iter]

    def modify_loc(self,ra,dec):
        '''Set the coordinates (ra,dec) of your source.'''
        self.__ra_src  = x
        self.__dec_src = y

    def fill_phaseogram(self, nbins=20, bin_range=[0.,2.], phase_shift=0.):
        '''Fill both phaseograms (TH1F, TH2F) and pulse_phase objects.

        parameters:
        -----------
        nbins       : number of bins in your histogram
        bin_range   : min and max phase interval
        phase_shift : add a phase shift in your histogram
        '''
        print "\n___pending___: loading histograms ..."
        self.nbins = nbins
        self.phase_shift = phase_shift
        self.binmin = bin_range[0]
        self.binmax = bin_range[1]

        # === initialization of histograms ===
        for item in self.phaseogram:
            item.Reset()
            item.SetBins(nbins*int(self.binmax-self.binmin),self.binmin,self.binmax)

        for item in self.phaseogram2D:
            item.Reset()
            item.SetBins(nbins,0.,1.,int((self.tmax-self.tmin)/(14*86400)),self.tmin,self.tmax)

        for i, item in enumerate(self.pulse_phase):
            self.pulse_phase[i] = numpy.array([])

        # ====================== LOOP OVER THE EVENTS  ======================
        #
        # ===================================================================
        for i, it_events in enumerate(self.eventlist):
            progressbar(len(self.eventlist),i)

            phase = it_events.phase + self.phase_shift
            if phase > 1.: phase -= 1.
            if phase < 0.: phase += 1.

            angsep = sepangle_deg(self.__ra_src,self.__dec_src,it_events.ra,it_events.dec)
            theta = set_psfcut(self.psf_selection,it_events.energy)

            basic_filter = self.get_basic_filter(it_events)
            radius_filter = ( angsep <= self.radius )
            psf_filter    = (angsep<theta or angsep<self.radmin) if self.psfcut else True

            # --- light curve > energy threshold ---
            if basic_filter and radius_filter and psf_filter:
                for j in range(len(self.phaseogram)):
                    if it_events.energy > self.energy_range[j][0] and it_events.energy < self.energy_range[j][1]:
                        self.phaseogram[j].Fill(phase); self.phaseogram[j].Fill(phase+1.)
                        self.pulse_phase[j] = numpy.append(self.pulse_phase[j],phase)
                        self.phaseogram2D[j].Fill(phase,it_events.time); self.phaseogram2D[j].Fill(phase+1,it_events.time)

	# to make a better display due to errorbars
        for histo in self.phaseogram:
            histo.SetMinimum(0.01)
            max_histo = histo.GetBinContent(histo.GetMaximumBin())
            factor = 1.12 + numpy.sqrt(max_histo)/max_histo
            histo.SetMaximum(int(max_histo*factor))
            if histo.GetMaximum()%5 == 0:
		histo.SetMaximum(histo.GetMaximum()+3)

    def get_basic_filter(self, event):
        '''Return True if the event goes through the cuts.'''
        basic_filter = ( event.evtclass >= self.eclsmin and event.evtclass <= self.eclsmax and
                         event.energy >= self.emin and event.energy <= self.emax and
                         event.time >= self.tmin and event.time <= self.tmax and
                         event.zenith_angle <= self.zenithcut )

        return basic_filter

    def get_background(self, which=0, phase_range=[0.,1.], ring_range=[1.,2.], psfcut=True):
        '''Return background/bin for a phaseogram.

        parameters:
        -----------
        which       : phaseogram number
        phase_range : phase selection
        ring_range  : selection of the ring [degrees]

        return:
        -------
        background/bin
        '''

        print "___pending___: evaluating background ..."

        if self.psf_selection is None and self.psfcut:
            print "Error, you have to fill histograms before. Exiting..."
            sys.exit()

        basic_events = 0.; psf_events = 0.; ring_events = 0.
        pmin = phase_range[0]; pmax = phase_range[1]

        for it_events in self.eventlist:

            phase = it_events.phase + self.phase_shift
            if phase > 1.: phase -= 1.
            if phase < 0.: phase += 1.

            angsep = sepangle_deg(self.__ra_src,self.__dec_src,it_events.ra,it_events.dec)
            theta = set_psfcut(self.psf_selection,it_events.energy)

            basic_filter  = self.get_basic_filter(it_events)
            energy_filter = it_events.energy >= self.energy_range[which][0] and it_events.energy <= self.energy_range[which][1]
            phase_filter  = (phase>=pmin and phase<=pmax) if pmin<pmax else (phase>=pmin or phase<=pmax)
            radius_filter = angsep <= self.radius
            psf_filter    = (angsep<theta or angsep<self.radmin) if self.psfcut else True
            ring_filter   = angsep >= ring_range[0] and angsep <= ring_range[1]

            if basic_filter and energy_filter and phase_filter:
                if radius_filter                : basic_events += 1
                if psf_filter and radius_filter : psf_events += 1
                if ring_filter                  : ring_events += 1

        phase_norm = (pmax-pmin) if pmin<pmax else (1.+pmax-pmin)
        surf_norm = (numpy.power(ring_range[1],2)-numpy.power(ring_range[0],2)) / numpy.power(self.radius,2)
        return ring_events * ((psf_events/basic_events) / phase_norm / float(self.nbins) / surf_norm)

    def print_psf(self):
        sys.stdout.write("\tPSF-cut: "),
        print self.psf_selection[0],
        sys.stdout.write(" x (E/100)^"),
        print self.psf_selection[1],
        sys.stdout.write(" ++ "),
        print self.psf_selection[2],
        print "(++: in quadrature)"

    def print_summary(self):
	print "___PARAMETERS___:"
        print "\t(tmin,tmax): (", self.tmin, ",", self.tmax, ") MET"
        print "\t(emin,emax): (", self.emin, ",", self.emax, ") MeV"
        print "\t(ra,dec): (", self.__ra_src, ",", self.__dec_src, ") deg"
        print "\t(radmin,radmax): (", self.radmin, ",", self.radius, ") deg"
        if self.psfcut: self.print_psf()
        print ""
        
        print "PARAMETERS | \t",
        for i in range(len(self.phaseogram)): print "PHASO_" +  str(i), "\t",
        print ""
        for x in range(37): print "- ",
        print ""
        print "ENERGY(MeV)| \t",
        for it, histo in enumerate(self.phaseogram):
            print '%(emin)d %(emax)d \t' %{'emin':self.energy_range[it][0], 'emax':self.energy_range[it][1]},
        print ""
        print "NB OF EVTS | \t",
        for histo in self.phaseogram:
            print "%.0f \t \t" %(histo.GetEntries()/2),
        print ""
        print "HTEST      | \t",
        for pulse in self.pulse_phase:
            if len(pulse) > 0.:
                print "%.2f \t \t" %(h_statistic(pulse)),
            else:
                print "0 \t \t",
        print ""
        print "HTEST(sig) | \t",
        for pulse in self.pulse_phase:
            if len(pulse) > 0.:
                print "%.2f \t \t" %(sig2sigma(h_sig(h_statistic(pulse)))),
            else:
                print "0 \t \t",
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
            factor = 1.12 + numpy.sqrt(max_histo)/max_histo
            template.SetMaximum(int(max_histo*factor))
        else:
            template = TGraph(2*len(phase_ext))
            for i, item in enumerate(phase_ext):
                template.SetPoint(i,float(i)/len(phase_ext),item)
                template.SetPoint(i+len(phase_ext),1.+float(i)/len(phase_ext),item)
                template.GetXaxis().SetRangeUser(self.binmin,self.binmax)

        return template, ylabel, comment

    
    def plot_lightcurve( self, nbands=4, xdim=1200, ydim=1500, outdir="", background=None,
                         background_sub=False, inset=False, profile=None, profile2=None ):
        '''ROOT function for ploting x gamma-ray bands.
        ___arguments___:
        nbands         : number of gamma-ray bands in the figure [2,..,5]
        profile        : radio or x-ray template (use load_profile)
        profile2       : radio or x-ray template (use load_profile)
        backgroung     : add a background level on the first panel
        background_sub : active the zero-suppress on the top panel [false]
        inset          : add an inset on the top panel [false]
        xdim, ydim     : plot dimensions
        outdir         : output directory
        '''
        root.initialization(batch=True)
        phaseogram = self.phaseogram

        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        pad = []
        text = TText()

        if nbands is 3:
            text.SetTextSize(0.08)
            ylow, yup, ystep = 0.68, 1., 0.31
            if profile is not None: ylow, yup, ystep = 0.64, 1., 0.34
            
        if nbands is 4:
            text.SetTextSize(0.09)
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
            
        # ==========> (0) TOP PANEL
	ecol = 0
        pad[ecol].cd(); pad[ecol].SetTopMargin(0.05); pad[ecol].SetBottomMargin(0.)

        # title size, title offset, label size
        if nbands == 3:
            root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.08,0.62,0.06)
        if nbands == 4:
            root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.073,0.69,0.053)
        if nbands == 5:
            root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.073,0.69,0.06)

        if background_sub:
            min_histo = phaseogram[ecol].GetBinContent(phaseogram[ecol].GetMinimumBin())
            factor = .8 - numpy.sqrt(min_histo)/min_histo
            phaseogram[ecol].SetMinimum(int(min_histo*factor))
            
        phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
        tt = "> %.2f GeV" %(self.energy_range[ecol][0]/1e3)
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.9,tt)
	text.DrawText(1.5,phaseogram[ecol].GetMaximum()*0.9,self.psrname)

        if background is not None:
            line = TLine(0,background,2,background)
            line.SetLineStyle(2)
            line.Draw()

        if inset:
            ecol = -1
            color = 14
            phaseogram[ecol].SetFillColor(color)
            phaseogram[ecol].Draw("same")
            text.SetTextColor(color)
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[0].GetMaximum()*0.82,tt)
            text.SetTextColor(1)

        # =========> (1) PANEL
        ecol = 1
        if nbands > ecol:
            pad[ecol].cd(); pad[ecol].SetTopMargin(0.); pad[ecol].SetBottomMargin(0.)

            if nbands is 3:
                root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.083,0.6,0.07)
                textsize = 0.09
            if nbands is 4:
                root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.11,0.42,0.09)
                textsize = 0.11
            if nbands is 5:
                root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.13,0.39,0.11)
                textsize = 0.13

            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            if self.energy_range[ecol][1] > 50000: tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            else: tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
            text.SetTextSize(textsize)
            phaseogram[ecol].SetFillColor(0)
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)

        # =========> (2) PANEL
        ecol = 2
        if nbands > ecol:
            pad[ecol].cd(); pad[ecol].SetTopMargin(0.); pad[ecol].SetBottomMargin(0.)

            if nbands is 3:
                pad[ecol].SetBottomMargin(0.2)
                root.set_axis_thx(phaseogram[ecol],"Pulsar Phase",0.09,0.95,0.06,"Counts",0.068,0.72,0.05)
                tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
                textsize = 0.07
            if nbands is 4:
                root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.11,0.42,0.09)
                tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
                textsize = 0.11
            if nbands is 5:
                root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.13,0.39,0.11)
                tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
                textsize = 0.13

            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            text.SetTextSize(textsize)
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)

        # ============> (3) PANEL
        ecol = 3
        if nbands > ecol:
            pad[ecol].cd(); pad[ecol].SetTopMargin(0.); pad[ecol].SetBottomMargin(0.)

            if nbands is 4:
                pad[ecol].SetBottomMargin(0.2)
                root.set_axis_thx(phaseogram[ecol],"Pulsar Phase",0.1,0.95,0.09,"Counts",0.092,0.55,0.08)
                textsize = 0.09
            if nbands is 5:
                root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.13,0.39,0.11)
                textsize = 0.13

            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            text.SetTextColor(1); text.SetTextSize(textsize)
            tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.86,tt)

        # ==========> (4) PANEL
        ecol = 4
        if nbands > ecol:
            pad[ecol].cd(); pad[ecol].SetTopMargin(0.); pad[ecol].SetBottomMargin(0.)
            
            if nbands is 5:
                pad[ecol].SetBottomMargin(0.3)
                root.set_axis_thx(phaseogram[ecol],"Pulsar Phase",0.11,1.1,0.09,"Counts",0.096,0.53,0.08)
                textsize = 0.095
                
            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
            text.SetTextSize(textsize)
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.85,tt)

        # ===========> TEMPLATE 1 (e.g. RADIO)
        if profile is not None:
            pad[-1].Clear()
            pad[-1].cd(); pad[-1].SetTopMargin(0.03); pad[-1].SetBottomMargin(0.)

            histo, ylabel, comment = profile
            # title size, title offset, label size
            if nbands is 3:
                pad[-1].SetTopMargin(0.02)
                pad[-1].SetBottomMargin(0.22)
                root.set_axis_thx(histo,"Pulsar Phase",0.095,1.,0.08,ylabel,0.09,0.53,0.075)
                textsize = 0.095
            if nbands is 4:
                pad[-1].SetBottomMargin(0.22)
                root.set_axis_thx(histo,"Pulsar Phase",0.1,0.95,0.09,ylabel,0.095,0.5,0.09)
                textsize = 0.11
            if nbands is 5:
                pad[-1].SetBottomMargin(0.28)
                root.set_axis_thx(histo,"Pulsar Phase",0.11,1.1,0.09,ylabel,0.09,0.55,0.08)
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
            if nbands is 3: root.set_axis_thx(histo,"Pulsar Phase",0.11,1.0,0.11,ylabel,0.1,0.44,0.09)
            if nbands is 4:
                root.set_axis_thx(histo,"Pulsar Phase",0.1,0.95,0.09,ylabel,0.11,0.43,0.08)
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
        outfilename = self.psrname + "_phaseogram_gamma_" + str(nbands)
        eps_outfilename = outfilename + "bands.eps"
        canvas.Print(os.path.join(outdir,eps_outfilename))
        root_outfilename = outfilename + "bands.root"
        canvas.Print(os.path.join(outdir,root_outfilename))

        
    def plot_phase_time(self, which=0, xdim=600, ydim=800, outdir="."):
        '''Plot phases as a function of the time.
        
        ___arguments___:
        which      : histogram number
        xdim, ydim : plot dimensions
        outdir     : output directory
        '''

        root.initialization()
        phaseogram2D = self.phaseogram2D
        
        num = which
        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        text = TText()
        pad1 = TPad("pad1","pad1",0,0,1,1); pad1.Draw()        

        pad1.cd()
        pad1.SetTopMargin(0.05); pad1.SetBottomMargin(0.1)
        pad1.SetLeftMargin(0.15); pad1.SetRightMargin(0.15)
        root.set_axis_thx(phaseogram2D[num],"Pulsar Phase",0.04,0.94,0.03,"Time [MET]",0.05,1.1,0.02)
        phaseogram2D[num].Draw("CONT4 Z")
        
        # =======> OUTPUT
        print "\t(emin,emax): (", self.energy_range[num][0], ",", self.energy_range[num][1], ") MeV"
        outfilename = self.psrname + "_phaseogram2D.eps"
        canvas.Print(os.path.join(outdir,outfilename))

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
        if htest.GetHistogram().GetMaximum() > 5: 
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
        
        filename = os.path.join(outdir,self.psrname+"_phaseogram.asc")
        outfile = open(filename,"w")
        outfile.write("# ------------------------------------------------------------------ \n")
        outfile.write("# Private results: LAT Collaboration and PTC ONLY \n")
        outfile.write("# Gamma-Ray Phase Histogram for %s \n" %self.psrname)
        outfile.write("# column: number of counts per energy band [MeV] \n")
        outfile.write("#%.0f-%.0f   %.0f-%.0f   %.0f-%.0f   %.0f-%.0f   %.0f-%.0f\n"
                      %(self.energy_range[0][0],self.energy_range[0][1],
                        self.energy_range[1][0],self.energy_range[1][1],
                        self.energy_range[2][0],self.energy_range[2][1],
                        self.energy_range[3][0],self.energy_range[3][1],
                        self.energy_range[4][0],self.energy_range[4][1])
                      )
        outfile.write("# ------------------------------------------------------------------ \n")

        for i in range(self.nbins):
            for histo in self.phaseogram:
                outfile.write("%.0f\t" %histo.GetBinContent(i+1))
            outfile.write("\n")

        print "INFO:", filename, "has been created."            

