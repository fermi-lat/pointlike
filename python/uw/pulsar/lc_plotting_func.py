'''
Module plotting pulse phase from a FT1 file.
'''
__author__  = 'Damien Parent' 
__version__ = 1.2

import os, sys
import numpy, pyfits
import itertools
from ROOT import TH1F, TH2F, TGraph, TText, TCanvas, TPad, TLegend, TLine

from uw.utilities.coords import sepangle_deg
import uw.utilities.fits_utils as utilfits
import uw.utilities.root_utils as root
from uw.utilities.phasetools import h_statistic, z2m, sig2sigma


def set_psfcut(par,energy):
    '''Apply a psf cut of the form (ROI x (e/100)**(alpha))**2 + beta**2'''
    theta = numpy.sqrt( numpy.power(par[0],2.) * numpy.power(energy/100.,par[1]*2.) + numpy.power(par[2],2.) )
    #theta = par[0] * numpy.power(energy/1000.,par[1])
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
        self.eclsmax       = 3
        self.zenithcut     = 105.
        self.phase_colname = 'PULSE_PHASE'
        self.psfcut        = True

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
        radmin        : Keep all the photons in x (deg). Default: 0.
        emin          : Minimum energy (MeV).            Default: FT1 header.
        emax          : Maximum energy (MeV).            Default: FT1 header.
        tmin          : Start time (MET).                Default: FT1 header.
        tmax          : End time (MET).                  Default: FT1 header.
        eclsmin       : Minimum event class              Default: 0
        eclsmax       : Maximum event class              Default: 3
        zenithcut     : Zenith Cut (deg)                 Default: 105.
        phase_colname : Phase column name.               Default: PULSE_PHASE
        psfcut        : PSF-cut selection [True/False]   Default: True 
        '''

        self.init()
        self.__dict__.update(kwargs)
        self.__ra_src  = ra_src
        self.__dec_src = dec_src
        
        self.psf_selection = None
        self.energy_range  = [ [100.,3e5] , [100.,300.] , [300.,1000.] , [1000.,3e5] , [3000.,3e5] ]

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

        # ===== declaration of histograms (ROOT ...) =====
        self.nbins = 0
        self.phaseogram = []
        self.pulse_phase = []
        self.phaseogram2D = []
        for i in range(5):
            histname = "phaseogram_" + str(i)
            self.phaseogram.append(TH1F())
            self.pulse_phase.append(numpy.array([]))
            self.phaseogram2D.append(TH2F())
                        
        # ===== Get the events from the FT1 file =====
        print "\n___pending___: loading photons from FT1 file ..."
        ft1file = pyfits.open(ft1name)
        # ft1hdr  = hdulist[1].header
        # ft1dat  = hdulist[1].data
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

    def fill_phaseogram(self, nbins=20, bin_range=[0.,2.], convtype=-1, phase_shift=0.):
        '''Fill both phaseograms (TH1F, TH2F) and pulse_phase objects.
        
        parameters:
        -----------
        nbins       : number of bins in your histogram
        bin_range   : Min and Max phase interval
        convtype    : -1=both, 0=Front, 1=Back
        phase_shift : add a phase shift in your histogram
        '''

        print "\n___pending___: loading histograms ..."
        
        self.nbins = nbins
        self.phase_shift = phase_shift

        # === initialization of histograms ===
        for item in self.phaseogram:
            item.Reset()
            item.SetBins(nbins*int(bin_range[1]-bin_range[0]),bin_range[0],bin_range[1])

        for item in self.phaseogram2D:
            item.Reset()
            item.SetBins(nbins,0.,1.,int((self.tmax-self.tmin)/(14*86400)),self.tmin,self.tmax)

        for i, item in enumerate(self.pulse_phase):
            self.pulse_phase[i] = numpy.array([])
            
        if self.psf_selection is None:
            if   convtype == -1: psf_selection = [ 5.12 , -0.8 , 0.07 ]
            elif convtype ==  0: psf_selection = [ 3.3 , -0.8 , 0.05 ]
            elif convtype ==  1: psf_selection = [ 5.5 , -0.8 , 0.1 ]
            else: print "convtype that you selected does not exist. Exiting ..."
            self.psf_selection = psf_selection
            
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
            histo.SetMaximum(max_histo*factor)

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

    def print_summary(self):
        print "___PARAMETERS___:"
        print "\t(tmin,tmax): (", self.tmin, ",", self.tmax, ") MET"
        print "\t(emin,emax): (", self.emin, ",", self.emax, ") MeV"
        print "\t(ra,dec): (", self.__ra_src, ",", self.__dec_src, ") deg"        
        print "\t(radmin,radmax): (", self.radmin, ",", self.radius, ") deg"        
        if self.psfcut:
            sys.stdout.write("\tPSF-cut: "),
            print self.psf_selection[0],
            sys.stdout.write(" x (E/100)^"),
            print self.psf_selection[1],
            sys.stdout.write(" ++ "),
            print self.psf_selection[2],
            print "(++: in quadrature) \n"

        print "PARAMETERS | \t PHASO_0 \t PHASO_1 \t PHASO_2 \t PHASO_3 \t PHASO_4"
        for x in range(30): print "- ",
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
        
    # ================================= PLOT =================================
    # The following functions are ...
    # ========================================================================
    def plot_lc_gamma_radio(self, nbands=3, xdim=1200, ydim=1700, outdir=".", profile=None, ncol=0,
                            fidpt=0., ylabel="Radio Flux (au)", comment="1.4 GHz", mean_sub=False, background=None):
        '''ROOT function for plotting [3,4,5] gamma-ray bands and 1 other external band.
        
        ___arguments___:
        nbands     : number of gamma-ray bands in the figure [3,4,5]
        xdim, ydim : plot dimensions
        outdir     : output directory (Default=.)
        profile    : radio or x-ray template file (ascii)
        ncol       : column number for the radio template file
        fidpt      : fiducial point
        ylabel     : ylabel title for the bottom panel
        comment    : add a comment about the radio observed frequency
        '''
        
        root.initialization()
        phaseogram = self.phaseogram

        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        text = TText()

        if nbands is 3:
            text.SetTextSize(0.08)
            pad1 = TPad("pad1","pad1",0,0.63,1.,1.00); pad1.Draw()
            pad3 = TPad("pad3","pad3",0,0.26,1.,0.62); pad3.Draw()
            pad5 = TPad("pad5","pad5",0,0.00,1.,0.25); pad5.Draw()

        if nbands is 4:
            text.SetTextSize(0.09)
            pad1 = TPad("pad1","pad1",0,0.69,1.,1.00); pad1.Draw()
            pad2 = TPad("pad2","pad2",0,0.47,1.,0.68); pad2.Draw()
            pad3 = TPad("pad3","pad3",0,0.26,1.,0.47); pad3.Draw()
            pad5 = TPad("pad5","pad5",0,0.00,1.,0.25); pad5.Draw()
                                            
        if nbands is 5:
            text.SetTextSize(0.09)
            pad1 = TPad("pad1","pad1",0,0.70,1.,1.00); pad1.Draw()
            pad2 = TPad("pad2","pad2",0,0.51,1.,0.69); pad2.Draw()
            pad3 = TPad("pad3","pad3",0,0.33,1.,0.51); pad3.Draw()
            pad4 = TPad("pad4","pad4",0,0.15,1.,0.33); pad4.Draw()
            pad5 = TPad("pad5","pad5",0,0.00,1.,0.14); pad5.Draw()
        
        # ===========> TOP PANEL
        pad1.cd()
        pad1.SetTopMargin(0.02); pad1.SetBottomMargin(0.)
        
        ecol = 0
        if nbands == 3: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.08,0.6,0.05)
        if nbands == 4: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.07,0.65,0.055)
        if nbands == 5: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.053,0.84,0.047)
        ### phaseogram[ecol].SetFillColor(14) #######
        phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
        tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.88,tt)
        text.DrawText(1.45,phaseogram[ecol].GetMaximum()*0.88,self.psrname)

        if background is not None:
              line = TLine(0,background,2,background)
              line.SetLineStyle(2)
              line.Draw()
        
        # ===========> (2) PANEL
        if nbands is not 3:
            pad2.cd()
            pad2.SetTopMargin(0.); pad2.SetBottomMargin(0.)
            
            ecol = 3
            if nbands is 4: root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.11,0.42,0.08)
            if nbands is 5: root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.09,0.5,0.08)
            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            text.SetTextSize(0.12)
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.88,tt)
            
            ecol = 4
            phaseogram[ecol].SetFillColor(14)
            phaseogram[ecol].Draw("same")
            text.SetTextColor(14)
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[3].GetMaximum()*0.78,tt)
            
        # ============> (3) PANEL
        pad3.cd()
        pad3.SetTopMargin(0.); pad3.SetBottomMargin(0.)
        
        ecol = 2
        if nbands is 3: root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.08,0.6,0.05)
        if nbands is 4: root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.11,0.42,0.08)
        if nbands is 5: root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.09,0.5,0.08)
        phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
        text.SetTextColor(1)
        tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.88,tt)
        
        # ============> (4) PANEL
        if nbands is 5:
            pad4.cd()
            pad4.SetTopMargin(0.); pad4.SetBottomMargin(0.)
            
            ecol = 1
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.09,0.5,0.08)
            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            text.SetTextColor(1)
            tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.88,tt)

        # ==========> TEMPLATE - BOTTOM PANEL
        pad5.cd()
        pad5.SetBorderSize(0); pad5.SetTopMargin(0.); pad5.SetBottomMargin(0.3)
        
        if profile is not None:
            phase_ext = numpy.loadtxt(profile,comments="#")[:,ncol-1]
            template = TGraph(2*len(phase_ext))
            phase_ext = numpy.roll(phase_ext,-int(len(phase_ext)*fidpt))  # fiducial point
            if mean_sub is True:
                phase_ext -= numpy.mean(phase_ext)
            phase_ext /= numpy.max(phase_ext)   # normalization

            for i, item in enumerate(phase_ext):
                template.SetPoint(i,float(i)/len(phase_ext),item)
                template.SetPoint(i+len(phase_ext),1.+float(i)/len(phase_ext),item)
                
            # title size, title offset, label size
            if nbands is 3: root.set_axis_thx(template,"Pulsar Phase",0.11,1.0,0.11,ylabel,0.1,0.44,0.09)
            if nbands is 4: root.set_axis_thx(template,"Pulsar Phase",0.11,1.0,0.11,ylabel,0.1,0.44,0.09)
            if nbands is 5: root.set_axis_thx(template,"Pulsar Phase",0.12,1.0,0.11,ylabel,0.1,0.44,0.12)

            template.GetXaxis().SetRangeUser(0,2)
            template.Draw("al");
            text.SetTextSize(.14);
            text.DrawText(0.1,0.8,comment)
            
        # =======> OUTPUT
        outfilename = self.psrname + "_phaseogram_gamma_radio_" + str(nbands) + "bands.eps"
        canvas.Print(os.path.join(outdir,outfilename))

    def plot_lc_gamma(self, nbands=4, xdim=1200, ydim=1500, outdir=".", background=None):
        '''ROOT function for ploting x gamma-ray bands.
            
        ___arguments___:
        nbands     : number of gamma-ray bands in the figure [2,..,5]
        xdim, ydim : plot dimensions
        outdir     : output directory
        '''
        
        root.initialization(batch=True)
        phaseogram = self.phaseogram
        
        # ==========> canvas and pad
        canvas = TCanvas("canvas","canvas",xdim,ydim);
        text = TText()
        
        if nbands is 3:
            text.SetTextSize(0.08)
            pad1 = TPad("pad1","pad1",0,0.63,1.,1.00); pad1.Draw()
            pad3 = TPad("pad3","pad3",0,0.34,1.,0.62); pad3.Draw()
            pad4 = TPad("pad4","pad4",0,0.00,1.,0.34); pad4.Draw()
        
        if nbands is 4:
            text.SetTextSize(0.09)
            pad1 = TPad("pad1","pad1",0,0.67,1.,1.00); pad1.Draw()
            pad3 = TPad("pad3","pad3",0,0.46,1.,0.66); pad3.Draw()
            pad4 = TPad("pad4","pad4",0,0.26,1.,0.46); pad4.Draw()
            pad5 = TPad("pad5","pad5",0,0.00,1.,0.26); pad5.Draw()
            
        if nbands is 5:
            text.SetTextSize(0.09)
            pad1 = TPad("pad1","pad1",0,0.71,1.,1.00); pad1.Draw()
            pad2 = TPad("pad2","pad2",0,0.54,1.,0.70); pad2.Draw()  
            pad3 = TPad("pad3","pad3",0,0.38,1.,0.54); pad3.Draw()
            pad4 = TPad("pad4","pad4",0,0.22,1.,0.38); pad4.Draw()
            pad5 = TPad("pad5","pad5",0,0.00,1.,0.22); pad5.Draw()
            
        # ===========> TOP PANEL
        pad1.cd()
        pad1.SetTopMargin(0.05); pad1.SetBottomMargin(0.)
        
        ecol = 0
        if nbands == 3: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.07,0.65,0.05)
        if nbands == 4: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.067,0.67,0.047)
        if nbands == 5: root.set_axis_thx(phaseogram[ecol],"",0.,0.,0.,"Counts",0.065,0.74,0.05)
        phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
        tt = "> %.2f GeV" %(self.energy_range[ecol][0]/1e3)
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)
        text.DrawText(1.45,phaseogram[ecol].GetMaximum()*0.87,self.psrname)

        if background is not None:
            line = TLine(0,background,2,background)
            line.SetLineStyle(2)
            line.Draw()

        # ===========> (2) PANEL
        if nbands is 5:
            pad2.cd(); pad2.SetTopMargin(0.); pad2.SetBottomMargin(0.)
            
            ecol = 4
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.12,0.39,0.08)
            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            text.SetTextSize(0.13)
            phaseogram[ecol].SetFillColor(0)
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)
            
        # ===========> (3) PANEL
        pad3.cd(); pad3.SetTopMargin(0.); pad3.SetBottomMargin(0.)
                                    
        ecol = 3
        if nbands is 3:
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.10,0.45,0.07)
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            textsize = 0.10
        if nbands is 4:
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.11,0.42,0.08)
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            textsize = 0.12
        if nbands is 5:
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.12,0.39,0.08)
            tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
            textsize = 0.13

        phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
        text.SetTextSize(textsize)
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.87,tt)

        if nbands is not 5:
            ecol = 4
            phaseogram[ecol].SetFillColor(14)
            phaseogram[ecol].Draw("same")
            text.SetTextColor(14)
            tt = "> " + str(self.energy_range[ecol][0]/1e3) + " GeV"
            text.DrawText(0.1,phaseogram[3].GetMaximum()*0.78,tt)
        
        # ===========> (4) PANEL
        pad4.cd()
        pad4.SetTopMargin(0.); pad4.SetBottomMargin(0.)
            
        ecol = 2
        if nbands is 3:
            pad4.SetBottomMargin(0.2)
            root.set_axis_thx(phaseogram[ecol],"Pulsar Phase",0.085,0.9,0.065,"Counts",0.085,0.55,0.055)
            textsize = 0.08
        if nbands is 4:
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.11,0.42,0.08)
            textsize = 0.13
        if nbands is 5:
            root.set_axis_thx(phaseogram[ecol],"",0.0,0.0,0.0,"Counts",0.12,0.39,0.08)
            textsize = 0.13
        
        phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
        text.SetTextColor(1); text.SetTextSize(textsize)
        tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
        text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.86,tt)

        # ===========> (5) BOTTOM PANEL
        if nbands is not 3:
            pad5.cd()
            pad5.SetTopMargin(0.); pad5.SetBottomMargin(0.2)
            
            ecol = 1
            if nbands is 4:
                root.set_axis_thx(phaseogram[ecol],"Pulsar Phase",0.09,1.0,0.07,"Counts",0.086,0.55,0.06)
            if nbands is 5:
                pad5.SetBottomMargin(0.3)
                root.set_axis_thx(phaseogram[ecol],"Pulsar Phase",0.11,1.1,0.08,"Counts",0.086,0.55,0.06)
                
            phaseogram[ecol].Draw(); phaseogram[ecol].Draw("Esame")
            tt = str(self.energy_range[ecol][0]/1e3) + " - " + str(self.energy_range[ecol][1]/1e3) + " GeV"
            text.SetTextSize(0.1)
            text.DrawText(0.1,phaseogram[ecol].GetMaximum()*0.85,tt)
                                                                                                                        
        # =======> OUTPUT
        outfilename = self.psrname + "_phaseogram_gamma_" + str(nbands) + "bands.eps"
        canvas.Print(os.path.join(outdir,outfilename))

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
        
        filename = os.path.join(output,self.psrname+"_phaseogram.asc")
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

