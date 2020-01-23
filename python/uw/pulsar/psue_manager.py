"""
   name: psue_manager.py
   author(s): Damien Parent <dmnparent@gmail.com>
   creation Date: July 2011

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/psue_manager.py,v 1.10 2012/11/29 00:30:37 kerrm Exp $
"""

__version__ = 1.9

from os import system, remove, F_OK, mkdir, access, environ, chdir, linesep, getcwd
from os.path import join, basename, abspath, isfile
from sys import exit
import cPickle, glob
from shutil import copy, copyfile, move
import time, random
import numpy as np
import astropy.io.fits as pyfits
from GtApp import GtApp

import matplotlib
matplotlib.use('Agg')    # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as pl
from uw.utilities.xml_utils import xml_manager
import uw.pulsar.lc_plotting_func as plot
from uw.pulsar.parfiles import ParFile
from uw.pulsar.stats import sf_h20_dj2010 as h_sig, hmw, sigma_trials, sig2sigma,h2sig

# this function should be moved into fits_utils.py
def SelectPhase(infile,outfile=None,phmin=0,phmax=1,phcolname='PULSE_PHASE'):
    '''Select a phase range into a FT1 file'''

    ft1file = pyfits.open(infile)    
    try:
        ft1file['EVENTS'].data.field(phcolname)

        if isinstance(phmin,float) and isinstance(phmax,float):
            phname = '_%s_%sph.fits' %(phmin,phmax)
            phfactor = phmax-phmin
            expr = 'PULSE_PHASE>%.2f && PULSE_PHASE<%.2f' %(phmin,phmax)
            if phmin>phmax:
                expr = expr.replace('&&','||')
                phfactor = (1-phmin)+phmax
            
        elif isinstance(phmin,list) and isinstance(phmax,list):
            phname = '_%s_%s_%s_%s_phase.fits' %(phmin[0],phmin[1],phmax[0],phmax[1])
            if phmax[0]<phmax[1]:
                expr = '(PULSE_PHASE>%.2f && PULSE_PHASE<%.2f) || (PULSE_PHASE>%.2f && PULSE_PHASE<%.2f)' %(phmin[0],phmin[1],phmax[0],phmax[1])
                phfactor = (phmin[1]-phmin[0]) + (phmax[1]-phmax[0])
            else:
                expr = 'PULSE_PHASE<%.2f || (PULSE_PHASE>%.2f && PULSE_PHASE<%.2f) || PULSE_PHASE>%.2f' %(phmax[1],phmin[0],phmin[1],phmax[0])
                phfactor = phmax[1] + (phmin[1]-phmin[0]) + (1.-phmax[0])

        outfile = infile.replace('.fits',phname) if outfile is None else outfile
        fselect = GtApp('fselect')
        fselect.run(infile=infile,outfile=outfile,expr=expr,clobber='yes')
        return outfile, phfactor
        
    except KeyError:
        from sys import exit
        print ("No %s column. Exiting ..." %phcolname); exit()

# =============================
# Class to manage the PSUE file
# =============================
class PSUEFile(dict):

    def _switch_lc_format(self):
        """ Temporary function wihle adding support for asymmetric peak widths. """
        for key,postfix in zip(['PEAK1','PEAK2','BRIDGE'],['P1','P2','B']):
            for ik,k in enumerate(self.ordered_keys):
                if k==key:
                    if (ik+2 < len(self.ordered_keys)) and (self.ordered_keys[ik+2] == 'FWHM_%s'%postfix):
                        old_key = self.ordered_keys.pop(ik+2)
                        old_val = self.pop(old_key)
                        self.ordered_keys.insert(ik+2,'FWHM_R_%s'%postfix)
                        self.ordered_keys.insert(ik+2,'FWHM_L_%s'%postfix)
                        self['FWHM_R_%s'%postfix] = old_val
                        self['FWHM_L_%s'%postfix] = old_val
        # what a kluge....
        while ('FWHM' in self.ordered_keys[-2]):
            self.ordered_keys.pop(-2)
        while ('FWHM' in self.ordered_keys[-1]):
            self.ordered_keys.pop()

    def __init__(self,psuefile):
        self.ordered_keys = []
        if isfile(psuefile):
            for line in file(psuefile):
                tok = line.strip().split()
                if len(tok)==0: continue
                self.ordered_keys.append(tok[0])
                self[tok[0]] = tok[1:] if (len(tok[1:]) > 1) else tok[1:][0]
            # NB -- to remove ASAP
            self._switch_lc_format()

    def get(self,key,first_elem=True,type=str):
        try:
            if type==float:
                type = lambda x: float(x.replace('D','E'))
            t = self[key]
            if hasattr(t,'__iter__'):
                if first_elem: return type(t[0])
                return [type(x) for x in t]
            return type(t)
        except KeyError:
            return None

    def set(self,key,val):
        if isinstance(val,list) or isinstance(val,tuple):
            val = [str(x) for x in val]
        if key not in self.keys():
            self.ordered_keys.append(key)
        self[key] = val

    def write(self,output):
        f = file(output,'w')
        for key in self.ordered_keys:
            val = self[key]
            if hasattr(val,'__iter__'):
                val = '  '.join(val)
            key += ' '*(18-len(key))
            f.write('%s%s\n'%(key,val))
        f.close()

# ==================================================
# Object
# ==================================================
class PSUEAnalysis():

    def init(self):
        self.data_version     = 'P7.6_P120_BASE'
        self.radius           = 10
        self.emin, self.emax  = 1e2, 3e5                 # MeV
        self.tmin, self.tmax  = 239557414, 334108802     # MET
        self.zmax             = 100
        self.lc_radius        = 2.
        self.fit_emin         = 1e2
        self.fit_emax         = 1e5
        self.w_radius        = 2.
        self.w_emin         = 1e2
        self.w_emax         = 1e5
        self.evtclass         = 'Source'
        self.phase_colname    = 'PULSE_PHASE'
        self.offpulse         = []
        self.irfs             = 'P7SOURCE_V6'
        self.psfpar           = [5.3,0.745,0.09]
        self.ft2file          = '/afs/slac/g/glast/groups/pulsar/data/WeeklyData/ft2/FT2.fits'
        self.xmlcat           = '/afs/slac/g/glast/groups/pulsar/share/catalog/gll_psc_v05.xml'
        self.galdiff          = '/afs/slac/g/glast/groups/diffuse/rings/2year/ring_2year_P76_v0.fits'
        self.isotropic        = '/afs/slac/g/glast/groups/diffuse/rings/2year/isotrop_2year_P76_source_v0.txt'
        self.limb_spec        = '/afs/slac/g/glast/groups/diffuse/rings/2year/smooth_limb/limb_2year_P76_source_v0_smooth.txt'
        self.limb_spatial     = '/afs/slac/g/glast/groups/diffuse/rings/2year/smooth_limb/limb_smooth.fits'
        self.template_dir     = '/afs/slac/g/glast/groups/pulsar/share/Templates'
        self.outdir           = None
        self.batch            = False
        self.pulsar_class     = 0

        assert(self.radius>=self.lc_radius)
        
    def __init__( self, parfile=None, **kwargs ):
        '''
        ============    ===============================================
        keyword         description
        ============    ===============================================
        data_version    LAT data version                 [P7.6_P120_BASE]
        radius          Maximum radius (deg)             [10]
        emin            FT1 minimum energy (MeV)         [100]
        emax            FT1 maximum energy (MeV)         [300000]
        tmin            Start time (MET)                 [239557514]
        tmax            End time (MET)                   [334108802]
        zmax            Maximum zenith angle value (deg) [100]
        lc_radius       Radius for pulse profile (deg)   [2]
        fit_emin        Minimum energy for spec (MeV)    [100]
        fit_emax        Maximum energy for spec (MeV)    [100000]
        w_emin          Minimum energy for weights (MeV) [100]
        w_emax          Maximum energy for weights (MeV) [100000]
        w_radius        Maximum radius for weights (deg) [2]
        evtclass        Event class selection            [SOURCE]
        phase_colname   Name of the phase column         [PULSE_PHASE]
        irfs            IRFs version                     [P7SOURCE_V6]
        psfpar          PSF parameters                   [5.3,0.745,0.09]
        pulsar_class    0=unid,1=RL,2=RL-MSP,3=RQ        [0]
        '''        

        self.init()
        self.__dict__.update(kwargs)

        # load parfile
        if isinstance(parfile,str): parfile = [parfile]
        parfile = [abspath(p) for p in parfile]
        self.parfile = parfile

        # get pulsar parameters
        ephem = ParFile(self.parfile[0])
        self.psrname = ephem.get_psrname()
        self.ra, self.dec = ephem.get_ra(), ephem.get_dec()
        self.period = ephem.p()

        # outdir
        if self.outdir is None: self.outdir = abspath(join(getcwd(),self.psrname))            
        self.outdir_local = self.outdir
        if not access(self.outdir_local,F_OK): mkdir(self.outdir_local)

        # ft1 and ft2 init and more
        self.ft1file         = join(self.outdir,'%s_%.0fDeg_%.0f_%.0fMET_%.0f_%.0fMeV_z%.0f_%s.fits') \
                               %(self.psrname,self.radius,self.tmin,self.tmax,self.emin,self.emax,self.zmax,self.evtclass)
        self.ft1file_gtis    = join(self.outdir,'%s_%.0fDeg_%.0f_%.0fMET_%.0f_%.0fMeV_z%.0f_%s_gtis.fits') \
                               %(self.psrname,self.radius,self.tmin,self.tmax,self.emin,self.emax,self.zmax,self.evtclass)
        self.ft1file_diffrsp = join(self.outdir,'%s_%.0fDeg_%.0f_%.0fMET_%.0f_%.0fMeV_z%.0f_%s_gtis_diffrsp.fits') \
                               %(self.psrname,self.w_radius,self.tmin,self.tmax,self.w_emin,self.w_emax,self.zmax,self.evtclass)
        self.ft1file_weights = join(self.outdir,'%s_%.0fDeg_%.0f_%.0fMET_%.0f_%.0fMeV_z%.0f_%s_gtis_weights.fits') \
                               %(self.psrname,self.w_radius,self.tmin,self.tmax,self.w_emin,self.w_emax,self.zmax,self.evtclass)
        self.gtltcube        = self.ft1file_gtis.replace('.fits','_gtltcube.fits')
        self.ft2file         = abspath(self.ft2file)
        self.srcmdl_in       = self.ft1file_gtis.replace('.fits','_srcmdl_input.xml')
        self.srcmdl_out      = self.ft1file_gtis.replace('.fits','_srcmdl_output.xml')

        if not isfile(self.ft2file):
            print ("Error. Cannot open %s! Exiting ..." %self.ft2file; exit())
        
        # BATCH mode
        if self.batch:
            self.outdir = '/scratch/psue' + str(random.randint(1,99999))
            try: mkdir(self.outdir)
            except OSError: raise Exception("You are trying the BATCH mode. Run with bsub. Exiting...")

            # setup pfiles
            gtselect  = join(environ['INST_DIR'],'syspfiles','gtselect.par'); copy(gtselect,self.outdir)
            gtmktime  = join(environ['INST_DIR'],'syspfiles','gtmktime.par'); copy(gtmktime,self.outdir)
            gtsrcprob = join(environ['INST_DIR'],'syspfiles','gtsrcprob.par'); copy(gtsrcprob,self.outdir) 
            gtdiffrsp = join(environ['INST_DIR'],'syspfiles','gtdiffrsp.par'); copy(gtdiffrsp,self.outdir)
            fselect   = join(environ['HEADAS'],'syspfiles','fselect.par'); copy(fselect,self.outdir)
            environ['PFILES'] = self.outdir
            chdir(self.outdir) # for tempo2            
                
        # PSUE output
        filename = join(self.outdir_local,'%s_psue.txt' %(self.psrname))
        self.PSUEoutfile = PSUEFile(filename)
        self.PSUEoutfile.set('PSRNAME',self.psrname)
        self.PSUEoutfile.set('RAJ',self.ra)
        self.PSUEoutfile.set('DECJ',self.dec)
        self.PSUEoutfile.set('TMIN',[self.tmin,'MET'])
        self.PSUEoutfile.set('TMAX',[self.tmax,'MET'])
        self.PSUEoutfile.set('PULSAR_CLASS',self.pulsar_class)
        [self.PSUEoutfile.set('PARFILE'+str(i),basename(p)) for i, p in enumerate(self.parfile)]
        
    def WritePSUE(self):
        self.PSUEoutfile.write(join(self.outdir_local,self.psrname+'_psue.txt'))

    def CleanFiles(self):
        '''Clean scratch directory'''
        from shutil import rmtree
        if self.batch: rmtree(self.outdir)        

    def CreateFT1file(self,fullft1file=None):
        '''Download a FT1 file from the astroserver or
           Select a region from an existing FT1 file [fullft1file].'''
        astro=False
        if not isfile(self.ft1file):
            ft1file = join(self.outdir,basename(self.ft1file))

            if fullft1file is not None and isfile(fullft1file):
                infile = fullft1file
            else:
                # download from the astro server at SLAC
                astro=True; infile = join(self.outdir,self.psrname+'.fits')
                system('~glast/astroserver/prod/astro -b -q --output-ft1 %s --event-sample %s --minEnergy %.2f --maxEnergy %.2f \
                --minTimestamp %.2f --maxTimestamp %.2f --ra %.6f --dec %.6f --radius %.2f --maxZenith %.2f --event-class-name "%s" store' \
                       %(infile,self.data_version,self.emin,self.emax,self.tmin,self.tmax,self.ra,self.dec,self.radius,self.zmax,self.evtclass))

            GtApp('gtselect').run(infile=infile,outfile=ft1file,rad=self.radius,ra=self.ra,dec=self.dec,
                                  emin=self.emin,emax=self.emax,tmin=self.tmin,tmax=self.tmax,zmax=self.zmax)
            
            if astro: remove(infile)  # clean astroft1file
            if self.batch: copy(ft1file,self.outdir_local) # copy the ft1 file to the local directory
                
    def FoldEvents(self,forceFolding=False):
        '''Fold gamma-rays using the TEMPO2 plugin
        ============    =========================================
        keyword         description
        ============    ========================================= 
        forceFolding    If True, the phases will be replaced '''

        if not isfile(self.ft1file):
            print ("Error. Cannot open %s! Must create FT1file. Exiting ..." %self.ft1file; exit())
            
        phase_exists = False
        file = pyfits.open(self.ft1file)
        for c in file['EVENTS'].get_coldefs():
            if c.name == self.phase_colname: phase_exists = True

        if not phase_exists or forceFolding:

            ft1file = join(self.outdir,basename(self.ft1file))
            ft2file = self.ft2file

            if self.batch:
                if not isfile(ft1file):
                    print ('Copying %s to %s'%(self.ft1file,self.outdir))
                    copy(self.ft1file,self.outdir)
                if not isfile(ft2file): copy(self.ft2file,self.outdir)
                
            # run TEMPO2
            for p in self.parfile:
                start, finish = ParFile(p).get("START",type=float), ParFile(p).get("FINISH",type=float)
                # remove line below for J2215 -- what other pulsars does it affect?
                #if len(self.parfile) == 1: start, finish = 54000, 58000
                cmd = "tempo2 -gr fermi -ft1 %s -ft2 %s -f %s -tmin %s -tmax %s -phase -graph 0" %(ft1file,ft2file,p,start,finish)
                print ('Issuing the folding command: \n%s'%cmd)
                system(cmd)

            # copy the ft1 file to the local directory
            if self.batch: 
                print ('Copying %s to %s'%(ft1file,self.outdir_local))
                copy(ft1file,self.outdir_local)

            # regenerate gtmktime-processed file
            ft1file_gtis = join(self.outdir,basename(self.ft1file_gtis))
            #if isfile(ft1file_gtis):
            filter = "DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52"
            GtApp('gtmktime').run(evfile=ft1file,outfile=ft1file_gtis,scfile=self.ft2file,filter=filter,roicut='no',chatter=4)
            if self.batch: 
                print ('Copying %s to %s'%(ft1file,self.outdir_local))
                copy(ft1file_gtis,self.outdir_local)

    def SearchPulsation(self):
        '''Run SearchPulsation (Philippe macro)'''

        ft1file = self.ft1file
        if not isfile(ft1file): print ("Cannot open %s!" %ft1file)
            
        self.PSUEoutfile.set('#SEARCH_PULSATION','................')        
        spfile = join(self.outdir_local,"searchpulsation_%s_%s.txt"%(self.psrname,self.evtclass))
        system("SearchPulsation -b fitsfile=%s psrtitle='%s' ra=%.6f dec=%.6f nbint=20 zenithcut=%.1f ephemfile=%s outfile=%s psfpar0=%.3f psfpar1=%.3f psfpar2=%.3f radialdistmax=5" \
               %(ft1file,self.psrname,self.ra,self.dec,self.zmax,self.parfile[0],spfile.split('.txt')[0],self.psfpar[0],self.psfpar[1],self.psfpar[2]) )

        for p in self.GetSearchPulsationParams(spfile):
            if p['name'] == 'sig_1D': val = float(p['value']); self.PSUEoutfile.set('1D',['%.1f'%val,'sig'])
            if p['name'] == 'sig_2D': val = float(p['value']); self.PSUEoutfile.set('2D',['%.1f'%val,'sig'])
            if p['name'] == 'sig_W' : val = float(p['value']); self.PSUEoutfile.set('W',['%.1f'%val,'sig'])
            if p['name'] == 'sig_allW' : val = float(p['value']); self.PSUEoutfile.set('allW',['%.1f'%val,'sig'])
            
    def GetSearchPulsationParams(self, filename):
        names, params = ['method','prob','sig','tref','loge','logr','alpha'], []
        if isfile(filename):
            for line in open(filename,'r'):
                tok = line.strip().split()
                if tok[0] == '#': continue
                else: params += [{'name':n+"_"+tok[0],'value':val} for (n,val) in zip(names,tok)]
            return params

    def TimingAnalysis( self, nbins=None, weight=False, offmin=0, offmax=1, background=None, 
                        profile=None, peakShape=None, peakPos=None, 
                        use_spw=False, dm_unc=0, phaseShift=0):
        '''Generate gamma-ray pulse profiles.
        ============    =============================================
        keyword         description
        ============    ============================================= 
        weight          if True, plot weighted counts       [False]
        nbins           Number of bins for phaso            [None]
        profile         Name of the radio template          [None]
        use_spw         Use Weights from SearchPulsation    [False]
        dm_unc          Uncertainty on dispersion (DM)      [0]
        phaseShift      A common shift to both radio/gamma  [0]
        '''

        if weight: ft1file = self.ft1file_weights
        else: ft1file = self.ft1file
        if not isfile(ft1file): print ("Cannot open %s!" %ft1file)
        
        self.PSUEoutfile.set("#PULSE_PROFILE",'................')
        self.PSUEoutfile.set('OFF_RANGE',['%.2f'%offmin,'%.2f'%offmax])
    
        # ============ PulsarLightCurve Object =========
        #   initiatilize the PulsarLightCurve object
        # ==============================================
        plc = plot.PulsarLightCurve(ft1file,psrname=self.psrname,radius=self.lc_radius,emin=self.emin,emax=self.emax)

        # default energy bands
        plc.set_energy_range(0,1e2,1e5)
        plc.set_energy_range(1,3e3,1e5)
        plc.set_energy_range(2,1e3,3e3)
        plc.set_energy_range(3,3e2,1e3)
        plc.set_energy_range(4,1e2,3e2)
        plc.set_energy_range(-1,1e4,1e5)

        if not weight:  # counts 
            plc.psfcut=True
            for i in range(plc.nhisto):
                emin, emax = plc.get_energy_range(which=i)
                ntrials, presig, postsig, bestE, bestR = plc.pulseopt(erange=[emin,emin],emax=emax,ebins=1,rrange=[0.2,self.lc_radius],rbins=40)
                plc.set_radius_range(which=i,radius=bestR)
            pulse_profile = ft1file.replace('.fits','_pulse_profile_c.eps')
            phaseogram2d = ft1file.replace('.fits','_phaseogram2D_c.eps')

        if weight:  # weighted counts from spectral analysis
            [plc.set_radius_range(which=i,radius=self.lc_radius) for i in range(plc.nhisto)]
            pulse_profile = ft1file.replace('.fits','_pulse_profile_wc.eps')
            phaseogram2d = ft1file.replace('.fits','_phaseogram2D_wc.eps')
            
        if weight and use_spw:  # weighted counts from Philippe's method
            logeref, logesig = 2., 0.5
            infile = join(self.outdir_local,"searchpulsation_%s_%s.txt"%(self.psrname,self.evtclass))
            for p in self.GetSearchPulsationParams(infile):
                if p['name'] == 'loge_W': logeref = float(p['value'])
                if p['name'] == 'logr_W': logesig = float(p['value'])
                plc.use_searchpulsation_weights(logeref=logeref,logesig=logesig)
            pulse_profile = ft1file.replace('.fits','_pulse_profile_spwc.eps')
            phaseogram2d = ft1file.replace('.fits','_phaseogram2D_spwc.eps')

        # Fill histograms
        phase_shift = plc.fill_phaseogram(nbins=nbins,weight=weight,profile=profile,phase_shift=phaseShift)
        self.PSUEoutfile.set('PHASE_SHIFT','%.5f'%phase_shift)
        plc.print_summary()

        # Write out new parfiles with correct TZR* parameters
        for p in self.parfile:
            outp = ParFile(p)
            outp.shift_tzr(phase_shift)
            outp.write(p+'.shifted')

        # background
        if background == 'ring':
            bkg = []
            for i in range(plc.nhisto):
                radius = plc.get_radius_range(which=i)
                bkg += [plc.get_background(which=i,phase_range=[offmin,offmax],ring_range=[radius-0.5,radius+0.5])[0]]
        elif background == 'weight':
            bkg = plc.get_background(method='weight')
        else: bkg = None
        
        # off-pulse region
        if offmin==0 and offmax==1: reg=None
        else: reg=[offmin,offmax]
                        
        # phase vs time
        plc.plot_phase_time(which=0,outfile=phaseogram2d)
        plc.plot_phase_time(which=0,outfile=phaseogram2d.replace('.eps','.png'))
        
        # ascii
        plc.toASCII(outdir=self.outdir_local)

        # H-test output
        phases  = plc.get_phases(which=0)
        weights = plc.get_weights(which=0)
        H = hmw(phases,weights); sig = h2sig(H)
        if weight: self.PSUEoutfile.set('WHTEST',['%.1f'%sig,'sig'])
        else: self.PSUEoutfile.set('HTEST',['%.1f'%sig,'sig'])

        # =========== FIT PULSE PROFILE =============
        best_template = None # "forward declaration"
        if (peakShape is not None) and (peakPos is not None):

            peakPos = np.mod([x+phase_shift for x in peakPos],1)
            
            import uw.pulsar.lcprimitives as lp
            import uw.pulsar.lcfitters as lf
            import uw.pulsar.lctemplate as lt

            loglike = -np.inf
            #consts = [lp.LCGaussian,lp.LCLorentzian,lp.LCGaussian2,lp.LCLorentzian2]
            consts = [lp.LCGaussian,lp.LCGaussian2,lp.LCLorentzian,lp.LCLorentzian2]

            # use previous peak shape if available
            """
            m_name = self.PSUEoutfile.get('PEAK1')
            if m_name is not None:
                if m_name == 'Gaussian': consts = [lp.LCGaussian]
                elif m_name == 'Lorentzian': consts = [lp.LCLorentzian]
                elif m_name == 'Gaussian2': consts = [lp.LCGaussian2]
                elif m_name == 'Lorentzian2': consts = [lp.LCLorentzian2]
                else: print ("Error. Unknown model!")
            """
            # use seed values for initial Gaussian fit
            npeaks = len(peakShape) - ('B' in peakShape)
            prims = []
            norms = [0.9/len(peakShape)]*len(peakShape)
            const = consts[0]
            for (shape,pos) in zip(peakShape,peakPos):
                print (peakPos)
                prim = const()
                if shape == 'P':
                    # normalization, width (standard dev), position
                    prim.p[:3] = [0.03,pos]
                    #norms.append(0.9/len(peakShape))
                    twoside = len(prim.p) > 2
                    if twoside:
                        prim.p[1:] = [0.03,pos]
                        prim.free[1] = False
                        prim.lock_widths = True
                    prim.peak = True # definition
                    prims += [prim]
                    #npeaks += 1
                elif shape == 'B':
                    prims += [lp.LCGaussian(p=[0.2,pos])]
                    prims[-1].peak = False
                    #norms.append(0.3)
                else:
                    print ("Do not recognize the shape!"; continue)

            for const in consts:
                dummy  = const()
                twoside = len(dummy.p) > 3
                # replace initial gaussian shapes with proper ones
                for iprim,prim in enumerate(prims):
                    if prim.peak:
                        prims[iprim] = lp.convert_primitive(prim,ptype=const)
                        prims[iprim].peak = True
                    else:
                        # make a copy of the bridge component
                        prims[iprim] = lp.LCGaussian(p=prim.p.copy())
                        prims[iprim].peak = False
                print ('Fitting with %s'%(dummy.name))

                lct = lt.LCTemplate(prims[:],norms[:]) # shallow copy
                binned_bins = 512 # perhaps overkill
                lcf = lf.LCFitter(lct,phases,weights,binned_bins=binned_bins)
                #unbinned = len(weights)<20000
                unbinned = False
                if not unbinned:
                    print ('We are using binned likelihood with %d bins.'%(binned_bins))
                else:
                    print ('We are using unbinned likelihood.')
                t1 = time.time()
                print ('Initial values')
                print (lcf)
                lcf.fit(quick_fit_first=False,unbinned=unbinned,estimate_errors=False)
                if twoside:
                    for prim in prims: prim.free[2] = True
                lcf.fit(quick_fit_first=False,unbinned=unbinned,estimate_errors=False)
                #print ('\nFitted values')
                #print (lcf)
                t2 = time.time()
                print ('Required %ds for light curve fitting.'%(t2-t1))
                print ('Loglikelihood for %s: %.2f'%(dummy.name,lcf.ll))

                def primitive_string(prim,postfix=''):
                    l = prim.get_location(error=True)
                    loc = ['%.4f'%l[0],'%.4f'%l[1]]
                    w = prim.get_width(error=True,hwhm=True,right=False)
                    w_l = ['%.4f'%w[0],'%.4f'%w[1]]
                    w = prim.get_width(error=True,hwhm=True,right=True)
                    w_r = ['%.4f'%w[0],'%.4f'%w[1]]
                    return ( ['%s_%s'%(x,postfix) for x in ['LOC','HWHM_L','HWHM_R']],
                             [loc, w_l, w_r] )
            
                if (loglike < lcf.ll) or (loglike==-np.inf): # always write at least one model
                    # only calculate errors for models we're saving
                    if lcf.hess_errors():
                        self.PSUEoutfile.set('BOOTSTRAP','False')
                        print ('Using hessian errors.')
                    else:
                        self.PSUEoutfile.set('BOOTSTRAP','True')
                        print ('Using bootstrap errors.')
                        fit_kwargs = dict(unbinned=unbinned,unbinned_refit=False)
                        try:
                            lcf.bootstrap_errors(set_errors=True,fit_kwargs=fit_kwargs)
                        except ValueError:
                            print ('Could not estimate bootstrap errors.')
                    norms[:] = lcf.template.norms()
                    print ('Saving values for %s'%(dummy.name),loglike,lcf.ll)
                    print (lcf)
                    if not np.isnan(lcf.ll): loglike = lcf.ll
                    best_template = lcf.template
                    self.PSUEoutfile.set('N_PEAKS','%0.f'%npeaks)

                    it_peak=1
                    for (shape,peak) in zip(peakShape,prims):
                        if shape == 'P':
                            label = 'PEAK%d'%it_peak; postfix = 'P%d'%it_peak
                            it_peak += 1
                        elif shape == 'B':
                            label = 'BRIDGE'; postfix = 'B'
                        else:
                            print ("Do not recognize the shape!"; continue)
                        self.PSUEoutfile.set(label,peak.name)
                        keys,vals = primitive_string(peak,postfix)
                        for k,v in zip(keys,vals):
                            self.PSUEoutfile.set(k,v)
                    
                    try: loc_p1, err_loc_p1 = self.PSUEoutfile.get('LOC_P1',first_elem=False,type=float)
                    except (KeyError,TypeError): loc_p1, err_loc_p1 = 0, 0
                    
                    try: loc_p2, err_loc_p2 = self.PSUEoutfile.get('LOC_P2',first_elem=False,type=float)
                    except (KeyError,TypeError): loc_p2, err_loc_p2 = 0, 0

                    # delta, DELTA
                    if self.pulsar_class&1:
                        delta,err_delta = -0.1,0 # radio quiet
                    else:
                        delta,err_delta = loc_p1,err_loc_p1
                    #if delta > 0.9: delta = 1-delta
                    self.PSUEoutfile.set('delta',[delta,err_delta])

                    if npeaks > 1:
                        # use P3 for \Delta if present
                        try:
                            loc_p2, err_loc_p2 = self.PSUEoutfile.get('LOC_P3',first_elem=False,type=float)
                        except (KeyError,TypeError):
                            pass
                        DELTA, err_DELTA = np.mod(loc_p2-loc_p1,1), np.sqrt(err_loc_p1**2+err_loc_p2**2)
                        self.PSUEoutfile.set('DELTA',['%.3f'%DELTA,'%.3f'%err_DELTA])

                    # uncertainty from DM
                    ephem = ParFile(self.parfile[0])
                    self.PSUEoutfile.set('DM_delta_UNC',
                        ephem.get_dm_alignment_uncertainty(dme=dm_unc))

                    # uncertainty from TRES
                    self.PSUEoutfile.set('TRES_UNC',
                        ephem.get_tres_uncertainty())


                    # figure
                    nbins = 25
                    if H > 150: nbins = 50
                    if H > 2000: nbins = 100
                    print ('Using %d bins'%nbins)
                    fig = pl.figure(1);
                    lcf.plot(fignum=1,nbins=nbins); pl.show()
                    fig.savefig(ft1file.replace('.fits','_pulse_profile_template.png'))
                    pl.clf()

                    # save object to a file
                    outfile = open(join(self.outdir_local,'%s_profile.pickle'%(self.psrname)),'w')
                    outdata = dict(lct=lct,lcf=lcf,plc=plc)
                    cPickle.dump(outdata,outfile)
                    outfile.close()

        d = self.PSUEoutfile.get('delta',first_elem=False,type=float)
        D = self.PSUEoutfile.get('DELTA',first_elem=False,type=float)
        delta = None
        if d is not None:
            delta = (d[0],d[1],-1,0)
            if D is not None:
                delta = (d[0],d[1],D[0],D[1])
        nbands = 5
        outascii = ft1file.replace('.fits','_profile.asc')
        offpulse = [ [x+phase_shift for x in y] for y in self.offpulse ]
        plc.plot_lightcurve(nbands=nbands, background=bkg, zero_sup=True,
            substitute=False, color='black', reg=reg, template=best_template,
            outfile=pulse_profile,outascii=outascii,
            period=self.period,inset=1,suppress_template_axis=True,
            htest=H,delta=delta,offpulse=offpulse,
            pulsar_class=self.pulsar_class)

        plc.plot_lightcurve(nbands=nbands, background=bkg, zero_sup=True,
            substitute=False, color='black', reg=reg, template=best_template,
            outfile=pulse_profile.replace('.eps','.png'),
            period=self.period,inset=1,suppress_template_axis=True,
            htest=H,delta=delta,offpulse=offpulse,
            pulsar_class=self.pulsar_class)

    def PointlikeSpectralAnalysis( self, free_roi=5, phase_range=None, phname="", plot=False, **kwargs):
        '''Perform a full spectral analysis using pointlike
        =============    ===========================================================
        keyword          description
        =============    ===========================================================
        fit_emin         Minimum energy to fit (MeV)                        [100]
        fit_emax         Maximum energy to fit (MeV)                        [100000]
        free_roi         Sources within roi (deg) have free spectral params [5]
        phase_range      Phase selection (e.g [0.2,0.4] or [[0.2,0.3],[0.8,0.1]])
        plot             If True: generate SED, TSmap, residuals            [False]
        '''
        
        from uw.like.pointspec import DataSpecification, SpectralAnalysis
        from skymaps import SkyDir, SkyImage, PySkyFunction
        from uw.like.Models import PLSuperExpCutoff, PowerLaw

        if not isfile(self.ft1file):
            print ("Error. Cannot open %s! Must create FT1file. Exiting ..." %self.ft1file; exit()                        )

        ft1file = join(self.outdir,basename(self.ft1file))
        ft1file_gtis = join(self.outdir,basename(self.ft1file_gtis))
        # ft2file = join(self.outdir,basename(self.ft2file))

        if self.batch:
            if not isfile(ft1file): copy(self.ft1file,self.outdir)
            # if not isfile(ft2file): copy(self.ft2file,self.outdir)

        # apply a mktime cut
        filter = "DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52"
        GtApp('gtmktime').run(evfile=ft1file,outfile=ft1file_gtis,scfile=self.ft2file,filter=filter,roicut='no',chatter=4)
        if self.batch: copy(ft1file_gtis,self.outdir_local)
        
        phase_factor = 1
        if phase_range is not None:
            ft1file_gtis, phase_factor = SelectPhase(ft1file_gtis,phmin=phase_range[0],phmax=phase_range[1])
            phname = "_OFF"
            print ("phase factor = %.2f" %phase_factor            )

        ltcube = ft1file_gtis.replace('.fits','_ltcube.fits')
        if self.batch and isfile(join(self.outdir_local,basename(ltcube))):
            copyfile(join(self.outdir_local,basename(ltcube)),ltcube)
        
        # binfile
        binfile = ft1file_gtis.replace('.fits','_binned.fits')
        if isfile(binfile): print ("remove binfile ..."; remove(binfile))

        # Create DataSpecification Object
        from uw.like.pointspec import DataSpecification 
        ds = DataSpecification( ft1files=ft1file_gtis, ft2files=self.ft2file, ltcube=ltcube, binfile=binfile )

        # Create a SpectralAnalysis Object
        from uw.like.pointspec import SpectralAnalysis
        roi_center = SkyDir(self.ra,self.dec)
        sa = SpectralAnalysis( ds,
                               binsperdec    = 8,
                               emin          = 100,             # emin and emax only affect the *binning*
                               emax          = 1e5,
                               roi_dir       = roi_center,
                               irf           = self.irfs,
                               zenithcut     = self.zmax,
                               exp_radius    = self.radius+5,   # radius for which to calculate the livetime
                               maxROI        = self.radius,     # ROI is energy-dependent; give max ROI
                               minROI        = self.radius )    # give min ROI; equal values implies no energy dependence

        # copy ltcube to the local directory
        if self.batch: copy(ltcube,self.outdir_local)
        
        # Generate a xml source model
        xml = xml_manager(self.xmlcat,template_dir=self.template_dir)

        self.srcmdl_in = ft1file_gtis.replace('.fits','_srcmdl_input.xml')
        xml.fill_srclist( ra=self.ra, dec=self.dec, max_roi=self.radius+5, free_roi=free_roi,
                          galactic=self.galdiff, isotropic=self.isotropic )    
        xml.write_srcmdl(filename=self.srcmdl_in)
        xml.print_srclist()
        if self.batch: copy(self.srcmdl_in,self.outdir_local)        
        
        # Create a ROI object
        roi = sa.roi_from_xml( roi_dir      = roi_center,
                               xmlfile      = self.srcmdl_in,
                               fit_emin     = self.fit_emin,
                               fit_emax     = self.fit_emax,
                               phase_factor = phase_factor)

        # ======> FIT <=======
        srcname = xml.get_srcname()
        
        for i, name in enumerate(roi.get_names()):
            if name == srcname: which = i

        model = roi.get_model(srcname)

        # fix the position to the timing position
        roi.modify_loc(skydir=SkyDir(self.ra,self.dec),which=srcname)
        
        if model.name != 'PLSuperExpCutoff':
            roi.modify_model(which=srcname,model=PLSuperExpCutoff())
            model = roi.get_model(srcname)
            model['b'] = 1; model.freeze('b')

        # fit using a power law
        print ("SPECTRUM: POWER LAW")
        model['Cutoff'] = 500000.; model.freeze('Cutoff')        # fix the energy cutoff
        roi.fit(use_gradient=True)
        ll_pl = -1*roi.logLikelihood(roi.parameters())        
        print (roi)
        
        # fit using a plexpcutoff
        print ("SPECTRUM: POWER LAW WITH EXPCUTOFF")
        model['Cutoff'] = 5000. ; model.freeze('Cutoff',False) 
        roi.fit(use_gradient=True)
        ll_plec = -1*roi.logLikelihood(roi.parameters())
        print (roi)

        # print the analysis summary
        roi.print_summary()

        # save parameters (flux, index, cutoff, ts)
        self.PSUEoutfile.set("#PTLIKE_SPECANA",'................')
        self.PSUEoutfile.set("SRCNAME",srcname)
        angsep = xml.get(which=0,key='angsep')
        self.PSUEoutfile.set("ANGSEP",['%.2f'%angsep,'deg'])
        phflux, err_phflux = xml.get(which=0,key='phflux'), xml.get(which=0,key='err_phflux')
        self.PSUEoutfile.set("2FGL_F100",['%.2e'%phflux,'%.2e'%err_phflux,'ph/cm2/s'])                       
        eflux, err_eflux = xml.get(which=0,key='eflux')/6.2415e5, xml.get(which=0,key='err_eflux')/6.2415e5
        self.PSUEoutfile.set("2FGL_G100",['%.2e'%eflux,'%.2e'%err_eflux,'erg/cm2/s'])
        self.PSUEoutfile.set("EMIN",['%.0f'%self.fit_emin,'MeV'])
        self.PSUEoutfile.set("EMAX",['%.0f'%self.fit_emax,'MeV'])
        self.PSUEoutfile.set("IRFS",self.irfs)
        self.PSUEoutfile.set("MODEL",model.name)

        if phase_range is None:
            self.PSUEoutfile.set("#PHASE-AVERAGED",'................')
        else:
            self.PSUEoutfile.set("#OFFPULSE",'................')        
        
        for i, pname in enumerate(model.param_names):
            val = model.statistical()[0][i]
            eval = val * model.statistical()[1][i]**0.5
            pname += phname
            self.PSUEoutfile.set(pname,['%.2e'%val,'%.2e'%eval])
            
        flux = model.i_flux
        phflux = flux(e_weight=0,error=True,emax=3e5)
        self.PSUEoutfile.set('F100'+phname,['%.2e'%(phflux[0]*phase_factor),'%.2e'%(phflux[1]*phase_factor),'ph/cm2/s'])        
        eflux = flux(e_weight=1,error=True,emax=3e5,cgs=True)
        self.PSUEoutfile.set('G100'+phname,['%.2e'%(eflux[0]*phase_factor),'%.2e'%(eflux[1]*phase_factor),'erg/cm2/s'])
        TS = roi.TS(which=srcname)
        self.PSUEoutfile.set('TSdc'+phname,'%.1f'%TS)
        TScutoff = 2*(ll_plec-ll_pl)
        self.PSUEoutfile.set('TScutoff'+phname,'%.1f'%TScutoff)
        # self.PSUEoutfile.set('LOG(LIKE)','%.1f' %ll_plec)

        # Upper limit
        if TS < 25 and phase_factor==1:
            model['Index'] = 1.5; model.freeze('Index')      # fix the spectral index
            model['Cutoff'] = 3000.; model.freeze('Cutoff')  # fix the energy cutoff
            ul = roi.upper_limit(which=0,confidence=.95)
            self.PSUEoutfile.set('95UL',['%.2e'%ul,'ph/cm2/s'])

        # XML output
        self.srcmdl_out = join(self.outdir_local,basename(ft1file_gtis.replace('.fits','_srcmdl_output.xml')))
        roi.toXML(self.srcmdl_out)

        if plot:
            outfile = join(self.outdir_local,basename(ft1file_gtis))

            # make a TS map
            print ("making a TS map ...")
            from uw.like.roi_localize import ROILocalizer
            from uw.utilities.image import TSplot
            fig = pl.figure(1)
            ax = fig.add_subplot(111)
            rl = ROILocalizer(roi,which=which)
            tsp = TSplot(PySkyFunction(rl),roi_center,0.3,pixelsize=0.03)
            tsp.cross(roi_center,size=0.01,label=self.psrname,color='blue')
            tsp.show()
            pl.savefig(outfile.replace('.fits','_tsmap.png')); pl.clf()
        
            # make a SED
            from uw.like.sed_plotter import plot_sed
            # from pointlike_utils import plot_all_seds
            print ("making a sed ...")
            fig = pl.figure(2,(9,7))
            plot_sed(roi,which=srcname,fignum=2,use_ergs=True,axis=(90,1.1e5,1e-13,eflux[0]))
            pl.savefig(outfile.replace('.fits','_sed.png'))
            # plot_all_seds(roi,fignum=3,filename=ft1file.split('.fits')[0] + '_allseds.png')
            
            # make counts map
            print ("making count map ...")
            roi.plot_counts_map(filename=outfile.replace('.fits','_count_map.png'))
            
            # make counts spectra
            print ("making count spectra ...")
            roi.plot_counts_spectra(filename=outfile.replace('.fits','_count_spectra.png'))
            
            # make a residual TS map
            print ("make a residual ts map ...")
            roi.plot_tsmap(filename=join(outfile.replace('.fits','_residual_tsmap.png')),figsize=(5.5,4.5))
            
            # smoothed counts
            print ("make a smoothed count map ...")
            roi.plot_source(which=srcname,filename=outfile.replace('.fits','_smoothed_counts.png'))

    def LivetimeCube(self,clobber=False):
        """ Generate a livetime cube with gtltcube."""
        if (not clobber) and isfile(self.gtltcube):
            print ('Livetime cube already exists; returning.')
        if not isfile(self.ft1file_gtis):
            raise IOError('Could not find gtmktimed FT1 file to make livetime cube.')

        ft1file_gtis = join(self.outdir,basename(self.ft1file_gtis)) # potentially scratch
        if self.batch and (not isfile(ft1file_gtis)):
            copy(self.ft1file_gtis,self.outdir)
        gtltcube = join(self.outdir,basename(self.gtltcube))
        GtApp('gtltcube','Likelihood').run(evfile=ft1file_gtis,scfile=self.ft2file,
            outfile=ltcube,dcostheta=0.025,binsz=1)
        if self.batch:
            copy(gtltcube,self.outdir_local) # copy the ft1 file to the local directory

    def WriteWeights(self, ozlem_xml=None, clobber=False, outdir=None):
        '''Add the weights into a FT1 file from a source model.
        =============    ===============================================
        keyword          description
        =============    ===============================================
        ozlem_xml        XML file from Ozlem's spectral fits [None]
        clobber          Overwrite files                     False
        '''
        if ozlem_xml is None:
            if not isfile(self.srcmdl_out):
                print ("Error. Cannot find %s! Exiting ..." %(self.srcmdl_out); exit())
        else:
            if not isfile(ozlem_xml):
                print ("Error. Cannot find %s! Exiting ..." %(ozlem_xml); exit())
            
        # Generate a new xml source model
        xml = xml_manager(ozlem_xml or self.srcmdl_out,template_dir=self.template_dir)        
        srcmdl = ozlem_xml or join(self.outdir_local,basename(self.ft1file_gtis.replace('.fits','_srcmdl_weights.xml')))
        xml.fill_srclist(ra=self.ra,dec=self.dec,max_roi=self.w_radius)
        if ozlem_xml is None:
            xml.write_srcmdl(filename=srcmdl)
            xml.print_srclist()
        srcname = xml.get_srcname(which=0)        

        print ('Using this for diffuse response')
        print (self.ft1file_diffrsp)
        # this is the local name of file (e.g. if running batch mode)
        ft1file_diffrsp = join(self.outdir,basename(self.ft1file_diffrsp))

        # check that gtmktimed file exists
        if not isfile(self.ft1file_gtis):
            ft1file = join(self.outdir,basename(self.ft1file))
            ft1file_gtis = join(self.outdir,basename(self.ft1file_gtis))
            if self.batch:
                if not isfile(ft1file): copy(self.ft1file,self.outdir)
            filter = "DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52"
            GtApp('gtmktime').run(evfile=ft1file,outfile=ft1file_gtis,scfile=self.ft2file,filter=filter,roicut='no',chatter=4)
            if self.batch:
                copy(ft1file_gtis,self.outdir_local)

        if (not isfile(self.ft1file_diffrsp)) or clobber:
            # select a sub-ft1 file
            GtApp('gtselect').run(infile=self.ft1file_gtis,outfile=ft1file_diffrsp,
                rad=self.w_radius,ra=self.ra,dec=self.dec,emin=self.w_emin,
                emax=self.w_emax,tmin=self.tmin,tmax=self.tmax,zmax=self.zmax)
        elif self.batch:
            copy(self.ft1file_diffrsp,ft1file_diffrsp)
            
        # check the diffuse response columns
        GtApp('gtdiffrsp').run(convert='yes',evfile=ft1file_diffrsp,evtable="EVENTS",
            scfile=self.ft2file,srcmdl=srcmdl,irfs=self.irfs,edisp='no')
        if self.batch:
            copy(ft1file_diffrsp,self.ft1file_diffrsp)
        
        # create a source list
        srclist = join(self.outdir_local,'srclist.txt')
        srclist_file = open(srclist,'w')
        srclist_file.write(srcname)
        srclist_file.close()
        
        print ('Using this for weights')
        print (self.ft1file_weights)
        # run gsrcprob
        print ("add weights into the FT1 file ...")
        GtApp('gtsrcprob').run(evfile=ft1file_diffrsp,scfile=self.ft2file,
            outfile=self.ft1file_weights,srcmdl=srcmdl,irfs=self.irfs,srclist=srclist)

        # rename the weight column to WEIGHT
        hdulist = pyfits.open(self.ft1file_weights,mode='update')
        col = hdulist[1].columns
        for c in col:
            if c.name == srcname: c.name = 'WEIGHT'
        hdulist[1].update()
        hdulist.writeto(self.ft1file_weights,clobber=True)
        hdulist.close()

        # clean
        #if batch: # think this will happen in CleanFiles
        #    remove(ft1file_diffrsp)

