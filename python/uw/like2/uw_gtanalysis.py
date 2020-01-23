"""
Provide interface to the fermipy GTAnalysis system
"""
import os,glob, healpy, yaml, logging
import numpy as np 
import pandas as pd 
from astropy.io import fits
from fermipy import gtanalysis 
from . import (configuration,)
from uw.data import binned_data
from skymaps import Band

config_template="""\
    data:
        evfile: {ft1}
        scfile: null
        ltcube: {ltcube}
        
    binning:
        projtype: 'HPX'
        enumbins: 1
        hpx_order: 6 #64 override in components
        coordsys: 'GAL'
        roiwidth : 10.0 #radius of 5 deg.

    selection:
        ra: {ra}
        dec: {dec}
        zmax: 100
        thmax: 66.4 # arccos(0.4) THB added

    gtlike:
        edisp: False
        use_external_srcmap: True
        irfs: {irf}
        bexpmap: {bexpmap} 
        ccube: {ccube} # added to standard
        
    model:
        galdiff: {galdiff}
        isodiff: # set in components
        # what about Sun/Moon?
        #catalogs: 
        #    - '3FGL'
        sources: []

    components: {components}

    fileio:
        outdir: 'fermipy'
    logging:
        prefix: ''
        chatter: 3
        verbosity: 3
    """   

# logging.config.dictConfig( dict(
#     version = 1,
#     formatters = {
#         'f': {'format':
#               '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}
#         },
#     handlers = {
#         'h': {'class': 'logging.StreamHandler',
#               'formatter': 'f',
#               'level': logging.DEBUG}
#         },
#     root = {
#         'handlers': ['h'],
#         'level': logging.DEBUG,
#         },
#     )
# )

class MyLogger(object):
    def __init__(self, uselogger=True):
        if uselogger:
            self.logger = logging.getLogger()
        else:
            self.logger = None
    def info(self, msg):
        if self.logger: self.logger.info(msg)
        else: print (msg)
    def warning(self, msg):
        if self.logger: self.logger.warning(msg)
        else: print ('WARNING', msg)
    def error(self, msg):
        if self.logger: self.logger.error(msg)
        else: print ('ERROR', msg)
    def log(self, level, msg, appname='', arg=None):
        if self.logger: self.logger.log(level, msg, appname, arg )
        else: 
            if arg is not None: 
                print ('INFO', msg % (appname, arg))
            else: print ('INFO', msg)

class UWtoGT(object):
    """ manage adaption of a UW model to gtanalysis
    """
    def __init__(self, path='.' ,uselogger=True ):
        self.uwcfg= cfg= configuration.Configuration(path, quiet=True, postpone=True)
        self.datapath= os.path.split(cfg.dataset.binfile)[0]
        self.logger = logger = MyLogger(uselogger)
        logger.info('start UW interface to fermipy in folder {}'.format(os.getcwd()))
        os.environ['CUSTOM_IRF_NAMES']='' # it wants this set
        if False:
            os.environ['CALDB']=cfg.caldb
            self.irf = cfg.irf
        else:
            self.irf='P8R2_SOURCE_V6' 
            logger.warning('CALDB set to {} for irf {}'.format(os.environ['CALDB'], self.irf))

        diffuse_dir =os.path.expandvars('$FERMI/diffuse/') 
        self.galdiff = diffuse_dir+'gll_iem_v06.fits' # cfg['diffuse']['ring']['filename']
        self.isodiff = [diffuse_dir+cfg['diffuse']['isotrop']['filename'].replace('**',x) for x in('FRONT','BACK')] 
        self.ltcube = cfg.dataset.ltcube

        # expect to find "ccube" files here, under folders with roi index values
        self.ccube_root=os.path.join(self.datapath,'8years_ccube')
        # the "bexpmap" files go herea
        self.bexpmap_dir=os.path.join(self.datapath, '8years_bexpmap')

        # extract component info from the FITS file
        t = fits.open(cfg.dataset.binfile)
        assert 'BANDS' in t, 'Expected new-format binned data file'
        try:
            cmp=[]
            for i, (nside, emin, emax, et) in enumerate(t['BANDS'].data):
                cmp.append(
                    dict(
                        name='{:02d}'.format(i),
                        selection=dict(evtype=1<<et,
                            emin=1e3*emin, emax=1e3*emax,), #convert to GeV??
                        binning=dict(hpx_order=int(np.log2(nside)), enumbins=1), ##allow to change
                        model=dict(
                            isodiff=self.isodiff[et],
                            galdiff = self.galdiff,
                            ),
                    )
                )  
            self.components=cmp
        except Exception as msg:
            logger.error('Failed to interpret file {}: {}'.format(cfg.dataset.binfile, msg))
            raise
        except:
            logger.error('Failed to interpret file {}'.format(cfg.dataset.binfile))

    
    def __call__(self, roi_index, ncomp=None, overwrite=False, **kwargs):
        """Setup gtanalysis with a pointlike ROI
        """
        self.roi_index = roi_index
        roidir = Band(12).dir(roi_index)
        logger=self.logger

        gt_config = yaml.load(config_template.format(ltcube=self.ltcube,
                ft1='', 
                irf=self.irf, 
                galdiff=self.galdiff,
                components= self.components if ncomp is None else self.components[:ncomp], # to 100 GeV 
                ra=roidir.ra(), dec=roidir.dec(),
                bexpmap=self.bexpmap_dir+'/bexpmap%s.fits',
                ccube=self.ccube_root+'/{:04d}/ccube%s.fits'.format(roi_index),
                #srcmdl='srcmdl.xml', # same for all components
                ))
        self.gtcfg = gt_config
        # Check and/or setup the ccube files
        assert os.path.exists(self.ccube_root), 'Expected to find folder {}'.format(self.cube_root)
        roi_root = os.path.join(self.ccube_root, str(roi_index))
        if not os.path.exists(roi_root):
            os.mkdir(roi_root)

        clist = range(len(gt_config['components']))
        ccfiles = [os.path.join(self.ccube_root,\
             '{:04d}'.format(roi_index), gt_config['gtlike']['ccube']% '_{:02d}'.format(j) )for j in clist]
        ok = np.all(map(os.path.exists, ccfiles))
        if not ok:
            logger.info( 'Need to generate ccfiles')
            bf = binned_data.BinFile(self.uwcfg.dataset.binfile)
            bf.generate_ccube_files(self.ccube_root, roi_index,   channels=clist, overwrite=overwrite)

        # Check for bexpmap
        assert os.path.exists(self.bexpmap_dir), 'expected to find folder {}'.format(self.bexpmap_dir)
        self.create_expcubes(gt_config['components'], overwrite=overwrite)
        
        logger.info('Create the gtanalysis')
        self.gta = gtanalysis.GTAnalysis(gt_config, **kwargs)

        # Check for srcmaps
        logger.info('Checking srcmap files')
        srcmap_ok = [os.path.exists(c.files['srcmap']) for c in self.gta]
        
        if np.any(np.logical_not(srcmap_ok)) or overwrite:
            logger.info('Will create srcmap files')
            try:
                self.create_srcmaps(overwrite)
            except Exception as msg:
                logger.error('Failed : {}. Returning anyway'.format(msg))
                return self.gta
        
        logger.info('GRanalysis ready to complete setup')
        return self.gta

    def file_check(self):
        checklist = ' ccube bexpmap srcmdl srcmap'.split()
        fstat = pd.DataFrame(
            np.array([[os.path.isfile(cmp.files[check]) for check in checklist] for  cmp in self.gta]),
            columns=checklist)
        return fstat

    def create_expcubes(self, components, overwrite=False):
        logger= self.logger
        cfg = self.gtcfg

        def create_expmap(cmp):
            cmap='none'
            appname = 'gtexpcube2'
            outfile = cfg['gtlike']['bexpmap'] % ('_'+cmp['name'])
            if os.path.isfile(outfile) and not overwrite:
                logger.info( 'skipping {}'.format(outfile))
                return
            kw = dict(infile=cfg['data']['ltcube'], cmap=cmap,
                    ebinalg='LOG',
                    emin=cmp['selection']['emin'],
                    emax=cmp['selection']['emax'],
                    enumbins=cmp['binning']['enumbins'],
                    outfile=outfile, 
                    proj='CAR',nxpix=360, nypix=180, binsz=1,
                    xref=0.0, yref=0.0,
                    evtype=cmp['selection']['evtype'],
                    thmax=cfg['selection']['thmax'], #THB
                    irfs=cfg['gtlike']['irfs'],
                    coordsys=cfg['binning']['coordsys'],
                    chatter=cfg['logging']['chatter']
                    )
            gtanalysis.run_gtapp(appname, logger, kw, loglevel=logging.INFO)
            if not os.path.isfile(kw['outfile']):
                logger.error( 'Failed to create output file')
                logger.error( '{} keywords:'.format(appname))
                for k,v in kw.items():
                    logger.error('  {:10s}: {}'.format(k,v))
            else:
                logger.info('created {}'.format(outfile))
        map(create_expmap, components)

    def create_srcmaps(self, overwrite=False):
        for cmp in self.gta:
            if not os.path.isfile(cmp.files['srcmdl'] or overwrite):
                cmp.roi.write_xml(cmp.files['srcmdl'], cmp.config['model'])
            cmp._create_srcmaps(overwrite=overwrite)

    def __getitem__(self, i):
        return self.components[i]
    def __len__(self):
        return len(self.components)

