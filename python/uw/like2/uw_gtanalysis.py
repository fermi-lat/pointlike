"""
Provide interface to the fermipy GTAnalysis system
"""
import os,glob, healpy
import numpy as np 
import pandas as pd 

from fermipy import gtanalysis
from . import (configuration,)

class GTconfig(object):
    def __init__(self, path='.'):
        self.uwcfg= cfg= configuration.Configuration(path, quiet=True, postpone=True)
        os.environ['CUSTOM_IRF_NAMES']='' # it wants this set
        if False:
            os.environ['CALDB']=cfg.caldb
            self.irf = cfg.irf
        else:
            self.irf='P8R2_SOURCE_V6' 
        print 'CALDB set to {} for irf {}'.format(os.environ['CALDB'], self.irf)

        diffuse_dir =os.path.expandvars('$FERMI/diffuse/') 
        self.galdiff = diffuse_dir+'gll_iem_v06.fits' # cfg['diffuse']['ring']['filename']
        self.isodiff = [diffuse_dir+cfg['diffuse']['isotrop']['filename'].replace('**',x) for x in('FRONT','BACK')] 

