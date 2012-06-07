import glob
import numpy as np
from CLHEP import Hep3Vector
rd = 180./np.pi

## Data specification class for stacked analysis
#  manages files locations 

basedir = '/phys/groups/tev/scratch1/users/Fermi/mar0'
datadir = basedir+'/data/7.3src/'
srcdir = basedir+'/pointlikedev/uw/stacklike/lists/'

class SDataSpec(object):

    def __init__(self,**kwargs):
        self.name=''
        self.filedir=''
        self.srcdir=''
        self.trange=[]
        self.ft1glob=''
        self.ft2glob=''
        self.phase=[]
        self.irf = ''
        self.rot=[0,0,0]
        self.binfile=None
        self.__dict__.update(kwargs)
        self.loadsrcs()
        self.useft2s = not self.ft2glob==''
        self.ft1s = np.sort(glob.glob(self.filedir+self.ft1glob))
        self.ft2s = np.sort(glob.glob(self.filedir+self.ft2glob))
    
    def loadsrcs(self):
        self.srcs=[]
        ff = open(self.srcdir+self.name+'.txt')
        header = ff.readline()
        for lines in ff:
            line = lines.split()
            self.srcs.append(Hep3Vector([float(line[1])/rd,float(line[2])/rd]))

alignspec = SDataSpec(name='agn-psf-study-bright',
                        filedir=datadir, 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='[0-9][0-9]*-ft1.fits',
                        ft2glob='[0-9][0-9]*-ft2.fits',
                        irf = 'P7SOURCE_V4MC')

alignspec_repro = SDataSpec(name='agn-psf-study-bright',
                        filedir=basedir+'/data/7.3srcrp/', 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='[0-9][0-9]*-ft1.fits',
                        ft2glob='[0-9][0-9]*-ft2.fits',
                        irf = 'P7SOURCE_V4MC')

velaon = SDataSpec(name='velaon',
                        filedir=basedir+'/data/pulsar7/source/', 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='vela-ft1.fits',
                        ft2glob='psr-ft2.fits',
                        irf = 'P7SOURCE_V11',
                        phase=[0.1,0.15,0.5,0.6]
                        )

velaoff = SDataSpec(name='velaoff',
                        filedir=basedir+'/data/pulsar7/source/', 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='vela-ft1.fits',
                        ft2glob='psr-ft2.fits',
                        irf = 'P7SOURCE_V11',
                        phase=[0.7,1.0]
                        )

gemon = SDataSpec(name='gemon',
                        filedir=basedir+'/data/pulsar7/source/', 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='gem-ft1.fits',
                        ft2glob='psr-ft2.fits',
                        irf = 'P7SOURCE_V4MC',
                        phase=[0.1,0.17,0.6,0.68]
                        )

gemoff = SDataSpec(name='gemoff',
                        filedir=basedir+'/data/pulsar7/source/', 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='gem-ft1.fits',
                        ft2glob='psr-ft2.fits',
                        irf = 'P7SOURCE_V11',
                        phase=[0.25,0.55]
                        )

agnlobzb = SDataSpec(name='agn_redshift2_lo_bzb',
                        filedir=datadir, 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='[0-9][0-9]*-ft1.fits',
                        ft2glob='[0-9][0-9]*-ft2.fits',
                        irf = 'P7SOURCE_V11')

agnhibzb = SDataSpec(name='agn_redshift2_hi_bzb',
                        filedir=datadir, 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='[0-9][0-9]*-ft1.fits',
                        ft2glob='[0-9][0-9]*-ft2.fits',
                        irf = 'P7SOURCE_V11')

1es0229 = SDataSpec(name='1es0229p200',
                        filedir=datadir, 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='[0-9][0-9]*-ft1.fits',
                        ft2glob='[0-9][0-9]*-ft2.fits',
                        irf = 'P7SOURCE_V11')

1es0347 = SDataSpec(name='1es0347-121',
                        filedir=datadir, 
                        srcdir = srcdir, 
                        trange = [239557417,356000000],
                        ft1glob='[0-9][0-9]*-ft1.fits',
                        ft2glob='[0-9][0-9]*-ft2.fits',
                        irf = 'P7SOURCE_V11')