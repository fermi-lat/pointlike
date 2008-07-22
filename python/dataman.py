from uw.sane import setup
setup()
from GtApp import GtApp

from sys import path
path.insert(0,'d:/users/kerrm/python/spectrum_dev4')
path.insert(0,'d:/users/kerrm/python/spectrum_dev4/uw')
import pointlike
path.pop(0)
reload(pointlike)

import glob, os
from pointlike import Data
import numpy as N
from numpy import array

"""
Alignment constants

"""

alignment = [
    # start     stop        x     y    z
(236511638, 	237149305, 	-161, 	-170, 	-475 ),
(237154733, 	237338067, 	-303, 	-218, 	-740 ),
(237343796, 	237777636, 	-170, 	-173, 	-491 ),
    ]

class DataManager(object):
    """Manage a list of runs and generate analysis products."""

    def __init__(self,
                 datapath = r'f:\glast\data\flight\leo',
                 analysispath=r'd:\common\first_light',
                 table=alignment,
                 ft2=None,
                 quality=None):

        #TODO -- incorporate class level cuts

        self.alignment = Alignment(table)
        self.ft1files = glob.glob('%s/r*ft1.fit'% datapath)

        print 'found %d files ' % len(self.ft1files)
        if len(self.ft1files)==0:
            print 'files not found in path %s '%datapath
            return

        if ft2:
            print 'using combined ft2 file %s ' % ft2
            Data.setHistoryFile(ft2)

        quality = quality if quality else datapath+'\qual_list.txt'
        try:
            good_quality = RunQuality(quality)
            self.ft1files = [f for f in self.ft1files if good_quality(self.run_number(f))]

        except: print 'Not applying run cuts!'
        
        self.datapath,self.analysispath=datapath,analysispath
        print 'done with initialization'


    def run_number(self,filename):
        runid = filename.replace('\\','/').split('/')[-1]
        return float(runid[1:11]) #10 digits, including inital '0'

    def in_range(self, runid, r):
        t = self.run_number(runid)
        return t>= r[0] and t<r[1]

    def set_bands(self, b = [100*10**(0.2*x) for x in xrange(18)]):
        Data.set_bands(b)
    
    def write_binned(self,sum_alignment = False):
        
        Data.set_rot(self.alignment[0][2:5])
        fl = [f for f in self.ft1files if self.in_range(f, self.alignment[0])]
        d0 = Data(fl)

        print d0.map().info()

        if len(self.alignment)>1:
            for s in self.alignment[1:]:
                Data.set_rot(s[2:5])
                fl = [f for f in self.ft1files if self.in_range(f, s)]
                if len(fl)==0: continue
                d = Data(fl)

                if sum_alignment: d0.map().add(d.map())
                else:
                    outfile = '%s/r%d_%d_bands.fits' % (self.analysispath, s[0], s[1])
                    d.map().write(outfile)
                    print 'wrote bands file to %s' %outfile

        final_index = -1 if sum_alignment else 0
        s = self.alignment
        outfile = '%s/r%d_%d_bands.fits' % (self.analysispath, s[0][0], s[final_index][1])
        d0.map().write(outfile)

    
    def livetime_cubes(self):
        
        ltfiles = glob('%s/r*ltcube.fits'%self.datapath)
        done_runs = [self.run_number(f) for f in ltfiles]
        to_do_ft1_files = [f for f in self.ft1files and self.run_number(f) not in done_runs]
         
        gt=GtApp('gtltcube', 'General')
        
        for r in to_do_ft1_files:

            gt.pars['evfile']=r
            gt.pars['scfile']=r.replace('ft1','ft2')
            gt.pars['outfile']=r.replace('ft1','ltcube')
            gt.pars['dcostheta']=0.25
            gt.pars['clobber']='yes'
            gt.run()

        f = open(self.datapath+'ltcube_list.txt','w')
        g = open(self.datapath+'ft1_list.txt','w')
        for r in to_do_ft1_files:
            f.write(r.replace('ft1','ltcube')+'\n')
            g.write(r+'\n')
        f.close()
        g.close()

    def sum_livetime_cubes(self):
    
        gt=GtApp('gtltsum','General')
        gt.pars['infile1']='@%s/ltcube_list.txt'%self.datapath
        gt.pars['outfile']='%s/summed_ltcube.fits'%self.datapath
        gt.pars['clobber']='yes'
        gt.run()

    def write_exposure_maps(self):

        #irfs_dict={'1':'TRANSIENT', '2':'SOURCE','3':'DIFFUSE}
        #irfs = 'P6_V1_'+irfs_dict(self.class_level)
        ifrs = 'P6_V1_DIFFUSE'

        enumbins = 100
        emin = 29.998
        emax = 300000.01

        gt=GtApp('gtexpcube', 'map_tools')
        gt.pars['infile']='%s/summed_ltcube.fits'%self.datapath
        gt.pars['evfile']='%s/ft1_list.txt'%self.datapath
        gt.pars['cmfile']='NONE'

        gt.pars['enumbins']=enumbins
        gt.pars['emin']=emin
        gt.pars['emax']=emax
        gt.pars['bincalc']='EDGE'
        gt.pars['nxpix']=1
        gt.pars['nypix']=1
        gt.pars['pixscale']=1
        gt.pars['proj']='CAR'
        gt.pars['coordsys']='CEL'
        gt.pars['clobber']='yes'
        gt.pars['xref']=0
        gt.pars['yref']=0

        enumbins,emin,emax=N.array([enumbins,emin,emax]).astype(int)

        gt.pars['outfile']='%s/%s_%ib_%i_%i_%s.fits'%(self.datapath,irfs,enumbins,emin,emax,'front')
        gt.pars['irfs']=irfs+'::FRONT'
        gt.run(CatchError=None)

        gt.pars['outfile']='%s/%s_%ib_%i_%i_%s.fits'%(self.datapath,irfs,enumbins,emin,emax,'back')
        gt.pars['irfs']=irfs+'::BACK'
        gt.run(CatchError=None)
       

class Alignment(list):
    def __init__(self, table=alignment):
        self += table

    def __call__(self, run):
        for t in self:
            if run>= t[0] and run<t[1]:
                return t[2:]
        return None

class RunQuality(object):

    def __init__(self,filename):
        try:
            f = open(filename)
            tokens = f.readlines()
            f.close()
        except: print 'Could not open %s'%filename

        import numpy as N

        for i in xrange(len(tokens)):
            tokens[i] = tokens[i].strip().split()

        runs = N.array([token[0] for token in tokens])

        mask_qual = N.asarray([token[1]=='Good' for token in tokens])
        mask_mode = N.asarray([N.any(['nomSciOps' in t for t in token]) for token in tokens])

        self.good_runs = N.array(runs[N.logical_and(mask_qual,mask_mode)]).astype(float)

    def __call__(self,run):
        return run in self.good_runs


    
if __name__=='__main__':
    d = DataManager(analysispath=r'd:\common\first_light\aligned_qual')

    #Test
    d.write_binned(sum_alignment = True)
    d.write_binned(sum_alignment = False)
