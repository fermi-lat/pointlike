"""
Seed processing code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/seeds.py,v 1.1 2015/08/16 01:13:19 burnett Exp $

"""
import os, sys, time, pickle, glob, pyfits, types
import numpy as np
import pandas as pd
from skymaps import SkyDir, Band
from uw.utilities import keyword_options
from uw.like2 import (tools, sedfuns, maps, sources, localization, roimodel,)


def read_seedfile(seedkey, 
    pgw_filename= '/nfs/farm/g/glast/g/catalog/transients/P302/PGW/1m_1mp15_PGW_ALL_NoSun.fits'):

    model_name = os.getcwd().split('/')[-1]

    if model_name.startswith('month') and seedkey=='pgw':
        #monthly mode, need to find and load PGW analysis with rouighly equivalent months
        month=int(model_name[5:]); 
        pgw_filename='/nfs/farm/g/glast/g/catalog/transients/TBIN_%d_all_pgw.txt'% (month-1)
        assert os.path.exists(pgw_filename), 'PGWAVE file %s not found'% pgw_filename  
        try:
            seeds = pd.read_table(pgw_filename, sep=' ', skipinitialspace=True, index_col=1,
                header=None,
                names='tbin ra dec k_signif pgw_roi fgl_seed fgl_ra fgl_dec fgl_assoc'.split())
        except Exception,msg:
            raise Exception('Failed to read file %s: %s' % (pgw_filename, msg))
        names=[]
        for i,s in seeds.iterrows():
            j = int(s.name[4:6]) if s.name[6]=='_' else int(s.name[4:5])
            names.append('PGW_%02d_%03d_%02d' % (month, int(s.pgw_roi), j))
        seeds['name'] = names    
    elif model_name.startswith('month') and seedkey=='PGW':
        # monthly mode, new format PGwave, in a single FITS file
        month=int(model_name[5:]); 

        
        assert os.path.exists(pgw_filename), 'PGWAVE file {} not found'.format( pgw_filename)  
        t = pyfits.open(pgw_filename)
        df=pd.DataFrame(t[1].data)
        selector = lambda month : (df.run=='1m   ') & (df.TBIN=='TBIN_{:<2d}'.format(month-1))
        cut = selector(month)
        assert sum(cut)>0, 'No seeds found for month {}'.format(month)
        print 'Found {} PGWave seeds'.format(sum(cut))
        ra = np.array(df.Ra[cut],float)
        dec = np.array(df.Dec[cut],float)
        prefix = 'PG{:02d} '.format(int(month))
        # note making it a string type
        name = np.array([prefix + n.split('_')[-1].strip() for n in 'TBIN_{}_'.format(month-1)+df.PGW_name[cut]])
        seeds = pd.DataFrame([name, ra,dec], index='name ra dec'.split()).T

    else:
        # reading a TS seeds file
        t = glob.glob('seeds_%s*' % seedkey)
        assert len(t)==1, 'Seed file search, using %s, failed to find one file' % seedkey
        seedfile=t[0]
        try:
            csv_format=seedfile.split('.')[-1]=='csv'
            if csv_format:
                seeds = pd.read_csv(seedfile)
            else:
                seeds = pd.read_table(seedfile)
        except Exception, msg:
            raise Exception('Failed to read file %s, perhaps empty: %s' %(seedfile, msg))
    seeds['skydir'] = map(SkyDir, seeds.ra, seeds.dec)
    seeds['hpindex'] = map( Band(12).index, seeds.skydir)
    # check for duplicated names
    dups = seeds.name.duplicated()
    if sum(dups)>0:
        print '\tRemoving {} duplicate entries'.format(sum(dups))
        return seeds[np.logical_not(dups)]
    return seeds

def add_seeds(roi, seedkey, model='PowerLaw(1e-14, 2.2)', 
            prefix=None,
            associator=None, tsmap_dir=None,
            tsmin=5, lqmax=8,
            **kwargs):
    """ add "seeds" from a text file the the current ROI
    
        roi : the ROI object
        seedkey : string
            Expect one of 'pgw' or 'ts' for now. Used by read_seedfile to find the list
        model : string
            model to use for the generated source
        prefix : None or string
            Name to 
        associator :
        tsmap_dir
        mints : float
            minimum TS to accept for addition to the model
        lqmax : float
            maximum localization quality for tentative source
    """
    seeds = read_seedfile(seedkey)
    inside = seeds.hpindex==Band(12).index(roi.roi_dir)
    seednames=seeds.name[inside]
    if sum(inside)==0:
        print 'no seeds in ROI'
        return
    if prefix is None:
        prefix = seeds.name[0][:4]
    srclist = []
    for i,s in seeds[inside].iterrows():
        try:
            srclist.append(roi.add_source(sources.PointSource(name=s['name'], skydir=s['skydir'], model=model)))
            print 'added %s at %s' % (s['name'], s['skydir'])
        except roimodel.ROImodelException:
            srclist.append(roi.get_source(s['name']))
            print 'updating existing %s at %s ' %(s['name'], s['skydir'])
    
    # Fit only fluxes for each seed first
    parnames = roi.sources.parameter_names
    seednorms = np.arange(len(parnames))[np.array([s.startswith(prefix) and s.endswith('_Norm') for s in parnames])]
    assert len(seednorms)>0, 'Did not find any seeds.'
    try:
        roi.fit(seednorms, tolerance=0.2, ignore_exception=False)
    except Exception, msg:
        print 'Failed to fit seed norms: \n\t{}\nTrying full fit'.format(msg)
        roi.fit(ignore_exception=True)
    #remove those with low TS
    print 'TS values'
    goodseeds = []
    for sname in seednames:
        ts = roi.TS(sname)
        print '%8s %5.1f ' %(sname, ts),
        if ts<tsmin:
            print '<-- remove'
            roi.del_source(sname)
        else:
            goodseeds.append(sname)
            print 'OK'
        
    # fit all parameter for each one, remove if ts< tsmin 
    seednames = goodseeds
    goodseeds = []
    for sname in seednames:
        roi.fit(sname, tolerance=0, ignore_exception=True)
        ts = roi.TS()
        print '  TS = %.1f' % ts
        if ts<tsmin:
            print ' TS<%.1f, removing from ROI' % tsmin
            roi.del_source(sname)
        else:
            # one iteration of pivot change
            s = roi.get_source(sname)
            roi.repivot([s])
            # and a localization: remove if fails or poor
            roi.localize(sname, update=True)
            ellipse = s.__dict__.get('ellipse', None) 
            if ellipse is None or ellipse[5]>lqmax:
                print '--> removing %s since did not localize' %s
                roi.del_source(sname)
                continue
            roi.get_sed(sname)
            goodseeds.append(sname)
            
    if len(goodseeds)>0:
        # finally, fit full ROI 
        roi.fit(ignore_exception=True, tolerance=0.2)
        return True
    else:
        return False
    
            
#def add_sources(self, csv_file='plots/seedcheck/good_seeds.csv'):
#    """Add new sources from a csv file
#    Any that are already in the model are ignored.
#    Default assumes analysis by seedcheck
#    """
#    assert os.path.exists(csv_file), 'csv file %s not found' % csv_file
#    good_seeds = pd.read_csv(csv_file, index_col=0)
#    print 'Check %d sources from file %s: ' % (len(good_seeds), csv_file),
#    myindex = Band(12).index(self.roi_dir)
#    inside = good_seeds['index']==myindex
#    ni = sum(inside)
#    if ni>0:
#        print '%d inside ROI' % ni
#    else:
#        print 'No sources in ROI %04d' % myindex
#        return 
#    snames = [ s.name for s  in self.free_sources]
#    for name, s in good_seeds[inside].iterrows():
#        if name in snames:
#            print 'Source %s already in model: skipping' % name
#            continue
#        e0 = s.get('e0', 1000)
#        #try: #in case already exists (debugging perhaps)
#        #    self.del_source(name)
#        #    print 'replacing source %s' % name 
#        #except: pass
#        source=self.add_source(name=name, skydir=SkyDir(s['ra'],s['dec']), 
#                    model=sources.LogParabola(s['eflux']/e0**2/1e6, s['pindex'], s.get('par2',0), e0, 
#                                              free=[True,True,False,False]))
#    # perform a preliminary fit including the new sources
#    self.fit(ignore_exception=True, tolerance=0.2)
#    # make sure they have SED and localization info (relocalizes any that were duplicated)
#    for name in good_seeds[inside].index:
#        self.get_sed(name)
#        self.localize(name, update=True)
#    return good_seeds[inside]   
#
#       
#def check_seeds(self, seedkey, model='PowerLaw(1e-14, 2.2)', 
#            prefix=None,
#            associator=None, tsmap_dir=None,
#            **kwargs):
#    """ check "seeds" from a text or csv file 
#        seedkey : string
#            If 'pgw', special check for special monthly
#        model : string
#            Model to use
#    """
#    repivot_flag = kwargs.pop('repivot', True)
#    seedcheck_dir = kwargs.get('seedcheck_dir', 'seedcheck')
#    if not os.path.exists(seedcheck_dir): os.mkdir(seedcheck_dir)
#    associator= kwargs.pop('associate', None)
#
#    model_name = os.getcwd().split('/')[-1]
#    if model_name.startswith('month') and seedkey=='pgw':
#        #monthly mode, need to find and load PGW analysis with rouighly equivalent months
#        month=int(model_name[5:]); 
#        pgw_filename='/nfs/farm/g/glast/g/catalog/transients/TBIN_%d_all_pgw.txt'% (month-1)
#        assert os.path.exists(pgw_filename), 'PGWAVE file %s not found'% pgw_filensme  
#        try:
#            seeds = pd.read_table(pgw_filename, sep=' ', skipinitialspace=True, index_col=1,
#                header=None,
#                names='tbin ra dec k_signif pgw_roi fgl_seed fgl_ra fgl_dec fgl_assoc'.split())
#        except Exception,msg:
#            raise Exception('Failed to read file %s: %s' % (pgw_filename, msg))
#        names=[]
#        for i,s in seeds.iterrows():
#            j = int(s.name[4:6]) if s.name[6]=='_' else int(s.name[4:5])
#            names.append('PGW_%02d_%03d_%02d' % (month, int(s.pgw_roi), j))
#        seeds['name'] = names    
#    else:
#        t = glob.glob('seeds_%s*' % seedkey)
#        assert len(t)==1, 'Seed file search, using %s, failed to find one file' % seedkey
#        seedfile=t[0]
#        csv_format=seedfile.split('.')[-1]=='csv'
#        if csv_format:
#            seeds = pd.read_csv(seedfile)
#        else:
#            seeds = pd.read_table(seedfile)
#    
#    assert len(seeds)>0, 'No seeds found in the file %s' % seedfile
#    seeds['skydir'] = map(SkyDir, seeds.ra, seeds.dec)
#    seeds['hpindex'] = map( Band(12).index, seeds.skydir)
#    inside = seeds.hpindex==Band(12).index(self.roi_dir)
#    seednames= seeds.name[inside]
#    if sum(inside)==0:
#        print 'no seeds in ROI'
#        return
#    if prefix is None:
#        prefix = seeds.name[0][:4]
#    srclist = []
#    for i,s in seeds[inside].iterrows():
#        try:
#            srclist.append(self.add_source(sources.PointSource(name=s['name'], skydir=s['skydir'], model=model)))
#            print 'added %s at %s' % (s['name'], s['skydir'])
#        except roimodel.ROImodelException:
#            srclist.append(self.get_source(s['name']))
#            print 'updating existing %s at %s ' %(s['name'], s['skydir'])
#    # Fit only fluxes for each seed first
#    parnames = self.sources.parameter_names
#    seednorms = np.arange(len(parnames))[np.array([s.startswith(prefix) and s.endswith('_Norm') for s in parnames])]
#    assert len(seednorms)>0, 'Did not find any seeds.'
#    try:
#        self.fit(seednorms, tolerance=0.2, ignore_exception=True)
#    except Exception, msg:
#        print 'Failed to fit seed norms %s' %msg
#        return
#    # now fit all parameter for each one 
#    deleted = []
#    for s in srclist:
#        self.fit(s.name, tolerance=0, ignore_exception=True)
#        s.ts = self.TS()
#        print '  TS = %.1f' % s.ts
#        if s.ts<5:
#            print ' TS<5, removing from ROI'
#            self.del_source(s.name)
#            deleted.append(s.name)
#        else:
#            # one iteration of pivot change
#            self.repivot([s])
#            
#    # now localize each, moving to new position
#    localization.localize_all(self, prefix=prefix, tsmap_dir=tsmap_dir, associator=associator, update=True, tsmin=10)
#    # fit full ROI 
#    self.fit(ignore_exception=True, tolerance=0.2)
#    # save pickled source object for each seed
#    for s in srclist:
#        if s.name in deleted: continue
#        sfile = os.path.join(seedcheck_dir, s.name+'.pickle')
#        pickle.dump(s, open(sfile, 'w'))
#        print 'wrote file %s' % sfile
            
