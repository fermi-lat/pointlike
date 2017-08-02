"""
Seed processing code
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/seeds.py,v 1.4 2016/10/28 21:21:51 burnett Exp $

"""
import os, sys, time, pickle, glob,  types
import numpy as np
import pandas as pd
from astropy.io import fits
from skymaps import SkyDir, Band
from uw.utilities import keyword_options
from uw.like2 import (tools, sedfuns, maps, sources, localization, roimodel,)


def read_seedfile(seedkey,  filename=None, config=None):

    model_name = os.getcwd().split('/')[-1]

    if model_name.startswith('month') and seedkey=='pgw':
        #monthly mode, need to find and load PGW analysis with rouighly equivalent months
        month=int(model_name[5:]); 
        filename='/nfs/farm/g/glast/g/catalog/transients/TBIN_%d_all_pgw.txt'% (month-1)
        assert os.path.exists(filename), 'PGWAVE file %s not found'% filename  
        try:
            seeds = pd.read_table(filename, sep=' ', skipinitialspace=True, index_col=1,
                header=None,
                names='tbin ra dec k_signif pgw_roi fgl_seed fgl_ra fgl_dec fgl_assoc'.split())
        except Exception,msg:
            raise Exception('Failed to read file %s: %s' % (filename, msg))
        names=[]
        for i,s in seeds.iterrows():
            j = int(s.name[4:6]) if s.name[6]=='_' else int(s.name[4:5])
            names.append('PGW_%02d_%03d_%02d' % (month, int(s.pgw_roi), j))
        seeds['name'] = names    
    elif model_name.startswith('month') and seedkey=='PGW':
        # monthly mode, new format PGwave, in a single FITS file
        month=int(model_name[5:]); 

        
        assert os.path.exists(filename), 'PGWAVE file {} not found'.format( filename)  
        t = fits.open(filename)
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

    elif filename is None and config is not None:
        # assume that config[seedkey] is the filename
        if seedkey in config:
            filename = config[seedkey]
        elif os.path.exists('seeds_{}.csv'.format(seedkey)):
            filename='seeds_{}.csv'.format(seedkey)
        else:
            raise Exception('seedkey {} not found in config, or filename'.format(seedkey))
        if os.path.splitext(filename)=='.fits':
            # a standard FITS catalog
            f = fits.open(os.path.expandvars(filename))
            name, ra, dec = [f[1].data.field(x) for x in 'Source_Name RAJ2000 DEJ2000'.split()]
            seeds = pd.DataFrame([name, np.array(ra,float),np.array(dec,float)],
            index='name ra dec'.split()).T
        else:
            seeds = pd.read_csv(filename)

    elif filename is not None:
        # file is cvs
        seeds = pd.read_csv(filename)
    else:
        # reading a TS seeds file
        t = glob.glob('seeds_%s*' % seedkey)
        assert len(t)==1, 'Seed file search, using key {}, failed to find one file\n\t{}'.format( seedkey,t)
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

def add_seeds(roi, seedkey, config=None,
            model='PowerLaw(1e-14, 2.2)', 
            prefix=None,
            associator=None, tsmap_dir='tsmap_fail',
            tsmin=5, lqmax=20,
            update_if_exists=False,
            **kwargs):
    """ add "seeds" from a text file the the current ROI
    
        roi : the ROI object
        seedkey : string
            Expect one of 'pgw' or 'ts' for now. Used by read_seedfile to find the list
        model : string
            model to use for the generated source, unless found in maps.table_info
        prefix : None or string
            Name to 
        associator :
        tsmap_dir
        mints : float
            minimum TS to accept for addition to the model
        lqmax : float
            maximum localization quality for tentative source
    """
    seeds = read_seedfile(seedkey, config=config)
    inside = seeds.hpindex==Band(12).index(roi.roi_dir)
    seednames=seeds.name[inside]
    if sum(inside)==0:
        print 'no seeds in ROI'
        return
    if prefix is None:
        prefix = seeds.name[0][:4]
    if seedkey in maps.table_info.keys():
        model = maps.table_info[seedkey][1]['model']
    print 'Using model {} for seed {}'.format(model, seedkey)
    # if seedkey=='hard':
    #     model='PowerLaw(1e-16, 1.7)'
    srclist = []
    for i,s in seeds[inside].iterrows():
        try:
            src=roi.add_source(sources.PointSource(name=s['name'], skydir=s['skydir'], model=model))
            if src.model.name=='LogParabola':
                roi.freeze('beta',src.name)
            srclist.append(src)
            print '%s: added at %s' % (s['name'], s['skydir'])
        except roimodel.ROImodelException, msg:
            if update_if_exists:
                srclist.append(roi.get_source(s['name']))
                print '{}: updating existing at {} '.format(s['name'], s['skydir'])
            else:
                print '{}: Fail to add "{}"'.format(s['name'], msg)
    # Fit only fluxes for each seed first
    seednames = [s.name for s in srclist]
    llbefore = roi.log_like() # for setup?

    # set normalizations with profile fits
    for src in srclist:
        prof= roi.profile(src.name, set_normalization=True)
        src.ts= prof['ts'] if prof is not None else 0
        print '\tTS={:.1f}'.format(src.ts)

    # now fit all norms at once
    seednorms = [s.name+"_Norm" for s in srclist if roi.get_source(s.name).ts>0]
    if len(seednorms)==0:
        print 'Did not find any seeds with prefix {}.'.format(prefix)
        return False
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
        elif tsmin>0:
            # one iteration of pivot change
            s = roi.get_source(sname)
            roi.repivot([s], min_ts=5)
            # and a localization: remove if fails or poor
            roi.localize(sname, update=True)
            ellipse = s.__dict__.get('ellipse', None) 
            if ellipse is None or ellipse[5]>lqmax:
                if tsmap_dir is not None and os.path.exists(tsmap_dir):
                    roi.plot_tsmap(sname, outdir=tsmap_dir)
                print '--> removing {} \n\t since did not localize'.format(s)
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
    
            
