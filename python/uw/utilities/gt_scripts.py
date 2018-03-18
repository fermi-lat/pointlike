""" Scripts for interfacing with the ScienceTools.

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/gt_scripts.py,v 1.3 2013/10/24 10:16:19 kerrm Exp $

    author: Matthew Kerr
"""

def get_coords(ft1):
    """Try to get ROI center from the FITS header."""
    possible_keys = ['DSVAL%d'%x for x in range(1,5)]
    for key in possible_keys:
        try:
            from pyfits import getheader
            val = getheader(ft1,'EVENTS').get(key)
            ra,dec,orig_rad = val.split('(')[-1].split(')')[0].split(',') # note strings
            return ra,dec,orig_rad
        except Exception:
            pass
    return [None]*3

def do_gtselect(ft1,outfile=None,ra=None,dec=None,rad=180,tmin=0,tmax=0,
                emin=100,emax=1e5,zmax=100,opt_strings=[]):
    if ((ra is None) or (dec is None)) and (rad != 180):
        ra,dec,orig_rad = get_coords(ft1)
    if ra is None and (rad != 180):
        raise Exception('Unable to determine RA/DEC from FT1 file.  Aborting.')
    if outfile is None: outfile = ft1[:-5] + '_gtselect.fits'
    command = ' '.join(['gtselect',
                        'infile=%s'%(ft1),
                        'outfile=%s'%(outfile),
                        'ra=%s'%(ra),
                        'dec=%s'%(dec),
                        'rad=%s'%(rad),
                        'tmin=%s'%(tmin),
                        'tmax=%s'%(tmax),
                        'emin=%s'%(emin),
                        'emax=%s'%(emax),
                        'zmax=%s'%(zmax),
                        ] + opt_strings)
    return command            
            
def do_gtmktime(ft1,ft2,filter_expr=None,outfile=None,roicut='yes',opt_strings=[]):
    if outfile is None: outfile = ft1[:-5] + '_gtmktime.fits'
    if filter_expr is None:
        filter_expr = '&&'.join(['(DATA_QUAL==1)',
                                 '(LAT_CONFIG==1)',
                                 '((ABS(ROCK_ANGLE)<52 && START > 273628805) || (ABS(ROCK_ANGLE)<43 && START < 273628805))',
                                ])
    command = ' '.join(['gtmktime',
                        'evfile=%s'%(ft1),
                        'scfile=%s'%(ft2),
                        'outfile=%s'%(outfile),
                        'roicut=%s'%(roicut),
                        'filter="(%s)"'%(filter_expr),
                        ] + opt_strings)
    return command

def do_gtltcube(ft1,ft2,outfile=None,dcostheta=0.025,pixelsize=1,opt_strings=[]):
    if outfile is None:
        outfile = ft1[:-5] + '-LTCUBE.fits'
    command = ' '.join(['gtltcube',
                        'evfile=%s'%(ft1),
                        'scfile=%s'%(ft2),
                        'outfile=%s'%(outfile),
                        'dcostheta=%s'%(str(dcostheta)),
                        'binsz=%s'%(str(pixelsize)),
                        ] + opt_strings)
    return command

def do_gtbin(ft1,ft2,outfile=None,algorithm='CCUBE',pixelsize=0.1,numpix=200,
             ra=None,dec=None,emin=100,emax=100000,enumbins=24,
             opt_strings=[]):
    if outfile is None:
        outfile = ft1[:-5] + '-%s.fits'%(algorithm)
    if (ra is None) or (dec is None):
        ra,dec,orig_rad = get_coords(ft1)
    if ra is None:
        raise Exception('Unable to determine RA/DEC from FT1 file.  Aborting.')
    command = ' '.join(['gtbin',
                        'evfile=%s'%(ft1),
                        'outfile=%s'%(outfile),
                        'scfile=%s'%(ft2),
                        'algorithm=%s'%(algorithm),
                        'nxpix=%d'%(numpix),
                        'nypix=%d'%(numpix),
                        'binsz=%s'%(str(pixelsize)),
                        'coordsys=CEL',
                        'xref=%s'%(ra),
                        'yref=%s'%(dec),
                        'proj=ZEA',
                        'axisrot=0',
                        'ebinalg=LOG',
                        'emin=%s'%(emin),
                        'emax=%s'%(emax),
                        'enumbins=%d'%(enumbins),
                        ] + opt_strings)
    return command

def do_gtexpcube(ltcube,ccube,ft1,outfile=None,irf='P6_V3_DIFFUSE',opt_strings=[]):
    if outfile is None:
        outfile = ft1[:-5] + '-EXPCUBE.fits'
    command = ' '.join(['gtexpcube',
                        'infile=%s'%(ltcube),
                        'cmfile=%s'%(ccube),
                        'outfile=%s'%(outfile),
                        'irfs=%s'%(irf),
                        'evfile=%s'%(ft1),
                        'bincalc="EDGE"',
                        ] + opt_strings)
    return command

def do_gtexpcube2(ltcube,ccube,ft1,outfile=None,irf='P6_V3_DIFFUSE',opt_strings=[]):
    if outfile is None: outfile = ft1[:-5] + '-EXPCUBE.fits'
    command = ' '.join(['gtexpcube2',
                        'infile=%s'%(ltcube),
                        'cmap=%s'%(ccube),
                        'outfile=%s'%(outfile),
                        'irfs=%s'%(irf),
                        ] + opt_strings)
    return command

def do_gtexpcube2_new(ltcube,ft1,outfile=None,
                      ra=None,dec=None,pixelsize=0.1,numpix=200,
                      emin=100,emax=100000,enumbins=25,
                      irf='P6_V3_DIFFUSE',opt_strings=[]):
    if outfile is None: outfile = ft1[:-5] + '-EXPCUBE.fits'
    command = ' '.join(['gtexpcube2',
                        'infile=%s'%(ltcube),
                        'cmap=none',
                        'outfile=%s'%(outfile),
                        'irfs=%s'%(irf),
                        'xref=%s'%(ra),
                        'yref=%s'%(dec),
                        'coordsys=CEL',
                        'axisrot=0',
                        'nxpix=%d'%(numpix),
                        'nypix=%d'%(numpix),
                        'binsz=%s'%(pixelsize),
                        'proj=ZEA',
                        'emin=%s'%(emin),
                        'emax=%s'%(emax),
                        'enumbins=%s'%(enumbins)
                        ] + opt_strings)
    return command

def do_gtsrcmaps(xml,ft2,ltcube,ccube,binned_expmap,outfile=None,
                 irf='P6_V3_DIFFUSE',do_ptsrcs=False,no_outer=False):
    if outfile is None:
        outfile = ccube[:-11] + '-SRCMAPS.fits' # fragile
    command = ' '.join(['gtsrcmaps',
                        'scfile=%s'%(ft2),
                        'expcube=%s'%(ltcube),
                        'cmap=%s'%(ccube),
                        'srcmdl=%s'%(xml),
                        'bexpmap=%s'%(binned_expmap),
                        'outfile=%s'%(outfile),
                        'irfs=%s'%(irf),
                        'ptsrc=%s'%('yes' if do_ptsrcs else 'no'),
                        'emapbnds=%s'%('no' if no_outer else 'yes')
                        ])
    return command

def do_gtmodel(outfile,xml,srcmaps,ltcube,binned_expmap,irf='P6_V3_DIFFUSE'):
    command = ' '.join(['gtmodel',
                        'srcmaps=%s'%(srcmaps),
                        'srcmdl=%s'%(xml),
                        'outfile=%s'%(outfile),
                        'irfs=%s'%(irf),
                        'expcube=%s'%(ltcube),
                        'bexpmap=%s'%(binned_expmap),
                        ])
    return command

def do_gtsrcprob(ft1,ft2,xml,outfile=None,irf='P7SOURCE_V6',srclist=None):
    if outfile is None: outfile = ft1[:-5]+'_gtsrcprob.fits'
    command = ' '.join(['gtsrcprob',
                        'evfile=%s'%(ft1),
                        'scfile=%s'%(ft2),
                        'outfile=%s'%(outfile),
                        'srcmdl=%s'%(xml),
                        'irfs=%s'%(irf),
                        'srclist=%s'%(srclist or 'none')
                        ])
    return command

def do_gtdiffrsp(ft1,ft2,xml,irf='P7SOURCE_V6'):
    command = ' '.join(['gtdiffrsp',
                        'evfile=%s'%(ft1),
                        'scfile=%s'%(ft2),
                        'srcmdl=%s'%(xml),
                        'irfs=%s'%(irf),
                        ])
    return command

def get_binned_obs(ltcube,binned_expmap,irf='P6_V3_DIFFUSE',srcmaps=None):
    import BinnedAnalysis
    bo = BinnedAnalysis.BinnedObs(srcMaps=srcmaps,expCube=ltcube,
                                  binnedExpMap=binned_expmap,irfs=irf)
    return bo

# N.B. -- for unbinned analysis
def make_resids(cmap,mmap,outfile=None,weighted=False,renormalize=False):
    if outfile is None:
        outfile = cmap[:-5]+'-RMAP.fits'
    import pyfits
    f1 = pyfits.open(cmap)
    f2 = pyfits.open(mmap)
    ocounts = f1[0].data.sum()
    mcounts = f2[0].data.sum()
    print 'Total Obs Counts: %.1f'%(ocounts)
    print 'Total Mod Counts: %.1f'%(mcounts)
    scale = float(ocounts-mcounts)/mcounts
    print 'Relative Difference (obs - mod)/mod: %.3f'%( scale )
    f1[0].data = f1[0].data - f2[0].data*float(ocounts)/mcounts
    if renormalize:
        print 'Renormalizing model counts for mean 0 residuals...(multiplying model by %.3f)'%(float(ocounts)/mcounts)
    f1.writeto(outfile,clobber=True)

def do_gtbary(ft1,ft2,outfile=None,ra=None,dec=None,tcorrect='BARY'):
    if ((ra is None) or (dec is None)):
        ra,dec,orig_rad = get_coords(ft1)
    if ra is None:
        raise Exception('Unable to determine RA/DEC from FT1 file.  Aborting.')
    if outfile is None: outfile = ft1[:-5] + '_gtbary.fits'
    command = ' '.join(['gtbary',
                        'evfile=%s'%(ft1),
                        'scfile=%s'%(ft2),
                        'outfile=%s'%(outfile),
                        'ra=%s'%(ra),
                        'dec=%s'%(dec),
                        'tcorrect=\'%s\''%(tcorrect)
                       ])
    return command

