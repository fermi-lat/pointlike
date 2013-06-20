"""
Application module, allowing command-line access to analysis/plotting tasks

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/app.py,v 1.3 2013/06/20 03:56:40 burnett Exp $

"""
import os, argparse, types
import numpy as np
import pylab as plt
import pandas as pd

from uw.like2.analyze import ( html, seedcheck, export, sourceinfo, )
#    associations,components, countplots, data, environment, 
#    flux_corr_iso, fluxcorr, frontbacksedplots,galactic, galacticspectra, gtlikecomparison,
#    hptables, isotropic, isotropicspectra, limb,limbrefit, localization, pgwseedcheck,
#    ptstable, pulsarseedcheck, roi_info, seedcheck,sourcecomparison, sourceinfo, sourcetotal, 
#    sunmoon, sunmoonrefit, uwsourcecomparison, diagnostics, find_peak, export)


opts = dict(
#        associations=  (associations.Associations,),
#        components= (components.Components,),
#        counts=     (countplots.CountPlots,),
#        data =      (data.Data,),
#        environment=   (environment.Environment,),
#        fluxcorriso=(flux_corr_iso.FluxCorrIso,),
#        fluxcorr=   (fluxcorr.FluxCorr,),
#        fb=         (frontbacksedplots.FrontBackSedPlots,),
#        galactic=   (galactic.Galactic,),
#        galspect =  (galacticspectra.GalacticSpectra,),
#        gtlike_comparison=(gtlikecomparison.GtlikeComparison,),
#        hptables =  (hptables.HPtables,),
#        isotropic=  (isotropic.Isotropic,),
#        isospect =  (isotropicspectra.IsotropicSpectra,),
#        limb=       (limb.Limb,),
#        limb_refit= (limbrefit.LimbRefit,),
#        localization=(localization.Localization,),
#        pgwseedcheck=(pgwseedcheck.PGWSeedCheck,),
#        pts=        (ptstable.PTStable,),
#        pseedcheck= (pulsarseedcheck.PulsarSeedCheck,),
#        roi=        (roi_info.ROIinfo,),
        seedcheck=  (seedcheck.SeedCheck,),
#        comparison= (sourcecomparison.SourceComparison,),
        sourceinfo= (sourceinfo.SourceInfo,),
#        sourcetotal=(sourcetotal.SourceTotal,),
#        sunmoon=    (sunmoon.SunMoon,),
#        sunmoon_refit = (sunmoonrefit.SunMoonRefit,),
#        uw_comparison=(uwsourcecomparison.UWsourceComparison,),
#        findpeak=   (find_peak.FindPeak,),
        export=     (export.Export,),
        ) 

        
        
def main(args, update_top=False , raise_exception=False):
    np.seterr(invalid='warn', divide='warn')
    success=True
    if type(args)==types.StringType: args = args.split()
    for arg in args:
        if arg=='all':
            cs = set(np.hstack(opts.values()))
            for cls in cs:
                if os.path.exists(cls.require):
                    print 'running %s' % cls.__name__
                    try:
                        cls('.').all_plots()
                        plt.close('all')
                    except Exception, msg:
                        print '=====failed====\n %s\n=============='% msg
                else:
                    print 'skipped %s, missing %s' % (cls.__name__, cls.require)
            break
        if arg=='menu': 
            update_top=True
            continue

        if arg not in opts.keys():
            print 'found %s; expect one of %s' % (arg, opts.keys())
            continue
            success = False
        try:
            for cls in opts[arg]:
                cls('.').all_plots()
                plt.close('all')
        except FloatingPointError, msg:
            print 'Floating point error running %s: "%s"' % (arg, msg)
            print 'seterr:', np.seterr()
            success=False
        except Exception, msg:
            print 'Exception running %s: "%s"' % (arg, msg)
            if raise_exception: raise
            success = False
    if success: 
        html.HTMLindex().create_menu()
        if update_top: html.HTMLindex().update_top()
        
    return success  
      
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='run an analysis/plotting job; must be in skymodel folder')
    parser.add_argument('args', nargs='+', help='processsor identifier: must be one of %s' %opts.keys())
    parser.add_argument('--update_top', action='store_true', help='Update the top level Web  menu')
    parser.add_argument('--raise_exception', action='store_true', help ='set to catch exceptions')
    args = parser.parse_args()
    if not main(args.args, update_top=args.update_top, raise_exception=args.raise_exception):
        raise Exception
