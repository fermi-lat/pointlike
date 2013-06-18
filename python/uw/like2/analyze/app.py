"""
Application module, allowing command-line access to analysis/plotting tasks

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/app.py,v 1.1 2013/06/17 21:48:45 burnett Exp $

"""

from uw.like2.pipeline import diagnostic_plots as dp 
from uw.like2.analyze import (find_peak, export )
import os, argparse, types
import numpy as np
import pylab as plt
import pandas as pd

opts = dict(
        environment=   (dp.Environment,),
        counts=  (dp.CountPlots,),
        sources= (dp.SourceInfo, dp.Localization, ), #SourceTotal,),
        sourceinfo=(dp.SourceInfo,),
        localization=(dp.Localization,),
        diffuse= (dp.Galactic, dp.Isotropic, dp.Limb, dp.SunMoon),
        isotropic=(dp.Isotropic,),
        galactic=(dp.Galactic,),
        limb=    (dp.Limb,),
        limb_refit=(dp.LimbRefit,),
        sunmoon= (dp.SunMoon,),
        sunmoon_refit = (dp.SunMoonRefit,),
        isospect =  (dp.IsotropicSpectra,),
        galspect =  (dp.GalacticSpectra,),
        fb=      (dp.FrontBackSedPlots,),
        fluxcorr=(dp.FluxCorr,),
        fluxcorriso=(dp.FluxCorrIso,),
        loc =    (dp.Localization,),
        loc1K =  (dp.Localization1K,),
        hptables = (dp.HPtables,),
        tables = (dp.HPtables,),
        sourcetotal=(dp.SourceTotal,),
        seedcheck=(dp.SeedCheck,),
        pseedcheck=(dp.PulsarSeedCheck,),
        pts=     (dp.PTStable,),
        data =   (dp.Data,),
        comparison=(dp.SourceComparison,),
        association=(dp.Associations,),
        gtlike_comparison=(dp.GtlikeComparison,),
        uw_comparison=(dp.UWsourceComparison,),
        findpeak= (find_peak.FindPeak,),
        export= (export.Export,),
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
        dp.HTMLindex().create_menu()
        if update_top: dp.HTMLindex().update_top()
        
    return success  
      
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='run an analysis/plotting job; must be in skymodel folder')
    parser.add_argument('args', nargs='+', help='processsor identifier: must be one of %s' %opts.keys())
    parser.add_argument('--update_top', action='store_true', help='Update the top level Web  menu')
    args = parser.parse_args()
    if not main(args.args, update_top=args.update_top):
        raise Exception
    
