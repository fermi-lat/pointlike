"""
Implementation of various roi printing
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/printing.py,v 1.9 2014/01/26 20:07:56 burnett Exp $
"""
import os
import numpy as np
import pandas as pd
from skymaps import SkyDir

def print_summary(roi, sdir=None, galactic=False, maxdist=5, title=None, print_all_ts=False, no_global=False):
    """ formatted table point sources positions and parameter in the ROI, 
        followed by summary of diffuse names, parameters values.
        values are followed by a character to indicate status:
            blank: fixed
            ?      sigma=0, not fit yet
            !      relative error<0.1
            *      otherwise
    Parameters
    ----------
         sdir : SkyDir, optional, default None for center
               for center: default will be the first source
         galactic : bool, optional, default False
             set true for l,b
         maxdist : float, optional, default 5
             radius in degrees
         title: None or a string
            if None, use the name of the ROI. 

    """
    self = roi
    def makefreeflag(free, sigma, goodfit=0.1):
        if not free: return ' '
        if sigma==0: return '?'
        if sigma<goodfit: return '!'
        return '*'
    def beta_string(beta):
        return 8*'' if beta<1e-2 else '%8.2f'%beta
    if sdir is None:
        center =roi.roi_dir
    else:
        center = sdir if isinstance( sdir,SkyDir,) else SkyDir(*sdir)
    if title is None: 
        title = self.name if hasattr(self,'name') else ''
    print 90*'-', '\n\t Nearby sources within %.1f degrees %s' % (maxdist,title)
    colstring = 'name dist ra dec TS eflux(eV) index energy beta/b'
    if galactic: colstring =colstring.replace('ra dec', 'l b')
    colnames = tuple(colstring.split())
    n = len(colnames)-1
    print ('%-13s'+4*'%10s'+' '+(n-4)*'%9s')% colnames
    local_sources =  filter(lambda s: s.skydir is not None, roi.sources[:])
    local_sources.sort( key=lambda s: s.skydir.difference(center))
    global_sources = filter(lambda s: s.skydir is None, roi.sources[:]) if not no_global else []

    for ps in local_sources:
        sdir = ps.skydir
        model = roi.get_model(ps.name)
        dist=np.degrees(sdir.difference(center))
        if maxdist and dist>maxdist:  continue
        loc = (sdir.l(),sdir.b()) if galactic else (sdir.ra(),sdir.dec())
        try:
            par, sigpar= model.statistical()
        except:
            par  =model.parameters
            sigpar = np.zeros(model.len())
            
        npar = len(par)
        index_order = range(1,npar) if model.name!='LogParabola' else (1,3,2)
        expcutoff = model.name=='ExpCutoff'
        ts = '%10.0f'% self.TS(ps.name) if (np.any(model.free) or print_all_ts) else 10*' '
        fmt = '%-18s%5.1f'+2*'%10.3f'+ '%10s'+ '%10.1f%1s'
        freeflag = map(makefreeflag, model.free, sigpar)
        values = (ps.name.strip(), dist) +loc+ (ts,)+( model.i_flux(e_weight=1, emax=1e5)*1e6, freeflag[0], )
        
        #index
        fmt += '%8.2f%1s' 
        j = index_order[0]
        values += (par[j], freeflag[j])
        #energy
        j = index_order[1]
        fmt += '%8.0f%1s'
        values += (par[j], freeflag[j])
        # beta or b
        j = index_order[2]
        if (par[j]==1) or abs(par[j])<0.01:
            fmt+='%5.0f'; values += (par[j],)
        else:
            fmt+='%8.2f%1s'
            values += (par[j], freeflag[j])
        
        print fmt % values
        
    if len(global_sources) > 0 :
        print 90*'-','\n\tDiffuse sources\n',90*'-'
    for source in global_sources:
        par, sigpar = source.spectral_model.statistical()
        n= len(par)
        freeflag = map(makefreeflag, source.spectral_model.free, sigpar)
        fmt ='%-22s' 
        values = (source.name.strip(),)
        for v,f in zip(par, freeflag):
            fmt +='%10.2f%1s'
            values +=(v,f)
        print fmt % values
    print 90*'-'


def gtlike(roi, sources=None, tsmin=25):
    sourcenames = [s.name for s in roi.sources if s.skydir is not None\
        and np.any(s.spectral_model.free) and s.ts>tsmin]
        

    for source in sorted(sourcenames):
        model = roi.get_model(source)#.model
        print source
        print '   Spectrum: %s' %model.name
        for i in np.arange(model.npar)[model.free]:
            value = model.get_all_parameters()[i]
            error = np.sqrt(max(0, model.get_cov_matrix()[i,i]))
            gtname =model.gtlike_param_names()[i]
            scale = 10**(int(np.log10(value))-1) if gtname=='norm' else 1.0
            print '%16s: %7.4f %7.4f ( %.1e)' % (gtname, value/scale, error/scale, scale)
            
def band_summary(roi, event_type=0):
    snames = ['data','model','residual', 'normresid']+[s.name for s in roi.sources]
    d = dict()
    for index,b in enumerate(roi.selected):
        if b.band.event_type != event_type: continue
        energy = b.band.energy
        data = sum(b.data)
        resid = data-b.counts
        models = [data,b.counts, resid, resid/np.sqrt(b.counts) ]+[ s.counts for s in b]
        d[index/2]= models
    df = pd.DataFrame(d, index =snames).T
    return df
        
        
