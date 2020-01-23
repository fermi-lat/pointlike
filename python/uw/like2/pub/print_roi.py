"""
Implementation of various roi printing
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pub/print_roi.py,v 1.1 2011/02/11 21:27:34 burnett Exp $
"""
import os, math
import numpy as np

def print_summary(roi, sdir=None, galactic=False, maxdist=5, title=None):
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
    if sdir is None: sdir = self.roi_dir
    if title is None: title = self.name
    print (90*'-', '\n\t Nearby point sources within %.1f degrees %s' % (maxdist,title))
    colstring = 'name dist ra dec TS flux8 index beta cutoff'
    if galactic: colstring =colstring.replace('ra dec', 'l b')
    colnames = tuple(colstring.split())
    n = len(colnames)-1
    print (('%-13s'+n*'%10s')% colnames)
    point    = list(self.psm.point_sources)
    extended = [s for s in self.bgm.diffuse_sources if 'spatial_model' in s.__dict__ ]
    for ps in point+extended:
        sdir = ps.skydir
        model = roi.get_model(ps.name)
        dist=math.degrees(sdir.difference(self.roi_dir))
        if maxdist and dist>maxdist:  continue
        loc = (sdir.l(),sdir.b()) if galactic else (sdir.ra(),sdir.dec())
        par, sigpar= model.statistical()
        expcutoff = model.name=='ExpCutoff'
        npar = len(par)
        ts = '%10.0f'% self.TS(which=ps.name) if np.any(model.free) else 10*' '
        fmt = '%-18s%5.1f'+2*'%10.3f'+ '%10s'+ '%10.2f%1s'
        freeflag = map(makefreeflag, model.free, sigpar)
        values = (ps.name.strip(), dist) +loc+ (ts,)+( model.fast_iflux()/1e-8, freeflag[0], )
        for i in range(1,npar): # parameters beyond flux
            if expcutoff and i==npar-1: fmt+=10*' '# gap if ExpCutoff to line up with cutoff 
            fmt    += '%9.2f%1s' 
            values += (par[i], freeflag[i]) 
        print (fmt % values)
        
    print (90*'-','\n\tDiffuse sources\n',90*'-')
    for source in self.bgm.diffuse_sources:
        if  'spatial_model' in source.__dict__: continue
        par, sigpar = source.model.statistical()
        n= len(par)
        freeflag = map(makefreeflag, source.model.free, sigpar)
        fmt ='%-22s' 
        values = (source.name.strip(),)
        for v,f in zip(par, freeflag):
            fmt +='%10.2f%1s'
            values +=(v,f)
        print (fmt % values)
    print (90*'-')
