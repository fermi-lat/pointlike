"""
Implementation of various roi printing
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_printing.py,v 1.10 2012/10/04 21:21:07 lande Exp $
"""
import os, math
import numpy as N

from . pointspec_helpers import PointSource

def print_summary(roi, sdir=None, galactic=False, maxdist=5, title=None, print_all_ts=False, indent=''):
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
    if title is None: 
        title = self.name if hasattr(self,'name') else ''
    print (indent,90*'-', '\n\t Nearby sources within %.1f degrees %s' % (maxdist,title))
    colstring = 'name dist ra dec TS flux8 index beta cutoff'
    if galactic: colstring =colstring.replace('ra dec', 'l b')
    colnames = tuple(colstring.split())
    n = len(colnames)-1
    print (indent,('%-13s'+n*'%10s')% colnames)
    sources = self.get_sources()
    sources.sort(key=lambda s:s.skydir.difference(self.roi_dir))
    for ps in sources:
        sdir = ps.skydir
        model = roi.get_model(ps.name)
        dist=math.degrees(sdir.difference(self.roi_dir))
        if maxdist and dist>maxdist:  continue
        loc = (sdir.l(),sdir.b()) if galactic else (sdir.ra(),sdir.dec())
        par, sigpar= model.statistical()
        expcutoff = model.name=='ExpCutoff'
        npar = len(par)
        ts = '%10.0f'% self.TS(which=ps.name) if (N.any(model.free) or print_all_ts) else 10*' '
        fmt = '%-18s%5.1f'+2*'%10.3f'+ '%10s'+ '%10.2f%1s'
        freeflag = map(makefreeflag, model.free, sigpar)
        values = (ps.name.strip(), dist) +loc+ (ts,)+( model.fast_iflux()/1e-8, freeflag[0], )
        for i in range(1,npar): # parameters beyond flux
            if expcutoff and i==npar-1: fmt+=10*' '# gap if ExpCutoff to line up with cutoff 
            fmt    += '%9.2f%1s' 
            values += (par[i], freeflag[i]) 
        print (indent,fmt % values)
        
    print (indent,90*'-')
    print (indent,'\tDiffuse sources')
    print (indent,90*'-')
    for source in self.bgm.diffuse_sources:
        if  'spatial_model' in source.__dict__: continue
        par, sigpar = source.smodel.statistical()
        n= len(par)
        freeflag = map(makefreeflag, source.smodel.free, sigpar)
        fmt ='%-22s' 
        values = (source.name.strip(),)
        for v,f in zip(par, freeflag):
            fmt +='%10.2f%1s'
            values +=(v,f)
        print (indent,fmt % values)
    print (indent,90*'-')
    print (indent,'logLikelihood = ',-roi.logLikelihood(roi.parameters()))
    print (indent,90*'-')

def print_resids(roi):
    """Print out (weighted) residuals for each energy range, both in
       separate front/back columns and in a joint column.

       Useful for identifying systematic effects that distinguish between
       front and back events.
    """
    self=roi

    d = dict()
    for b in self.bands:
        key = (-1 if b.ct==1 else 1)*int(b.e)
        d[key] = b
    ens = N.sort(list(set([b.e for b in self.bands]))).astype(int)
    print ('')
    print ('        -------CT=0--------     -------CT=1--------     ------CT=0+1-------')
    print ('Energy  Mod     Obs     Res     Mod     Obs     Res     Mod     Obs     Res')
    print ('        -------------------     -------------------     -------------------')
    for en in ens:
        s1 = '%-6.0f'%(en)
        tm = 0; to = 0
        for key in [en,-en]:
            if key in d.keys():
                b  = d[key]
                m = b.ps_all_counts + b.bg_all_counts
                o = b.photons
            else:
                m = o = 0
            tm += m; to += o
            wres = (o-m)/m**0.5 if m>0 else 0
            s1 = '\t'.join([s1,'%-6.0f\t%-6d\t%.1f'%(m,o,wres)])
        wres = (to-tm)/tm**0.5 if tm>0 else 0
        s1 = '\t'.join([s1,'%-6.0f\t%-6d\t%.1f'%(tm,to,(to-tm)/tm**0.5)])
        print (s1)

def printSpectrum(roi,sources=None):
    """Print total counts and estimated signal in each band for a list of sources.

    Sources can be specified as PointSource objects, source names, or integers
    to be interpreted as indices for the list of point sources in the roi. If
    only one source is desired, it needn't be specified as a list. If no sources
    are specified, all sources with free fit parameters will be used."""

    self=roi

    if sources is None:
        sources = [s for s in self.psm.point_sources if N.any(s.model.free)]
    elif type(sources) != type([]):
        sources = [sources]
    if sources == []: return # No point sources in ROI
    bad_sources = []
    for i,s in enumerate(sources):
        if type(s) == PointSource:
            if not s in self.psm.point_sources:
                print ('Source not found in source list:\n%s\n'%s)
                bad_sources += [s]
        elif type(s) == int:
            try:
                sources[i] = self.psm.point_sources[s]
            except IndexError:
                print ('No source #%i. Only %i source(s) specified.'\
                        %(s,len(self.psm.point_sources)))
                bad_sources += [s]
        elif type(s) == type(''):
            names = [ps.name for ps in self.psm.point_sources]
            try:
                sources[i] = self.psm.point_sources[names.index(s)]
            except ValueError:
                print ('No source named %s'%s)
                bad_sources += [s]
        else:
            print ('Unrecognized source specification:', s)
            bad_sources += [s]
    sources = set([s for s in sources if not s in bad_sources])
    indices = [list(self.psm.point_sources).index(s) for s in sources]
    self.setup_energy_bands()

    fields = ['  Emin',' f_ROI',' b_ROI' ,' Events','Galactic','Isotropic']\
                 +[' '*15+'Signal']*len(sources)
    outstring = 'Spectra of sources in ROI about %s at ra = %.2f, dec = %.2f\n'\
                      %(self.psm.point_sources[0].name, self.roi_dir.ra(), self.roi_dir.dec())
    outstring += ' '*54+'  '.join(['%21s'%s.name for s in sources])+'\n'
    outstring += '  '.join(fields)+'\n'
    print (outstring)
    for eb in self.energy_bands:
        print (eb.spectralString(which=indices))
