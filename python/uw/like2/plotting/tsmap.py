"""
Code to plot TS maps

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/plotting/tsmap.py,v 1.1 2011/09/18 17:43:57 burnett Exp $

"""
import math, os
import numpy as np
from uw.utilities import image
import pylab as plt

def plot(localizer, name=None, center=None, size=0.5, pixelsize=None, outdir=None, 
        which=0, catsig=99, axes=None, fignum=99, 
        bandfits=True,
        galmap=True, galactic=False,
        assoc = None,
        notitle = False,
        nolegend = False,
        markercolor='blue', markersize=12,
        primary_markercolor='green', primary_markersize=14,
         **kwargs):
    """ create a TS map for the source. These are localization style
        TS maps (where the source is not in the background model) and
        are useful for source localization.

    localizer: a Localization object

        Optional keyword arguments:

      =========   =======================================================
      Keyword     Description
      =========   =======================================================
      name        [None]  -- provide name for title, and to save figure if outdir set
      center      [None] -- center, default the location
      outdir       [None] if set, save sed into <outdir>/<source_name>_tsmap.png if outdir is a directory, save into filename=<outdir> if not.  
      catsig      [99]  -- if set and less than 1.0, draw cross with this size (degrees)
      size        [0.5]  -- width=height (deg)
      pixelsize   [None] -- if not set, will be 20 x20 pixels
      galmap      [True] -- if set, draw a galactic coordinate inset image with the source position shown
      galactic    [False] -- plot using galactic coordinates
      assoc       [None] -- if set, a list of tuple of associated sources 
      notitle     [False] -- set to turn off (allows setting the current Axes object title)
      nolegend    [False]
      markersize  [12]   -- set 0 to not plot nearby sources in the model
      markercolor [blue]
      =========   =======================================================

    returns the image.TSplot object for plotting positions, for example
    """
    if name is None: name=localizer.source.name
    sdir = center if center is not None else localizer.skydir
    if axes is None: 
        plt.figure(fignum,figsize=(5,5)); plt.clf()
    
    tsp = image.TSplot(localizer.TSmap, sdir, size, 
                pixelsize=pixelsize if pixelsize is not None else size/20. , 
                axes=axes, galactic=galactic, galmap=galmap, **kwargs)
    if hasattr(localizer, 'ellipse'): 
        loc = localizer.ellipse
        sigma = np.sqrt(loc['a']*loc['b']) #?? scale factor needed?
        qual = loc['qual']
        if sigma<1 and qual <50:
            #tsp.overplot([loc['ra'],loc['dec'],loc['a'],loc['b'],loc['ang'], qual], sigma)
            tsp.overplot(localizer.qform, sigma)
        else:
            print 'bad fit sigma %g, >1 or qual %.1f >50' % (sigma, qual)
    tsp.show(colorbar=False)
    if catsig<1:
        tsp.cross(sdir, catsig, lw=2, color='grey')
        
    # plot the primary source, any nearby from the fit
    x,y = tsp.zea.pixel(sdir)
    tsp.zea.axes.plot([x],[y], '*', color=primary_markercolor, label=name, markersize=primary_markersize)
    marker = 'ov^<>1234sphH'; i=k=0
    if markersize!=0:
        pass
        # plot nearby sources in the ROI -- disable for now
        #for ps in roi.psm.point_sources: # skip 
        #    x,y = tsp.zea.pixel(ps.skydir)
        #    if ps.name==name or x<0 or x>tsp.zea.nx or y<0 or y>tsp.zea.ny: continue
        #    tsp.zea.axes.plot([x],[y], marker[k%12], color=markercolor, label=ps.name, markersize=markersize)
        #    k+=1
    
    tsp.plot(tsp.tsmaxpos, symbol='+', color='k') # at the maximum
    if not notitle: plt.title( name, fontsize=24)

    if assoc is not None:
        # eventually move this to image.TSplot
        last_loc,i=SkyDir(0,90),0
        for aname, loc, prob, catid in zip(assoc['name'],assoc['dir'],assoc['prob'],assoc['cat']):
            #print 'associate with %s, prob=%.2f' % (aname.strip(),prob)
            if catid in ('ibis',): 
                print '---skip gamma cat %s' % catid
                continue
            if i>8:
                print '---skip because too many for display'
                continue
            x,y = tsp.zea.pixel(loc)
            diff = np.degrees(loc.difference(last_loc)); last_loc=loc
            if diff>1e-3: k+=1 # new marker only if changed place
            tsp.zea.axes.plot([x], [y], marker=marker[k%12], color='green', linestyle='None',
                label='%s[%s] %.2f'%(aname.strip(), catid, prob ), markersize=markersize)
            i+=1
    
    fs = plt.rcParams['font.size']
    plt.rcParams.update({'legend.fontsize':7, 'font.size':7})
    # put legend on left.
    if not nolegend: tsp.zea.axes.legend(loc=2, numpoints=1, bbox_to_anchor=(-0.15,1.0))
    plt.rcParams['font.size'] = fs

    #if outdir is not None: 
    #  if os.path.isdir(outdir):
    #    plt.savefig(os.path.join(outdir,'%s_tsmap.png'%name.strip()))
    #  else :
    #    plt.savefig(outdir)
    #return tsp

    if hasattr(localizer,'qform') and localizer.qform is not None:
        tsp.overplot(localizer.qform)
    tsp.zea.axes.set_title('%s'% name, fontsize=16)  # big title
    if outdir is not None:
        filename = name.replace(' ','_').replace('+','p')
        fout = os.path.join(outdir, ('%s_tsmap.png'%filename) )
        plt.savefig(fout)
        print 'saved tsplot to %s' % fout 
        if kwargs.get('tsfits', False): 
            fitsname = os.path.join(outdir, '%s_tsmap.fits' % filename)
            tsp.zea.skyimage.reimage(tsm.zea.center,fitsname , pixelsize, tsize)
            print 'saved fits format to %s' % fitsname
    return tsp
