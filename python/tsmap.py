"""
Special function to make TS maps
--------------------------------

Currently coded as a function, but which could be added as a member function to roi_analysis.ROIAnalysis
$Header:$

"""
import numpy as np
import pylab as pl
import os,  math

# imports from pointlike
import  roi_localize 
from pypsf     import Psf,OldPsf,NewPsf

# from skymaps
from skymaps import SkyDir,  PySkyFunction
import image

def plot_tsmap(self, name=None, center=None, size=0.5, pixelsize=None, outdir=None, 
        which=0, catsig=99, axes=None, fignum=99, 
        bandfits=True,
        galmap=True, galactic=False,
        ): 
    """ create a TS map for the source

    Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  name        [None]  -- provide name for title, and to save figure if outdir set
  center      [None] -- center, default the roi center
  outdir      [None] -- folder name to save the image in, outdir/<name>_tsmaps.png
  catsig      [99]  -- if set and less than 1.0, draw cross with this size (degrees)
  size        [0.5]  -- half width=height (deg)
  pixelsize   [None] -- if not set, will be 20 x20 pixels
  galmap      [True] -- if set, draw a galactic coordinate image with the source position shown
  which       [0]    -- chose a different source in the ROI to plot
  =========   =======================================================

    returns the image.TSplot object for plotting positions, for example
    """
    roi = self

    name = self.psm.point_sources[0].name if name is None else name
    self.localizer = roi_localize.ROILocalizer(self, which, bandfits)

    tsm = PySkyFunction(self.localizer)
    sdir = center if center is not None else self.sa.roi_dir
    if axes is None: 
        pl.figure(fignum,figsize=(5,5)); pl.clf()
        axes= pl.axes([0.18, 0.1, 0.8,0.8])
    
    tsp = image.TSplot(tsm, sdir, size, pixelsize =pixelsize if pixelsize is not None else size/20. , 
                axes=axes, galactic=galactic)
    if 'qform' in roi.__dict__ and roi.qform is not None:
        sigma = math.sqrt(roi.qform.par[3]*roi.qform.par[4]) # why do I need this?
        qual = roi.qform.par[6]
        if sigma<1 and qual <20:
            tsp.overplot(roi.qform, sigma)
        else:
            print 'bad fit sigma %g, >1 or qual %.1f >20' % (sigma, qual)
    tsp.show(colorbar=False)
    if catsig<1:
        tsp.cross(sdir, catsig, lw=2)
    if galmap:
        axi = pl.gcf().add_axes((0.75, 0.80, 0.20, 0.10))
        ait_insert=image.AIT_grid(axes=axi, labels=False, color='w')
        ait_insert.plot([sdir], 'sr')

    if 'tsmax' in self.__dict__ and self.tsmax is not None:
        tsp.plot(self.tsmax, symbol='*')
    pl.figtext(0.5,0.92, name, fontsize=12, ha='center')
    if outdir is not None: pl.savefig(os.path.join(outdir,'%s_tsmap.png'%name.strip()))
    return tsp
