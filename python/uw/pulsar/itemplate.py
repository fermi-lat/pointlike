#!/usr/bin/env python
"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/itemplate.py,v 1.1 2011/04/27 18:32:03 kerrm Exp $

Provide a method for interactively fitting a multi-gaussian template to data.

Authors: Paul S. Ray <paul.ray@nrl.navy.mil>
         Matthew Kerr <matthew.kerr@gmail.com>
"""
import numpy as np
import pylab as pl
import pyfits
from lcprimitives import LCGaussian
from lcfitters import LCTemplate,LCFitter
from lcspeclike import light_curve
from optparse import OptionParser

def get_phases(ft1file,get_weights=False,weightcol='WEIGHT'):
    f = pyfits.open(ft1file)
    phases = np.asarray(f['EVENTS'].data.field('PULSE_PHASE'),dtype=float)
    if get_weights:
        weights = np.asarray(f['EVENTS'].data.field(weightcol),dtype=float)
    else: weights = None
    f.close()
    return phases,weights

class InteractiveFitter(object):

    def init(self):
        self.nbins = 50
        self.weights = None
        self.fignum = 1

    def welcome(self):
        print 'Welcome to the interactive unbinned template fitter!'
        print 'Displaying the profile... now, we will specify where to put Gaussians.'
        print 'For each peak, drag a horizontal line'
        print '         AT THE HIGHEST POINT ON THE PEAK'
        print '         from HALF-MAX to HALF-MAX'
        print 'After each drag, the plot will refresh with the current template.'
        print 'After all Gaussians are specified, close the plot, and fitting will start.'
        print '(Note -- if using interactively, you will start the fit with do_fit; but do close the plot!'
        
    def __init__(self,phases,**kwargs):
        self.init()
        self.__dict__.update(**kwargs)
        self.phases = phases
        self.primitives = []
        self.dom = np.linspace(0,1,100)
        self.welcome()
        pl.close(self.fignum)
        self.fig = pl.figure(self.fignum)
        self.ax  = pl.gca()
        self.connect()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax)
        pl.show()

    def do_fit(self):
        print 'Fitting the template with unbinned likelihood...'
        template = LCTemplate(self.primitives)
        fitter   = LCFitter(template,self.phases,weights=self.weights)
        fitter.fit()
        print 'Fitting finished!'
        print fitter
        print 'Overlaying fitted template...'
        self.fig = pl.figure(self.fignum)
        self.ax = pl.gca()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax,template=template)
        pl.show()
        self.fitter = fitter

    def connect(self):
        self.cidpress  = self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cidrelese = self.fig.canvas.mpl_connect('button_release_event',self.on_release)

    def on_press(self,event):
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self,event):
        x1 = event.xdata
        y1 = event.ydata

        fwhm  = x1 - self.x0
        peak  = (y1 + self.y0)/2.
        phase = (x1 + self.x0)/2.

        # just Gaussian for now
        sigma = fwhm/(8 * np.log(2))**0.5
        ampl  = peak * sigma * (2*np.pi)**0.5

        self.primitives += [LCGaussian(p=[ampl,sigma,phase])]
        template = LCTemplate(self.primitives)
        self.ax.clear()
        light_curve(self.phases,weights=self.weights,nbins=self.nbins,axes=self.ax,template=template)
        pl.draw()

    def write_template(self,outfile):
        if not hasattr(self,'fitter'):
            print 'Must do fit first!'; return
        self.fitter.write_template(outfile)


if __name__ == '__main__':

    desc="Read an FT1 file containing PULSE_PHASE info and interactively fit a template."""
    parser=OptionParser(usage=" %prog [options] [FT1_FILENAME]", description=desc)
    parser.add_option('-n','--nbins',type='int',default=50,help="Number of bins to use in phase histogram.")
    parser.add_option('-w','--weights',action='store_true',default=False,help='Use weighted light curve')
    parser.add_option('-c','--weightcol',type='string',default='WEIGHT',help='Column in FT1 file that holds the weight')
    parser.add_option('-p','--prof',type='string',default=None,help='Output name for profile')
    parser.add_option('-m','--min_weight',type='float',default=0.1,help='Minimum weight to include in fit.')
    
    ## Parse arguments
    (options,args) = parser.parse_args()

    phases,weights = get_phases(args[0],get_weights=options.weights,weightcol=options.weightcol)

    if options.weights:
        phases = phases[weights > options.min_weight]
        print '%d of %d photons survive weight cut'%(len(phases),len(weights))
        weights = weights[weights > options.min_weight]

    intf = InteractiveFitter(phases,nbins=options.nbins,weights=weights)
    intf.do_fit()

    if options.prof is not None: out = options.prof
    else:
        out = ''
        out = raw_input('Enter filename for gaussian profile output file, or just hit ENTER to exit...:  ')
    if len(out) > 0:
        print 'Writing Gaussian-style template to %s...'%(out)
        intf.write_template(out)
    print 'Goodbye!'
    
    

