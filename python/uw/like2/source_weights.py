"""
module that calculates weights for a source, given a model and a FT1 file
$Header$
"""

import os, argparse, pyfits
import numpy as np
import skymaps
from uw.like2 import main

class PulsarWeights(object):
    def __init__(self, filename,  keys = ['ENERGY','CONVERSION_TYPE','RA','DEC', 'TIME', 'PULSE_PHASE']):
        """ Set up a dataframe with selected columns
        """
        self.filename= filename
        self.fits = pyfits.open(filename)
        self.frec = self.fits[1].data # the FITS_rec 
        
    def calc_weights(self, roi, source_name, emin=100):
        """ calculate the weights, using model from the roi
        """
        self.weight = roi.weights(source_name, self.frec, emin)
        
    def add_weight_columns(self, colnames=['WEIGHT']):
        def update_or_add(cols, colname, values):
            if (colname in dnames):
                self.fits[1].data.field(colname)[:] = values
            else:
                cols += [pyfits.Column(name=colname,format='E',array=values)]
        dnames = self.fits[1].data.names
        cols = []
        update_or_add(cols, colnames[0], self.weight)
        if len(cols) > 0:
            tbhdu = pyfits.new_table(self.fits[1].columns.data+cols,header=self.fits[1].header)
            self.fits[1] = tbhdu
        self.fits.writeto(self.filename,clobber=True)
        self.fits.close()
        
def doit(modeldir, sourcename, filename, emin=100):
    modeldir = os.path.expandvars(modeldir)
    assert os.path.exists(os.path.join(modeldir, 'config.txt')), 'No "config.txt" found in %s' % modeldir
    assert os.path.exists(os.path.expandvars(filename)), 'file %s does not exist' %filename
    factory = main.factory(modeldir, quiet=True)
    roi = factory(sourcename)
    pw = PulsarWeights(filename)
    pw.calc_weights(roi, sourcename, emin=emin)
    pw.add_weight_columns()

if __name__=='__main__':
    parser = argparse.ArgumentParser( description="""Add a column the the FT1 FTIS file with the weights for the given source
    Note that modeldir must be a folder defining the pointlike all-sky model""")
    parser.add_argument('args', nargs=3, help='modeldir sourcename filename')

    parser.add_argument('--emin', default=100, help='minimum energy')
    #parser.add_argument('--cuts',  default='(sources.ts>10)',
    #       help='selection cuts')
    args = parser.parse_args()
    doit(*args.args, emin=args.emin)
