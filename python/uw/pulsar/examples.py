"""
Code to provide examples and also basic tests.

$Header:  $

author: M. Kerr <matthew.kerr@gmail.com>

"""

# demonstrate simulation and fitting of light curves with primitive components

import lcprimitives as lp; reload(lp)
import lcfitters as lf; reload(lf)

# ----------------------------------------------------------------------------#
# a simple light curve with 2 (wrapped) gaussian peaks
# ----------------------------------------------------------------------------#

# make objects for the two gaussians -- parameters are normalization,
# position, and width (standard deviation)
# NB there is an "unpulsed" component equal to 20% of the total emission
g1 = lp.LCGaussian(p=[0.4,0.05,0.05])
g2 = lp.LCGaussian(p=[0.4,0.02,0.55])

# make a light curve template with a profile consisting of the two peaks
lct = lf.LCTemplate(primitives=[g1,g2])

# simulate 200 photons from this light curve for fitting
photons = lct.random(200)

# now, do a maximum likelihood fit for the parameters
lcf = lf.LCFitter(lct,photons)
print ('Performing ML fit -- may take a few moments.')
lcf.fit()

# print the results
print (lcf)

## now, compare with the "quick" (chi^2) fit
#print ('Doing a quick fit -- note no errors are reported yet.')
#lcf.quick_fit()
#print lcf

# get the FWHM of the peaks
print (g1.fwhm(),g2.fwhm())

