#  setup for point fit test
# $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointfit_setup.py,v 1.3 2007/07/19 14:45:12 burnett Exp $

# data selection parameters

radius = 7.0   # radius in degrees for initial data selection
event_type = -1 # 0, select front only; -1 no selection
source_id =-1  # -1: all sources -- select according to Monte Carlo source id, if present


# HEALpix level range for energy band  fits

minlevel=8   # minimum level to use for fits (>-6)  
maxlevel=13  # maximum level for fits  (<=13)
minalpha=0.15# minimum value for signal fraction to use in TS total

# parameters governing iteration cycle of bands/position

TSmin= 5     # minimum TS value to allow during iteration
skip1=1      # inital number of layers to skip in localization fit
skip2=4      # don't skip beyond this
itermax=1    # maximum number of iterations

verbose = 0  # set non-zero to get lots of output

# diffuse input file image file
diffusefile = ''
import os
if 'EXTFILESSYS' in os.environ:
    diffusefile = os.path.join(os.environ['EXTFILESSYS'],'galdiffuse','GP_gamma_v0r0p1.fits')
    diffusefile = os.path.join(os.environ['EXTFILESSYS'],'galdiffuse','GP_gamma.fits')

#  specify files with FT1 data and points to fit. pixelfile for PhotonMap, files for FT1 or merit
# files = [] 
pixelfile = r'F:\glast\data\SC2\obssim\allsky_noGRBs.fits'
#files = [r'F:\glast\data\octobertest\merit\r0252672900_t0213642152_merit.root']

# near gal equator, confused by diffuse
name = ['UW_J266p622', 'UW_J3587p604']
ra   = [  266.55  , 358.662]
dec  = [ 31.534 , 60.4437]

# high latitude, strong local source
name = ['UW_J83m142', 'UW_J78m151']
ra   = [ 8.3398,   7.8164]
dec  = [-14.225, -15.1306]

# the troublesome triplet: the 3EG is very strong, affects the HLCloud
name = ['DC2_3EGJ1635m1751', 'HLCloud_SC1_05']#, 'Giommi_blazar_1237']
ra   = [248.788,    248.4804] #, 248.34]
dec  = [-17.861,   -18.294] #,  -18.71]
verbose = 1  # set non-zero to get lots of output
 
factor=40
exposure = factor* 3.16e7 * 1e4 
print 'using exposure factor %f' %factor