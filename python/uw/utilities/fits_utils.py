from __future__ import division
import os, sys
import astropy.io.fits as pyfits

"""
NAME: fits_utils.py

VERSION: 1.0

AUTHORS: Tyrel, Aous, Damien

  downloadFT1 : download a FT1 from the astroserver

  get_header_position = return the Position keyword in the FT1 header and
                        extract the ra, dec, and radius values

  get_header_erange = return the energy range

"""

def downloadFT1( ft1name, emin=100, emax=100000, tmin=0., tmax=0., ra=0., dec=0., rad=1.,
                 data_version='P6_public_v2', zmax=100., evclsmin=0, evtclass='Diffuse' ):
    '''Download FT1 file from the astroserver'''
    cmd = '~glast/astroserver/prod/astro -b -q --output-ft1 %s --event-sample %s \
--minEnergy %.2f --maxEnergy %.2f --minTimestamp %.2f --maxTimestamp %.2f --ra %.6f --dec %.6f \
--radius %.2f --maxZenith %.2f --event-class-name "%s" store' \
    %(ft1name,data_version,emin,emax,tmin,tmax,ra,dec,rad,zmax,evtclass)
    return cmd

def check_file(fname,display=False):
    '''Check if the file exists'''
    if not os.access(fname,os.F_OK):
        print "\033[1;31m The file %s does not exist => exit \033[0m" %(fname)
        sys.exit()
    if os.access(fname,os.F_OK):
        if display:
            print '\tReading file: ', fname
        return True
                              

def add_angsep_column(filename,ra=0.,dec=0.):
  '''Add a angular separation column between the photon direction and (ra,dec).
  ___arguments___:
  filename : name of your FT1 file
  ra       : Right Ascension (deg)
  dec      : Declination (deg)
  '''
  file = pyfits.open(filename)
  RA   = numpy.asarray(file['EVENTS'].data.field('RA'))
  DEC  = numpy.asarray(file['EVENTS'].data.field('DEC'))

  angsep=[]
  for it in range(len(RA)):
    angsep.append(angular_separation(ra,dec,RA[it],DEC[it]))

  colname = 'ANGSEP'
  # clobber old values
  try:
    file['EVENTS'].data.field(colname)[:] = angsep
    print 'Clobbered old ANGSEP column.'
  except:
    cols = file['EVENTS'].columns
    newcols = cols.add_col(pyfits.Column(name=colname,format='D',array=angsep))
    
    table = pyfits.new_table(newcols,header=file['EVENTS'].header)
    table.name = 'EVENTS'
    file[1] = table
    
  file.writeto(filename,clobber=True)

     
def get_header_position(filename):
  '''return the Position keyword in the FT1 header and extract the ra, dec, and radius values.
  ___arguments___: filename (name of FITS file)
  '''
  
  file   = pyfits.open(filename)
  num    = file[1].header['NDSKEYS']; header = file[1].header
  right  = 'POS(RA,DEC)'
  i=1; keynum=0
  
  while i<=num:  #this step is necessary since it is not clear that the POS key word will have the same number always
    word='DSTYP%i' %i
    test=file[1].header[word]
    if test == right:
      keynum=i; i=num
    i+=1

  if keynum == 0:  #DSKEYS start numbering at 1, if this value hasn't been updated, KEYword doesn't exist
    print 'Error: No position keyword found in fits header (assuming position is RA and DEC.  Exiting...'
    exit()
                
  keyword='DSVAL%i' %keynum

  try:
      ra,dec,rad=header[keyword].strip('circle()').split(',') #gets rid of the circle and parenthesis part and splits around the comma
      return float(ra),float(dec),float(rad)
  except ValueError:
      ra,dec,rad=header[keyword].strip('CIRCLE()').split(',') #gets rid of the circle and parenthesis part and splits around the comma
      return float(ra),float(dec),float(rad)


def get_header_erange(filename):
  '''Return the FT1 energy range from the HEADER: emin, emax.
  ___arguments___: filename (name of FITS file)
  '''
  
  file  = pyfits.open(filename)
  num   = file[1].header['NDSKEYS']; header=file[1].header
  right = 'ENERGY'
  i=1; keynum=0

  while i<=num:
    word='DSTYP%i' %i
    test=file[1].header[word]
    if test == right:
      keynum=i; i=num
    i+=1

  if keynum==0:
    print 'Error: No energy keyword found in fits header.  Exiting...'
    exit()

  keyword='DSVAL%i' %keynum
  emin,emax=header[keyword].split(':')
  return float(emin),float(emax)

def ft1_time_range( ft1filelistname ):
  '''Return the FT1 time range.'''
  ft1start = sys.maxint
  ft1stop = 0
  if ft1filelistname.startswith( '@' ):
    ft1filelist=open( ft1filelistname[1:], 'r' )
    while ft1filelist:
      ft1filename=ft1filelist.readline()
      if len(ft1filename) == 0:
        break
      ft1file=pyfits.open(ft1filename)
      ft1filestart=float(ft1file[0].header['TSTART'])
      if ft1filestart < ft1start:
        ft1start=ft1filestart
      ft1filestop=float(ft1file[0].header['TSTOP'])
      if ft1filestop > ft1stop:
        ft1stop=ft1filestop
  else:
    ft1file=pyfits.open(ft1filelistname)
    ft1filestart=float(ft1file[0].header['TSTART'])
    if ft1filestart < ft1start:
      ft1start=ft1filestart
    ft1filestop=float(ft1file[0].header['TSTOP'])
    if ft1filestop > ft1stop:
      ft1stop=ft1filestop

  ft1range = {'START':ft1start, 'STOP':ft1stop}
  return ft1range

# ========================================
# This function returns the FT2 time range
# ========================================
def ft2_time_range( ft2filename ):
  ft2file = pyfits.open( ft2filename )
  # NOTE: header information unreliable if using ftmerge, so we need to look into the data for FT2.
  #ft2start=ft2file[0].header['TSTART']
  #ft2stop=ft2file[0].header['TSTOP']
  # NOTE: using Aous's method of getting information from the data fields for FT2       
  sc_data = ft2file[1].data
  ft2start = sc_data.field('START')[0]
  ft2stop = sc_data.field('STOP')[-1]
  ft2range = {'START':ft2start, 'STOP':ft2stop}
  return ft2range

# =======================================================
# This function checks the time between FT1 and FT2 files
# =======================================================
def time_check( ft1range, ft2range ):
  ft1start = ft1range['START']
  ft1stop  = ft1range['STOP']

  ft2start = ft2range['START']
  ft2stop  = ft2range['STOP']
    
  print "\tFT1: ", ft1start, " ", ft1stop
  print "\tFT2: ", ft2start, " ", ft2stop

  if ft2start < (ft1start+100):
    if ft2stop>(ft1stop-100):
      return True
    else:
      print 'WARNING: you are losing', (ft1stop-ft2stop-100.0)/60.0, ' minutes of data at the end of your FT1 because of a short FT2 file. Override this check with -T 0'
  else:
    print 'WARNING: you are losing', (ft2start-ft1start-100.0)/60.0, ' minutes of data at the beginning of your FT1 because of a short FT2 file. Override this check with -T 0'
    return False

# ==================================================
# This function gets the photon time from a ft1 file
# ==================================================
def load_photon_mjds( ft1name ):
  # Load photon times (presumed to be barycentered) from FT1 file
  hdulist = pyfits.open(ft1name)
  ft1hdr=hdulist[1].header
  ft1dat=hdulist[1].data
  MJDREF = -1
  TIMEZERO = -1

  if ft1hdr['TIMESYS'] != 'TDB':
    print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "# !!!!!!!!! WARNING !!!!!!!!!! TIMESYS is NOT TDB!  You probably want to barycenter your photons!!!"
    print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    time.sleep(5)
    # Collect TIMEZERO and MJDREF
    try:
      TIMEZERO = ft1hdr['TIMEZERO']
    except KeyError:
      TIMEZERO = ft1hdr['TIMEZERI'] + ft1hdr['TIMEZERF']
    #print >>outfile, "# TIMEZERO = ",TIMEZERO
    try:
      MJDREF = ft1hdr['MJDREF']
    except KeyError:
      # Here I have to work around an issue where the MJDREFF key is stored
      # as a string in the header and uses the "1.234D-5" syntax for floats, which
      # is not supported by Python
      MJDREF = ft1hdr['MJDREFI'] + float(ft1hdr['MJDREFF'].replace('D','E'))
    #print >>outfile, "# MJDREF = ",MJDREF
  mjds = ft1dat.field('TIME')/86400.0 + MJDREF + TIMEZERO

  return mjds
    

