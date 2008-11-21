"""A suite of tools for processing FITS files.

   $Header$

   author: Matthew Kerr

"""

import pyfits as pf
import numpy as N
from types import ListType
from math import cos,sin,pi

#TODO: GTI
#e.g. counts_plot(file_list,coordsys='galactic',cuts=['L  > 50','L < 100','B > -60'])
def counts_plot(ft1files,scale='log',pixels=256,coordsys='equatorial',cuts = None):
   ft1 = merge_flight_data(ft1files,cuts=cuts)
   events = ft1[1]
   if coordsys == 'equatorial':
      lon,lat = events.data.field('RA'),events.data.field('DEC')
   else:
      lon,lat = events.data.field('L'),events.data.field('B')
   img,x,y = N.histogram2d(lon,lat,bins=pixels)
   if scale == 'log':
      img = N.where(img > 0,N.log10(img),0)
   from pylab import pcolor,xlabel,ylabel
   pcolor(x,y,img.transpose())
   xlabel( ('RA'  if coordsys=='equatorial' else 'L') + ' (deg)')
   ylabel( ('DEC' if coordsys=='equatorial' else 'B') + ' (deg)')

def merge_flight_data(files, outputfile = None, cuts = None, fields = None):
   """Merge FT1 or FT2 files and make cuts on the columns."""

   handles = __get_handles__(files)

   try:
      handles[0]['EVENTS']
      table_name = 'EVENTS'
   except:
      table_name = 'SC_DATA'

   event_table = __merge_events__(handles, table_name = table_name)

   if cuts is not None:
      __arbitrary_cuts__(event_table,cuts)
      interval = [event_table.data.field('TIME').min(),event_table.data.field('TIME').max()]
   else: interval = None

   #Overwrite data in dummy table and write it to file
   handles[0][table_name].data = event_table.data
   if table_name == 'EVENTS':
      handles[0]['GTI'].data = __merge_gti__(handles,interval=interval).data
   
   if outputfile is not None: handles[0].writeto(outputfile,clobber=True)
   for x in handles: x.close()

   return handles[0]

def FT1_to_GTI(files):
   """Convert FT1 files to a single GTI object."""
   handles = __get_handles__(files)
   from skymaps import Gti
   g = Gti(files[0])
   starts,stops = __merge_gti__(handles[1:],no_table = True)
   if len(starts) == 0: return g
   for i in xrange(len(starts)):
      g.insertInterval(starts[i],stops[i])
   return g

def sum_ltcubes(files,outputfile = 'summed_ltcube.fits'):
   """Pass either a name of an ASCII file with other FITS file names or a list of filenames."""

   handles = __get_handles__(files)

   exp_table = __merge_exposures__(handles)
   gti_table = __merge_gti__(handles)

   #Overwrite data in dummy table and write it to file
   handles[0]['EXPOSURE'].data = exp_table.data
   handles[0]['GTI'].data = gti_table.data
   handles[0].writeto(outputfile,clobber=True)
   
   for x in handles: x.close()

#EVERYTHING BELOW IS AN INTERNAL CALL.

def __FITS_parse__(files):
   """Parse input and return a list of (FITS) filenames.

      files -- a glob-style wildcard, an ASCII file containing a list of filenames, a list of filenames, or a single FITS file."""
   if type(files) is ListType: return files
   try: #See if it's a FITS file
      f = pf.open(files)
      f[0]
      f.close()
      return [files]
   except:
      pass
   if files[0] == '@':
      return [line.strip() for line in file(ft1files[1:]) if len(line)>0 and line[0]!='#']
   from glob import glob
   return glob(files)        

def __get_handles__(files):
   files = __FITS_parse__(files)
   handles = [pf.open(x,memmap=1) for x in files]
   return handles

def __merge_exposures__(handles):
   """Return a FITS table of merged exposure cubes."""
   
   exposures = [x['EXPOSURE'].data.field('COSBINS') for x in handles]

   for exp in exposures:
      if exp.shape != exposures[0].shape:
         print 'Livetime cube binning inconsistent -- bailing out!'
         return

   #Sum the exposures and put in a new table
   summed_exposure = N.array(exposures).sum(axis=0)
   exp_table = pf.new_table(handles[0]['EXPOSURE'].columns,nrows=exposures[0].shape[0])
   exp_table.data.field('COSBINS')[:] = summed_exposure
   
   return exp_table

def __merge_events__(handles,table_name = 'EVENTS'):
   """Return a FITS table of merged FT1 events.  Now works for FT2 too!"""
   
   num_events = [x[table_name].data.shape[0] for x in handles]
   columns = handles[0][table_name].columns
   event_table = pf.new_table(columns,nrows=sum(num_events))
   previous_loc = 0
   for i,handle in enumerate(handles):
      for j in xrange(len(columns)):
         event_table.data.field(j)[previous_loc:previous_loc+num_events[i]] = handle[table_name].data.field(j)[:]
      previous_loc += num_events[i]
   return event_table

def __merge_gti__(handles,no_table=False,interval=None):
   """Return a FITS table of merged GTI."""

   if len(handles) == 0: return ([],[])

   #Merge the gti and sort the results
   starts = N.concatenate([x['GTI'].data.field('START') for x in handles])
   stops  = N.concatenate([x['GTI'].data.field('STOP')  for x in handles])
   sorting = N.argsort(starts)
   starts = starts[sorting]
   stops = stops[sorting]

   if interval is not None:
      mask   = N.logical_and(starts > interval[0], stops < interval[1])
      starts = N.append(interval[0],starts[mask])
      stops  = N.append(stops[mask],interval[1])

   if no_table: return (starts,stops)

   #Put GTI in new table
   gti_table = pf.new_table(handles[0]['GTI'].columns,nrows = len(starts))
   gti_table.data.field('START')[:] = starts
   gti_table.data.field('STOP')[:] = stops

   return gti_table

def __arbitrary_cuts__(events,cuts):
   """Perform the cuts provided as arguments.

      events -- a major table in FITS format such as might be returned by __merge_events__
      cuts -- a list of cuts; the final result is the intersection of all cuts.  e.g,
              ['ENERGY>100','MET > 5000','MET < 10000']
      columns -- a list of the columns involved in the cuts to save me some parsing effort"""
   
   from numarray import array,logical_and

   for cut in cuts: #put cut columns in namespace
      for name in events.columns.names:
         if name in cut:          
            exec('%s = events.data.field(\'%s\')'%(name,name))
   mask = array([True]*events.data.shape[0])
   for cut in cuts: #build mask cumulatively
      exec('mask = logical_and(mask,%s)'%cut)
   
   new_table = pf.new_table(events.columns,nrows=len(mask[mask]))
   for i in xrange(len(events.columns)):
      new_table.data.field(i)[:] = events.data.field(i)[mask]
   events.data = new_table.data
   