"""A suite of tools for processing FITS files.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/fitstools.py,v 1.11 2009/06/22 23:32:51 kerrm Exp $

   author: Matthew Kerr

"""

import pyfits as pf
import numpy as N
from types import ListType,FunctionType,MethodType
from math import cos,sin,pi

def rect_mask(lons,lats,cut_lon,cut_lat,lon_hwidth,lat_hwidth):
   mask = N.abs( lons - cut_lon)/N.cos(lats * N.pi / 180.) < lon_hwidth
   return N.logical_and(mask, N.abs( lats - cut_lat) < lat_hwidth)

def trap_mask(ras,decs,cut_dir,radius):
   """Make a conservative, trapezoid cut as a precursor to a radius cut."""

   c1 = N.cos( (cut_dir.dec() - radius)*N.pi/180. )
   c2 = N.cos( (cut_dir.dec() + radius)*N.pi/180. )
   ra_radius =  radius/N.minimum(c1,c2) #conservative

   mask = N.abs(ras - cut_dir.ra()) <= ra_radius
   mask = N.logical_and(mask,N.abs(decs - cut_dir.dec()) <= radius)
   return mask

def rad_mask(ras,decs,cut_dir,radius):
   """Make a slower, exact cut on radius."""
   from skymaps import SkyDir
   diffs   = N.asarray([cut_dir.difference(SkyDir(ra,dec)) for ra,dec in zip(ras,decs)])   
   mask    = diffs*(180/N.pi) < radius
   return mask,diffs[mask]

def rad_extract(eventfiles,center,radius_function,return_cols=['PULSE_PHASE'],cuts=None):
   """Extract events with a radial cut.  Return specified columns and perform additional boolean cuts.

      Return is in form of a dictionary whose keys are column names (and 'DIFFERENCES') and values are
      numpy arrays with the column values.  These will have been concatenated.

Arguments:

  =========   =======================================================
  Argument    Description
  =========   =======================================================
  eventfiles  -- a list of FT1 filenames
  center      -- a SkyDir giving the center of the radial cut
  radius_function -- can be either a float specifying a cookier cutter radial cut, or
              a function taking as arguments the energy and event_class and speciying
              the radius in degrees, e.g.

                  def radius(energy,event_class):
                     return numpy.where(event_class,2,1)*(energy/1000)**-0.75
  =========   =======================================================

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  return_cols ['RA','DEC','ENERGY','EVENT_CLASS','PULSE_PHASE'] - a list of FT1 column names to return
  cuts        None - an optional list of boolean cuts to apply, e.g., ['ENERGY > 100']
              NB -- cuts not yet implemented!!
  =========   =======================================================
   """
   if not (type(radius_function) is FunctionType or type(radius_function) is MethodType):
      simple_scalar = True
      rval = radius_function
      radius_function = lambda e,event_class: rval
   else:
      simple_scalar = False

   eventfiles = __FITS_parse__(eventfiles)

   from collections import defaultdict,deque
   coldict = defaultdict(deque)
   cols = {}
   keys = list(set(['RA','DEC','ENERGY','EVENT_CLASS','ZENITH_ANGLE']+return_cols))

   for eventfile in eventfiles:
      from pyfits import open
      e = open(eventfile,memmap=1)

      for key in keys: cols[key] = N.asarray(e['EVENTS'].data.field(key)).astype(float)

      rad   = radius_function(cols['ENERGY'],cols['EVENT_CLASS'])
      tmask = N.logical_and(trap_mask(cols['RA'],cols['DEC'],center,rad),cols['ZENITH_ANGLE'] < 105.)
      if simple_scalar:
         rmask,diffs = rad_mask(cols['RA'][tmask],cols['DEC'][tmask],center,rad)
      else:
         rmask,diffs = rad_mask(cols['RA'][tmask],cols['DEC'][tmask],center,rad[tmask])

      for key in keys: coldict[key].append(cols[key][tmask][rmask])
      coldict['DIFFERENCES'].append(diffs)
      e.close()
   
   for key in coldict.keys():
      if key == 'ZENITH_ANGLE' and 'ZENITH_ANGLE' not in return_cols: continue
      cols[key] = N.concatenate([x for x in coldict[key]])
   return cols

#TODO: GTI
#e.g. counts_plot(file_list,coordsys='galactic',cuts=['L  > 50','L < 100','B > -60'])
def counts_plot(ft1files,center,fov=10,scale='log',pixels=256,coordsys='equatorial',
                cuts = None, print_locs = False):
   
   ft1 = merge_flight_data(ft1files,cuts=cuts)
   events = ft1[1]
   if coordsys == 'equatorial':
      lon,lat = N.asarray(events.data.field('RA')),N.asarray(events.data.field('DEC'))
      clon,clat = center.ra(),center.dec()
   else:
      lon,lat = N.asarray(events.data.field('L')),N.asarray(events.data.field('B'))
      clon,clat = center.l(),center.b()
   mask = rect_mask(lon,lat,clon,clat,fov/2.,fov/2.)

   lon = lon[mask]
   lat = lat[mask]   
   
   img,x,y = N.histogram2d(lon,lat,bins=pixels)
   if scale == 'log':
      img = N.where(img > 0,N.log10(img),-1)
   from pylab import pcolor,xlabel,ylabel,imshow
 
   #pcolor(x,y,img.transpose()) #TODO -- make RA go the right way!
   imshow(img.transpose())
   xlabel( ('RA'  if coordsys=='equatorial' else 'L') + ' (deg)')
   ylabel( ('DEC' if coordsys=='equatorial' else 'B') + ' (deg)')
   if print_locs:
      if coordsys == 'equatorial': print 'RA     DEC      ENERGY      TIME     EVENT_CLASS'
      else: print 'L     B      ENERGY       TIME      EVENT_CLASS'
      en = events.data.field('ENERGY')
      time = events.data.field('TIME')
      ec = events.data.field('EVENT_CLASS')
      for i in xrange(len(lon)):
         
         print '%.2f  %.2f  %.2g  %.10g  %d'%(lon[i],lat[i],en[i],time[i],ec[i])
   return img

def merge_flight_data(files, outputfile = None, cuts = None, fields = None):
   """Merge FT1 or FT2 files and make cuts on the columns.
   
      If a cut is made on MET, the GTI are pruned and updated so that the exposure
      calculation will accurately reflect the MET cuts.

      outputfile -- [None] string argument giving output file name; if None, return the
                        pyfits Table instance
      cuts       -- [None] a list of logical cuts to apply to FITS columns, e.g. 'ENERGY > 100'
      fields     -- not implemented

      NOTA BENE: the headers will not be treated correctly!!!
   """

   handles = __get_handles__(files)
   table_name = handles[0][1].name #set to 'EVENTS' for FT1 or 'SC_DATA' for FT2

   event_table = __merge_events__(handles, table_name = table_name)

   if cuts is not None:
      __arbitrary_cuts__(event_table,cuts)
      interval = [event_table.data.field('TIME').min(),event_table.data.field('TIME').max()]
   else: interval = None

   #Overwrite data in dummy table and write it to file
   #handles[0][table_name].data = event_table.data
   #handles[0][table_name].columns = event_table.columns
   handles[0][table_name] = event_table

   if table_name == 'EVENTS':
      handles[0]['GTI'].data = __merge_gti__(handles,interval=interval).data

   if outputfile is not None: handles[0].writeto(outputfile,clobber=True)
   for x in handles: x.close()

   return handles[0]

def get_fields(files, fields, cuts = None):
   """A lightweight version to get only certain fields of flight data."""

   #from collections import defaultdict
   #data = defaultdict(list)
   data = dict()
   files = __FITS_parse__(files)

   #Process the cuts
   f = pf.open(files[0],memmap=1)
   if cuts is not None:
      cut_fields = set()
      for i,cut in enumerate(cuts): #put cut columns in namespace
         tokens = [tok.strip() for tok in cut.split()]
         for name in f['EVENTS'].columns.names:
            if name in tokens:
               cut_fields.add(name)
               cuts[i] = cuts[i].replace(name,'data[\'%s\']'%name)
      cut_fields = list(cut_fields)
   else: cut_fields = []
   f.close()

   counts = N.empty(len(files),dtype=int)

   for nfi,fi in enumerate(files):
      f = pf.open(fi,memmap=1)
      counts[nfi] = len(f['EVENTS'].data) #f['EVENTS'].data.getshape()[0]
      f.close()

   for field in fields + cut_fields:
      data[field] = N.empty(counts.sum(),dtype=N.float32)
   
   #Get fields
   counter = 0
   for fi,counts in zip(files,counts):
      f = pf.open(fi,memmap=1)
      for field in fields + cut_fields:
         #data[field] += [N.array(f['EVENTS'].data.field(field),dtype=N.float32)]
         data[field][counter:counter+counts] = f['EVENTS'].data.field(field)
      counter += counts
      f.close()

   import gc
   gc.collect()
   gc.collect()

   """
   #Concatenate results
   for field in fields + cut_fields:
      data[field] = N.concatenate(data[field]).astype(float) #note hard case
      gc.collect()
   """

   #Apply cuts
   if cuts is not None:
      mask = N.asarray([True]*len(data[data.keys()[0]]))
      for cut in cuts:
         mask = N.logical_and(mask,eval(cut))
      for field in fields:
         data[field] = data[field][mask]

   #Remove "helper" data for cuts
   for cut in cut_fields:
      if cut_fields not in fields: data.pop(cut)

   return data


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
   handles = [pf.open(x,memmap=1) for x in files] #some versions may need to disable memmap (=0)
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
   columns,header = __common_columns__(handles,table_name)
   event_table = pf.new_table(columns,header=header,nrows=sum(num_events))
   previous_loc = 0
   for i,handle in enumerate(handles):
      for j in xrange(len(columns)):
         name = columns[j].name
         event_table.data.field(name)[previous_loc:previous_loc+num_events[i]] = handle[table_name].data.field(name)[:]
      previous_loc += num_events[i]
   return event_table

def __common_columns__(handles,table_name = 'EVENTS'):
   """Find the columns common to all files."""
   #Quick kluge - return the shortest one...
   all_cols = [handle[table_name].columns for handle in handles]
   arg = N.argmin([len(col) for col in all_cols])
   return all_cols[arg],handles[arg][table_name].header
   
   """   
   selector = {}
   for col in all_cols:
      for n in col: selector[n.name] = n
   for col in all_cols:
      names = [n.name for n in col]
      for key in selector.keys():
         if not (key in names): selector.pop(key)
   print selector
   """
         
   


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
   """

   if cuts is None: return
   
   from numarray import array,logical_and #some installations may need the numpy version of these calls

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
   