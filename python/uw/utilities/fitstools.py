"""A suite of tools for processing FITS files.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/fitstools.py,v 1.21 2013/10/24 10:16:19 kerrm Exp $

   author: Matthew Kerr

"""

INT_TYPES  = ['EVENT_CLASS','CONVERSION_TYPE']


from astropy.io import fits as pyfits; pf= pyfits
import numpy as N; import numpy as np
from types import ListType,FunctionType,MethodType
from math import cos,sin,pi
from skymaps import SkyDir,Gti,BinnedPhotonData,PythonUtilities

def rect_mask(lons,lats,cut_lon,cut_lat,lon_hwidth,lat_hwidth):
    """ lons -- longitude coordinate (deg.)
        lats -- latitude coordinate (deg.)
    """
    mask = np.abs( lons - cut_lon)/np.cos(np.radians(lats)) < lon_hwidth
    return mask & (np.abs(lats - cut_lat) < lat_hwidth)

def trap_mask(ras,decs,cut_dir,radius):
    """ Make a conservative, trapezoidal cut as a precursor to a radius cut."""
    try:
        rad_test = radius[0]
    except:
        rad_test = radius
    if rad_test > 30:
        return np.asarray([True] * len(ras)) # cut no good for large radii

    c1 = np.cos( np.radians(cut_dir.dec() - radius) )
    c2 = np.cos( np.radians(cut_dir.dec() + radius) )
    ra_radius =  radius/np.minimum(c1,c2) #conservative

    mask = np.abs(ras - cut_dir.ra()) <= ra_radius
    mask &= np.abs(decs - cut_dir.dec()) <= radius
    #If cut radius overlaps the pole, keep a polar cap above cut_dir.dec()
    if np.any(cut_dir.dec()+radius > np.pi/2.) or np.any(cut_dir.dec() - radius < -N.pi/2.):
        mask |= np.abs(decs) > np.abs(cut_dir.dec())
    mask[radius > 20] = True # cut doesn't work for large radii?
    return mask

def rad_mask(ras,decs,cut_dir,radius,mask_only=False):
    """Make a slower, exact cut on radius."""
    ra0,dec0 = np.radians(cut_dir.ra()),np.radians(cut_dir.dec())
    ras,decs = np.radians(ras),np.radians(decs)
    cos_diffs = np.sin(decs)*np.sin(dec0)+np.cos(decs)*np.cos(dec0)*np.cos(ras-ra0)
    mask = cos_diffs > np.cos(np.radians(radius))
    if mask_only:
        return mask
    else:
        return mask,np.arccos(cos_diffs)[mask]

def get_gti_mask(ft1file,times):
    gti = Gti(ft1file)
    gti_starts,gti_stops = \
        np.asarray([(x.minValue(),x.maxValue()) for x in gti]).transpose()
    a = np.argsort(gti_stops)
    gti_starts = gti_starts[a]; gti_stops = gti_stops[a]
    indices = np.searchsorted(gti_stops,times)
    accept = (times > gti_starts[indices]) & (times <= gti_stops[indices])
    return accept

def rad_extract(eventfiles,center,radius_function,return_cols=['PULSE_PHASE'],cuts=None,apply_GTI=True,theta_cut=66.4,zenith_cut=105,return_indices=False):
    """ Extract events with a radial cut.  
        Return specified columns and perform additional boolean cuts.

        Return is in form of a dictionary whose keys are column names 
        (and 'DIFFERENCES') and values are numpy arrays with the column 
        values.  These will have been concatenated if there are multiple
        FT1 files.

    =========   =======================================================
    Argument    Description
    =========   =======================================================
    eventfiles  -- a list of FT1 filenames
    center      -- a SkyDir giving the center of the radial cut
    radius_function -- can be either a float specifying a cookier cutter 
                radial cut, or a function taking as arguments the energy 
                and event_class and speciying the radius in degrees, e.g.

              def radius(energy,event_class):
                 return numpy.where(event_class,2,1)*(energy/1000)**-0.75

    =========   =======================================================
    Keyword     Description
    =========   =======================================================
    return_cols ['RA','DEC','ENERGY','EVENT_CLASS','PULSE_PHASE'] - 
                a list of FT1 column names to return
    cuts        None - an optional list of boolean cuts to apply, 
                e.g., ['ENERGY > 100']
                NB -- cuts not yet implemented!!
    no_cuts     [False] do not apply default zenith and incidence angle cuts
    apply_GTI   [True] accept or reject an event based on GTI if True; 
                else ignore GTI
    return_indices [False] if True, return an array giving the index in the
                original file of each event; obviously only useful in the 
                case of a single event file
    =========   =======================================================
    """
    if not hasattr(radius_function,'__call__'):
        simple_scalar = True
        rval = radius_function
        radius_function = lambda e,event_class: rval
    else:
        simple_scalar = False

    eventfiles = __FITS_parse__(eventfiles)

    from collections import defaultdict,deque
    coldict = defaultdict(deque)
    cols = {}
    cut_cols = ['ZENITH_ANGLE','THETA','TIME']
    keys = list(set(['RA','DEC','ENERGY','CONVERSION_TYPE']+cut_cols+return_cols))
    accepted = 0
    total = 0

    for eventfile in eventfiles:
        #e = pf.open(eventfile,memmap=1)
        #nrows = e[1].data.shape[0]
        #e.close()
        nrows = pyfits.getheader(eventfile,'EVENTS')['NAXIS2']

        for key in keys:
            cols[key] = np.empty(nrows,dtype=float)
            PythonUtilities.get_float_col(cols[key],eventfile,'EVENTS',key)

        rad   = radius_function(cols['ENERGY'],cols['CONVERSION_TYPE'])
        tmask = trap_mask(cols['RA'],cols['DEC'],center,rad)
        tmask &= (cols['ZENITH_ANGLE'] < zenith_cut) & (cols['THETA'] < theta_cut)
        if apply_GTI:
            tmask &= get_gti_mask(eventfile,cols['TIME'])
            print 'GTI will remove %d of %d photons.'%((~tmask).sum(),len(tmask))
        if simple_scalar:
            rmask,diffs = rad_mask(cols['RA'][tmask],cols['DEC'][tmask],center,rad)
        else:
            rmask,diffs = rad_mask(cols['RA'][tmask],cols['DEC'][tmask],center,rad[tmask])

        for key in keys:
            coldict[key].append(cols[key][tmask][rmask])
        if return_indices:
            if 'EVENT_INDICES' not in return_cols:
                return_cols.append('EVENT_INDICES')
            coldict['EVENT_INDICES'].append(np.arange(len(tmask))[tmask][rmask])
        coldict['DIFFERENCES'].append(diffs)
        accepted += tmask.sum()
        total += len(tmask)

    for key in coldict.keys():
        if (key in cut_cols) and not (key in return_cols):
            cols.pop(key)
            continue
        cols[key] = np.concatenate([x for x in coldict[key]])
        if key in INT_TYPES: cols[key] = cols[key].astype(int)

    print 'Cuts removed %d of %d photons.'%(total-accepted,total)
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

def get_fields(files, fields, cuts = None, memmap = False):
   """A lightweight version to get only certain fields of flight data."""

   data = dict()
   files = __FITS_parse__(files)

   #Process the cuts
   f = pf.open(files[0],memmap=memmap)
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
   del(f)

   counts = N.empty(len(files),dtype=int)

   for nfi,fi in enumerate(files):
      f = pf.open(fi,memmap=memmap)
      counts[nfi] = len(f['EVENTS'].data) #f['EVENTS'].data.getshape()[0]
      f.close()
      del(f)

   for field in fields + cut_fields:
      data[field] = N.empty(counts.sum(),dtype=float)

   #Get fields
   counter = 0
   for fi,counts in zip(files,counts):
      f = pf.open(fi,memmap=memmap)
      for field in fields + cut_fields:
         #data[field] += [N.array(f['EVENTS'].data.field(field),dtype=N.float32)]
         data[field][counter:counter+counts] = f['EVENTS'].data.field(field)
      counter += counts
      f.close()
      del(f)

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
      if cut not in fields: data.pop(cut)

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

    files = __FITS_parse__(files)

    try:
        f = pf.open(files[0])
        header = f['EXPOSURE'].header
        exposures = [f['EXPOSURE'].data.field('COSBINS')] 
    except KeyError:
        print('file %s has no EXPOSURE table: aborting'%files[0])
        return
    finally:
        f.close() 

    for file in files[1:]:
        try:
            f = pf.open(file)
            h = f['EXPOSURE'].header
            for key in header.keys():
                if key not in ['DATE','DATE-OBS','DATE-END','TSTART','TSTOP']:
                   assert(h[key]==header[key])
            exposures+=[f['EXPOSURE'].data.field('COSBINS')]
        except AssertionError:
            print('Inconsistent header values, file %s, keyword %s'%(file,key))
            return
        except KeyError:
            print('File %s has no table "EXPOSURE": aborting'%file)
            return
        finally:
            f.close()

    summed_exposure = N.array(exposures).sum(axis=0)
    null = header.get('TNULL1', None) #THB: kluge fix.
    coldef = pf.ColDefs([pf.Column(name='COSBINS',format=header['TFORM1'],null=null,array=summed_exposure)])
    exp_table = pf.new_table(coldef,header=header,nrows=exposures[0].shape[0])
    gti = merge_gti(files)
    exp_table.writeto(outputfile,clobber=True)
    gti.writeExtension(outputfile)

def merge_gti(files,table_name = 'GTI',interval = None):
    """Return the union of Gtis specified in files, with an optional time range cut."""
    gtis = [Gti(f,table_name) for f in files]
    new_gti = Gti(gtis[0])
    for gti in gtis:
        new_gti.combine(gti)

    if interval is not None:
        new_gti = new_gti.applyTimeRangeCut(*interval)

    return new_gti

def merge_bpd(bpd_files,outfile = None):
    """Merge a set of BinnedPhotonData files.
    
    outfile: File to write the merged BinnedPhotonData to.  If None, don't save it."""
    bpds = [BinnedPhotonData(bf) for bf in bpd_files]
    bpds = [bpd for bpd in bpds if bpd.gti().computeOntime()>0] #Ignore entries with empty GTIs
    new_bpd = bpds[0]
    for b in bpds[1:]:
        new_bpd.add(b)
    if outfile:
        new_bpd.write(outfile)

    return new_bpd

def merge_lt(lt_files,outfile = 'merged_lt.fits',weighted = True):
    """Merge a list of LivetimeCube files, handling the WEIGHTED_EXPOSURE extension.

    kwargs:
        outfile: File to save the merged LivetimeCube to
        weighted: If true, merge WEIGHTED_EXPOSURE tables, in addition to
                  EXPOSURE.
    """

    files = __FITS_parse__(lt_files)

    try:
        f = pf.open(files[0])
        p_header = f['PRIMARY'].header
        err = 'EXPOSURE'
        header = f['EXPOSURE'].header
        exposures = [f['EXPOSURE'].data.COSBINS]
        err = 'WEIGHTED_EXPOSURE'
        if weighted:
            w_header = f['WEIGHTED_EXPOSURE'].header
            w_exposures = [f['WEIGHTED_EXPOSURE'].data.COSBINS]
    except KeyError:
        print('file %s has no %s table: aborting'%(files[0],ext))
        return
    finally:
        f.close() 

    for file in files[1:]:
        try:
            f = pf.open(file)
            ext = 'EXPOSURE'
            h = f['EXPOSURE'].header
            for key in header.keys():
                assert(h[key]==header[key])
            exposures+=[f['EXPOSURE'].data.COSBINS]

            if weighted:
                ext = 'WEIGHTED_EXPOSURE'
                hw = f['WEIGHTED_EXPOSURE'].header
                for key in w_header.keys():
                    assert(hw[key]==w_header[key])
                w_exposures += [f['WEIGHTED_EXPOSURE'].data.COSBINS]

        except AssertionError:
            print('Inconsistent header values, file %s, extension %s, keyword %s'%(file,ext,key))
            return
        except KeyError:
            print('File %s has no table "%s": aborting'%(file,ext))
            return
        finally:
            f.close()
    primary = pf.PrimaryHDU(data=None,header = p_header)
    hdulist = pf.HDUList([primary])
    summed_exposure = N.array(exposures).sum(axis=0)
    null = header.get('TNULL1',None)
    coldef = pf.ColDefs([pf.Column(name='COSBINS',format=header['TFORM1'],null=null,array=summed_exposure)])
    exp_table = pf.new_table(coldef,header=header,nrows=exposures[0].shape[0])
    hdulist.append(exp_table)
    if weighted:
        summed_w_exposure = N.array(w_exposures).sum(axis=0)
        null = w_header.get('TNULL1',None)
        coldef_w = pf.ColDefs([pf.Column(name='COSBINS',format=w_header['TFORM1'],null=null,array=summed_w_exposure)])
        w_exp_table = pf.new_table(coldef_w,header=w_header,nrows=w_exposures[0].shape[0])
        hdulist.append(w_exp_table)
    hdulist.writeto(outfile,clobber=True)
    gti = merge_gti(files)
    gti.writeExtension(outfile)

#EVERYTHING BELOW IS AN INTERNAL CALL.

def __FITS_parse__(files):
   """Parse input and return a list of (FITS) filenames.

      files -- a glob-style wildcard, an ASCII file containing a list of filenames, a list of filenames, or a single FITS file."""
   #This doesn't detect ndarrays.
   #if type(files) is ListType: return files

   #Try testing for non-string iterability instead.
   if getattr(files,'__iter__',False):
      return files

   try: #See if it's a FITS file
      f = pf.open(files)
      f[0]
      f.close()
      return [files]
   except:
      pass
   if files[0] == '@':
      return [line.strip() for line in file(files[1:]) if len(line)>0 and line[0]!='#']
   from glob import glob
   return glob(files)

def __get_handles__(files):
   files = __FITS_parse__(files)
   handles = [pf.open(x,memmap=1) for x in files] #some versions may need to disable memmap (=0)
   return handles

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

