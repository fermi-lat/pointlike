import pyfits as pf
import numpy as N
from types import ListType
from math import cos,sin,pi
from scipy.interpolate.interpolate import interp2d

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

def merge_ft1files(files,outputfile = 'merged_ft1.fits',cuts = None, columns = None):
   """Pass either a name of an ASCII file listing the FT1 files,a list of FT1 filenames, or a single FT1 file.
   
      See __arbitrary_cuts__ for explanation of how to supply cuts."""
   
   handles = __get_handles__(files)

   event_table = __merge_events__(handles)      
   gti_table   = __merge_gti__(handles)

   if cuts is not None and columns is not none:
      __arbitrary_cuts__(event_table,cuts,columns)

   #Overwrite data in dummy table and write it to file
   handles[0]['EVENTS'].data = event_table.data
   handles[0]['GTI'].data = gti_table.data
   handles[0].writeto(outputfile,clobber=True)
   
   for x in handles: x.close()

class Livetime(object):

   def __init__(self,ft2file,**kwargs):
      self.init()
      self.__dict__.update(kwargs)
      if self.ft1files is not None: self.__get_gti__()

      ft2 = pf.open(ft2file,memmap=1)
      livetimes,b_mask = self.__process_gtis__(ft2)
      z_mask = livetimes > 0.

      self.livetimes = livetimes[z_mask]
      self.ra_scz = N.asarray(ft2['SC_DATA'].data.field('RA_SCZ'))[b_mask][z_mask]*pi/180.
      self.dec_scz = N.asarray(ft2['SC_DATA'].data.field('DEC_SCZ'))[b_mask][z_mask]*pi/180.
      self.ra_zenith = N.asarray(ft2['SC_DATA'].data.field('RA_ZENITH'))[b_mask][z_mask]*pi/180
      self.dec_zenith = N.asarray(ft2['SC_DATA'].data.field('DEC_ZENITH'))[b_mask][z_mask]*pi/180

      ft2.close()

   def init(self):
      self.ft1files = None
      self.gti = None
      self.zenithcut = 180.#105.
      self.fovcut = 180.#66.
      #self.bins = N.linspace(.4,1.,20)
      self.bins = N.linspace(.2,1.,20)

   def __get_gti__(self):
      handles = __get_handles__(self.ft1files)
      starts,stops = __merge_gti__(handles,no_table=True)
      self.gti = [[starts[i],stops[i]] for i in xrange(len(starts))]
      for handle in handles: handle.close()

   def __process_gtis__(self,ft2):

      livetimes = N.asarray(ft2['SC_DATA'].data.field('LIVETIME'))

      if self.gti is not None:

         ft2_starts = N.asarray(ft2['SC_DATA'].data.field('START'))
         ft2_stops = N.asarray(ft2['SC_DATA'].data.field('STOP'))

         min_gti,max_gti = self.gti[0][0],self.gti[-1][1]
         mask = N.logical_and(ft2_starts < max_gti, ft2_stops > min_gti)

         ft2_starts = ft2_starts[mask]
         ft2_stops = ft2_stops[mask]
         livetimes = livetimes[mask]            
         overlaps = N.empty_like(livetimes)

         print ft2_starts[0],self.gti[0][0]

         #A perhaps faster implementation would loop over FT2 and GTI simultaneously,
         #at the cost of O(a few seconds*GTI_intervals) of inaccuracy
         gti_counter = 0
         for i in xrange(len(ft2_starts)):
            start,stop = ft2_starts[i],ft2_stops[i]
            gstart,gstop = self.gti[gti_counter]
            
            if start >= gstart:
               if stop < gstop: continue
               else: lt_frac_so_far = max(0,gstop - start)
            else: continue #Must have skipped ahead a GTI, advance the FT2 file

            gti_counter +=1
            gstart,gstop = self.gti[gti_counter]
            if stop < gstart: pass
            else: lt_frac_so_far += max(0,stop-gstart)
            livetimes[i]*=(lt_frac_so_far)/(stop-start)              
               
         #This is a very clean & concise approach, also totally accurate
         #for i,gti_interval in enumerate(self.gti):
         #   maxi = N.maximum(gti_interval[0],ft2_starts)
         #   mini = N.minimum(gti_interval[1],ft2_stops)
         #   overlaps += N.maximum(0,mini - maxi)
         #livetimes = overlaps/(ft2_stops-ft2_starts)*livetimes
      
         return livetimes,mask

      return livetimes,N.asarray([True]*len(livetimes))

   def livetime(self,skydir):

      local_vars = ['zenithcut','fovcut','livetimes','ra_scz','dec_scz','ra_zenith','dec_zenith','bins']
      for var in local_vars: exec('%s=self.%s'%(var,var))         

      #Cut on FOV
      ra,dec = float(skydir.ra())*pi/180.,float(skydir.dec())*pi/180.
      cosines = N.cos(dec_scz)*cos(dec)*N.cos(ra-ra_scz) + N.sin(dec_scz)*sin(dec)
      mask = cosines >= cos(fovcut*pi/180.)

      #Cut on zenith
      zenith_cosines = N.cos(dec_zenith)*cos(dec)*N.cos(ra-ra_zenith) + N.sin(dec_zenith)*sin(dec)
      mask = N.logical_and(mask,zenith_cosines >= cos(zenithcut*pi/180.))

      total_cut_rate = float(len(mask[mask]))/len(mask)
      print 'Total cut rate: %f'%total_cut_rate

      #Get rid of any values outside of bins
      mask = N.logical_and(mask,N.logical_and(cosines>=bins[0],cosines<=bins[-1]))
      livetimes = livetimes[mask]
      cosines = cosines[mask]

      #Bin in cosine and put in the appropriate livetime
      sorting = N.argsort(cosines)
      livetimes = livetimes[sorting]
      cosines = cosines[sorting]
      summed_livetimes = N.zeros(len(bins)-1)

      counter = 0
      current_bin = bins[counter+1]
      low = 0
      for i,c in enumerate(cosines):
         if c <= current_bin: continue
         summed_livetimes[counter] = livetimes[low:i].sum()
         counter += 1
         current_bin = bins[counter]
         low = i
      summed_livetimes[-1] = livetimes[low:].sum()
      return (summed_livetimes,bins)


class EffectiveArea:

   def __init__(self):

      path = r'f:\glast\packages\ScienceTools-v9r7p1\irfs\caldb\v0r7p1\CALDB\data\glast\lat\bcf\ea'
      front_file = r'\aeff_P6_v1_diff_front.fits'
      back_file = r'\aeff_P6_v1_diff_back.fits'
      ea = pf.open(path+front_file)
      cbins = N.append(ea['EFFECTIVE AREA'].data.field('CTHETA_LO')[0],ea['EFFECTIVE AREA'].data.field('CTHETA_HI')[0][-1])
      ebins = N.append(ea['EFFECTIVE AREA'].data.field('ENERG_LO')[0],ea['EFFECTIVE AREA'].data.field('ENERG_HI')[0][-1])
      feffarea = N.array(ea['EFFECTIVE AREA'].data.field('EFFAREA')[0])
      ea.close()
      ea = pf.open(path+back_file)
      beffarea = N.array(ea['EFFECTIVE AREA'].data.field('EFFAREA')[0])
      ea.close()
      self.cbins,self.ebins = cbins,ebins
      nc,ne = len(self.cbins),len(self.ebins)
      self.feffarea,self.beffarea = feffarea.reshape(nc-1,ne-1),beffarea.reshape(nc-1,ne-1)
      self.i_ebins,self.i_cbins = N.log((ebins[:-1]*ebins[1:])**0.5),(cbins[1:]+cbins[:-1])/2.

   def image(self,event_class=-1,interpolate = False):

      if event_class < 0: effarea = self.feffarea + self.beffarea
      elif event_class == 0: effarea = self.feffarea
      else: effarea = self.beffarea
      ebins,cbins = self.ebins,self.cbins

      import pylab as P
      
      P.figure(2)
      P.gca().set_xscale('log')
      P.pcolor((ebins[:-1]*ebins[1:])**0.5,(cbins[:-1]+cbins[1:])/2.,effarea.reshape(len(cbins)-1,len(ebins)-1))
      P.figure(3)
      P.imshow(effarea,origin='lower',interpolation=None)

      energies = N.logspace(N.log10(ebins[0]),N.log10(ebins[-1]),120)
      #cthetas = N.array([0.975]*60)#N.linspace(cbins[0],cbins[-1],60)
      vals = N.array([self(e,.99,interpolate=True) for e in energies])
      P.figure(4)
      P.gca().set_xscale('log')
      P.plot(energies,vals)
      vals = N.array([self(e,.99,interpolate=False) for e in energies])
      P.plot(energies,vals)
      #P.pcolor(energies,cthetas,vals.reshape(len(cthetas),len(energies)))
      #P.figure(5)
      #P.imshow(vals.reshape(len(cthetas),len(energies)).transpose(),origin='lower',interpolation=None)

   def __call__(self,e,c,event_class=-1,interpolate = True):
      """Return bilinear (or nearest-neighbour) interpolation."""
      eb,cb = self.i_ebins,self.i_cbins
      e = N.log(e)
      ne,nc = len(eb),len(cb)

      if e < eb[0]: e = eb[0]
      if e > eb[-1]: e = eb[-1]
      if c < cb[0]: c = cb[0]
      if c > cb[-1]: c = cb[-1]

      #for i in xrange(ne-2):
      #      if eb[i] <= e and eb[i+1] >= e: break
      #for j in xrange(nc-2):
      #      if cb[j] <= c and cb[j+1] >= c: break
      
      #Cute way to find nearest neighbour
      i = N.argmin(N.abs(eb-e))
      j = N.argmin(N.abs(cb-c))

      # effarea[:,-1] increasing effective area
      # effarea[-1,:] increasing effective area
      # effarea[:,-1].shape 32
      # effarea[-1,:].shape 64

      if not interpolate:
         #if (eb[i+1]-e) < (e - eb[i]) and i < ne-1: i+=1
         #if (cb[j+1]-c) < (c - cb[j]) and j < nc-1: j+=1
         if event_class < 0: return (1e4*self.feffarea[i,j],1e4*self.beffarea[i,j])
         elif event_class == 0: return 1e4*self.feffarea[i,j]
         else: return 1e4*self.beffarea[i,j]
         return 1e4*effarea[j,i]

      i = i if eb[i]<=e and i < ne-1 else i-1 #adjust nearest neighbor to lower bin
      j = j if cb[j]<=c and j < nc-1 else j-1

      def bilinear(effarea):
      
         c2,c1 = cb[j+1],cb[j]
         e2,e1 = eb[i+1],eb[i]
         f00 = effarea[j,i]
         f11 = effarea[j+1,i+1]
         f01 = effarea[j+1,i]
         f10 = effarea[j,i+1]

         return 1e4/(e2-e1)/(c2-c1)*( (e2-e)*(f00*(c2-c) + f01*(c-c1)) + (e-e1)*(f10*(c2-c) + f11*(c-c1)) )

      if event_class < 0: return (bilinear(self.feffarea),bilinear(self.beffarea))
      elif event_class == 0: return bilinear(self.feffarea)
      else: return bilinear(self.beffarea)

   def dual_avg(self,e_range,c_range,event_class=-1,steps=10):
      return N.array([self(x,y,event_class=event_class) for x in N.linspace(e_range[0],e_range[1],steps) for y in N.linspace(c_range[0],c_range[1],steps)]).sum()/steps**2

   def theta_avg(self,e,c_range,event_class=-1,steps=10):
      return N.array([self(e,ctheta,event_class=event_class) for ctheta in N.linspace(c_range[0],c_range[1],steps)]).sum()/steps

class Exposure:

   def __init__(self,ft2file,ft1files=None):
      self.ft2file = ft2file
      self.ea = EffectiveArea()
      self.lt = Livetime(ft2file,ft1files=ft1files,bins=self.ea.cbins)      

   def __call__(self,skydir,energies,event_class=-1):

      e_centers = N.asarray(energies)
      lt = self.lt.livetime(skydir)
      c_centers = (lt[1][1:]+lt[1][:-1])/2.
      
      vals = N.array([[self.ea(e_center,c_center,event_class=event_class) for c_center in c_centers] for e_center in e_centers])
      if event_class == -1:
         return N.append(N.sum(vals[:,:,0]*lt[0],axis=1),N.sum(vals[:,:,1]*lt[0],axis=1))
      else: return N.sum(vals*lt[0],axis=1)

      
class BinnedPicture(object):
   
   def __init__(self,ft1files,bins=None):
      self.ft1files = ft1files
      self.bins = 10**N.arange(2,5.01,0.1) if bins is None else bins

      handles = __get_handles__(ft1files)
      energies = N.concatenate([x['EVENTS'].data.field('ENERGY') for x in handles])
      sorting = N.argsort(energies)
      energies = energies[sorting]
      mask = N.logical_and(energies>=self.bins[0],energies<=self.bins[-1])
      ras = N.concatenate([x['EVENTS'].data.field('L') for x in handles])
      decs = N.concatenate([x['EVENTS'].data.field('B') for x in handles])

      self.energies,self.ras,self.decs = energies[mask],ras[sorting][mask],decs[sorting][mask]
      for x in handles: x.close()

      import psf
      self.psf = psf.PSF()

   def __smoothed__(self,in_pairs,e):

      #Find bin
      bins = self.bins
      if e < bins[0] or e > bins[-1]: return
      for i in xrange(len(bins)-1):
         if bins[i] <= e and bins[i+1] >=e: break
      bin_cent = (bins[i]*bins[i+1])**0.5

      sigma = self.__smoothing_factor__(bin_cent)/180.*N.pi

      mask = N.logical_and(self.energies>=bins[i],self.energies<bins[i+1])
      ras,decs = self.ras[mask].astype(float)*N.pi/180.,self.decs[mask].astype(float)*N.pi/180.

      sigma / float(len(ras))**(1./6) #factor for KDE -- dunno if should be used or not!

      from math import cos,sin,pi
      in_ras,in_decs = (N.asarray(in_pairs)*N.pi/180.).transpose()
      diffs = N.empty_like(in_ras)
      for i in xrange(len(in_ras)):
         mask = N.abs(in_decs[i]- decs) < 3*sigma
         diffs[i] = N.exp( (N.cos(decs[mask])*cos(in_decs[i])*N.cos(ras[mask]-in_ras[i]) + N.sin(decs[mask])*sin(in_decs[i])-1)/sigma**2).sum()
      return diffs

   def __smoothing_factor__(self,e):
      return self.psf.sigma(e)
      
   def quick_grid(self):
      ras = list(N.linspace(180,540,720))
      ras.reverse()
      ras = N.asarray(ras)
      decs = N.linspace(-30,30,240)
      pairs = [[x,y] for x in ras for y in decs]
      return self.__smoothed__(pairs,1010)


   def psf_smooth(self):
      handles = __get_handles__(self.ft1files)

      #Assume just one file for now
      ft1 = handles[0]

      #First, determine energy bins


#EVERYTHING BELOW IS AN INTERNAL CALL.

def __get_handles__(files):
   if type(files) is not ListType:
      try:
         f = open(files)
         f.close()
         files = f.readlines()
         files = [x.strip() for x in files]
         files = [x for x in files if len(x) > 0]
      except:
         #files is a single filename
         files = [files]

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

def __merge_events__(handles):
   """Return a FITS table of merged FT1 events."""
   
   num_events = [x['EVENTS'].data.shape[0] for x in handles]
   columns = handles[0]['EVENTS'].columns
   event_table = pf.new_table(columns,nrows=sum(num_events))
   previous_loc = 0
   for i,handle in enumerate(handles):
      for j in xrange(len(columns)):
         event_table.data.field(j)[previous_loc:previous_loc+num_events[i]] = handle['EVENTS'].data.field(j)[:]
      previous_loc += num_events[i]
   return event_table

def __merge_gti__(handles,no_table=False):
   """Return a FITS table of merged GTI."""

   #Merge the gti and sort the results
   starts = N.concatenate([x['GTI'].data.field('START') for x in handles])
   stops  = N.concatenate([x['GTI'].data.field('STOP')  for x in handles])
   sorting = N.argsort(starts)
   starts = starts[sorting]
   stops = stops[sorting]

   if no_table: return (starts,stops)

   #Put GTI in new table
   gti_table = pf.new_table(handles[0]['GTI'].columns,nrows = len(starts))
   gti_table.data.field('START')[:] = starts
   gti_table.data.field('STOP')[:] = stops

   return gti_table

def __arbitrary_cuts__(events,cuts,columns):
   """Perform the cuts provided as arguments.

      events -- a major table in FITS format such as might be returned by __merge_events__
      cuts -- a list of cuts; the final result is the intersection of all cuts.  e.g,
              ['ENERGY>100','MET > 5000','MET < 10000']
      columns -- a list of the columns involved in the cuts to save me some parsing effort"""

   from numarray import array as narray
   mask = N.asarray([True]*events.data.shape[0])
   for entry in columns:
      exec('%s = events.data.field(\'%s\')'%(entry,entry))
   for cut in cuts:
      exec('mask = N.logical_and(mask,%s)'%cut)
   print int(sum(mask))
   new_table = pf.new_table(events.columns,nrows=int(sum(mask)))
   mask = narray(mask)
   for i in xrange(len(columns)):
      new_table.data.field(i)[:] = events.data.field(i)[mask]
   events.data = new_table.data
   