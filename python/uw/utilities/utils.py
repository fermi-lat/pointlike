"""A collection of convenience functions, included mostly for the sake of get_data, and its dependencies.
Those should probably be moved to their own module at some point, and the sexigesimal-decimal conversions
subsumed into a more sensible module, but for now this works.

Author: E. Wallace
"""
import os
from datetime import datetime,timedelta
import calendar
from glob import glob
import numpy as np
from skymaps import Gti
from uw.utilities.fitstools import merge_gti
from uw.utilities.fermitime import MET,utc_to_met


def hms_to_deg(hr,min,sec):
    """Convert a (hr,min,sec) to decimal degrees"""
    return hr*15+min/4.+sec/240.

def dms_to_deg(deg,min,sec):
    """Convert a (deg,arcmin,arcsec) to decimal degrees"""
    return deg+min/60.+sec/3600.

def deg_to_dms(deg):
    """Convert decimal degrees to (deg,arcmin,arcsec)"""
    d = int(deg)
    deg-=d
    m = int(deg*60.)
    s=deg-m//60 
    return d,m,s 

def deg_to_hms(deg):
    """Convert decimal degrees to (hr,min,sec)"""
    h = int(deg)//15
    deg -= h*15
    m = int(deg*4)
    s = (deg-m//4)/15.
    return h,m,s
    
def get_data(tstart,tstop,binsperdec = 4,data_dir = '/phys/groups/tev/scratch1/users/Fermi/data'):
    """Find saved data products for the time range tstart to tstop (MET).
    
    Directory structure and file names are assumed to be as described in the docstring
    for uw.like.pointspec.SavedData.
    """
    start_date = MyDate(*MET(tstart).time.timetuple()[:6])
    stop_date = MyDate(*MET(tstop).time.timetuple()[:6])
    files = dict(monthly=dict(bpd=None,lt=None),weekly=dict(bpd=None,lt=None),daily=dict(bpd=None,lt=None))
    for t in ['monthly','weekly','daily']:
            files[t]['bpd'] = np.asarray(sorted(glob(os.path.join(data_dir,t,'bpd','%s_*_%ibpd.fits'%(t[:-2].replace('i','y'),binsperdec)))))
            files[t]['lt'] = np.asarray(sorted(glob(os.path.join(data_dir,t,'lt','%s_*_lt.fits'%(t[:-2].replace('i','y'))))))
    month_mask,gti = accept_files(files['monthly']['bpd'],start_date,stop_date,months = True)
    week_mask,gti = accept_files(files['weekly']['bpd'],start_date,stop_date,gti=gti)
    day_mask,gti = accept_files(files['daily']['bpd'],start_date,stop_date,gti=gti)
    bpds = np.append(files['monthly']['bpd'][month_mask],np.append(files['weekly']['bpd'][week_mask],files['daily']['bpd'][day_mask]))
    lts = np.append(files['monthly']['lt'][month_mask],np.append(files['weekly']['lt'][week_mask],files['daily']['lt'][day_mask]))
    return bpds,lts

class MyDate(datetime):
    """Extend datetime with convenience functions for adding days,weeks, and months."""
    def add_day(self):
        new = self+timedelta(1,0,0)
        return MyDate(*(new.timetuple()[:6]))
    def add_week(self):
        new = self+timedelta(7,0,0)
        return MyDate(*(new.timetuple()[:6]))
    def add_month(self):
        new = self+timedelta(calendar.monthrange(self.year,self.month)[1],0,0)
        return MyDate(*(new.timetuple()[:6]))


def accept_files(files,start_date,stop_date,gti=None,months = False):
    """Decide what subset of files should be accepted based on start_date,stop_date and gti.

    Intended for use by get_data; there should be no reason to use it on its own, and it may
    have counterintuitive results if you do.  start_date and stop_date are datetime objects. 
    gti is a Gti describing time ranges that should NOT be accepted.
    """
    if not hasattr(files,'__iter__'):
        files = [files]
    if gti is None:
        gti = Gti()
    files = np.asarray(files)
    mask = np.array([True]*len(files))
    for i,f in enumerate(files):
        toks = os.path.basename(f).split('_')
        l = toks[0]
        mult = 100 if months else 1
        d = MyDate(*yyyymmdd_to_date(int(toks[1])*mult + months).timetuple()[:6])
        if d<start_date or (gti.minValue()<(utc_to_met(*(d+timedelta(0,1,0)).timetuple()[:6])) and
		            gti.maxValue()>(utc_to_met(*(d+timedelta(0,1,0)).timetuple()[:6]))):
            mask[i] = False
        d2 = eval('d.add_%s()'%l)
        if d2>stop_date or (gti.minValue()<(utc_to_met(*(d2-timedelta(0,1,0)).timetuple()[:6])) and
		            gti.maxValue()>(utc_to_met(*(d2-timedelta(0,1,0)).timetuple()[:6]))):
            mask[i] = False
    try:
        gti.combine(merge_gti(files[mask]))
    except IndexError:
        pass
    return mask,gti


def yyyymmdd_to_date(i):
    """Take a date in the form of an int with format yyyymmdd, and return a datetime"""
    y = i//10000
    m = (i-y*10000)//100
    d = i-y*10000-m*100
    return datetime(y,m,d)
