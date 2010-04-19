""" manage times 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/fermitime.py,v 1.1 2010/01/13 22:09:52 burnett Exp $

"""
import datetime

class MET(object):
    """ convert time in MET to a datetime object"""
    mission_start = datetime.datetime(2001,1,1)
    def __init__(self, met):
        if met>252460801: met=met-1 # 2008 leap second
        if met>157766400: met=met-1 # 2005 leap second
        self.time = MET.mission_start + datetime.timedelta(0,met)
    def __str__(self):
        return str(self.time)

def utc_to_met(year,month,day,hour = 0,min = 0,sec =0):
    """Convert a datetime in utc to MET."""
    utc = datetime.datetime(year,month,day,hour,min,sec)
    diff = utc-datetime.datetime(2001,1,1)
    leap_secs = 0
    if utc.year>2005: leap_secs+=1
    if utc.year>2008: leap_secs+=1
    return diff.days*86400+diff.seconds+leap_secs


def date_tag():
    """ useful to tag plots"""
    import  pylab
    pylab.figtext(0.04, 0.02, str(datetime.datetime.today())[:16], size=8)

