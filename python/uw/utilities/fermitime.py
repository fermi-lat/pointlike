""" manage times 
$Header$

"""
import datetime

class MET(object):
    """ convert time in MET to a datetime object"""
    mission_start = datetime.datetime(2001,1,1)
    def __init__(self, met):
        self.time = MET.mission_start + datetime.timedelta(0,met)
    def __str__(self):
        return str(self.time)

def date_tag():
    """ useful to tag plots"""
    import  pylab
    pylab.figtext(0.04, 0.02, str(datetime.datetime.today())[:16], size=8)

