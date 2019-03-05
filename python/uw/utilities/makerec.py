""" Various useful utilities for creating, dumping numpy recarry objects
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/makerec.py,v 1.9 2013/09/14 05:39:01 burnett Exp $



"""
import os, pickle
from astropy.io import fits as pyfits
from pylab import mlab
import numpy as np

def makefits(r, filename=None, **kwargs):
    """ convert a recarray to a pyfits object, write to filename if present.
        **kwargs can be used to pass clobber to the pyfits writeto function.
    """
    def convertformat(dtype):
        if dtype[:2]=='|S': return dtype[2:]+'A'
        try:
            return {'<f8':'D', '<f4':'E', '<i4':'J', '|b1':'L', '|O4':'|O4','<i8':'I',}[dtype]
        except KeyError:
            print 'recarry type %s not recognized' %dtype
            raise
    def column(n,f):
        if f!='|O4': return pyfits.Column(name=n, format=f, array=r[n])
        # treat an object as a string, expand to match 
        maxsize = np.array([len(str(e)) for e in r[n]]).max()
        fmt = '%%-%ds' % maxsize
        data = [fmt % str(e) for e in r[n]] 
        return pyfits.Column(name=n, format='%dA'%maxsize, array=data)
    names = r.dtype.names
    formats = [convertformat(t[1]) for t in r.dtype.descr]
    columns = [column(n,f ) for n,f in zip(names, formats)]
    thdulist = pyfits.HDUList([ pyfits.PrimaryHDU(), 
                                pyfits.new_table(pyfits.ColDefs(columns))])
    if filename is not None:
        thdulist.writeto(filename, **kwargs)
    return thdulist
    

def fitsrec(filename, HDU=1, quiet=False):
    " make and return a record array from a FITS file"

    d=pyfits.open(filename)[HDU].data
    if not quiet: print 'loaded file, %s, found %d entries' % (filename, len(d))
    names = [name.lower() for name in d.names]

    fields = [np.asarray(d.field(name)).astype(float if format.find('A')==-1 else str )
                 for name,format in zip(names,d.formats)]
    trunc_fields = [f if len(f.shape)==1 else f[:,0] for f in fields]
    r = np.rec.fromarrays(trunc_fields, names=names)
    if not quiet: print r.dtype
    return r

def textrec(filename, quiet=False, insertname=False, delimiter=' '):
    nameline = open(filename).readline()
    if nameline[0]=='#':nameline=nameline[1:]
    delim=None if nameline.find(',')<0 else ','
    names = [name.lower().strip() for name in nameline.split(delim)]
    if insertname: names.insert(0,'name')
    r = mlab.csv2rec(filename, skiprows=1, delimiter=delimiter, names=names)
    if not quiet: 
        print 'loaded file %s, found %d entries' %(filename, len(r))
        print r.dtype

    return r


def multifitsrec(filenames, HDU=1, quiet=False):
    tables = [pyfits.open(filename)[HDU].data for filename in filenames]
    recs = sum([t.shape[0] for t in tables]) #total length
    names = [name.lower() for name in tables[0].names]
    formats = tables[0].formats
    fields = [np.zeros(recs, float if format.find('A')==-1 else str) for  format in formats]
    i=0
    for t in tables:
        n = t.shape[0]
        for f, name in zip(fields, names):
            f[i: i+n] = t.field(name)
        i += n
    r = np.rec.fromarrays(fields, names=names)

class RecArray(object):
    def __init__(self, names, fmt=None):
        self.names = names
        self.n = len(names)
        self.fields=list([list() for i in range(self.n)])

    def append(self,*arg):
        for i,x in enumerate(arg):
            self.fields[i].append(x)
    def __call__(self):
        """ return finished recarray"""
        return np.rec.fromarrays(self.fields, names=self.names)

      
def load(filename):
    ext = os.path.splitext(filename)[1]
    if ext=='.txt':
        return textrec(filename)
    elif ext=='.pickle' or ext=='.rec':
        return pickle.load(open(filename))
    elif ext=='.csv':
        return textrec(filename, delimiter=',')
    raise Exception('extension %s not recognized: expect txt, rec or pickle' %ext)
