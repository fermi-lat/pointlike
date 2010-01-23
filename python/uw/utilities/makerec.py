import matplotlib
import numpy as np
import pyfits

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
    r = matplotlib.mlab.csv2rec(filename, skiprows=1, delimiter=delimiter, names=names)
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

        