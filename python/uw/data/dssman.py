"""
Module to handle reading/writing/consistency of DSS keywords
in FITS files.

author(s): M. Kerr, T. Burnett
"""

__version__ = '$Revision: 1.8 $'
#$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/data/dssman.py,v 1.8 2016/06/22 17:02:50 wallacee Exp $

from astropy.io import fits as pyfits
from collections import deque
import numpy as np

#TODO class(es) for other DSS types, e.g. Region type
DSSKeys = ['TYP','UNI','VAL','REF']


def isfloat(string):
    try:
        float(string); return True
    except (ValueError,TypeError):
        return False

def isint(string):
    try:
        a = int(string)
        b = float(string)
        return a==b
    except (ValueError,TypeError):
        return False

def isnum(string):
    return isfloat(string)

def num(string):
    if isint(string): return int(string)
    if isfloat(string): return float(string)
    #raise ValueError,'Could not cast %s to float or int'%(string)
    return None

class DSSEntry(dict):
    def __init__(self,d):
        self.update(d)
        if 'index' not in self.keys(): self['index'] = -1
    def __copy__(self):
        return DSSEntry(d) # do we need copy?
    def __eq__(self,other):
        for key in DSSKeys:
            if self[key] != other[key]:
                return False
        return True
    def __ne__(self,other):
        return not self.__eq__(other)
    def __str__(self):
        return '\n'.join(['DS{0}{1}: {2}'.format(key,self['index'],self[key]) for key in DSSKeys])
    def to_header_entries(self,index=None):
        """ Return tuples of header key/val, e.g. to write to an FT1 file."""
        index = index or self['index']
        tup = []
        for key in DSSKeys:
            if self[key] is None: continue
            tup.append(['DS{0}{1}'.format(key,index),'{0}'.format(self[key])])
        return tup
    def roi_info(self):
        """ Check if specifies an extraction region."""
        if (('POS' in self['TYP']) and
            (self['UNI'] == 'deg') and
            ('CIRCLE' in self['VAL'])):
            try:
                tup = self['VAL'].split('(')[1][:-1].split(',')
                val = [float(x) for x in tup]
                assert(len(val)==3)
                return val
            except Exception:
                pass
        return None

class DSSSimpleRange(DSSEntry):
    """ A DSS entry specifying a single range on a column.
        E.g., ZENITH_ANGLE < 100 deg., 1e2 < ENERGY < 1e5."""
    def get_bounds(self):
        return self['lower'],self['upper']

    @staticmethod
    def bounds_to_val(lower,upper):
        if isnum(lower):
            if isnum(upper): val = '{0}:{1}'.format(lower,upper)
            else: val = '{0}:'.format(lower)
        elif isnum(upper): val = ':{0}'.format(upper)
        else: val = None
        return val

    def set_bounds(self,lower,upper):
        val = DSSSimpleRange.bounds_to_val(lower,upper)
        if val is None:
            raise ValueError('Both lower and upper were empty!')
        self['lower'] = num(lower); self['upper'] = num(upper)

    def intersection(self,other):
        """ Compare with another DSSSimpleRange and update bounds if other
            is more restrictive."""
        print(other.get_bounds())
        print(self.get_bounds())
        if (other['lower'] is not None):
            if self['lower'] is None: self['lower'] = other['upper']
            else: self['lower'] = max(self['lower'],other['lower'])
        if (other['upper'] is not None):
            if self['upper'] is None: self['upper'] = other['upper']
            else: self['upper'] = min(self['upper'],other['upper'])
            print(min(self['upper'],other['upper']))
            print(self['upper'])

class DSSBitMask(DSSSimpleRange):
    """ The only current use for a bit mask is selecting on EVENT_CLASS
        in Pass7+ data.  As such, only need is to keep DSS keywords
        consistent and pass the single masked bit into Data.
        
        The single bit is stored in both the lower and upper field."""
        # NOTES on EVENT_CLASS DSS format
        # PASS 7 style
        # DSTYP2 = 'BIT_MASK(EVENT_CLASS,2)'
        # DSUNI2 = 'DIMENSIONLESS'
        # DSVAL2 = '1:1 '
        # Pass 6 style
        # DSTYP3 = 'EVENT_CLASS'
        # DSUNI3 = 'dimensionless'
        # DSVAL3 = '3:10 ' 

    def __init__(self,d):
        super(DSSBitMask,self).__init__(d)
        t = self['TYP'][1:-1] #split off parens
        self['lower'] = self['upper'] = int(t.split(',')[1])

    def intersection(self,other):
        """ A bit mask intersection at the level we require is trivial --
            the masked bit MUST be the same."""
        if (other['lower'] is None) or (other['lower'] is None):
            raise ValueError('Not DSSBitMask format')
        if (self['lower'] is None) or (self['lower'] is None):
            raise ValueError('Not DSSBitMask format')
        if self['lower']!=other['lower'] or self['upper']!=other['upper']:
            raise ValueError('Bit masks are not compatible!')
        return

def DSSFactory(keys,vals):
    """ Return the correct type of DSS subclass for header entries."""

    # first, parse the values present in the set of header entries
    d = dict(TYP=None,UNI=None,VAL=None,REF=None)
    for k,v in zip(keys,vals):
        d['index'] = int(k[-1])
        for l in d.keys():
            if l in k : d[l] = v
        assert d[l] is not None, 'DSS key {}, value {} not valid'.format(k,v)

    # determine if cut is a bitmask
    if 'BIT_MASK' in d['TYP']:
        return DSSBitMask(d)

    # determine if cut is of a simple type
    assert d['VAL'] is not None, 'DSS key {}, bad,  value {}'.format(keys,vals)
    toks = d['VAL'].split(':') 
    if len(toks)==2:
        a,b = toks
        # conditions to accept: empty low and numeric high,
        # complement, and both are numbers
        if ((len(a)==0 and isnum(b)) or
            (len(b)==0 and isnum(a)) or
            (isnum(a) and isnum(b))):
            d['lower'] = num(a); d['upper'] = num(b)
            return DSSSimpleRange(d)

    return DSSEntry(d)

class DSSEntries(list):
    # NB -- may need to handle duplicate entries?

    def __init__(self,fits_name,header_key='EVENTS'):
        try:
            h = pyfits.getheader(fits_name,header_key)
        except IOError:
            print('Could not find file {0}'.format(fits_name))
            return
        except IndexError, msg: 
            print('Invalid header index for fits file %s: %s'% (fits_name,msg))
            return
        except KeyError, msg:
            print('Invalid key for fits file %s: %s'% (fits_name,msg))
            return
        keys = [x for x in h.keys() if x.startswith('DS')]
        if len(keys)==0: return
        vals = [h[k] for k in keys]
        # list of the indices, not necewssarily in numeric order (THB change 02/23/16)
        indeces = sorted(list(set([int(k[-1]) for k in keys])))
        kdeque,vdeque = deque(),deque()
        counter = 0 # index of list of DSS indeces
        for i in xrange(len(keys)):
            if int(keys[i][-1])!=indeces[counter]:
                self.append(DSSFactory(kdeque,vdeque)) 
                kdeque.clear(); vdeque.clear()
                counter += 1
            kdeque.append(keys[i]); vdeque.append(vals[i])
        self.append(DSSFactory(kdeque,vdeque))
        self._remove_duplicates()

    def __str__(self):
        return '\n'.join([str(x) for x in self])

    def __eq__(self,other):
        """Test for equality of DSS entries, independent of order."""
        return sorted(self)==sorted(other)
    def __ne__(self,other):
        return not self.__eq__(other)

    def get_simple_dss(self,colname):
        """ Return a DSS entry corresponding to a simple cut on colname."""
        for idss,dss in enumerate(self):
            if dss['TYP'] == colname:
                return dss,idss
            #Slight kludge to handle weird formatting of EVENT_CLASS bitmasks
            if colname=='EVENT_CLASS' and (colname in dss['TYP']):
                return dss,idss
        return None,None


    def delete(self,index):  
        """ Delete a DSS entry and re-index the remaining ones."""
        ret = self.pop(index)
        if index < len(self)-1:
            for i in xrange(index,len(self)):
                self[i]['index'] = i+1
        return ret

    def write(self,fits_name,header_key='EVENTS'):
        f = pyfits.open(fits_name,uint=False)
        h = f[header_key]._header
        for d in self:
            tup = d.to_header_entries()
            for t in tup:
                h[t[0]] = t[1]
        # convert unsigned ints to ints -- this is a kluge but perhaps necessary
        for hdu in f:
            if not isinstance(hdu,pyfits.BinTableHDU): continue
            for icol,col in enumerate(hdu.columns):
                if col.format=='1J':
                    #print 'update %s'%col.name
                    data = hdu.data.field(col.name).astype(np.int32) # apply transform
                    # not sure why above line must be done -- delayed, perhaps?
                    hdu.columns.change_attrib(col.name,'bzero',0)
                    # above line is crucial
        f.writeto(fits_name,clobber=True)
        #f.writeto('/tmp/test.fits',clobber=True)

    def roi_info(self,tol=0.01,delete_duplicates=False):
        """ Return any information about the extraction region.
            Will check for duplicates, which can arise from a slight 
            mismatch of astroserver and gtselect coordinates, and
            optionally delete them."""
        roi_info = None
        offset = 0
        for i in xrange(len(self)):
            i += offset
            d = self[i]
            r = d.roi_info()
            if (r is not None):
                if (roi_info is None):
                    roi_info = r
                else:
                    # need to check for compatibility
                    agree = all((abs(x-y)<tol for x,y in zip(r,roi_info)))
                    if not agree:
                        raise ValueError('Found multiple DSS keywords for ROI differing by greater than {0} degrees'.format(tol))
                    elif delete_duplicates:
                        self.delete(i)
                        offset += 1
        return roi_info
    
    def _remove_duplicates(self):
        duplicates = []
        for i, keyword in enumerate(self):
            if keyword in self[i+1:]: duplicates.append(i)

        offset = 0
        for dup in duplicates:
            self.delete(dup - offset)
            offset += 1

def make_simple_dss(colname,unit,low,high,index=1):
    """ Return a DSSSimpleRange object with bounds specified by low/high.
        index gives the order of the DSS keyword, e.g. DSVAL1, DSVAL2, etc. """
    val = DSSSimpleRange.bounds_to_val(low,high)
    if val is None:
        raise ValueError('Could not interpret arguments as a simple cut.')
    d = dict(REF=None,TYP=colname,UNI=unit,VAL=val,
        lower=num(low),upper=num(high),index=index)
    return DSSSimpleRange(d)

"""
def process_pixeldata(pd):
    # for now, assume all FT1 files have same DSS keywords...
    dsse = get_dss_entries(pd.ft1files[0])

    colnames = ['ZENITH_ANGLE','THETA','ENERGY','ENERGY']
    units = ['deg','deg','MeV','MeV']
    ptlvars = ['zenithcut','thetacut','emin','emax']
    indices = [1,1,0,1]

    for i in xrange(len(ptlvars)):
        # check the appropriate simple cut
        ptl_var = pd.__dict__[ptlvars[i]]
        dss = dsse.get_simple_entry(colnames[i])
        if dss is None:
            if not pd.quiet:
                print 'DSS keywords specified no %s cut.  Applying the specified pointlike cut %s = %.1f '%(colnames[i],ptlvars[i],ptl_var)
            if indices[i]:
                # NB only works for cut variables with range [0,max]...
                # TODO -- make work for all cuts (necessary?)
                dsse.append(make_simple_dss(colnames[i],units[i],0,ptl_var,index=len(dsse)+1))
            continue
        dss_var = dss.get_simple_bound(indices[i])
        dss_more_stringent = \
            (indices[i]==0 and dss_var>ptl_var) or \
            (indices[i]==1 and dss_var<ptl_var)
        ptl_more_stringent = (not dss_more_stringent) and (ptl_var != dss_var)
        sign_string = 'upper bound' if indices[i] else 'lower bound'
        if ptl_more_stringent:
            dss.set_simple_bound(ptl_var,indices[i])
            if not pd.quiet:
                print 'Applying more stringent pointlike %s %s cut (%s=%.1f) over that found in DSS keywords (%.1f)'%(colnames[i],sign_string,ptlvars[i],ptl_var,dss_var)
        elif dss_more_stringent:
            pd.__dict__[ptlvars[i]] = dss_var
            if not pd.quiet:
                print 'Applying more stringent DSS %s %s cut (%.1f) over that found in pointlike (%s=%.1f)'%(colnames[i],sign_string,dss_var,ptlvars[i],ptl_var)
        else:
            if not pd.quiet:
                print 'Identical cuts on %s %s (%s=%.1f)'%(colnames[i],sign_string,ptlvars[i],ptl_var)
    pd.dsse = dsse
"""

class DummyPD(object):
    """ For testing."""

    def __init__(self):
        self.ft1files = ['j1838m0536_p7_07062011-ft1_gtselect.fits']
        self.quiet = False
        self.zenithcut = 100.
        self.thetacut = 66.4
        self.emin = 1e2
        self.emax = 1e6

