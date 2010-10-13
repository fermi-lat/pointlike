"""
do pivot stuff 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/pivot.py,v 1.7 2010/08/06 18:15:11 burnett Exp $

"""
version='$Revision: 1.7 $'.split()[1]
from uw.utilities import collection
from uw.thb_roi import pipeline, findsource
from skymaps import SkyDir
import numpy as np
import os, sys, exceptions

class InvalidParameter(exceptions.Exception):
    pass

def near(s, cat, tol=1.0):
    sdir = SkyDir(s.ra, s.dec)
    rtol = np.radians(tol)
    return [c for c in cat if sdir.difference(SkyDir(c.ra,c.dec))<rtol and c != s]

def selection(indir, select_all=False, mindist=0.25, other_keys=None):
    r0 = pipeline.load_rec_from_pickles(indir, other_keys=other_keys)
     #--- 
    if select_all: return r0
    cut = (r0.band_ts>9)* (r0.qual<25) \
        * (r0.delta_ts<9) #* (r0.pnorm>1e-15) * (r0.pindex>0.5) * (r0.pindex<4)
    r1 = r0[cut]
    print 'selected %d with standard cuts' %len(r1)
    if mindist==0: return r1
    # need to prune
    prunecut = findsource.prune(r1)
    print 'selected %d after prune within %.2f deg' % (prunecut.sum(), mindist)    
    return r1[prunecut]

class Pivot(object):

    def __init__(self,indir, outdir, name, 
                select_all=False, 
                dzc='dzc.xml',
                other_keys = None,
                **kwargs
        ):
        
        """
         indir: where to find the pickle files describing the sources, 
         outdir: targe where to setup the pivot files: expect to find a dzc file here with the names
         name:  name for the Pivot Collection
         dzc:  xml file with the Deepzoom collection, already set up in the target folder.
                If 'dzc.xml', the default, there must also be a folder 'dzc_files'
         select_all: if True, do not filter the sources
         others: a list of keys to pass to 
        """
        if not os.path.exists(outdir): 
            raise InvalidParameter('folder %s not found')
        self.outdir = outdir
        
        # get a recarray with all the data
        self.rec = selection(indir, select_all, other_keys=other_keys)

        # create the Collection
        self.collection=collection.Collection(name, outdir,  self.rec.name, dzc, **kwargs) 
        
        self.fill_default()
    def limit(self, v, a, b, nan=None):
        r = v[:]
        if nan is not None: r[np.isnan(v)]=nan
        if a is not None: r[v<a]=a
        if b is not None: r[v>b]=b
        return r

    def fill_default(self):
        uwc = self.rec
        col = self.collection
        # now fill it with default Facets
        sdir = map(SkyDir, uwc.ra, uwc.dec)
        glon = np.asarray([s.l() for s in sdir])
        cglon = glon.copy()
        cglon[glon>180] = cglon[glon>180]-360
        col.add_facet('ra',  'Number', 'F3', uwc.ra) 
        col.add_facet('dec', 'Number', 'F3', uwc.dec) 
        col.add_facet('glon', 'Number', 'F3', np.round(cglon,4)) 
        glat = np.asarray([s.b() for s in sdir])
        glat[np.abs(glat)<1e-3] = 0
        col.add_facet('glat', 'Number', 'F3',glat) 
        uwc.band_ts[uwc.band_ts<1]=0
        uwc.ts[uwc.ts<1]=0
        uwc.ts[uwc.ts>1e5]=1e5
        col.add_facet('band_ts', 'Number', 'F1', uwc.band_ts)
        col.add_facet('TS', 'Number', 'F1', self.limit(uwc.ts, 0, 1e5))  # should fix this!
        
        if 'isonorm' in uwc.dtype.names:
            # has diffuse normalization info (old style)
            col.add_facet('isonorm', 'Number', 'F2', uwc.isonorm)
            uwc.galnorm[uwc.galnorm<1e-3]=0
            uwc.isonorm[uwc.isonorm<1e-3]=0
            col.add_facet('galnorm', 'Number', 'F2', uwc.galnorm)
        
        col.add_facet('a', 'Number', 'F3', self.limit(uwc.a,0, 1, nan=2))
        col.add_facet('b', 'Number', 'F3', self.limit(uwc.b, 0,1, nan=2))
        col.add_facet('qual', 'Number', 'F1', self.limit(uwc.qual,0, 100, 100))
        col.add_facet('e_pivot', 'Number', 'F1', uwc.pivot_energy)
        col.add_facet('lognorm', 'Number', 'F3', self.limit(np.log10(uwc.pnorm), -15, -8, nan=-15) )
        col.add_facet('pindex',  'Number', 'F3', self.limit(uwc.pindex, 0,5) ) 

        col.add_facet('delta_ts', 'Number', 'F1', self.limit(uwc.delta_ts,0, 100, 100))
        
        if 'id_prob' in uwc.dtype.names:
            # has association data
            col.add_facet('id_prob', 'Number', 'F3', uwc.id_prob)
            col.add_facet('class', 'String', 'C', ['%3s'% s if s!='' else 'None' for s in uwc.aclass])
        inridge = (abs(cglon)<60)*(abs(glat)<1.)
        col.add_facet('ridge', 'String', 'C', inridge)
        col.add_facet('high lat', 'String', 'C', (abs(glat)>10))

        ## nearby sources 
        ## (commented out since takes a long time with large collections: need to fix)
        #nearlist = []
        #nearmax=1.0
        #nearcount = []
        #for s in uwc:
        #    near_s = near(s, uwc, nearmax)
        #    nearlist.append( None if len(near_s)==0 else \
        #        (','.join([t.name for t in near_s]), '#'+'&amp;'.join(['$SEARCH$=FL.%s'%t.name for t in [s]+near_s])))
        #    nearcount.append(len(near_s))
        #col.add_facet('Near', 'Link', 'C', nearlist)
        #col.add_facet('nearby sources', 'Number', 'F', nearcount)
        #
      
    def add_facet(self,  name, type, format, data, **kwargs):
        """ add  a facet """
        self.collection.add_facet( name, type, format, data, **kwargs)
        
    def add_related(self, related):
        self.collection.add_related(related)
        
    def write(self, outfile='pivot.cxml', id_offset=0):
        fulloutfile = os.path.join(self.outdir, outfile)

        print 'writing collection file with %d Items, %d Facets to "%s" ...' % \
            (self.collection.n, len(self.collection.facets), fulloutfile), 
        self.collection.write(fulloutfile, id_offset)
        print 'done!'

class MultiPivot(Pivot):
    """
    Subclass of Pivot implementing multiple Deepzoom collections
    (note that the Silverlight control does not support this feature :-( )
    
    """
    def __init__(self, rec,  dzc_list, outdir, name,  fill_default_facets=True, **kwargs):
        
        """
         rec:      record array with all info, expect it to have a name field
         outdir:   target where to set up the pivot files
         name:     name for the Pivot Collection
         dzc_list: xml files with the Deepzoom collections, already set with paths relative to the target folder.
        """
        if not os.path.exists(outdir): 
            raise InvalidParameter('folder %s not found')
        self.outdir = outdir
        self.rec = rec

        # create the Collection, actually a subclass that knows about mulitple dzc's
        self.collection=collection.MultiCollection(name, outdir,  self.rec.name, dzc_list,  **kwargs) 
        
        # add default Facets
        if fill_default_facets: self.fill_default()


if __name__=='__main__':
    from optparse import OptionParser
    usage = """\n
usage: %prog [options] indir outdir name\n
Generate pivot collection
    indir: where to find the pickle files describing the sources
    outdir: where to setup the pivot files: expect to find a dzc file here with the names
    name:  name for the Pivot Collection (surround by double quotes if more than one token) """
    parser = OptionParser(usage, version=version)
    parser.add_option('-a', '--all', help='keep all sources (otherwise apply filter)', action='store_true', dest='select_all',default=False)
    parser.add_option( '--noorigin', help='do not check for standard origin from name', action='store_true', dest='noorigin', default=False)
    options, args = parser.parse_args()
    if len(args)!=3: 
        parser.print_usage()
        sys.exit(-1)
    origin = None if options.noorigin else {'1F':'1FGL', 'PG':'PGW', 'MR':'MRF', 'UW':'UW'}
    Pivot(args[0],args[1],args[2], select_all=options.select_all, origin=origin).write()
 
