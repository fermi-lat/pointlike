"""
do pivot stuff 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/pivot.py,v 1.4 2010/05/11 19:05:00 burnett Exp $

"""
version='$Revision: 1.4 $'.split()[1]
from uw.utilities import collection
from uw.thb_roi import pipeline, catalog
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
    prunecut = catalog.prune(r1)
    print 'selected %d after prune within %.2f deg' % (prunecut.sum(), mindist)    
    return r1[prunecut]

class Pivot(object):

    def __init__(self,indir, outdir, name, 
                select_all=False, 
                dzc='dzc.xml',
                other_keys = None,
                icon=None,
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
        self.collection=collection.Collection(name, outdir,  self.rec.name, dzc,  icon=icon) 
        
        self.fill_default()

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
        col.add_facet('ts', 'Number', 'F1', uwc.ts)
        col.add_facet('isonorm', 'Number', 'F2', uwc.isonorm)
        uwc.galnorm[uwc.galnorm<1e-3]=0
        uwc.isonorm[uwc.isonorm<1e-3]=0
        col.add_facet('galnorm', 'Number', 'F2', uwc.galnorm)
        col.add_facet('a', 'Number', 'F3', uwc.a)
        col.add_facet('b', 'Number', 'F3', uwc.b)
        col.add_facet('qual', 'Number', 'F1', uwc.qual)
        col.add_facet('pivot_energy', 'Number', 'F1', uwc.pivot_energy)
        pnorm = uwc.pnorm
        pnorm[np.isinf(pnorm) * ( pnorm< 1e-20)] = 1e-20
        col.add_facet('lognorm', 'Number', 'F3', np.log10(pnorm) )
        col.add_facet('pindex',  'Number', 'F3', uwc.pindex) 

        dts = uwc.delta_ts
        dts[dts<1e-3]=0
        col.add_facet('delta_ts', 'Number', 'F1', dts)
        col.add_facet('id_prob', 'Number', 'F3', uwc.id_prob)
        inridge = (abs(cglon)<60)*(abs(glat)<1.)
        col.add_facet('class', 'String', 'C', ['%3s'% s if s!='' else 'None' for s in uwc.aclass])
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
        
    def write(self, outfile='pivot.cxml', id_offset=0):
        fulloutfile = os.path.join(self.outdir, outfile)

        print 'writing collection file to "%s" ...' % fulloutfile, 
        self.collection.write(fulloutfile, id_offset)
        print 'done!'

class MultiPivot(Pivot):
    """
    Subclass of Pivot implementing multiple Deepzoom collections
    
    """
    def __init__(self, rec,  dzc_list, outdir, name,  icon=None, ):
        
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
        self.collection=collection.MultiCollection(name, outdir,  self.rec.name, dzc_list,  icon=icon) 
        
        # add default Facets
        self.fill_default()


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
 
