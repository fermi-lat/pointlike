"""
supplemental setup of ROI
----------------------------------
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/thb_roi/roi_setup.py,v 1.11 2010/08/29 20:19:58 burnett Exp $

These are near-duplicates of the classes with the same name in uw.like, but modifed for the interactive access

"""
import numpy as np
import os, pyfits, types, pickle

from uw.utilities import makerec, keyword_options, xml_parsers
from uw.like import Models,  pointspec_helpers
from skymaps import SkyDir

     
class CatalogManager(object):
    """Read a  catalogue and use it to set up a source list for a ROI."""
    defaults = (
        ('prune_radius',  0.10,  'deg; in a merge, consider sources closer than this duplicates\n'
                                ' or to exclude all sources within a larger region for separate specification'),
        ('free_radius',   1,      'deg; sources within this distance have free spectral parameters'),
        ('min_flux',     2e-9,   'ph/cm2/s; minimum flux for sources beyond a certain proximity'),
        ('max_distance',  5,      'deg; distance inside which sources are returned regardless of flux'),
        ('min_ts',       25,     ''),
        ('pulsar_dict',   None,   ' optional pulsar dictionary with ExpCutoff parameters to use'),
        ('point_source',  None,   ' class that creates point source objects: if None, use pointspec_helpers.PointSource'),
        ('quiet',         False,   ''),
        ('verbose',       False,   ''),
        )

    @keyword_options.decorate(defaults)
    def __init__(self,catalog_file, *args, **kwargs):
        """ Create the catalog: 
            catalog_file: a Fermi-LAT FITS format catalog
            """
        keyword_options.process(self, kwargs)
        if self.point_source is None: self.point_source = pointspec_helpers.PointSource
        if not self.quiet: 
            print 'creating a CatalogManager from %s...' %catalog_file
        cdata = pyfits.open(catalog_file)[1].data
        try:
            ts   = np.asarray(cdata.field('Test_Statistic'),float)
        except KeyError:
            ts   = np.asarray(cdata.field('Signif_Avg'))**2
        good = ts>self.min_ts
             
        ras  = np.asarray(cdata.field('RA'),float)[good]
        decs = np.asarray(cdata.field('DEC'),float)[good]
        pens = cdata.field('PIVOT_ENERGY')[good]
        n0s  = cdata.field('FLUX_DENSITY')[good]
        inds = cdata.field('SPECTRAL_INDEX')[good]
        try:
            cuts = cdata.field('Cutoff_Energy')[good]
        except KeyError:
            cuts = good.sum()*[np.nan]
        inds = np.where(inds > 0, inds, -inds)
        try:
            self.names  = np.asarray(cdata.field('Source_Name'))[good]
        except KeyError:
            self.names  = np.asarray(cdata.field('NickName'))[good]

        self.dirs   = map(SkyDir,ras,decs)
        
        pdkeys = [] 
        if self.pulsar_dict is not None:
            if type(self.pulsar_dict) == types.StringType:
                self.pulsar_dict = pickle.load(open(self.pulsar_dict))
            pdkeys = self.pulsar_dict.keys()
        def load_model(name, n0,ind, pen, cut):
            #if name[:3]=='PSR': assert False, 'breakpoint'
            if name not in pdkeys:
                return Models.PowerLaw(p=[n0,ind],e0=pen) if np.isnan(cut) or cut==0 else\
                       Models.ExpCutoff(p=[n0,ind,cut],e0=pen)
            # Use pulsar_dict to override catalog           
            psr = self.pulsar_dict[name]
            if psr['TS']<100:return Models.PowerLaw(p=[n0,ind],e0=pen)
            stat = psr['stat'][0]
            if self.verbose: 
                print ('replacing catalog fits for source %s, par='+3*'%12.2e')\
                    % ((name,) + tuple(stat) )
            return Models.ExpCutoff(p=stat)
            
        self.models = np.asarray([load_model(name.strip(),n0,ind,pen,cut) for name,n0,ind,pen,cut in zip(self.names,n0s,inds,pens,cuts)])
        if not self.quiet: 
            print 'Loaded %d sources  for roi backgrounds' % (len(cdata),)
            if np.sum(-good)>0:
                print '\tNot using %d entries with TS<25' % np.sum(-good)

            if self.pulsar_dict is not None:
                print '\tUsing a pulsar dictionary with %d entries' % len(pdkeys)
    
    def append(self, acat):
        """ 
        append sources found in an auxilliary list of sources
        
        acat: a recarray with fields name, ra, dec, pnorm, pindex
        """
        names  = acat.name
        dirs   = map(SkyDir,acat.ra,acat.dec)
        models = [Models.PowerLaw(p=[n0,ind],e0=1000.) for n0,ind in zip(acat.pnorm,acat.pindex)]
        self.names = np.hstack((self.names, names))
        self.dirs= self.dirs+dirs
        self.models = np.hstack((self.models, models))
        if not self.quiet: print 'Added %d sources to roi background catalog, total now: %d'\
            % (len(acat), len(self.names))
 
        
    def __call__(self,skydir, radius, **kwargs):
        """ return sources, as a list of PointSource objects, within radius of skydir
            set as free those within free_radius
            exclude those within prune_radius
        """
        quiet = kwargs.pop('quiet', self.quiet)
        free_radius = kwargs.pop('free_radius', self.free_radius)
        prune_radius = kwargs.pop('prune_radius', self.prune_radius)
        diffs   = np.degrees(np.asarray([skydir.difference(d) for d in self.dirs]))
        #mask    = ((diffs < radius)&(self.ts > self.min_ts)) & \
        #          ((self.fluxes > self.min_flux)|(diffs < self.max_distance))         
        mask    = (diffs < radius) 
        diffs   = diffs[mask]
        sorting = np.argsort(diffs)

        #sort SkyDirs -- pesky vector behavior...
        dirs    = [x for i,x in enumerate(self.dirs) if mask[i]]
        dirs    = [dirs[x] for x in sorting]

        names   = self.names[mask][sorting]
        models  = [x.copy() for x in self.models[mask][sorting]]

        # now select those in annulus between prune_radius and free_radius for inclusion in the model
        exclude = (diffs[sorting] < prune_radius) 
        numex = exclude.sum()
        makefree = (diffs[sorting] < free_radius) 
        point_sources = map(self.point_source,dirs,names,models,makefree)
        if not quiet:
            print '...selected %d sources within %.1f deg for roi'   % (len(point_sources), radius)
            print '...selected %d sources within %.1f deg for refit' % (makefree.sum(), free_radius)
            if numex>0:
                print '...excluded %d sources within %.1f deg : %s ' %\
            (numex, prune_radius, ', '.join([s.name.strip() for s in point_sources[:numex]]))
        
        self.exclude= point_sources[:numex]
        return point_sources[numex:]

    def merge_lists(self,skydir,radius=15,user_list=None):
        """Get a list of catalog sources and merge it with an (optional) list of PointSource objects
         provided by the user.  In case of duplicates (closer than prune_radius), the user object
         takes precedence.
         """
        cat_list = self(skydir,radius)
        if user_list is None: return cat_list

        from collections import deque
        merged_list = deque(user_list)

        for ncps,cps in enumerate(cat_list):
            merged_list.append(cps)
            for ups in user_list:
                if np.degrees(ups.skydir.difference(cps.skydir)) < self.prune_radius:
                   merged_list.pop(); break

        merged_list = list(merged_list)
        merged_list.sort(key = lambda ps: ps.skydir.difference(skydir))

        return merged_list
        
    def source_recarray(self):
        """ return a recarry of the full list of sources used for models
        """
        return np.rec.fromarrays(
                [   [name.strip() for name in self.names], 
                    [s.ra() for s in self.dirs], 
                    [s.dec() for s in self.dirs], 
                    [10**m.p[0] for m in self.models],
                    [10**m.p[1] for m in self.models],
                    [10**m.p[2] if len(m.p)==3 else np.nan for m in self.models ],
                    [m.e0 for m in self.models],
                ],
                names = 'name ra dec pnorm pindex cutoff pivot'.split(),
           )

    def write_reg_file(self, filename, color='green'):
        """ generate a 'reg' file from the catalog, write to outfile
        """
        catrec = self.source_recarray()
        have_ellipse = 'Conf_95_SemiMajor' in catrec.dtype.names
        out = open(filename, 'w')
        print >>out, "# Region file format: DS9 version 4.0 global color=%s" % color
        for s in catrec:
            if have_ellipse:
                print >>out, "fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) #text={%s}" % \
                                (s.ra,s,dec,
                                  s.Conf_95_SemiMinor,Conf_95_SemiMajor,Conf_95_PosAng,
                                  s.name)
            else:
                print >> out, "fk5; point(%.4f, %.4f) # point=cross text={%s}" %\
                                (s.ra, s.dec, s.name)
        out.close()

    def write_xml_file(self, filename, title='source_library'):
        
        point_sources = map(self.point_source,self.dirs,self.names,self.models)
        stacks= xml_parsers.unparse_point_sources(point_sources)
        xml_parsers.writeXML(stacks, filename, title=title)

if __name__=='__main__':
    pass
