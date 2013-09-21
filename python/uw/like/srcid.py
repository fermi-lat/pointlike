"""
Python support for source association, equivalent to the Fermi Science Tool gtsrcid
author:  Eric Wallace <wallacee@uw.edu>
"""

__version__ = "$Revision: 1.33 $"
#$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/srcid.py,v 1.33 2012/07/16 16:44:09 lande Exp $

import os
import sys
import math
import re
import operator
import glob

import numpy as np
import pyfits as pf

import skymaps
from uw.utilities import fitstools,keyword_options,path


class SrcidError(Exception): pass

conv95 = (-2*np.log(0.05))**.5
class SourceAssociation(object):
    """A class to associate LAT sources with multiwavelength counterparts."""

    defaults = (('catalog_dir',None,
                 'Path to counterpart catalogs, overriding "srcid_dir/cat"'),
                ('class_dir',None,
                 '''Path to directory containing class modules, overriding
                    "srcid_dir/classes"'''),
                ('verbosity',1,
                 '''How verbose to be: 0 for no output, 1 for normal output,
                    2 for extra output'''),
                ('quiet',False,'Set verbosity=0. DEPRECATED'))

    @keyword_options.decorate(defaults)
    def __init__(self,srcid_dir='$FERMI/catalog/srcid',**kwargs):
        """
        Arguments:
           srcid_dir : Path to a directory containing a subfolders 'cat',
                       containing counterpart catalogs, and 'classes',
                       containing python modules describing the catalogs.
        Keyword Arguments:
           catalog_dir: Path to counterpart catalogs, overriding 'srcid_dir/cat'
           class_dir: Path to directory containing class modules,
                      overriding 'srcid_dir/classes'
           verbosity : int, 0 for no output, 1 for normal output, 2 for extra output
           quiet: bool, deprecated, set verbosity = 0
        """
        keyword_options.process(self,kwargs)
        self.srcid_dir = path.expand(srcid_dir)
        if self.catalog_dir is not None:
            self.catalog_dir = path.expand(self.catalog_dir)
        else:
            self.catalog_dir = os.path.join(self.srcid_dir,'cat')
        if self.class_dir is not None:
            self.class_dir = path.expand(self.class_dir)
        else:
            self.class_dir = os.path.join(self.srcid_dir,'cat')
        if self.quiet: self.verbosity=0
        if not self.srcid_dir in sys.path:
            sys.path.insert(0,self.srcid_dir)
        self.class_list=None
        try:
            global classes
            import classes
            self.class_list = classes.__all__
            if self.verbosity:
                print('Available counterpart classes are:\n%s'%(' '.join(self.class_list or [''])))
        except ImportError:
            raise SrcidError('Cannot find class modules to import.')
        self.catalogs = {}
        self.sources = {}

    def id(self,position,error,**kwargs):
        """Find associations for each class in classes.

        Arguments:
            position:   A SkyDir corresponding to the source to be associated
            error:      Either the radius of a 1-sigma error circle, or a list
                        of major and minor axis and position angle for a 1-sigma
                        error ellipse
        Keyword Arguments:
            cpt_class[None]: string or [strings], Counterpart class(es) to use for
                          associations. Can be a single name or a list. If None,
                          use all availabile classes.
            name[None]: string, key under which to save association information
                        in self.sources. If None, don't save.
            unique[False]: bool, if True, compute probabilities for unique 
                           associations (one counterpart per LAT source)
            trap_mask[False]: bool, if True, use a trapezoidal mask for
                              computing local counterpart densities.
                              This is only included to get the same numbers
                              as gtsrcid, no reason to use it in general.
            accept_in_r95[False]: bool, if True, accept sources within r95 as
                                  associations, regardless of formal probability
        Return: a dictionary with keys=class, value = a list of (name, prob. skydir) tuples sorted by prob.
        """
        kw = dict(cpt_class = None,name=None,trap_mask=False,unique = False,accept_in_r95=True)
        for k,v in kw.items():
            kw[k] = kwargs.pop(k,v)
        if kwargs:
            raise SrcidError('\n'.join(['Unrecognized kwargs for method "id":']+
                            ["\t%s"%k for k in kwargs.keys()]+['']))
        if kw['cpt_class'] is None or kw['cpt_class']=='all':
            kw['cpt_class'] = self.class_list
        if hasattr(kw['cpt_class'],'__iter__'):
            associations = dict()
            id_kw = kw.copy()
            for c in kw['cpt_class']:
                id_kw['cpt_class'] = c
                ###################### recursive call. ugh #####################
                ##ass = SourceAssociation.id(self,position,error,**id_kw)
                ass = self.id(position,error,**id_kw)
                if ass:
                    associations[c] = ass
            return associations
        try:
            classes = __import__('classes',globals(),locals(),[kw['cpt_class']])
            class_module = getattr(classes,kw['cpt_class'])
        except ImportError, AttributeError:
            raise SrcidError("Counterpart class %s not found."%kw['cpt_class'])
        if not self.catalogs.has_key(kw['cpt_class']):
            self.catalogs[kw['cpt_class']] = Catalog(class_module,self.catalog_dir,verbosity=self.verbosity)
            #print 'loaded catalog %s' %class_module
        these = self.catalogs[kw['cpt_class']].associate(position,error,
                                                         trap_mask=kw['trap_mask'],
                                                         unique = kw['unique'],
                                                         accept_in_r95 = kw['accept_in_r95'])
        associations = [(a[0].name, a[1], a[0].skydir) for a in these if these]
        if kw['name'] is not None and associations!=[]:
            if not self.sources.has_key(kw['name']):
                self.sources[kw['name']] = dict()
            self.sources[kw['name']][kw['cpt_class']] = associations
        return associations

    def id_list(self,r,class_list = None,trap_mask=False,unique = False):
        """Perform associations on a recarray of sources.
        r: a recarray with columns name,ra,dec,a,b,ang
        class_list: list of counterpart classes
        return: dict with key = name, value = return from id(SkyDir(ra,dec),(a,b,ang))"""

        associations = {}
        for s in r:
            associations[s.name] = self.id(skymaps.SkyDir(s.ra,s.dec),(s.a,s.b,s.ang),
                                           cpt_class = class_list,name = s.name,
                                           trap_mask=trap_mask,unique = False)
        return associations

    def associate_catalog(self,catalog,**kwargs):
        """Cross-correlate a Fermi catalog with a given counterpart catalog

        Arguments:
            catalog: string, path to a LAT catalog for which to perform
                     associations
        Keyword Arguments:
            cpt_class[None]: string or [strings], Counterpart class(es) to use for
                          associations. Can be a single name or a list. If None,
                          use all availabile classes.
            name[None]: string, key under which to save association information
                        in self.sources. If None, don't save.
            unique[False]: bool, if True, compute probabilities for unique 
                           associations (one counterpart per LAT source)
            trap_mask[False]: bool, if True, use a trapezoidal mask for
                              computing local counterpart densities.
                              This is only included to get the same numbers
                              as gtsrcid, no reason to use it in general.
            elliptical_errors[True]: bool, if True, use the full elliptical
                                     error paraemeters, otherwise use a circle
                                     with radius equal to the semi-major axis
                                     of the ellipse. Provided for compatibility
                                     with older versions of gtsrcid, no reason
                                     to use it otherwise.
            accept_in_r95[False]: bool, if True, accept sources within r95 as
                                  associations, regardless of formal probability
        Return: a dictionary with keys=names, values=dicts returned by id()
        """
        #Read in Fermi catalog
        kw = dict(cpt_class = None,unique=True,
                  trap_mask=False,elliptical_errors=True,accept_in_r95=False)
        for k,v in kw.items():
            kw[k] = kwargs.pop(k,v)
        if kwargs:
            raise SrcidError('\n'.join(['Unrecognized kwargs for method "id":']+
                            ["\t%s"%k for k in kwargs.keys()]+['']))
        _cat = pf.open(catalog)
        dat = _cat['LAT_POINT_SOURCE_CATALOG'].data
        names = dat.Source_Name if 'Source_Name' in dat.names else dat.NickName
        try:
            ras,decs = dat.RA.astype('float'),dat.DEC.astype('float')
        except AttributeError:
            ras,decs = dat.RAJ2000.astype('float'),dat.DEJ2000.astype('float')
        if kw.pop('elliptical_errors'):
            maj_axes,min_axes = dat.Conf_95_SemiMajor/conv95,dat.Conf_95_SemiMinor/conv95
            angles = dat.Conf_95_PosAng
        else:
            maj_axes=min_axes=angles = dat.Conf_95_SemiMajor/conv95
        return dict(((name,self.id(skymaps.SkyDir(ra,dec),(a,b,ang),**kw))
                      for name,ra,dec,a,b,ang in zip(names,ras,decs,maj_axes,min_axes,angles)))

    def __str__(self):
        n = 0
        a = 0
        for v in self.sources.values():
            n+=1
            a += 1 if v else 0
        return 'SourceAssociation: %i sources, %i associated'%(n,a)

    def get_class(self,source):
        """Given a dictionary like self.sources, decide what class to ascribe the source to.  Partly guesswork!"""
        cat_list = source.keys()
        priority = '''bllac bzcat cgrabs crates crates_fom seyfert seyfert_rl qso agn
                    vcs galaxies pulsar_lat snr snr_ext pulsar_high pulsar_low pulsar_fom
                    msp pwn hmxb lmxb globular
                   '''.split()
        ass_class = ['bzb','bzcat']+['bzq']*3+['agn']*6+['LAT psr']+['snr']*2 + ['psr']*4 + ['pwn'] + ['hmxb'] + ['lmxb']+ ['glc'] 
        cls = None
        for c,a in zip(priority,ass_class):
            if c in cat_list:
                cls = a
                break
        if cls == 'bzcat':
            cls = source[cls][0][0][:3].lower()
        return cls

class Catalog(object):
    """A class to manage the relevant information from a FITS catalog."""
    def __new__(cls,class_module,catalog_dir,verbosity=1):
        gamma_catalogs = 'agile egr cosb eg3 fermi_bsl fermi_1fgl'.split()
        extended_catalogs = 'dwarfs snr_ext'.split()
        modname = class_module.__name__.split('.')[-1]
        if modname in gamma_catalogs:
            obj = object.__new__(GammaCatalog)
        elif modname in extended_catalogs:
            obj = object.__new__(ExtendedCatalog)
        else:
            obj = object.__new__(Catalog)
        return obj

    def init(self,class_module,catalog_dir,verbosity=1):
        self.verbosity = verbosity
        #self.class_module = self._get_class_module(class_file)
        if isinstance(class_module,str):
            self.class_module = __import__(class_module)
        else:
            self.class_module = class_module
        ## allows wild card in file specification
        cat_files = sorted(glob.glob(os.path.join(catalog_dir,self.class_module.catname)))
        #self.cat_file = os.path.join(catalog_dir,self.class_module.catname)
        assert len(cat_files)>0, 'File "%s" does not exist; module=%s' % (self.cat_file, class_module)
        self.cat_file = cat_files[-1]
        if self.verbosity > 1:
            print('Setting up catalog for source class "%s" from file "%s"'%(self.class_module.catid,self.cat_file))
        if hasattr(self.class_module,'name_prefix'):
            self.name_prefix = self.class_module.name_prefix
        else:
            self.name_prefix = ''
        self.coords = skymaps.SkyDir.EQUATORIAL
        self.prior = self.class_module.prob_prior
        self.prob_threshold = self.class_module.prob_thres
        self.max_counterparts = self.class_module.max_counterparts
        self.source_mask_radius = None #For selection of subset for association
        try:
            fits_cat = pf.open(self.cat_file)
        except IOError:
            raise CatalogError(self.cat_file,'Error opening catalog file. Not a FITS file?')
        #Find hdu with catalog information
        self.hdu = self._get_hdu(fits_cat)
        if self.hdu is None:
            raise CatalogError(self.cat_file,'No catalog information found.')
        names = [' '.join([self.name_prefix,x]).strip() for x in self._get_ids()]
        if names is None:
            raise CatalogError(self.cat_file,'Could not find column with source names')
        lons,lats = self._get_positions()
        if lons is None:
            raise CatalogError(self.cat_file,'Could not find columns with source positions')
        return names,lons,lats

    def __init__(self,class_module,catalog_dir,verbosity=1):
        self.names,lons,lats = self.init(class_module,catalog_dir,verbosity=verbosity)
        self.mask = self._make_selection()
        self.sources = np.array([CatalogSource(self,name,skymaps.SkyDir(lon,lat,self.coords))
                        for name,lon,lat in zip(self.names,lons,lats)])[self.mask]
        if self.coords == skymaps.SkyDir.GALACTIC:
            self.ras = np.array([source.skydir.ra() for source in self.sources])
            self.decs = np.array([source.skydir.dec() for source in self.sources])
        else:
            self.ras,self.decs = lons[self.mask],lats[self.mask]
        self._get_foms()

    def __iter__(self):
        return self.sources

    def _get_class_module(self,class_module):
        """Import and return module class_file"""
        return __import__(class_module)

    def _get_hdu(self,fits_cat):
        """Find and return HDU with catalog information."""
        #First check for HDUS with CAT-NAME or EXTNAME in header.
        for hdu in fits_cat:
            cards = hdu.header.ascardlist()
            try:
                self.cat_name = cards['CAT-NAME'].value
                return hdu
            except KeyError:
                try:
                    self.cat_name = cards['EXTNAME'].value
                    return hdu
                except KeyError:
                    self.cat_name = self.class_module.catid
                    pass
        #No CAT-NAME or EXTNAME found, just return second HDU
        if len(fits_cat)>=2:
            return fits_cat[1]
        #Only one HDU, no name info, return None
        return

    def _get_ids(self):
        """Find source name information and return as a list."""
        name_key = ''
        cards = self.hdu.header.ascardlist()
        #First check for UCD in header
        for card in cards:
            if card.key[:5]=='TBUCD' and card.value in ['ID_MAIN','meta.id;meta.main']:
                name_key = cards['TTYPE'+card.key[5:8]].value
                break
            #Sometimes UCDs are declared in comments
            #May be fragile - depends on specific format for comments as in gamma-egr catalog
            value_comment = card.ascardimage().split('/')
            if len(value_comment)>1:
                comment = value_comment[1]
                ucd_string = comment[comment.find('UCD'):].split()
                if ucd_string:
                    try:
                        if ucd_string[0].split('=')[1].strip('.')=='ID_MAIN':
                            name_key = cards[''.join(['TTYPE',card.key[5:8]])].value
                            break
                    except IndexError:
                        pass
            if card.key[:5]=='TTYPE' and card.value.upper() in ['NAME','ID','PSR_NAME','SOURCE_NAME']:
                name_key = card.value
                break
        try:
            return self.hdu.data.field(name_key)
        except KeyError:
            return

    def _get_positions(self):
        """Find columns containing position info and return a list of SkyDirs"""

        cards = self.hdu.header.ascardlist()
        ucds = cards.filterList('TBUCD*')
        ttypes = cards.filterList('TTYPE*')
        lon_key = lat_key = ''
        if not lon_key:
            if 'POS_EQ_RA_MAIN' in ucds.values():
                ucd = ucds.keys()[ucds.values().index('POS_EQ_RA_MAIN')]
                lon_key = ttypes[''.join(['TTYPE',ucd[5:8]])].value
                #Assumes that if POS_EQ_RA_MAIN exists, POS_EQ_DEC_MAIN does too.
                ucd = ucds.keys()[ucds.values().index('POS_EQ_DEC_MAIN')]
                lat_key = ttypes[''.join(['TTYPE',ucd[5:8]])].value
            elif 'RAdeg' in ttypes.values():
                lon_key = ttypes[ttypes.keys()[ttypes.values().index('RAdeg')]].value
                lat_key = ttypes[ttypes.keys()[ttypes.values().index('DEdeg')]].value
            elif '_RAJ2000' in ttypes.values():
                lon_key = ttypes[ttypes.keys()[ttypes.values().index('_RAJ2000')]].value
                lat_key = ttypes[ttypes.keys()[ttypes.values().index('_DEJ2000')]].value
            elif 'RAJ2000' in ttypes.values():
                lon_key = ttypes[ttypes.keys()[ttypes.values().index('RAJ2000')]].value
                lat_key = ttypes[ttypes.keys()[ttypes.values().index('DEJ2000')]].value
            elif 'RA' in ttypes.values():
                lon_key = ttypes[ttypes.keys()[ttypes.values().index('RA')]].value
                try:
                    lat_key = ttypes[ttypes.keys()[ttypes.values().index('DE')]].value
                except ValueError:
                    lat_key = ttypes[ttypes.keys()[ttypes.values().index('DEC')]].value
        if not lon_key:
            self.coords = skymaps.SkyDir.GALACTIC
            if 'POS_GAL_LON' in ucds.values():
                lon_key = ucds.keys()[ucds.values().index('POS_GAL_LON')]
                lat_key = ucds.keys()[ucds.values().index('POS_GAL_LAT')]
            elif '_GLON' in ttypes.values():
                lon_key = ttypes[ttypes.keys()[ttypes.values().index('_GLON')]].value
                lat_key = ttypes[ttypes.keys()[ttypes.values().index('_GLAT')]].value
            elif 'GLON' in ttypes.values():
                lon_key = ttypes[ttypes.keys()[ttypes.values().index('GLON')]].value
                lat_key = ttypes[ttypes.keys()[ttypes.values().index('GLAT')]].value
        if lon_key:
            return (self.hdu.data.field(lon_key).astype('float'),
                    self.hdu.data.field(lat_key).astype('float'))
        else:
            return

    def _make_selection(self):
        """Make selections specified in class module."""
        selections = [x for x in self.class_module.selection if (x != '' and 'ANGSEP' not in x)]
        dat = self.hdu.data
        mask = np.array([True]*len(dat))
        if not selections: return mask
        catid_pattern = re.compile('@%s_([A-Za-z0-9]+)'%self.class_module.catid.upper())
        fields = {}
        for sel in selections:
            field_names = catid_pattern.findall(sel)
            for f in field_names:
                field = dat.field(f)
                sel = catid_pattern.sub(f,sel,1)
                defnull_pattern = re.compile('DEFNULL\(%s,([0-9e\+\.]+)\)'%f)
                dn_matches = defnull_pattern.findall(sel)
                for dnm in dn_matches:
                    field[np.isnan(field)]=dnm
                fields[f]=field
                sel = defnull_pattern.sub(f,sel)
                sel = sel.replace(f,'fields["%s"]'%f)
            if '||' in sel:
                sel = 'np.logical_or(%s)'%(sel.replace('||',',').strip('()'))
            elif '&&' in sel:
                sel = 'np.logical_and(%s)'%(sel.replace('&&',',').strip('()'))
            mask = np.logical_and(mask,eval(sel))
        return mask


    def _get_foms(self):
        """Compute figure of merit for sources, as specified in class module."""
        fom = self.class_module.figure_of_merit
        dat = self.hdu.data
        if not fom:
            return
        catid_pattern = re.compile('@%s_([A-Za-z0-9_]+)'%self.class_module.catid.upper())
        fields = {'Name':self._get_ids()}
        field_names = catid_pattern.findall(fom)
        for f in field_names:
            field = dat.field(f)
            fom = catid_pattern.sub(f,fom,1)
            #Treat DEFNULL in a log slightly differently
            #This is an awful kludge to deal with the pulsar_fom catalog.
            #Hopefully it won't break anything in the future.
            defnull_pattern = re.compile('(?<!LOG10\()DEFNULL\(%s,([0-9e\.\+]+)\)'%f)
            dn_matches = defnull_pattern.findall(fom)
            for dnm in dn_matches:
                field[np.isnan(field)] = dnm
            defnull_pattern2 = re.compile('(?<=LOG10\()DEFNULL\(%s,([0-9e\.\+]+)\)'%f)
            dn_matches2 = defnull_pattern2.findall(fom)
            for dnm in dn_matches2:
                field[np.logical_or(np.isnan(field),field==0)] = dnm
            fields[f] = field
            fom = defnull_pattern.sub(f,fom)
            fom = defnull_pattern2.sub(f,fom)
            fom = fom.replace(f,'fields["%s"]'%f)
        fom = fom.replace('exp','np.exp')
        fom = fom.replace('LOG10','np.log10')
        fom_dict = dict(zip([' '.join([self.name_prefix,x]).strip() for x in fields['Name']],eval(fom)))
        for source in self.sources:
            source.fom = fom_dict[source.name]

    def select_circle(self,position,radius,trapezoid=False):
        """Return an array of CatalogSources within radius degrees of position.

        Arguments:
            position    : SkyDir for center of selection region.
            radius      : radius of selection region.
        """
        ras,decs = self.ras,self.decs
        if trapezoid:
            tmask = trap_mask(self.ras,self.decs,position,radius)
            sources = self.sources[tmask]
            rmask = fitstools.rad_mask(self.ras[tmask],self.decs[tmask],position,radius,mask_only=True)
            return sources[rmask]
        else:
            rmask = fitstools.rad_mask(self.ras,self.decs,position,radius,mask_only=True)
            return self.sources[rmask]

    def local_density(self,position,radius=4,fom=1.0,trap_mask=False):
        """Return the local density of catalog sources in a radius-degree region about position.

        Only counts sources with figures of merit >= fom. The default fom for CatalogSources should be 1.,
        so the default fom=0 should not cut anything out.  However, this method ought to be independent
        of the implementation of the fom in CatalogSource, so this should get refactored at some point."""
        n_sources = sum([1. for source in self.select_circle(position,radius,trapezoid = trap_mask) if
                         source.fom >= fom])
        #If no sources within radius, set n_sources = 1 to give lower limit on density
        #Maybe better to expand the radius in this case?
        if n_sources < 1 : n_sources = 1
        solid_angle = (1-math.cos(np.radians(radius)))*2*math.pi
        solid_angle = np.degrees(np.degrees(solid_angle))
        return n_sources/solid_angle

    def associate(self,position,error_ellipse,unique = True,trap_mask=False,accept_in_r95=False):
        """Given a skydir and error ellipse, return associations.

        Arguments:
            position      : A SkyDir representing the location of the source to be associated.
            error_ellipse : Sequence of length 3 representing major and minor axes and position angle of an
                          error ellipse in degrees
            accept_in_r95 : If True, accept anything within r95 as an association

        Returns all sources with posterior probability greater than prob_threshold,
        sorted by probability.
        """

        try:
            assert(len(error_ellipse)==3)
        except TypeError:
            if self.verbosity > 0:
                print("Got scalar instead of sequence for error ellipse: Assuming this is 1-sigma error radius in degrees")
            error_ellipse = [error_ellipse]*2+[0]
        except AssertionError:
            print("Wrong length for error_ellipse: Needed length 3, got %i"%len(error_ellipse))
            return
        if self.source_mask_radius is None:
            self.source_mask_radius = error_ellipse[0]*conv95*5
        if np.isnan(self.source_mask_radius) or self.source_mask_radius == 0.:
            self.source_mask_radius = 1.
        #filter sources by position, ~5-sigma radius
        sources = self.select_circle(position,self.source_mask_radius,trapezoid=trap_mask)
        post_probs = np.array([source.posterior_probability(position,error_ellipse,trap_mask=trap_mask)
                               for source in sources])
        #If desired, require no more than 1 counterpart per LAT source.
        if unique:
            inv_probs = 1-post_probs
            phk = np.zeros(post_probs.shape) #P(Hk) as defined in 1FGL paper
            norm = inv_probs.prod()
            for i in xrange(post_probs.shape[0]):
                mask = np.zeros(post_probs.shape,dtype='bool')
                mask[i] = True
                phk[i] = post_probs[mask]*inv_probs[~mask].prod()
                norm += phk[i]
            post_probs = phk/norm
        #return sources above threshold with posterior probability, sorted by posterior probability
        if accept_in_r95:
            source_list = [(source,prob) for prob,source in zip(post_probs,sources) if prob > self.prob_threshold or np.degrees(source.skydir.difference(position))<error_ellipse[0]*conv95]
        else:
            source_list = [(source,prob) for prob,source in zip(post_probs,sources) if prob > self.prob_threshold]
        source_list.sort(key = lambda x:-x[1])
        #source_dict = dict((el[1].name,el) for el in source_list[:self.max_counterparts])
        return source_list

class GammaCatalog(Catalog):
    """A catalog of gamma-ray sources (i.e. sources with error circles comparable to LAT)"""

    def __init__(self,class_module,catalog_dir,verbosity = 1):
        self.names,lons,lats = self.init(class_module,catalog_dir,verbosity = verbosity)
        errors = self.get_position_errors()
        self.source_mask_radius = 3*max(errors)
        self.mask = self._make_selection()
        self.sources = np.array([GammaRaySource(self,name,skymaps.SkyDir(lon,lat,self.coords),error)
                        for name,lon,lat,error in zip(self.names,lons,lats,errors)])[self.mask]
        #Save ras and decs for use in select_circle. MUCH faster than using the SkyDirs.
        if self.coords == skymaps.SkyDir.GALACTIC:
            self.ras = np.array([source.skydir.ra() for source in self.sources])
            self.decs = np.array([source.skydir.dec() for source in self.sources])
        else:
            self.ras,self.decs = lons[self.mask],lats[self.mask]

    def get_position_errors(self):
        q = [x for x in self.class_module.new_quantity if
                    (('LAT' not in x) and ('ECOM' not in x))][0]
        patt = re.compile('@%s_([A-Za-z0-9_]+)'%self.class_module.catid.upper())
        lhs,rhs = q.split('=')
        match = patt.search(rhs)
        if match:
            error_field = match.groups()[0]
            rhs = patt.sub('self.hdu.data.field("%s")'%error_field,rhs)
            return eval(rhs.split('*')[0])
        else:
            raise CatalogError(self.cat_file,'Could not find position uncertainties.')

class ExtendedCatalog(Catalog):
    """A catalog of extended sources"""

    def __init__(self,class_module,catalog_dir,verbosity = 1):
        self.names,lons,lats= self.init(class_module,catalog_dir,verbosity = verbosity)
        radii = self.get_radii()
        self.source_mask_radius = max(radii)*3
        self.mask = self._make_selection()
        self.sources = np.array([ExtendedSource(self,name,skymaps.SkyDir(lon,lat,self.coords),radius)
                                 for name,lon,lat,radius in zip(self.names,lons,lats,radii)])[self.mask]
        if self.coords == skymaps.SkyDir.GALACTIC:
            self.ras = np.array([source.skydir.ra() for source in self.sources])
            self.decs = np.array([source.skydir.dec() for source in self.sources])
        else:
            self.ras,self.decs = lons[self.mask],lats[self.mask]

    def get_radii(self):
        q = self.class_module.new_quantity[0]
        lhs,rhs = q.split('=')
        terms = rhs.split('+')
        rad_term = [term.strip() for term in terms if '@%s'%self.class_module.catid.upper() in term][0]
        radius = rad_term.replace('@%s_'%self.class_module.catid.upper(),'')
        num,denom = radius.split('/')
        num = 'self.hdu.data.field("%s")'%num
        radius = '/'.join([num,denom])
        return eval(radius)

class CatalogSource(object):
    """A class representing a catalog source."""
    def __init__(self,catalog,name,skydir):
        self.catalog = catalog
        self.name = name
        self.skydir = skydir
        self.fom = 1.

    def __str__(self):
        return '\t'.join([self.catalog.cat_name,self.name,str(self.skydir.ra()),str(self.skydir.dec())])

    def angular_separation(self,other_skydir):
        """Calculate angular separation between this source and other_skydir"""
        return np.degrees(self.skydir.difference(other_skydir))

    def position_angle(self,other_skydir):
        """Calculate other_skydir angle _from_ other_skydir _to_ this source."""
        ra1,dec1,ra2,dec2 = np.radians([other_skydir.ra(),other_skydir.dec(),
                            self.skydir.ra(),self.skydir.dec()])
        denom = math.sin(dec1)*math.cos(ra1-ra2) - math.cos(dec1)*math.tan(dec2)
        return np.degrees(math.atan2(math.sin(ra1-ra2), denom))

    def delta_logl(self,position,error_ellipse):
        """Compute Delta(logl) for this association with source at position, with error_ellipse

        Arguments:
            position : SkyDir representing position of other source
            error_ellipse : 3-tuple representing error ellipse of other source in degrees"""

        #error_ellipse = np.radians(error_ellipse)
        phi = np.radians(self.position_angle(position)-error_ellipse[2])
        angsep = self.angular_separation(position)
        return .5*angsep**2*((math.cos(phi)/error_ellipse[0])**2 +
                               (math.sin(phi)/error_ellipse[1])**2)

    def positional_likelihood(self,position,error_ellipse):
        """Likelihood for point source at position with error_ellipse to be associated with this source.

        Arguments:
            position : SkyDir representing position of other source
            error_ellipse : 3-tuple representing error ellipse of other source in degrees"""

        norm = 2.*math.pi*operator.mul(*error_ellipse[:2])
        return math.exp(-self.delta_logl(position,error_ellipse))/norm

    def chance_probability(self,position,radius = 4,trap_mask=False):
        """Probability of chance association"""
        return self.catalog.local_density(position,radius=radius,
                                          fom=self.fom,trap_mask=trap_mask)

    def prior_probability(self):
        """Prior probability for association"""
        return self.catalog.prior

    def posterior_probability(self,position,error_ellipse,trap_mask=False):
        """Posterior probability for association with given position and error ellipse"""
        #Kludge to associate lat pulsars, which have no errors in catalog
        if np.isnan(error_ellipse)[0] and position.difference(self.skydir)<=1e-5:
            return self.catalog.prob_threshold+1e-5
        ring_rad = error_ellipse[0]*(-2*np.log(1-.95))**.5*5/2
        if ring_rad <4:
            ring_rad = 4
        denom = self.positional_likelihood(position,error_ellipse)*self.prior_probability()
        if denom==0: return 0
        arg = (self.chance_probability(position,radius=ring_rad,trap_mask=trap_mask)*
               (1-self.prior_probability())/denom)
        prob = self.fom/(1.+arg)
        if prob > 1e-8:
            return prob
        else:
            return 0.

class GammaRaySource(CatalogSource):
    """A source from a gamma-ray catalog."""

    def __init__(self,catalog,name,skydir,error):
        super(GammaRaySource,self).__init__(catalog,name,skydir)
        self.error = error

    def posterior_probability(self,position,error_ellipse,trap_mask = True):
        error_ellipse = np.where(np.isnan(error_ellipse),0,error_ellipse)
        if self.combined_error(error_ellipse) >= self.angular_separation(position):
            return self.catalog.prob_threshold + 1e-5
        else:
            return 0.0

    def combined_error(self,error_ellipse):
        #note 2.45 factor - conversion from 1-sigma ellipse to 95%
        return ((error_ellipse[0]*conv95)**2 + self.error**2)**.5

class ExtendedSource(CatalogSource):
    """An extended catalog source"""

    def __init__(self,catalog,name,skydir,radius):
        super(ExtendedSource,self).__init__(catalog,name,skydir)
        self.radius = radius

    def posterior_probability(self,position,error_ellipse,trap_mask=True):
        #note 2.45 factor - conversion from 1-sigma ellipse to 95%
        error_ellipse = np.where(np.isnan(error_ellipse),0,error_ellipse)
        if self.angular_separation(position)<=(error_ellipse[0]*conv95+self.radius):
            return self.catalog.prob_threshold + 1e-5
        else:
            return 0.0

class SrcidError(Exception):
    """Exception class for general problems."""
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return self.message

class CatalogError(Exception):
    """Exception class for problems with a catalog."""
    def __init__(self,catalog,message):
        self.catalog = catalog
        self.message = message
    def __str__(self):
        return 'In catalog %s:\n\t%s'%(self.catalog,self.message)

def trap_mask(ras,decs,cut_dir,radius):
    dec_min = max(cut_dir.dec()-radius,-90)
    dec_max = min(cut_dir.dec()+radius,90)
    ra_rad = radius/np.cos(np.radians(cut_dir.dec()))
    ra_min,ra_max = cut_dir.ra()-ra_rad,cut_dir.ra()+ra_rad
    ra_min = ra_min if ra_min>0 else ra_min+360
    ra_max = ra_max if ra_max<360 else ra_max-360
    dec_mask = np.logical_and(decs>=dec_min,decs<=dec_max)
    ra_mask = (np.logical_and(ras>=ra_min,ras<=ra_max) if ra_max>ra_min else
              np.logical_or(ras<min(ra_min,ra_max),ras>max(ra_min,ra_max)))
    return np.logical_and(dec_mask,ra_mask)

def test():

    assoc = SourceAssociation()
    #3C 454.3
    pos, error = skymaps.SkyDir(343.495,16.149), .016/2.45*1.51
    associations = assoc.id(pos,[error,error,0],name = '3C 454.3')
    print(associations)
    #for cat,ass in associations.items():
    #    print 'Associations in %s:'%cat
    #    print ass
    #print('\n'.join([str(x[1]) for x in assoc.id(pos,error,'obj-blazar-crates',.33,.8)]))
    #Couldn't find elliptical errors, but want to test input for error.
    #error = (error,error,0.0)
    #print('\n'.join([str(x[1]) for x in assoc.id(pos,error,'obj-blazar-crates',.33,.8)]))
    #Test for correct failure for wrong length list
    #error = [.5]*4
    #print(assoc.id(pos,error,'obj-agn',.039,.9735))
    
if __name__=='__main__':
    test()
