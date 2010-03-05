import os
import sys
import math
import re
from operator import mul
from glob import glob
import numpy as np
import pyfits as pf
from skymaps import SkyDir
from uw.utilities.fitstools import trap_mask, rad_mask


class SourceAssociation(object):
    """A class to find association probabilities for sources with a given set of counterpart catalogs."""

    def __init__(self,catdir):
        self.catdir = catdir
        cat_files = glob(os.path.join(self.catdir,'*.fits'))
        self.catalogs = {}

    def id(self,position,error,classes,class_dir = None):
        """Find associations for each class in classes.

        Arguments:
            position:   A SkyDir corresponding to the source to be associated
            error:      Either the radius of a 1-sigma error circle, or a list
                        of major and minor axis and position angle for a 1-sigma
                        error ellipse
            classes:    List of source classes to be associated.  Names should correspond
                        to python modules in class_dir
            class_dir:  Directory to look in for source class modules. Defaults to
                        self.catdir/../classes.

        For now, just gets the prior and threshold for each catalog and calls single_catalog_id.
        This means that all catalogs are treated in the same way, which will not properly handle
        catalogs of extended sources, or point sources with error circles comparable to those of
        LAT sources.
        """
        if class_dir is None:
            class_dir = os.path.join(self.catdir,os.path.pardir,'classes')
        sys.path.insert(0,class_dir)
        #If classes is not a list, make it one
        if not hasattr(classes,'__iter__'): classes = [classes]
        class_files = [os.path.join(class_dir,cls) for cls in classes]
        associations = {}
        for cls,cf in zip(classes,class_files):
            if not self.catalogs.has_key(cls):
                self.catalogs[cls] = Catalog(cf)
            these = self.single_catalog_id(position,error,self.catalogs[cls])
            associations[cls] = these
        return associations

    def single_catalog_id(self,position,error_ellipse,catalog):
        """Given a skydir and error ellipse, return associations from catalogs.

        Arguments:
            position      : A SkyDir representing the location of the source to be associated.
            error_ellipse : Sequence of length 3 representing major and minor axes and position angle of an
                          error ellipse in degrees
            catalog       : A Catalog object in which to look for associations

        Returns all sources with posterior probability greater than prob_threshold,
        sorted by probability.

        Maybe make this a method on Catalog?
        """

        try:
            assert(len(error_ellipse)==3)
        except TypeError:
            print("Got scalar instead of sequence for error ellipse: Assuming this is r68 in degrees")
            error_ellipse = [error_ellipse]*2+[0]
        except AssertionError:
            "Wrong length for error_ellipse: Needed length 3, got %i"%len(error_ellipse)
            return

        #filter sources by position, ~5-sigma radius
        sources = catalog.select_circle(position,error_ellipse[0]*5)
        post_probs = [source.posterior_probability(position,error_ellipse) for source in sources]
        #return sources above threshold with posterior probability, sorted by posterior probability
        source_list = [(prob,source) for prob,source in zip(post_probs,sources) if prob > catalog.prob_threshold]
        source_list.sort()
        source_dict = dict((el[1].name,el) for el in source_list[:catalog.max_counterparts])
        return source_dict


class Catalog(object):
    """A class to manage the relevant information from a FITS catalog."""

    def __new__(cls,class_file):
        gamma_catalogs = ['%s'%g for g in 'agile egr cosb eg3'.split()]
        extended_catalogs = ['%s'%e for e in 'dwarf snr_ext'.split()]
        cf = os.path.splitext(os.path.basename(class_file))[0]
        if cf in gamma_catalogs:
            obj = object.__new__(GammaCatalog)
        elif cf in extended_catalogs:
            obj = object.__new__(ExtendedCatalog)
        else:
            obj = object.__new__(Catalog)
        return obj

    def init(self,class_file):
        self.class_module = self._get_class_module(class_file)
        self.cat_file = os.path.join(os.path.dirname(class_file),os.path.pardir,'cat',self.class_module.catname)
        print('Setting up catalog for source class "%s" from file "%s"'%(self.class_module.catid,self.cat_file))
        self.coords = SkyDir.EQUATORIAL
        self.prior = self.class_module.prob_prior
        self.prob_threshold = self.class_module.prob_thres
        self.max_counterparts = self.class_module.max_counterparts
        try:
            fits_cat = pf.open(self.cat_file)
        except IOError:
            raise CatalogError(self.cat_file,'Error opening catalog file. Not a FITS file?')
        #Find hdu with catalog information
        self.hdu = self._get_hdu(fits_cat)
        if self.hdu is None:
            raise CatalogException(self.cat_file,'No catalog information found.')
        names = self._get_ids()
        if names is None:
            raise CatalogException(self.cat_file,'Could not find column with source names')
        lons,lats = self._get_positions()
        if lons is None:
            raise CatalogException(self.cat_file,'Could not find columns with source positions')
        return names,lons,lats

    def __init__(self,class_file):
        names,lons,lats = self.init(class_file)
        self.mask = self._make_selection()
        self.sources = np.array([CatalogSource(self,name,SkyDir(lon,lat,self.coords))
                        for name,lon,lat in zip(names,lons,lats)])[self.mask]

    def __iter__(self):
        return self.sources

    def _get_class_module(self,class_file):
        """Import and return module class_file"""
        if not os.path.dirname(class_file) in sys.path:
            sys.path.insert(0,os.path.dirname(class_file))
        cls = os.path.basename(class_file).split('.')[0]
        exec('import %s'%cls)
        return eval(cls)

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
            if card.key[:5]=='TBUCD' and card.value == 'ID_MAIN':
                name_key = cards['TTYPE'+card.key[5:8]].value
                break
            #Sometimes UCDs are declared in comments
            #May be fragile - depends on specific format for comments as in gamma-egr catalog
            value_comment = card._getValueCommentString().split('/')
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
            if card.key[:5]=='TTYPE' and card.value.upper() in ['NAME','ID']:
                name_key = card.value
        try:
            return self.hdu.data.field(name_key)
        except KeyError:
            return
        names = self.hdu.data.field(name_key)

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
                lat_key = ttypes[ttypes.keys()[ttypes.values().index('DE')]].value
        if not lon_key:
            self.coords = SkyDir.GALACTIC
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
        """Make selections specified in class module. Not really implemented yet."""
        selections = [x for x in self.class_module.selection if (x != '' and 'ANGSEP' not in x)]
        dat = self.hdu.data
        mask = np.array([True]*len(dat))
        if not selections: return mask
        catid_pattern = re.compile('@%s_([A-Za-z0-9]+)'%self.class_module.catid.upper())
        defnull_pattern = re.compile('DEFNULL\(([A-Z0-9a-z]+),([0-9e\.]+)\)')
        for sel in selections:
            field_names = catid_pattern.findall(sel)
            fields = {}
            for f in field_names:
                field = dat.field(f)
                sel = catid_pattern.sub(f,sel,1)
                defnull_pattern = re.compile('DEFNULL\(%s,([0-9e\.]+)\)'%f)
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


    def select_circle(self,position,radius):
        """Return an array of CatalogSources within radius degrees of position.

        Arguments:
            position    : SkyDir for center of selection region.
            radius      : radius of selection region.
        """

        ras = np.array([source.skydir.ra() for source in self.sources])
        decs = np.array([source.skydir.dec() for source in self.sources])
        tmask = trap_mask(ras,decs,position,radius)
        ras,decs = ras[tmask],decs[tmask]
        rmask = rad_mask(ras,decs,position,radius)[0]
        return self.sources[tmask][rmask]

    def local_density(self,position,radius=4):
        """Return the local density of catalog sources in a radius-degree region about position."""

        n_sources = len(self.select_circle(position,radius))
        #If no sources within radius, set n_sources = 1 to give lower limit on density
        #Maybe better to expand the radius in this case?
        if n_sources < 1 : n_sources = 1
        solid_angle = deg2rad(radius)**2*math.pi
        return n_sources/solid_angle

class GammaCatalog(Catalog):
    """A catalog of gamma-ray sources (i.e. sources with error circles comparable to LAT)"""

    def __init__(self,class_file):
        names,lons,lats = self.init(class_file)
        errors = self.get_position_errors()
        self.mask = self._make_selection()
        self.sources = np.array([GammaRaySource(self,name,SkyDir(lon,lat,self.coords),error)
                        for name,lon,lat,error in zip(names,lons,lats,errors)])[self.mask]

    def get_position_errors(self):
        q = [x for x in self.class_module.new_quantity if
                    (('LAT' not in x) and ('ECOM' not in x))][0]
        patt = re.compile('@%s_([A-Za-z0-9_]+)'%self.class_module.catid.upper())
        match = patt.search(q)
        if match:
            error_field = match.groups()[0]
            return self.hdu.data.field(error_field)
        else:
            raise CatalogException(self.cat_file,'Could not find position uncertainties.')

class ExtendedCatalog(Catalog):
    """A catalog of extended sources"""

    def __init__(self,class_file):
        names,lons,lats= self.init(class_file)
        radii = self.get_radii()
        self.mask = self._make_selection()
        self.sources = np.array([ExtendedSource(self,name,SkyDir(lon,lat,self.coords),radius)
                                 for name,lon,lat,radius in zip(names,lons,lats,radii)])[self.mask]

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

    def __str__(self):
        return '\t'.join([self.catalog.cat_name,self.name,str(self.skydir.ra()),str(self.skydir.dec())])

    def angular_separation(self,other_skydir):
        """Calculate angular separation between this source and other_skydir"""
        return self.skydir.difference(other_skydir)

    def position_angle(self,other_skydir):
        """Calculate other_skydir angle _from_ other_skydir _to_ this source."""
        ra1,dec1,ra2,dec2 = deg2rad([other_skydir.ra(),other_skydir.dec(),
                            self.skydir.ra(),self.skydir.dec()])
        denom = math.sin(dec1)*math.cos(ra1-ra2) - math.cos(dec1)*math.tan(dec2)
        return math.atan2(math.sin(ra1-ra2), denom)

    def delta_logl(self,position,error_ellipse):
        """Compute Delta(logl) for this association with source at position, with error_ellipse

        Arguments:
            position : SkyDir representing position of other source
            error_ellipse : 3-tuple representing error ellipse of other source in degrees"""

        error_ellipse = deg2rad(error_ellipse)
        phi = self.position_angle(position)-error_ellipse[2]
        angsep = self.angular_separation(position)
        return .5*angsep**2*((math.cos(phi)/error_ellipse[0])**2 +
                               (math.sin(phi)/error_ellipse[1])**2)

    def positional_likelihood(self,position,error_ellipse):
        """Likelihood for point source at position with error_ellipse to be associated with this source.

        Arguments:
            position : SkyDir representing position of other source
            error_ellipse : 3-tuple representing error ellipse of other source in degrees"""

        norm = 2.*math.pi*mul(*deg2rad(error_ellipse)[:2])
        return math.exp(-self.delta_logl(position,error_ellipse))/norm

    def chance_probability(self,position,radius = 4):
        """Probability of chance association"""
        return self.catalog.local_density(position,radius=radius)

    def prior_probability(self):
        """Prior probability for association"""
        return self.catalog.prior

    def posterior_probability(self,position,error_ellipse):
        """Posterior probability for association with given position and error ellipse"""

        arg = (self.chance_probability(position)*(1-self.prior_probability())/
               (self.positional_likelihood(position,error_ellipse)*self.prior_probability()))
        return 1./(1.+arg)

class GammaRaySource(CatalogSource):
    """A source from a gamma-ray catalog."""
    
    def __init__(self,catalog,name,skydir,error):
        super(GammaRaySource,self).__init__(catalog,name,skydir)
        self.error = error

    def posterior_probability(self,position,error_ellipse):
        if self.combined_error(error_ellipse) >= self.angular_separation(position):
            return 1.0
        else:
            return 0.0

    def combined_error(self,error_ellipse):
        return (error_ellipse[0]**2 + self.error**2)**.5

class ExtendedSource(CatalogSource):
    """An extended catalog source"""

    def __init__(self,catalog,name,skydir,radius):
        super(ExtendedSource,self).__init__(catalog,name,skydir)
        self.radius = radius

    def posterior_probability(self,position,error_ellipse):
        if self.angular_separation(position)<=(error_ellipse[0]+self.radius):
            return 1.0
        else:
            return 0.0

class CatalogError(Exception):
    """Exception class for problems with a catalog."""
    def __init__(self,catalog,message):
        self.catalog = catalog
        self.message = message
    def __str__(self):
        return 'In catalog %s:\n\t%s'%(self.catalog,self.message)

def deg2rad(angles):
    try:
        return [angle*math.pi/180. for angle in angles]
    except TypeError:
        return angles*math.pi/180.

def rad2deg(angles):
    try:
        return [angle*180./math.pi for angle in angles]
    except TypeError:
        return angles*180./math.pi

if __name__=='__main__':
    assoc = SourceAssociation('/home/eric/research/catalog/srcid/cat')
    #3C 454.3
    pos, error = SkyDir(343.495,16.149), .016/2.45*1.51
    associations = assoc.id(pos,error,['agn','bzcat','cgrabs','crates'])
    for cat,ass in associations.items():
        print 'Associations in %s:'%cat
        for name,data in ass.items():
            print '\t'.join([name,'prob:',str(data[0])])
    #print('\n'.join([str(x[1]) for x in assoc.id(pos,error,'obj-blazar-crates',.33,.8)]))
    #Couldn't find elliptical errors, but want to test input for error.
    #error = (error,error,0.0)
    #print('\n'.join([str(x[1]) for x in assoc.id(pos,error,'obj-blazar-crates',.33,.8)]))
    #Test for correct failure for wrong length list
    #error = [.5]*4
    #print(assoc.id(pos,error,'obj-agn',.039,.9735))
