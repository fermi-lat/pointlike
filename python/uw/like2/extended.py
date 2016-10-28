"""
Extended source code
Much of this adapts and utilizes 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/extended.py,v 1.14 2016/05/31 17:29:50 burnett Exp $

"""
import os, copy, glob
import numpy as np
from astropy.io import fits as pyfits
from skymaps import SkyDir
from uw.like import  Models
from uw.like.SpatialModels import Disk,EllipticalDisk,Gaussian,EllipticalGaussian,SpatialMap

from . import sources, response


class ExtendedSourceCatalog(object):
    """ Object to read in spatially extended sources from Elizabeth Ferrara suggested
        table of extended sources.

        The input to this object is the directory contaning the extended source table
        and a folder for xml descriptions and templates of extended sources.

        For a description of the format of these sources, see:

            https://confluence.slac.stanford.edu/x/Qw2JBQ
        Copied from like.roi_catalogs.py since no longer supported there

    '"""

    def __init__(self,archive_directory, force_map=False):
        self.force_map=force_map
        self.__open_catalog__(archive_directory)

    def __open_catalog__(self,archive_directory):
        """ Parses the LAT_extended_sources.fit table 
            to get a list of the extended sources. """
        self.archive_directory = archive_directory

        filename=os.path.join(self.archive_directory,"LAT_extended_sources*.fit*")
        filename=glob.glob(filename)
        if len(filename)!=1: raise Exception("Unable to find LAT_extended_sources.fit archive file.")
        filename=filename[0]
        f = pyfits.open(filename)
        self.names = np.array([n.strip() for n in f[1].data.field('Source_Name')])
        ras   = f[1].data.field('RAJ2000')
        decs  = f[1].data.field('DEJ2000')
        form = [x.strip() for x in f[1].data.field('Model_Form')]
        major = f[1].data.field('Model_SemiMajor')
        minor = f[1].data.field('Model_SemiMinor')
        posang = f[1].data.field('Model_PosAng')

        # The xml filename for the extended sources.
        self.xmls      = f[1].data.field('Spectral_Filename').astype(str)
        self.templates = f[1].data.field('Spatial_Filename').astype(str)

        self.dirs    = map(SkyDir,np.asarray(ras).astype(float),np.asarray(decs).astype(float))

        if self.archive_directory in [ "Extended_archive_v01", "Extended_archive_v02", 
                                       "Extended_archive_jbb02", "Extended_archive_jbb03" ]:

            # old style archives
            os.environ["LATEXTDIR"]=os.path.join(self.archive_directory,'Templates')
        else:
            os.environ["LATEXTDIR"]=self.archive_directory


        # build up a list of the analytic extended source shapes (when applicable)
        self.spatial_models = []
        fail = False
        for i in range(len(self.names)):
            #print 'unpacking source "{}", form: "{}"'.format(self.names[i], form[i]), 
            if self.force_map or form[i] in ('Map','Ring'):
                self.spatial_models.append(
                    SpatialMap(file=self.templates[i].replace(' ', '')))
            elif form[i] in ('Disk', 'RadialDisk'):
                if major[i] == minor[i] and posang[i] == 0:
                    self.spatial_models.append(Disk(p=[major[i]],center=self.dirs[i]))
                else:
                    self.spatial_models.append(EllipticalDisk(p=[major[i],minor[i],posang[i]],center=self.dirs[i]))
            elif form[i] in ('2D Gaussian', 'RadialGaussian'):
                if major[i] == minor[i] and posang[i] == 0:
                    self.spatial_models.append(Gaussian(p=[major[i]/Gaussian.x68],center=self.dirs[i]))
                else:
                    self.spatial_models.append(
                        EllipticalGaussian(p=[major[i]/Gaussian.x68,minor[i]/Gaussian.x68,posang[i]],
                                           center=self.dirs[i]))
            else:
                print 'Unrecognized spatial model: {}, assuming a Map'.format(form[i])
                try:
                    self.spatial_models.append(
                        SpatialMap(file=self.templates[i])
                    )
                except Exception, msg:
                    print 'Failure: {}'.format(msg)
                    fail=True
                    
            assert not fail, 'Quit due to errors'
            #remember the fits file template in case the XML needs to be saved out.
            self.spatial_models[-1].original_template = self.templates[i]
            self.spatial_models[-1].original_parameters = self.spatial_models[-1].p.copy()


        self.spatial_models = np.asarray(self.spatial_models)

    def get_sources(self,skydir,radius=15):
        """ Returns a list of ExtendedSource objects from the extended source
            catalog that have a center withing a distance radius of the
            position skydir. 
           
        Note that if there are none, it returns an empty list.   
        """

        from uw.utilities.xml_parsers import parse_sources

        diffs    = np.degrees(np.asarray([skydir.difference(d) for d in self.dirs]))
        mask     = diffs < radius
        if sum(mask)==0: return []
        diffs    = diffs[mask]
        sorting = np.argsort(diffs)

        names     = self.names[mask][sorting]
        xmls      = self.xmls[mask][sorting]
        spatials  = self.spatial_models[mask][sorting]

        sources = []
        for name,xml,spatial in zip(names,xmls,spatials):

            full_xml=os.path.join(self.archive_directory,'XML',os.path.basename(xml))

            # Use the built in xml parser to load the extended source.
            ps,ds=parse_sources(xmlfile=full_xml.replace(' ', ''))
            if len(ps) > 0: 
                raise Exception("A point source was found in the extended source file %s" % xmlfile)
            if len(ds) > 1: 
                raise Exception("No diffuse sources were found in the extended soruce file %s" % xmlfile)
            if len(ds) < 1: 
                raise Exception("More than one diffuse source was found in the extended source file %s" % xmlfile)

            source=ds[0]
            if False: # don't support analytic maps  spatial is not None:
                # replace the SpatialMap extended source with an analytic one.
                analytic_source = ExtendedSource(name=source.name,model=source.model,
                                                 spatial_model=spatial)
                analytic_source.original_template = source.spatial_model.file # for reference
                sources.append(analytic_source)
            else:
                sources.append(source)

        return sources

    def merge_lists(self,skydir,radius=15,user_point_list=[],user_diffuse_list=[]):
        """ Get a list of extended sources from the  catalog and merge it with 
            and already existin glist of point and diffuse sources.

            Unlike the FermiCatalog, no source pruning is done. """

        cat_list = self.get_sources(skydir,radius)

        return user_point_list,user_diffuse_list+cat_list


class ExtendedCatalog( ExtendedSourceCatalog):
    """ subclass to add this lookup function 
    TODO: merge, keeping only needed features
    """

    def __init__(self, extended_catalog_name, force_map=False,  **kwargs):
        """ initialize by also filling an array with all source spectral models"""
        self.alias = kwargs.pop('alias', dict())
        self.quiet = kwargs.pop('quiet', True)
        self.catname = extended_catalog_name
        if os.path.isabs(self.catname):
            extended_catalog_name = self.catname
        else:
            extended_catalog_name = \
                os.path.expandvars(os.path.join('$FERMI','catalog',extended_catalog_name))
        if not os.path.exists(extended_catalog_name):
            raise Exception('extended source folder "%s" not found' % extended_catalog_name)
        if not self.quiet:
            print 'Loaded extended catalog %s' % extended_catalog_name
        
        super(ExtendedCatalog,self).__init__(extended_catalog_name, force_map=force_map)#**kwargs)
        
        # create list of sources using superclass, for lookup by name
        self.sources = [self.get_sources(self.dirs[i], 0.1)[0] for i in range(len(self.names))]
        for source, spatial_model in zip(self.sources, self.spatial_models):
            model = source.model
            if model.mappers[0].__class__.__name__== 'LimitMapper':
                #print 'converting mappers for model for source %s, model %s' % (source.name, model.name)
                source.model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))

            # Replace the spatial model with the one specified by the FITS catalog list    
            dmodel = source.dmodel
            if spatial_model.name != 'SpatialMap':
                assert spatial_model.center.difference(dmodel[0].center)<0.001,\
                     'Center discrepance for {}: {},{}'.format(source.name, spatial_model.center, dmodel[0].center)
                source.dmodel[0]=source.spatial_model=spatial_model

    def __repr__(self):
        return '%s.%s: %s' % (self.__module__, self.__class__.__name__, self.catname)
    def realname(self, cname):
        """ cname was truncated"""
        if cname in self.names: return cname
        for name in self.names:
            assert name is not None, 'bad name'
            t = name.replace(' ','')
            if t==cname: return name
        assert 'compressed name %s not found in list of names, %s' %(cname,self.names)

    def lookup(self, name):
        """ return an ExtendedSource object by name, None if not found """
        aname = self.alias.get(name,name) #alias will be the new name
        try:
            i = list(self.names).index(name)
        except ValueError:
            return None
        source = self.sources[i]
        # make a new object copied from original
        if source.model.name=='BrokenPowerLaw': #convert this
            model = Models.LogParabola()
        else: model = source.model
        if model.name=='LogParabola': 
           model.free[-1]=False # E_break ne free

        ### seems to be necessary for some models created from 
        if model.mappers[0].__class__.__name__== 'LimitMapper':
            print 'wrong mappers: converting model for source %s, model %s' % (name, model.name)
            model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))
        extsource = sources.ExtendedSource(name=self.realname(aname), 
            skydir=source.skydir,
            model = model, 
            dmodel= source.spatial_model
            )
        if extsource.model.name=='LogParabola': 
            assert sum( extsource.model.free)<4, 'E_break is free?'
        return extsource  
    def __getitem__(self, name): return self.lookup(name)

def make_pictures(config_dir='.', image_folder = 'extended_images'):
    """ make a folder with images for all extended sources """
    from uw.like2 import configuration
    from skymaps import Band
    if not os.path.exists(image_folder):
        os.mkdir(image_folder)
    cf = configuration.Configuration(config_dir, quiet=True, postpone=True)
    class Bandlite(object):
        def __init__(self, roi_dir, event_type=1, energy =133.352):
            self.event_type=event_type
            self.energy=energy
            self.sd=roi_dir
            self.psf=cf.psfman(event_type,energy)
            self.exposure = cf.exposureman(event_type,energy)
            self.radius =5

    ecat = ExtendedCatalog(cf.extended)
    for name in ecat.names:
        print 'processing %s ...' % name ,
        source = ecat.lookup(name)
        b12 = Band(12); roi_index = b12.index(source.skydir); 
        roi_dir = b12.dir(roi_index)
        conv = convolution.ExtendedConvolver(source, Bandlite(roi_dir=roi_dir))
        conv.create_grid()
        fig = conv.show_source()
        fig.savefig('%s/%s_map.png' %(image_folder,name.replace(' ','_')))
        print 



