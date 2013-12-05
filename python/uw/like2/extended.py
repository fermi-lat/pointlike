"""
Extended source code
Much of this adapts and utilizes 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/extended.py,v 1.7 2013/12/04 05:21:21 burnett Exp $

"""
import os, copy
import numpy as np
from uw.like import roi_catalogs, Models
from . import sources, response


class ExtendedCatalog( roi_catalogs.ExtendedSourceCatalog):
    """ subclass to add this lookup function """

    def __init__(self, extended_catalog_name, **kwargs):
        """ initialize by also filling an array with all source spectral models"""
        self.alias = kwargs.pop('alias', dict())
        self.quiet = kwargs.pop('quiet', True)
        self.catname = extended_catalog_name
        extended_catalog_name = \
            os.path.expandvars(os.path.join('$FERMI','catalog',extended_catalog_name))
        if not os.path.exists(extended_catalog_name):
            raise Exception('extended source folder "%s" not found' % extended_catalog_name)
        if not self.quiet:
            print 'Loaded extended catalog %s' % extended_catalog_name
        
        super(ExtendedCatalog,self).__init__(extended_catalog_name, force_map=True)#**kwargs)
        
        # create list of sources using superclass, for lookup by name
        self.sources = [self.get_sources(self.dirs[i], 0.1)[0] for i in range(len(self.names))]
        for source in self.sources:
            model = source.model
            if model.mappers[0].__class__.__name__== 'LimitMapper':
                #print 'converting mappers for model for source %s, model %s' % (source.name, model.name)
                source.model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))

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
        if model.name=='LogParabola': model.free[-1]=False # E_break ne free
        ### seems to be necessary for some models created from 
        if model.mappers[0].__class__.__name__== 'LimitMapper':
            print 'wrong mappers: converting model for source %s, model %s' % (name, model.name)
            model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))
        extsource= sources.ExtendedSource(name=self.realname(aname), 
            skydir=source.skydir,
            model = model, 
            dmodel= source.spatial_model
            )
        if extsource.model.name=='LogParabola': extsource.free[-1]=False # E_break not free
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



