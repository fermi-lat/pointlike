"""
Extended source code
Much of this adapts and utilize3s
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/extended.py,v 1.1 2013/10/29 03:23:52 burnett Exp $

"""
import os, copy
import numpy as np
from uw.like import pointspec_helpers, Models
from . import sources, convolution


class ExtendedSource(sources.Source):

    def __str__(self):
        return self.name + ' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')  
  
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)
        
    def copy(self):
        """ return a new ExtendSource object, with a copy of the model object"""
        ret = ExtendedSource(**self.__dict__)
        ret.model = self.model.copy()
        if ret.model.name=='LogParabola':
            ret.model.free[-1]=False # make sure Ebreak is frozen
        return ret
 
    def grid_generator(self,  band, pixelsize):
        """convolvable grid for this object
            psf : PSF (from configuration) for given event_type
            exposure : Exposure (from configuration) for event type        """
        grid = convolution.ExtendedGridGenerator(band, self.dmodel, pixelsize=pixelsize)
        return grid
                
    def grid_for_band(self, band, pixelsize=0.1):
        """return convolved grid for a Band object, which has the event_type, position, psf, esposure, and energy
        """
        assert hasattr(band, 'sd'), 'Not a Band-like object? %s' %band
        grid = self.grid_generator( band, pixelsize=pixelsize)
        return grid(band.energy)


class ExtendedCatalog( pointspec_helpers.ExtendedSourceCatalog):
    """ subclass to add this lookup function """

    def __init__(self, extended_catalog_name, **kwargs):
        """ initialize by also filling an array with all source spectral models"""
        self.alias = kwargs.pop('alias', dict())
        extended_catalog_name = \
            os.path.expandvars(os.path.join('$FERMI','catalog',extended_catalog_name))
        if not os.path.exists(extended_catalog_name):
            raise Exception('extended source folder "%s" not found' % extended_catalog_name)
        print 'Loaded extended catalog %s' % extended_catalog_name
        
        super(ExtendedCatalog,self).__init__(extended_catalog_name, **kwargs)
        
        # create list of sources using superclass, for lookup by name
        self.sources = [self.get_sources(self.dirs[i], 0.1)[0] for i in range(len(self.names))]
        for source in self.sources:
            model = source.model
            if model.mappers[0].__class__.__name__== 'LimitMapper':
                #print 'converting mappers for model for source %s, model %s' % (source.name, model.name)
                source.model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))

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
        extsource= ExtendedSource(name=self.realname(aname), 
            skydir=source.skydir,
            model = model, 
            #spatial_model = source.spatial_model,
            #smodel= model,      # these reference copies needed
            dmodel= source.spatial_model
            )
        if extsource.model.name=='LogParabola': extsource.free[-1]=False # E_break not free
        return extsource  

def test(config_dir='.', source_name='LMC', sdir=(78.032,-70.216), event_type=1, energy=133 ):
    """Test and demonstrate interface
    """
    from uw.like2 import configuration
    import skymaps
    cf = configuration.Configuration(config_dir, quiet=True, postpone=True)
    ecat = ExtendedCatalog(cf.extended)
    lmc =ecat.lookup(source_name)
    roi_dir = skymaps.SkyDir(*sdir) #ROI that contains LMC
    class Bandlite(object):
        def __init__(self, roi_dir=roi_dir, event_type=event_type, energy =energy):
            self.event_type=event_type
            self.energy=energy
            self.sd=roi_dir
            self.psf=cf.psfman(event_type,energy)
            self.exposure = cf.exposureman(event_type,energy)
            self.radius = 5
    band = Bandlite()
    print 'Testing source "%s" with band parameters\n' %lmc
    for item in band.__dict__.items():
        print '\t%-10s %s' % item
    grid = lmc.grid_for_band(band, pixelsize=0.2)
    print 'overlap: %.3f,  exposure_ratio: %.3f' %( grid.overlap(),grid.exposure_ratio())
    return grid.show() # return a figure showing convolution
    