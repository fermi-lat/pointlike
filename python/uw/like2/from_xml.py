"""
$Header$
"""
import numpy as np
from uw.utilities import  xml_parsers
from . import (roimodel, sources, diffuse, extended)

class ROImodelFromXML(roimodel.ROImodel):
    
    def load_sources(self, xmlfile):
        self.filename = xmlfile
        handler = xml_parsers.parse_sourcelib(xmlfile)
        self.source_library = handler.outerElements[0]
        # get the HEALPix index from the source_library element
        self.index = int(self.source_library['index'])
        
        if self.ecat is None: #speed up if already loaded
            self.ecat = extended.ExtendedCatalog(self.config.extended, quiet=self.quiet)
        
        xtm = xml_parsers.XML_to_Model()
        for src in handler.sources:
            name, stype = str(src['name']), src['type']
            
            # get the gtlike-style model constructed by Josn's code and create a default version
            spec = src.getChild('spectrum')
            mx  = xtm.get_model(spec,name)
            p, free = mx.get_all_parameters(), mx.free
            classname = mx.__class__.__name__
            if classname=='FrontBackConstant': 
                # needs special attention, not translated properly
                # also the only one that does not have e0 property
                model = mx.__class__( f=p[0], b=p[1])
                model.free = np.array([k['free']==1 for k in spec.children], bool)
            elif classname=='LogParabola' or not hasattr(mx, 'e0'):
                model = mx.__class__(p=p, free=free)
            else: # could be PowerLaw, need to convert I think
                model = mx.__class__(p=p, free=free, e0=mx.e0)
            model.set_external_cov_matrix( mx.get_cov_matrix() )
            sources.set_default_bounds(model)
            
            if stype == 'ExtendedSource': 
                sm = src.getChild('spatialModel')
                # could check that this is consistent with the extended cat
                esrc = self.ecat.lookup(name)
                esrc.model = model
                self.append( esrc )
            elif stype=='PointSource':    
                sd = xml_parsers.get_skydir(src.getChild('spatialModel'))
                self.append(
                    sources.PointSource(name=name, skydir=sd, model=model, 
                        associate=src.get('associate', None),
                        ts = src.get('ts', None)
                    )
                )
            elif stype=='GlobalSource':
                model.background = True # to suppress flux integrals on __str__
                self.append(
                    sources.GlobalSource(name=name, skydir=None, model=model,
                        dmodel = diffuse.diffuse_factory(self.config.diffuse[name], 
                            self.config.event_type_names,),
                    )
                )
            else:
                raise Exception('unrecognized type %s for source %s' % (stype, name))
                
    def __repr__(self):
        return '%s.%s: file %s' % (self.__module__, self.__class__.__name__, self.filename)
    