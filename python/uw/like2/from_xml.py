"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/from_xml.py,v 1.5 2014/04/14 17:42:21 burnett Exp $
"""
import os, numpy as np
from uw.utilities import  xml_parsers
from . import (roimodel, sources, diffuse, extended)

class XMLparser(dict):
    """ manage the XML file 
    Assume an enclosing element 'source_library' containing 'source' elements.
    This class is a subclass of dict, containing the attributes of source_library.
    The attribute 'sources' is a list of source elements
    """
    def __init__(self, xmlfile):
        if not os.path.isabs(xmlfile):
            xmlfile = os.path.expandvars(xmlfile)
            
        self.handler = xml_parsers.parse_sourcelib(xmlfile)
        self.update(self.handler.outerElements[0])
        self.sources = self.handler.sources


class ROImodelFromXML(roimodel.ROImodel):
    
    def load_sources(self, roi_spec=None):
        """
        """
        if roi_spec is not None:
            handler = XMLparser(roi_spec)
            self.input_xml=roi_spec
        else:
            handler = self.config.xml_parser
            self.input_xml=self.config.input_xml
        self.index = int(handler.get('index', -1))

        xtm = xml_parsers.XML_to_Model()
        for src in handler.sources:
            name, stype = str(src['name']), src['type']
            
            # get the gtlike-style model constructed by Josh's code and create a default version
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
        return '%s.%s: file %s' % (self.__module__, self.__class__.__name__, self.input_xml)
    
