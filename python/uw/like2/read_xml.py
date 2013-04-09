"""
Copied from skymodel, not tested
"""
import os
from xml import sax
from uw.like2 import skymodel
from uw.utilities import keyword_options

from uw.utilities import xml_parsers

class XMLSkyModel(skymodel.SkyModel):
    """A SkyModel initialized from a stored XML representation."""

    defaults= (
        ('auxcat', None, 'name of auxilliary catalog of point sources to append or names to remove',),
        ('newmodel', None, 'if not None, a string to eval\ndefault new model to apply to appended sources'),
        ('filter',   lambda s: True,   'selection filter: see examples at the end.'), 
        ('global_check', lambda s: None, 'check global sources: can modify parameters'),
        ('closeness_tolerance', 0., 'if>0, check each point source for being too close to another, print warning'),
        ('quiet',  False,  'make quiet' ),
    )
    @keyword_options.decorate(defaults)
    def __init__(self,xml,**kwargs):
        keyword_options.process(self,kwargs)
        self._parse_xml(xml)
        self.nside = self.handler.nside
        self._load_sources()
        self._load_globals()
        #self.nside = int(np.sqrt(len(self.global_sources)/12))
        self.load_auxcat()

    def _parse_xml(self,xml):
        self.parser = sax.make_parser()
        self.handler = SkyModelHandler()
        self.parser.setContentHandler(self.handler)
        self.parser.parse(xml)

    def _parse_global_sources(self):
        pass

    def _load_sources(self):
        self.point_sources = xml_parsers.parse_point_sources(self.handler,SkyDir(0,0),180)
        #parse diffuse sources checks the sources list, so won't grab the globals
        self.extended_sources = xml_parsers.parse_diffuse_sources(self.handler)

    def _load_globals(self):
        gds = pointspec_helpers.get_diffuse_source
        xtm = xml_parsers.XML_to_Model()
        self.global_sources = []
        self.diffuse = []
        self.diffuse_dict = {}
        for roi in self.handler.rois:
            index = int(roi['index'])
            gss = []
            for source in roi.children:
                spatial = source.getChild("spatialModel")
                spectral = source.getChild("spectrum")
                name = str(source['name'])
                if spatial['type'] == 'ConstantValue':
                    if spectral['type'] == 'FileFunction':
                        diffdir,fname = os.path.split(str(os.path.expandvars(spectral['file'])))
                        mo = xtm.get_model(spectral,name)
                        if not self.diffuse_dict.has_key(name):
                            self.diffuse_dict[name] = [gds('ConstantValue',None,mo,fname,name,diffdir=diffdir)]
                            self.diffuse += [fname]
                        gss += [sources.GlobalSource(model=mo,index=index,skydir=None,name=name)]
                    elif (spectral['type'] == 'PowerLaw' ) or (spectral['type'] == 'PowerLaw2'):
                        mo = xtm.get_model(spectral,name)
                        if not self.diffuse_dict.has_key(name):
                            self.diffuse_dict[name] = [gds('ConstantValue',None,mo,None,name)]
                        gss += [sources.GlobalSource(model=mo,index=index,skydir=None,name=name)]
                    elif spectral['type']=='CompositeSpectrum':
                        dss = []
                        fnames = []
                        for i,sp in enumerate(spectral.children):
                            diffdir,fname = os.path.split(str(os.path.expandvars(sp['file'])))
                            dss+=[gds('ConstantValue',None,mo,fname,name,diffdir=diffdir)]
                            fnames += [fname]
                            if i==0:
                                mo = xtm.get_model(sp,name)
                                gss += [sources.GlobalSource(model=mo,index=index,skydir=None,name=name)]
                        if not self.diffuse_dict.has_key(name):
                            self.diffuse_dict[name] = [dss]
                            self.diffuse += [fnames]
                    else:
                        raise Exception,'Isotropic model not implemented'
                elif spatial['type'] == 'MapCubeFunction':
                    diffdir,fname = os.path.split(str(os.path.expandvars(spatial['file'])))
                    if spectral['type'] == 'ConstantValue' or spectral['type'] == 'FrontBackConstant':
                        mo = xtm.get_model(spectral,name)
                        gss += [sources.GlobalSource(model=mo,index=index,skydir=None,name=name)]
                    elif spectral['type'] == 'PowerLaw' or spectral['type'] == 'PowerLaw2':
                        mo = xtm.get_model(spectral,name,index_offset=1)
                        gss += [sources.GlobalSource(model=mo,index=index,skydir=None,name=name)]
                    else:
                        raise Exception('Non-isotropic model "%s" not implemented' % spatial['type'])
                    if not self.diffuse_dict.has_key(name):
                        self.diffuse+=[os.path.split(fname)[1]]
                        self.diffuse_dict[name] = [gds('MapCubeFunction',fname,mo,None,name,diffdir=diffdir)]
                else:
                    raise Exception('Diffuse spatial model "%s" not recognized' % spatial['type'])
            self.global_sources += [gss]


class SkyModelHandler(sax.handler.ContentHandler, xml_parsers.Stack):
    """ContentHandler for parsing the XML representation of a SkyModel
    """
    def __init__(self):
        self.outerElements = collections.deque()
        self.sources       = collections.deque()
        self.rois = collections.deque()

    def __call__(self): return self.lastOff

    def startElement(self,name,attrs):
        self.push(xml_parsers.XMLElement(name,attrs))
        if name=='roi_info':
            self.nside = int(attrs.get('nside',12))

    def endElement(self,name):
        t = self.pop()
        l = self.peek()
        if l is not None:
            l.addChild(t)
            if l.name!='roi' and t.name == 'source':
                self.sources.append(t)
            elif t.name == 'roi':
                self.rois.append(t)
        else:
            self.outerElements.append(t)

