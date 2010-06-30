"""Class for parsing and writing gtlike-style source libraries.
   Barebones implementation; add additional capabilities as users need.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pointspec_helpers.py,v 1.8 2010/05/24 08:10:30 lande Exp $

   author: Matthew Kerr
"""

import numpy as N
import xml.sax as x
from collections import deque
from uw.like.pointspec_helpers import PointSource,get_diffuse_source
from uw.like.roi_diffuse import DiffuseSource
from uw.like.Models import *
from skymaps import SkyDir,DiffuseFunction,IsotropicSpectrum,IsotropicPowerLaw
import os

JAC = N.log10(N.exp(1))

def decorate(instring,tablevel=1):
	return ''.join(['\t'*tablevel,instring,'\n'])

class Stack(deque):

    def push(self,element):
        self.append(element)

    def pop(self):
        if len(self) == 0: return
        return super(Stack,self).pop()

    def peek(self):
        if len(self) == 0: return
        return self[-1]

class XMLElement(dict):

    def __init__(self,name,*args,**kwargs):
        self.name = name
        self.children = deque()
        super(XMLElement,self).__init__(*args,**kwargs)

    def addChild(self,child):
        self.children.append(child)

    def getChild(self,name):
        for child in self.children:
            if child.name == name: return child

    def __str__(self):
        n  = 'Element Name: %s'%(self.name)
        a1 = 'Attributes:\n\tKey:\t\tValue:\n\t---\t\t-----'
        astrings = '\n'.join(['\t%s\t\t%s'%(k,v) for k,v in zip(self.keys(),self.values())])
        c  = 'Child Element Names:'
        cstrings = '\n'.join([x.name for x in self.children])
        return '\n'.join([n,a1,astrings,c,cstrings])

class SourceHandler(x.ContentHandler,Stack):
    
    def __init__(self):
        self.outerElements = deque()
        self.sources       = deque()

    def __call__(self): return self.lastOff

    def startElement(self,name,attrs):
        self.push(XMLElement(name,attrs))

    def endElement(self,name):
        t = self.pop()
        l = self.peek()
        if l is not None:
            l.addChild(t)
        else:
            self.outerElements.append(t)
        if t.name == 'source':
            self.sources.append(t)

class XML_to_Model(object):
    """Map a parsed XML model onto the Model classes.
    
       This class should be extended as more use cases arise."""

    def __init__(self):

        self.modict   = dict({
            'PowerLaw'  : 'PowerLaw',
            'PowerLaw2' : 'PowerLawFlux',
            'PLSuperExpCutoff' : 'PLSuperExpCutoff',
            'Constant'  : 'Constant',
            'FileFunction' : 'Constant' # a big klugey
            })
        
        self.specdict = dict({
            'PowerLaw'  : ['Prefactor','Index'],
            'PowerLaw2' : ['Integral','Index'],
            'PLSuperExpCutoff' : ['Prefactor','Index1','Cutoff','Index2'],
            'ConstantValue' : ['Value'],
            'FileFunction'  : ['Normalization']
            })

        self.kwargdict = dict({
            'PowerLaw'  : [ ['Scale','e0' ] ],
            'PowerLaw2' : [ ['LowerLimit','emin'], ['UpperLimit','emax'] ],
            'PLSuperExpCutoff' : [ ['Scale','e0' ] ],
            'ConstantValue' : [],
            'FileFunction' : []
            })

    def get_model(self,xml_dict,index_offset=0):

        specname = xml_dict['type']
        params   = xml_dict.children

        d = dict()
        for p in params:
            d[p['name']] = p

        model = eval('%s()'%(self.modict[specname]))

        for ip,p in enumerate(self.specdict[specname]):
            pdict = d[p]
            scale = float(pdict['scale'])
            value = float(pdict['value'])
            if (p == 'Index') or (p == 'Index1'):
                # gtlike uses a neg. index internally so scale > 0
                # means we need to take the negative of the value
                if scale > 0: value = -value
                else:         scale = -scale
                value += index_offset
                model.index_offset = index_offset
            model.p[ip] = value*scale
            model.free[ip] = (pdict['free'] == '1')
            if 'error' in pdict.keys():
                err = float(pdict['error'])
                model.cov_matrix[ip,ip] = (err/value*JAC)**2

        model.p = N.log10(model.p)

        for p in self.kwargdict[specname]:
            pdict = d[p[0]]
            model.__dict__[p[1]] = float(pdict['value'])

        return model


class Model_to_XML(object):
    """Convert a Model instance to an XML entry.
       
       This class should be extended as more use cases arise."""
    
    def __init__(self,debug=False):
        self.x2m = XML_to_Model()
        self.debug = debug

    def update(self,name,scaling=False):
        self.name = name

        if name == 'PowerLaw2':
            self.pname  = ['Integral','Index','LowerLimit','UpperLimit']
            self.pfree  = [1,1,0,0]
            self.pscale = [1e-7,-1,1,1]
            self.pmin   = [1e-4,0,30,30]
            self.pmax   = [1e4,5,5e5,5e5]
            self.pval   = [2,1,1e3,1e5]
            self.perr   = [0,0,-1,-1]

        elif name == 'PowerLaw':
            self.pname  = ['Prefactor','Index','Scale']
            self.pfree  = [1,1,0]
            self.perr   = [0,0,-1]
            if scaling:
                self.pscale = [1,-1,1] ; self.pmin  = [0.1,-1,30]
                self.pmax  = [10,1,5e5]; self.pval = [1,0,1e3]
            else:
                self.pscale = [1e-9,-1,1]
                self.pmin   = [1e-5,0,30]
                self.pmax   = [1e4,5,5e5]
                self.pval   = [1,2,1e3]

        elif name == 'PLSuperExpCutoff':
            self.pname = ['Prefactor','Index1','Cutoff','Index2','Scale']
            self.pfree  = [1,1,1,0,0]
            self.pscale = [1e-9,-1,1000,1,1]
            self.pmin   = [1e-4,0,0.1,0,30]
            self.pmax   = [1e4,5,3e5,5,5e5]
            self.pval   = [1,2,1,1,1000]
            self.perr   = [0,0,0,0,-1]

        elif name == 'FileFunction':
            self.pname = ['Normalization']
            self.pfree = [1]
            self.pscale = [1]
            self.pmin = [0.1]
            self.pmax = [10]
            self.pval = [1]
            self.perr = [0]

        else:
            raise Exception,'Unrecognized model %s'%(name)

    def find_param(self,pname):
        for iparam,param in enumerate(self.pname):
            if param == pname: return iparam
        return None

    def param_strings(self):
        err_strings = ['error = "%s"'%(e) if (e>0) else '' for e in self.perr]
        return ['\t<parameter %s free="%d" max="%s" min="%s" name="%s" scale="%s" value="%s"/>'%(e,f,m1,m2,n,s,v)
            for e,f,m1,m2,n,s,v in zip(err_strings,self.pfree,self.pmax,self.pmin,self.pname,self.pscale,self.pval)]

    def getXML(self,tablevel=1,spec_attrs='',comment_string=None):
        cstring = [] if comment_string is None else [comment_string]
        strings = ['<spectrum %s type="%s">'%(spec_attrs,self.name)] + \
                   cstring + self.param_strings() + ['</spectrum>']
        return ''.join([decorate(s,tablevel=tablevel) for s in strings])

    def process_photon_index(self,model,index,val,err):
        # note this appears to be inconsistent in xml
        self.perr[index]   = abs(self.perr[index])
        if (not hasattr(model,'index_offset')) or (model.index_offset == 0):
            self.pval[index] = abs(self.pval[index])
            self.pscale[index] = -abs(self.pscale[index])
        else:
            self.pval[index] = val - model.index_offset
            self.pscale[index] = -1

    def process_model(self,model,scaling=False):
        """model an instance of Model"""

        if model.name == 'ExpCutoff':
            model = convert_exp_cutoff(model)

        # map the Model instance onto an xml-style model
        for xml_name,v in self.x2m.modict.iteritems():
            if v == model.name: break
        if v != model.name:
            raise Exception,'Unable to find an XML model for %s'%(model.name)
        self.update(xml_name,scaling=scaling)

        # replace spectral parameters
        specparams = self.x2m.specdict[xml_name]
        vals,errs  = model.statistical(absolute=True)
        for iparam,param in enumerate(specparams):
            if self.debug: print 'Processing %s'%(param)
            index = self.find_param(param)
            if index is None:
                raise Exception,'Unrecognized parameter %s'%(param)
            self.pfree[index] = model.free[iparam]
            self.pval[index] = vals[iparam] / self.pscale[index]
            self.perr[index] = errs[iparam] / self.pscale[index]
            if (param == 'Index') or (param == 'Index1'):
                self.process_photon_index(model,index,vals[iparam],errs[iparam])
            if self.pval[index] < self.pmin[index]:
                print 'WARNING: Found %s for minimum value %s (%s)'%(str(self.pval[index]),str(self.pmin[index]),param)
                print 'Setting parameter value to minimum.'
                self.pval[index] = self.pmin[index]
            if self.pval[index] > self.pmax[index]:
                print 'Found %s for maximum value %s (%s)'%(str(self.pval[index]),str(self.pmax[index]),param)
                print 'Setting parameter value to maximum.'
                self.pval[index] = self.pmax[index]

        # replace non-variable parameters
        kwargparams = self.x2m.kwargdict[xml_name]
        for xml_key,mod_key in kwargparams:
            if self.debug: print 'Processing %s'%(xml_key)
            index = self.find_param(xml_key)
            if index is None:
                raise Exception,'Unrecognized parameter %s'%(param)
            self.pfree[index] = False
            self.pval[index]  = model.__dict__[mod_key]
            self.pscale[index] = 1
            self.pmin[index] = self.pval[index]
            self.pmax[index] = self.pval[index]
            self.perr[index] = -1

def get_skydir(elem):
    """convert a parsed SpatialModel into a SkyDir."""
    d = dict()
    for p in elem.children:
        d[p['name']] = float(p['value'])
    return SkyDir(d['RA'],d['DEC'])

def makePSSpatialModel(ra,dec,tablevel=1):
    """Encode the spatial model in XML for a point source at ra,dec."""

    strings = [
        '<spatialModel type="SkyDirFunction">',
        '\t<parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="%s"/>' %(ra),
        '\t<parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="%s"/>' %(dec),
        '</spatialModel>'
        ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeDSConstantSpatialModel(tablevel=1):
    """Encode an isotropic model."""
    
    strings = [
        '<spatialModel type="ConstantValue">',
        '\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>',
        '</spatialModel>'
	]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeDSMapcubeSpatialModel(filename='ERROR',tablevel=1):
    """Encode a mapcube model."""
    strings = [
        '<spatialModel file="%s" type="MapCubeFunction">'%(filename),
        '\t<parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>',
	    '</spatialModel>'
    ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def parse_sourcelib(xml):
    parser = x.make_parser()
    handler = SourceHandler()
    parser.setContentHandler(handler)
    parser.parse(xml)
    return handler

def parse_point_sources(handler):
    point_sources = deque()
    xtm = XML_to_Model()
    for src in handler.sources:
        if src['type'] != 'PointSource': continue
        sd = get_skydir(src.getChild('spatialModel'))
        mo = xtm.get_model(src.getChild('spectrum'))
        point_sources.append(PointSource(sd,src['name'],mo,leave_parameters=True))
    return point_sources

def parse_diffuse_sources(handler):
    gds = get_diffuse_source
    ds = deque()
    xtm = XML_to_Model()
    for src in handler.sources:
        if src['type'] != 'DiffuseSource': continue
        spatial  = src.getChild('spatialModel')
        spectral = src.getChild('spectrum')
        if spatial['type'] == 'ConstantValue':
            if spectral['type'] == 'FileFunction':
                fname = str(os.path.expandvars(spectral['file']))
                mo = xtm.get_model(spectral)
                ds.append(gds('ConstantValue',None,mo,fname))
            elif (spectral['type'] == 'PowerLaw' ) or (spectral['type'] == 'PowerLaw2'):
                mo = xtm.get_model(spectral)
                ds.append(gds('ConstantValue',None,mo,None))
            else:
                raise Exception,'Isotropic model not implemented'
    
        elif spatial['type'] == 'MapCubeFunction':
            fname = str(os.path.expandvars(spatial['file']))
            if spectral['type'] == 'ConstantValue':
                mo = xtm.get_model(spectral)
            elif spectral['type'] == 'PowerLaw':
                mo = xtm.get_model(spectral,index_offset=1)
                print mo.index_offset
            else:
                raise Exception,'Non-isotropic model not implemented'
            ds.append(gds('MapCubeFunction',fname,mo,None))
            
        else:
            raise Exception,'Spatial model not recognized'
        ds[-1].name = src['name']
    return ds

def parse_sources(xmlfile):
    """Convenience function to parse an entire file into
       point sources and diffuse sources."""
    handler = parse_sourcelib(xmlfile)
    ps = parse_point_sources(handler)
    ds = parse_diffuse_sources(handler)
    return ps,ds

def unparse_point_sources(point_sources):
    """Convert a list (or other iterable) of PointSource objects into XML."""
    xml_blurbs = Stack()
    m2x = Model_to_XML()
    for ps in point_sources:
        skyxml = makePSSpatialModel(ps.skydir.ra(),ps.skydir.dec())
        m2x.process_model(ps.model)
        specxml = m2x.getXML()
        s1 = '\n<source name="%s" type="PointSource">\n'%(ps.name)
        s2 = '</source>'
        xml_blurbs.push(''.join([s1,specxml,skyxml,s2]))
    return xml_blurbs

def process_diffuse_source(ds):
    """Convert an instance of DiffuseSource into an XML blurb."""
    m2x = Model_to_XML()
    dm = ds.dmodel
    if hasattr(dm,'__len__'):  dm = dm[0]

    if isinstance(dm,DiffuseFunction):
        filename = os.path.abspath(dm.name())
        skyxml = makeDSMapcubeSpatialModel(filename=filename)
        m2x.process_model(ds.smodel,scaling=True)
        specxml = m2x.getXML()
    else:
        skyxml = makeDSConstantSpatialModel()
        if isinstance(dm,IsotropicSpectrum):
            filename = os.path.abspath(dm.name())
            m2x.process_model(ds.smodel,scaling=True)
            specxml = m2x.getXML(spec_attrs='file=\"%s\"'%(filename))
        elif isinstance(dm,IsotropicPowerLaw):
            sd = SkyDir()
            iflux = dm.integral(sd,100,3e5)
            spind = dm.value(sd,100)*100/iflux + 1
            pl = PowerLawFlux()
            pl.p[0] = ds.smodel.p[0] + N.log10(iflux)
            p1.p[1] = N.log10(spind+10**(ds.smodel.p[1]-ds.smodel.index_offset))
            pl.cov_matrix = ds.smodel.cov_matrix.copy() #correct?
            m2x.process_model(pl,scaling=False)
            specxml = m2x.getXML()
        else:
            raise Exception,'Did not recognize %s'%(ds.name)
    s1 = '<source name="%s" type="DiffuseSource">\n'%(ds.name)
    s2 = '</source>'
    return ''.join([s1,specxml,skyxml,s2])
    
def unparse_diffuse_sources(diffuse_sources):
    """Convert a list of DiffuseSources into XML blurbs."""
    xml_blurbs = Stack()
    for ds in diffuse_sources:
        xml_blurbs.push(process_diffuse_source(ds))
    return xml_blurbs

def writeXML(stacks,filename):
    """Write XML blurbs to a gtlike-style XML file."""
    f = open(filename,'wb')
    f.write('<source_library title="source_library">\n')
    for stack in stacks:
        for elem in stack:
            f.write(elem)
    f.write('</source_library>')

def writeROI(roi,filename):
    """Write out the contents of an ROIAnalysis source model
       to a gtlike XML file."""
    dxml = unparse_diffuse_sources(roi.dsm.diffuse_sources)
    pxml = unparse_point_sources(roi.psm.point_sources)
    writeXML([dxml,pxml],filename)
            
