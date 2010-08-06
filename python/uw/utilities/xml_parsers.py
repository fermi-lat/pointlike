"""Class for parsing and writing gtlike-style source libraries.
   Barebones implementation; add additional capabilities as users need.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/xml_parsers.py,v 1.8 2010/08/04 23:27:11 kerrm Exp $

   author: Matthew Kerr
"""

import numpy as N
import xml.sax as x
from collections import deque
from uw.like.pointspec_helpers import PointSource,get_diffuse_source
from uw.like.roi_diffuse import DiffuseSource
from uw.like.roi_extended import ExtendedSource
from uw.like.Models import *
from uw.like.SpatialModels import *
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
            'BrokenPowerLaw' : 'BrokenPowerLaw',
            'BrokenPowerLaw2' : 'BrokenPowerLawFlux',
            'PLSuperExpCutoff' : 'PLSuperExpCutoff',
            'Constant'  : 'Constant', # should this have been ConstantValue key?
            'ConstantValue' : 'Constant',
            'FileFunction' : 'Constant' # a big klugey
            })
        
        self.specdict = dict({
            'PowerLaw'  : ['Prefactor','Index'],
            'PowerLaw2' : ['Integral','Index'],
            'BrokenPowerLaw' : ['Prefactor', 'Index1', 'Index2', 'BreakValue'],
            'BrokenPowerLaw2' : ['Integral', 'Index1', 'Index2', 'BreakValue'],
            'PLSuperExpCutoff' : ['Prefactor','Index1','Cutoff','Index2'],
            'ConstantValue' : ['Value'],
            'FileFunction'  : ['Normalization']
            })

        self.kwargdict = dict({
            'PowerLaw'  : [ ['Scale','e0' ] ],
            'PowerLaw2' : [ ['LowerLimit','emin'], ['UpperLimit','emax'] ],
            'BrokenPowerLaw' : [],
            'BrokenPowerLaw2' : [ ['LowerLimit','emin'], ['UpperLimit','emax'] ],
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
            if (p == 'Index') or (p == 'Index1') or (p == 'Index2'):
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

class XML_to_SpatialModel(object):

    def __init__(self):

        # exclude position from spatialdict
        self.spatialdict = dict({
            'Gaussian'           : ['Sigma'],
            'PseudoGaussian'     : [],
            'Disk'               : ['Sigma'],
            'PseudoDisk'         : [],
            'NFW'                : ['Sigma'],
            'PseudoNFW'          : [],
            'EllipticalGaussian' : ['MajorAxis','MinorAxis','PositionAngle'],
            })

    def get_spatial_model(self,xml_dict):
        """ Kind of like getting a spectral model, but
            instead returns a SpatialModel object. 
            
            Since extended ources don't exist in gtlike, no need
            to keep backwards compatability. Just trying
            to keep a farmilliar interface and make the
            naming convention similar to gtobssim. """

        spatialname = xml_dict['type']
        params      = xml_dict.children

        d = dict()
        for p in params: d[p['name']] = p

        if d.has_key('RA') and d.has_key('DEC') and not (d.has_key('L') or d.has_key('B')):
            coordsystem="SkyDir.EQUATORIAL"
            self.spatialdict[spatialname]=['RA','DEC']+self.spatialdict[spatialname]
        elif d.has_key('L') and d.has_key('B') and not (d.has_key('RA') or d.has_key('DEC')):
            coordsystem="SkyDir.GALACTIC"
            self.spatialdict[spatialname]=['L','B']+self.spatialdict[spatialname]
        else: raise Exception("Unable to parse spatial model %s. Either RA and Dec or L and B must be parameters." % spatialname)

        spatial_model = eval('%s(coordsystem=%s)' % (spatialname,coordsystem))

        for ip,p in enumerate(self.spatialdict[spatialname]):
            pdict = d[p]

            scale = float(pdict['scale'])
            value = float(pdict['value'])

            if spatial_model.log[ip]:
                spatial_model.p[ip] = N.log10(value*scale) 
                if pdict.has_key('min'):
                    spatial_model.limits[ip,0] = N.log10(float(pdict['min'])*scale) 
                if pdict.has_key('max'):
                    spatial_model.limits[ip,1] = N.log10(float(pdict['max'])*scale) 

            else:
                spatial_model.p[ip] = value*scale
                if pdict.has_key('min'):
                    spatial_model.limits[ip,0] = float(pdict['min'])*scale
                if pdict.has_key('max'):
                    spatial_model.limits[ip,1] = float(pdict['max'])*scale

            spatial_model.free[ip] = (pdict['free'] == '1')

            if 'error' in pdict.keys():
                err = float(pdict['error'])
                if spatial_model.log[ip]: 
                    spatial_model.cov_matrix[ip,ip] = (err/value*JAC)**2
                else:
                    spatial_model.cov_matrix[ip,ip] = err**2

        spatial_model.cache()
        return spatial_model


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

        elif name == 'BrokenPowerLaw':
            self.pname  = ['Prefactor', 'Index1', 'Index2', 'BreakValue']
            self.pfree  = [1,1,1,1]
            self.pscale = [1e-9,-1,-1,1]
            self.pmin   = [1e-5,0,0,30]
            self.pmax   = [1e4,5,5,5e5]
            self.pval   = [1,2,2,1e3]
            self.perr   = [0,0,0,0]

        elif name == 'BrokenPowerLaw2':
            self.pname  = ['Integral','Index1','Index2','BreakValue','LowerLimit','UpperLimit']
            self.pfree  = [1,1,1,1,0,0]
            self.pscale = [1e-7,-1,-1,1,1,1]
            self.pmin   = [1e-4,0,0,30,30,30]
            self.pmax   = [1e4,5,5,5e5,5e5,5e5]
            self.pval   = [2,1,1,1e3,1e2,1e5]
            self.perr   = [0,0,0,0,-1,-1]

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

        elif name == 'ConstantValue':
            self.pname  = ['Value']
            self.pfree  = [1]
            self.pscale = [1]
            self.pmin   = [0.1]
            self.pmax   = [10]
            self.pval   = [1]
            self.perr   = [0]

        else:
            raise Exception,'Unrecognized model %s'%(name)

    def find_param(self,pname):
        for iparam,param in enumerate(self.pname):
            if param == pname: return iparam
        return None

    def param_strings(self):
        err_strings = ['error="%s"'%(e) if (e>0) else '' for e in self.perr]
        return ['\t<parameter name="%s" value="%s" %s free="%d" max="%s" min="%s" scale="%s" />'%(e,f,m1,m2,n,s,v)
            for e,f,m1,m2,n,s,v in zip(self.pname,self.pval, err_strings,self.pfree,self.pmax,self.pmin,self.pscale,)]

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

    def process_model(self,model,xml_name=None,scaling=False):
        """model an instance of Model
        
           Allow an override of the xml_name to map the model on
           when there is ambiguity (i.e. for ConstantValue vs.
           FileFunction). """

        my_xml_name = xml_name # fix scope issue
        if model.name == 'ExpCutoff':
            model = convert_exp_cutoff(model)

        # map the Model instance onto an xml-style model
        if xml_name == None:
            if model.name == 'Constant': my_xml_name='ConstantValue'
            else:
                for l_xml_name,v in self.x2m.modict.iteritems():
                    if v == model.name:
                        my_xml_name = l_xml_name;break
                if v != model.name:
                    raise Exception,'Unable to find an XML model for %s'%(model.name)
        self.update(my_xml_name,scaling=scaling)

        # replace spectral parameters
        specparams = self.x2m.specdict[my_xml_name]
        vals,errs  = model.statistical(absolute=True)
        for iparam,param in enumerate(specparams):
            if self.debug: print 'Processing %s'%(param)
            index = self.find_param(param)
            if index is None:
                raise Exception,'Unrecognized parameter %s'%(param)
            self.pfree[index] = model.free[iparam]
            self.pval[index] = vals[iparam] / self.pscale[index]
            self.perr[index] = errs[iparam] / self.pscale[index]
            if (param == 'Index') or (param == 'Index1') or (param == 'Index2'):
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
        kwargparams = self.x2m.kwargdict[my_xml_name]
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
    if elem['type'] != 'SkyDirFunction':
        raise Exception("""The PointSource's SpatialModel must have type="PointSource".""")
    for p in elem.children:
        d[p['name']] = float(p['value'])
    return SkyDir(d['RA'],d['DEC'])

def makePSSpatialModel(ra,dec,tablevel=1):
    """Encode the spatial model in XML for a point source at ra,dec."""

    strings = [
        '<spatialModel type="SkyDirFunction">',
        '\t<parameter name="RA"  value="%s" free="0" max="360.0" min="-360.0" scale="1.0" />' %(ra),
        '\t<parameter name="DEC" value="%s" free="0" max="90" min="-90" scale="1.0" />' %(dec),
        '</spatialModel>'
        ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeDSConstantSpatialModel(tablevel=1):
    """Encode an isotropic model."""
    
    strings = [
        '<spatialModel type="ConstantValue">',
        '\t<parameter  name="Value" value="1.0" free="0" max="10.0" min="0.0"scale="1.0" />',
        '</spatialModel>'
	]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeDSMapcubeSpatialModel(filename='ERROR',tablevel=1):
    """Encode a mapcube model."""
    strings = [
        '<spatialModel file="%s" type="MapCubeFunction">'%(filename),
        '\t<parameter name="Normalization" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />',
	    '</spatialModel>'
    ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeExtendedSourceSpatialModel(es,tablevel=1):
    """Encode a spatial model."""
    strings = [
        '<spatialModel type="%s">' % es.pretty_name,
    ]

    xtsm = XML_to_SpatialModel()
    param_names = xtsm.spatialdict[es.pretty_name]
    if es.coordsystem == SkyDir.GALACTIC:
        param_names=['L','B']+param_names
    if es.coordsystem == SkyDir.EQUATORIAL:
        param_names=['RA','DEC']+param_names

    params,param_errors=es.statistical(absolute=True)
    err_strings = ['error="%s" '%(e) if (e>0) else '' for e in param_errors]
    min_params,max_params=N.transpose(es.get_limits(absolute=True))
    min_params[0],max_params[0]=[-360,360]
    min_params[1],max_params[1]=[-90,90]
    for param,err,free,min,max,name in zip(params,err_strings,
                                       es.free,min_params,max_params,
                                       param_names):
        strings.append('\t<parameter name="%s" value="%g" error="%g" %sfree="%d" max="%g" min="%g" scale="1.0" />' % \
                       (name,param,err,free,max,min))
    strings.append('</spatialModel>')
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

def parse_diffuse_sources(handler,diffdir=None):
    gds = get_diffuse_source
    ds = deque()
    xtm = XML_to_Model()
    for src in handler.sources:
        if src['type'] != 'DiffuseSource': continue
        spatial  = src.getChild('spatialModel')
        spectral = src.getChild('spectrum')
        name = src['name']
        if spatial['type'] == 'ConstantValue':
            if spectral['type'] == 'FileFunction':
                fname = str(os.path.expandvars(spectral['file']))
                mo = xtm.get_model(spectral)
                ds.append(gds('ConstantValue',None,mo,fname,name,diffdir=diffdir))
            elif (spectral['type'] == 'PowerLaw' ) or (spectral['type'] == 'PowerLaw2'):
                mo = xtm.get_model(spectral)
                ds.append(gds('ConstantValue',None,mo,None,name))
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
            ds.append(gds('MapCubeFunction',fname,mo,None,name,diffdir=diffdir))
            
        elif spatial['type'] in [ 'PseudoGaussian', 'Gaussian', 'Disk', 'PseudoDisk', 'NFW', 'PseudoNFW']:
            xtsm = XML_to_SpatialModel()
            spatial_model=xtsm.get_spatial_model(spatial)
            spectral_model=xtm.get_model(spectral)
            ds.append(ExtendedSource(name=name,
                                     model=spectral_model,
                                     spatial_model=spatial_model,
                                     leave_parameters=True))
        else:
            raise Exception('Diffuse spatial model "%s" not recognized' % spatial['type'])
    return ds

def parse_sources(xmlfile,diffdir=None):
    """Convenience function to parse an entire file into
       point sources and diffuse sources."""
    handler = parse_sourcelib(xmlfile)
    ps = parse_point_sources(handler)
    ds = parse_diffuse_sources(handler,diffdir=diffdir)
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

    if isinstance(ds,ExtendedSource):
        m2x.process_model(ds.smodel,scaling=True)
        specxml = m2x.getXML()
        skyxml = makeExtendedSourceSpatialModel(ds.spatial_model)
    elif isinstance(dm,DiffuseFunction):
        filename = os.path.abspath(dm.name())
        skyxml = makeDSMapcubeSpatialModel(filename=filename)
        m2x.process_model(ds.smodel,scaling=True)
        specxml = m2x.getXML()
    else:
        skyxml = makeDSConstantSpatialModel()
        if isinstance(dm,IsotropicSpectrum):
            filename = os.path.abspath(dm.name())
            m2x.process_model(ds.smodel,xml_name='FileFunction',scaling=True)
            specxml = m2x.getXML(spec_attrs='file=\"%s\"'%(filename.replace('\\','/')))
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
    s1 = '\n<source name="%s" type="DiffuseSource">\n'%(ds.name)
    s2 = '</source>'
    return ''.join([s1,specxml,skyxml,s2])
    
def unparse_diffuse_sources(diffuse_sources):
    """Convert a list of DiffuseSources into XML blurbs."""
    xml_blurbs = Stack()
    for ds in diffuse_sources:
        xml_blurbs.push(process_diffuse_source(ds))
    return xml_blurbs

def writeXML(stacks,filename, title='source_library'):
    """Write XML blurbs to a gtlike-style XML file."""
    f = open(filename,'wb')
    f.write('<source_library title="%s">'% title)
    for stack in stacks:
        for elem in stack:
            f.write(elem)
    f.write('\n</source_library>')

def writeROI(roi,filename):
    """ out the contents of an ROIAnalysis source model
       to a gtlike XML file."""
    source_xml = [unparse_point_sources(roi.psm.point_sources)]
    try:
        source_xml.append( unparse_diffuse_sources(roi.dsm.diffuse_sources))
    except AttributeError: 
        print 'warning: no diffuse sources found to write to xml'
    writeXML(source_xml,filename)
            
