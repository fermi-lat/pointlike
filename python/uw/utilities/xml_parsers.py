"""Class for parsing and writing gtlike-style source libraries.
   Barebones implementation; add additional capabilities as users need.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/xml_parsers.py,v 1.54 2012/01/30 22:44:35 wallacee Exp $

   author: Matthew Kerr
"""

import numpy as N
import xml.sax as x
from collections import deque
from uw.like.pointspec_helpers import PointSource,get_diffuse_source
from uw.like.roi_diffuse import DiffuseSource
from uw.like.roi_extended import ExtendedSource
from uw.like.Models import *
from uw.like2.models import *
from uw.like.SpatialModels import *
from uw.like.dark_matter import *
from skymaps import SkyDir,DiffuseFunction,IsotropicSpectrum,IsotropicPowerLaw
from os.path import join
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
        self.parent = None
        super(XMLElement,self).__init__(*args,**kwargs)

    def addChild(self,child):
        child.parent = self
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

        self.modict = dict(
                PowerLaw             = PowerLaw,
                PowerLaw2            = PowerLawFlux,
                BrokenPowerLaw       = BrokenPowerLaw,
                BrokenPowerLaw2      = BrokenPowerLawFlux,
                SmoothBrokenPowerLaw = SmoothBrokenPowerLaw,
                PLSuperExpCutoff     = PLSuperExpCutoff,
                Constant             = Constant, # should this have been ConstantValue key?
                ConstantValue        = Constant,
                FrontBackConstant    = FrontBackConstant,
                FileFunction         = Constant, # a big klugey
                LogParabola          = LogParabola,
                DMFitFunction        = DMFitFunction,
                )

        self.specdict = dict(
                PowerLaw             = ['Prefactor','Index'],
                PowerLaw2            = ['Integral','Index'],
                BrokenPowerLaw       = ['Prefactor', 'Index1', 'Index2', 'BreakValue'],
                BrokenPowerLaw2      = ['Integral', 'Index1', 'Index2', 'BreakValue'],
                SmoothBrokenPowerLaw = ['Prefactor', 'Index1', 'Index2', 'BreakValue'],
                PLSuperExpCutoff     = ['Prefactor','Index1','Cutoff','Index2'],
                ConstantValue        = ['Value'],
                FrontBackConstant    = ['Scale_front','Scale_back'],
                FileFunction         = ['Normalization'],
                LogParabola          = ['norm', 'alpha', 'beta', 'Eb'],
                DMFitFunction        = ['sigmav','mass'],
                )

        self.kwargdict = dict(
                PowerLaw             = [ ['Scale','e0' ] ],
                PowerLaw2            = [ ['LowerLimit','emin'], ['UpperLimit','emax'] ],
                BrokenPowerLaw       = [],
                BrokenPowerLaw2      = [ ['LowerLimit','emin'], ['UpperLimit','emax'] ],
                SmoothBrokenPowerLaw = [ ['Scale','e0' ], ['Beta','beta'] ],
                PLSuperExpCutoff     = [ ['Scale','e0' ] ],
                ConstantValue        = [],
                FrontBackConstant    = [],
                FileFunction         = [],
                LogParabola          = [],
                DMFitFunction        = [['norm','norm'], ['bratio','bratio'],
                                        ['channel0','channel0'], ['channel1','channel1']],
                )

    def get_model(self,xml_dict,source_name,index_offset=0):
        """ source_name is used for better error message printing. """

        specname = xml_dict['type']
        params   = xml_dict.children

        d = dict()
        for p in params:
            d[p['name']] = p

        # certain spectral models require the file to be set
        kwargs = {}
        if specname == 'DMFitFunction' and xml_dict.has_key('file'):
            kwargs['file']=str(xml_dict['file'])

        model = self.modict[specname](**kwargs)

        for ip,p in enumerate(self.specdict[specname]):
            try:
                pdict = d[p]
            except:
                raise Exception("For source %s, %s parameter %s not found in xml file." % (source_name,specname,p))

            scale = float(pdict['scale'])
            value = float(pdict['value'])
            if (p == 'Index') or (p == 'Index1') or (p == 'Index2' and self.modict[specname]!=PLSuperExpCutoff):
                # gtlike uses a neg. index internally so scale > 0
                # means we need to take the negative of the value
                if scale > 0: value = -value
                else:         scale = -scale
                value += index_offset
                model.index_offset = index_offset
	    if value==0 : value=1.e-5
            if value*scale<=0: raise Exception('For source %s, %s parameter %s cannot be negative' % (source_name,specname,p))
            if N.isnan(value*scale): raise Exception('For source %s, %s parameter %s is NaN' % (source_name,specname,p))
            model.setp(ip, value*scale)
            model.free[ip] = (pdict['free'] == '1')
            if 'error' in pdict.keys():
                err = float(pdict['error'])*scale
                model.cov_matrix[ip,ip] = (err/value*JAC)**2

        for p in self.kwargdict[specname]:
            pdict = d[p[0]]

            if pdict['free'] != '0':
                # Sanity check on validity of xml 
                raise Exception('For source %s, %s parameter %s cannot be fit (must be free="0")' % (source_name,specname,p[0]))

            model.__dict__[p[1]] = float(pdict['value'])*float(pdict['scale'])

        return model

class XML_to_SpatialModel(object):

    def __init__(self):

        # exclude position from spatialdict
        self.spatialdict = dict({
            'Gaussian'           : ['Sigma'],
            'PseudoGaussian'     : [],
            'Disk'               : ['Sigma'],
            'PseudoDisk'         : [],
            'Ring'               : ['Sigma','Fraction'],
            'NFW'                : ['Sigma'],
            'PseudoNFW'          : [],
            'EllipticalGaussian' : ['MajorAxis','MinorAxis','PositionAngle'],
            'EllipticalDisk'     : ['MajorAxis','MinorAxis','PositionAngle'],
            'EllipticalRing'     : ['MajorAxis','MinorAxis','PositionAngle','Fraction'],
            'SpatialMap'         : [],
            })

        self.spatialdict['PseudoEllipticalGaussian'] = self.spatialdict['PseudoGaussian']
        self.spatialdict['RadiallySymmetricEllipticalGaussian'] = self.spatialdict['Gaussian']
        self.spatialdict['PseudoEllipticalDisk'] = self.spatialdict['PseudoDisk']
        self.spatialdict['RadiallySymmetricEllipticalDisk'] = self.spatialdict['Disk']

    def get_spatial_model(self,xml_dict,diffdir=None):
        """ Kind of like getting a spectral model, but
            instead returns a SpatialModel object. 
            
            Since extended ources don't exist in gtlike, no need
            to keep backwards compatability. Just trying
            to keep a farmilliar interface and make the
            naming convention similar to gtobssim. """

        spatialname = xml_dict['type']
        params      = xml_dict.children

        if spatialname == 'SpatialMap':
            # For spatial maps, ignore any of the parameters.
            if diffdir:
                file = join(diffdir,os.path.basename(str(xml_dict['file'])))
            else:
                file = str(xml_dict['file'])
            return SpatialMap(file=file)

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
            value = float(pdict['value'])*scale

            if spatial_model.log[ip]:
                if value<=0: 
                    raise Exception("Unable to parse %s's parameter %s. Value must be greater than 0" % (spatialname,pdict['name']))
                spatial_model.p[ip] = N.log10(value) 
                if pdict.has_key('min'):
                    minimum=float(pdict['min'])*scale
                    if minimum<=0: 
                        raise Exception("Unable to parse %s's parameter %s. Minimum fit range must be greater than 0" % (spatialname,pdict['name']))
                    spatial_model.limits[ip,0] = N.log10(minimum) 
                if pdict.has_key('max'):
                    maximum=float(pdict['max'])*scale
                    if maximum<=0: 
                        raise Exception("Unable to parse %s's parameter %s. Maximum fit range must be greater than 0" % (spatialname,pdict['name']))
                    spatial_model.limits[ip,1] = N.log10(maximum) 

            else:
                spatial_model.p[ip] = value
                if pdict.has_key('min'):
                    spatial_model.limits[ip,0] = float(pdict['min'])*scale
                if pdict.has_key('max'):
                    spatial_model.limits[ip,1] = float(pdict['max'])*scale

            if spatial_model.limits[ip][0] >= spatial_model.limits[ip][1]:
                raise Exception("Unable to parse %s's parameter %s. Maximum fit range must be greater minimum fit range." % (spatialname,pdict['name']))

            spatial_model.free[ip] = (pdict['free'] == '1')

            if 'error' in pdict.keys():
                err = float(pdict['error'])*scale
                if spatial_model.log[ip]: 
                    spatial_model.cov_matrix[ip,ip] = (err/value*JAC)**2
                else:
                    spatial_model.cov_matrix[ip,ip] = err**2

        spatial_model.cache()
        return spatial_model


class Model_to_XML(object):
    """Convert a Model instance to an XML entry.
       
       This class should be extended as more use cases arise."""
    
    def __init__(self,debug=False, strict=False):
        self.x2m = XML_to_Model()
        self.debug = debug
        self.strict = strict

    def update(self,name,scaling=False):
        self.name = name

        if name == 'PowerLaw2':
            self.pname  = ['Integral','Index','LowerLimit','UpperLimit']
            self.pfree  = [1,1,0,0]
            self.pscale = [1e-10,-1,1,1]
            self.pmin   = [1e-4,0,30,30]
            self.pmax   = [1e4,5,5e5,5e5]
            self.pval   = [2,1,1e3,1e5]
            self.perr   = [0,0,-1,-1]
            self.oomp   = [1,0,0,0]

        elif name == 'PowerLaw':
            self.pname  = ['Prefactor','Index','Scale']
            self.pfree  = [1,1,0]
            self.perr   = [0,0,-1]
            if scaling:
                self.pscale = [1,-1,1] ; self.pmin  = [0.1,-1,30]
                self.pmax  = [10,1,5e5]; self.pval = [1,0,1e3]
                self.oomp   = [0,0,0]
            else:
                self.pscale = [1e-9,-1,1]
                self.pmin   = [1e-6,0,30]
                self.pmax   = [1e4,5,5e5]
                self.pval   = [1,2,1e3]
                self.oomp   = [1,0,0]

        elif name == 'BrokenPowerLaw':
            self.pname  = ['Prefactor', 'Index1', 'Index2', 'BreakValue']
            self.pfree  = [1,1,1,1]
            self.pscale = [1e-9,-1,-1,1]
            self.pmin   = [1e-5,0,0,30]
            self.pmax   = [1e4,5,5,5e5]
            self.pval   = [1,2,2,1e3]
            self.perr   = [0,0,0,0]
            self.oomp   = [1,0,0,1]

        elif name == 'BrokenPowerLaw2':
            self.pname  = ['Integral','Index1','Index2','BreakValue','LowerLimit','UpperLimit']
            self.pfree  = [1,1,1,1,0,0]
            self.pscale = [1e-7,-1,-1,1,1,1]
            self.pmin   = [1e-4,0,0,30,30,30]
            self.pmax   = [1e4,5,5,5e5,5e5,5e5]
            self.pval   = [2,1,1,1e3,1e2,1e5]
            self.perr   = [0,0,0,0,-1,-1]
            self.oomp   = [1,0,0,1,0,0]

        elif name == 'PLSuperExpCutoff':
            self.pname = ['Prefactor','Index1','Cutoff','Index2','Scale']
            self.pfree  = [1,1,1,0,0]
            self.pscale = [1e-9,-1,1000.,1,1]
            self.pmin   = [1e-4,0,0.1,0,30]
            self.pmax   = [1e4,5,3e5,5,5e5]
            self.pval   = [1,2,1,1,1000]
            self.perr   = [0,0,0,0,-1]
            self.oomp   = [1,0,1,0,0]

        elif name == 'SmoothBrokenPowerLaw':
            self.pname  = ['Prefactor', 'Index1', 'Index2', 'BreakValue', 'Scale', 'Beta']
            self.pfree  = [1,1,1,1,0,0]
            self.pscale = [1e-9,-1,-1,1,1,1]
            self.pmin   = [1e-5,0,0,30,30,-10]
            self.pmax   = [1e4,5,5,5e5,5e5,10]
            self.pval   = [1,2,2,1e3,1e3,0.1]
            self.perr   = [0,0,0,0,-1,-1]
            self.oomp   = [1,0,0,1,0,0]

        elif name == 'FileFunction':
            self.pname = ['Normalization']
            self.pfree = [1]
            self.pscale = [1]
            self.pmin = [0.1]
            self.pmax = [10]
            self.pval = [1]
            self.perr = [0]
            self.oomp = [0]

        elif name == 'ConstantValue':
            self.pname  = ['Value']
            self.pfree  = [1]
            self.pscale = [1]
            self.pmin   = [0.1]
            self.pmax   = [10]
            self.pval   = [1]
            self.perr   = [0]
            self.oomp   = [0]

        elif name == 'FrontBackConstant':
            self.pname  = ['Scale_front','Scale_back']
            self.pfree  = [1,1]
            self.pscale = [1,1]
            self.pmin   = [0,0]
            self.pmax   = [10,10]
            self.pval   = [1,1]
            self.perr   = [0,0]
            self.oomp   = [0,0]

        elif name == 'LogParabola':
            self.pname  = ['norm', 'alpha', 'beta', 'Eb' ]
            self.pfree  = [1,1,1,1]
            self.perr   = [0,0,0,0]
            self.pscale = [1e-9,1,1,1]
            self.pmin   = [1e-6,0,0,30]
            self.pmax   = [1e4,5,5,5e5]
            self.pval   = [1,2,0,1e3]
            self.oomp   = [1,0,0,1]

        elif name == 'DMFitFunction':
            self.pname  = ['sigmav','mass', 'norm', 'bratio', 'channel0', 'channel1']
            self.pfree  = [       1,     1,      0,        0,          0,          0]
            self.perr   = [       0,     0,     -1,       -1,         -1,         -1]
            self.pscale = [   1e-25,     1,      1,     1e17,          1,          1]
            self.pmin   = [       0,     1,   1e-5,        0,          0,          0]
            self.pmax   = [     1e6,   1e4,    1e5,        1,         10,         10]
            self.pval   = [       1,   100,      1,        1,          1,          1]
            self.oomp   = [       1,     0,      1,        0,          0,          0]

        else:
            raise Exception,'Unrecognized model %s'%(name)

    def find_param(self,pname):
        for iparam,param in enumerate(self.pname):
            if param == pname: return iparam
        return None

    def param_strings(self):
        err_strings = ['error="%s" '%(e) if (e>0) else '' for e in self.perr]
        return ['\t<parameter name="%s" value="%s" %sfree="%d" max="%s" min="%s" scale="%s" />'%(e,f,m1,m2,n,s,v)
            for e,f,m1,m2,n,s,v in zip(self.pname,self.pval, err_strings,self.pfree,self.pmax,self.pmin,self.pscale,)]

    def getXML(self,tablevel=1,spec_attrs='',comment_string=None):
        cstring = [] if comment_string is None else [comment_string]
        strings = ['<spectrum %s type="%s">'%(" ".join([spec_attrs,self.extra_attrs]),self.name)] + \
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

    def process_scale(self,val,index):
        if not self.oomp[index]: return
        # else adjust scale to nearest order of magnitude
        scale = round(N.log10(val))

        if scale == float('-inf'): return # happens if a parameter is too close to 0.

        self.pscale[index] = 10**scale
        self.pmin[index] = 1e-2
        self.pmax[index] = 1e2

    def process_model(self,model,xml_name=None,scaling=False):
        """model an instance of Model
        
           Allow an override of the xml_name to map the model on
           when there is ambiguity (i.e. for ConstantValue vs.
           FileFunction). """

        my_xml_name = xml_name # fix scope issue
        if model.name == 'ExpCutoff':
            model = convert_exp_cutoff(model)
        elif model.name == 'LogParabola' and model[2]<= 2e-3 and not model.free[2]:
            # convert to equivalent powerlaw
            model = model.create_powerlaw()

        # map the Model instance onto an xml-style model
        if xml_name == None:
            if model.name == 'Constant': my_xml_name='ConstantValue'
            else:
                for l_xml_name,v in self.x2m.modict.iteritems():
                    if v == type(model):
                        my_xml_name = l_xml_name;break
                if v != type(model):
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
            self.process_scale(vals[iparam],index)
            self.pfree[index] = model.free[iparam]
            self.pval[index] = vals[iparam] / self.pscale[index]
            self.perr[index] = errs[iparam] / self.pscale[index]
            # Note that logParabola is already consistent with gtlike
            if (param == 'Index') or (param == 'Index1') or (param == 'Index2' and model.name!='PLSuperExpCutoff'):
                self.process_photon_index(model,index,vals[iparam],errs[iparam])
            if self.pval[index] < self.pmin[index]:
                msg = 'Found %s=%s < %s, minimum allowed value'%(param, str(self.pval[index]),str(self.pmin[index]))
                if self.strict: raise Exception(msg)
                print 'WARNING: %s, \n\tSetting parameter value to minimum.' %msg
                self.pval[index] = self.pmin[index]
            if self.pval[index] > self.pmax[index]:
                msg = 'Found %s=%s > %s, maximum allowed value'%(param, str(self.pval[index]),str(self.pmax[index]))
                if self.strict: raise Exception(msg)
                print 'Warning %s, \n\tSetting parameter value to maximum.'%msg
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

        # this will also be necessary for FileSpectrum objects.
        if isinstance(model,DMFitFunction): 
            self.extra_attrs='file="%s"' % model.file
        else:
            self.extra_attrs=''

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
        '\t<parameter  name="Value" value="1.0" free="0" max="10.0" min="0.0" scale="1.0" />',
        '</spatialModel>'
    ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeDSMapcubeSpatialModel(filename='ERROR',tablevel=1):
    """Encode a mapcube model."""
    strings = [
        '<spatialModel file="%s" type="MapCubeFunction">'%(filename.replace('\\','/')),
        '\t<parameter name="Normalization" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />',
        '</spatialModel>'
    ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeExtendedSourceSpatialModel(es,expand_env_vars,tablevel=1):
    """Encode a spatial model."""
    if es.name != 'SpatialMap':
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
            strings.append('\t<parameter name="%s" value="%g" %sfree="%d" max="%g" min="%g" scale="1.0" />' % \
                           (name,param,err,free,max,min))
    else:
        file=os.path.expandvars(es.file) if expand_env_vars else es.file
        strings = [
            '<spatialModel type="%s" file="%s" >' % (es.pretty_name,file),
            '\t<parameter name="Prefactor" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />'
        ]

    strings.append('</spatialModel>')
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def parse_sourcelib(xml):
    parser = x.make_parser()
    handler = SourceHandler()
    parser.setContentHandler(handler)
    if isinstance(xml,list):
	    [parser.parse(xmlfile) for xmlfile in xml]
    else:
	    parser.parse(xml)
    return handler

def parse_point_sources(handler,roi_dir,max_roi):
    point_sources = deque()
    xtm = XML_to_Model()
    for src in handler.sources:
        if src['type'] != 'PointSource': continue
        name = str(src['name'])
        sd = get_skydir(src.getChild('spatialModel'))
        mo = xtm.get_model(src.getChild('spectrum'),name)
        if None in [roi_dir,max_roi] or N.degrees(sd.difference(roi_dir))<max_roi:
            point_sources.append(PointSource(sd,name,mo,leave_parameters=True))
    return list(point_sources)

def parse_diffuse_sources(handler,diffdir=None):
    gds = get_diffuse_source
    ds = deque()
    xtm = XML_to_Model()
    for src in handler.sources:
        #if src['type'] == 'CompositeDiffuseSource':
        #    handler = SourceHandler()
        if src['type'] != 'DiffuseSource': continue
        spatial  = src.getChild('spatialModel')
        spectral = src.getChild('spectrum')
        name = str(src['name'])
        if spatial['type'] == 'ConstantValue':
            if spectral['type'] == 'FileFunction':
                fname = str(os.path.expandvars(spectral['file']))
                mo = xtm.get_model(spectral,name)
                ds.append(gds('ConstantValue',None,mo,fname,name,diffdir=diffdir))
            elif (spectral['type'] == 'PowerLaw' ) or (spectral['type'] == 'PowerLaw2'):
                mo = xtm.get_model(spectral,name)
                ds.append(gds('ConstantValue',None,mo,None,name))
            else:
                raise Exception,'Isotropic model not implemented'
    
        elif spatial['type'] == 'MapCubeFunction':
            fname = str(os.path.expandvars(spatial['file']))
            if spectral['type'] == 'ConstantValue':
                mo = xtm.get_model(spectral,name)
            elif spectral['type'] == 'PowerLaw' or spectral['type'] == 'PowerLaw2':
                mo = xtm.get_model(spectral,name,index_offset=1)
            else:
                raise Exception('Non-isotropic model "%s" not implemented' % spatial['type'])
            ds.append(gds('MapCubeFunction',fname,mo,None,name,diffdir=diffdir))
            
        else:
            xtsm = XML_to_SpatialModel()

            if spatial['type'] in xtsm.spatialdict.keys():

                spatial_model=xtsm.get_spatial_model(spatial,diffdir=diffdir)
                spectral_model=xtm.get_model(spectral,name)
                ds.append(ExtendedSource(name=name,
                                         model=spectral_model,
                                         spatial_model=spatial_model))
            else:
                raise Exception('Diffuse spatial model "%s" not recognized' % spatial['type'])
    return list(ds)

def parse_sources(xmlfile,diffdir=None,roi_dir=None,max_roi=None):
    """Convenience function to parse an entire file into
       point sources and diffuse sources."""
    handler = parse_sourcelib(xmlfile)
    ps = parse_point_sources(handler,roi_dir,max_roi)
    ds = parse_diffuse_sources(handler,diffdir=diffdir)
    return ps,ds

def unparse_point_sources(point_sources, strict=False, properties=lambda x:''):
    """Convert a list (or other iterable) of PointSource objects into XML.
        strict : bool
            set True to generate exception, error message identifying offending source, reason
        properties : a function
            the function, if specified, returns a string for the source element with properties, like TS
    """
    xml_blurbs = Stack()
    m2x = Model_to_XML(strict=strict)
    for ps in point_sources:
        skyxml = makePSSpatialModel(ps.skydir.ra(),ps.skydir.dec())
        try:
            m2x.process_model(ps.model)
        except Exception, emsg:
            print 'Failed to process source %s: %s' %(ps.name, emsg)
        specxml = m2x.getXML()
        s1 = '\n<source name="%s" type="PointSource" %s >\n'%(ps.name, properties(ps))
        s2 = '</source>'
        xml_blurbs.push(''.join([s1,specxml,skyxml,s2]))
    return xml_blurbs


def process_diffuse_source(ds,convert_extended=False,expand_env_vars=True,filename=None,ctype=None):
    """Convert an instance of DiffuseSource into an XML blurb."""
    m2x = Model_to_XML()
    dm = ds.dmodel
    if hasattr(dm,'__len__') and len(dm)==1: dm = dm[0]

    if isinstance(ds,ExtendedSource) or hasattr(ds,'spatial_model'):
        m2x.process_model(ds.smodel,scaling=False)
        specxml  = m2x.getXML()
        spatial  = ds.spatial_model
        spectral = ds.smodel
        if convert_extended and not isinstance(spatial,SpatialMap): 
            if spatial.__dict__.has_key('original_template') and N.all(spatial.original_parameters == spatial.p):
                # Kludge! this is incase the xml was read in from the 
                # pointspec_helpers.ExtendedSourceArchive and should be saved 
                # out with the original template.
                spatial=SpatialMap(file=spatial.original_template)
            else:
                folder=os.path.dirname(filename or os.getcwd())
                template_name=folder+os.sep if folder != '' else ''
                template_name+='template_%s_%s_%s.fits' % (ds.name.replace(' ','_'),
                                                           spatial.pretty_name, 
                                                           spectral.pretty_name)
                spatial = convert_spatial_map(spatial,template_name)
                spatial.file = template_name
        skyxml = makeExtendedSourceSpatialModel(spatial,expand_env_vars)
        if isinstance(spatial,SpatialMap) and not N.all(spatial.p==spatial.init_p):
            print 'Warning: When saving out SpatialMap object which has been localized, the original unmoved template is saved in the xml model.'

    elif isinstance(dm,DiffuseFunction):
        filename = os.path.abspath(dm.name())
        skyxml = makeDSMapcubeSpatialModel(filename=filename)
        m2x.process_model(ds.smodel,scaling=True)
        specxml = m2x.getXML()
    else:
        skyxml = makeDSConstantSpatialModel()
        #Handle the case where the isotropic is specified by separate front and back components
        if not hasattr(dm,'__len__'): dm = [dm]
        if len(dm)==1:
            ctypes=[-1]
            specxml=''
        elif len(dm)==2:
            ctypes = [0,1]
            specxml='\t<spectrum type="CompositeSpectrum">\n'
        else:
            raise Exception("Don't know how to handle isotropic model with >2 components")
        for m,ct in zip(dm,ctypes):
            if isinstance(m,IsotropicSpectrum):
                filename = os.path.abspath(m.name())
                m2x.process_model(ds.smodel,xml_name='FileFunction',scaling=True)
                m2x.extra_attrs+='ctype="{0}"'.format(ct)
                specxml += m2x.getXML(spec_attrs='file=\"%s\"'%(filename.replace('\\','/')),tablevel=len(dm))
            elif isinstance(m,IsotropicPowerLaw):
                sd = SkyDir()
                iflux = m.integral(sd,100,3e5)
                spind = m.value(sd,100)*100/iflux + 1
                pl = PowerLawFlux()
                pl.p[0] = ds.smodel.p[0] + N.log10(iflux)
                pl.p[1] = N.log10(spind+10**(ds.smodel.p[1]-ds.smodel.index_offset))
                pl.cov_matrix = ds.smodel.cov_matrix.copy() #correct?
                m2x.process_model(pl,scaling=False)
                specxml += m2x.getXML(tablevel=len(dm))
            else:
                raise Exception('Did not recognize %s'%(ds.name))
        if len(dm)==2:
            specxml+='\t</spectrum>\n'
    s1 = '\n<source name="%s" type="DiffuseSource">\n'%(ds.name)
    s2 = '</source>'
    return ''.join([s1,specxml,skyxml,s2])
    
def unparse_diffuse_sources(diffuse_sources,convert_extended=False,expand_env_vars=True,filename=None):
    """Convert a list of DiffuseSources into XML blurbs."""
    xml_blurbs = Stack()
    for ds in diffuse_sources:
        xml_blurbs.push(process_diffuse_source(ds,
                                               convert_extended=convert_extended,
                                               expand_env_vars=expand_env_vars,
                                               filename=filename))
    return xml_blurbs

def writeXML(stacks,filename, title='source_library'):
    """Write XML blurbs to a gtlike-style XML file."""
    f = open(filename,'wb') if type(filename)==str else filename
    f.write('<source_library title="%s">'% title)
    for stack in stacks:
        for elem in stack:
            f.write(elem)
    f.write('\n</source_library>')
    f.close()

def writeROI(roi,filename,strict=False,convert_extended=False,expand_env_vars=False):
    """ Out the contents of an ROIAnalysis source model
        to a gtlike XML file.

        the strict flag raises an exception if any of the spectral
        parameters are outside of the fit range.

        The convert_extended flag converts extended sources to be fully
        compliant with gtlike.  By default, extended sources are saved
        out with a format incompatable with gtlike. This is because
        each diffuse source has its spatial parameters saved out in
        the xml file. This is useful for saving errors on the fit
        spatial parameters as well as for easibly being read back into
        pointlike. OTOH, it is sometimes desirable for extended source
        output to be stirctly compatable with gtlike. This can be done
        with the convert_extended flag, which converts all extended
        sources cto SpatialMap objects before the xml is created. 
        
        Currently, expand_env_vars only applies for extended
        sources. There really isn't a good way in pointlike to have
        diffuse sources mantain a memory of environment varaibles in their
        pathname, but this would be a nice addition in the future. """
    source_xml = [unparse_point_sources(roi.psm.point_sources, strict=strict)]
    if len(roi.dsm.diffuse_sources)>0:
        source_xml.append(unparse_diffuse_sources(roi.dsm.diffuse_sources,
                                                  convert_extended=convert_extended,
                                                  expand_env_vars=expand_env_vars,
                                                  filename=filename))
    writeXML(source_xml,filename)
