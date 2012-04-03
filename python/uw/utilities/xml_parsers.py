"""Class for parsing and writing gtlike-style sourceEQUATORIAL libraries.
   Barebones implementation; add additional capabilities as users need.

   $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/utilities/xml_parsers.py,v 1.64 2012/03/27 16:41:34 lande Exp $

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
from uw.darkmatter.spatial import *
from skymaps import SkyDir,DiffuseFunction,IsotropicSpectrum,IsotropicPowerLaw,IsotropicConstant
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
    
    def __init__(self, pattern = 'source'):
        self.outerElements = deque()
        self.sources       = deque()
        self.pattern       = pattern

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
        if t.name == self.pattern:
            self.sources.append(t)

class XML_to_Model(object):
    """ Map a parsed XML model onto the Model classes.
    
        This class should be extended as more use cases arise.

        Here is a simple example to load in a point source.
        First, we have to setup the xml parser:

            >>> def xml2model(xml):
            ...     from StringIO import StringIO
            ...     parser = x.make_parser()
            ...     xtm = XML_to_Model()
            ...     handler = SourceHandler(pattern='spectrum')
            ...     parser.setContentHandler(handler)
            ...     parser.parse(StringIO(xml))
            ...     return xtm.get_model(handler.sources[0],'source')

        Now, we can parse spectral models:

            >>> model=xml2model('''
            ...     <spectrum   type="PowerLaw">
            ...         <parameter name="Prefactor" value="1.0" free="1" max="100.0" min="0.01" scale="1e-11" />
            ...         <parameter name="Index" value="2.0" free="1" max="5" min="0" scale="-1" />
            ...         <parameter name="Scale" value="1000.0" free="0" max="1000.0" min="1000.0" scale="1" />
            ...     </spectrum>''')
            >>> real_model=PowerLaw()
            >>> np.allclose(model.get_all_parameters(), real_model.get_all_parameters())
            True

"""

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

    extended_sources = [Gaussian, PseudoGaussian, Disk,
                             PseudoDisk, Ring, NFW, PingNFW,
                             PseudoPingNFW,
                             EllipticalGaussian, EllipticalDisk, 
                             EllipticalRing, SpatialMap,
                             RadialProfile]

    spatialdict = {i.__name__:i for i in extended_sources}

    @staticmethod
    def get_spatial_model(xml_dict,diffdir=None):
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

        input_kwargs = dict()

        if spatialname == 'SpatialMap' or spatialname == 'RadialProfile':
            # For spatial maps, ignore any of the parameters.
            if diffdir:
                file = join(diffdir,os.path.basename(str(xml_dict['file'])))
            else:
                file = str(xml_dict['file'])
            input_kwargs['file'] = file

        if d.has_key('RA') and d.has_key('DEC'):
            input_kwargs['coordsystem']=SkyDir.EQUATORIAL
        elif d.has_key('L') and d.has_key('B'):
            input_kwargs['coordsystem']=SkyDir.GALACTIC
    
        if not spatialname in XML_to_SpatialModel.spatialdict:                                
            raise Excception("Unrecognized spatial model %s" % spatialname)

        spatial_model = XML_to_SpatialModel.spatialdict[spatialname](**input_kwargs)

        if spatialname == 'SpatialMap':
            return spatial_model

        for p in spatial_model.param_names:
            pdict =d[p]

            scale = float(pdict['scale'])
            value = float(pdict['value'])


            spatial_model[p] = scale*value

            free = pdict['free'] == '1'
            spatial_model.set_free(p, free)

        return spatial_model


class Model_to_XML(object):
    """Convert a Model instance to an XML entry.
       
       This class should be extended as more use cases arise.
       
        Here is a simple test making a PowerLaw model:

            >>> m2x = Model_to_XML()
            >>> m2x.process_model(PowerLaw())
            >>> print m2x.getXML(tablevel=0).replace('\\t',' '*4).strip()
            <spectrum   type="PowerLaw">
                <parameter name="Prefactor" value="1.0" free="1" max="100.0" min="0.01" scale="1e-11" />
                <parameter name="Index" value="2.0" free="1" max="5" min="0" scale="-1" />
                <parameter name="Scale" value="1000.0" free="0" max="1000.0" min="1000.0" scale="1" />
            </spectrum>

   """
    
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
            self.pmin   = [0.001]
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
                    #Check the class name instead of the class, to account for FrontBackConstant moving from like2.models to like.Models
                    if v.__name__ == model.__class__.__name__: 
                        my_xml_name = l_xml_name;break
                if v.__name__ != model.__class__.__name__:
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

def make_position_params(skydir, ra_free=False, dec_free=False):
    ra = skydir.ra()
    dec = skydir.dec()
    return ['\t<parameter name="RA"  value="%s" free="%s" max="360.0" min="-360.0" scale="1.0" />' %(ra, int(ra_free)),
            '\t<parameter name="DEC" value="%s" free="%s" max="90" min="-90" scale="1.0" />' %(dec, int(ra_free))]

def makePSSpatialModel(skydir, tablevel=1):
    """Encode the spatial model in XML for a point source at ra,dec."""

    strings = [ '<spatialModel type="SkyDirFunction">', ] + \
            make_position_params(skydir) + \
            [ '</spatialModel>' ]
    return ''.join([decorate(st,tablevel=tablevel) for st in strings])

def makeDSConstantSpatialModel(value=1.0, tablevel=1):
    """Encode an isotropic model."""
    
    strings = [
        '<spatialModel type="ConstantValue">',
        '\t<parameter  name="Value" value="%s" free="0" max="10.0" min="0.0" scale="1.0" />' % value,
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
    """ Create a spatial model.

        Test saving out extended source:

            >>> def make(ds):
            ...     # remove newlines and convert tabs to spaces.
            ...     m=makeExtendedSourceSpatialModel(ds,expand_env_vars=False, tablevel=0)
            ...     return m.replace('\\t',' '*4).strip()

        Save Disk source:
    
            >>> disk = Disk(sigma=1, ra=244, dec=57)
                    
            >>> print make(disk)
            <spatialModel type="Disk">
                <parameter name="RA" value="244" free="1" max="360" min="-360" scale="1.0" />
                <parameter name="DEC" value="57" free="1" max="180" min="-180" scale="1.0" />
                <parameter name="Sigma" value="1" free="1" max="3" min="1e-10" scale="1.0" />
            </spatialModel>

            >>> disk.set_cov_matrix(np.asarray([[1,0,0],[0,1,0],[0,0,1]]))
            >>> print disk.error('sigma')
            2.30258509299
            >>> print make(disk)
            <spatialModel type="Disk">
                <parameter name="RA" value="244" error="1.0" free="1" max="360" min="-360" scale="1.0" />
                <parameter name="DEC" value="57" error="1.0" free="1" max="180" min="-180" scale="1.0" />
                <parameter name="Sigma" value="1" error="2.30258509299" free="1" max="3" min="1e-10" scale="1.0" />
            </spatialModel>

        When the coordsystem is GALACTIC, the output will be in L & B:

            >>> disk = Disk(sigma=1, l=244, b=57)
            >>> print make(disk)
            <spatialModel type="Disk">
                <parameter name="L" value="244" free="1" max="360" min="-360" scale="1.0" />
                <parameter name="B" value="57" free="1" max="180" min="-180" scale="1.0" />
                <parameter name="Sigma" value="1" free="1" max="3" min="1e-10" scale="1.0" />
            </spatialModel>


        Save SpatialModel source:

            >>> from tempfile import NamedTemporaryFile
            >>> temp = NamedTemporaryFile(delete=True)
            >>> disk.save_template(temp.name)
            >>> map=SpatialMap(file=temp.name)
            >>> print make(map).replace(temp.name, '<FILENAME>')
            <spatialModel type="SpatialMap" file="<FILENAME>" >
                <parameter name="Prefactor" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />
            </spatialModel>

        Save RadialProfile source:

            >>> temp = NamedTemporaryFile(delete=True)
            >>> disk.save_profile(temp.name)
            >>> profile=RadialProfile(file=temp.name)
            >>> print make(profile).replace(temp.name, '<FILENAME>')
            <spatialModel type="RadialProfile" file="<FILENAME>" >
                <parameter name="Normalization" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />
                <parameter name="RA" value="0" free="1" max="360" min="-360" scale="1.0" />
                <parameter name="DEC" value="0" free="1" max="180" min="-180" scale="1.0" />
            </spatialModel>

    """
    if es.name == 'SpatialMap' or es.name == 'RadialProfile':
        file=os.path.expandvars(es.file) if expand_env_vars else es.file
        strings = [
            '<spatialModel type="%s" file="%s" >' % (es.pretty_name,file),
        ]

    else:
        strings = [ '<spatialModel type="%s">' % es.pretty_name, ]

    if es.name == 'SpatialMap':
        strings += [ '\t<parameter name="Prefactor" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />' ]

    if es.name == 'RadialProfile':
        strings += [ '\t<parameter name="Normalization" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />' ]

    if es.name != 'SpatialMap':
        # All spatialm models but SpatialMap have RA & DEC

        for name in es.param_names:
            param = es[name]
            err = es.error(name)
            err = 'error="%s" ' % err if err > 0 else ''
            free = es.get_free(name)

            if name == 'L' or name == 'RA':
                min,max=-360,360
            elif name == 'B' or name == 'DEC':
                min,max=-180,180
            else:
                min,max=es.get_limits(absolute=True)[es.mapper(name)]

            strings.append('\t<parameter name="%s" value="%g" %sfree="%d" max="%g" min="%g" scale="1.0" />' % \
                           (name,param,err,free,max,min))

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
    """ Some simple testing. 

        In pointlike, we can load in the example point-like source.

        (XML from http://www.slac.stanford.edu/exp/glast/wb/prod/pages/sciTools_likelihoodTutorial/xmlModelDefinitions.htm)

            >>> from StringIO import StringIO
            >>> def l(xml): 
            ...     ps = parse_point_sources(parse_sourcelib(StringIO(xml)),roi_dir=None,max_roi=None)
            ...     assert len(ps) == 1
            ...     return ps[0]

            >>> ps=l('''<source name="PowerLaw_source" type="PointSource">  
            ...           <spectrum type="PowerLaw">
            ...             <parameter free="1" max="1000.0" min="0.001" name="Prefactor" scale="1e-09" value="1"/>
            ...             <parameter free="1" max="-1.0" min="-5." name="Index" scale="1.0" value="-2.1"/>
            ...             <parameter free="0" max="2000.0" min="30.0" name="Scale" scale="1.0" value="100.0"/>
            ...           </spectrum>      
            ...           <spatialModel type="SkyDirFunction">
            ...             <parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="83.45"/>
            ...             <parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="21.72"/>
            ...           </spatialModel>    
            ...         </source>''')
            >>> print ps.name
            PowerLaw_source
            >>> print ps.skydir.ra(), ps.skydir.dec()
            83.45 21.72
            >>> print ps.model.name
            PowerLaw
            >>> print ps.model['Norm'], ps.model['Index']
            1e-09 2.1

    """
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
    """ 
        Some simple testing:

        First, we define a simpler helper function

            >>> from StringIO import StringIO
            >>> def l(xml): 
            ...     ds = parse_diffuse_sources(parse_sourcelib(StringIO(xml)))
            ...     assert len(ds) == 1
            ...     return ds[0]


        Example loading in the galactic diffuse emission

        (XML from http://www.slac.stanford.edu/exp/glast/wb/prod/pages/sciTools_binnedLikelihoodTutorial/binnedLikelihood_v01archived.htm)
    
            >>> ds=l('''
            ...   <source name="Galpro Diffuse" type="DiffuseSource">
            ...     <spectrum type="ConstantValue">
            ...       <parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>
            ...     </spectrum>
            ...     <spatialModel file="$GLAST_EXT/diffuseModels/v2r0p1/ring_2year_P76_v0.fits" type="MapCubeFunction">
            ...       <parameter free="0" max="1000.0" min="0.001" name="Normalization" scale="1.0" value="1.0"/>
            ...     </spatialModel>
            ...   </source>''')
            >>> print ds.smodel.name
            Constant
            >>> print ds.smodel['scale']
            1.0
            >>> type(ds.dmodel) == list and len(ds.dmodel) == 1
            True
            >>> dm=ds.dmodel[0]
            >>> print type(dm)
            <class 'skymaps.DiffuseFunction'>

        Example loading an isotropic powerLaw source. This code is a klugey, and could be improved

        (XML from http://www.slac.stanford.edu/exp/glast/wb/prod/pages/sciTools_binnedLikelihoodTutorial/binnedLikelihood_v01archived.htm)

            >>> ds=l('''
            ...   <source name="Extragalactic Diffuse" type="DiffuseSource">
            ...     <spectrum type="PowerLaw">
            ...       <parameter free="1" max="100.0" min="1e-05" name="Prefactor" scale="1e-07" value="1.6"/>
            ...       <parameter free="0" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>
            ...       <parameter free="0" max="200.0" min="50.0" name="Scale" scale="1.0" value="100.0"/>
            ...     </spectrum>
            ...     <spatialModel type="ConstantValue">
            ...       <parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>
            ...     </spatialModel>
            ...   </source>''')
            >>> print ds.smodel.name
            PowerLaw
            >>> print ds.smodel['norm']
            1.0
            >>> print ds.smodel['index']
            1.0
            >>> type(ds.dmodel) == list and len(ds.dmodel) == 1
            True
            >>> dm=ds.dmodel[0]
            >>> print type(dm)
            <class 'skymaps.IsotropicPowerLaw'>
            >>> print dm.index()
            2.1
            >>> np.allclose(PowerLaw(norm=1.6e-7, index=2.1, e0=1e2).i_flux(1e2, np.inf), dm.flux())
            True

        Example loading an isotropic source from a file
    
        (xml from http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/likelihood_tutorial.html)

            >>> ds=l('''
            ...   <source name="iso_p7v6source" type="DiffuseSource">
            ...     <spectrum file="$GLAST_EXT/diffuseModels/v2r0p1/isotrop_2year_P76_source_v1.txt" type="FileFunction">
            ...       <parameter free="1" max="1000" min="1e-05" name="Normalization" scale="1" value="1" />
            ...     </spectrum>
            ...     <spatialModel type="ConstantValue">
            ...       <parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>
            ...     </spatialModel>
            ...   </source>''')
            >>> print ds.smodel.name
            Constant
            >>> print ds.smodel['scale']
            1.0
            >>> type(ds.dmodel) == list and len(ds.dmodel) == 1
            True
            >>> dm=ds.dmodel[0]
            >>> print type(dm)
            <class 'skymaps.IsotropicSpectrum'>

            
        Now, we are going to test loading extended sources:

        First, create the profile from a simple Gaussian

            >>> from tempfile import NamedTemporaryFile
            >>> gaussian = Gaussian(center=SkyDir(82.73, 13.38))

        First, create & load a new Disk extended soure

            >>> ds=l('''
            ...   <source name="test_map" type="DiffuseSource">
            ...      <spectrum type="ConstantValue">
            ...        <parameter free="1" max="1000.0" min="0.001" name="Value" scale="1" value="1.0"/>
            ...      </spectrum>
            ...      <spatialModel type="Gaussian">
            ...        <parameter free="0" name="RA" scale="1.0" value="%s"/>
            ...        <parameter free="1" name="DEC" scale="1.0" value="%s"/>
            ...        <parameter free="1" name="Sigma" scale="1.0" value="%s"/>
            ...      </spatialModel>
            ...    </source>''' % (gaussian['ra'], gaussian['dec'],gaussian['sigma']))
            >>> new_gauss =  ds.spatial_model
            >>> np.allclose(new_gauss.p, gaussian.p)
            True
            >>> print new_gauss.free
            [False  True  True]

        Now, create & load an xml file for a SpatialMap

            >>> temp = NamedTemporaryFile(delete=True)
            >>> gaussian.save_template(temp.name)

            >>> ds=l('''
            ...   <source name="test_map" type="DiffuseSource">
            ...      <spectrum type="ConstantValue">
            ...        <parameter free="1" name="Value" scale="1" value="1.0"/>
            ...      </spectrum>
            ...      <spatialModel type="SpatialMap" file="%s">
            ...        <parameter free="0" name="Prefactor" scale="1.0" value="1"/>
            ...      </spatialModel>
            ...    </source>''' % temp.name)
            >>> print ds.model['scale']
            1.0
            >>> map=ds.spatial_model
            >>> print map.name
            SpatialMap
            >>> np.allclose([map.center.ra(), map.center.dec()], [gaussian.center.ra(), gaussian.center.dec()], atol=1e-2, rtol=1e-2)
            True
            >>> np.allclose(map(map.center),gaussian(map.center), atol=1e-2, rtol=1e-2)
            True

        Now, create & load a an xml file for a RadialProfile

            >>> temp = NamedTemporaryFile(delete=True)
            >>> gaussian.save_profile(temp.name, numpoints=1000)

            >>> ds=l('''
            ...   <source name="test_profile" type="DiffuseSource">
            ...      <spectrum type="ConstantValue">
            ...        <parameter free="1" name="Value" scale="1" value="1.0"/>
            ...      </spectrum>
            ...      <spatialModel type="RadialProfile" file="%s">
            ...        <parameter free="0" name="Normalization" scale="1.0" value="1"/>
            ...        <parameter free="0" name="RA" scale="1.0" value="%s"/>
            ...        <parameter free="0" name="DEC" scale="1.0" value="%s"/>
            ...      </spatialModel>
            ...    </source>''' % (temp.name, gaussian.center.ra(), gaussian.center.dec()))
            >>> profile=ds.spatial_model
            >>> print profile.name
            RadialProfile
            >>> np.allclose([profile.center.ra(), profile.center.dec()], [gaussian.center.ra(), gaussian.center.dec()])
            True
            >>> np.allclose(profile(profile.center),gaussian(profile.center), atol=1e-5, rtol=1e-5)
            True

        You should also be able to load in extended soruces in galactic coordinates

            >>> ds=l('''
            ...   <source name="test_map" type="DiffuseSource">
            ...      <spectrum type="ConstantValue">
            ...        <parameter free="1" max="1000.0" min="0.001" name="Value" scale="1" value="1.0"/>
            ...      </spectrum>
            ...      <spatialModel type="Gaussian">
            ...        <parameter free="0" name="L" scale="1.0" value="0"/>
            ...        <parameter free="1" name="B" scale="1.0" value="0"/>
            ...        <parameter free="1" name="Sigma" scale="1.0" value="1"/>
            ...      </spatialModel>
            ...    </source>''')
            >>> ds.spatial_model['l']
            0.0
            >>> ds.spatial_model['b']
            0.0
            >>> ds.spatial_model.coordsystem == SkyDir.GALACTIC
            True
    """
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

            if spatial['type'] in XML_to_SpatialModel.spatialdict.keys():

                spatial_model=XML_to_SpatialModel.get_spatial_model(spatial,diffdir=diffdir)
                spectral_model=xtm.get_model(spectral,name)
                ds.append(ExtendedSource(name=name,
                                         model=spectral_model,
                                         spatial_model=spatial_model))
            else:
                raise Exception('Diffuse spatial model "%s" not recognized' % spatial['type'])
    return list(ds)

def parse_sources(xmlfile,diffdir=None,roi_dir=None,max_roi=None):
    """ Convenience function to parse an entire file into
        point sources and diffuse sources. """
    handler = parse_sourcelib(xmlfile)
    ps = parse_point_sources(handler,roi_dir,max_roi)
    ds = parse_diffuse_sources(handler,diffdir=diffdir)
    return ps,ds

def unparse_point_sources(point_sources, strict=False, properties=lambda x:''):
    """ Convert a list (or other iterable) of PointSource objects into XML.
        strict : bool
            set True to generate exception, error message identifying offending source, reason
        properties : a function
            the function, if specified, returns a string for the source element with properties, like TS

            >>> ps = PointSource(name='test', model=Constant(), skydir=SkyDir(-30,30))
            >>> ret=unparse_point_sources([ps])
            >>> print len(ret)
            1
            >>> print ret[0].strip().replace('\\t', ' '*4)
            <source name="test" type="PointSource"  >
                <spectrum   type="ConstantValue">
                    <parameter name="Value" value="1.0" free="1" max="10" min="0.001" scale="1" />
                </spectrum>
                <spatialModel type="SkyDirFunction">
                    <parameter name="RA"  value="330.0" free="0" max="360.0" min="-360.0" scale="1.0" />
                    <parameter name="DEC" value="30.0" free="0" max="90" min="-90" scale="1.0" />
                </spatialModel>
            </source>
    """
    xml_blurbs = Stack()
    m2x = Model_to_XML(strict=strict)
    for ps in point_sources:
        skyxml = makePSSpatialModel(ps.skydir)
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
    """ Convert an instance of DiffuseSource into an XML blurb.
        
        Some simple testing of saving out diffuse sources:

            >>> from uw.like.pointspec_helpers import get_default_diffuse
            >>> from os.path import expandvars
            >>> diffdir=expandvars('$GLAST_EXT/diffuseModels/v2r0p1/')
            >>> ds = get_default_diffuse(diffdir=diffdir,
            ...     gfile="ring_2year_P76_v0.fits",
            ...     ifile="isotrop_2year_P76_source_v1.txt")
            >>> p = lambda d: process_diffuse_source(d).replace(diffdir,'').replace('\\t',' '*4).strip()
            >>> print p(ds[0])
            <source name="Galactic Diffuse (ring_2year_P76_v0.fits)" type="DiffuseSource">
                <spectrum   type="PowerLaw">
                    <parameter name="Prefactor" value="1.0" free="1" max="10" min="0.1" scale="1" />
                    <parameter name="Index" value="0.0" free="1" max="1" min="-1" scale="-1" />
                    <parameter name="Scale" value="1000.0" free="0" max="1000.0" min="1000.0" scale="1" />
                </spectrum>
                <spatialModel file="ring_2year_P76_v0.fits" type="MapCubeFunction">
                    <parameter name="Normalization" value="1.0" free="0" max="1e3" min="1e-3" scale="1.0" />
                </spatialModel>
            </source>
            >>> print p(ds[1])
            <source name="Isotropic Diffuse (isotrop_2year_P76_source_v1.txt)" type="DiffuseSource">
                <spectrum file="isotrop_2year_P76_source_v1.txt" ctype="-1" type="FileFunction">
                    <parameter name="Normalization" value="1.0" free="1" max="10" min="0.1" scale="1" />
                </spectrum>
                <spatialModel type="ConstantValue">
                    <parameter  name="Value" value="1.0" free="0" max="10.0" min="0.0" scale="1.0" />
                </spatialModel>
            </source>
    """
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
    elif isinstance(dm,IsotropicConstant):

        value = dm.constant()
        skyxml = makeDSConstantSpatialModel(value)

        m2x.process_model(ds.smodel,scaling=False)
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
                flux,index=m.flux(),m.index()
                pl=PowerLawFlux(index=index)
                pl.set_flux(flux,100,np.inf)

                if isinstance(ds.smodel,Constant):
                    pl['Int_Flux'] *= ds.smodel['Scale']
                else:
                    raise Exception("...")
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

def write_sources(point_sources, diffuse_sources, filename, strict=False,convert_extended=False,expand_env_vars=False):
    source_xml = [unparse_point_sources(point_sources, strict=strict)]
    if len(diffuse_sources)>0:
        source_xml.append(unparse_diffuse_sources(diffuse_sources,
                                                  convert_extended=convert_extended,
                                                  expand_env_vars=expand_env_vars,
                                                  filename=filename))
    writeXML(source_xml,filename)

def writeROI(roi,*args, **kwargs):
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
    write_sources(roi.psm.point_sources, roi.dsm.diffuse_sources, *args, **kwargs)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
