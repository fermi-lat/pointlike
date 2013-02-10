"""
Generate the XML representation of a skymodel
$Header$

"""
import os, collections, argparse
from xml import sax
from . import skymodel
from uw.utilities import keyword_options

from uw.utilities import xml_parsers

class ToXML(object):

    def __init__(self, skymodel, filename, ts_min=None, title=None, source_filter=lambda x:True, strict=False, gtlike = False):
        """ generate a file with the XML version of the sources in the model
        parameters
        ----------
        skymodel : SkyModel object to convert
        filename : string
            name of file to write to
            ts_min : None or float
                if set, only select sources with ts>ts_min
            title : string
                set title property of the source_library
            source_filter : function
                if set, function of a source that returns bool
            strict : bool
                set to True to apply strict rules
            gtlike :bool
                set to True to generate only a list of sources, no ROI information
        """
        self.skymodel = skymodel
        point_sources = self.skymodel.point_sources if ts_min is None\
            else filter(lambda s: s.ts>ts_min, self.skymodel.point_sources)
        print 'SkyModel: writing XML representations of %d point sources %s and %d extended sources to %s' \
            %(len(point_sources), ('' if ts_min is None else '(with TS>%.f)'%ts_min), len(self.skymodel.extended_sources), filename)
        def pointsource_properties(s):
            if hasattr(s.model,'e0'): e0 = s.model.e0 # don't think this is needed now
            else: e0 = 10**s.model.getp(3)
            return 'Pivot_Energy="%.1f" TS="%.1f"' % (e0, s.ts)
        stacks= [
            xml_parsers.unparse_diffuse_sources(self.skymodel.extended_sources,convert_extended=True,filename=filename),
            xml_parsers.unparse_point_sources(point_sources,strict=strict, properties=pointsource_properties),
        ]
        gs_xml = self._global_sources_to_xml(filename)
        with open(filename,'wb') as f:
            if not gtlike:
                f.write('<skymodel>\n')
            f.write('<source_library title="%s">'% title)
            for stack in stacks:
                for elem in stack:
                    f.write(elem)
            f.write('\n</source_library>')
            if not gtlike:
                f.write('\n'.join(['\n<roi_info nside="{0}">'.format(self.nside),
                                   gs_xml,
                                   '</roi_info>']))
                f.write('\n</skymodel>')

    def _global_sources_to_xml(self,filename):
        stacks = []
        bad =0
        for i in xrange(1728):
            stack = xml_parsers.Stack()
            s1 = '<roi index="{0}">'.format(i)
            s2 = '</roi>'
            globals = self.skymodel.global_sources[i]
            for s in globals:
                prefix = s.name.split('_')[0]
                s.name, s.dmodel = prefix, self.skymodel.diffuse_dict[prefix]
                s.smodel = s.model
            try:
                diffuse_xml = xml_parsers.unparse_diffuse_sources(globals,filename=filename)
            except:
                bad +=1
                continue
            for x in diffuse_xml:
                x = '\t'+x
            diffuse_xml.appendleft(s1)
            diffuse_xml.append(s2)
            stacks+=['\n'.join(diffuse_xml)]
        if bad>0:
            print 'Failed to convert %d ROIs' % bad
        return '\n'.join(stacks)


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

def main( args ):
    sm = skymodel.SkyModel('.')
    ToXML(sm, args.filename, ts_min=args.ts_min, strict=args.strict, gtlike=args.gtlike)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser( description=""" Convert the skymodel in the current folder to XML""")
    parser.add_argument('filename', nargs=1, help='filename to write to')
    parser.add_argument('--ts_min',   default=10, help='minimum TS')
    parser.add_argument('--strict', action='store_true', help='apply strict' )
    parser.add_argument('--gtlike', action='store_true', help='make compatible with gtlike' )
    args = parser.parse_args()
    main(args)
