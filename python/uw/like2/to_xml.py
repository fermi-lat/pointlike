"""
Generate the XML representation of a skymodel
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/to_xml.py,v 1.2 2013/03/11 18:20:39 burnett Exp $

"""
import os, collections, argparse, types
from uw.like2 import skymodel
from uw.utilities import keyword_options

from uw.utilities import xml_parsers

class XMLelement(object):
    level=0
    stream=None
    def __init__(self, name, **kw):
        self.name, self.kw = name, kw
    def __str__(self):
        return ' '.join('%s="%s"'%item  for item in self.kw.items()) 
    def text(self, text):
        print >> self.stream ,'  '*XMLelement.level + text
    def __enter__(self):
        self.text('<%s %s>'% (self.name, self))
        XMLelement.level += 1
        return self
    def __exit__(self, type, value, traceback):
        XMLelement.level -= 1
        self.text('</%s>'%self.name )

class ToXML(object):

    def __init__(self, skymodel, filename, ts_min=10, a_max=0.25, title=None, source_filter=lambda x:True, strict=False, gtlike = False):
        """ generate a file with the XML version of the sources in the model
        parameters
        ----------
        skymodel : SkyModel object to convert
        filename : string or None
            name of file to write to; if None make it up from local folder
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
        if filename is None:
            filename = '_'.join(os.path.abspath('.').split('/')[-2:])+'.xml'
            # for example, 'P202_uw10.xml'
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
                f.write('\n'.join(['\n<roi_info nside="12">',
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


def main( args ):
    sm = skymodel.SkyModel('.')
    filename = args.filename[0] if len(args.filename)>0 else None
    ToXML(sm, filename, ts_min=args.ts_min, a_max=args.a_max, strict=args.strict, gtlike=args.gtlike)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser( description=""" Convert the skymodel in the current folder to XML""")
    parser.add_argument('filename', nargs='*', help='filename to write to (default: make it up)')
    parser.add_argument('--cuts',  default='(sources.ts>10)*(sources.a<0.25)*(sources.locqual<10)', help='selection cuts')
    parser.add_argument('--ts_min', default=16, help='minimum TS')
    parser.add_argument('--a_max', default=0.25, help='maximum major radius')
    parser.add_argument('--strict', default=False, action='store_true', help='apply strict' )
    parser.add_argument('--gtlike', default=False, action='store_true', help='make compatible with gtlike' )
    args = parser.parse_args()
    main(args)
