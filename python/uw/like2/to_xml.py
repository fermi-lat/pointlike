"""
Generate the XML representation of a skymodel
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/to_xml.py,v 1.10 2013/07/09 16:05:45 burnett Exp $

"""
import os, collections, argparse, types, glob, pyfits
import numpy as np
import pandas as pd
from uw.like2 import skymodel
from uw.utilities import keyword_options
from uw.like import Models
from skymaps import SkyDir

from uw.utilities import xml_parsers
from collections import OrderedDict

class Element(object):
    """ Manage creation of an  XML document, containing nested elements
    """
    level=0
    stream=None
    def __init__(self, element_name,  **kw):
        """
        To write to a file, set Element.stream to the open file. Otherwise to standard out
        
        Example:
        
        with Element('A', Aprop=2, b=3) as t:
            t.text('A string')
            with Element('B', u=99) as s:
                s.text('test string')
                SimpleElement('C', ctest=9):
        
        """
        # if there is a name in the keywords, put if first for readability
        name = kw.pop('name', None)
        self.kw=OrderedDict( kw.items() if name is None else [('name',name)]+kw.items() )
        self.name =element_name
    def __str__(self):
        return ' '.join('%s="%s"'%item  for item in self.kw.items()) 
    def text(self, text):
        offset = '  '*Element.level
        for line in text.split('\n'):
            if line=='': continue
            print >> self.stream , offset + line
    def __enter__(self):
        self.text('<%s %s>'% (self.name, self))
        Element.level += 1
        return self
    def __exit__(self, type, value, traceback):
        Element.level -= 1
        self.text('</%s>'%self.name )

class SimpleElement(Element):
    """ An element not containing text."""
    def __init__(self, element_name, **kw):
        self.kw=OrderedDict(**kw)
        self.text('<%s %s/>'% (element_name, self))


def pmodel(source):
    """ create a pointlike model from a Series object from DataFrame row
    """
    expandbits = lambda x : [(x & 2**i)!=0 for i in range(4)]
    
    modelname, freebits, e0, norm, norm_unc, pindex, pindex_unc, index2,index2_unc, cutoff, cutoff_unc=\
        [source[x] for x in 'modelname freebits e0 flux flux_unc pindex pindex_unc index2 index2_unc cutoff cutoff_unc'.split()]
    free = expandbits(freebits)
    errors = [norm_unc, pindex_unc]
    if modelname=='LogParabola':
        if np.abs(index2)<2e-3:
            modelname='PowerLaw'
            model = Models.PowerLaw(p=[norm, pindex ], e0=e0)
        else:
            if index2<0: 
                print  'Source %s has beta (%.2f) <0: setting to 0.' % (  source.name, index2, )
                index2=0
            model =Models.LogParabola(p= [norm, pindex, index2, e0])
            model.free[-1]=False
            errors.append(index2_unc)
    elif modelname=='PLSuperExpCutoff':
        prefactor = np.exp(-(e0/cutoff)**index2)
        model = Models.PLSuperExpCutoff(p = [norm/prefactor, pindex, cutoff, index2], e0=e0)
        errors[0] /= prefactor
        errors += [cutoff_unc, index2_unc]
    elif modelname=='ExpCutoff':
        prefactor = np.exp(-(e0/cutoff))
        model = Models.PLSuperExpCutoff(p = [norm/prefactor, pindex, cutoff, 0.], e0=e0)
        model.free[3]==False
        errors[0] /= prefactor
        errors += [cutoff_unc]
    else:
        raise Exception('model "%s" not recognized' % modelname)
    map(model.set_error, range(len(errors)), errors)
    model.free[:]=free[:model.len()]    
    return model
        

def source_library(source_list, title='sources', stream=None, strict=False, maxi=None):
    """ Generate sources in XML format 
    source_list : pandas.DataFrame
    
    """
    Element.stream = stream
    m2x = xml_parsers.Model_to_XML(strict=True)
    ns=ne=0
    try:
        extended = pd.DataFrame(pyfits.open(glob.glob(
            os.path.expandvars('$FERMI/catalog/Extended_archive*/LAT_extended_sources*.fit'))[-1]
            )[1].data)
        extended.index= [x.strip() for x in extended['Source_Name']]
    except Exception, msg:
        raise Exception('Failed to find the Extended archive: %s' %msg)

    with Element('source_library', title=title) as sl:
        for i,source in source_list.iterrows():
            stype = 'DiffuseSource' if np.isnan(source['locqual']) else 'PointSource'
            with Element('source', type=stype, **source) as src:
                m2x = xml_parsers.Model_to_XML(strict=strict)
                m2x.process_model(pmodel(source))
                src.text(m2x.getXML(tablevel=0))
                if stype=='PointSource':
                    src.text(xml_parsers.makePSSpatialModel(SkyDir(source['ra'],source['dec']),tablevel=0))
                    ns +=1
                else:
                    with Element('spatialModel', type='SpatialMap', 
                            file=extended.ix[source['name']]['Spatial_Filename'].strip() ) as sm:
                        SimpleElement('parameter', name='Prefactor', value=1.0, free=0, max=1e3,min=1e-3, scale=1.0)
                    ne += 1
            if maxi is not None and i>maxi: break
    return ns,ne

def main( filename=[], sources='sources*.csv', cuts='(sources.ts>10)' ):
    t = sorted(glob.glob(sources))[-1] #shold get the one we want
    sources = pd.read_csv(t)
    print 'read %d sources from %s' %(len(sources), t)
    cut_sources = sources[eval(cuts)]
    print 'applied cut "%s", %d remain' % (cuts, len(cut_sources))
    modelname = '_'.join(os.path.abspath('.').split('/')[-2:])
    filename = filename[0] if len(filename)>0 else None
    if filename is None:
        filename = modelname+'.xml'
        # for example, 'P202_uw10.xml'
    with open(filename,'w') as stream:
        ns,ne= source_library(cut_sources, title=modelname, stream=stream)
    print 'wrote file %d point sources, and %d extended sources to file %s' % (ns,ne,filename)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser( description=""" Convert the skymodel in the current folder to XML""")
    parser.add_argument('filename', nargs='*', help='filename to write to (default: make it up)')

    parser.add_argument('--sources', default='sources*.csv', help='input table')
    parser.add_argument('--cuts',  default='(sources.ts>10)',
            help='selection cuts')
    args = parser.parse_args()
    main(args.filename, args.sources, args.cuts)
