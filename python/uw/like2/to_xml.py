"""
Generate the XML representation of a list of sources
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/to_xml.py,v 1.20 2017/08/02 23:07:03 burnett Exp $

"""
import os, collections, argparse, types, glob
from astropy.io import fits as  pyfits
import numpy as np
import pandas as pd
from uw.like import Models
from skymaps import SkyDir
from . import ( sources, extended)
from uw.utilities import xml_parsers
from collections import OrderedDict

class Element(object):
    """ Manage creation of an  XML document, containing nested elements
    """
    level=0
    stream=None
    def __init__(self, element_name,  ignore=(), **kw):
        """
        To write to a file, set Element.stream to the open file. Otherwise to standard out
        
        Example:
        
        with Element('A', Aprop=2, b=3) as t:
            t.text('A string')
            with Element('B', u=99) as s:
                s.text('test string')
                SimpleElement('C', ctest=9):
        
        """
        # if there is a name and/or type in the keywords, put them first for readability
        items = []
        for first in ('title', 'name', 'type'):
            t = kw.pop(first, None)
            if t is not None:
                items.append((first, t))
        self.kw=OrderedDict( items + [k for k in kw.items() if (k[0] not in ignore) and (k[0][0]!='_') ])
        self.name =element_name
        self.finish()
        
    def finish(self): # for SimpleElement to override
        pass
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
    def finish(self):
        self.text('<%s %s/>'% (self.name, self))


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
            # if index2<0: 
            #     print  'Source %s has beta (%.2f) <0: setting to 0.' % (  source.name, index2, )
            #     index2=0
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
        

def from_roi(roimodel, title=None, stream=None, strict=True, maxi=None):
    """
    Create an XML file describing the complete data model for an ROI
    
    roimodel : roimodel.ROImodel object
        Expect to be a list of sources.Source objects
    """
    Element.stream = stream
    m2x = xml_parsers.Model_to_XML(strict=strict)
    if title is None:
        title = 'Sources from roi %04d' % roimodel.index
        
    with Element('source_library', title=title,  
        ignore=('selected_source','quiet'), **roimodel.__dict__) as sl:
    
        for source in roimodel:
            stype = source.__class__.__name__

            with Element('source',  type=stype, ignore=('model','sedrec', 'free', 'changed', 'spatial_model'),
                        **source.__dict__) as src:
                # insert model, a spectrum element. Uses class from xml_parsers
                m2x.process_model(source.model)
                src.text(m2x.getXML(tablevel=0))
                
                # insert spatial model depending on source type
                if stype=='PointSource':
                    src.text(xml_parsers.makePSSpatialModel(source.skydir, tablevel=0))
                elif stype=='ExtendedSource':
                    with Element('spatialModel', type=source.dmodel.name, 
                            file="%s"%source.dmodel.file) as sm:
                        SimpleElement('parameter', name='Prefactor', value=1.0, free=0, max=1e3,min=1e-3, scale=1.0)
                elif stype=='GlobalSource':
                    SimpleElement('spatialModel', type=source.dmodel[0].__class__.__name__,
                        ignore=('spectral_function','files', 'opener','limbfun','energy', 'loaded'),
                            **source.dmodel[0].__dict__)
            if maxi is not None and i>maxi: break


def source_library(source_list, title='sources', stream=None, strict=False, maxi=None):
    """ Generate sources in XML format 
    
    source_list : pandas.DataFrame
    stream : output stream
        if None, to console
    
    """
    Element.stream = stream
    m2x = xml_parsers.Model_to_XML(strict=True)
    ns=ne=0
    try:
        extended = pd.DataFrame(pyfits.open(glob.glob(
            os.path.expandvars('$FERMI/catalog/Extended_archive*/LAT_extended_sources*.fit*'))[-1]
            )[1].data)
        extended.index= [x.strip() for x in extended['Source_Name']]
    except Exception, msg:
        raise Exception('Failed to find the Extended archive: %s' %msg)

    with Element('source_library', title=title) as sl:
        for i,source in source_list.iterrows():
            stype = 'ExtendedSource' if np.isnan(source['locqual']) else 'PointSource'
            with Element('source', type=stype,  ignore=('model','sedrec', 'free'),
                        **source) as src:
                m2x = xml_parsers.Model_to_XML(strict=strict)
                m2x.process_model(pmodel(source))
                src.text(m2x.getXML(tablevel=0))
                if stype=='PointSource':
                    src.text(xml_parsers.makePSSpatialModel(SkyDir(source['ra'],source['dec']),tablevel=0))
                    ns +=1
                else:
                    sname = source['name']
                    if sname not in extended.index:
                        print 'Failed to find %s, with no localizaiton, in extended sources' %sname
                        print 'Skipping ...'
                        continue
                    with Element('spatialModel', type='SpatialMap', 
                            file="%s" % extended.ix[sname]['Spatial_Filename'].strip() ) as sm:
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
