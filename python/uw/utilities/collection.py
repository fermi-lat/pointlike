"""
generate collection file for the  Pivot viewer

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/collection.py,v 1.9 2010/12/09 23:31:50 burnett Exp $
See <http://getpivot.com>
Author: Toby Burnett <tburnett@uw.edu>
"""

import os, types, exceptions
import xml.sax

import numpy as np
class InvalidParameter(exceptions.Exception):
    pass

def get_image_ids( dzc):
    """ look for the names in the DZC file 
    
    """
    class Handler(xml.sax.ContentHandler):
        def __init__(self):
            self.id_dict={}
        def startElement(self, name, attrs):
            if name == "I":
                id = attrs['Id']
                source =attrs['Source']
                fname =os.path.split(source)[1]
                self.id_dict[fname]=id
    parser = xml.sax.make_parser()
    h=Handler()
    parser.setContentHandler(h)
    parser.parse(open(dzc,'r'))
    return h.id_dict

def image_list(item_names, full_dzc):
    """
     search the DZC for names that correspond to the names in item_names
    # the names are source names, so file names have had blanks and + signs replaced
    """
    if not os.path.exists(full_dzc):
        raise InvalidParameter('DZC file %s does not exist' % full_zdc)
    ids = get_image_ids(full_dzc) # this is a dictionary, key the source name, value the index of the dzi
    # 
    names = [n.replace('+', 'p').replace(' ','_')  for n in item_names]
    try:
        return [ids[n+'.xml'] for n in names]
    except:
        for n in names:
            if n+'.xml' not in ids:
                print ('Item %s not found in list of images %s' %  ( n, full_dzc))
        raise


class Collection(object):
    """ manage a collection object.
    
    """
    header="""<?xml version="1.0"?>\n<Collection 
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
    xmlns:xsd="http://www.w3.org/2001/XMLSchema" 
    Name="%s" 
    SchemaVersion="1" 
    xmlns:d1p1="http://schemas.microsoft.com/livelabs/pivot/collection/2009" 
    xmlns="http://schemas.microsoft.com/collection/metadata/2009" %s>"""

    trailer = """\n</Collection> """ 


    def __init__(self, cname, folder, item_names, dzc,  icon=None, 
        href=None, 
        ignore_exception=True):
        """
        parameters description
        ---------- -----------
        cname       name for the collection
        folder     folder where files exist, the result will be put
        item_names list of names for each item
        dzc        filename for the Deep Zoom (DZ) Collection of images, relative to folder
        href       list of references, same length as item_names if specified
                    
        """
        self.name = cname
        self.folder = folder
        self.ignore_exception=ignore_exception
        if not os.path.exists(folder):
            raise InvalidParameter('folder %f does not exist' % folder)
        self.facets =[]
        self.item_names = item_names
        self.n = len(item_names)
        if self.n==0:
            raise InvalidParameter('List of Item names is empty')
        self.icon = icon
        self.href =href
        self.related = None
        self.dzc  = dzc
        

    class Facet(object):
        def __init__(self, name, type, format, data, filter):
            self.name,self.type,self.format, self.data, self.filter= name,type,format,data,filter
                
    def add_facet(self, name, type, format, data, filter='true'):
        """add a Facet to the collection.
        
        parameters description
        ---------- -----------
        name       facet name
        type       type, one of String, LongString, Number, DateTime, or Link
        format     format to use for info ('C', 'F', 'D')
        data       array of values. For String, each value can be an array itself
        
        """
        assert( type in 'String Number LongString DateTime Link'.split())
        assert(len(data)==self.n) # make sure all columns the same length
        self.facets.append(Collection.Facet(name, type, format, data, filter))
    
    def add_related(self, related):
        """
                related    None, or a list of n lists of pairs (name, href)
                    href is relative path to a related collection
        """

        self.related=related
        assert len(related)==self.n, 'related has wrong length'

    def write(self, outfile='pivot.cxml', last_id=-1):
        """ The Facet objects have all been added
            
            outfile ['pivot.cxml'], the full file name of the cxml, position relative to base
            last_id [-1]            used to sequence: add one for first element, etc.
        """
        if len(self.facets)==0: # fail if no Facets were added
            raise InvalidParameter('No Facets were added')

        out= open(outfile, 'w')
        out.write(Collection.header %  (self.name, ('' if self.icon is None else '\n    d1p1:Icon="%s"'%self.icon) ))

        self.writeFacetCategories(out)
        
        # write out one or more Items elements
        
        self.writeAllItems(out, last_id)
        
        out.write(Collection.trailer)
        out.close()
   
    def writeFacetCategories(self, out):
        out.write('\n<FacetCategories>' )
        for facet in self.facets:
            #IsFilterVisible="%(filter)s"
            out.write('\n  <FacetCategory Name="%(name)s" Type="%(type)s" Format="%(format)s" d1p1:IsFilterVisible="%(filter)s"/>'\
                %(facet.__dict__))
        out.write('\n</FacetCategories>')
    
    def writeAllItems(self, out,last_id):
        "write one or more Items elements: only one in this class"
        self.last_id = last_id
        img_list = image_list(self.item_names, os.path.join(self.folder,self.dzc))
        self.writeItems(out, self.dzc, self.item_names, img_list, self.facets)
  

    def write_related(self, out, related):
        out.write(
            '\n  <Extension><Related xmlns="http://schemas.microsoft.com/livelabs/pivot/collection/2009">')
        for name, href in related:
            out.write('\n   <Link Name="%s" Href="%s"/>' % (name,href) )
        out.write(    '\n  </Related></Extension>')
  
 
    def writeItems(self,out, imagebase, item_names, img_list, facets, select=None):
        """ write out an Items element
            out: open File object
            imgbase: the DZC file
            item_names: list of names for each item
            img_list: list of corresponding dzi indices
            facets: list of Facet objects with arrays of data (attributes name, data)
            select: if not None, an array of bools to select the set of Item's for this DZC
            
            """
        out.write('\n<Items ImgBase="%s">' %  imagebase)
        for i,name in enumerate(item_names):
            if select is not None and not select[i]: continue
            img = img_list[i]
            if img is None: continue # image was not found
            href= ' Href="%s"' % self.href[i] if self.href is not None else ""
            self.last_id+=1
            out.write('\n<Item Id="%d" Img="#%d" Name="%s" %s>' % (self.last_id, int(img), name, href))
            out.write('\n <Facets>')
            for facet in facets:
                if facet.type == 'Number':
                    datum = facet.data[i]
                    out.write('\n  <Facet Name="%s"> <%s Value="%.4f"/> </Facet>' % (facet.name, facet.type, datum))
                elif facet.type == 'Link':
                    t = facet.data[i] # has name, href
                    if t is not None:
                        out.write('\n  <Facet Name="%s"> <Link Name="%s" Href="%s"/> </Facet>' % (facet.name, t[0], t[1]))
                else:
                    out.write('\n  <Facet Name="%s"> '% facet.name)
                    d = facet.data[i]
                    if getattr(d,'__iter__',False):
                        for x in d: out.write('<%s Value="%s"/>' %( facet.type, x))
                    else: 
                        out.write('<%s Value="%s"/>' %( facet.type, d))
                    out.write(' </Facet>' ) 
                
            out.write('\n </Facets>')
            if self.related is not None:
                self.write_related(out, self.related[i])
            out.write('\n</Item>') 
        out.write('\n</Items>')
    
class MultiCollection(Collection):
    """ Subclass of Collection to support multiple DZCs
        difference is that the DZC argument is a list. 
        For each name in the item_names argument, it searches for the appropriate DZI, and organizes an Items element to 
        contain each set.
    
    """

    def _getid(self,name):
        fname = name.replace(' ','_').replace('+','p')+'.xml'
        for i, dict in enumerate(self.imageid_dicts):
            try:
                return i, dict[fname]
            except: 
                print ('name %s not found in DeepZoom collections'%fname)
                if self.ignore_exception: return i, None
        raise InvalidParameter('name %s not found in DeepZoom collections'%fname)
        
    def writeAllItems(self, out,last_id):
        """ override base class to manage multiple Items, according to where the individual DZI images are found
        """
        self.last_id = last_id
        self.imageid_dicts = [get_image_ids(os.path.join(self.folder,dzc)) for dzc in self.dzc]
        
        # make a list of pairs of the dzc index, and the dzi index for each name
        found = np.array([ self._getid(name) for name in self.item_names])
        #z = np.rec.fromarrays([found[:,0], found[:,1], self.item_names], names='dzc i name'.split())
        for i, dzc in enumerate(self.dzc):
            g = found[:,0]
            select = np.asarray([int(f)==i for f in g])
            if select.sum()==0:
                print ('warning: no Items with images in DZC %s' % dzc)
                continue
            self.writeItems(out, dzc, self.item_names, found[:,1], self.facets, select)
        
     
if __name__=='__main__':
    pass
    


