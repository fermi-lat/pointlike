"""
generate collection file for LiveLabs Pivot viewer

$Header$
See <http://getpivot.com>
Author: Toby Burnett <tburnett@uw.edu>
"""

import os, pickle, pyfits
import numpy as np

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


    def __init__(self, name, item_names, dzc, img_list=None, icon=None, href=None):
        """
        parameters description
        ---------- -----------
        name       name for the collection
        item_names list of names for each item
        dzc        filename for the Deep Zoom (DZ) Collection of images, relative to base
        img_list  if specified, a list of integers for the images in the DZ collection corresponding to the names
        href      list of references
       
        """
        self.name = name
        self.facets =[]
        self.item_names = item_names
        self.n = len(item_names)
        assert(self.n>0)
        self.dzc = dzc
        if img_list is not None:
            self.img_list = img_list
            assert(len(img_list)==self.n)
        else:
            self.img_list = range(self.n) 
        self.icon = icon
        self.href=href

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
        data       array of values
        
        """
        assert( type in 'String Number LongString DateTime Link'.split())
        assert(len(data)==self.n) # make sure all columns the same length
        self.facets.append(Collection.Facet(name, type, format, data, filter))
        
    def write(self, outfile):
        assert(len(self.facets)>0) # fail if no Facets were added
        out= open(outfile, 'w')
        out.write(Collection.header %  (self.name, ('' if self.icon is None else ' d1p1:Icon="%s"'%self.icon) ))

        out.write('\n<FacetCategories>' )
        for facet in self.facets:
            #IsFilterVisible="%(filter)s"
            out.write('\n  <FacetCategory Name="%(name)s" Type="%(type)s" Format="%(format)s" d1p1:IsFilterVisible="%(filter)s"/>' %(facet.__dict__))
        out.write('\n</FacetCategories>')
        
        out.write('\n<Items ImgBase="%s">' %  self.dzc)
        for i,name in enumerate(self.item_names):
            img = self.img_list[i]
            href= ' Href="%s"' % self.href[i] if self.href is not None else ""
            out.write('\n<Item Id="%d" Img="#%d" Name="%s" %s>' % (i, img, name, href))
            out.write('\n <Facets>')
            for facet in self.facets:
                if facet.type != 'Link':
                    out.write('\n  <Facet Name="%s"> <%s Value="%s"/> </Facet>' % (facet.name, facet.type, facet.data[i]))
                else:
                    href = facet.data[i]
                    out.write('\n  <Facet Name="%s"> <Link Name="%s" Href="%s"/> </Facet>' % (facet.name,  href,href))
                
            out.write('\n </Facets>\n</Item>')    
        out.write('\n</Items>')
        
        out.write(Collection.trailer)
        out.close()
        
if __name__=='__main__':
    pass
    


