"""
HTML clickable map stuff
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/tsmap_grid.py,v 1.4 2010/06/20 14:16:30 burnett Exp $

author: Toby Burnett <tburnett@uw.edu>
"""
from PIL import Image
import glob, os, sys

def rename(indir='tsmap'):
    """ fix filenames with + """
    filenames = glob.glob(os.path.join(indir, '*_tsmap.png'))
    for file in filenames:
        if '+' in file:
            cmd = 'ren "%s" "%s"' % (file, (os.path.split(file)[1]).replace('+','p'))
            print cmd,  os.system(cmd)
 
html_header="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
  <title>%s</title>
</head>
<body><div align="center"><h2>%s</h2></div>\n"""
html_trailer="""\n</body></html>\n"""

class Holder(object):
    def __init__(self, size=(1536,1280)):
        self.image = Image.new('RGB', size, (255,255,255)) 
    def insert(self, filename, size=None, pos=(0,0)):
        """
        """
        im = im2= Image.open(filename)
        if size is not None:
            im2 = im.resize(size)
        self.image.paste(im2, pos)
        return im2 
        
    def save(self, filename, crop = None):
        """ crop is a box
        """
        if crop is not None:
            self.image.crop(crop).save(filename)
        else:
            self.image.save(filename)
    
        

def make_map(images, names=None,
            title = 'TS map',
            html_file='tsmap_grid.htm',  img_file=None, 
            size=50, columns=None, 
            map_name=None,
            border=1,
            html_page=True,
            ):
    """ create a grid of thumbnails from the list of images
        and a corresponding html map
        
        images: list of image file names
        title: title for web page, and heading for the image
        names : [None] corresponding titles for tool tips (default: the filenames)
        title:  ['TSmap'] title to use for the web page.
        html_file: name of the HTML file to generate or append to (if mtml_page is False)
        img_file: [None] if None, use the same name as the html_file. (must specify for multiple maps)
        size:   [50] square thumnail size
        columns: [None] number of columns, default 1000/size, or 20 if size=50
        html_page: [True] generate entire HTML page, with proper head. If False, will append the map to a presumably existing file
    """
    map_header ="""\n<div align="center">\n<map name="%s"><h3>%s</h3>"""
    map_entry  ="""\n<area href="%s" title="%s" coords="%d,%d,%d,%d"/>"""
    map_trailer="""\n</map>\n<img alt="map grid" src="%s" border="%d" width="%d" height="%d" usemap="#%s"/>\n</div>\n"""
    
    n = len(images)
    assert(n>0)
    if names is None: 
        names = images
    else:
        assert(n==len(names))
    if columns is None: columns = 1000//size 
    if img_file is None: img_file=html_file.replace('.htm','.jpg')
    if map_name is None: map_name = img_file.split('.')[0]
    rows = n//columns+1
    width,height = size*columns, size*rows
    thumbnails = Image.new('RGB', (width,height), (255,255,255)) 
    out=open(html_file, 'w' if html_page else 'a')
    if html_page: out.write(html_header%(title,title))
    out.write(map_header% (map_name,title))
    for i, pair in enumerate(zip(images, names)):
        img,name = pair
        left,upper = i%columns*size, i//columns*size 
        # open the image in PIL, convert to a thumbnail, and insert into the map
        im = Image.open(img)
        im.thumbnail((size,size), Image.ANTIALIAS)
        thumbnails.paste(im, (left,upper))
        # add line with area setup
        out.write(map_entry % (img.replace('\\','/'), name, left,upper,left+size,upper+size) ) 
    thumbnails.save(img_file) 
    out.write(map_trailer% (img_file,border,width,height,map_name)) 
    if html_page: out.write(html_trailer)
    out.close()

def main(image_dir='tsmap', title='1FGL TS maps from 11 month data set', html_file='tsmap_grid.htm', **kwargs):
    """create a grid of thumbnails from the files in image_dir, in alphabetical order
    """
    assert(os.path.exists(image_dir))
    images = glob.glob(os.path.join(image_dir,'*.png'))
    assert(len(images)>0)
    images.sort()
    names = ['%d: %s' %(i+1,(os.path.split(img)[1][:-3]).replace('p','+')) for (i,img) in enumerate(images)]
    make_map(images, names,title=title, html_file=html_file, **kwargs)
    
if __name__=='__main__':
    main()
     
  
    