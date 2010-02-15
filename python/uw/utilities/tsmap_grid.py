"""
HTML clickable map stuff
$Header$

author: T Burnett <tburnett@uw.edu>
"""
from PIL import Image
import glob, os

def rename(indir='tsmap'):
    """ fix filenames with + """
    filenames = glob.glob(os.path.join(indir, '*_tsmap.png'))
    for file in filenames:
        if '+' in file:
            cmd = 'ren "%s" "%s"' % (file, (os.path.split(file)[1]).replace('+','p'))
            print cmd,  os.system(cmd)
    

  
def make_map(image_dir='tsmap', 
            img_file='tmap.jpg', 
            html_file='tmap.htm',
            size=50, columns=25, 
            map_name='tsmap_grid'
            ):
    """ create a grid of thumbnails from the files in image_dir, 
        and a corresponding html map
        assumes tsmap files, could be more general
    """
    map_header="""\n<div align="center">\n<map name="%s">"""
    map_entry="""\n<area href="%s" title="%s" coords="%d,%d,%d,%d">"""
    map_trailer="""\n</map>\n<img alt="grid" src="%s" border="0" width="%d" height="%d" usemap="#%s"/>\n</div>"""
    assert(os.path.exists(image_dir))
    images = glob.glob(os.path.join(image_dir,'*_tsmap.png'))
    images.sort()
    n = len(images)
    assert(n>0)
    rows = n//columns+1
    lx,ly = size*columns, size*rows
    cim = Image.new('RGB', (lx,ly), (255,255,255)) 
    out=open(html_file, 'w')
    out.write(map_header% map_name)
    for i, img in enumerate(images):
        # open the image in PIL
        im = Image.open(img)
        name = (os.path.split(img)[1][:-10]).replace('p','+')
        x,y = i%columns*size, i//columns*size
        # convert to a thumbnail and insert into the map
        im.thumbnail((size,size), Image.ANTIALIAS)
        cim.paste(im, (x,y))
        # add line with area setup
        out.write(map_entry % (img, '%d: %s'% (i+1,name), x,y,x+size,y+size) ) 
    cim.save(img_file)  #write out the combined thumbnail and html file
    out.write(map_trailer% (img_file,lx,ly,map_name)) 
    out.close()


    