""" 
Make combinded figures for Pivot, perhaps
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/makefig.py,v 1.5 2010/08/03 22:29:42 burnett Exp $

""" 

from PIL import Image
import glob, os, sys, exceptions
import pylab  as plt
version='$Revision: 1.5 $'.split()[1]

class InvalidParameter(exceptions.Exception):
    pass
    
def convert_log_to_png(outdir):
    """ convert txt files to png """
    path = os.path.join(outdir, 'log')
    files = glob.glob(os.path.join(path, '*.txt'))
    if len(files)==0: raise InvalidParameter('no .txt files found in folder "%s"' % path)
    plt.figure(figsize=(5,10))
    for file in files:
        plt.clf()
        t = open(file).read()
        plt.figtext(0.02,0.98,t,fontsize=8,va='top', fontname='monospace')
        plt.savefig(file.replace('.txt', '.png'))

def convert_one_log(filename):
    plt.figure(figsize=(5,10))
    plt.clf()
    t = open(filename).read()
    plt.figtext(0.02,0.98,t,fontsize=8,va='top', fontname='monospace')
    plt.savefig(filename.replace('.txt', '.png'))

class Holder(object):
    def __init__(self, size):
        self.image = Image.new('RGB', size, (255,255,255)) 
    def insert(self, filename, size=None, pos=(0,0)):
        """
        """
        im = im2= Image.open(filename)
        if size is not None:
            im2 = im.resize(size)
        self.image.paste(im2, pos)
        
    def save(self, filename, crop = None):
        """ crop is a box
        """
        if crop is not None:
            self.image.crop(crop).save(filename)
        else:
            self.image.save(filename)

def combine_images( sname, names, 
                    hsize=(2048,1600), 
                    sizes=     [(1050,1050) , None, None, None], 
                    positions= [(0,0), (70,1000), (1050,80), (570, 1000)],
                    crop=(0,32, 1568, 1470),  
                    outdir=None,
                    overwrite = False,
                    title = None,
                    ):
    if outdir is not None:
        outfile = os.path.join(outdir, sname )
    else:
        outfile =  sname
    if os.path.exists(outfile) and not overwrite: return False
    c = Holder(size=hsize)
    for i,(name, size,pos) in enumerate(zip(names, sizes, positions)):
        c.insert(name, size, pos)
    if title is not None: # expect a tuple (function, (l,t))
        title[0](sname, 'temp.png') # write a temporary file
        c.insert('temp.png', None, title[1])
        os.remove('temp.png')
    c.save(outfile, crop=crop)
    return True

def combine_image(name, path, subpaths=('tsmap', 'sedfig', 'log', 'light_curves'), outfolder='combined'):
    """ special to combine a create a single image
    """
    try:
        names = [glob.glob(os.path.join(path,subpath, '%s*.png'%name))[0] \
                    for subpath in subpaths ]
    except IndexError:
        raise InvalidParameter('Source name %s not found  in one of the folders %s'%(name, subpaths))
    combine_images(names, outdir=os.path.join(path, outfolder))

def make_dzc(infolder, outfolder,
        imagefolder='dzi', #'healpipe1_dzimages', 
        collection_name='dzc', #'healpipe1_dzcollection' 
        type = 'jpg',  
        ):
    """ infolder   :  path to a folder containing a bunch of images to convert
        outfolder  :   where to set up the dzc and dzi
    keyword parameters
        imagefolder: ['dzi'] name of a folder for the created dzi images
        collection_name: ['dzc'] name to apply to the deep zoom collection: will be the name of a xml file.
        type       : ['jpg'] graphics type, perhaps 'png'
    """
    imagefolder = imagefolder or '%s_dzi' % infolder
    collection_name = collection_name or '%s_dzc' % infolder
    t = '%s\\*.%s' % (infolder, type)
    n = len(glob.glob(t))
    if n ==0: raise InvalidParameter('no %s files found in folder "%s"' % (type,infolder))
    print 'Converting %d jpg images from %s, dzi to %s, collection %s' \
        % (n, infolder, imagefolder, collection_name)
    def cmd(c): 
        print c
        os.system(c)
    curdir = os.getcwd()
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    os.chdir(outfolder)
    print 'cd %s' %outfolder
    cmd('DZconvert %s %s' % (os.path.join(curdir,t), imagefolder) ) 
    cmd('DZcollection %s %s' % (imagefolder, collection_name))
    os.chdir(curdir)


def main(
        names,  figpath ,
        pivot_dir  = None, 
        nmax=0,
        force_log_gen = False,
        subpaths=('tsmap', 'sedfig', 'log', 'light_curves'),
        allow_missing=True,
    ):
    """
    Args:
    names:  list of source names to combine for Pivot
        assume that filename is the name with spaces converted to '_', plus sign to 'p'
    figpath: where to find the figures, assuming tsmap, sedfig are subfolders
    pivot_dir [None]: folder to save merged figures
    nmax [0]: if >0, limit number to convert (for testing)
    
    Will save combined figures in 'combined' folder
    """
    
    if force_log_gen or not os.path.exists(os.path.join(figpath, 'log', names[0])):
        convert_log_to_png(figpath)
        
    
    outdir = os.path.join(figpath, 'combined')
    if not os.path.exists(outdir): os.mkdir(outdir)
    print 'Creating combined images in folder %s' % outdir

    for name in names:
        try:
            combine_image(name.replace(' ','_').replace('+', 'p'), figpath, subpaths=subpaths)
        except:
            print 'fail to convert source "%s"' % name
            if not allow_missing: raise
            
    if pivot_dir is not None: make_dzc(outdir, pivot_dir)


if __name__=='__main__':
    from optparse import OptionParser
    usage = """\n
usage: %prog [options] figpath pivot_dir\n
Generate deep zoom collection for Pivot
    figpath: where to find the figures, assuming tsmap, sedfig, log are subfolders
    pivot_dir: folder to save merged figures
    
    Will save combined figures in 'combined' folder under figpath """
    parser = OptionParser(usage, version=version)
    parser.add_option('-n', '--nmax', help='if set, limit number to process', 
            action='store', dest='nmax',default=0, type='int')
    options, args = parser.parse_args()
    if len(args)!=2: 
        parser.print_usage()
        sys.exit(-1)
    main(args[0],args[1], nmax = options.nmax)
     
