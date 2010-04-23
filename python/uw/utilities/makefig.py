""" 
Make combinded figures for Pivot, perhaps
""" 

from PIL import Image
import glob, os, sys, exceptions
import pylab  as plt

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

w= 1536 # default initial width -- numbers below tuned by hand for tsmap, sed, log
def combine_images(names, hsize=(2304,w), 
                    sizes=     [(w,w) , None, None], 
                    positions= [(0,0), (w,75), (w,w-1015)],
                    crop=(50,50,2050, 1500),
                    outdir=None,
                    ):
    c = Holder(size=hsize)
    for i,(name, size,pos) in enumerate(zip(names, sizes, positions)):
        c.insert(name, size, pos)

    #sname = os.path.split(names[0])[1].split('_')[-2]
    sname = os.path.split(names[0])[1].replace('_tsmap.png','.jpg')
    print sname
    if outdir is not None:
        outfile = os.path.join(outdir, sname )
    else:
        outfile =  sname
    c.save(outfile, crop=crop)


def make_dzc(infolder, outfolder,
        imagefolder='dzi', #'healpipe1_dzimages', 
        collection_name='dzc', #'healpipe1_dzcollection' 
        ):
    """ infolder: path to a folder containing a bunch of jpg images to convert
        outfolder: where to set up the dzc and dzi
        imagefolder: name of a folder for the created dzi images
        collection_name: name to apply to the deep zoom collection: will be the name of a xml file.
    """
    imagefolder = imagefolder or '%s_dzi' % infolder
    collection_name = collection_name or '%s_dzc' % infolder
    t = '%s\\*.jpg' % infolder
    n = len(glob.glob(t))
    if n ==0: raise InvalidParameter('no jpg files found in folder "%s"' % infolder)
    print 'Converting %d jpg images from %s, dzi to %s, collection %s' \
        % (n, infolder, imagefolder, collection_name)
    def cmd(c): 
        print c
        os.system(c)
    curdir = os.getcwd()
    if not os.path.exists(outfolder): os.mkdir(outfolder)
    os.chdir(outfolder)
    cmd('DZconvert %s %s' % (os.path.join(curdir,t), imagefolder) ) 
    cmd('DZcollection %s %s' % (imagefolder, collection_name))
    os.chdir(curdir)


def main(
        figpath ,
        pivot_dir  = None, 
        nmax=0,
    ):
    """
    Args:
    figpath: where to find the figures, assuming tsmap, sedfig are subfolders
    pivot_dir [None]: folder to save merged figures
    nmax [0]: if >0, limit number to convert (for testing)
    
    Will save combined figures in 'combined' folder
    """
    
    tsmaps = glob.glob(os.path.join(figpath,'tsmap','*.png'))
    tsmaps.sort()
    sedfigs= glob.glob(os.path.join(figpath,'sedfig','*.png'))
    sedfigs.sort()
    n = len(tsmaps)
    logfigs = glob.glob(os.path.join(figpath,'log', '*.png'))
    logfigs.sort()
    if n>0 and len(logfigs)!=n:
        print 'generating png versions of log files...'
        convert_log_to_png(figpath)
        logfigs = glob.glob(os.path.join(figpath,'log', '*.png'))
        logfigs.sort()

        
    print 'found %d tsmap, %d sedfig images, and %d logfigs in %s' % (n, len(sedfigs), len(logfigs), figpath)
    if (n==0 or n!=len(sedfigs) or n!=len(logfigs)):
        raise InvalidParameter('expected equal numbers of files in each folder to combine')
    outdir = os.path.join(figpath, 'combined')
    if not os.path.exists(outdir): os.mkdir(outdir)
    print 'Creating combined images in folder %s' % outdir

    for i,names in enumerate(zip(tsmaps, sedfigs,logfigs)):
        if nmax>0 and i>nmax: break
        
        combine_images(names, outdir=outdir)
        if i%500==0: print '.',
        
    if pivot_dir is not None: make_dzc(outdir, pivot_dir)


if __name__=='__main__':
    argv = sys.argv
    if len(argv)<3: raise InvalidParameter('expect 3 parameters: figpath, outdir')
    figpath = argv[1]
    pivotdir= argv[2]
    collection_name = argv[3]
    main(figpath, pivotdir)
