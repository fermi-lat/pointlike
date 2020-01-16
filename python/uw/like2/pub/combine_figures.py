"""
Utility for making composite images for pivot
$Header$
"""
import os, glob
from uw.utilities import makefig
  

  
def get_files(fname, skymodel_dir, subdirs, allow_missing=1):
    """ helper function: return set of filenames associated with given source or ROI
    """

    def find_file(subfolder, filename):
        assert os.path.exists(os.path.join(skymodel_dir, subfolder)), 'folder %s not found' %subfolder
        filepat = os.path.join(skymodel_dir, subfolder, filename)
        #print (filepat,)
        t = glob.glob(filepat)
        if len(t)==0: 
            # not here: recursively check redirect folder
            rf =os.path.join(skymodel_dir, subfolder, 'redirect.txt')
            if os.path.exists(rf):
                redir = open(rf).read()
                return find_file(redir, filename)
            else: 
                f = None
        else: f = t[0]
        #print ('-->', f)
        return f
    
    ret = []
    missing = 0
    for sd in subdirs:
        filename = find_file(sd,'%s*.png' % fname)
        if filename is None :
            if missing<allow_missing: 
                missing +=1
            else:
                raise Exception, 'file for source %s not found in %s [%d missing]'\
                % (fname, skymodel_dir+'/'+sd, missing)
        ret.append(filename) # this might be the wrong one
    return ret

def make_composite_image(name, outdir, subdirs, target, overwrite, layout_kw):
       fname = name.replace(' ','_').replace('+','p')
       files = get_files(fname, outdir, subdirs)
       return makefig.combine_images(fname+'.jpg',files, outdir=target, overwrite=overwrite, **layout_kw)

def combine_images(names, outdir, subdirs, layout_kw,  outfolder, overwrite=False, mec=None):
    target = os.path.join(outdir,outfolder)
    if not os.path.exists(target): os.makedirs(target)
    missing = []
    for name in names:
        try: 
            ret= make_composite_image(name, outdir, subdirs, target, overwrite, layout_kw)
        except Exception, msg:
            raise
            missing.append(name)
    if len(missing)>0: 
        print ('Missing one or more figures in %d sources' %len(missing))

class CombineFigues(object):

    def __init__(self, skymodel_dir, outfolder='combined',
            overwrite=False, mec=None):
        self.skymodel =skymodel_dir
        assert os.path.exists(skymodel_dir), 'folder %s not found'
        self.overwrite=overwrite
        self.outfolder=outfolder
        self.mec = mec
    
    def __call__(self,sources=None, rois=None, log=True):
        if log and len(glob.glob(os.path.join(self.skymodel,'pnglog', '*.png')))==0:
            print ('converting log files...')
            makefig.convert_log_to_png(self.skymodel)
        if rois is not None:
            combine_images(rois, outdir=self.skymodel, 
                subdirs=('kde_maps','ts_maps', 'counts_dir', 'pnglog'), 
                layout_kw=dict(
                    hsize=(1600,1100), 
                    sizes=     [None, None, None, None], 
                    positions= [(0,0),(0,500), (600, 50), (900,0)],
                    crop=None,  
                    #title = (make_title, (800,0)), # pass function to make title, position
                ) ,
                outfolder=self.outfolder, overwrite = self.overwrite, mec=self.mec)
        if sources is not None:  
            # updated with stuff from Eric's light curves
            combine_images(sources, outdir=self.skymodel, 
                subdirs=('tsmap', 'sedfig', 'light_curves'),
                layout_kw=dict(
                   hsize=(900,700), 
                   sizes=     [None, None, None], #(500,500) , (500,500) , (800,600)], 
                   positions= [(100,0), (500, 0), (0,400)],
                   crop= None,
                   ),
                outfolder=self.outfolder, overwrite=self.overwrite,mec=self.mec)
            

