"""
Utility for making composite images for pivot
"""
import os, glob
from uw.utilities import makefig
    
def get_files(fname, outdir, subdirs):
    ret = []
    for sd in subdirs:
        t = glob.glob(os.path.join(outdir, sd, '%s*.png' % fname))
        if len(t)==0: raise Exception, 'file for source %s not found in %s' % (fname, sd)
        elif len(t)>1:
            ret.append(t[0]) # this might be the wrong one
        else:    
            ret.append(t[0])
    return ret

def make_composite_image(name, outdir, subdirs, target, overwrite, layout_kw):
       fname = name.replace(' ','_').replace('+','p')
       files = get_files(fname, outdir, subdirs)
       return makefig.combine_images(fname+'.jpg',files, outdir=target, overwrite=overwrite, **layout_kw)

def combine_images(names, outdir, subdirs, layout_kw,  outfolder, overwrite=False, mec=None):
    target = os.path.join(outdir,outfolder)
    if not os.path.exists(target): os.makedirs(target)
    for name in names:
        ret= make_composite_image(name, outdir, subdirs, target, overwrite, layout_kw)

class CombineFigues(object):

    def __init__(self, outdir=None, outfolder='combined', overwrite=False, mec=None):
        if outdir is None:
            outdir = 'uw%2d' % int(open('version.txt').read())
            assert os.path.exists(outdir), 'folder %s not found'
        self.outdir =outdir
        self.overwrite=overwrite
        self.outfolder=outfolder
        self.mec = mec
    
    def __call__(self,sources=None, rois=None):
        if len(glob.glob(os.path.join(self.outdir,'log', '*.png')))==0:
            print 'converting log files...'
            makefig.convert_log_to_png(self.outdir)
        if rois:
            combine_images(rois, outdir=self.outdir, 
                subdirs=('kde_figs','residual_ts', 'counts_dir', 'log'), 
                layout_kw=dict(
                    hsize=(1800,1100), 
                    sizes=     [(600,600), (600,600), None, None], 
                    positions= [(0,0),(0,500), (575, 50), (1175,0)],
                    crop=None, #(0,32, 1568, 1470), 
                    #title = (make_title, (800,0)), # pass function to make title, position
                ) ,
                outfolder=self.outfolder, overwrite = self.overwrite, mec=self.mec)
        if sources:        
            combine_images(sources, outdir=self.outdir, 
                subdirs=('tsmap', 'sedfig', ),
                layout_kw=dict(
                    hsize=(1100,600), 
                    sizes=     [(600,600) , (500,500) , ], 
                    positions= [(0,0), (600, 0)],
                    crop= None,#(0,32, 1000, 1470), 
                ) ,
                outfolder=self.outfolder, overwrite=self.overwrite,mec=self.mec)

