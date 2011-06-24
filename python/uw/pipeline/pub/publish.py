"""
manage publishing 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/pub/publish.py,v 1.6 2011/04/09 14:18:37 burnett Exp $
"""
import sys, os, pickle, glob, types, time
import PIL
import pyfits
import numpy as np
import pylab as plt
from uw.utilities import makefig, makerec, colormaps
from . import healpix_map, dz_collection 
from . import source_pivot, roi_pivot, makecat, display_map, healpix_map
from ..analysis import residuals
from .. import skymodel

from skymaps import Band, SkyDir 

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

def make_title(title, filename=None):
    plt.close(100);
    plt.figure(100,figsize=(5,0.5));
    plt.figtext(0.1,0.8, title, va='top',ha='left', fontsize=24)    
    if filename is not None: plt.savefig(filename)
    
def combine_images(names, outdir, subdirs, layout_kw,  outfolder, overwrite=False,):
    if not np.all(map(lambda subdir: os.path.exists(os.path.join(outdir,subdir)), subdirs)): 
        print 'WARNING not all subdirs, %s, exist: not combining images' % list(subdirs)
        return
    target = os.path.join(outdir,outfolder)
    if not os.path.exists(target): os.makedirs(target)
    for name in names:
        fname = name.replace(' ','_').replace('+','p')
        files = get_files(fname, outdir, subdirs)
        if makefig.combine_images(fname+'.jpg',  files,  outdir=target, overwrite=overwrite, **layout_kw):
            print '.',

def make_thumbnail(outfile):
        im = PIL.Image.open(outfile)
        im.thumbnail((200,100))
        im.save(outfile.replace('.png', '_thumbnail.png'))

class Publish(object):
    """ manage generating all the products from a build
    """
    def __init__(self, outdir=None, 
            pivot_root = r'D:\common\pivot' if os.name=='nt' else '/phys/users/tburnett/pivot',
            ts_min = 10,
            catalog_name = None,
            dzc = 'dzc.xml',
            overwrite=False,
            mec=None,
            **kwargs
            ):
            
        """
            outdir : 
            ts_min : float,or None
                apply cut to sources
            catalog_name : string, or None
                use as name, or create with root '24M_' + outdir
            mec : None or IPython.kernel.engine.client.MultiEngineClient object, to use for multprocessing
        """
        self.source_pivot_kw = kwargs.pop('source_pivot_kw', {})
        if outdir is None:
            outdir='uw%02d' % int(open('version.txt').read())
        if type(outdir)==types.StringType:
            self.outdir=outdir
            self.refdir=''
        else:
            self.outdir,self.refdir=outdir
        self.mec=mec
        self.get_config()
        self.skymodel = skymodel.SkyModel(outdir, extended_catalog_name=self.config['extended'])
        self.rois = self.skymodel.roi_rec()
        self.sources = self.skymodel.source_rec()
        if ts_min is not None:
            self.sources = self.sources[self.sources.ts>ts_min]
            print 'Selected subset: %d sources with ts>%.0f' % (len(self.sources), ts_min)
        self.ts_min = ts_min
        self.name=catalog_name if catalog_name is not None else os.path.split(os.getcwd())[-1]+'_'+self.outdir
        print 'using name %s for files' % self.name

        self.pivot_dir = os.path.join(pivot_root, self.name)
        if self.refdir!='':
            self.ref_pivot_dir = os.path.join(pivot_root, self.name.replace(self.outdir, self.refdir))
        self.dzc = dzc
        if not os.path.exists(self.pivot_dir): os.mkdir(self.pivot_dir)
        print 'writing files to %s' % self.pivot_dir
        files = glob.glob(os.path.join(self.pivot_dir, '*'))
        print 'existing files:\n\t%s' % '\n\t'.join(files)
        self.overwrite=overwrite
        self.zips =     ( 'sedfig',      'tsmap',        'ts_maps',      'kde_maps')
        self.zipnames = ( 'source seds', 'source tsmaps','roi residual ts','roi kdemaps')
        for z in self.zips:
            if not os.path.exists(os.path.join(self.outdir,self.ref(z))):
                print 'WARNING: missing folder in outdir: %s' % z
    
    def get_config(self, fn = 'config.txt'):
        """ parse the items in the configuration file into a dictionary
        """
        file = open(os.path.join(self.outdir, fn))
        self.config={}
        for line in file:
            item = line.split(':')
            if len(item)>1:
                self.config[item[0].strip()]=item[1].strip()
        
        if 'extended' not in self.config: self.config['extended']=None
        
    def ref(self,x):
        if self.refdir =='': return x
        return '../%s/%s' %( self.refdir,x)
    
    def roi_plots(self):
        r = self.rois
        for field, vmin, vmax in (
             ('galnorm', 0.5,1.4),
             ('isonorm', 0,  2.0),
             ('limbnorm', 0,  2.0),
             ('chisq',   0,  50 )):
            dm = display_map.DisplayMap(r.field(field))
            dm.fill_ait(show_kw=dict(vmin=vmin, vmax=vmax))
            plt.title( '%s for %s' % (field, self.outdir))
            outfile = os.path.join(self.pivot_dir,'roi_%s.png'% field)
            plt.savefig(outfile, bbox_inches='tight')
            make_thumbnail(outfile)

    def ait_plots(self,   **kwargs):
        """ images from the high-resolution HEALPix, extract from fields of aladin512.fits
            kwargs: can specify pixelsize, dpi
        """
        show_kw_dict=dict( 
            kde=dict(nocolorbar=True, scale='log',vmin=4.5,vmax=7.5, cmap=colormaps.sls),
            ts = dict(nocolorbar=False, vmin=10, vmax=25),
            galactic = dict(nocolorbar=True, scale='log', cmap=colormaps.sls),
            counts = dict(scale='log'),
            )
        dpi = kwargs.pop('dpi', 120)
        t = pyfits.open(os.path.join(self.pivot_dir,'aladin512.fits'))[1].data
        for table in t.dtype.names:
            outfile = self._check_exist(table+'_ait.png')
            if outfile is None: continue
            
            dm = display_map.DisplayMap(t.field(table))
            show_kw=show_kw_dict.get(table, {})
            dm.fill_ait(show_kw=show_kw, **kwargs)
            #plt.title( '%s for %s' % (field, self.outdir))
            plt.savefig(outfile, bbox_inches='tight', dpi=dpi)
            make_thumbnail(outfile)
            print 'wrote %s and its thumbndail' % outfile

    def roi_fit_html(self):
        def entry(name):
            return """<b>%(title)s</b> <br/> <a href="roi_%(name)s.png"> <img alt="%(name)s_ait.png"  src="roi_%(name)s_thumbnail.png" /></a> """\
                %dict(name=name, 
                    title = dict(isonorm='Isotropic nomalization', 
                        galnorm='Galactic normalization', 
                        limbnorm='Limb normalization',
                        chisq='Chi Squared', kde='photon density')[name])
        files = glob.glob(os.path.join(self.pivot_dir,'roi_*.png'))
        if len(files)==0: return ''
        return '<table cellspacing="0" cellpadding="5"><tr>\n' +'\n'.join(['\t<td>%s</td>'%entry(f) for f in 'galnorm isonorm limbnorm chisq'.split()])+'</tr></table>'

    def source_plot(self, field='kde', healpix_file='aladin512.fits'):
        s = self.sources
        try:
            kde =pyfits.open(os.path.join(self.outdir,healpix_file))[1].data.field(field)
        except:
            print 'cannot create a source plot over the photon density: '\
                    'field %s not found in %s'% (field, healpix_file)
            return
        sm = display_map.SourceMap(kde,s)
        sm.fill_ait()
        if len(s)>0: sm.add_sources()
        outfile = os.path.join(self.pivot_dir,'sources_ait.png')
        plt.savefig(outfile, bbox_inches='tight')
        make_thumbnail(outfile)
    
        
    def combine_plots(self, sources=True, rois=True, dzc=True):
        if   len(glob.glob(os.path.join(self.outdir,'log', '*.png'))) \
            <len(glob.glob(os.path.join(self.outdir,'log', '*.txt'))):
            print 'converting log files...'
            makefig.convert_log_to_png(self.outdir)
        if rois:
            combine_images(self.rois.name, outdir=self.outdir, 
                #subdirs=('countsmap','residual_ts', 'counts_dir', 'log'), #gtsmap_dir='gtsmap', 
                subdirs=(self.ref('kde_maps'),'ts_maps', 'counts_dir', 'log'), #gtsmap_dir='gtsmap', 
                layout_kw=dict(
                    hsize=(1800,1100), 
                    sizes=     [(600,600), (600,600), None, None], 
                    positions= [(0,0),(0,500), (600, 50), (1250,0)],
                    crop=None, #(0,32, 1568, 1470), 
                    #title = (make_title, (800,0)), # pass function to make title, position
                ) ,
                outfolder='combined', overwrite = self.overwrite, )
        if sources:        
            combine_images(self.sources.name, outdir=self.outdir, 
                subdirs=('tsmap', 'sedfig', ),
                layout_kw=dict(
                    hsize=(1100,600), 
                    sizes=     [(600,600) , (500,500) , ], 
                    positions= [(0,0), (600, 0)],
                    crop= None,#(0,32, 1000, 1470), 
                ) ,
                outfolder='combined', overwrite=self.overwrite)
        if dzc: 
            self.make_dzc()

    def make_dzc(self):
        if not os.path.exists(os.path.join(self.pivot_dir,'dzc.xml')) or self.overwrite:
            mc= dz_collection.MakeCollection(os.path.join(self.outdir, 'combined'), self.pivot_dir)
            mc.convert(self.mec) # use the MEC if provided!
            mc.collect()
        
    def make_pivot(self):
        roi_pivot.make_pivot(self.rois, outdir=self.outdir,
            pivot_dir=self.pivot_dir,  
            pivot_name='%s - rois'%self.outdir , 
            pivot_file='rois.cxml',
            dzc = self.dzc,
            )
        self.write_source_pivot()
        
    def write_source_pivot(self):
        source_pivot.make_pivot(self.sources, outdir=self.outdir,
            pivot_dir=self.pivot_dir, 
            pivot_name='%s - sources'%self.outdir,
            pivot_file='sources.cxml',  
            dzc = self.dzc,
            ** self.source_pivot_kw
            )
            
    def write_xml(self, catname=None):
        """ generate the xml a of the catalog in the folder with the pivots"""
        fn = os.path.join(self.pivot_dir, self.name+'.xml')
        self.skymodel.toXML(fn, ts_min=self.ts_min, title='catalog %s sources'%self.name)
    
    def _check_exist(self, filename):
        """ return full filename, or None if the file exists and overwrite is not set"""
        fn = os.path.join(self.pivot_dir, filename)
        if os.path.exists(fn):
            if self.overwrite: os.remove(fn)
            else: return None
        return fn
    
    def write_reg(self):
        """ generate the reg file """
        fn = self._check_exist(self.name+'.reg.txt')
        if fn is None: return
        self.skymodel.write_reg_file(fn, ts_min=self.ts_min)
        print 'wrote reg file %s' %fn
        
    def write_FITS(self, TSmin=None):
        """ write out the Catalog format FITS version, and the rois """
        catname = os.path.join(self.pivot_dir,self.name+'.fits') 
        if os.path.exists(catname):
            if self.overwrite: os.remove(catname)
            else: 
                return
        cat = makecat.MakeCat(self.sources, TScut=0 if TSmin is None else TSmin)
        cat(catname)
        for rec, rectype in ((self.sources, 'sources'), (self.rois,'rois')):
            outfile = os.path.join(self.pivot_dir,self.name+'_%s.fits'%rectype)
            if os.path.exists(outfile):
                if self.overwrite: os.remove(outfile)
                else:
                    print 'file %s exists: set overwrite to replace it' % catname
                    continue
            makerec.makefits(rec, outfile)
            print 'wrote %s' %outfile
        #self.write_xml_and_reg(catname)
    
    def write_map(self, name, fun=lambda x:x, title=None, vmax=None, table=None):
        """ DISABLED for now 
        write out image files for a single table
        name : string 
            name of the table, to be loaded by roi_maps.load_tables if table is None 
        fun  : ufunc function
            ufunc to apply to each pixel when generating the image files
        title : string
            title to apply to plots
        table : None, or a HEALpix array
            """
        if table is None:
            table = healpix_map.load_tables(name, outdir=self.outdir)
        nside = int(np.sqrt(len(table)/12))
        assert len(table)==12*nside**2, 'table length, %d, not HEALpix compatible' % len(table)

        filename = '%s_ait.fits' % name
        fitsfile = os.path.join(self.pivot_dir,filename)
        if os.path.exists(fitsfile):
            if self.overwrite:  os.remove(fitsfile)
            else: fitsfile = ''

        outfile=os.path.join(self.pivot_dir,'%s_ait.png'%name)
        if fitsfile != '':
            # first without the usual scaling function for the FITS version
            print 'generating FITS image %s' %fitsfile
            q=display_map.skyplot(table, title,
                ait_kw=dict(fitsfile=fitsfile,pixelsize=0.1))
            del q
        if not os.path.exists(outfile) or self.overwrite:
            print 'generating %s and its thumbnail...' % outfile
            # now with scaling function to generate png
            plt.close(30)
            fig = plt.figure(30, figsize=(16,8))
            q=display_map.skyplot(fun(table), title, axes=fig.gca(),
                ait_kw=dict(fitsfile='', pixelsize=0.1), vmax=vmax)
            plt.savefig(outfile, bbox_inches='tight', pad_inches=0)
            im = PIL.Image.open(outfile)
            im.thumbnail((200,100))
            im.save(outfile.replace('.png', '_thumbnail.png'))
            del q
     
        outfile = os.path.join(self.pivot_dir,'%s_map.fits' % name)
        if os.path.exists(outfile):
            if self.overwrite: os.remove(outfile)
            else: return
        print 'generating %s' % outfile
        dirs = map(Band(nside).dir, xrange(len(table)))
        ra  = np.array(map(lambda s: s.ra(), dirs), np.float32)
        dec = np.array(map(lambda s: s.dec(),dirs), np.float32)
        outrec = np.rec.fromarrays([ra,dec,np.array(table,np.float32)], 
            names = ['ra','dec', 'value'])
        makerec.makefits(outrec, outfile)

    
    def write_maps(self):
        print 'write_maps is disabled for now'
        return
        for pars in ( 
            #('data', np.log10, 'log10(counts)', 4),
            ('ts',     np.sqrt,  'sqrt(ts)', 10),
            ('kde',    np.sqrt,  'sqrt(kde)', 5000),
            ('ring', lambda x: np.log10(x), eval(self.config['diffuse'])[0], None, 
                     pickle.load(open(os.path.join(self.outdir, 'galactic.pickle')))),
            ):
            if pars[0]=='kde' and self.refdir!='': continue
            try:
                self.write_map(*pars)
            except:
                print 'failed to process %s' % pars
    
    def write_hpmaps(self, outfile='aladin512.fits', nside=512):
        outfile = os.path.join(self.pivot_dir, outfile)
        if os.path.exists(outfile):
            if self.overwrite: 
                os.remove(outfile)
            else: return
        # look for tables to include
        t = glob.glob(os.path.join(self.outdir, '*_table'))
        cols = healpix_map.tables(self.outdir,names=map(lambda x: os.path.split(x)[1][:-6], t))
        # add the diffuse directly
        gal = eval(self.config['diffuse'])[0]
        cols.append(healpix_map.diffusefun(gal, nside=nside))
        healpix_map.HEALPixFITS(cols, nside=nside).write(outfile)
        print 'wrote file %s with %d columns' %(outfile, len(cols))
        
    def write_residualmaps(self, outfile='aladin64.fits',nside=64):
        """ generate set of maps with residuals, chisq, etc"""
        outfile = os.path.join(self.pivot_dir, outfile)
        if os.path.exists(outfile):
            if self.overwrite: 
                os.remove(outfile)
            else: return
        residuals.tofits(self.outdir, outfile, )
        
    
    
    def write_zips(self):
        """ generate zip files containing the contents of a folder
        """
        for dirname in self.zips:
            fulldir = os.path.join(self.outdir, dirname)
            if not os.path.exists(fulldir):
                print 'did not find folder %s: continue' % fulldir
                continue
            tozip = os.path.join(self.pivot_dir, '%s_images.zip' % dirname)
            if os.path.exists(tozip) and not self.overwrite: continue
            cmd = 'zip -r %s %s' %(tozip, fulldir)
            print cmd
            os.system(cmd)
    
    def image_html(self):
        """ generate a list of images: look for 'name_ait.png' files generated by write_map where name is in the dict below.
        """
        template=""" <td>%(title)s <br/> 
            <!--<a href="%(path)s_ait.fits">[FITS AIT Image]</a> 
            <a href="%(path)s_map.fits">[FITS table]</a><br/> -->
            <a href="%(path)s_ait.png"> 
            <img alt="%(path)s_ait.png"  
                 src="%(path)s_ait_thumbnail.png" /></a> <br/>
        </td>"""
        q = glob.glob(os.path.join(self.pivot_dir,'*_ait.png'))
        if self.refdir !='':
            q.append(os.path.join(self.ref_pivot_dir,'kde_ait.png'))
        ret = '<ul>\n'
        if os.path.exists(os.path.join(self.pivot_dir, 'aladin512.fits')):
            ret += '<li><a href="aladin512.fits">HEALPix-format file for aladin</a> contains columns %s</li>'\
                % pyfits.open(os.path.join(self.pivot_dir, 'aladin512.fits'))[1].data.names
        
        if os.path.exists(os.path.join(self.pivot_dir, 'aladin64.fits')):
            ret += '<li><a href="aladin64.fits">HEALPix-format file for aladin</a> contains columns %s</li>'\
                % pyfits.open(os.path.join(self.pivot_dir, 'aladin64.fits'))[1].data.names
                
        if len(q)==0: 
            print  'WARNING: no image files found'
            return ret +'\n</ul>\n'
        heads = [os.path.split(t)[0] for t in q]
        names = [ os.path.split(t)[-1].split('_')[0] for t in q]
        paths = ['../%s/'% os.path.split(t)[-1] for t in heads] # relative path to the file
        titles = [dict( kde='<b>Photon density</b> (photons/sr)<br/>using a <a href="../images/kde_reference.PNG">Kernel Density Estimator</a>', 
                        ts=  '<a href="../notes/residual_ts.htm"><b>Residual TS</b></a><br/> assuming powerlaw index=2.0',
                        ts15='<a href="../notes/residual_ts.htm">Residual TS</a> assuming powerlaw index=1.5',
                        ts25='<a href="../notes/residual_ts.htm">Residual TS</a> assuming powerlaw index=2.5',
                        counts='<b>Counts map</b>, E>1GeV',
                        galactic='<b>Galactic diffuse</b>, at 1 GeV',
                        sources='<b>Source locations</b>',
                    )[name] for name in names]
        return ret+'\n<table cellpadding="5" cellspacing="0"><tr>'\
            + '\n  '.join([template%dict(title=a, path=c+b) for a,b,c in zip(titles,names,paths)])\
            +'\n</tr></table>\n</ul>\n' 

    def write_html(self, tables=False):
        """ write out a little overview to go in the folder"""
        
        d = dict(catname=self.name,
            zipfiles='\n'.join(['<li><a href="%s_images.zip">%s <a/></li>'
                % (file,name) for file,name in zip(self.zips, self.zipnames)\
                    if os.path.exists(os.path.join(self.pivot_dir,file+'_images.zip'))]),
            diffuse = eval(self.config['diffuse']),
            irf = self.config['irf'],
            extended= self.config['extended'],
            roi_fit = self.roi_fit_html(),
            ts10=np.sum(self.sources.ts>10), ts25=np.sum(self.sources.ts>25),
            imagefiles='',
        )
        if not tables: d.update( imagefiles=self.image_html())

        html_basic="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Language" content="en-us" />
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Fermi-LAT sky model %(catname)s</title></head>
<body><h2>Fermi-LAT sky model version %(catname)s</h2>
    <dl>
    <dt>IRF: <dd>%(irf)s
    <dt>Diffuse components: <dd>%(diffuse)s
    <dt>Extended definition: <dd>%(extended)s
    </dl>

<h3>Generated files:</h3>
<ul>
 <li><a href="%(catname)s.fits">Catalog-format source fits file</a> (Standard list for external. 
    <br/>There is an internal version with more information<a href="%(catname)s_sources.fits">(FITS)</a>)
    <br/>Total sources with TS>10: %(ts10)d; TS>25: %(ts25)d
 </li>
 <li>ROI information<a href="%(catname)s_rois.fits"> (FITS)</a> <br/>
    ROI fit parameters, chisquared: %(roi_fit)s
 </li>
 <li><a href="%(catname)s.xml">XML-format file for gtlike</a></li>
 <li><a href="%(catname)s.reg.txt">Region file for DS9 </a>(note, should be renamed back to '.reg')</li>
 <li>ZIP files, if any, of ROI- or source-based images. Note: TS cut is not applied, all sources in the model are included.
 </h3> <ul> %(zipfiles)s</li></ul>
 </ul>
"""%d
        html_images="""<h3>All-sky HEALPix Images:</h3>
  %(imagefiles)s
<h3><a href="http://www.silverlight.net/learn/pivotviewer/">Pivot</a> Collections</h3>
<ul><li><a href="http://fermi-ts.phys.washington.edu/pivot/viewer/?uw=%(catname)s/sources.cxml">sources</a>
    All sources in the model.
   </li>
	<li><a href="http://fermi-ts.phys.washington.edu/pivot/viewer/?uw=%(catname)s/rois.cxml">rois</a>
    Entries for all 1728 ROIs
    </li>
</ul>"""%d
        html_tail="<hr>Last update: %s</body></html>"% time.asctime()    
        open(os.path.join(self.pivot_dir, 'default.htm'),'w').write(html_basic + ( html_images if not tables else '') + html_tail)
        print 'wrote HTML file %s' % os.path.join(self.pivot_dir, 'default.htm')

    def doall(self):
        self.combine_plots()
        self.source_plot()
        self.roi_plots()
        self.make_pivot()
        self.write_FITS()
        self.write_maps() 
        self.write_hpmaps()
        self.ait_plots()
        self.write_residualmaps()
        self.write_zips()
        self.write_xml()
        self.write_reg()
        self.write_html()
    def tables(self):
        self.write_FITS()
        self.write_zips()
        self.write_xml()
        self.write_reg()
        self.roi_plots()
        self.write_html(tables=True)
        
def main(outdir=None, **kwargs):
    return Publish(outdir, **kwargs) 
 
def doall(outdir, **kwargs):
    pub = Publish(outdir, **kwargs) 
    pub.combine_plots()
    pub.source_plot()
    pub.roi_plots()
    pub.make_pivot()
    pub.write_FITS()
    pub.write_maps() 
    pub.write_zips()
    pub.write_xml()
    pub.write_reg()
    pub.write_html()
    
    
def robocopy(fromdir, todir=None, path=r'P:\ '):
    if todir is None: todir = fromdir
    cmd = r'C:\Windows\System32\Robocopy.exe %s d:\common\pivot\%s /MIR' %(fromdir, todir)
    os.system('start /path=%s %s'% (path, cmd))
    
if __name__=='__main__':
    try:
        doall(sys.argv[1])
    except:
        print 'expected one arg, the name of the release, like "uw64"'
        raise