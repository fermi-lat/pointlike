"""
Code to generate a standard Fermi-LAT catalog FITS file
also, see to_xml, to generate XML for the sources
$Header$
"""
import os
import pyfits
from skymaps import SkyDir

#catdir  = config.catalog_path
#catalog = config.default_catalog
#assoc   = 'gll_psc18month_uw8_assoc.fits' # file from Jean with associations
##newcat  = '24M_%s.fits' %outdir

#def get_rec(outdir,name='sources'):    
#    return pickle.load(open('%s_%s.rec'%(name, outdir))) 
#def full_path(fn):
#    return os.path.join(catdir, fn)
#
#
#def to_xml(outdir):
#    recfile ='sources_%s.rec'%outdir 
#    assert os.path.exists(recfile), 'pickled rec file %s not found' %recfile
#    cat = CatalogManager(recfile)
#    ps = map(cat.point_source, cat.dirs, cat.names, cat.models)
#    stack=xml_parsers.unparse_point_sources(ps)
#    xmlfile = '24M_%s.xml'%outdir
#    xml_parsers.writeXML(stack, xmlfile, '24M_%s source library'%outdir)
#    print 'wrote out %d source entries XML file %s ' % (len(stack), xmlfile)
#
#class Assoc(object):
#    def __init__(self,cat, curcat, cat_ref=assoc):
#        self.hdu = pyfits.open(full_path(curcat))
#        self.data = self.hdu[1].data
#        self.cat_ref = pyfits.open(full_path(cat_ref))[2]
#    def __call__(self, name):
#        return self.data[self.data.field('NickName')==name]
#
# default column definitions 
coldata ="""\
    NickName               20A None
    RA                       E deg
    DEC                      E deg
    GLON                     E deg
    GLAT                     E deg
    Conf_95_SemiMajor        E deg
    Conf_95_SemiMinor        E deg
    Conf_95_PosAng           E deg
    Test_Statistic           E None
    Pivot_Energy             E MeV
    Cutoff_Energy            E MeV
    Flux_Density             E photon/cm**2/MeV/s
    Unc_Flux_Density         E photon/cm**2/MeV/s
    Spectral_Index           E None
    Unc_Spectral_Index       E None
    Flux1000                 E photon/cm**2/s
    Unc_Flux1000             E photon/cm**2/s
    Energy_Flux              E erg/cm**2/s
    Unc_Energy_Flux          E erg/cm**2/s
    SpectralFitQuality       E None
    ID_Number                I None
    ID_Name               520A None
    ID_Probability         26E None
    ID_RA                  26E deg
    ID_DEC                 26E deg
    ID_Angsep              26E deg
    ID_Catalog             26I None""".split('\n')
coldict = {}
for name,format,unit in [c.split() for c in coldata]:
    coldict[name]= dict(format=format, unit = unit)

def makecol(name, array):
    if name in coldict:
        format = coldict['name']['format']
        unit = coldict['name']['unit']
        if unit=='None': unit=''
    elif name[:4] == 'Flux' or name[:8]=='Unc_Flux':
        format, unit = 'E', 'photon/cm**2/s'

    else:
        format, unit = 'E', ''
    return pyfits.Column(name=name, format=format, unit=unit, array=array)
        
def band_cols():
    elist = (100,300,1000,3000,10000,100000)
    for i in range(len(elist)-1):
        e1,e2 = elist[i:i+2]
        print '%d_%d' % (e1,e2)
  
def get_data():
    return pipeline.load_rec_from_pickles(outdir) 
   
def source_class(z):
    source = [{'1F':'1FGL', 'PG':'PGW', 'MR':'MRF', 'UW':'UW', 'MS':'MST', 'SE':'SEED', '18':'18M',
        'Cy':'bin', 'LS':'bin','PS':'PSR',}[n[:2]] for n in z.name]
    return source
    
def move(z):
    source = [{'1F':'1FGL', 'PG':'PGW', 'MR':'MRF', 'UW':'UW', 'MS':'MST', 'SE':'SEED', '18':'18M',
        'Cy':'bin', 'LS':'bin','PS':'PSR',}[n[:2]] for n in z.name]
    fgl = np.array(source)=='1FGL'
    em  = np.array(source)=='18M'
    canmove = em+fgl
    print 'moving %d sources' % sum(canmove)
    z.ra[canmove]=z.fit_ra[canmove]
    z.dec[canmove]=z.fit_dec[canmove]

    
class MakeCat(object):
    
    def __init__(self, z,  canmove=None, TScut=0, add_assoc=False):
        self.z = z  
        self.canmove = canmove
        self.TScut = TScut
        self.add_assoc = add_assoc
        
    def add(self, name, array, fill=0):
        print ' %s ' % name ,
        if name in coldict:
            format = coldict[name]['format']
            unit = coldict[name]['unit']
        else:
            format, unit = 'E', ''
        t = array
        if self.check:
            t = array.copy()
            t[self.bad] = fill
        self.cols.append(pyfits.Column(name=name, format=format, unit=unit, array=t))
        
        
    def __call__(self, outfile):
        self.cols = []
        z = self.z[self.z.ts>self.TScut] # limit for now
        z.sort(order=('ra'))
        #z.ts = z.ts2 #kluge for now
        self.check=False
        self.bad = z.ts<9
        self.add('NickName', z.name)
        self.add('RA', z.ra)
        self.add('DEC', z.dec)
        sdir = map(SkyDir, z.ra, z.dec)
        self.add('GLON', [s.l() for s in sdir])
        self.add('GLAT', [s.b() for s in sdir])
        f95 = 2.45*1.1 # from 
        self.add('Conf_95_SemiMajor', f95*z.a)
        self.add('Conf_95_SemiMinor', f95*z.b)
        try:
            self.add('Conf_95_PosAng',    z.ang)
        except:
            self.add('Conf_95_PosAng',    z.phi)
        self.add('Test_Statistic',    z.ts)
        #self.check = True # rows for bad sources set to zero
        self.add('Pivot_Energy',      z.e0)  # note that pivot_energy is the measured value
        self.add('Flux_Density',      z.pnorm)
        self.add('Unc_Flux_Density',  z.pnorm_unc)
        self.add('Spectral_Index',    z.pindex)
        self.add('Unc_Spectral_Index',z.pindex_unc)
        try:
            self.add('beta',              z.beta)
            self.add('Unc_beta',          z.beta_unc)
        except:
            print 'attempt to add beta, Unc_beta failed'
        self.add('Cutoff_Energy',     z.cutoff) ## need to get this from info
        self.add('Cutoff_Energy_Unc', z.cutoff_unc) ## need to get this from info
        self.add('SpectralFitQuality',    z.band_ts-z.ts)  # sort of an approximation?
        #if self.add_assoc:
        #    assoc = Assoc()
        #    for idcol in 'Number Name Probability RA DEC Angsep Catalog'.split():
        #        h = 'ID_'+idcol
        #        adata = np.array([assoc(name).field(h)[0] for name in z.name])
        #        self.add(h, adata)
        
        # make the FITS stuff
        table = pyfits.new_table(self.cols)
        table.name = 'LAT_Point_Source_Catalog' 
        if os.path.exists(outfile):
            os.remove(outfile)
        self.hdus =  [pyfits.PrimaryHDU(header=None),  #primary
                 table,      # this table
                ]
        if self.add_assoc:
            self.hdus += [assoc.cat_ref,]    # the catalog reference (copied)
            
        self.finish(outfile)
        
    def finish(self, outfile):
        pyfits.HDUList(self.hdus).writeto(outfile)
        print '\nwrote FITS file to %s' % outfile
 
if __name__=='__main__':
    m = MakeCat()
            