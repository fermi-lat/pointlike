"""
Manage Fermi catalogs
$Header:$

"""
import os, cStringIO, glob
import numpy as np
import pandas as pd
from astropy.io import fits
from skymaps import SkyDir, Band
from uw.like import Models
from uw.like2 import sedfuns

# this gets set to most recent Fermi Calalog 
energy_bounds=None

class GLL_PSC(object):
    """Parse a Fermi Catalog
    """
    def __init__(self, filename='$FERMI/catalog/gll_psc7year10GeV_v7r1.fit'):
        global energy_bounds
        if not os.path.exists(os.path.expandvars(filename)):
            check = glob.glob(os.path.expandvars('$FERMI/catalog/')+filename+'_*.fit*')
            if len(check)>0: 
                filename = check[-1]
            else:
                raise Exception('could not resolve filename or catalog name {}'.format(filename))
        self.version = filename.split('_')[-1].split('.')[0]
        print 'Loaded {}'.format(filename)
        self.hdulist = fits.open(os.path.expandvars(filename))
        energy_bounds = pd.DataFrame(self.hdulist['EnergyBounds'].data)
        self.pscdata = self.hdulist['LAT_Point_Source_Catalog'].data
        field = self.pscdata.field
        self.colnames = [x.name for x in self.pscdata.columns]
        if 'Source_Name' in self.colnames:
            self.names = field('Source_Name')
        else:
            self.names=[]
        self.nicknames=field('NickName')
        self.skydirs = map(SkyDir, np.array(field('RAJ2000'),float),
              np.array(field('DEJ2000'),float))
        
    def __repr__(self):
        return self.hdulist.info()

    def __getitem__(self, name):
        """Return an entry as a pandas Series object
        name : int or string
            if int, the row number. If string, the source or nick name
        """
        if type(name)==int:
            i = name
        else:
            if name in self.names:# 
                i = list(self.names).index(name)
            elif name in self.nicknames:
                i = list(self.nicknames).index(name)
            else:
                raise Exception(
                    'Name "{}" not found in {}'.format(name, (self.names,self.nicknames)))
        return pd.Series(data=self.pscdata[i], index=self.colnames)
        
    def table(self):
        """
        This takes a very long time, since each Series is expensive
        """
        u = []
        for i in range(len(self.names)):
            print i
            u.append(self[i])
        return pd.DataFrame(u)

    def data_frame(self):
        """ Return a DataFrame with a simplified subset, including pointlike Model objects

        """
        field_names = ('RAJ2000 DEJ2000 GLAT GLON Test_Statistic Npred Conf_95_SemiMajor Pivot_Energy'
            +' Flux_Density Spectral_Index beta').split()
        # make columns, with type either float, or str
        col_data = [np.array(self.pscdata.field(fname),float) for fname in field_names]
        col_data.append(np.array(self.pscdata.field('SpectrumType'),str)) 

        # now put it into a DataFrame with simple names, index with NickName
        cols = 'ra dec glat glon ts npred r95 pivot_energy flux pindex beta spectral_type'.split()
        df = pd.DataFrame(col_data, index=cols).T
        df.index=self.pscdata.field('NickName'); 
        df.index.name='name'
        df['cutoff'] = np.nan # make dependent on catalog type
 
        # add Pointlike stuff
        df['skydir'] = map(lambda ra,dec:SkyDir(float(ra),float(dec)), df.ra, df.dec)
        df['roi'] = map(Band(12).index, df.skydir)
        models = []

        def model(st, flux, index, beta, pivot, cutoff):
            if st=='PowerLaw':
               return Models.PowerLaw(p=[flux, index], e0=pivot)
            elif st=='PLSuperExpCutoff':
                prefactor = flux*np.exp((pivot/cutoff)**b)
                return Models.PLSuperExpCutoff(p=[prefactor, index, cutoff,b], e0=pivot)
            elif st=='LogParabola':
                return Models.LogParabola(p=[flux, index, beta, pivot])
            elif st=='PowerLaw2':   ### same as PowerLaw in table
                return Models.PowerLaw(p=[flux, index], e0=pivot)
            else:
                raise Exception('unexpected spectrum type {}, source name {}'.format(st,name))

        for name, row in df.iterrows():
            models.append(
                model( *[row[q] for q in 'spectral_type flux pindex beta pivot_energy cutoff'.split()])
                )
        df['model']=models
        return df
     
 
    def band_flux(self, name):
        """
        """
        return BandFlux(self[name])

    def __repr__(self):
        output = cStringIO.StringIO()
        self.hdulist.info(output)
        return 'FITS catalog:\n{}'.format(output.getvalue())

class BandFlux(object):
    
    @staticmethod
    def conversion_factors(e1,e2,gamma):
        """ Return array of conversion factors from photon flux per energy bin to 
            differential energy flux, in eV units

            e1,e2: array of float
                upper, lower band limits in MeV
            gamma : array of float
                spectral index 
        """
        assert len(gamma)==len(e1), 'Wrong list of energies'
        c = (gamma-1)/(e1**(1-gamma)-e2**(1-gamma))
        return c*(e1*e2)**(1-gamma/2) * 1e6
 

    def __init__(self, cat_entry):
        """return a Pandas.DataFrame for the band flux information
        cat_entry: pandas.Series object
        
        Adds columns with conversion to differential energy flux, in eV units
        """
        cols = 'Flux Unc_Flux Index Npred Sqrt_TS'.split()
        t=[cat_entry[cname+'_Band'] for cname in cols]  
        bf =pd.DataFrame(t, index=cols).T
        
        # add columns with energy flux
        # evaluate energy flux = e**2 * df/de at geometric mean, assuming powerlaw
        # convert to units eV/(cm**2 * s)
        eb=energy_bounds
        r = BandFlux.conversion_factors(eb.LowerEnergy, eb.UpperEnergy, bf.Index)
        bf['eflux'] = bf.Flux * r 
        bf['eflux_unc'] = (bf.Unc_Flux) * r
        self.bf=bf
       # construct the spectral model object if info is there
       # TODO: deal with exponential
        s=cat_entry
        #self.name = cat_entry['Source_Name']
        self.model=None
        self.nickname = cat_entry['NickName']
        if 'Flux_Density' in cat_entry:
            if s['SpectrumType']=='LogParabola':
                self.model = Models.LogParabola(
                    [s['Flux_Density'], s['Spectral_Index'], s['beta'], s['Pivot_Energy']],
                    free=[True,True,False,False])

    def plot(self, ax=None, ylim=(1e-2,1e2), **kwargs):
        """Make a plot with data and spectral function
        """
        from matplotlib import pyplot as plt
        ul_kwargs = kwargs.copy()
        ul_kwargs['color']='gray'
        if 'color' not in kwargs:
            kwargs['color'] = 'k'
        ul_kwargs['color']='gray'

        bf = self.bf
        eb = energy_bounds # from global set above by class
        if ax is None:
            fig, ax = plt.subplots(1,1, figsize=(5,5))
        else: fig = ax.figure
        xc = np.sqrt(eb.LowerEnergy * eb.UpperEnergy)
        yc = np.array(bf.eflux, float)
        yerr = np.array([abs(t) for t in bf.eflux_unc]).T
        xerr =(xc-eb.LowerEnergy, eb.UpperEnergy-xc)
        for i in range(len(eb)):
            xl,xh = eb.LowerEnergy[i], eb.UpperEnergy[i]
            bc = xc[i]
            f, df = bf.eflux[i], bf.eflux_unc[i]
            if f>1e-2: #essentially zero
                ax.plot([xl,xh], [f,f],  **kwargs)
                ax.plot([bc,bc], f+df,  **kwargs)
            else:
                x,y = bc, 2*(f+df[1])
                ax.plot([xl,xh], [y,y] , **ul_kwargs) # bar at upper limit
                # plot arrow 0.6 long by 0.4 wide, triangular head (in log coords)
                ax.plot([x, x,     x*1.2, x,     x/1.2, x],
                        [y, y*0.6, y*0.6, y*0.4, y*0.6, y*0.6], **ul_kwargs)
    
        # overplot the function
        dom = np.logspace(np.log10(eb.LowerEnergy[0]),np.log10(list(eb.UpperEnergy)[-1]))
        if self.model is not None:
            ax.plot(dom, dom**2 * self.model(dom)*1e6, color='red', lw=2, ) 

        ax.set_title(self.nickname)
        plt.setp(ax, xlabel='Energy [MeV]', xscale='log', 
                yscale='log', ylim=(1e-1,None) if ylim is None else ylim,
                xlim=(None, 2e6));
        ax.set_ylabel(r'$\mathsf{Energy\ Flux\ (%s\ cm^{-2}\ s^{-1})}$' % 'eV', labelpad=0)

        plt.grid(alpha=0.5)
        fig.set_facecolor('white')
        return fig
       
    def table(self):
        """return a simplified DataFrame
        """
        bf = self.bf.copy()
        bf['emin']=energy_bounds.LowerEnergy
        bf['emax']=energy_bounds.UpperEnergy
        bf['TS']= bf.Sqrt_TS**2
        return bf['emin emax Flux Npred TS eflux eflux_unc'.split()]
    
class CreateFermiFITS(object):
    def __init__(self, df, names = None, filename='sed_info.fits'):
        """
        df : DataFrame
            contains all analysis
            Will extract the SED info for each source in the dataframe
        names : list of string
            If set and filename is present, use this subset of the df entries
        filename : string | None
            File to write to. If filename is None
        """
        self.df=df
        # get an entry from the source list to extract elow, ehigh
        # (assume that all are the same!)
        sr = df.ix[0]['sedrec']
        self.elow = sr['elow']; 
        self.ehigh = sr['ehigh']

        if filename is None: return

        # making a stand-alone FITS file. Add the names
        names = list(df.index) if names is None else names
        dcols = [fits.Column(name='NickName', format = '18A', array=names)]

        self.hdulist =fits.HDUList(hdus=[
            fits.PrimaryHDU(header=None),
            fits.BinTableHDU.from_columns(dcols+self.make_bands(names),
                 name='LAT_Point_Source_Catalog'),
            fits.BinTableHDU.from_columns(self.make_ebounds(), name='EnergyBounds'),
            ], 
            ) 
        self.hdulist.writeto(open(filename,'w'))
        print 'Created FITS file {} with {} sources'.format(filename, len(names))

 
    def make_Band(self,sr):
        Index_Band = sr['pindex']
        Npred_Band = sr['npred']
        r = BandFlux.conversion_factors(self.elow, self.ehigh, Index_Band)
        f, u, l = sr['flux'], sr['uflux'], sr['lflux'] 
        nuFnu = f * 1e-6  # convert to MeV
        Flux_Band = f / r
        Unc_Flux_Band = np.array([(l-f)/r,(u-f)/r]).T 
        Sqrt_TS_Band = np.sqrt(sr['ts'])
        return pd.DataFrame([nuFnu, Flux_Band, Unc_Flux_Band, Npred_Band, Index_Band, Sqrt_TS_Band],
                    index='nuFnu Flux_Band Unc_Flux_Band Npred_Band Index_Band Sqrt_TS_Band'.split()).T
    
    def make_bands(self, names):
        """ return a list of fits.Columns for the bands
        """

        nufnu = [];  fb = []; ufb = [];  npb =[];  ib = [];   sts=[]
        for i,name in enumerate(names):
            if name not in self.df.index:
                print 'Warning: name {} not found'.format(name)
                continue
            sr = self.df.ix[name]['sedrec']
            if sr is None: 
                print 'No sedrec for {} at with TS={:.0f}'.format(name,  row['roiname'],row['ts'],)
                continue
            sr=self.make_Band(sr)
            nufnu.append(np.array(sr['nuFnu'],float))
            fb.append(np.array(sr['Flux_Band'],float))
            ufb.append(np.array(list(sr['Unc_Flux_Band']),float))
            npb.append(np.array(sr['Npred_Band'],float))
            ib.append(np.array(sr['Index_Band'],float))
            sts.append(np.array(sr['Sqrt_TS_Band'],float))

        nb = len(self.elow) # number of bands
        nbE = '{}E'.format(nb)
        dcols = [
            fits.Column(name='Flux_Band',   format=nbE, unit='photon/cm**2/s', array=fb),
            fits.Column(name='Unc_Flux_Band', format='{}E'.format(2*nb), unit='photon/cm**2/s', 
                            dim = '(2,{})'.format(nb),    array=ufb),
            fits.Column(name='nuFnu',       format=nbE, unit='MeV/cm**2/s',array=nufnu),
            fits.Column(name='Index_Band',  format=nbE, array=ib),
            fits.Column(name='Npred_Band',  format=nbE, array=npb),
            fits.Column(name='Sqrt_TS_Band',format=nbE, array=sts),
            fits.Column(name='Spectral_Fit_Quality', format='E'),
            ]
        return dcols

    def make_ebounds(self):
        ecols = [
            fits.Column(name='LowerEnergy', format = 'E', unit = 'MeV', array=self.elow),
            fits.Column(name='UpperEnergy', format = 'E', unit = 'MeV', array=self.ehigh),
            fits.Column(name='SysRel',      format = 'E')
            ]
        return ecols
    
        