"""
Association analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/associations.py,v 1.15 2014/02/15 23:06:33 burnett Exp $

"""
import os, glob, sys, pyfits
import numpy as np
import pylab as plt
import pandas as pd

from skymaps import SkyDir, Band
from . import sourceinfo
from . analysis_base import FloatFormat, html_table


class Associations(sourceinfo.SourceInfo):
    """Analysis of associations
    <p>
    A feature of the UW pipeline is the application of the LAT association algorithm, with the same catalogs as
    were used for 2FGL, except that the latest LAT pulsar catalog is used. Note that priors are not recalculated,
    we use the values determined for 2FGL. The method is explained in the <a href="http://arxiv.org/pdf/1002.2280v1.pdf">1FGL paper</a>, see Appendix B.
    <p>
    <p>Note that Jurgen's pipeline now has the Planck and WISE catalogs.
    """
    
    def setup(self, **kw):
        self.args = kw.pop('args', None)
        super(Associations, self).setup(**kw)
        self.plotfolder='associations'
        self.load_assoc(self.args)
    
    def load_assoc(self, fromdf=None):
        if fromdf is not None:
            print 'loading associations from file %s' %fromdf
            self.df['altassoc']=pd.load(fromdf)
        else: print 'using associations found in sourceinfo'
        associations = self.df.associations if fromdf is None else self.df.altassoc
        probfun = lambda x: x['prob'][0] if not pd.isnull(x) else 0
        self.df['aprob'] = np.array([ probfun(assoc) for  assoc in associations])
        self.df['acat']  = np.array([ assoc['cat'][0] if not pd.isnull(assoc) else 'unid' for  assoc in associations])
        self.df['aname'] = np.array([ assoc['name'][0] if not pd.isnull(assoc) else 'unid' for  assoc in associations])
        self.df['aang']  = np.array([ assoc['ang'][0] if not pd.isnull(assoc) else np.nan for  assoc in associations])

        self.df['adeltats'] = np.array([assoc['deltats'][0] if not pd.isnull(assoc) else np.nan for assoc in associations])
        
        self.df10 = self.df.ix[self.df.ts>10]
        print 'associated: %d/%d' % (sum(self.df10.aprob>0.8), len(self.df10))
        
    def association_vs_ts(self, aprob_min=0.5):
        """ Associations vs. TS
        
        <br>Left: histogram of TS, showing the fraction that have associations.
        <br>Right: The fractions themselves.
        """
        ts = self.df10.ts
        lowlat = np.abs(self.df10.glat)<5
        assoc = self.df.aprob>aprob_min
        def plota(ax, bins=np.logspace(1,5,41) ):
            ax.hist(ts, bins, label='all sources')
            ax.hist(ts[assoc], bins, color='orange', label='associated')
            plt.setp(ax, xscale='log', xlabel='TS', xlim=(10,1e5))
            ax.legend(prop=dict(size=10)); ax.grid()
        def plotb(ax, bins=np.logspace(1,4.5,8)):
            for tsvals, label,color in zip( (self.df10[~lowlat].ts, self.df10[lowlat].ts), 
                    ('|b|>5', '|b|<5'), ('blue','red')):
                all = np.array(np.histogram(tsvals, bins)[0],float)
                subset = np.histogram(tsvals[assoc], bins)[0]
                fraction = (subset/all)
                x = np.sqrt(bins[:-1]*bins[1:])
                yerr = np.sqrt(subset*(all-subset)/all )/all
                xerr = [x-bins[:-1], bins[1:]-x]
                ax.errorbar(x=x, y=fraction,xerr=xerr, yerr=yerr, fmt= 'o', color=color, label=label)
            plt.setp(ax, xscale='log', xlim=(bins[0],bins[-1]), ylim=(0,1), xlabel='TS', ylabel='associated fraction')
            ax.grid()
            ax.legend(loc='upper left', prop=dict(size=10))
            
        fig, axx = plt.subplots(1,2, figsize=(12,5))
        plt.subplots_adjust(left=0.1)
        for f,ax in zip((plota,plotb), axx.flatten()): f(ax) 
        return fig
            
    def summary(self):
        """Summary
        %(summary_html)s
        """
        # load the catalogs used
        srcid_path=sys.path[0] = os.path.expandvars('$FERMI/catalog/srcid/')
        import classes
        cats = dict()
        for module_name in classes.__all__:
            cd = dict()
            # this is equivalent to "from classes import module_name"
            module = __import__('classes.'+module_name,  fromlist=['classes'])

            for var in 'catid catname prob_prior prob_thres figure_of_merit max_counterparts new_quantity selection'.split():
                cd[var] = module.__dict__[var] 
            cats[module_name]= cd
            if cd['catname'].find('*')>0:
                try:
                    t = glob.glob(srcid_path+'/cat/'+cd['catname'])[-1]
                    cd['catname']=os.path.split(t)[-1]
                except:
                    print 'File %s not found' %s
            
        self.catdf = pd.DataFrame(cats).T
        self.catdf['objects']=[len(pyfits.open(srcid_path+'cat/'+fname)[1].data) for fname in self.catdf.catname]
        # make dict of catnames with values the number of associations
        t = dict()
        for c in self.df10.acat:
            if c not in t: t[c]=1
            else: t[c]+=1
        
        self.catdf['associations'] = pd.DataFrame(t.items(), index=t.keys(), columns='name associations'.split())['associations']
        
        self.summary_html = html_table(
            self.catdf[~pd.isnull(self.catdf.associations)]['objects associations prob_prior prob_thres catname'.split()],
            dict(objects='Objects,number of entries in table',
                 prob_prior='Prior,prior probability', prob_thres='Threshold,threshold probability',
                 catname='Catalog,name of the FITS file in the folder %s'%(srcid_path+'/cat')),
            heading='<br>Table of catalogs with associations', name=self.plotfolder+'/associated_catalogs', 
            href=False, maxlines=100)
                 
    def pulsar_check(self):
        """LAT pulsar check
        %(atable)s
        """
        # compare with LAT pulsar catalog     
        tt = set(self.df.name[self.df.psr])
        pulsar_lat_catname = sorted(glob.glob(os.path.expandvars('$FERMI/catalog/srcid/cat/obj-pulsar-lat_*')))[-1]
        print 'opening LAT catalog file %s' %pulsar_lat_catname
        pp = pyfits.open(pulsar_lat_catname)[1].data
        lat = pd.DataFrame(pp, index=[n.strip() for n in pp.Source_Name])
        lat['ts'] = self.df[self.df.psr]['ts']
        lat['ROI_index'] = [Band(12).index(SkyDir(float(ra),float(dec))) for ra,dec in zip(lat.RAJ2000,lat.DEJ2000)]
        
        lat['skydir'] = [SkyDir(float(ra),float(dec)) for ra,dec in zip(lat.RAJ2000, lat.DEJ2000)]
        #map(SkyDir, np.array(lat.RAJ2000,float), np.array(lat.DEJ2000,float))
        lat['sourcedir'] = self.df.skydir[self.df.psr]
        lat['delta'] = [np.degrees(s.difference(t)) if not type(t)==float else np.nan for s,t in zip(lat.skydir,lat.sourcedir)]
        if sum(lat.delta>0.25)>0:
            print 'LAT pulsar catalog entries found more than 0.25 deg from catalog:'
            #print lat[lat.delta>0.25]['ts delta'.split()]
        dc2names =set(pp.Source_Name)
        print 'sources with exp cutoff not in LAT catalog:', list(tt.difference(dc2names))
        print 'Catalog entries not found:', list(dc2names.difference(tt))
        missing = np.array([ np.isnan(x) or x<10. for x in lat.ts])
        
        # this used to work but now generates 'endian' message
        #latsel = lat[missing]['RAJ2000 DEJ2000 ts ROI_index'.split()]
        missing_names = lat.index[missing]
        cols = 'RAJ2000 DEJ2000 ts ROI_index'.split()
        latsel = pd.DataFrame( np.array([lat[id][missing] for id in cols]), index=cols, columns=missing_names).T

        psrx = np.array([x in 'pulsar_fom pulsar_low msp pulsar_big'.split() for x in self.df.acat])
        print '%d sources found in other pulsar catalogs' % sum(psrx)

        self.atable = '<h4>Compare with LAT pulsar catalog: %s</h4>' % os.path.split(pulsar_lat_catname)[-1]
        self.atable += '<p>Sources fit with exponential cutoff not in catalog %s' %list(tt.difference(dc2names))
        self.atable += html_table(latsel,
                    dict(ts='TS,Test Statistic', ROI_index='ROI Index,Index of the ROI, a HEALPix ring index'),
                    heading = '<p>%d LAT catalog entries not in the model (TS shown as NaN), or too weak.' % sum(missing),
                    name=self.plotfolder+'/missing', maxlines=20,
                    float_format=(FloatFormat(2)))
        far = lat.delta>0.25
        if sum(far)>0:
            far_names = lat.index[far]
            cols = 'ts delta'.split()
            latfar = pd.DataFrame( np.array([lat[id][far] for id in cols]),index=cols, columns=far_names).T
            self.atable += '<p>Pulsars located > 0.25 deg from nominal'\
                    + latfar.to_html(float_format=FloatFormat(2))
        if sum(psrx)>0:
            self.atable+= html_table(self.df[psrx]['aprob acat aname aang ts delta_ts locqual'.split()],
                          dict(name='Source Name,click for link to SED',
                          ts='TS,Test Statistic for the source', 
                          acat='catalog,Catalog nickname',
                          aprob='Probability,Association probability',
                          aname='Source Name,Catlog name for the source',
                          aang='Angle,distance to the catalog source (deg)',
                          delta_ts='Delta TS,change in TS to the catalog source\n'
                                          'should be positive negative means peak of TS map was not at source',
                          locqual='Localization quality,measure of the goodness of the localization fit\n greater than 5 is questionable',
                          ),
                          float_format=FloatFormat(2), 
                          heading = """<p>%d sources with pulsar association not in LAT pulsar catalog. Note, no cut on 
                            association probability.""" % sum(psrx),
                          name=self.plotfolder+'/atable',
                          maxlines=50)        
        
    def localization_check(self, tsmin=10, dtsmax=9, qualmax=5):
        r"""Localization resolution test
        
        The association procedure records the likelihood ratio for consistency of the associated location with the 
        fit location, expressed as a TS, or the difference in the TS for source at the maximum, and at the associated
        source. The distribution in this quantity should be an exponential, $\exp(-TS/2/f^2)$, where $f$ is a scale factor
        to be measured from the distribution. If the PSF is a faithful representation of the distribution of photons
        from a point source, $f=1$. For 1FGL and 2FGL we assumed 1.1. The plots show the results for AGN, LAT pulsars, and
        all other associations. They are cut off at 9, corresponding to 95 percent containment.
        """
        
        t = self.df.acat
        agn = np.array([x in 'crates bzcat agn bllac'.split() for x in t])
        psr = np.array([x in 'pulsar_lat'.split() for x in t])
        unid= np.array([x in 'unid'.split() for x in t])
        otherid=~(agn | psr | unid)

        def select(sel,  df = self.df, tsmin=tsmin, qualmax=qualmax):
            cut = sel * (df.aprob>0.8) * (df.ts>tsmin) * (df.locqual<qualmax)
            return df[cut].adeltats
            
        fig, axx = plt.subplots(1,3, figsize=(14,5))

        for sel, name, ax in zip((agn, psr,otherid), ('AGN','LAT pulsars', 'other ids'), axx):
            z = FitExponential(select(sel), name, vmax=dtsmax)
            z.plot(ax, xlabel=r'$\Delta TS$')
            print '%s: localization factor=%.2f' %(name, z.factor)

    def all_plots(self):    
        self.runfigures([self.summary, self.pulsar_check, self.association_vs_ts, self.localization_check,])


class FitExponential(object):
    """manage the fit to an exponential distribution, with no background
    """
    
    def __init__(self, v, label, vmax=9,  binsize=0.25):
        from scipy import optimize

        self.vmax, self.binsize, self.label = vmax, binsize, label
        self.vcut=vcut = v[v<vmax]
        self.vmean = vmean = vcut.mean() 
        # find factor that has same average over the interval
        self.factor = optimize.brentq( lambda x : self.cfactors(x)[3]-self.vmean, 1.0, 1.5)
        beta, c0, c1, r = self.cfactors(self.factor)
        self.alpha = len(vcut) / c0 * binsize
        self.beta=beta
        
    def cfactors(self, f):
        """for factor f, return beta, c0, c1, and c1/c0
        where c0 and c1 are integrals(o,vmax) of exp(-x/beta) and x*exp(-x/beta)
        Thus c1/c0 is <x> over the interval
        """
        beta = 2 * f**2
        u = self.vmax/beta
        c0 = beta * (1 - np.exp(-u))
        c1 = beta**2 * (1 - np.exp(-u) * (1+u)) 
        return beta, c0, c1, c1/c0
    
    def plot(self, ax=None, xlabel=''):
        if ax is None:
            fig,ax = plt.subplots(figsize=(5,5))
        else:fig = ax.figure
        x = np.linspace(0, self.vmax, int(self.vmax/self.binsize)+1) 
        ax.hist( self.vcut, x, label='%d %s'%(len(self.vcut),self.label))
        ax.set_ylim(ymin=0)
        ax.plot(x, self(x), '-r', lw=2,  label='factor=%.2f'% self.factor)
        ax.grid(); ax.legend()
        plt.setp(ax, ylim=(0, 1.1*self.alpha), xlabel=xlabel)

    def __call__(self, x):
        return self.alpha * np.exp(-x/self.beta)

