"""   Analyze localization 

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/localization.py,v 1.17 2018/01/27 15:39:29 burnett Exp $

"""
import os, pickle, collections
import numpy as np
import pylab as plt
import pandas as pd
from uw.utilities import makepivot
from . import (sourceinfo, _html)
from . analysis_base import FloatFormat, html_table
from skymaps import Band, SkyDir
#from . _html import HTMLindex

class Localization(sourceinfo.SourceInfo):
    """Localization summary
    
    <br>Plots summarizing the localization of all sources: precision, quality, closeness, and confusion.
    <br><a href="../sources/index.html?skipDecoration">Back to sources</a>
    """
    require='pickles.zip'
    def setup(self, **kw):
        super(Localization, self).setup(**kw)
        self.plotfolder = 'localization'
        df = kw.get('df', None)
        if df is not None:
            self.df = df
        # unpack the ellipse info into a new DataFrame
        self.ebox = pd.DataFrame([x if x is not None else [np.nan]*7 for x in self.df.ellipse], index=self.df.index)
        self.ebox.columns = 'fit_ra fit_dec a b ang locqual delta_ts'.split()
        self.ebox['roiname']=self.df.roiname
        self.ebox['ts'] = self.df.ts
        
        self.tscut = kw.get('tscut', 10.)
        self.acut =  kw.get('acut', 0.25)
        self.qualcut=kw.get('qualcut', 8.0)
        self.delta_tscut = kw.get('delta_tscut', 2.0)
        poorcut=((self.ebox.locqual>self.qualcut) | (self.ebox.a>self.acut) | 
            (abs(self.ebox.delta_ts)>self.delta_tscut))&(self.df.ts>self.tscut)
        self.poorloc=self.ebox[poorcut] ['ts a locqual delta_ts roiname'.split()].sort_index()
        if len(self.poorloc)>0:
            print '%d poorly localized (locqual>%.1f or a>%.2f or delta_ts>%.2f) '%\
                (len(self.poorloc), self.qualcut,self.acut, self.delta_tscut)
            self.poorloc.to_csv('poorly_localized.csv')
            print 'wrote file "poorly_localized.csv"'
        self.unloc = unloc = np.array(self.df.unloc & (self.df.ts>self.tscut), bool)
        if sum(unloc)>0:
            print '%d point sources (TS>10) without localization information' % sum(unloc)

    def localization(self, maxdelta=9, mints=10, maxqual=5):
        """Localization plots
        The 'finish' stage of creating a model runs the localization code to check that the current position is 
        still appropriate. This is measured by the change in the value of the TS at the best fit position. The position is only 
        updated based on this information at the start of a new series of interations.
           <br> <b>Left</b>: histogram of the square root of the TS difference from current position to
           the fit; corresponds the number of sigmas. <br>
            <b>Right</b>: scatter plot of this vs. TS
            """
        bins=np.linspace(0,np.sqrt(maxdelta),26)
        fig, axx = plt.subplots(1,2,figsize=(13,5)); 
        plt.subplots_adjust(wspace=0.4, left=0.1)
        wp = self.ebox
        cut = self.df.ts>mints
        ax=axx[0]
        for tcut,color in zip((mints, 100), ('blue','orange')):
            t = np.sqrt(wp.delta_ts[(self.df.ts>tcut) & (self.df.locqual<maxqual)].clip(0,maxdelta))
            if  len(t)==0: continue
            ax.hist(t, bins, log=True, color=color, histtype='step', lw=2, label='ts>%d: mean=%.2f'%(tcut, t.mean()) )
        #ax.hist(np.sqrt(wp.delta_ts[self.df.ts>100].clip(0,maxdelta)), bins,label='TS>100\nmean:%f.1'%wp.delta)
        ax.legend(prop=dict(size=10))
        ax.grid()
        plt.setp(ax, xlabel='sqrt(delta TS)', ylim=(0.8,None))
        ax=axx[1]
        ax.plot( self.df.ts[cut],np.sqrt(wp.delta_ts[cut].clip(0,maxdelta)), '.')
        ax.grid()
        plt.setp(ax, xscale='log', xlabel='TS', ylabel='sqrt(delta TS)')
        return fig
        
    def locqual_hist(self, ax=None,  maxqual=10, mints=10, tscut=25, grid_flag=False):
        """Localization quality
            <br>histogram of the fit quality. This is a measure of the difference between the sampled
            TS map points and the prediction of the quadratic model. 

        
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(5,5))
        else:
            fig = ax.figure
        wp = self.ebox
        bins=np.linspace(0,maxqual,26)

        for x, color in zip((mints, tscut), ('blue', 'orange')):
            ax.hist(wp.locqual[self.df.ts>x].clip(0,maxqual), bins,
                histtype='step', lw=2, log=True, color=color, label='TS>%d'%x)
        ax.legend(prop=dict(size=10))
        ax.grid(grid_flag)
        ax.set_ylim(ymin=1)
        plt.setp(ax, xlabel='localization fit quality')
        return fig

    
    def localization_quality(self, maxqual=10, mints=10, tscut=25):
        """Localization quality plots
            <br><b>Left</b>: histogram of the fit quality. This is a measure of the difference between the sampled
            TS map points and the prediction of the quadratic model. <br>
            <br><b>Center</b>: scatter plot of the quality vs. TS. <br>
            <br><b>Right</b>: locations of poorly-fit sources, see the <a href="poorly_localized_table.html?skipDecoration">table</a>.
        """
        bins=np.linspace(0,maxqual,26)
        fig, axxx = self.subplot_array( hsize=(1.0, 0.6, 1.0, 0.2, 2.0, 0.5), figsize=(13,5))
        axx = axxx[0]
        plt.subplots_adjust(wspace=0.4)
        wp = self.ebox
        badfit = pd.isnull(wp.locqual)
        if sum(badfit)>0:
            print 'Warning: {} sources with no fits'.format(sum(badfit))
        cut = (self.df.ts>mints) & (~ badfit)
        ax=axx[0]
        for x in (mints, tscut):
            ax.hist(np.array(wp.locqual[cut & (self.df.ts>x)],float).clip(0,maxqual), bins,label='TS>%d'%x)
        ax.legend(prop=dict(size=10))
        ax.grid()
        plt.setp(ax, xlabel='localization fit quality')
        ax=axx[1]
        ax.plot( self.df.ts[cut],wp.locqual[cut].clip(0,maxqual), '.')
        ax.grid()
        plt.setp(ax, xscale='log', xlim=(10,1e5), xlabel='TS', ylabel='localization fit quality')
        ax=axx[2]
        self.skyplot(self.poorloc.locqual, ax=ax, s=50, vmin=10, vmax=100, cbtext='localization quality')
        return fig 
        
    def r95(self, qualmax=5):
        """ Error circle radius
        R95 is the 95 %%%% containment radius. Here I show the semi-major axis.
        Applying cut quality < %(qualmax)d.
        """
        fig, ax = plt.subplots(1,2, figsize=(13,5))
        plt.subplots_adjust(left=0.1)
        r95 = 60* 2.6 * self.ebox.a[self.ebox.locqual<qualmax]
        ts = self.df.ts[self.ebox.locqual<qualmax]
        self.qualmax = qualmax
        def hist(ax, rmax=30):
            bins = np.logspace(-1,2,31)
            ax.hist(r95, bins, label='all')
            for tsmin in (25,1000):
                ax.hist(r95[self.df.ts>tsmin].clip(0,rmax), bins, label='TS>%d' % tsmin)
            plt.setp(ax, xlabel='R95 (arcmin)', xscale='log')
            ax.grid(); ax.legend(prop=dict(size=10))
        def scat(ax):
            #ax.plot( ts, r95, '.')
            #plt.setp(ax, xscale='log', xlim=(10,1000), ylim=(0,30), xlabel='TS', ylabel='r95')
            ax.plot( ts**-0.5, r95.clip(0,30), '.')
            plt.setp(ax, xlim=(0,0.30),  ylim=(0,30.5), xlabel='1/sqrt(TS)', ylabel='R95 (arcmin)')
            ax.grid()
        for f,a in zip((hist,scat), ax): f(a)
        return fig
        
    def check_closeness(self, tol=0.15, bmin=0, tsmin=10):
        """ Closeness check
            <p>Table of pairs closer than %(close_tol).2f degrees or the sum of the r95 values for each source, but no more than 0.5 deg.
            %(close_table)s
        """
        cut = (np.abs(self.df.glat)>bmin) & (self.df.ts>tsmin)
        indeces = self.df.index[cut] 
        sdirs = self.df[ cut ]['skydir'].values
        r95 = 2.5 * self.df[cut]['a'].values
        r95[r95>0.5]=0.5
        ts = self.df[cut]['ts'].values
        name1=[]; name2=[]; distance=[]; tolerance=[]; roi_index=[]
        ts1=[]; ts2=[]
        for i in range(len(sdirs)-1): 
            a = sdirs[i]
            a95 = r95[i]
            for j in range(i+1, len(sdirs)):
                dist = np.degrees(a.difference(sdirs[j]))
                t = max(a95+r95[j], tol)
                if dist< t:
                    name1.append(indeces[i])
                    name2.append(indeces[j])
                    ts1.append(ts[i])
                    ts2.append(ts[j])
                    roi = Band(12).index(a)
                    roi_index.append(roi)
                    tolerance.append(t)
                    distance.append(dist.round(2))
                    print 'Closer than tolerance: sources %s[%d], %s, %.2f deg < %.2f; ts=%.0f,%.0f' \
                        % (indeces[i], roi, indeces[j], dist,t, ts[i],ts[j])
        self.close_tol = tol
#        def hreftag(name):
#           fn = 'sedfig/' + name.replace(' ','_').replace('+','p') + '_sed_%s.jpg' % self.skymodel
#           if not os.path.exists(fn): return name
#           return '<a href="../../%s">%s</a>' %(fn,name)
#
#        name1href = map(hreftag, name1)
#        name2href = map(hreftag, name2)
        
        self.close_df = tdf= pd.DataFrame(
            dict(source1=name1, source2=name2, distance=distance, 
                    tolerance=tolerance, roi=roi_index, ts1=ts1, ts2=ts2),
                columns = 'source1 ts1 source2 ts2 distance tolerance roi'.split(),
            ).sort_values(by='roi')
        self.close_table = html_table(tdf, 
            name=self.plotfolder+'/close',
            heading='<h4>Table of %d pairs of close sources</h4>'%len(name1),
            float_format=FloatFormat(2), href=False, href_cols=['source1','source2']) 
        return None
        
    def source_confusion(self, bmin=10, dtheta=0.1, nbins=50, deficit_angle=1.0, tsmin=10):
        """ Source Confusion
        Distribution of the distances to the nearest neighbors of all detected sources with |b|> %(bmin)s degrees.
        <br> Left:The number of entries per angular bin divided by the bin's solid angle. The overlaid curve is the expected
        distribution for a uniform distribution of sources with no confusion.
        <br> Right: ratio of measured to expected. 
        <br> Estimated loss: %(loss)s.
        """
        sdirs = self.df[ (np.abs(self.df.glat)>bmin) & (self.df.ts>tsmin) ]['skydir'].values
        closest = ([sorted(map(sdirs[i].difference, sdirs))[1] for i in range(len(sdirs))])
        z = np.degrees(closest)
        self.closest = z
        n = len(z)
        rho = n /(4*np.pi*np.degrees(1)**2) / (1-np.sin(np.radians(bmin)))
        f = lambda x : n*rho * np.exp( -np.pi*rho*x**2)
        
        def deficit(theta):
            n1, n1t =n-sum(z>theta), n*(1-np.exp(-np.pi * rho*theta**2))
            return n1t-n1, n, 100*(n1t-n1)/n
        loss = deficit(deficit_angle)
        self.loss ='%d/%d, or %.1f percent' % loss
        self.bmin = bmin
        print 'lost: ', self.loss
        n += loss[0]
        rho *=(1+loss[0]/n)
        fig, axx = plt.subplots(1,2, figsize=(12,5))
        plt.subplots_adjust(wspace=0.3)
        bins = np.arange(nbins+1)*dtheta
        cumarea=np.cos(np.radians(bins))*2*np.pi*np.degrees(1)**2
        dA = cumarea[:-1]-cumarea[1:]
        h = np.histogram(z,bins)[0]
        x = bins[:-1]+dtheta/2
        
        ax = axx[0]
        ax.errorbar(x, h/dA ,yerr=np.sqrt(h)/dA,  fmt='.', label='TS>%d: %d sources above |b|=%d'%(tsmin, len(sdirs),bmin))
        ax.plot(bins,f(bins), '--g')
        ax.set( yscale='log', ylim=(1,None), 
            ylabel='Number of sources per square degree', xlabel='closest distance (deg)')
        ax.legend(prop=dict(size=10))
        ax.grid()
        
        ax=axx[1]
        ax.errorbar(x, h/dA/f(x), yerr=np.sqrt(h)/dA/f(x),  fmt='o', )
        ax.axhline(1.0, color='k')
        ax.grid()
        ax.set( xlim=(0,3), ylim=(0,1.5),  xlabel='closest distance (deg)',
            ylabel='ratio of detected to expected')
        return fig
    
    def unlocalized(self):
        """Sources without localization
        Check for sources which completely failed attempt at localization.<br>
            %(unlocalized_sources)s
            
        """
        unloc = self.unloc
        if sum(unloc)>0:
            unloc_table = self.df.ix[unloc]['ra dec ts roiname'.split()].sort_values(by='roiname')
            self.unlocalized_sources =html_table(unloc_table,
                name=self.plotfolder+'/unlocalized',
                heading='<h4>Table of %d unlocalized sources</h4>'%len(unloc_table),
                href_pattern='tsmap_fail/%s*_tsmap*jpg',
                float_format=FloatFormat(2))
        else:
            self.unlocalized_sources='<p>No unlocalized sources'
        
    def poor_loc(self):
        """ Poorly localized sources
                %(poorly_localized_table_check)s
        """
        if len(self.poorloc)>0:
            self.poorly_localized_table_check  = html_table(self.poorloc, 
                name=self.plotfolder+'/poorly_localized',
                heading='<h4>Table of %d poorly localized (a>%.2f deg, or qual>%.1f with TS>%d) sources</h4>'\
                        % ( len(self.poorloc),self.acut,self.qualcut, self.tscut),
                float_format=FloatFormat(2))
                
            #poorly_localized_tablepath = os.path.join(self.plotfolder,'poorly_localized_table.html')
            #open('poorly_localized_table.html','w').write(tohtml)
            #print 'Wrote poorly_localized_table.html'
            #open(os.path.join(poorly_localized_tablepath),'w').write(
            #    '<head>\n'  + _html.style + '</head>\n<body>\n<h3>Poorly Localized Source Table</h3>'\
            #                +  tohtml+'\n</body>')
            #print 'saved html doc to %s' % os.path.join(poorly_localized_tablepath)
            #self.poorly_localized_table_check =\
            #            '<p><a href="%s?skipDecoration"> Table of %d poorly localized (a>%.2f deg, or qual>%.1f with TS>%d) sources</a>'\
            #            % ( 'poorly_localized_table.html',len(self.poorloc),self.acut,self.qualcut, self.tscut)
            try:
                version = os.path.split(os.getcwd())[-1]
                pv = makepivot.MakeCollection('poor localizations %s'%version, 'tsmap_fail', 'poorly_localized.csv',refresh=True)
                self.poorly_localized_table_check +=\
                 '<br>A  <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">pivot collection </a>of TS maps for these sources can be examined.'%pv.cId 
            except Exception, msg:
                self.poorly_localized_table_check += '<br>(No pivot table: %s)' %msg
                print '**** Failed to create pivot: %s' % msg
                        
        else:
            self.poorly_localized_table_check ='<p>No poorly localized sources!'

        
    def load_moment_analysis(self, make_collection=True):
        """ check results of moment analysis
        """
        m =self.df.moment
        has_moment = [x is not None for x in m]
        print 'Found %d sources with moment analysls' % sum(has_moment)
        mdf = pd.DataFrame(m[has_moment]) 
        u = np.array([list(x) for x in mdf.moment])
        self.dfm=md= pd.DataFrame(u,  index=mdf.index, columns='rax decx ax bx angx size peak_fract'.split())
        md['locqual'] = self.df.locqual
        md['ts'] = self.df.ts
        md['delta_ts'] = self.df.delta_ts
        md['roiname'] = self.df.roiname
        md['a'] = self.df.a
        
        # generate the angular difference of fit vs. current positions 
        delta=[]
        for n,t in self.df.iterrows():
            ellipse = t['ellipse']
            delta.append(np.degrees(SkyDir(*ellipse[:2]).difference(t['skydir'])) if ellipse is not None else None )
        self.df['delta'] = delta
        md['delta']= self.df.delta
        
        filename = 'moment_localizations.csv'
        md.to_csv(filename)
        print 'Write file %s' % filename
        if not make_collection: return md
        
        # now make a collection with images and data
        self.moment_collection_html=''
        try:
            t = makepivot.MakeCollection('moment analysis localizations %s'%self.skymodel, 'tsmap_fail', 'moment_localizations.csv', 
                refresh=True) 
            makepivot.set_format(t.cId)
            self.moment_collection_html="""
                <p>The images and associated values can be examined with a 
                <a href="http://deeptalk.phys.washington.edu/PivotWeb/SLViewer.html?cID=%d">Pivot browser</a>,
                which requires Silverlight."""  % t.cId
        except Exception, msg: 
            print "**** Failed to make moment pivot table: %s" % msg

        

    def moment_plots(self):
        """Plots of properties of the moment analysis
        This analysis was done on localizations that had  locqual>5, or a>0.25, or delta_ts>2. 
        It invovled creating a 15x15 TS plot with size estimated from the major axis, but limited to 2 degrees.
        
        The values shown are:
        <ol><li>peak fraction: fraction of the total weight for the largest bin. This is a check on the scale.
            Small means a uniform distribution, large means scale is too large </li>
            <li>major: Size of major axis (1 sigma)</li>
            <li>major/size : this ratio should be neither small, not >1 for the fit to make sense</li>
            <li>minor/major: Ratio of minor to major axis size</li>
            <li>size: Size of image. The analysis that generated these images used the peak fitting analysis to set 
            the scale, but limited to 2 degrees.</li>
            <li>locqual : The peak-finding localization quality. This being too large is the primary reason these sources
            were selected. </li>
         </ol> 
        <p>Subsets correspond to peak fraction betweek 0.1 and 0.5.
         %(moment_collection_html)s
        """
        
        fig, axx = plt.subplots(2,3, figsize=(10,8))
        plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.1)
        if not hasattr(self,'dfm'):
            self.load_moment_analysis()
        dfs = self.dfm
        goodfrac= (dfs.peak_fract<0.5) & (dfs.peak_fract>0.1)
        for ix, ax in enumerate(axx.flatten()):
            if ix==0:                       
                ax.hist(dfs.peak_fract, np.linspace(0,1,26))
                ax.hist(dfs.peak_fract[goodfrac], np.linspace(0,1,26))
                plt.setp(ax, xlabel='peak fraction', xlim=(0,1))
            elif ix==2:
                ax.hist(dfs.ax/dfs['size'], np.linspace(0,0.5,26))
                ax.hist((dfs.ax/dfs['size'])[goodfrac], np.linspace(0,0.5,26))
                plt.setp(ax, xlabel='major/size')
            elif ix==1:
                ax.hist(dfs.ax.clip_upper(0.5), np.linspace(0,0.5,26))
                ax.hist(dfs.ax[goodfrac].clip_upper(0.5), np.linspace(0,0.5,26))
                plt.setp(ax, xlabel='major')
            elif ix==3:
                ax.hist(dfs.bx/dfs.ax, np.linspace(0,1,26))
                ax.hist((dfs.bx/dfs.ax)[goodfrac], np.linspace(0,1,26))
                plt.setp(ax, xlabel='minor/major')
            elif ix==4:
                ax.hist(dfs['size'], np.linspace(0,2,26))
                ax.hist(dfs['size'][goodfrac], np.linspace(0,2,26))
                plt.setp(ax, xlabel='size')
            elif ix==5:
                ax.hist(np.array(dfs.locqual.clip_upper(8),float), np.linspace(0,8))
                plt.setp(ax, xlabel='locqual')
        return fig
    
    
    def all_plots(self):
        return self.runfigures([self.r95, self.localization,self.localization_quality,
            self.unlocalized, 
            #self.poor_loc,
            self.moment_plots,
            self.check_closeness,self.source_confusion,
            ])
