"""
Create a table of source fits as a confluence table

"""
class RegFile(object):
##    class Source(object):
##        def __init__(self, source):
##            tokens = line.split()
##            self.name = tokens[0]
##            self.sigma, self.ra, self.dec = [float(t) for t in tokens[2:5] ]
    def __init__(self, data):
        self.sourcelist = data
        print 'found %d sources' % len(self.sourcelist)
    def __getitem__(self, i): return self.sourcelist[i]
    def __len__(self): return len(self.sourcelist)

    def write_regfile(self, filename='test.reg', color='red'):
        header='global color=%s font="helvetica 10 normal" select=1 edit=1 move=0 delete=0 include=1 fixed=0 width=2;fk5;'%color
        out = file(filename, 'w')
        print >> out, header
        for s in self.sourcelist:
            print >> out,\
                  'point(%5.3f, %5.3f) # point=cross 20  text={TS=%5.0f};'\
                      % ( s.ra, s.dec, s.TS)
        out.close()



class Table(object):
    def sinbad_ref(self, ra, dec):
        sd = '%2b' if dec>0 else '%2d'
        ad = abs(dec)
        coord = '%.4f%%20%s%.4f%%20'%(ra,sd, ad)
        
        return 'http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord=%s&submit=submit%%20query&Radius.unit=deg&CooEqui=2000&CooFrame=FK5&Radius=0.2'%coord
        #1.6362%20%2b73.0255%20
    def link_ref(self, ra):
        return '^RA%03d.png' % int(ra)
    
    def __init__(self, path, version, verbose=False, link_sed=True):
        self.version=version
        self.verbose=verbose
        self.link_sed = link_sed
        self.infilename=r'%spointfit_%s.txt' %(path, version)
        self.outfilename='%sconfluence_%s.txt' %(path, version)
        infile =file(self.infilename)
        self.outfile = file(self.outfilename, 'w')
        titles = infile.readline()[2:] # split off initial #?
        self.cols = len(titles.split())
        self.title(titles)
        for line in infile:
            self.entry(line)
        self.outfile.close()
        self.histogram()

    def title(self, line):
        q='|| '
        for i,t in enumerate(line.split()):
            if t in ('ra','dec'): t+=r'\\(deg)'
            elif t in ('localization','moved'): t+=r'\\(deg)'
            elif t=='[neighbor]': break
            q+= t+' || '
        if self.verbose: print q
        self.outfile.write(q+'\n')
        
    def entry(self, line):
        toks = line.split()
        ra,dec = [float(t) for t in toks[1:3]]
        q = '|'
        for i,t in enumerate(toks):
            if i==0 :
                if self.link_sed:
                  q+= '[%s| %s] [S|%s] | ' % (t, self.link_ref(ra), self.sinbad_ref(ra,dec))
                else:
                  q+= '%s [S|%s] | ' % (t, self.sinbad_ref(ra,dec))
            else: q+= t + ' | '
        #for i in range(self.cols-len(toks)+1): q+= ' |' # note extra
            
        if self.verbose: print q
        self.outfile.write(q+'\n')

    def histogram(self):
        import matplotlib
        from pylab import figure, hist, grid, plot, axvline, title, xlabel,ylabel, text, savefig
        from numpy import arange,array,exp
        t =matplotlib.mlab.csv2rec(self.infilename, delimiter=' ')
        figure(figsize=(4,4))
        ratio=[s[5] for s in t if s[4]<0.25] # ratio if localization < 0.25
        binsize=0.25
        xmax=5.0
        hist(ratio  , arange(0,xmax,binsize))
        grid()
        x = arange(0,xmax,0.1)
        y = len(ratio)*x*exp(-x**2/2)*binsize
        plot(x,y, 'r-', lw=2)
        ylabel('sources/%3.2f'%binsize)
        title('resolution check')
        xlabel('difference/sigma')
        axvline(2.45); text(2.5, 4.5, '95% error circle')
        outfile = self.infilename.replace('.txt', '.png')
        savefig(outfile)
        print 'wrote png histogram to %s' % outfile


if __name__=='__main__':
    analysis_path =r'D:/common/first_light/'
    suffix='v2d'
    Table(analysis_path, suffix, link_sed=True, verbose=False)
        
        
            
        
        
