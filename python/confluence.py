"""
Create a table of source fits as a confluence table

"""
class RegFile(list):
    """ manage reg file showing initial and final positions, the latter with an error circle"""

    def __init__(self, data):
        for s in data: self.append(s) # must be better way
        print 'found %d sources' % len(self)
        
    def write_regfile(self, filename='test.reg', color='red'):
        header='global color=%s font="helvetica 10 normal" select=1 edit=1 move=0 delete=1 include=1 fixed=0 width=2;fk5;'%color
        out = file(filename, 'w')
        print >> out, header
        for s in self:
            print >> out,\
                  'point(%5.3f, %5.3f) # point=cross 20 text={%s}; circle(%5.3f,%5.3f,%5.3f) #  text={TS=%5.0f};'\
                      % ( s.ra, s.dec, s.name, s.ra_fit,s.dec_fit, (2.5*s.sigma), s.TS)
        out.close()

class FitSources(list):
    
    class Source(object):
        def __init__(self, name, ra, dec,  TS, ra_fit, dec_fit, sigma):
            self.name, self.ra, self.dec, self.TS, self.ra_fit, self.dec_fit, self.sigma=\
            ( name,    ra,      dec,      TS,      ra_fit,      dec_fit,     sigma)
            #print name,    ra,      dec,      TS,      ra_fit,      dec_fit,     sigma
            
    def __init__(self, t, select):
        """ table was created by Table; select is a selection array, like (moved<4) """
        from numpy import sum
        name =   t.table.field(0)
        ra_fit = t.table.field(1)
        dec_fit= t.table.field(2)
        TS =     t.table.field(3)
        sigma =  t.table.field(4)
        inname=t.intable.field(0)
        ra =   t.intable.field(1)
        dec =  t.intable.field(2)
        check = sum(name!=inname)
        if  check>0:
            print name; print inname
            raise Exception('the lists are not consistent, differ in %d places' %check)
        for i,s in enumerate(select):
            if s:
               self.append(FitSources.Source(name[i], ra[i], dec[i], TS[i], ra_fit[i], dec_fit[i], sigma[i]))

    

class Table(object):
    def sinbad_ref(self, ra, dec):
        sd = '%2b' if dec>0 else '%2d'
        ad = abs(dec)
        coord = '%.4f%%20%s%.4f%%20'%(ra,sd, ad)
        
        return 'http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord=%s&submit=submit%%20query&Radius.unit=deg&CooEqui=2000&CooFrame=FK5&Radius=0.2'%coord
        #1.6362%20%2b73.0255%20
    def link_ref(self, ra):
        return '^RA%03d.png' % int(ra)
    
    def __init__(self, sourcefile, resultfile, verbose=False, link_sed=True):
        import matplotlib
        self.verbose=verbose
        self.link_sed = link_sed
        self.table  =matplotlib.mlab.csv2rec(resultfile, delimiter=' ')
        self.intable=matplotlib.mlab.csv2rec(sourcefile, delimiter=' ')

        self.TS = self.table.field(3)
        self.moved = self.table.field(4)
        self.good = (self.TS>10)* (self.moved<4)


        #self.write_table
        #self.histogram()

    def write_table(self):        
        self.outfilename='%sconfluence_%s.txt' %(path, version)
        infile =file(self.infilename)
        self.outfile = file(self.outfilename, 'w')
        titles = infile.readline()[2:] # split off initial #?
        self.cols = len(titles.split())
        self.title(titles)
        for line in infile:
            self.entry(line)
        self.outfile.close()

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
        from pylab import figure, hist, grid, axis, plot, axvline, title, xlabel,ylabel, text, savefig
        from numpy import arange,array,exp
        figure(figsize=(4,4))
        ratio=[s[5] for s in self.table if s[4]<0.25] # ratio if localization < 0.25
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
        axvline(2.45); text(2.5, 0.75*axis()[3], '95% error circle')
        outfile = self.infilename.replace('.txt', '.png')
        savefig(outfile)
        print 'wrote png histogram to %s' % outfile


if __name__=='__main__':
#    analysis_path =r'D:/common/first_light/'
#    suffix='v2d'
#    Table(analysis_path, suffix, link_sed=True, verbose=False)
    t = Table('cgrabs_sorted.txt', 'pointfit_02.txt')     
        
            
        
        
