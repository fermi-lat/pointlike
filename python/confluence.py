"""
Create a table of source fits as a confluence table

"""


class Table(object):
    def sinbad_ref(self, ra, dec):
        sd = '%2b' if dec>0 else '%2d'
        ad = abs(dec)
        coord = '%.4f%%20%s%.4f%%20'%(ra,sd, ad)
        
        return 'http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord=%s&submit=submit%%20query&Radius.unit=deg&CooEqui=2000&CooFrame=FK5&Radius=0.2'%coord
        #1.6362%20%2b73.0255%20
    def __init__(self, path, version, verbose=False):
        self.version=version
        self.verbose=verbose
        filename=r'%spointfit_%s.txt' %(path, version)
        self.outfilename='%sconfluence_%s.txt' %(path, version)
        infile =file(filename)
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
            q+= t+' || '
        q+= ' tentative association||'
        if self.verbose: print q
        self.outfile.write(q+'\n')
        
    def entry(self, line):
        toks = line.split()
        ra,dec = [float(t) for t in toks[1:3]]
        q = '|'
        for i,t in enumerate(toks):
            link = t.replace('+','_') # confluence attached file cannot have +
            if i==0: q+= '[%s|^%s_SED_%s.png] [S|%s] | ' % (t,link,self.version, self.sinbad_ref(ra,dec))
            else: q+= t + ' | '
        for i in range(self.cols-len(toks)+1): q+= ' |' # note extra
            
        if self.verbose: print q
        self.outfile.write(q+'\n')

if __name__=='__main__':
    analysis_path =r'D:/common/first_light/'

    suffix='06'
    Table(analysis_path, suffix, True)
        
        
            
        
        
