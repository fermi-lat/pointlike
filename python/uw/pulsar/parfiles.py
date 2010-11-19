
def sex2dec(s,mode='ra'):
    toks = s.split(':')
    multi = (-1 if '-' in s else 1) * (15 if mode == 'ra' else 1)
    return multi*sum( (60**(-i) * abs(float(toks[i])) for i in xrange(len(toks)) ) )

class ParFile(dict):

    def __init__(self,parfile):

        for line in file(parfile):
            tok = line.strip().split()
            if len(tok)==0: continue
            self[tok[0]] = tok[1:] if (len(tok[1:]) > 1) else tok[1:][0]
    
    def get(self,key):
        t = self[key]
        if type(t) is type([]):
            return t[0]
        return t

    def get_skydir(self):
        from skymaps import SkyDir
        ra = sex2dec(self.get('RAJ'))
        dec = sex2dec(self.get('DECJ'),mode='dec')
        return SkyDir(ra,dec)

    def p(self):
        f0 = self.get('F0')
        return 1./float(f0)
 
    def pdot(self):
        f1 = self.get('F1')
        return -float(f1)*self.p()**2

    def has_waves(self):
        return np.any(['WAVE' in key for key in self.keys()])

    def is_msp(self):
        return (self.p() < 0.03) and (self.pdot < 1e-18) and (not self.has_waves())

    def get_time_cuts(self):
        if self.is_msp():
            return 0, 1e10 # "infinity" in MET
        t0 = float(self.get('START'))
        t1 = float(self.get('FINISH'))
        mjd_0 = 51910 + 7.428703703703703e-4
        met_0 = (t0 - mjd_0) * 86400
        met_1 = (t1 - mjd_0) * 86400
        return max(met_0,0), max(met_1,0)