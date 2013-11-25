"""
 Manage associations
 $Header$
 author: T. Burnett <tburnett@uw.edu>
"""
import os, glob
import numpy as np
from skymaps import SkyDir
from uw.like import srcid  # gtsrcid work-alike from Eric


class SrcId(srcid.SourceAssociation):
    """
    adapter to simplify source association
    """
    def __init__(self, catalog_path='$FERMI/catalog/', classes='all_but_gammas', quiet=True):
        """ 
        catalog_path : string
            path to the catalogs, expect to find srcid/classes under it
        clases: ['all' | 'all_but_gammas' | list ]
            list of classes to apply,
        
        """
        self.classes = classes
        catalog_path = os.path.expandvars(catalog_path)
        d = os.path.join(catalog_path, 'srcid', 'classes')
        q = glob.glob(os.path.join(d, '*.py'))
        assert len(q)>0, 'no association classes found in folder %s' % d
        allclasses =[os.path.split(u)[-1].split('.')[0] for u in q if '__init__' not in u]
        if classes=='all': 
            # special tag to really get everything
            q = glob.glob(os.path.join(d, '*.py'))
            self.classes = allclasses
        elif self.classes=='all_but_gammas':
            self.classes = ['agn', 'bllac', 'bzcat', 'cgrabs', 'crates', 'crates_fom', 'dwarfs', 
            'galaxies', 'globular', 'hmxb', 'ibis', 'lbv', 'lmxb',  'ocl', 'ostar', 
             #'pulsar_fom',
            'pulsar_lat', 'pulsar_big', #'msp', 'pulsar_high',  'pulsar_low', 'pulsar_nonATNF', 
            'pwn', 'qso', 'seyfert', 'seyfert_rl', 'snr', 'snr_ext', 'starbursts', 'tev']
        else:
            self.classes=classes
        for c in self.classes:
            if c not in allclasses:
                txt = 'class %s not in set classes: %s' % (c, allclasses)
                raise Exception(txt)
        super(SrcId, self).__init__(os.path.join(catalog_path, 'srcid'),quiet=quiet)
        self.class_list = self.classes # will be used by the id method
     
    def __str__(self):
        return 'SrcId(%s)' %self.classes
    
    def __repr__(self):
        return '%s.%s: %s' % (self.__module__, self.__class__.__name__, self.classes)
        
    def __call__(self, name, pos, error):
        """ name: source name, ignored, for reference
            pos: a SkyDir object
            error: [float | tuple]
                a or tuple: (a,b,ang) or, (ra,dec,a,b,ang,)...
                a,b: 1-sigma elipse in deg; ang orientation degrees
            returns None, or a dictionary consistent with Association above. (elements sorted by prob.)
        """
        if not hasattr(error, '__iter__'):
            error = (error,error,0)
        else:
            if len(error)==7: 
                error = error[3:6] 
            assert len(error)==3, 'wrong length for error ellipse specification'
        if not isinstance(pos, SkyDir):
            pos = SkyDir(*pos)
        source_ass = self.id(pos,error)
        # select first association per catalog, rearrange to sort on over-all prob.
        items =  source_ass.items()
        candidates = [(v[0][1], v[0][0], v[0][2], key, v[0][3]) for key,v in items if v!={}]
        if len(candidates)==0: return None
        candidates.sort(reverse=True)
        # format as a dict to conform to Association above (TODO: a recarray would be cleaner)
        d = {
            'name':  [a[1] for a in candidates],
            'cat':   [a[3] for a in candidates],
            'dir':   [a[2] for a in candidates],
            'prob':  [a[0] for a in candidates],
            'ra':    [a[2].ra() for a in candidates],
            'dec':   [a[2].dec() for a in candidates],
            'ang':   [np.degrees(a[2].difference(pos)) for a in candidates],
            'deltats':[a[4] for a in candidates],
            }
        # stuff depending on the catalog 
        cats = [self.catalogs[a[3]] for a in candidates]
        d['prior'] =   [cat.prior for cat in cats]
        d['density'] = [cat.local_density(pos) for cat in cats]
        #d['fom'] =     [cat.fom for cat in cats]
        return d
        
    def positional_likelihood(self, name, pos, error):
        a,b,ang = error

