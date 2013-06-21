"""
Application module, allowing command-line access to analysis/plotting tasks

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/app.py,v 1.5 2013/06/20 23:53:47 burnett Exp $

"""
_locs = locals().keys()
from uw.like2.analyze import *
newlocs = locals().keys()
module_names = filter(lambda x:not x[0]=='_',set(newlocs).difference(_locs))
from uw.like2.analyze import html
import numpy as np
import sys, types

class AppMenu(dict):
    
    def __init__(self, names):
        super(AppMenu,self).__init__(zip(names, map(self._create_entry, names)))
        
    def _create_entry(self, name):
        pack = globals()[name]
        doc = pack.__doc__.split('\n')
        title = (doc[0] if doc[0]!='' else doc[1]).strip()
        classname = self._check_for_class(pack)
        if classname is None:
            print 'no class found in %s' %name
            return None
        classobj = eval(name+'.'+classname)
        return dict(title=title, classname=classname, classobj=classobj, require=getattr(classobj,'require',None),)

    def _check_for_class(self,pack):
        for name, value in pack.__dict__.items():
            if name[0]=='_': continue
            if 'all_plots' in value.__dict__: return name
        return None

    def __call__(self, name):
        return self[name]['classobj']()

    def __str__(self):
        s = '%-15s %s\n' % ('name', 'description')
        return s+ '\n'.join(['%-15s %s' % (key, self[key]['title']) for key in sorted(self.keys())])  

menu = AppMenu(module_names)
        
        
def main(args, update_top=False , raise_exception=False):
    np.seterr(invalid='warn', divide='warn')
    success=True
    if type(args)==types.StringType: args = args.split()
    for arg in args:
        # ## does not work, needs fixing
        #if arg=='all':
        #    cs = set(np.hstack(menu.values()))
        #    for cls in cs:
        #        if os.path.exists(cls.require):
        #            print 'running %s' % cls.__name__
        #            try:
        #                cls('.').all_plots()
        #                plt.close('all')
        #            except Exception, msg:
        #                print '=====failed====\n %s\n=============='% msg
        #        else:
        #            print 'skipped %s, missing %s' % (cls.__name__, cls.require)
        #    break
        if arg=='menu': 
            update_top=True
            continue

        if arg not in menu.keys():
            print 'found %s; expect one of %s' % (arg, menu.keys())
            success = False
            continue
        try:
            menu(arg).all_plots()
        except FloatingPointError, msg:
            print 'Floating point error running %s: "%s"' % (arg, msg)
            print 'seterr:', np.seterr()
            success=False
        except Exception, msg:
            print 'Exception running %s: "%s"' % (arg, msg)
            if raise_exception: raise
            success = False
    if success: 
        html.HTMLindex().create_menu()
        if update_top: html.HTMLindex().update_top()
        
    return success  
      
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=""" Run an analysis application \n """+menu.__str__()
    )
    parser.add_argument('args', nargs='+', help='processsor identifier: must be one of %s' %menu.keys())
    parser.add_argument('--update_top', action='store_true', help='Update the top level Web  menu')
    parser.add_argument('--raise_exception', action='store_true', help ='set to catch exceptions')
    args = parser.parse_args()
    if not main(args.args, update_top=args.update_top, raise_exception=args.raise_exception):
        sys.exit(1)
