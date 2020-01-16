"""
Application module, allowing command-line access to analysis/plotting tasks

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/analyze/app.py,v 1.17 2016/10/28 20:48:13 burnett Exp $

"""
_locs = locals().keys()
from uw.like2.analyze import *

newlocs = locals().keys()
module_names = filter(lambda x:not x[0]=='_',set(newlocs).difference(_locs))
assert len(module_names)>0, 'No modules found'

#module_names = filter(lambda x:not x[0]=='_', uw.like2.analyze.__all__ ) #set(newlocs).difference(_locs))
#
from uw.like2.analyze import _html
import numpy as np
import sys, types

class AppMenu(dict):
    
    def __init__(self, names):
        assert len(names)>0, 'AppMenu failed: no module names available'
        super(AppMenu,self).__init__(zip(names, map(self._create_entry, names)))
        
    def _create_entry(self, name):
        pack = globals()[name]
        doc = pack.__doc__.split('\n')
        title = (doc[0] if doc[0]!='' else doc[1]).strip()
        try:
            classname = self._check_for_class(pack)
        except Exception as msg:
            print ('*** Failed to check module {}:{}'.format(name, msg))
            return None
        if classname is None:
            print ('***no class with all_plots found in module %s' %name)
            return None
        classobj = eval(name+'.'+classname)
        return dict(title=title, classname=classname, classobj=classobj, 
                    require=getattr(classobj,'require',None), 
                    module=pack,)

    def _check_for_class(self,pack):
        for name, value in pack.__dict__.items():
            if name[0]=='_': continue
            if 'all_plots' in value.__dict__: return name
        return None

    def __call__(self, name, args=None, reloadit=False):
        """ return an object for the class in this module, optionally reloading it"""
        try:
            pk = self[name]
        except:
            print ('KeyError: "%s" not in %s' % (name, sorted(self.keys())))
            raise
        if reloadit:
            print (reload(pk['module']))
            pk = self[name]= self._create_entry(name)
        return pk['classobj'](args=args)

        
    def __str__(self):
        s = '%-15s %s\n' % ('name', 'description')
        return s+ '\n'.join(['  %-15s %s' % (key, self[key]['title']) for key in sorted(self.keys())])  

menu = AppMenu(module_names)
        
        
def main(procs, args=None, update_top=False , raise_exception=False):
    np.seterr(invalid='warn', divide='warn')
    success=True
    if type(procs)==types.StringType: procs = procs.split()
    for arg in procs:
        # ## does not work, needs fixing
        #if arg=='all':
        #    cs = set(np.hstack(menu.values()))
        #    for cls in cs:
        #        if os.path.exists(cls.require):
        #            print ('running %s' % cls.__name__)
        #            try:
        #                cls('.').all_plots()
        #                plt.close('all')
        #            except Exception, msg:
        #                print ('=====failed====\n %s\n=============='% msg)
        #        else:
        #            print ('skipped %s, missing %s' % (cls.__name__, cls.require))
        #    break
        if arg=='menu': 
            update_top=True
            continue

        if arg not in menu.keys():
            print ('found %s; expect one of %s' % (arg, sorted(menu.keys())))
            success = False
            continue
        try:
            menu(arg,args=args).all_plots()
        except FloatingPointError as msg:
            print ('Floating point error running %s: "%s"' % (arg, msg))
            print ('seterr:', np.seterr())
            success=False
        except Exception as msg:
            print ('Exception running %s: "%s"' % (arg, msg))
            if raise_exception: raise
            success = False
    if success: 
        if update_top: _html.HTMLindex().update_top()
        
    return success  
      
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=""" Run one or more analysis applications for the sky model defined in the current folder. \nValid module names:\n """+menu.__str__()
    )
    parser.add_argument('module', nargs='+', help='module name: must be in the list of names in the table above')
    parser.add_argument('--args', default=None, help='argments for the processor')
    parser.add_argument('--update_top', action='store_true', help='Update the top level Web  menu')
    parser.add_argument('--raise_exception', action='store_true', help ='set to catch exceptions, default is to ignore.')
    args = parser.parse_args()
    if not main(args.module, args=args.args, update_top=args.update_top, raise_exception=args.raise_exception):
        sys.exit(1)
