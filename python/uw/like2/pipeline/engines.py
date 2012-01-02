"""
Support for running multiple IPEngines on a cluster or machines

$Header$
"""

import time, os, sys, types, subprocess
import numpy as np

from IPython import parallel
from IPython.parallel.util import interactive # decorator for function to run in an engine

version = '$Revision: 1.1 $'.split()[1]

class Engines(object):
    """ manage a set of IPengines
    """

    def __init__(self, default='default', **kwargs):
        """
        """
        self.logstream = kwargs.pop('log', None)
        self.quiet = kwargs.pop('quiet', False)
        self.sleep_interval = kwargs.pop('sleep_interval', 60.0)
        self.maxerrors = kwargs.pop('maxerrors', 5)
        self.usercallback = kwargs.pop('callback', None)
        self.elapsed=dict()
        self.results=dict()

        try:
            self.rc = parallel.Client(default=default)
        except:
            raise Exception( 'No engines available: you must run IPcluster')
        if not self.quiet: print '%d engines available' %len(self.rc)
        self.lview = self.rc.load_balanced_view()
        ncomp, nqueue, ntasks = self.status()
        if nqueue+ntasks>0: print 'warning: %d running tasks' % ( nqueue+ntasks)
   
    def __getitem__(self, i): return self.rc[i]
    def len(self): return len(self.rc)
    
    def clear(self, **kwargs):
        """
        Definition: self.rc.clear(self, targets=None, block=None)
        Docstring:  Clear the namespace in target(s).
        """
        self.rc.clear(**kwargs)
    
    def shutdown(self, **kwargs):
        """
        Definition: self.rc.shutdown(self, targets=None, restart=False, hub=False, block=None)
        Docstring:  Terminates one or more engine processes, optionally including the hub.
        """
        self.rc.shutdown(**kwargs)
        
    def execute(self, code):
        """ code: string to execute on all engines """
        self._wait_to_finish(self.rc[:].execute(code, block=False))

    def submit(self, code, *pars):
        """ submit jobs to the engines
        code : string to execute, should evaluate to a function that will be executed on the target machine
        pars : one or more list of arguments
        """
        self.log('submitting  %d tasks '%len(pars[0]))
        starttime=time.time()
        @interactive
        def fn(f,x):
            return eval(f,globals())(x)
        self._wait_to_finish(self.lview.map( fn, len(pars[0])*[code], *pars))
        self.log('Done: elapsed, total times=%4.0f, %.0f sec'\
                    %((time.time()-starttime), sum(self.elapsed.values())))

    def status(self):   
        stat = self.lview.queue_status()
        ncomp ,nqueue, ntasks  = [sum([q[key] for q in stat.values() if type(q)==types.DictType])\
                for key in ('completed', 'queue','tasks')]
        return ncomp, nqueue, ntasks
 
    def _wait_to_finish(self, amr):
        """ monitor progress for a set of taskss """
        rc = self.rc
        self.elapsed=dict()
        self.results=dict()

        errorcount=0
        nczero ,nqueue, ntasks  = self.status()

        # code adapted from custom written by MinRK <benjaminrk@gmail.com>
        pending = set(amr.msg_ids)
        while pending:
            try:
                rc.wait(pending, self.sleep_interval) #1e-3)
                loops = 0
            except parallel.TimeoutError:
                # ignore timeouterrors, since they only mean that at least one isn't done
                loops +=1
                if loops*self.sleep_interval >timeout:
                    raise parallel.TimeoutError
            except Exception, e:
                log('Exception %s'%e)
                errorcount +=1; 
                if errorcount > self.maxerrors:
                    log('Maximum error count: will leave the loop')
                    return
            # finished is the set of msg_ids that are complete
            finished = pending.difference(rc.outstanding)
            # update pending to exclude those that just finished
            pending = pending.difference(finished)
            for msg_id in finished:
                # we know these are done, so don't worry about blocking
                ar = rc.get_result(msg_id)
                self.results[msg_id] = ar.result[0] if ar.result is not None else None
                self.elapsed[msg_id] = (ar.completed-ar.started).total_seconds()
                if self.usercallback is not None:
                    # this needs some way to associate msg_id with the input parameters.
                    self.usercallback(msg_id, ar.stdout.replace('\n', '\n    ').rstrip())
            ncomp ,nqueue, ntasks  = self.status()
            self.log('completed: %4d, queued: %4d tasks: %4d  pending: %4d' %(ncomp-nczero,nqueue, ntasks, len(pending)))
    
    def log(self, message):
        if not self.quiet:
            print >>self.logstream,\
                '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (message,))
            if self.logstream is not None: self.logstream.flush()
            else: sys.stdout.flush()
