"""
Support for running multiple IPEngines on a cluster or machines

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/engines.py,v 1.3 2012/06/24 13:50:00 burnett Exp $
"""

import time, os, sys, types, subprocess
import numpy as np

from IPython import parallel
from IPython.parallel.util import interactive # decorator for function to run in an engine

version = '$Revision: 1.3 $'.split()[1]

class Engines(object):
    """ manage a set of IPEngines to simplify parallel processing
    Expect that an IPController has been set up with a number of IPEngines
    
    Example:
    mec = Engines()
    mec.execute('f=lambda x: x')
    n=1000
    mec.submit('f', range(n) )
    assert sum(mec.results.values()) ==n*(n-1)/2, 'failed test'
    """

    def __init__(self, profile='default', **kwargs):
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
            self.rc = parallel.Client(profile=profile, **kwargs)
        except:
            raise Exception( 'No controller available: you must run ipcluster')
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
        self.log('executing setup code')
        self._wait_to_finish(self.rc[:].execute(code, block=False))

    def submit(self, f, *pars):
        """ submit jobs to the engines
        f : string  or function 
            if string, should by the name of a function that will be executed on the target machine
        pars : one or more list of arguments
        """
        assert len(pars)>0, 'no parameters'
        self.log('submitting  %d tasks to %d engines'% (len(pars[0]), self.len()))
        starttime=time.time()
        self.results = dict()
        
        # suggestion from MinRK 
        if isinstance(f, basestring):
            name = f
            f = parallel.Reference(name)
            f.__name__ = name # workaround: not needed in next ipython version
        self._wait_to_finish(self.lview.map( f, *pars))
        
        self.log('Done: elapsed, total times=%4.0f, %.0f sec'\
                    %((time.time()-starttime), sum(self.elapsed.values())))

    def status(self):   
        stat = self.lview.queue_status()
        ncomp ,nqueue, ntasks  = [sum([q[key] for q in stat.values() if type(q)==types.DictType])\
                for key in ('completed', 'queue','tasks')]
        return ncomp, nqueue, ntasks
 
    def _wait_to_finish(self, amr):
        """ monitor progress for a set of tasks """
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
                self.log('Exception other than Timeout: %s'%e)
                errorcount +=1; 
                if errorcount > self.maxerrors:
                    self.log('Maximum error count: will leave the loop')
                    return
            # finished is the set of msg_ids that are complete
            finished = pending.difference(rc.outstanding)
            # update pending to exclude those that just finished
            pending = pending.difference(finished)
            for msg_id in finished:
                # we know these are done, so don't worry about blocking
                ar = rc.get_result(msg_id)
                try:
                    self.results[msg_id] = ar.result[0] if ar.result is not None else None
                except Exception,errmsg:
                    self.log('Exception %s for id %s' % (errmsg, msg_id))
                    errorcount +=1; 
                    if errorcount > self.maxerrors:
                        self.log('Maximum error count: will leave the loop')
                        return
                    continue
                self.elapsed[msg_id] = (ar.completed-ar.started).total_seconds()
                if self.usercallback is not None:
                    # this needs some way to associate msg_id with the input parameters.
                    self.usercallback(msg_id, ar.stdout.replace('\n', '\n    ').rstrip())
            ncomp, nqueue, ntasks  = self.status()
            self.log('completed: %4d, queued: %4d tasks: %4d  pending: %4d' %(ncomp-nczero,nqueue, ntasks, len(pending)))
    
    def log(self, message):
        if not self.quiet:
            print >>self.logstream,\
                '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (message,))
            if self.logstream is not None: self.logstream.flush()
            else: sys.stdout.flush()

### TODO: utility functions to start/stop/kill a cluster
# ssh -n -f tev10 " sh -c \"nohup ipcluster start --n=24 > /dev/null 2>&1 &\""
# ssh -n -f tev09 " sh -c \"nohup ipcluster engines --n-24 > /dev/null 2>&1 &\""
# ...
#ssh -n -f user@host "sh -c \"cd /whereever; nohup ./whatever > /dev/null 2>&1 &\""