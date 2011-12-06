"""
  Assign a set of tasks to multiengine clients

  $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/assigntasks.py,v 1.20 2011/06/24 04:34:19 burnett Exp $

"""
from IPython import parallel
import time, os, sys, types, pickle, subprocess
import numpy as np

version = '$Revision: 1.20 $'.split()[1]

#def execute(string, id):
#    try:
#        exec(string,  globals(), locals())
#        return (id,'ok')
#    except Exception, e: return (id,e)


class AssignTasks(object):
    """
        Schedule a set of tasks using a set of IPython's MultiEngineClient "engines"
        """

    def __init__(self,fn, taskids, 
            local=False, 
            quiet=False, log=None, timelimit=60,
            callback = None,
            ignore_exception = True,
            maxerrors = 2,
            except_count=5
            ):
        """
        parameters
        ----------
        fn : function to run
        taskids: list of integers
        local: bool
            [False] if True, run a single stream in the local processor
                very useful for debugging
        timelimit: int
            [60] limit (s) for a single task
        callback: None or function
            if set, will call as callback(task_index, result) 
        ignore_exception: bool
            [True] Determine what to do if an engine has an exception: continue, or throw it
        """
        self.fn = fn
        self.local = local
        self.timelimit = timelimit
        self.except_count=except_count
        self.usercallback = callback
        self.quiet = quiet
        if not local:
            try:
                self.mec = parallel.Client()
                dv = self.mec[:]
                dv.scatter('id', self.mec.ids, flatten=True)
                self.lview = self.mec.load_balanced_view()
                self.ids = self.mec.ids
                self.lview.block = False
            except:
                print 'No connection available: you must run ipcluster'
                raise
        else:
            self.ids = [-1]
        self.taskids = taskids
        assert len(taskids)>0, 'no tasks!'
        self.logstream = log
        self.index = 0
        self.ignore_exception = ignore_exception
        self.maxerrors = maxerrors
        self.starttime=time.time()
        self.elapsed=dict()
        for task in self.taskids:
            self.elapsed[task] = 0.
        
        self.log('Start AssignTasks, with %d tasks and  %d engines' % ( len(self.taskids), len(self.ids) ) )

    def _wait_to_finish(self, sleep_interval, timeout=100):
        amr,rc = self.async, self.mec
        # code adapted from custom written by MinRK <benjaminrk@gmail.com>
        self.pending = set(amr.msg_ids)
        first = True; nczero = errorcount=0
        while self.pending:
            try:
                rc.wait(self.pending, sleep_interval) #1e-3)
                loops = 0
            except parallel.TimeoutError:
                # ignore timeouterrors, since they only mean that at least one isn't done
                loops +=1
                if loops*sleep_interval >timeout:
                    raise parallel.TimeoutError
            except Exception, e:
                log('Exception %s'%e)
                errorcount +=1; 
                if errorcount > self.maxerrors:
                    log('Maximum error count: will leave the loop')
                    raise
            # finished is the set of msg_ids that are complete
            finished = self.pending.difference(rc.outstanding)
            # update pending to exclude those that just finished
            self.pending = self.pending.difference(finished)
            for msg_id in finished:
                # we know these are done, so don't worry about blocking
                ar = rc.get_result(msg_id)
                taskid,status = ar.result[0] # the output from execute
                #print 'taskid, time started, completed', taskid, ar.started, ar.completed
                self.elapsed[taskid] = (ar.completed-ar.started).total_seconds()
                if self.usercallback is not None:
                    self.usercallback(taskid, ar.stdout.replace('\n', '\n    ').rstrip())
                #if status == 'ok':
                #    pass #self.log('engine %2d finished task %4d (%s)'%(ar.engine_id, taskid, self.tasks[taskid]))
                #else:
                #    self.log('engine %2d  raised exception task %4d (%s): %s'\
                #            %(ar.engine_id, taskid, status))
                #    if not self.ignore_exception:
                #        raise RuntimeError(status)
            stat = self.lview.queue_status()
            ncomp ,nqueue, ntasks  = [sum([q[key] for q in stat.values() if type(q)==types.DictType])\
                for key in ('completed', 'queue','tasks')]
            if first: 
                nczero = ncomp
                first = False
            self.log('completed: %4d, queued: %4d tasks: %4d  pending: %4d' %(ncomp-nczero,nqueue, ntasks, len(self.pending)))

              
    def __call__(self, sleep_interval=5):
        """
        process all tasks
        """
        running = True
        if not self.local:
            self.lview.wait()
            self.log('submitting  %d tasks '%len(self.taskids))
            self.async = self.lview.map( self.fn, self.taskids )
            while self.except_count>0:
                try:
                    self._wait_to_finish(sleep_interval)
                    break
                except Exception, e:
                    self.log('Exception raised: %s'%e)
                    self.except_count-=1
                
                #sself.log('failed to commplete, %d unfinished' % len(self.pending))
                
        else:
            for taskid in self.taskids:
                self.fn(taskid)
            if self.post is not None: execute(self.post, -2)
        self.log('Done: elapsed, total times=%4.0f, %.0f sec' %((time.time()-self.starttime), sum(self.elapsed.values())))
    
    def log(self, message):
        if not self.quiet:
            print >>self.logstream,\
                '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (message,))
            if self.logstream is not None: self.logstream.flush()
            else: sys.stdout.flush()
 
    
    
def get_mec():
    return parallel.Client()

def kill_mec():
    try:
        get_mec().shutdown(True)
    except:
        print 'No mecs to kill'

def get_machines(machines=None, clusterfile='clusterfile.txt'):
    if machines is None:
        exec(open(clusterfile if os.path.exists(clusterfile) else '../'+clusterfile))
        machines = sorted(engines.keys())
    return machines

def free(machines=None, clusterfile='clusterfile.txt', cmd='free'):
    """ run the command on the machines specified, or get from the clusterfile 
    cmd : string
        default is 'free', which has special formating also suggest 'pkill -9 ipengine' to clean up after disaster.
    """
    if cmd=='free':
        print '\t             total       used       free     shared    buffers     cached'
    else: print 'running command "%s"' % cmd
    for m in get_machines(machines, clusterfile):
        print m, ':\t'
        t=subprocess.Popen('ssh %s %s'% (m,cmd) , shell=True, stdout=subprocess.PIPE).communicate()
        if cmd=='free':
            for line in  t[0].split('\n')[1:]: print '\t'+line
        else: print t
        
    
def test(tasks=9, **kwargs):
    def callback(id,result): 
        print 'callback from task %d: %s' % (id, result)
    atsk = AssignTasks('import time,random',
            ['t=round(10.*random.random()); print %d,t; time.sleep(t)'%i for i in range(tasks)], 
            callback = callback,
            **kwargs)
    t = atsk()
    return t, atsk

if __name__=='__main__':
    at=test()
