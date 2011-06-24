"""
  Assign a set of tasks to multiengine clients

  $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/assigntasks.py,v 1.19 2011/03/30 17:36:58 wallacee Exp $

"""
from IPython.kernel import client
import time, os, pickle, subprocess
version = '$Revision: 1.19 $'.split()[1]

class ProgressBar:
    def __init__(self, total=60, width=40):
        self.total = total
        self.prog_bar = '[]'
        self.fill_char = '#'
        self.width = width
        self.__update_amount(0)
    
    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) / 2) - len(str(percent_done))
        pct_string = '%3d%%' % min(percent_done,99)
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])
        
    def update(self, current):
        if current>=self.total:
            self.prog_bar = '['+(self.width-1)*self.fill_char+']'
        else:
            self.__update_amount((current / float(self.total)) * 100.0)
        self.prog_bar += '  %d/%s' % (current, self.total)
        
    def __str__(self):
        return str(self.prog_bar)
    def display(self):
        print self()+'\r'

def get_mec():
    return client.MultiEngineClient()

class AssignTasks(object):
    """
        Schedule a set of tasks using a set of IPython's MultiEngineClient "engines"
        """

    def __init__(self, setup_string, tasks, post=None,
            mec=None, local=False, 
            quiet=False, log=None, timelimit=60,
            callback = None,
            ignore_exception = True,
            progress_bar = False,
            ):
        """
        setup_string: python code to setup the clients for the tasks
        tasks: a list of strings of python code for the individual tasks
        post: string to execute on each engine when done (e.g., delete an object)
        mec: [None] if specified, a  MultiEngineClient object to use
        local: [False] if True, run a single stream in the local processor
                very useful for debugging
        quiet: [False] control output
        log:   [None]  if set, an open output stream for log info
        timelimit: [60] limit (s) for a single task
        callback: [None] if set, will call as callback(task_index, result) 
        ignore_exception: [True] Determine what to do if an engine has an exception: continue, or throw it
        progress_bar : [False] If set, change quiet to false, log to open('assign_tasks.log', 'w')
        """
        self.local = local
        self.timelimit = timelimit
        self.setup_string = setup_string
        self.post = post
        self.usercallback = callback
        if not local:
            try:
                self.mec = mec if mec is not None else get_mec()
            except:
                print 'No connection available: you must run ipcluster'
                raise
            time.sleep(0.5) 
            assert len(self.mec.get_ids())>0, 'no engines found'
        self.assigned = {}
        self.tasks = tasks
        assert(len(tasks)>0)
        self.quiet = quiet
        self.logstream = log
        self.index = 0
        self.pending  = {}
        self.result   = {}
        self.time     = {}  #dictionary of time per task (id=-1 is startup)
        self.engine_time={} #dict. of time per engine
        self.lost     = set() # list of dead ids
        self.ignore_exception = ignore_exception
        self.progress_bar = progress_bar
        if progress_bar:
            if self.logstream is None:
                self.logstream = open('assigntasks.log', 'w')
            self.progress_bar=ProgressBar(len(tasks))

        self.log('Start AssignTasks, with %d tasks and  %d engines'\
               %( len(self.tasks), len(self.get_ids())))

    def get_ids(self):
        """ return a list of available engine ids (but no more than the number of tasks) 
            if local mode, the id number is -1
        """
        if self.local: return [-1]
        ids = set(self.mec.get_ids())-set(self.lost) # note ignoring "lost" ids
        ntasks = len(self.tasks)
        return list(ids) if ntasks>len(ids) else list(ids)[:ntasks]

    def execute(self, string, id, block=True):
        """ execute the string on engine id  
        """
        if id<0:
            exec(string, globals())
            return string
        else:
            return self.mec.execute(string, id, block)

    def assign(self, id=0):
        """ assign next task to engine id
            return done status
        """
        self.log('Assigning task# %d (%s) to Engine %d'% ( self.index, self.tasks[self.index], id) )
        self.assigned[id]=self.index
        self.time[self.index]=self.engine_time[id]=time.time()
        self.pending[id] = self.execute(self.tasks[self.index], id, block=False)
        self.index+=1
        if self.index==len(self.tasks): # check if done
            #self.index=0 # reset index for another run
            return True  
        return False

    def log(self, message):
        if not self.quiet:
            print >>self.logstream,\
                '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (message,))
            if self.logstream is not None: self.logstream.flush()
        if self.progress_bar:
            self.progress_bar.update(self.index+1)
            print self.progress_bar.__str__()+ ' '+message+ ('\r' if os.name=='nt' else chr(27)+'[A'),

    def check_result(self, id, block=False):
        """ 
            check result from engine id
            if not initialized, run the setup code on it and return pending result
            if ready, save result in self.result and release the engine
            if not ready, just return False (except if block is True, will always wait
            return ready status in any case
        """
        if id not in self.get_ids():
            self.log("Engine %d dropped out: lost task %d" %(id, self.assigned[id]))
            return True

        if id not in self.assigned:
            self.log( "Engine %d added to pool" % (id))
            self.pending[id]=self.execute(self.setup_string, id, block=False)
            self.assigned[id]=-1 #flag that running setup
            self.engine_time[id]=time.time()
            self.time[-1] = time.time()

        index = self.assigned[id]
        if index is None: return True # not assigned, ready

        if id>=0:
            try:
                result = self.pending[id].get_result(block=block)
                if result is None:
                    # waiting: check time limit
                    if time.time()-self.engine_time[id]>self.timelimit:
                        self.log('Engine %d exceeded time limit (%f s) --lost task %d'\
                            %(id, self.timelimit, index))
                        self.lost.add(index)
                        self.mec.kill(targets=id, block=False)
                        # todo: if this does not work, need to remove from list, since get_ids will keep getting it
                    return False
                #if not self.quiet: print >>self.log, '%4d --> %s' %(index, result)
                if self.usercallback is not None and index>=0: self.usercallback(index, result)
                self.result[index]=result
            except:
                self.log("Engine %d raised exception executing task %d ('%s')" %(id, index, 'setup' if index<0 else self.tasks[index]) )
                if not self.ignore_exception:                
                    raise
                self.lost.add(index)
                #raise #TODO: option to allow dropping the engine and continue? Probably a serious issue.
        else:
            self.result[index]=self.tasks[index]

        self.time[index]=time.time()-self.time[index]
        self.assigned[id]=None
        return True

    def __call__(self, sleep_interval=1):
        """
        process all tasks
        save results in dict self.result, execution times in dict self.time
        sleep_interval should be on the order of the time it takes a task to complete
        """

        # loop through the tasks, assigning as engines are available
        loop_iters = wait_iters =sleepcount=0
        starttime= time.time()

        done=busy = False

        while not done:
            if busy and sleep_interval>0: #sleep if nothing ready
                sleepcount +=1
                time.sleep(sleep_interval)

            for id in self.get_ids(): 
                # loop thru engines
                if self.check_result(id):
                    done = self.assign(id) # start next task, flag if ran out 
                    if done: break 
                    busy=True  # busy set if any task is started: will force a sleep
            loop_iters += 1

        # all tasks have been assigned: wait for all engines to finish
        still = [id for id in self.get_ids() if not self.check_result(id)]
        curtime = time.time()
        while len(still)>0:
            self.log('Still running: %d engines %d ...%d tasks %d...%d'%\
                    (len(still), still[0], still[-1], self.assigned[still[0]],self.assigned[still[-1]]))
            if time.time()-curtime>self.timelimit:
                self.log('quitting, exceeded time limit: %f' %self.timelimit)
                break
            if wait_iters> 1000:
                self.log('quitting, in apparent loop')
                for t in still: self.lost.add(self.assigned[t])
                break
            time.sleep(sleep_interval)
            still = [id for id in still if not self.check_result(id)]
            wait_iters+=1
            
        if self.post is not None: 
            self.execute(self.post, self.get_ids(), True)

        if self.lost:
            self.log('possibly failed tasks: %s' % list(self.lost).sort() )
        self.log( 'Cycled through main loop %d times, slept %d times; elapsed, total times: %.1f, %.1f s'\
                %( loop_iters, sleepcount, time.time()-starttime, sum(self.time.values())-self.time[-1]) )
        if self.logstream is not None: self.logstream.close()


    def dump(self, filename='summary.pickle'):
        """ save times, logs to a pickle file """
        pickle.dump({'time':self.time, 'result':self.result}, open(filename, 'wr'))


def setup_mec(engines=None, machines=None, clusterfile='clusterfile.txt',clobber=False):
    """ On windows:start cluster and 4 engines on the local machine, in the current directory
        On linux: (our tev cluster) start controller on local machine, 16 engines/machine in all
            If clusterfile.txt exists in the path: use it
    """
    if os.name=='nt':
        engines = engines or 4
        os.system(r'start /MIN /D %s cmd /K ipcluster local -xy -n %d'% (os.getcwd(),engines))
    else:
        # on a tev machine
        cf = clusterfile if os.path.exists(clusterfile) else '../'+clusterfile
        assert os.path.exists(cf), 'clusterfile %s not found in local or containing folder' % clusterfile
        if not clobber:
            os.system('ipcluster ssh -xy --clusterfile %s &' %cf)
            return
        else:
            os.remove(cf)
        #engines = engines or 16 #default on a tev machine!
        if engines is None:
            engines = 16
        if type(engines) == type([]):
            clusterfile_data = '\n'.join(['send_furl=False',
                                         'engines={',
                                         '\n'.join(['\t"%s":%d,'%(machines[i],engines[i]) for i in xrange(len(machines))]),
                                         '}'])
        else:
            clusterfile_data = '\n'.join(['send_furl = False',
                                          'engines={',
                                          '\n'.join(['\t"%s" : %d,'%(m,engines) for m in machines]),
                                          '}'])

        print 'cluster info:\n%s' % clusterfile_data
        ofile=open(clusterfile, 'w')
        ofile.writelines(clusterfile_data)
        ofile.close()
        os.system('ipcluster ssh -xy --clusterfile %s &' %clusterfile)

def kill_mec():
    get_mec().kill(True)

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
        print 'callback: %d, %s' % (id, result)
    at = AssignTasks('import time,random', 
            ['t=round(10.*random.random()); print %d,t; time.sleep(t)'%i for i in range(tasks)], 
            #callback = callback,
            **kwargs)
    at()
    return at

if __name__=='__main__':
    at=test()
