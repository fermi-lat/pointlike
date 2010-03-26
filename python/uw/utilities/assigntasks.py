"""
  Assign a set of tasks to multiengine clients

  $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/assigntasks.py,v 1.10 2010/02/24 21:20:01 wallacee Exp $

"""
from IPython.kernel import client
import time, os, pickle
version = '$Revision: 1.10 $'.split()[1]


def get_mec():
    return client.MultiEngineClient()

class AssignTasks(object):
    """
        Schedule a set of tasks using a set of IPython's MultiEngineClient "engines"
        """

    def __init__(self, setup_string, tasks, post=None,
            mec=None, local=False, 
            quiet=False, log=None, timelimit=60,
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
        """
        self.local = local
        self.timelimit = timelimit
        self.setup_string = setup_string
        self.post = post
        if not local:
            try:
                self.mec = mec if mec is not None else get_mec()
            except:
                print 'No connection available: you must run ipcluster'
                raise
            assert(len(self.mec.get_ids())>0)
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

        self.log('Start AssignTasks, with %d tasks and  %d engines'\
               %( len(self.tasks), len(self.get_ids())))

    def callback(self, result, id):
        # still experimental--would like to to be notified
        print 'call back, id=%d, result=%s' % (id,result)
        index = self.assigned[id]
        self.log('%4d --> %s' %(index, result))
        self.result[index]=result
        self.time[index]=time.clock()-self.time[index]
        self.assigned[id]=None
        self.assign(id)

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
        self.log('Assigning task# %d to Engine %d'% ( self.index, id) )
        self.assigned[id]=self.index
        self.time[self.index]=self.engine_time[id]=time.clock()
        self.pending[id] = self.execute(self.tasks[self.index], id, block=False)
        self.index+=1
        if self.index==len(self.tasks): # check if done
            self.index=0 # reset index for another run
            return True  
        return False

    def log(self, message):
        if not self.quiet:
            print >>self.logstream,\
                '%4d-%02d-%02d %02d:%02d:%02d - %s' %(time.localtime()[:6]+ (message,))

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
            self.engine_time[id]=time.clock()
            self.time[-1] = time.clock()

        index = self.assigned[id]
        if index is None: return True # not assigned, ready

        if id>=0:
            try:
                result = self.pending[id].get_result(block=block)
                if result is None:
                    # waiting: check time limit
                    if time.clock()-self.engine_time[id]>self.timelimit:
                        self.log('Engine %d exceeded time limit (%f s) --lost task %d'\
                            %(id, self.timelimit, index))
                        self.lost.add(index)
                        self.mec.kill(targets=id, block=False)
                        # todo: if this does not work, need to remove from list, since get_ids will keep getting it
                    return False
                #if not self.quiet: print >>self.log, '%4d --> %s' %(index, result)
                self.result[index]=result
            except:
                self.log("Engine %d raised exception executing task %d" %(id, index))
                self.lost.add(index)
                #raise #TODO: option to allow dropping the engine and continue? Probably a serious issue.
        else:
            self.result[index]=self.tasks[index]

        self.time[index]=time.clock()-self.time[index]
        self.assigned[id]=None
        return True

    def __call__(self, sleep_interval=1):
        """
        process all tasks
        save results in dict self.result, execution times in dict self.time
        sleep_interval should be on the order of the time it takes a task to complete
        """

        # loop through the tasks, assigning as engines are available
        loop_iters = sleepcount=0
        starttime= time.clock()

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
        still = self.get_ids()
        curtime = time.clock()
        while True:
            if time.clock()-curtime>self.timelimit:
                self.log('quitting, exceeded time limit: %f' %self.timelimit)
                break
            time.sleep(sleep_interval)
            still = [id for id in still if not self.check_result(id)]
            if len(still)==0: break
            self.log('Still running: %d engines %d ... %d'% (len(still), still[0], still[-1]))
        if self.post is not None: 
            self.execute(self.post, self.get_ids(), True)

        if self.lost:
            self.log('possibly failed tasks: %s' % list(self.lost).sort() )
        self.log( 'Cycled through main loop %d times, slept %d times; elapsed, total times: %.1f, %.1f s'\
                %( loop_iters, sleepcount, time.clock()-starttime, sum(self.time.values())-self.time[-1]) )


    def dump(self, filename='summary.pickle'):
        """ save info to a pickle file """
        pickle.dump({'time':self.time, 'result':self.result}, open(filename, 'w'))


def setup_mec(engines=None, machines='tev1 tev2 tev3 tev4'.split()):
    """ On windows:start cluster and 4 engines on the local machine, in the current directory
        On linux: (our tev cluster) start controller on local machine, 16 engines/machine in all
    """
    if os.name=='nt':
        engines = engines or 4
        os.system(r'start /MIN /D %s cmd /K python C:\python25\scripts\ipcluster local -xy -n %d'% (os.getcwd(),engines))
    else:
        # on a tev machine
        engines = engines or 16 #default on a tev machine!
        
        #clusterfile_data='send_furl = False'\
        #    + '\nengines={'\
        #    + '\n'.join(['\t"%s" : %d,'%(m,engines) for m in machines])\
        #    + '\n}'
        #print 'cluster info:\n%s' % clusterfile_data
        #ofile=open('clusterfile', 'w')
        #ofile.writelines(clusterfile_data)
        #ofile.close()
        #os.system('ipcluster ssh -xy --clusterfile clusterfile &')

        # old, klugy way
        os.system('ipcontroller local -xy&')  
        for m in machines:
            for i in range(engines):
                time.sleep(0.1) # make sure the controller is started ?
                os.system('ssh %s ipengine&'% m) # this assumes that the environment is setup with non-interactive login

def kill_mec():
    get_mec().kill(True)


def test(tasks=9, **kwargs):
    at = AssignTasks('import time,random', 
            ['t=round(10.*random.random()); print %d,t; time.sleep(t)'%i for i in range(tasks)], 
            **kwargs)
    at()
    return at

if __name__=='__main__':
    at=test()
