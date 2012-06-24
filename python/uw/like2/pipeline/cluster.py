"""
Support for managing a cluster

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/cluster.py,v 1.1 2012/02/26 23:43:49 burnett Exp $
"""

import time, os, sys, types, subprocess
import numpy as np

from IPython import parallel
version = '$Revision: 1.1 $'.split()[1]

# setup
master='tev10'
slaves=['tev09', 'tev08']
profile='default'
npe=24 # number per engine to start
print 'master: %s, slaves: %s' % (master, slaves)
def clear():
    # clear 
    cmd = 'ssh %s ipcluster stop'%master
    print cmd, ':', os.system(cmd)
    #for machine in [master]+slaves: 
    #    cmd='ssh %s killall --user $USER python' % machine
    #    print cmd, ':', os.system(cmd)
        
def start_master(ipcmd='start', n=npe):
    # start
    # ssh -n -f user@host "sh -c \"cd /whereever; nohup ./whatever > /dev/null 2>&1 &\""
    cmdt="""ssh -n -f %(host)s "sh -c \\" nohup %(what)s > /dev/null 2>&1 &\\"" """
    cmd =cmdt % dict(host=master, where='.', what='ipcluster %s --n=%d'% (ipcmd, n))
    print cmd
    return subprocess.Popen(cmd , shell=True, stdout=subprocess.PIPE).communicate()
    
def engines( slave, n=npe):
    cmdt="""ssh -n -f %(host)s "sh -c \\" nohup %(what)s > /dev/null 2>&1 &\\"" """
    cmd =cmdt % dict(host=slave, where='.', what='ipcluster %s --n=%d'% ('engines', n))
    print cmd
    return subprocess.Popen(cmd , shell=True, stdout=subprocess.PIPE).communicate()

    
def system(host, cmd):
    return subprocess.Popen('ssh %s %s'% (host,cmd) , shell=True, stdout=subprocess.PIPE).communicate()

def start_engines(n, delay=10):
    for slave in slaves:
        if n<=0: return
        time.sleep(delay)
        engines(slave, n=min(n,npe))
        n-=npe

def check(profile='default'):
    try: 
        rc = parallel.Client(profile=profile)
        print '%d engines' %len(rc)
    except:
        print 'no cluster!'
        
def free():
    for machine in [master]+slaves:
        print machine
        for line in system(machine, 'free')[0].split('\n'): 
            print '\t%s' % line
            
def startup(n=None, interval=5, retry=5):
    #clear()
    if n is None:
        n = npe*(len(slaves)+1)
    start_master(n=min(n,npe))
    start_engines(n-npe)
    last=-1
    for i in range(retry):
        time.sleep( interval if i>0 else 10)
        try: 
            rc = parallel.Client(profile=profile)
            m =len(rc)
        except Exception, msg:
            m=0
            print msg
        print  m,
        if m==n: 
            print 'done' 
            return n
        last=m
    print 'timed out, missing %d engine(s)'%(n-m)
    return m