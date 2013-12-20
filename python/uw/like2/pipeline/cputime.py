"""
Analyze CPU times for batch execution, make diagnostic plots

$Header$
"""
import os, glob, re,  pandas as pd, numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec


class CPUtime(object):
    def __init__(self, stream, path='.'):
        """Parse, for display, the cpu execution times found in the stremlogs
        
        stream : int
            stream number
        path : str
            location of the streamlogs file
        """
        self.stream=stream
        streamlogs=sorted(glob.glob(os.path.join(path,'streamlogs','stream%s.*.log'%stream)))
        if len(streamlogs)==0:
            print 'no logs found for stream %s' %stream
            streamlogs=glob.glob(os.path.join(path,'streamlogs','*.log'))
            streams = set([re.search(r'logs/stream(.*)\.\d', s).group(1) for s in streamlogs])
            print 'found streams %s' % list(streams)
            return    
            
        print 'found %d logs for stream %s' % (len(streamlogs), stream)
        self.model_version = '_'.join(os.path.realpath(path).split('/')[-2:])
        print 'model: %s' % self.model_version
        roi_info = dict()
        job_info = dict()
        def parse_logfile(fn):
            text = open(fn).read()
            lines = text.split('\n')
            hc = lines[0].split()[-1][:4]
            j = 0
            while lines[j].find('Start setup')>0 and lines[j+1].find('elapsed=')<0: j+=1
            setup = re.search(r'elapsed=(.*) \(',lines[j+1]).group(1)
            total = re.search(r'total (.*)\)', lines[-2]).group(1)
            job_info[fn]=dict(host=hc, setup=float(setup), total=float(total))
            n=0
            for j,line in enumerate(lines[:-1]):
                if line.find('Start roi')>0 and lines[j+1].find('Finish: elapsed=')>0:
                    next = lines[j+1]
                    roi = re.search(r'Start roi (.*)$', line).group(1)
                    elapsed = re.search(r'elapsed=(.*) \(',next).group(1)
                    roi_info[int(roi)] = dict(host=hc, time=float(elapsed), n=n)
                    n+=1
        for fn in streamlogs:
            try:
                parse_logfile(fn) 
            except Exception,msg:
                print 'Failed to parse file %s: %s' % (fn, msg)
        self.tdf, self.jdf = pd.DataFrame(roi_info).T, pd.DataFrame(job_info).T
        print 'Times: number  mean  minimum maximum'
        fmt = ' %-5s'+ '%6d'+'%8.1f'*3
        for label, ser in zip(('setup','ROI', 'Job'), (self.jdf.setup,self.tdf.time, self.jdf.total)):
            print fmt % (label, len(ser), ser.mean(), ser.min(), ser.max())
        
    def hists(self, tsmax=50, tmax=50, jobmax=2000, axx=None):
        tdf, jdf = self.tdf, self.jdf
        hosts = set(tdf.host)
        jcuts =[jdf.host==hc for hc in hosts]
        tcuts =[tdf.host==hc for hc in hosts]
        
        def setup_time(ax):
            ax.hist([jdf.setup[hc].clip(0,tsmax) for hc in jcuts],
                    np.linspace(0,tsmax,26), stacked=True, label=hosts)
            plt.setp(ax, xlabel='Job setup time')
            ax.legend(prop=dict(size=10)); ax.grid()
        def roi_time(ax):
            ax.hist([tdf.time[hc].clip(0,tmax) for hc in tcuts], np.linspace(10,tmax,26), stacked=True, label=hosts)
            plt.setp(ax, xlabel='ROI process time')
            ax.legend(prop=dict(size=10)); ax.grid()
        def total_time(ax):
            ax.hist( [jdf.total[hc].clip(0,jobmax) for hc in jcuts], 
                np.linspace(0,jobmax,26), stacked=True, label=hosts)
            plt.setp(ax, xlabel='total CPU time per job')
            ax.grid(); ax.legend(prop=dict(size=10))
        
        if axx is None:
            fig, axx = plt.subplots(1,3, figsize=(12,5))
        for fun, ax in zip((setup_time, roi_time, total_time), axx):
            fun(ax)
        return ax.figure
    
    def scats(self, tmax=50, ax=None):
        tdf = self.tdf
        if ax is None:
            fig, ax = plt.subplots(figsize=(12,5))
        for ht in set(tdf.host):
            cut = (tdf.host==ht)# & (tdf.n>0)
            ax.plot(tdf.index[cut], tdf.time[cut].clip(0,tmax), '.', label=ht)
        plt.setp(ax, ylim=(0,tmax+1), ylabel='time (s)', xlabel='ROI number')
        ax.legend(); ax.grid()
        ax.set_title('Stream %d'% self.stream, size=14)
        return ax.figure

    def all(self,tsmax=50, tmax=100, jobmax=2000,):
        """combine all plots in a 2x3 grid"""
        gs = gridspec.GridSpec(2,3)
        gs.update(wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(12,8))
        axx = [plt.subplot(gs[1,col]) for col in range(3)]
        axx.append(plt.subplot(gs[0,:]) )
        self.hists(tsmax, tmax, jobmax, axx=axx)
        self.scats(tmax, ax=axx[3])
        fig.text(0.05, 0.02, self.model_version, size=8)
        return fig

