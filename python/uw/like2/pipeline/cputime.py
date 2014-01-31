"""
Analyze CPU times for batch execution, make diagnostic plots

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/cputime.py,v 1.4 2013/12/31 04:49:30 burnett Exp $
"""
import os, glob, re, argparse, pandas as pd, numpy as np
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
        if not os.path.exists(os.path.join(path, 'streamlogs')):
            print 'No streamlogs folder in path %s' % path
            return
        streamlogs=sorted(glob.glob(os.path.join(path,'streamlogs','stream%s.*.log'%stream)))
        if len(streamlogs)==0:
            print 'no logs found for stream %s' %stream
            streamlogs=glob.glob(os.path.join(path,'streamlogs','*.log'))
            if len(streamlogs)==0:
                print 'No stream logs in %s' % path
                return
            streams =[]
            for s in streamlogs:
                t = re.search(r'logs/stream(.*)\.\d', s)
                if t is not None:
                    streams.append(t.group(1))
                    
            if len(streams)==0:
                print 'Did not find any streams'
            else:
                print 'found streams %s' % sorted(map(int,list(set(streams))))
            return    
            
        print 'found %d logs for stream %s' % (len(streamlogs), stream)
        self.model_version = '_'.join(os.path.realpath(path).split('/')[-2:])
        print 'model: %s' % self.model_version
        roi_info = dict()
        job_info = dict()
        
            
        def parse_logfile(fn):
 
            def search( pattern, text):
                x = re.search(pattern, text)
                if x is None:
                   print 'Processing %s:\n\t Failed to find "%s" in "%s"' % (fn, pattern, text)
                   return 0
                return x.group(1)

            text = open(fn).read()
            lines = text.split('\n')
            hc = lines[0].split()[-1][:-4]
            j = 0
            while lines[j].find('Start setup')>0 and lines[j+1].find('elapsed=')<0: j+=1
            setup = search(r'elapsed=(.*) \(',lines[j+1])
            total = search(r'total (.*)\)', lines[-2])
            job_info[fn]=dict(host=hc, setup=float(setup), total=float(total))
            n=0
            for j,line in enumerate(lines[:-1]):
                if line.find('Start roi')>0 and lines[j+1].find('Finish: elapsed=')>0:
                    next = lines[j+1]
                    roi = search(r'Start roi (.*)$', line)
                    elapsed = search(r'elapsed=(.*) \(',next)
                    roi_info[int(roi)] = dict(host=hc, time=float(elapsed), n=n)
                    n+=1
        for fn in streamlogs:
            try:
                parse_logfile(fn) 
            except Exception,msg:
                print 'Failed to parse file %s: %s' % (fn, msg)
        self.tdf, self.jdf = pd.DataFrame(roi_info).T, pd.DataFrame(job_info).T
        # correct for appent setup in first ROI per job
        self.tdf['tcorr'] = self.tdf.time
        self.tdf.tcorr[self.tdf.n==0] -= 15
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

    def plots(self,tsmax=50, tmax=100, jobmax=2000,):
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
        
    def make_list(self, cmax=250, outfile='joblist.txt'):
        i = 0
        cum= 999
        cums= []
        start=[]
        end = []
        while i<1728:
            t = self.tdf.tcorr[i]
            if cum+t < cmax:
                cum += t
            else:
                start.append(i)
                if cum!=999: 
                    cums.append(cum)
                    end.append(i)
                cum = 45 + t
            i+=1
        cums.append(cum)
        end.append(1728)
        if outfile is not None:
            with open(outfile,'w') as f:
                f.write('\n'.join(['%04d %04d'% t for t in zip(start,end)]))
        return start, cums


def main(stream=None, path='.', plotto=None):
    ct = CPUtime(stream, path)
    if plotto is not None:
        fig = ct.plots()
        if plotto.find('.')<0:
            plotto +='.png'
        plt.savefig(plotto)


if __name__=='__main__':
    parser = argparse.ArgumentParser(
            description="""Analyze a streaminfo folder
    """)
    parser.add_argument('stream', nargs='*', default=None, help='stream to process')
    parser.add_argument('--path', default=os.getcwd(), help='path to skymodel folder,  default: "%(default)s"')
    parser.add_argument('--plotto', default=None, help='name for plot, default "%(default)s"')
    args = parser.parse_args()
    main(int(args.stream[0]),args.path, args.plotto)
    

