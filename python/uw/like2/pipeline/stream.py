"""
manage creating new PipelineII stream

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/stream.py,v 1.6 2018/01/27 15:38:26 burnett Exp $
"""
import os, datetime, glob, argparse
import numpy as np
import pandas as pd


class PipelineStream(object):
    """ manage starting streams
    Assume that POINTLIKE_DIR is defined, and that current directory is a skymodel
    """
    def __init__(self, fn='summary_log.txt'):
        self.pointlike_dir=os.path.expandvars('$POINTLIKE_DIR')
        self.summary= os.path.join(self.pointlike_dir,fn)
        self.fullskymodel = os.getcwd()
        assert os.path.exists(self.summary), 'File "%s" not found' % self.summary
        assert os.path.exists('config.txt') or os.path.exists('../config.txt'), \
            'File config.txt not found in %s or its parent' %self.fullskymodel
        t =self.fullskymodel.split('/')
        self.model = '/'.join(t[t.index('skymodels')+1:])
        with open(self.summary, 'r') as slog:
            lines = slog.read().split('\n')
            last = lines[-1] if lines[-1] != '' else lines[-2] 
            try:
                self.stream_number = int(last.split()[0])
            except:
                print 'failed to interpret stream number from file %s: last=%s' % (self.summary,self.last)
                raise
            lastsp = last.split()
            self.last_job=dict(stream=self.stream_number, skymodel=lastsp[3], 
                stage=lastsp[4], model=lastsp[3], job_list=lastsp[5])
            
    def __call__(self,  stage, job_list, test=True):
        """Submit a command to the Pipeline to create a stream
        """
        pipeline='/afs/slac/g/glast/ground/bin/pipeline -m PROD createStream --stream %d -D "%s" UWpipeline '

        self.stream_number += 1
        time = str(datetime.datetime.today())[:16]
        line = '%5d  %s %-15s %-10s %s' % (self.stream_number, time, self.model, stage, job_list)
        print line
        if not test:
            with open(self.summary, 'a') as slog:
                slog.write('\n'+line)
        os.environ['SKYMODEL_SUBDIR']=self.fullskymodel
        cmd=pipeline % (self.stream_number, 
                "stage=%s, POINTLIKE_DIR=%s, SKYMODEL_SUBDIR=%s, job_list=%s"
                    % (stage, self.pointlike_dir, self.fullskymodel, os.path.expandvars(job_list))
                )
        if not test:
            print '-->' , cmd
            os.system(cmd)
        else:
            print 'Test mode: would have submitted \n\t%s'%cmd
            self.stream_number -=1

    def restart(self, jobs=None, test=True):
        """Restart the presumably hung jobs, by starting a new stream, the same stage as the current one,
        but with only the specified (hung) substreamss
        Do this by creating a file with the substream info, with a file name the same as the hung stream
        Its presence should cause any jobs in the current stream to abort when starting.

        jobs: list of int
            the substream ids. 
        """
        if jobs is None:
            jobs = SubStreamStats().notdone
            print 'Found {} jobs to restart'.format(len(jobs))
        if len(jobs)==0:
            print 'no jobs to restart'
            return
        assert self.model==self.last_job['skymodel'], 'Must run in same folder as last job'
        # get the current job_list
        substreams=t = list(substreamids(self.last_job['job_list']))+[1728]
        assert np.all([job in substreams for job in jobs]), 'job problem'
        last_stream = str(self.last_job['stream'])
        with open(last_stream, 'w') as out:
            for job in jobs:
                i = t.index(job)
                j,k = t[i:i+2]
                out.write('{:5d}{:5d}\n'.format( j,k) )
        if test: 
            print 'Testing: Created file "{}" to rerun substreams: {}'.format(last_stream, jobs)
            print 'Now starting test of restart'
        self(self.last_job['stage'], '$SKYMODEL_SUBDIR/'+last_stream, test)                
            

def substreamids(job_list):
    # check the substream ids
    jl = os.path.expandvars(job_list)
    if not os.path.exists(jl):
        jl = os.path.expandvars('$POINTLIKE_DIR/infrastructure/'+ job_list)
    assert os.path.exists(jl), 'Did not find file {}'.format(jl)
    lines = open(jl).read().split('\n')
    sslist= []
    assert len(lines)>1, 'Bad file {}'.format(jl)
    for line in lines:
        if line[0]=='#':continue
        sslist.append(line.split()[0])
    return np.array(sslist,int)


class StreamInfo(dict):
    """A dictionary of the stream info
    """
    def __init__(self, model=None, fn='summary_log.txt'):
        """
        model : string, default None
            If set, a filter. 
            the name of a model, e.g. "P301_6years/uw972", or a group with ending wildcard
        
        """
        self.pointlike_dir=os.path.expandvars('$POINTLIKE_DIR')
        self.summary= os.path.join(self.pointlike_dir,fn)
        t = open(self.summary,'r').read().split('\n')
        for line in t[2:]:
            if len(line)==0: 
                continue
            tokens = line.split()
            if len(tokens)<6: 
                print 'bad line, %s: %d tokens' % (line, len(tokens))
                continue
            if model is None or model==tokens[3] or model[-1]=='*' and tokens[3].startswith(model[:-1]):
                self[int(tokens[0])] = dict(stage=tokens[4], date=tokens[1]+' '+tokens[2], model=tokens[3],
                job_list=tokens[5],)
    def __call__(self, stream ):
        return self[stream]

def recent_stream(model_name=None, filter=None):
    """ return a dict, key the model name of the most recent stream, with stage and date

    model_name : str | None
        the full name of the model, e.g. P302_8years/uw8000. If None, use the current folder path
    """
    if model_name is None: model_name='/'.join(os.getcwd().split('/')[-2:])
    sinfo = StreamInfo(model_name)
    sdf = pd.DataFrame(sinfo).T
    # select last one for each model
    recent = dict()
    for model,s in zip(sdf.model, sdf.index):
        m = model.split('/')[-1]
        if filter is not None and not filter(m): continue
        date = sdf.ix[s].date
        stage = sdf.ix[s].stage
        job_list= sdf.ix[s].job_list
        if m not in recent: recent[m]=dict(stream=s, date=date, stage=stage, job_list=job_list)
        else: recent[m].update(stream=s, date=date,stage=stage, job_list=job_list)
    return recent

class SubStreamStats(object):

    def __init__(self, stream_id=None, path='.'):
        """Print info for substreams of the specified stream, or the most recent one
        Class has members with text, and a time dataframe for further analysis
        """
        if stream_id is None:
            sinfo = recent_stream()
            model,v = sinfo.items()[0]
            stream_id=v['stream']
        else:
            v = StreamInfo()[int(stream_id)]
        v['stream']=stream_id
        print v
        self.info = v
        search_string= os.path.join(path,'streamlogs','stream%s.*.log'%stream_id)
        streamlogs=sorted(glob.glob(search_string))
        if len(streamlogs)==0:
            print 'cwd', os.getcwd()
            print 'Did not find stream {} logfiles with search {}'.format( stream_id, search_string)
            return
        
        # check the substream ids
        jl = os.path.expandvars(v['job_list'])
        if not os.path.exists(jl):
            jl = os.path.expandvars('$POINTLIKE_DIR/infrastructure/'+v['job_list'])
        assert os.path.exists(jl), 'Did not find file {}'.format(jl)
        lines = open(jl).read().split('\n')
        sslist= []
        assert len(lines)>1, 'Bad file {}'.format(jl)
        for line in filter(lambda line: len(line)>0 and line[0]!='#', lines):
            sslist.append(line.split()[0])
        self.ssidlist = np.array(sslist,int)

        etimes = []
        substream = []
        nex = []
        rois=[]
        self.fulltext=dict()
        notdone_dict=dict()
        def parse_it(ft):
            lines= ft[0].split()
            a,b,host = lines[5], lines[7],lines[-1]
            t = filter(lambda line: line.find('Start roi')>0, ft)
            if len(t)==0:#file exists, but no roi's started
                current=-1
            else:
                try:
                    current = t[-1].split()[-1]
                except:
                    print 'Fail parse for current: "{}"'.format(t)
                    current = -2
            return dict(host=host, first=int(a), last=int(b), current=int(current))

        for f in streamlogs:
            text=open(f).read()
            
            lines = text.split('\n')
            if lines[-2].find('Finish')<0:
                ss = f[-13:-4]
                #print 'Stream {} not finished:\n{}'.format(ss, text)
                notdone_dict[ss] = parse_it(lines)
                continue
            try:
                etimes.append(float(text.split('\n')[-2].split()[-1][:-1]))
            except:
                print 'fail to parse: stream {}, file {}: \n{}'.format(stream_id,f,text)
                continue #raise
            nex.append( sum(np.array(['setup' in line for line in text.split('\n') ] )))
            id = int(f.split('.')[-2]) 
            substream.append(id )
            rois.append(len(set([line[line.find('roi')+4:]
                     for line in lines if line.find('Start roi')>0 ])))

            self.fulltext[id]=lines
        found = set(sorted(self.fulltext.keys()))
        expected = set(self.ssidlist)
        self.notdone =  [v['first'] for v in notdone_dict.itervalues()]
        self.notstarted = sorted(filter(lambda x: int(x) not in self.notdone, list(expected.difference(found)))) 
        if len(self.notstarted)>0:
            missing=''
            t=list(self.ssidlist)+[1728]
            for j in sorted(self.notstarted):
                i = t.index(j)
                j,k = t[i:i+2] 
                missing += '{}..{}, '.format(j,k-1) if k>j+1 else '{}, '.format(j)
            print '**** missing ROIs: {}'.format(missing)
   
        if len(self.notdone)>0:
            print '************* Not finished: *************'
            print pd.DataFrame(notdone_dict, index='host first last current'.split()).T
            print '*****************************************\n'
        print 'Found %d streamlog files for stream %s ...' %(len(streamlogs), stream_id), 
        self.times =times = pd.DataFrame(dict(etimes=etimes,nex=nex, rois=rois), 
            index=pd.Index(substream, name='substream'))
     
        etimes = times.etimes
        print '\texecution times: mean=%.1f, max=%.1f'% (etimes.mean(), etimes.max())
        print '\ttotal execution time: %.1f hours' % (etimes.sum()/3600.)
        print '\tnumber of executions: %s' % np.histogram(times.nex, range(7))[0][1:]
        self.tail_check()

    def tail_check(self, lim=2.0):
        times = self.times
        bad = times[times.etimes>lim*times.etimes.mean()]
        if len(bad)>0:
            print '\tsubstreams with times>{}*mean: \n{}'.format(lim,bad)


def main( stream=None ):

    ss=SubStreamStats(stream)

if __name__=='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description=""" Check stats for current or designated stream
    """)
    parser.add_argument('stream', nargs='*', default=None, help='optional Strem number')
    args = parser.parse_args()
    stream =(args.stream[0] if len(args.stream)>0 else None)
    os.environ['POINTLIKE_DIR']='/afs/slac/g/glast/groups/catalog/pointlike'
    os.environ['SKYMODEL_SUBDIR']=os.getcwd()
    main(stream)
    