"""
manage creating new PipelineII stream

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/pipeline/stream.py,v 1.4 2015/08/16 01:12:11 burnett Exp $
"""
import os, datetime
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
            
    def __call__(self,  stage, job_list, test=True):
        """Submit a command to the Pipeline to create a stream
        """
        pipeline='/afs/slac/g/glast/ground/bin/pipeline -m PROD createStream --stream %d -D "%s" UWpipeline '

        self.stream_number += 1
        time = str(datetime.datetime.today())[:16]
        line = '%5d  %s %-15s %-10s %s' % (self.stream_number, time, self.model, stage, job_list)
        print line
        with open(self.summary, 'a') as slog:
            slog.write('\n'+line)
        cmd=pipeline % (self.stream_number, 
                "stage=%s, POINTLIKE_DIR=%s, SKYMODEL_SUBDIR=%s, job_list=%s"
                    % (stage, self.pointlike_dir, self.fullskymodel, job_list)
                )
        if not test:
            print '-->' , cmd
            os.system(cmd)
        else:
            print 'Test mode: would have submitted \n\t%s'%cmd

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
                self[int(tokens[0])] = dict(stage=tokens[4], date=tokens[1]+' '+tokens[2], model=tokens[3])
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
        if m not in recent: recent[m]=dict(stream=s, date=date, stage=stage)
        else: recent[m].update(stream=s, date=date,stage=stage)
    return recent