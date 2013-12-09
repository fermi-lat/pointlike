"""
Output the ROI info to as a pickle file.
$Header$
"""
import os, pickle, time
import numpy as np

def pickle_dump(roi,  pickle_dir, dampen=1, **kwargs):
    """ dump the source information from an ROI constructed from the sources here
    """
    assert os.path.exists(pickle_dir), 'output folder not found: %s' %pickle_dir
    name = roi.name.strip()
    fname = kwargs.pop('fname', name)
    filename=os.path.join(pickle_dir,fname+'.pickle')
    if os.path.exists(filename):
        # if the file exists
        output=pickle.load(open(filename))
        # update history of previous runs
        prev_logl = output.get('prev_logl', [])
        if prev_logl is None or len(prev_logl)==0 or prev_logl[0] is None: prev_logl = []
        if dampen>0 : # update list only if a new fit
            last_logl = output.get('logl')
            if last_logl is not None: prev_logl.append(last_logl)
            output['prev_logl'] = prev_logl
        oldsrc = output.get('sources', dict())
        print 'updating pickle file: log likelihood history:', \
             ''.join(map(lambda x: '%.1f, '%x, prev_logl)) 
    else:
        output = dict()
        oldsrc = dict()
    diffuse_sources =  [s for s in roi.sources if s.isglobal or s.isextended]
    
    output['name'] = name
    output['skydir']  = roi.roi_dir
    # add extended sources to diffuse for backwards compatibilty: extended are also in the sources dict.
    output['diffuse'] = [s.spectral_model for s in diffuse_sources] #roi.dsm.models
    output['diffuse_names'] = [s.name for s in diffuse_sources]  #roi.dsm.names
    output['parameters'] = roi.sources.parameters[:] 
    output['gradient'] = np.asarray(roi.gradient(), np.float32)
    output['time'] = time.asctime()        
    sources=dict()
    output['sources']= sources
    def getit(s, key, savekey=None):
        """ key: current key
            savekey: key to save, expect to find in saved pickle
        """
        if savekey is None: savekey=key
        t = s.__dict__.get(key, None)
        if t is not None: return t #use a new version
        if s.name in oldsrc: #was it measured before?
            t = oldsrc[s.name].get(savekey, None)
            if t is not None: #yes, make a note and return it
                print 'ROI pickle: keeping previous calculation of %s for %s' % (key, s.name)
        return t
    
    for s in roi.free_sources:
        if s.isglobal: continue # skip globals
        try:
            pivot_energy = s.spectral_model.pivot_energy()
        except: # if not fit
            pivot_energy = None 
        if pivot_energy is None: pivot_energy=s.spectral_model.e0 # could be a pulsar?
        
        if 'sedrec' not in s.__dict__ or s.sedrec is None:
            print 'warning: no sedrec in %s' %s.name
            s.sedrec = None
        sedrec = s.sedrec
        loc = s.__dict__.get('loc', None)
        sources[s.name]=dict(
            skydir=s.skydir, 
            model=s.spectral_model,
            isextended=s.isextended,
            ts = s.ts if hasattr(s, 'ts') else roi.TS(s.name),
            sedrec = sedrec,
            band_ts=0 if sedrec is None else sedrec.ts.sum(),
            pivot_energy = pivot_energy,
            # if ellipse or adict not done, but already in pickle, keep them
            ellipse= getit(s, 'ellipse'), #s.__dict__.get('ellipse', None), 
            associations = getit(s, 'adict', 'associations'), #s.__dict__.get('adict',None),
            )
    output.update(kwargs) # add additional entries from kwargs
    with open(filename,'wb') as f:  #perhaps overwrite
        pickle.dump(output,f)
    print 'saved pickle file to %s' % filename
        
