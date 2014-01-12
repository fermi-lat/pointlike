"""
Output the ROI info to as a pickle file.
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/to_healpix.py,v 1.5 2013/12/22 01:17:25 burnett Exp $
"""
import os, pickle, time
import numpy as np

def pickle_dump(roi,  pickle_dir, dampen, **kwargs):
    """ dump the source information from an ROI constructed from the sources here
    """
    assert os.path.exists(pickle_dir), 'output folder not found: %s' %pickle_dir
    name = roi.name.strip()
    fname = kwargs.pop('fname', name)
    filename=os.path.join(pickle_dir,fname+'.pickle')
    
    # start output record, a dict
    output = dict()
    output['name'] = name
    output['skydir']  = roi.roi_dir
    output['parameters'] = roi.sources.parameters[:] 
    output['gradient'] = np.asarray(roi.gradient(), np.float32)

    # add a record to the history read in originally
    history = roi.sources.history if hasattr(roi.sources, 'history') else []
    history.append( 
        dict(stream=kwargs.get('stream', '-1'),
            time=time.asctime(),
            logl=roi.log_like(),
            dampen= dampen,
            )
        )
    output['history'] = history
    
    # add diffuse and extended info
    diffuse_sources =  [s for s in roi.sources if s.isglobal or s.isextended]
    
    output['diffuse'] = [s.spectral_model for s in diffuse_sources] #roi.dsm.models
    output['diffuse_names'] = [s.name for s in diffuse_sources]  #roi.dsm.names

    # a dict for all variable sources
    sources=dict()
    output['sources']= sources
    
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
            ellipse= s.__dict__.get('ellipse', None), 
            associations = s.__dict__.get('associations',None),
            )
    output.update(kwargs) # add additional entries from kwargs
    with open(filename,'wb') as f:  #perhaps overwrite
        pickle.dump(output,f)
    print 'saved pickle file to %s' % filename
        