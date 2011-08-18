"""
Manage sources

$Header$
Author: T.Burnett <tburnett@uw.edu>
"""

import numpy as np

class SourceList(list):
    """ manage properties of the list of sources
    """

    def __init__(self, roi):
        """ initialization: get sources from the ROIAnalysis object
        """
        class ConvolvedSourceFactory(object):
            """ klugy way to create new object from extended or diffuse source
                encapsulates how to get such from an ROI
            """
            def __init__(self, roi):   self.roi = roi
            def __call__(self, source_index, **kwargs):
                cm = self.roi.dsm.bgmodels[source_index]
                if hasattr(cm,'diffuse_source'):
                    return cm.__class__(cm.sa, cm.diffuse_source, self.roi.roi_dir, **kwargs)
                else:
                    return cm.__class__(cm.sa, cm.extended_source, self.roi.roi_dir, **kwargs)

        class PointSourceFactory(object):
            def __init__(self, roi):    self.roi =roi
            def __call__(self, source_index):
                return roi.psm.point_sources[source_index], self.roi.roi_dir
                
        # note that sources are added in the order diffuse, point to agree with ROIAnalysis
        # also, add two attributes to ecach source object 
        for i,source in enumerate(roi.bgm.diffuse_sources):
            source.factory = ConvolvedSourceFactory(roi)
            source.manager_index = i
            self.append(source)
        for i,source in enumerate(roi.psm.point_sources):
            source.manager_index = i
            source.factory = PointSourceFactory(roi)
            self.append(source)
        self.source_names = np.array([s.name for s in self])
        self.models= np.array([s.model for s in self])
        self.free  = np.array([np.any(model.free) for model in self.models])

    def __str__(self):
        return 'SoruceList, with %d sources, %d free parameters' %(len(self), sum(self.free))
    @property
    def parameter_names(self):
        """ array of free parameter names """
        names = []
        for source_name, model in zip(self.source_names, self.models):
            for pname in np.array(model.param_names)[model.free]:
                names.append(source_name.strip()+'_'+pname)
        return np.array(names)
    
    def get_parameters(self):
        """ array of free parameters (internal rep)"""
        if len(self.models)==0: return []
        return np.concatenate([m.get_parameters() for m in self.models])
    
    def set_parameters(self,parameters):
        """ set the (internal rep) parameters"""
        current_position=0
        for m in self.models:
            cp,np = current_position, current_position+ sum(m.free)
            m.set_parameters(parameters[cp:np])
            current_position += np-cp

    def set_covariance_matrix(self,cov_matrix,):
        """ take over-all covariance matrix from a fit, give block-diagonal pieces to model objects"""
        current_position=0
        for m in self.models:
            cp,np = current_position,current_position+len(m.get_parameters())
            m.set_cov_matrix(cov_matrix[cp:np,cp:np])
            current_position += np-cp

    
    def get_external(self):
        return 10**self.get_parameters()
    
    def set_external(self, par):
        self.set_parameters(np.log10(par))
    parameters = property(get_external, set_external, doc='array of free parameters')

