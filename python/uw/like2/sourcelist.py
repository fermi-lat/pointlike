"""
Manage sources: single class SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sourcelist.py,v 1.4 2011/08/31 23:12:42 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""

import numpy as np

class SourceList(list):
    """ manage properties of the list of sources
        
    """

    def __init__(self, roi):
        """ initialization: get sources from the ROIAnalysis object: creates a polymorphic list of the 
            three types: point, global, and extended
            
            roi : a ROIAnalsis object:
                expect to find attributes:
                    dsm : a ROIDiffuseManager object, for the list of global and extended sources
                    psm : a ROIPointSourceManager object, with a list of point sources
                    roi_dir: the center of the ROI
        """
        class DiffuseSourceFactory(object):
            """ klugy way to create new object from diffuse source
            """
            def __init__(self, roi):   self.roi = roi
            def __call__(self, source_index, **kwargs):
                cm = self.roi.dsm.bgmodels[source_index]
                return cm.__class__(cm.sa, cm.diffuse_source, self.roi.roi_dir, **kwargs)
        class ExtendedSourceFactory(object):
            def __init__(self, roi):   self.roi = roi
            def __call__(self, source_index, **kwargs):
                cm = self.roi.dsm.bgmodels[source_index]
                return cm.__class__(cm.sa, cm.extended_source, self.roi.roi_dir, **kwargs)
        class PointSourceFactory(object):
            def __init__(self, roi):    self.roi =roi
            def __call__(self, source_index):
                return roi.psm.point_sources[source_index], self.roi.roi_dir
                
        # note that sources are added in the order diffuse, point to agree with ROIAnalysis
        # also, add two attributes to ecach source object and define 'spectral_model'
        for i,source in enumerate(roi.dsm.diffuse_sources):
            spatialclassname = roi.dsm.bgmodels[i].__class__.__name__
            source.factory = ExtendedSourceFactory(roi) if spatialclassname=='ROIExtendedModel'\
                else DiffuseSourceFactory(roi)
            source.manager_index = i
            source.spectral_model = source.smodel
            self.append(source)
        for i,source in enumerate(roi.psm.point_sources):
            source.manager_index = i
            source.factory = PointSourceFactory(roi)
            source.spectral_model = source.model
            self.append(source)
        self.source_names = np.array([s.name for s in self])
        self.models= np.array([s.model for s in self])
        self.free  = np.array([np.any(model.free) for model in self.models])
        
        ## temporary kluge to ensure that bands are setup for current use by diffuse
        ## this is needed if the ROI is generated without its setup, which adds these
        ## attributes to the bands
        nm  = len(roi.dsm.bgmodels)
        for band in roi.bands:
            band.bg_counts = np.empty(nm)
            band.bg_pix_counts = np.empty([len(band.wsdl),nm]) if band.has_pixels else 0

    def __str__(self):
        return 'SourceList, with %d sources, %d free parameters' %(len(self), sum(self.free))
        
    def find_source(self, source_name):
        names = [s.name for s in self]
        try:
            return self[names.index(source_name)]
        except:
            raise Exception('source %s not found' %source_name)

    @property
    def parameter_names(self):
        """ array of free parameter names """
        names = []
        for source_name, model in zip(self.source_names, self.models):
            for pname in np.array(model.param_names)[model.free]:
                names.append(source_name.strip()+'_'+pname)
        return np.array(names)
    
    def get_parameters(self):
        """ array of free parameters (fitter rep)"""
        if len(self.models)==0: return []
        return np.concatenate([m.get_parameters() for m in self.models])
    
    def set_parameters(self,parameters):
        """ set the (fitter rep) parameters"""
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

    def get_model_parameters(self):
        return 10**self.get_parameters()
    
    def set_model_parameters(self, par):
        self.set_parameters(np.log10(par))
    parameters = property(get_model_parameters, set_model_parameters,
                doc='array of free model parameters')

