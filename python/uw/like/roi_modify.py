"""
Provides functions to modify an ROIAnalysis object.

$Header:

author: Matthew Kerr, Joshua Lande
"""
import numpy as np

from . Models import Model
from . SpatialModels import SpatialModel
from . roi_extended import ExtendedSource
from . pointspec_helpers import PointSource,get_default_diffuse_mapper
from . import roi_localize

from skymaps import SkyDir

def modify_loc(roi,skydir,which):
    """Move point source given by which to new location given by skydir."""
    manager,index=roi.mapper(which)
    if manager==roi.psm:
        rl = roi_localize.ROILocalizer(roi,which=index,update=True,bandfits=False)
        rl.spatialLikelihood(skydir,update=True)
    elif manager==roi.dsm:
        if isinstance(roi.dsm.diffuse_sources[index],ExtendedSource):
            roi.dsm.bgmodels[index].modify_loc(roi.bands,skydir)
        else:
            raise Exception("Unable to modify_loc of diffuse source %s" % which)

def modify_spatial_model(roi,which,spatial_model,keep_old_center=True):
    """ Modify a source's spatial model.

        If spatial_model is a SpatialModel object and which is a point
        source, the source will be converted into an extended source.'

        If spatial_model is a skydir object and which is an extended
        source, the source will be converted to a point source. If
        spatial_model is a skydir object and which is a PointSource,
        this function will simply call modify_loc. 
        
        keep_old_center => keep the center from the old spatial model.
        """
    manager,index = roi.mapper(which)
    source = roi.get_source(which)
    model = roi.get_model(which)

    if isinstance(spatial_model,SkyDir):
        if isinstance(source,PointSource):
            roi.modify_loc(which,spatial_model)
        else:
            roi.del_source(which)
            ps=PointSource(name=source.name,
                           model=model.copy(),skydir=spatial_model)
            roi.add_source(ps)

    elif isinstance(spatial_model,SpatialModel):

        if keep_old_center: 
            spatial_model.modify_loc(source.skydir)

        if manager==roi.psm:
            roi.del_source(which)
            es=ExtendedSource(name=source.name,
                              model=model.copy(),spatial_model=spatial_model)
            roi.add_source(es)
        else:
            if not isinstance(manager.diffuse_sources[index],ExtendedSource):
                raise Exception("Can only modify the spatial model of point and extended sources.")

            roi.dsm.diffuse_sources[index].spatial_model=spatial_model

            diffuse_mapper = get_default_diffuse_mapper(roi.sa,roi.roi_dir,roi.quiet)
            bgmodel = diffuse_mapper(roi.dsm.diffuse_sources[index])

            roi.dsm.bgmodels[index]=bgmodel
            roi.dsm.bgmodels[index].initialize_counts(roi.bands)
            roi.dsm.update_counts(roi.bands)
    else:
        raise Exception("spatial must be either a skydir or a SpatialModel")
    roi.__update_state__()

def modify_model(roi,which,model,free=None,keep_old_flux=True):
    """ Modify a spectral model for the source specified by which.
        
        The parameter free can be used to modify the free/frozen state
        of the source's spectral paremters. If both the free list and
        model are passed, the free list will override those in the
        spectral model.

        By default, the overal normalization of the spectral model
        will be ignored and the spectral model will be renormalized so
        that it contains the same flux (in the fitted energy range) as
        the previous model. This can be overridden with the keep_old_flux flag. """

    manager,index = roi.mapper(which)
    source = roi.get_source(which)

    if model is not None:
        if not isinstance(model,Model):
            raise Exception("model must be of type Models.Model")

        if keep_old_flux: 
            emin,emax=roi.bin_edges[[0,-1]]
            model.set_flux(roi.get_model(which).i_flux(emin=emin,emax=emax),emin=emin,emax=emax)

        manager.models[index]=model

        if hasattr(manager,'bgmodels') and hasattr(manager.bgmodels[index],'smodel'):
            manager.bgmodels[index].smodel=model

        if hasattr(source,'model'): source.model=model
        if hasattr(source,'smodel'): source.smodel=model

    if free is not None: 
        model=roi.get_model(which)

        if isinstance(free,bool):
            free=np.asarray([free]*len(model.get_all_parameters()))

        assert(len(free)==len(model.get_all_parameters()))
        for i in range(len(free)):
            model.freeze(i,freeze=not free[i])

    roi.__update_state__()

def modify_name(roi,which,name):
    """ Modify the name of a source in the ROI."""
    manager,index = roi.mapper(which)
    source = roi.get_source(which)

    manager.names[index]=name
    source.name=name

def modify_spectral_kwargs(roi,which,keep_old_flux,kwargs):
    """ Modify in the spectral model all of the parameters 
        specified by kwargs. """

    if len(kwargs) == 0: return

    source = roi.get_source(which)
    model = roi.get_model(which)

    if len([i for i in kwargs.keys() if i in model]) == 0:
        return

    for key in kwargs.keys():

        if key in model: model[key] = kwargs.pop(key)

    modify_model(roi,which,model,keep_old_flux=keep_old_flux)

def modify_spatial_kwargs(roi,which,keep_old_center,kwargs):
    """ Modify in the spatial model all of the parameters 
        specified by kwargs. """

    if len(kwargs) == 0: return

    source = roi.get_source(which)

    if hasattr(source,'spatial_model'):

        spatial_model = source.spatial_model
       
        if len([i for i in kwargs.keys() if i in spatial_model]) == 0:
            return

        for key in kwargs.keys():

            if key in spatial_model: 
                spatial_model[key] = kwargs.pop(key)

        modify_spatial_model(roi,which,spatial_model,keep_old_center)

def modify(roi,which, name=None, skydir=None,model=None,spatial_model=None,
        keep_old_flux=True,keep_old_center=True,free=None,**kwargs):
    """ This is a just a glue function wich will call all of the required
        modification functions to fully modify the source.

        See function modify_loc, modify_spatial_model, and modify_model
        to see what the various parameters do. """

    if skydir is None and model is None and \
            spatial_model is None and free is None and name is None and kwargs is {}:
        raise Exception("Parameter to modify must be specified.")

    if skydir is not None and spatial_model is not None:
        raise Exception("Both skydir and spatial_model cannot be specified.")

    if skydir is not None:
        modify_loc(roi,which=which,skydir=skydir)

    if spatial_model is not None:
        modify_spatial_model(roi,which,spatial_model,keep_old_center)

    if model is not None or free is not None:
        modify_model(roi,which,model,free,keep_old_flux)

    modify_spectral_kwargs(roi,which,keep_old_flux,kwargs)
    modify_spatial_kwargs(roi,which,keep_old_center,kwargs)

    if name is not None:
        modify_name(roi,which,name)

    if kwargs != {}: raise Exception("Unable to parse the kwargs=%s" % kwargs)

