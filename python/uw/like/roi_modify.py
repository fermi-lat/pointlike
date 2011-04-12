"""
Provides functions to modify an ROIAnalysis object.

$Header:

author: Matthew Kerr, Joshua Lande
"""
import numpy as N

from . Models import Model
from . SpatialModels import SpatialModel
from . roi_extended import ExtendedSource
from . pointspec_helpers import PointSource
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

def modify_spatial_model(roi,which,spatial_model,preserve_center=False):
    """ Modify a source's spatial model.

        If spatial_model is a SpatialModel object and which is a point
        source, the source will be converted into an extended source.'

        If spatial_model is a skydir object and which is an extended
        source, the source will be converted to a point source. If
        spatial_model is a skydir object and which is a PointSource,
        this function will simply call modify_loc. """
    manager,index = roi.mapper(which)
    source = roi.get_source(which)

    if isinstance(spatial_model,SkyDir):
        if isinstance(source,PointSource):
            roi.modify_loc(which,spatial_model)
        else:
            roi.del_source(which)
            ps=PointSource(name=source.name,model=source.model.copy(),skydir=spatial_model)
            roi.add_source(ps)

    elif isinstance(spatial_model,SpatialModel):

        if not preserve_center: spatial_model.modify_loc(source.skydir)

        if manager==roi.psm:
            roi.del_source(which)
            es=ExtendedSource(name=source.name,model=source.model.copy(),spatial_model=spatial_model)
            roi.add_source(es)
        else:
            if not isinstance(manager.diffuse_sources[index],ExtendedSource):
                raise Exception("Can only modify the spatial model of point and extended sources.")

            roi.dsm.diffuse_sources[index].spatial_model=spatial_model
            roi.dsm.bgmodels[index].initialize_counts(roi.bands)
            roi.dsm.update_counts(roi.bands)
    else:
        raise Exception("spatial must be either a skydir or a SpatialModel")
    roi.__update_state__()

def modify_model(roi,which,model,free=None,preserve_flux=False):
    """ Modify a spectral model for the source specified by which.
        
        The parameter free can be used to modify the free/frozen state
        of the source's spectral paremters. If both the free list and
        model are passed, the free list will override those in the
        spectral model.

        By default, the overal normalization of the spectral model
        will be ignored and the spectral model will be renormalized so
        that it contains the same flux (in the fitted energy range) as
        the previous model. This can be overridden with the preserve_flux flag. """

    manager,index = roi.mapper(which)
    source = roi.get_source(which)

    if model is not None:
        if not isinstance(model,Model):
            raise Exception("model must be of type Models.Model")

        if not preserve_flux: 
            emin,emax=roi.bin_edges[[0,-1]]
            model.set_flux(source.model.i_flux(emin=emin,emax=emax),emin=emin,emax=emax)

        manager.models[index]=model

        if manager==roi.psm:
            source.model=model
        else:
            source.smodel=model
            if isinstance(roi.dsm.diffuse_sources,ExtendedSource):
                source.model=model

    if free is not None: 
        model=roi.get_model(which)

        if isinstance(free,bool):
            free=N.asarray([free]*len(model.get_all_parameters()))

        assert(len(free)==len(model.get_all_parameters()))
        for i in xrange(len(free)):
            model.freeze(i,freeze=not free[i])
    roi.__update_state__()


def modify(roi,which=0,skydir=None,model=None,spatial_model=None,
        preserve_flux=False,preserve_center=False,free=None):
    """ This is a just a glue function wich will call all of the required
        modification functions to fully modify the source. """

    if skydir is None and model is None and spatial_model is None and free is None:
        raise Exception("Parameter to modify must be specified.")

    if skydir is not None and spatial_model is not None:
        raise Exception("Both skydir and spatial_model cannot be specified.")

    if skydir is not None:
        roi.modify_loc(which=which,skydir=skydir)

    if spatial_model is not None:
        roi.modify_spatial_model(which,spatial_model,preserve_center)

    if model is not None or free is not None:
        roi.modify_model(which,model,free,preserve_flux)

