""" Class to write out gtlike-style results files. 

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/results_writer.py,v 1.10 2012/07/09 19:16:53 lande Exp $

author: Joshua Lande
"""
from pprint import pformat
import numpy as N
from uw.like.roi_extended import ExtendedSource
from skymaps import IsotropicSpectrum

def unparse_spectral(model,**kwargs):
    """ Convert a Model object to a gtlike style dictionary. """
    return {n:'%g +/- %g' % (model[n],model.error(n)) \
            if model.has_errors() else '%g' % model[n] \
            for n in model.param_names}

def unparse_spatial(model):
    """ Convert a SpatialModel object to a gtlike-inspired dictionary. """
    return dict(('%s' % name,'%g +/- %g' % (p,perr) if perr != 0 else '%g' % p)
                for name,p,perr in zip(model.param_names,
                                       *model.statistical(absolute=True,two_sided=False)))

def unparse_point_sources(roi,point_sources,emin,emax,**kwargs):
    """ Convert a PointSource object into a gtlike style dictionary. """
    point_dict = {}
    for ps in point_sources:

        name = str(ps.name)
        point_dict[name]={'TS value':'%g' % roi.TS(which=ps,**kwargs)}
        point_dict[name].update(unparse_spectral(ps.model,scaling=False))
        if ps.model.has_errors():
            point_dict[name]['Flux']='%g +/- %g' % ps.model.i_flux(emin,emax,cgs=True,two_sided=False,error=True)
        else:
            point_dict[name]['Flux']='%g' % ps.model.i_flux(emin,emax,cgs=True,two_sided=False,error=False)
    return point_dict

def unparse_diffuse_sources(roi,diffuse_sources,emin,emax,**kwargs):
    diffuse_dict= {}
    for ds in diffuse_sources:
        dm = ds.dmodel
        if hasattr(dm,'__len__'):  dm = dm[0]

        name = str(ds.name)

        diffuse_dict[name]={}

        diffuse_dict[name].update(unparse_spectral(ds.smodel))

        if isinstance(ds,ExtendedSource):
            diffuse_dict[name].update(unparse_spatial(ds.spatial_model))

            if ds.smodel.has_errors():
                diffuse_dict[name]['Flux']='%g +/- %g' % ds.smodel.i_flux(emin,emax,cgs=True,two_sided=False,error=True)
            else:
                diffuse_dict[name]['Flux']='%g' % ds.smodel.i_flux(emin,emax,cgs=True,two_sided=False,error=False)

            diffuse_dict[name]['TS value']='%g' % roi.TS(which=ds,**kwargs)

    return diffuse_dict

def writeResults(roi,filename=None,**kwargs):
    """ Saves out an ROI to a gtlike style results file. """
    emin,emax=roi.bin_edges[[0,-1]]
    if not roi.quiet:
        if filename is not None: print ("\nSaving ROI to results file %s" % filename)

        print ("\nPhoton fluxes are computed for the energy range %d to %d" % (emin,emax))

    source_dict={}
    source_dict.update(unparse_point_sources(roi,roi.psm.point_sources,emin,emax,**kwargs))
    source_dict.update(unparse_diffuse_sources(roi,roi.dsm.diffuse_sources,emin,emax,**kwargs))

    if filename is not None:
        file=open(filename,'w')
        file.write(pformat(source_dict))
        file.close()

        if not roi.quiet:
            print ("\nDone Saving ROI to results file %s" % filename)
    
    return source_dict
