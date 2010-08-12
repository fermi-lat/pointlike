""" Class to write out gtlike-style results files. 

$Header:$

author: Joshua Lande
"""
from pprint import pformat
import numpy as N
from uw.utilities.xml_parsers import Model_to_XML
from uw.like.roi_extended import ExtendedSource
from skymaps import IsotropicSpectrum

def unparse_spectral(model,**kwargs):
    """ Convert a Model object to a gtlike style dictionary. """
    m2x=Model_to_XML()
    m2x.process_model(model,**kwargs)
    names,vals,errs=m2x.pname,m2x.pval,m2x.perr

    return dict((name,'%g +/- %g' % (p,perr) if perr>0 else '%g' % p)
                for name,p,perr in zip(names,vals,errs))

def unparse_spatial(model):
    """ Convert a SpatialModel object to a gtlike-inspired dictionary. """
    return dict((name,'%g +/- %g' % (p,perr) if perr != 0 else '%g' % p)
                for name,p,perr in zip(model.param_names,
                                       *model.statistical(absolute=True,two_sided=False)))

def unparse_point_sources(roi,point_sources,emin,emax,**kwargs):
    """ Convert a PointSource object into a gtlike style dictionary. """
    point_dict = {}
    for ps in point_sources:
        point_dict[ps.name]={'TS value':'%g' % roi.TS(which=ps,**kwargs)}
        point_dict[ps.name].update(unparse_spectral(ps.model,scaling=False))
        if not N.all(ps.model.cov_matrix==0):
            point_dict[ps.name]['Flux']='%g +/- %g' % ps.model.i_flux(emin,emax,cgs=True,two_sided=False,error=True)
        else:
            point_dict[ps.name]['Flux']='%g' % ps.model.i_flux(emin,emax,cgs=True,two_sided=False,error=False)
    return point_dict

def unparse_diffuse_sources(roi,diffuse_sources,emin,emax,**kwargs):
    diffuse_dict= {}
    for ds in diffuse_sources:
        dm = ds.dmodel
        if hasattr(dm,'__len__'):  dm = dm[0]

        diffuse_dict[ds.name]={}

        diffuse_dict[ds.name].update(
            unparse_spectral(ds.smodel,scaling=False if isinstance(ds,ExtendedSource) else True,
                             xml_name = 'FileFunction'\
                             if isinstance(dm,IsotropicSpectrum) else None))

        if isinstance(ds,ExtendedSource):
            diffuse_dict[ds.name].update(unparse_spatial(ds.spatial_model))

            if not N.all(ds.smodel.cov_matrix==0):
                diffuse_dict[ds.name]['Flux']='%g +/- %g' % ds.smodel.i_flux(emin,emax,cgs=True,two_sided=False,error=True)
            else:
                diffuse_dict[ds.name]['Flux']='%g' % ds.smodel.i_flux(emin,emax,cgs=True,two_sided=False,error=False)

            diffuse_dict[ds.name]['TS value']='%g' % roi.TS(which=ds,**kwargs)

    return diffuse_dict

def writeResults(roi,filename,**kwargs):
    emin,emax=roi.bin_edges[[0,-1]]
    if not roi.quiet:
        print "\nPhoton fluxes are computed for the energy range %d to %d" % (emin,emax)

    source_dict={}
    source_dict.update(unparse_point_sources(roi,roi.psm.point_sources,emin,emax,**kwargs))
    source_dict.update(unparse_diffuse_sources(roi,roi.dsm.diffuse_sources,emin,emax,**kwargs))

    file=open(filename,'w')
    file.write(pformat(source_dict))
    file.close()
