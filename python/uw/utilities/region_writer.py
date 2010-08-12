""" Class to write out region files compatable with ds9. 

$Header:$

author: Joshua Lande
"""
from uw.like.roi_extended import ExtendedSource
from uw.like.SpatialModels import RadiallySymmetricModel,SpatialMap,PseudoSpatialModel
from math import degrees

def unparse_point_sources(point_sources):

    return ["fk5; point(%.4f, %.4f) # point=cross text={%s}" % \
            (ps.skydir.ra(),ps.skydir.dec(),ps.name) \
            for ps in point_sources]

def unparse_diffuse_sources(diffuse_sources):

    lines = []
    for ds in diffuse_sources:
        if isinstance(ds,ExtendedSource):

            lines.append("fk5; point(%.4f, %.4f) # point=cross text={%s}" % \
                          (ds.spatial_model.center.ra(),
                           ds.spatial_model.center.dec(),
                           ds.name))

            if isinstance(ds.spatial_model,RadiallySymmetricModel) and not \
                (isinstance(ds.spatial_model,PseudoSpatialModel) or \
                 type(ds.spatial_model) == SpatialMap):
                    lines.append("fk5; circle(%.4f, %.4f, %.4f)" % \
                                  (ds.spatial_model.center.ra(),
                                  ds.spatial_model.center.dec(),
                                  degrees(ds.spatial_model.r68())))
    return lines

def writeRegion(roi,filename,color='green'):

    lines = [
        "# Region file format: DS9 version 4.0",
        "global color=%s" % color,
    ]

    lines += unparse_diffuse_sources(roi.dsm.diffuse_sources)
    lines += unparse_point_sources(roi.psm.point_sources)

    file=open(filename,'w')
    file.write('\n'.join(lines))
    file.close()
