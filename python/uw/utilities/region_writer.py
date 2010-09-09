""" Class to write out region files compatable with ds9. 

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/utilities/region_writer.py,v 1.1 2010/08/12 23:09:22 lande Exp $

author: Joshua Lande
"""
from uw.like.roi_extended import ExtendedSource
from uw.like.SpatialModels import *
from math import degrees

def unparse_point_sources(point_sources):

    return ["fk5; point(%.4f, %.4f) # point=cross text={%s}" % \
            (ps.skydir.ra(),ps.skydir.dec(),ps.name) \
            for ps in point_sources]

def unparse_diffuse_sources(diffuse_sources):

    lines = []
    for ds in diffuse_sources:
        if isinstance(ds,ExtendedSource):
            sm = ds.spatial_model

            ra,dec=sm.center.ra(),sm.center.dec()

            lines.append("fk5; point(%.4f, %.4f) # point=cross text={%s}" % \
                          (ra,dec,ds.name))

            if isinstance(sm,SpatialModel):

                if isinstance(sm,PseudoSpatialModel) or type(sm) == SpatialMap:
                    continue

                if isinstance(sm,RadiallySymmetricModel):
                    sigma=sm.sigma
                    if isinstance(sm,Disk):
                        lines.append("fk5; circle(%.4f, %.4f, %.4f) # Circle encloses all of the disk." % \
                                      (ra,dec,degrees(sigma)))
                    elif isinstance(sm,Ring):
                        frac=sm.frac
                        lines += ["fk5; circle(%.4f, %.4f, %.4f)" % \
                                      (ra,dec,degrees(_)) for _ in [frac*sigma,sigma]]
                    else:    
                        lines.append("fk5; circle(%.4f, %.4f, %.4f) # Circle contaning 68 percent of the source." % \
                                      (ra,dec,degrees(sm.r68())))

                elif isinstance(sm,EllipticalSpatialModel):
                    sigma_x, sigma_y, theta = sm.sigma_x, sm.sigma_y, sm.theta
                    if isinstance(sm,EllipticalDisk):
                        lines.append("fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                                (ra,dec,degrees(sigma_x),degrees(sigma_y),degrees(sm.theta)+90))

                    elif isinstance(sm,EllipticalRing):
                        frac = sm.frac
                        lines += ["fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                                (ra,dec,degrees(_*sigma_x),degrees(_*sigma_y),degrees(sm.theta)+90) \
                                for _ in [frac,1]]
                    else:
                        a,b,c=sm.ellipse_68()
                        lines.append("fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                                (ra,dec,degrees(a),degrees(b),degrees(c)+90))


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
