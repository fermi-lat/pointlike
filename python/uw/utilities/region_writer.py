""" Class to write out region files compatable with ds9. 

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/utilities/region_writer.py,v 1.3 2010/10/02 23:41:56 lande Exp $

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
    """ There is the same inconsistency in ellipse definition 
        between extended sources ellipses and ds9 ellipses as
        is discussed in the docstring for unparse_localization,
        resulting in the same switch from maj,min <-> min,maj. """

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
                                      (ra,dec,sigma))
                    elif isinstance(sm,Ring):
                        frac=sm.frac
                        lines += ["fk5; circle(%.4f, %.4f, %.4f)" % \
                                      (ra,dec,_) for _ in [frac*sigma,sigma]]
                    else:    
                        lines.append("fk5; circle(%.4f, %.4f, %.4f) # Circle contaning 68 percent of the source." % \
                                      (ra,dec,sm.r68()))

                elif isinstance(sm,EllipticalSpatialModel):
                    sigma_x, sigma_y, theta = sm.sigma_x, sm.sigma_y, sm.theta
                    if isinstance(sm,EllipticalDisk):
                        lines.append("fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                                (ra,dec,sigma_y,sigma_x,sm.theta))

                    elif isinstance(sm,EllipticalRing):
                        frac = sm.frac
                        lines += ["fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                                (ra,dec,_*sigma_y,_*sigma_x,sm.theta) \
                                for _ in [frac,1]]
                    else:
                        a,b,c=sm.ellipse_68()
                        lines.append("fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                                (ra,dec,b,a,c))

    return lines

def unparse_localization(roi):
    """ Note that maj and min are switched. This is caused by a discrepancy between
        pointlike's localization code and ds9's ellipse definition. It was 
        discussed in http://confluence.slac.stanford.edu/x/JSCJBQ
        The same link is: https://confluence.slac.stanford.edu/display/SCIGRPS/LAT+Catalog+ds9+regions
        """

    if roi.__dict__.has_key('qform'):
        ra,dec,a,b,ang=roi.qform.par[0:5]
        return ["# The next line is the localization error",
                "fk5; ellipse(%.4f, %.4f, %.4f %.4f, %.4f)" % \
                (ra,dec,b,a,ang)]
    else:
        return []

def writeRegion(roi,filename,color='green'):

    lines = [
        "# Region file format: DS9 version 4.0",
        "global color=%s" % color,
    ]

    lines += unparse_diffuse_sources(roi.dsm.diffuse_sources)
    lines += unparse_point_sources(roi.psm.point_sources)
    lines += unparse_localization(roi)

    file=open(filename,'w')
    file.write('\n'.join(lines))
    file.close()
