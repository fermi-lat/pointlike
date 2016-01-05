""" Class to write out region files compatable with ds9. 

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/region_writer.py,v 1.15 2011/06/25 03:54:50 lande Exp $

author: Joshua Lande
"""
from uw.like.roi_extended import ExtendedSource
from uw.like.SpatialModels import *
from math import degrees

def unparse_point(ps,label_sources):
    string="fk5; point(%.4f, %.4f) # point=cross" % \
            (ps.skydir.ra(),ps.skydir.dec())
    if label_sources: string += " text={%s}" % ps.name
    return string

def unparse_point_sources(point_sources,show_sources,label_sources):

    return [unparse_point(ps,label_sources) for ps in point_sources] if show_sources else []


def unparse_extension(spatial_model,extension_color=None,r68=False):
    """ By default, return edge for disk source,
        inner and outer edge for ring sources,
        and otherwise the 68% edge.
        
        If r68 flag is true, always return 68% edge
        of extended source. """
    sm = spatial_model

    extra = 'color=%s' % extension_color if extension_color is not None else '' 

    if isinstance(sm,PseudoSpatialModel) or type(sm) == SpatialMap:
        return []

    ra,dec=sm.center.ra(),sm.center.dec()

    if isinstance(sm,RadiallySymmetricModel):
        sigma=sm.sigma

        if isinstance(sm,Disk) and r68 is False:
            return ["fk5; circle(%.4f, %.4f, %.4f) # %s Circle encloses all of the disk." % \
                          (ra,dec,sigma,extra)]
        elif isinstance(sm,Ring) and r68 is False:
            frac=sm.frac
            return ["fk5; circle(%.4f, %.4f, %.4f) # %s" % \
                          (ra,dec,_,extra) for _ in [frac*sigma,sigma]]
        else:    
            return ["fk5; circle(%.4f, %.4f, %.4f) # %s Circle containing 68 percent of the source." % \
                          (ra,dec,sm.r68(),extra)]

    elif isinstance(sm,EllipticalSpatialModel):
        sigma_x, sigma_y, theta = sm.sigma_x, sm.sigma_y, sm.theta
        if isinstance(sm,EllipticalDisk) and r68 is False:
            return ["fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # %s" % \
                    (ra,dec,sigma_y,sigma_x,sm.theta,extra)]

        elif isinstance(sm,EllipticalRing) and r68 is False:
            frac = sm.frac
            return ["fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # %s" % \
                    (ra,dec,_*sigma_y,_*sigma_x,sm.theta,extra) \
                    for _ in [frac,1]]
        else:
            a,b,c=sm.ellipse_68()
            return ["fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f) # %s" % \
                    (ra,dec,b,a,c,extra)]
    else:
        raise Exception("Unable to Parse Spatial Model %s" % type(sm))

def unparse_diffuse_sources(diffuse_sources,show_sources,label_sources,show_extension,extension_color):
    """ There is the same inconsistency in ellipse definition 
        between extended sources ellipses and ds9 ellipses as
        is discussed in the docstring for unparse_localization,
        resulting in the same switch from maj,min <-> min,maj. """

    lines = []
    for ds in diffuse_sources:
        if isinstance(ds,ExtendedSource):
            sm = ds.spatial_model

            if show_sources: 
                lines.append(unparse_point(ds,label_sources))

            if show_extension:
                lines += unparse_extension(sm,extension_color)


    return lines

def unparse_localization(roi):
    """ Note that maj and min are switched. This is caused by a discrepancy between
        pointlike's localization code and ds9's ellipse definition. It was 
        discussed in http://confluence.slac.stanford.edu/x/JSCJBQ
        The same link is: https://confluence.slac.stanford.edu/display/SCIGRPS/LAT+Catalog+ds9+regions
        """

    if roi.__dict__.has_key('qform'):
        ra,dec=roi.qform.par[0:2]
        a,b,ang=roi.qform.par[3:6]
        return ["# The next line is the localization error",
                "fk5; ellipse(%.4f, %.4f, %.4f, %.4f, %.4f)" % \
                (ra,dec,b,a,ang)]
    else:
        return []

def get_region(roi,color,show_sources=True, label_sources=True,show_localization=True,show_extension=True,extension_color=None):
    lines = [
        "# Region file format: DS9 version 4.0",
        "global color=%s" % color,
    ]

    lines += unparse_diffuse_sources(roi.dsm.diffuse_sources,show_sources,label_sources,show_extension,extension_color)
    lines += unparse_point_sources(roi.psm.point_sources,show_sources,label_sources)
    if show_localization:
        lines += unparse_localization(roi)
    if len(lines)==2: return None
    return '\n'.join(lines)

def writeRegion(roi,filename,color='green',show_sources=True, label_sources=True,show_localization=True, show_extension=True):
    """ Saves out an ROI to a ds9 style region file.
        
        The size of simple exended sources is saved to the region file
        as are elliptical localization errors if they exist. """

    file=open(filename,'w')
    file.write(get_region(roi,color,show_sources,label_sources,show_localization,show_extension))
    file.close()
