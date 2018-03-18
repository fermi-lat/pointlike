"""
Module to save an ROIAnalysis object to a file and to load it back in.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/roi_save.py,v 1.10 2013/05/07 04:18:08 lande Exp $

author: Joshua Lande
"""
import os
import cPickle
import collections
import numpy as N

from uw.utilities import path

class Empty: pass

def save(roi,filename):
    """ Save ROI to file. 
    
        This implementation has many limitations and is suitable only for quick hacks.

            * To keep the saved file small, this none of the model predictions are saved. After
              reloading the ROI, all of these time consuming calculations have to be redone.
            * Is fragile to the exact path of things. So you probably could not reload the ROI
              with a new version of the science tools, or after moving the location of your
              diffuse sources or anything like that.
            * Can only pickle certain kinds of diffuse sources.
            * Does not respect if a special diffuse_mapper was used to create the ROI.

        Nevertheless, it is very useful to be able to temporarily save and load an ROI. """
    self=roi

    d=collections.defaultdict(dict)

    from . pointspec import DataSpecification,SpectralAnalysis

    for i in DataSpecification.defaults:
        if len(i)==3:
            j=i[0]
            d['DataSpecification'][j]=self.sa.dataspec.__dict__[j]

    for i in SpectralAnalysis.defaults:
        if len(i)==3:
            j=i[0]
            d['SpectralAnalysis'][j]=self.sa.__dict__[j]

    d['roi_dir']=self.roi_dir
    d['point_sources']=self.psm.point_sources.tolist()
    d['diffuse_sources']=self.dsm.diffuse_sources.tolist()

    # allow storing extra stuff for custom analysis
    if hasattr(roi,'extra'): d['extra'] = roi.extra

    if self.__dict__.has_key('qform') and \
            self.__dict__.has_key('ldir') and \
            self.__dict__.has_key('lsigma') and \
            self.__dict__.has_key('delta_loc_logl'):

        # can't save out qform, but at least save qform.par
        empty_qform=Empty()
        empty_qform.par=self.qform.par

        # add in localization stuff, kinda ugly
        d['localization']={'qform':empty_qform,
                           'ldir':self.ldir,
                           'lsigma':self.lsigma,
                           'delta_loc_logl':self.delta_loc_logl}

    from . roi_analysis import ROIAnalysis

    for i in ROIAnalysis.defaults:
        if len(i)==3:
            j=i[0]
            d['ROIAnalysis'][j]=self.__dict__[j]

    d['LATEXTDIR']=os.environ['LATEXTDIR'] if os.environ.has_key('LATEXTDIR') else None

    cPickle.dump(d,open(path.expand(filename),'w'))

def load(filename,**kwargs):
    """ Factory method to return a ROIAnalysis object
        that has been saved to a file. 
        
        Any additional kwargs is used to modify DataSpecification, SpectralAnalysis,
        and ROIAnalysis objects."""
    if isinstance(filename, basestring):
        d=cPickle.load(open(path.expand(filename),'r'))
    elif isinstance(filename, dict):
        d=filename
    else:
        raise Exception("Unknown ROI file %s" % filename)

    # restore previous LATEXTDIR if it is not already set
    if not os.environ.has_key('LATEXTDIR') and d['LATEXTDIR'] not in [None,{}]:
        os.environ['LATEXTDIR']=d['LATEXTDIR'] 

    # SpatialMap objects may fail to load if LATEXTDIR wasn't previously
    # set. If so, try again.
    if N.any([hasattr(ds,'spatial_model') and \
              hasattr(ds.spatial_model,'skyfun') \
              and ds.spatial_model.skyfun is None 
              for ds in d['diffuse_sources']]):
        d=cPickle.load(open(path.expand(filename),'r'))

    from . pointspec import DataSpecification,SpectralAnalysis
    from . roi_analysis import ROIAnalysis

    xmlfile=None
    point_sources = d['point_sources']
    diffuse_sources = d['diffuse_sources']

    for k,v in kwargs.items():

        keys=lambda x: [ i[0] for i in x]

        if k in keys(DataSpecification.defaults):
            d['DataSpecification'][k]=v
        elif k in keys(SpectralAnalysis.defaults):
            # allow changing the ROIAnalysis parameters when
            d['SpectralAnalysis'][k]=v
        elif k in keys(ROIAnalysis.defaults):
            d['ROIAnalysis'][k]=v
        elif k == 'xmlfile':
            xmlfile = v
        elif k == 'point_sources':
            point_sources = v
        elif k == 'diffuse_sources':
            diffuse_sources = v
        else:
            raise Exception("Unknown argument %s to function load" % k)

    # backwards compatability
    for ps in point_sources: 
        if hasattr(ps.model,'p'): ps.model._p=ps.model.p
    for ds in diffuse_sources: 
        if hasattr(ds.smodel,'p'): ds.smodel._p=ds.smodel.p

    ds=DataSpecification(**d['DataSpecification'])
    sa=SpectralAnalysis(ds,**d['SpectralAnalysis'])
    roi=sa.roi(roi_dir=d['roi_dir'],
               xmlfile=xmlfile,
               point_sources=point_sources,
               diffuse_sources=diffuse_sources,
               **d['ROIAnalysis'])

    # add in localization stuff, kinda ugly
    if d.has_key('localization'):
        roi.qform=d['localization']['qform']
        roi.ldir=d['localization']['ldir']
        roi.lsigma=d['localization']['lsigma']
        roi.delta_loc_logl=d['localization']['delta_loc_logl']

    # load back any potential custom stuff
    if d.has_key('extra'): roi.extra = d['extra']

    # just to be safe
    roi.__update_state__()

    return roi

