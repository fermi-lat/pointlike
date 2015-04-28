"""  A module to handle finding irfs
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/pycaldb.py,v 1.11 2015/03/10 16:59:18 burnett Exp $

    author: Joshua Lande """
import os
from os.path import join

import pyfits
import numpy as N

from uw.utilities import keyword_options

class CALDBManager(object):
    """ This object abstracts access to IRF files which can be
        found in the file system. The intention is to provide a
        consistency between the pointlike and gtlike IRF loading.

        In particular, the caldb.indx is read so that a consistent
        naming convention can be used. Additionally, the CUSTOM_IRF_DIR
        environment varaible can used to define new IRFs. """

    defaults=(
        ('psf_irf',None,'specify a different IRF to use for the PSF'),
        ('CALDB',None,'override the CALDB specified by the env. variable or by py_facilities'),
        ('custom_irf_dir',None,'override the irf dir specified by the env. variable CUSTOM_IRF_DIR'),
        ('quiet', True, 'set false for output'),
   )

    @keyword_options.decorate(defaults)
    def __init__(self,irf,**kwargs):

        keyword_options.process(self, kwargs)
        self.irf = irf

        if self.CALDB is None:
            try:
                self.CALDB=os.environ['CALDB']
            except:
                try:
                    from facilities import py_facilities
                    os_environ = py_facilities.commonUtilities_getEnvironment
                    self.CALDB=os_environ('CALDB')
                except:
                    raise Exception('Environment variable CALDB must be set, or findable by py_facilities package')

        if self.custom_irf_dir is not None:
            if not os.path.exists(self.custom_irf_dir):
                raise Exception("custom_irf_dir %s does not exist" % self.custom_irf_dir)
        else:
            self.custom_irf_dir = os.environ.get('CUSTOM_IRF_DIR', None)
        if self.custom_irf_dir is not None and self.custom_irf_dir!='':
            if not self.quiet: print 'CALDBManager: using custom irf: "%s"' % self.custom_irf_dir
            
        self.bcf = join(self.CALDB,'bcf')

        if not os.path.exists(self.bcf):
            self.CALDB = join(self.CALDB,'data', 'glast', 'lat')
            self.bcf = join(self.CALDB,'bcf')

            if not os.path.exists(self.bcf):
                raise Exception('Invalid CALDB directory %s.' % self.bcf)

        self.CALDB_index=os.path.join(self.CALDB,'caldb.indx')

        if not os.path.exists(self.CALDB_index):
            raise Exception("caldb.indx file %s does not exist." % self.CALDB_index)

        self.load_caldb_indx()
        self.construct_psf()
        if not self.quiet: print 'PSF: %s' % self.psf_files
        self.construct_aeff()
        #self.construct_edisp()

    def load_caldb_indx(self):
        # parse CALDB_index
        index_fits = pyfits.open(self.CALDB_index)
        self.irf_names=index_fits[1].data.field('CAL_CBD') # name of the irfs (includes other junk in the string)
        self.conv_types=index_fits[1].data.field('DETNAM') # Whether irf is front/back
        self.irf_files=index_fits[1].data.field('CAL_FILE') # name of fits file for the irf
        self.irf_types=index_fits[1].data.field('CAL_CNAM') # type of irf "RPSF", "EFF_AREA", "EDISP"
        index_fits.close()

    def construct_psf(self):

        # see if the psf should be overloaded
        if self.psf_irf is None:
            irf=self.irf
        else:
            if not self.quiet: print 'Overriding default PSF; using %s'%(self.psf_irf)
            irf=self.psf_irf

        # see if the psf name is directly in the filename.
        self.psf_files = [join(self.bcf,'psf','psf_%s_%s.fits'%(irf,i)) for i in ['front','back']]

        if os.path.exists(self.psf_files[0]) and os.path.exists(self.psf_files[1]):
            return

        # try to read form caldb index
        front=N.where((N.chararray.find(self.irf_names,irf)!=-1)&(self.irf_types=='RPSF')&(self.conv_types=='FRONT'))[0]
        back=N.where((N.chararray.find(self.irf_names,irf)!=-1)&(self.irf_types=='RPSF')&(self.conv_types=='BACK'))[0]

        # if front & back exist in the caldb.indx
        if len(front)==1 and len(back)==1:
            psfs=self.irf_files[front[0]],self.irf_files[back[0]]
            self.psf_files = [join(self.bcf,'psf',x) for x in psfs]

            if os.path.exists(self.psf_files[0]) and os.path.exists(self.psf_files[1]):
                return
            print 'caldb.indx: did not find both %s' %self.psf_files
            
        # try the cusom_irf_dir
        if self.custom_irf_dir is not None and self.custom_irf_dir!='':
            if not os.path.exists(self.custom_irf_dir):
                raise Exception("custom_irf_dir '%s' does not exist." % self.custom_irf_dir)
            self.psf_files = [join(self.custom_irf_dir,'psf_%s_%s.fits' % (irf,i)) for i in ['front','back']]

            if os.path.exists(self.psf_files[0]) and os.path.exists(self.psf_files[1]):
                return
            print 'Custom irf_dir: did not find both %s ' %self.psf_files

        raise Exception("Unable to find the irf %s\n\tcustom_irf_dir: %s\n\tlooking for %s" \
              % (irf, self.custom_irf_dir, self.psf_files) )
    
    def construct_aeff(self):
        irf = self.irf

        self.aeff_files = [join(self.bcf,'ea','aeff_%s_%s.fits'%(irf,i)) for i in ['front','back']]

        if os.path.exists(self.aeff_files[0]) and os.path.exists(self.aeff_files[1]):
            return

        # try to read form caldb index
        front=N.where((N.chararray.find(self.irf_names,irf)!=-1)&(self.irf_types=='EFF_AREA')&(self.conv_types=='FRONT'))[0]
        back=N.where((N.chararray.find(self.irf_names,irf)!=-1)&(self.irf_types=='EFF_AREA')&(self.conv_types=='BACK'))[0]

        # if front & back exist in the caldb.indx
        if len(front)==1 and len(back)==1:
            aeffs=self.irf_files[front[0]],self.irf_files[back[0]]
            self.aeff_files = [join(self.bcf,'ea',x) for x in aeffs]

            if os.path.exists(self.aeff_files[0]) and os.path.exists(self.aeff_files[0]):
                return
            
        # try the cusom_irf_dir
        if self.custom_irf_dir is not None:
            self.aeff_files = [join(self.custom_irf_dir,'aeff_%s_%s.fits' % (irf,i)) for i in ['front','back']]

            if os.path.exists(self.aeff_files[0]) and os.path.exists(self.aeff_files[0]):
                return

        raise Exception("Unable to find the irf %s." % irf)


    def construct_edisp(self):
        raise NotImplementedException("No current need for this function.")

    def get_psf(self): return self.psf_files

    def get_aeff(self): 
        return self.aeff_files

    def get_edisp(self):
        raise NotImplementedException("No current need for this function.")
