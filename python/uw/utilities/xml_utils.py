""" xml_utils.py
    author: Damien Parent <dmnparent@gmail.com>

    parameter_element:
      Create a XML document parameter description element
"""
import os, sys, types
from xml.dom import minidom
from optparse import OptionParser
from uw.utilities.coords import sepangle_deg

def parameter_element(free, name, maximum, minimum, scale, value):
    """Create an XML document parameter description element"""
    impl = minidom.getDOMImplementation()
    xmldoc_out = impl.createDocument(None,None,None)
    parameter = xmldoc_out.createElement('parameter')
    parameter.setAttribute('free', str(free))
    parameter.setAttribute('name', str(name))
    parameter.setAttribute('max', str(maximum))
    parameter.setAttribute('min', str(minimum))
    parameter.setAttribute('scale', str(scale))
    parameter.setAttribute('value', str(value))
    return parameter

class xml_manager():

    def init(self):
        self.galactic         = None
        self.galactic_name    = 'gal_2yearp7v6_v0'
        self.isotropic        = None
        self.isotropic_name   = 'iso_p7v6source'
        self.template_dir     = '/afs/slac/g/glast/groups/pulsar/share/Templates'
       
    def __init__( self, catalogs, **kwargs ):        
        """
        ===========     ===============================================
        keyword         description
        ===========     ===============================================
        catalogs        XML catalog (single file or list)
        galactic        Galactic Diffuse model                   [None]
        isotropic       Isotropic Model                          [None]
        galactic_name   [gal_2yearp7v6_v0]
        isotropic_name  [iso_p7v6source]
        template_dir    [/afs/slac/g/glast/groups/pulsar/share/Templates]
        """

        self.catalogs = catalogs
        self.init()
        self.__dict__.update(kwargs)

        if isinstance(catalogs,types.StringType):
            if not os.path.isfile(catalogs): print ("Cannot open the catalog! Exiting ..."); sys.exit()
            self.sourcelist = minidom.parse(catalogs).getElementsByTagName('source')
        elif isinstance(catalogs,types.ListType):
            self.sourcelist = []
            for catalog in catalogs:
                if not os.path.isfile(catalog): print ("Cannot open %s! Exiting ..." %catalog)
                self.sourcelist += minidom.parse(catalog).getElementsByTagName('source')
        else:
            raise Exception("Do not recognize the format! Exiting...")

        # exception for diffuse sources
        self.DiffuseSourceList = []
        self.DiffuseSourceList += [{'name':'CenALobes.fits','ra':201.0,'dec':-43.5}]
        self.DiffuseSourceList += [{'name':'CygnusLoop.fits','ra':312.75,'dec':30.67}]
        self.DiffuseSourceList += [{'name':'IC443.fits','ra':94.31,'dec':22.58}]
        self.DiffuseSourceList += [{'name':'HESSJ1825-137.fits','ra':276.13,'dec':-13.8521}]
        self.DiffuseSourceList += [{'name':'MSH15-52.fits','ra':228.507,'dec':-59.256}]
        self.DiffuseSourceList += [{'name':'LMC.fits','ra':81.65,'dec':-68.42}]
        self.DiffuseSourceList += [{'name':'SMC.fits','ra':14.75,'dec':-72.7}]
        self.DiffuseSourceList += [{'name':'VelaX.fits','ra':128.287,'dec':-45.1901}]
        self.DiffuseSourceList += [{'name':'W28.fits','ra':270.34,'dec':-23.44}]
        self.DiffuseSourceList += [{'name':'W30.fits','ra':271.408,'dec':-21.6117}]
        self.DiffuseSourceList += [{'name':'W44.fits','ra':283.99,'dec':1.355}]
        self.DiffuseSourceList += [{'name':'W51C.fits','ra':290.818,'dec':14.145}]
                                
    def fill_srclist( self, ra=None, dec=None, max_roi=20, free_roi=5, tsmin=0, **kwargs ):
        """
        Create a xml source model from a xml catalog
        ===========     ===============================================
        keyword         description
        ===========     ===============================================
        ra              Right Ascension (deg)                    [None]
        dec             Declination (deg)                        [None]
        max_roi         Maximum radius (deg)                     [20]
        free_roi        ROI or which param are free (deg)        [5]
        tsmin           TS threshold                             [0]
        """

        self.__dict__.update(kwargs)
        self.srclist = []

        impl = minidom.getDOMImplementation()
        xmldoc_out = impl.createDocument(None, "source_library", None)
        xmldoc_out.documentElement.setAttribute('title', "source library")

        # ================================================================  
        # ====================== DIFFUSE BACKGROUND ======================
        # ================================================================
        # GALACTIC
        if self.galactic is not None:
            source_out = xmldoc_out.createElement('source')
            source_out.setAttribute('name', self.galactic_name)
            source_out.setAttribute('type', "DiffuseSource")
            
            spectrum_out = xmldoc_out.createElement('spectrum')
            spectrum_out.setAttribute('type', "PowerLaw")
            spectrum_out.appendChild(parameter_element("1", "Prefactor", "100.0", "1e-5", "1", "1."))
            spectrum_out.appendChild(parameter_element("1", "Index", "1", "-1", "1", "0."))
            spectrum_out.appendChild(parameter_element("0", "Scale", "2000", "50", "1", "500"))
            source_out.appendChild(spectrum_out)
            
            spatial_out = xmldoc_out.createElement('spatialModel')
            spatial_out.setAttribute('type', "MapCubeFunction")
            spatial_out.setAttribute('file', self.galactic)
            spatial_out.appendChild(parameter_element("0", "Normalization", "1000.0", "0.001", "1", "1"))
            source_out.appendChild(spatial_out)
            
            xmldoc_out.documentElement.appendChild(source_out)

        # ISOTROPIC
        if self.isotropic is not None:
            source_out = xmldoc_out.createElement('source')
            source_out.setAttribute('name', self.isotropic_name)
            source_out.setAttribute('type', "DiffuseSource")
            
            spectrum_out = xmldoc_out.createElement('spectrum')
            spectrum_out.setAttribute('type', "FileFunction")
            spectrum_out.setAttribute('file', self.isotropic)
            spectrum_out.appendChild(parameter_element("1", "Normalization", "1000.0", "1e-05", "1.0", "1.0"))
            source_out.appendChild(spectrum_out)
            
            spatial_out = xmldoc_out.createElement('spatialModel')
            spatial_out.setAttribute('type', "ConstantValue")
            spatial_out.appendChild(parameter_element("0", "Value", "10.0", "0.0", "1.0", "1.0"))
            source_out.appendChild(spatial_out)
            
            xmldoc_out.documentElement.appendChild(source_out)

        # Limb smooth
        '''
        if self.limb is not None:
            source_out = xmldoc_out.createElement('source')
            source_out.setAttribute('name', "LIMB")
            source_out.setAttribute('type', "DiffuseSource")
        
            spectrum_out = xmldoc_out.createElement('spectrum')
            spectrum_out.setAttribute('type', "ConstantValue")
            # spectrum_out.setAttribute('file', limb_spec_filename)
            spectrum_out.appendChild(parameter_element("0", "Value", "10", "0.1", "1.0", "1"))
            source_out.appendChild(spectrum_out)
        
            spatial_out = xmldoc_out.createElement('spatialModel')
            spatial_out.setAttribute('type', "MapCubeFunction")
            spatial_out.setAttribute('file', limb_filename)
            spatial_out.appendChild(parameter_element("0", "Normalization", "1e3", "1e-3", "1", "1"))
            source_out.appendChild(spatial_out)
        
            xmldoc_out.documentElement.appendChild(source_out)    
        '''

        # ====================== SOURCE LIST ======================
        #    add the nearby sources from the LAT xml catalog
        # =========================================================        
        nbsrc = 0

        for source in self.sourcelist:

            phflux, err_phflux = 0, 0
            eflux, err_eflux, ts_value = 0, 0, 0
            ra_xml, dec_xml = 0, 0
            
            srcname = source.getAttribute('name')
            if source.getAttribute('type'): type = source.getAttribute('type')
            if source.getAttribute('TS_value'): ts_value = float(source.getAttribute('TS_value'))
            if source.getAttribute('Energy_Flux100'): eflux = source.getAttribute('Energy_Flux100')
            if source.getAttribute('Unc_Energy_Flux100'): err_eflux = source.getAttribute('Unc_Energy_Flux100')
            if source.getAttribute('Flux100'): phflux = source.getAttribute('Flux100')
            if source.getAttribute('Unc_Flux100'): err_phflux = source.getAttribute('Unc_Flux100')
            if source.getAttribute('RA'): ra_xml = float(source.getAttribute('RA'))
            if source.getAttribute('DEC'): dec_xml = float(source.getAttribute('DEC'))
            
            spectrumList = source.getElementsByTagName('spectrum')
            specparamlist = spectrumList[0].getElementsByTagName('parameter')
            
            spatialList = source.getElementsByTagName('spatialModel')
            spatialParamList = spatialList[0].getElementsByTagName('parameter')
            
            for spparam in spatialParamList:
                spparam_name = str(spparam.getAttribute('name'))
                spparam_value = float(spparam.getAttribute('value'))
                if spparam_name == 'RA': ra_xml = spparam_value
                if spparam_name == 'DEC': dec_xml = spparam_value
      
            if type == 'DiffuseSource':
                template = spatialList[0].getAttribute('file')
                template_filename = os.path.join(self.template_dir,template)
                if os.path.isfile(template_filename):
                    spatialList[0].setAttribute('file',template_filename)
                    for src in self.DiffuseSourceList:
                        if os.path.basename(template) == src['name']:
                            ra_xml, dec_xml = src['ra'], src['dec']
                else: ra_xml, dec_xml = -1, -1
                if srcname == self.galactic_name: ra_xml, dec_xml = 0, 0
                if srcname == self.isotropic_name: ra_xml, dec_xml = 0, 0
                            
            # =====================================
            # check if the source is within MAX_ROI
            # =====================================
            angsep = sepangle_deg(ra,dec,ra_xml,dec_xml)
            if (angsep < max_roi or (ra_xml==0 and dec_xml==0)) and ts_value >= tsmin:

                self.srclist += [{'name':srcname,'ra':ra_xml,'dec':dec_xml,'angsep':angsep,'ts':ts_value,
                                  'type':type,'phflux':float(phflux),'err_phflux':float(err_phflux),
                                  'eflux':float(eflux),'err_eflux':float(err_eflux)}]

                for specparam in specparamlist:
                    specparam_name = specparam.getAttribute('name')

                    if specparam.getAttribute('error'): specparam.removeAttribute("error")
                    
                    # if the source is out of ROI_INT, the param are not free = 0
                    if angsep > free_roi: specparam.setAttribute('free', "0")

                    # if the source is within FREE_ROI all the param are free = 1                
                    else:                                                                        
                        mask = ( specparam_name == 'Integral' or specparam_name == 'Prefactor' or specparam_name == 'Index' or
                                 specparam_name == 'norm' or specparam_name == 'alpha' or specparam_name == 'Index1' or
                                 specparam_name == 'Cutoff' )                    
                        if mask: specparam.setAttribute('free', "1")

                        # exception: fit looks more stable
                        if specparam_name == 'beta': specparam.setAttribute('free', "0")

                xmldoc_out.documentElement.appendChild(source)

        # sort srclist by angsep
        self.srclist.sort()

        # save
        self.xmldoc_out = xmldoc_out

    def write_srcmdl(self,filename='srcmdl.xml'):
        outfile = open(filename,'w')
        xmlStr = self.xmldoc_out.toprettyxml('  ').splitlines(True)
        outStr = filter(lambda xmlStr: len(xmlStr) and not xmlStr.isspace(), xmlStr)
        outfile.write(''.join(outStr))
        outfile.close()

    def write_region(self,filename='ds9.reg'):        
        file = open(filename,'w')
        file.write("global color=white \n")
        file.write("fk5 \n")
        for src in self.get_srclist():
            file.write("point(%.2f,%.2f) # point=cross text={%s} \n" %(src['ra'],src['dec'],src['name']))
        file.close()    
	
    def get_srclist(self):
        return self.srclist

    def get_srcname(self,which=0):
        return self.srclist[which]['name']

    def get(self,which=0,key='name'):
        return self.srclist[which][key]

    def print_srclist(self):
        print ("====================================")
        print ("Number of sources =", len(self.srclist))
        print ("====================================")        
        print ("# NAME \t ANGSEP(deg) \t TS \t type")
        for src in self.srclist:
            print ("%s | %.2f | %.0f | %s" %(src['name'],src['angsep'],src['ts'],src['type']))
                            
