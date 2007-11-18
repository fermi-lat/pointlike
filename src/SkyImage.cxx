/** @file SkyImage.cxx

@brief implement the class SkyImage
$Header: /nfs/slac/g/glast/ground/cvs/map_tools/src/SkyImage.cxx,v 1.55 2007/05/07 18:29:40 burnett Exp $
*/

#include "pointlike/SkyImage.h"
#include "astro/SkyProj.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"

#include <cctype>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <errno.h> // to test result of std::remove()

namespace {
    static unsigned long lnan[2]={0xffffffff, 0x7fffffff};
    static double& dnan = *( double* )lnan;
}
using namespace pointlike;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SkyImage::SkyImage(const astro::SkyDir& center,  
                   const std::string& outputFile, 
                   double pixel_size, double fov, int layers, 
                   const std::string& ptype,
                   bool galactic)
: m_naxis3(layers)  
, m_image(0)
, m_save(true)
, m_layer(0)
, m_interpolate(false)
{

    if( fov>90) {
        std::string types[]={"" ,"CAR","AIT","ZEA"};
        int xsize[] =       {360, 360,  325,  230}; 
        int ysize[] =       {180, 180,  162,  230}; 
        for( unsigned int i = 0; i< sizeof(types)/sizeof(std::string); ++i){
            if( ptype == types[i]) {
                m_naxis1 = static_cast<int>(xsize[i]/pixel_size);
                m_naxis2 = static_cast<int>(ysize[i]/pixel_size);
                break;
            }
        }

        if( m_naxis1==0) {
            throw std::invalid_argument("SkyImage::SkyImage -- projection type " 
                +ptype +" does not have default image size");
        }
    }else{

        m_naxis1=m_naxis2 = static_cast<int>(fov/pixel_size + 0.5);
    }

    double crval[2] = { galactic?center.l():center.ra(),galactic? center.b(): center.dec()};
    double cdelt[2] = { -pixel_size, pixel_size };
    double crpix[2] = { (m_naxis1+1)/2.0, (m_naxis2+1)/2.0};

    m_wcs = new astro::SkyProj(ptype, crpix, crval, cdelt, 0., galactic);
    this->setupImage(outputFile);
}
#if 0
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Note: this constructor was stolen on 2/7/2006 by James Peachey from
// SkyImage::SkyImage(const map_tools::MapParameters&) and modified to use
// ScienceTools-compliant parameters for the map geometry. This constructor
// is thus redundant to the one above, but represents a migration path toward
// making map_tools more like the other ScienceTools.
// TODO: migrate all tools in map_tools to use this constructor, then
// remove Parameters and MapParameters classes and the constructor above.
SkyImage::SkyImage(const hoops::IParGroup& pars)
: m_naxis1(1)
, m_naxis2(1)
, m_naxis3(1)
, m_total(0)
, m_image(0)
, m_imageData()
, m_save(true)
, m_layer(0)
, m_wcs(0)
{
    using namespace astro;

    // see if there is an input count map
    std::string cm_file = pars["cmfile"];
    std::string uc_cm_file = cm_file;
    for ( std::string::iterator itor = uc_cm_file.begin(); itor != uc_cm_file.end(); ++itor) *itor = std::toupper(*itor);
        
    // determine how energy values are computed from bins using bincalc parameter.
    std::string layer_calc = pars["bincalc"];
    for ( std::string::iterator itor = layer_calc.begin(); itor != layer_calc.end(); ++itor) *itor = std::toupper(*itor);

    if ( "NONE" != uc_cm_file){
        // get as much info as possible from the count map
        // determine the projection
        m_wcs = new astro::SkyProj(cm_file,1);

        // read the count map
        std::auto_ptr<const tip::Image> count_map(tip::IFileSvc::instance().readImage(cm_file, ""));
        const tip::Header& header = count_map->getHeader();

        // first two dimensions come from count map image
        header["NAXIS1"].get(m_naxis1);
        header["NAXIS2"].get(m_naxis2);

        // read energies associated with layers from ebounds extension of count map.
        std::auto_ptr<const tip::Table> ebounds(tip::IFileSvc::instance().readTable(cm_file, "EBOUNDS"));

        static double s_MeV_per_keV = .001;

        if ( layer_calc == "CENTER"){
            // compute energies using N central bin values from ebounds
            m_naxis3 = ebounds->getNumRecords();
            m_energy.resize(m_naxis3);
            std::vector<double>::iterator out_itor = m_energy.begin();
            for ( tip::Table::ConstIterator in_itor = ebounds->begin(); in_itor != ebounds->end(); ++in_itor, ++out_itor){
                double e_min = (*in_itor)["E_MIN"].get();
                double e_max = (*in_itor)["E_MAX"].get();
                *out_itor = .5 * (e_max + e_min) * s_MeV_per_keV;
            }
        } else {
            // compute energies using N + 1 edge bin values from ebounds
            m_naxis3 = ebounds->getNumRecords() + 1;
            m_energy.resize(m_naxis3);
            tip::Table::ConstIterator last_in = ebounds->begin();
            std::vector<double>::iterator out_itor = m_energy.begin();
            for ( tip::Table::ConstIterator in_itor = ebounds->begin(); in_itor != ebounds->end(); ++in_itor, ++out_itor){
                *out_itor = (*in_itor)["E_MIN"].get() * s_MeV_per_keV;
                last_in = in_itor;
            }
            // get last bin edge from e_max column.
            *out_itor = (*last_in)["E_MAX"].get() * s_MeV_per_keV;
        }

        // size of image is now known, so initialize it.
        m_pixelCount = m_naxis1*m_naxis2*m_naxis3;
        m_imageData.resize(m_pixelCount, 0.);
    }else{
        //
        // the input is a livetime cube: get display from par file parameters
        //
        m_naxis1 = pars["numxpix"];
        m_naxis2 = pars["numypix"];
        int enumbins = pars["enumbins"];
        std::string ptype = pars["proj"];
        double pixelsize = pars["pixscale"];

        if( m_naxis1<=1){
            // special code to determine all-sky limits based on scale factor and transformation
            std::string types[]={"" ,"CAR","AIT","ZEA"};
            int xsize[] =       {360, 360,  325,  230}; 
            int ysize[] =       {180, 180,  162,  230}; 
            for( unsigned int i = 0; i< sizeof(types)/sizeof(std::string); ++i){
                if( ptype == types[i]) {
                    m_naxis1 = static_cast<int>(xsize[i]/pixelsize);
                    m_naxis2 = static_cast<int>(ysize[i]/pixelsize);
                    break;
                }
            }
            if( m_naxis1<=1) {
                throw std::invalid_argument("SkyImage::SkyImage -- projection type " 
                    +ptype +" does not have default image size");
            }
        }
        if( m_naxis2==0) m_naxis2=m_naxis1; // default square image
        std::string coord_sys = pars["coordsys"];
        for (std::string::iterator itor = coord_sys.begin(); itor != coord_sys.end(); ++itor) *itor = std::toupper(*itor);
        bool galactic = (coord_sys == "GAL");

        /// arrays describing transformation: assume reference in the center
        double          //lon            lat
            crval[2]={ pars["xref"],      pars["yref"]},
            crpix[2]={ (m_naxis1+1)/2.0, (m_naxis2+1)/2.0},
            cdelt[2]={ -pixelsize,       pixelsize },
            crota2=pars["axisrot"];
        m_wcs = new astro::SkyProj( pars["proj"], crpix, crval, cdelt, crota2, galactic);
  
        double emin = pars["emin"], emax = pars["emax"];

        // compute logarithmic bin ratio for edges
        std::vector<double> edge(enumbins+1);
        const double eratio = std::exp(std::log(emax / emin) / enumbins);
        double energy = emin;
        for ( int ii = 0; ii != enumbins; ++ii, energy *= eratio){
            edge[ii] = energy;
        }
        // prevent annoying round-off in the last bin
        edge[enumbins] = emax;

        // handle different styles of energy output
        if ( layer_calc == "CENTER"){
          // energies are taken at centers of bins
          m_naxis3=enumbins;
          m_energy.resize(m_naxis3);
          for( int ii = 0; ii != enumbins; ++ii){
            m_energy[ii] = .5 * (edge[ii] + edge[ii+1]);
          }
        }else{
          // energies are taken at edges of bins
          m_naxis3=enumbins+1;
          m_energy=edge;
        }
    }
    setupImage(pars["outfile"],  pars["clobber"]);
}
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::setupImage(const std::string& outputFile,  bool clobber)
{
    std::string extension("skyimage"); // maybe a parameter?

    if( clobber ){
        int rc = std::remove(outputFile.c_str());
        if( rc==-1 && errno ==EACCES ) throw std::runtime_error(
            std::string("SkyImage: cannot remove file "+outputFile)
            );
    }

    // setup the image: it needs an axis dimension array and the file name to write to
    std::vector<long> naxes(3);
    naxes[0]=m_naxis1;
    naxes[1]=m_naxis2;
    naxes[2]=m_naxis3;

    // now add an image to the file
    tip::IFileSvc::instance().appendImage(outputFile, extension, naxes);
    // create a float image
    m_image = tip::IFileSvc::instance().editImageFlt(outputFile, extension);

    m_pixelCount = m_naxis1*m_naxis2*m_naxis3;
    m_imageData.resize(m_pixelCount);

    // fill the boundaries with NaN
    //if( pars.projType()!="CAR") clear();

    m_wcs->setKeywords(m_image->getHeader());
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SkyImage::SkyImage(const std::string& fits_file, const std::string& extension, bool interpolate)
: m_save(false)
, m_layer(0)
, m_wcs(0)
, m_interpolate(interpolate)
{
    // note expect the image to be float
    m_image = tip::IFileSvc::instance().editImageFlt(fits_file, extension);
    tip::Header& header = m_image->getHeader();

    // standard ordering for ra, dec, cos(theta).
    header["NAXIS1"].get(m_naxis1);
    header["NAXIS2"].get(m_naxis2);
    header["NAXIS3"].get(m_naxis3);
    m_pixelCount = m_naxis1*m_naxis2*m_naxis3;

    m_wcs = new astro::SkyProj(fits_file,1);
    // finally, read in the image: assume it is float
    dynamic_cast<tip::TypedImage<float>*>(m_image)->get(m_imageData);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unsigned int SkyImage::setLayer(unsigned int newlayer)
{
    checkLayer(newlayer);
    unsigned int t = m_layer;
    m_layer = newlayer;
    return t;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool SkyImage::addPoint(const astro::SkyDir& dir, double delta, unsigned int layer)
{
    std::pair<double,double> p= dir.project(*m_wcs);
    // ignore if not in the image.
    if( p.first<0 || p.first >= m_naxis1 || p.second<0 || p.second>=m_naxis2) return false;
    unsigned int 
        i = static_cast<unsigned int>(p.first),
        j = static_cast<unsigned int>(p.second),
        k = i+m_naxis1*(j + layer*m_naxis2);
    
    if(  k< m_pixelCount){
        m_imageData[k] += delta;
        m_total += delta;
    }
    return true;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::checkLayer(unsigned int layer)const
{
    if( layer >= (unsigned int)m_naxis3){
        std::stringstream errmsg;
        errmsg << "SkyImage: requested layer " << layer 
            << " not compatible with axis3: " << m_naxis3 << std::endl;
        throw std::out_of_range(errmsg.str());
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::fill(const astro::SkyFunction& req, unsigned int layer)
{
    checkLayer(layer);
    m_total=m_count=m_sumsq=0;
    m_min=1e10;m_max=-1e10;
    int offset = m_naxis1* m_naxis2 * layer;
    for( size_t k = 0; k< (unsigned int)(m_naxis1)*(m_naxis2); ++k){
        double 
            x = static_cast<int>(k%m_naxis1)+1.0, 
            y = static_cast<int>(k/m_naxis1)+1.0;
        if( m_wcs->testpix2sph(x,y)==0) {
            astro::SkyDir dir(x,y, *m_wcs);
            double t= req(dir);
            m_imageData[k+offset] = t;
            m_total += t;
            ++m_count;
            m_sumsq += t*t;
            m_min = t<m_min? t:m_min;
            m_max = t>m_max? t:m_max;
        }else{
            // not valid (off the edge, perhaps)
            m_imageData[k+offset]=dnan; 
        }
    }
    return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::clear()
{
    size_t s = m_imageData.size();
    for( size_t k = 0; k< s; ++k){
        m_imageData[k]=0; 
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SkyImage::~SkyImage()
{
    if( m_save) {
        dynamic_cast<tip::TypedImage<float>*>(m_image)->set(m_imageData);
    }
    delete m_image; 
    delete m_wcs;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SkyImage::pixelValue(const astro::SkyDir& pos,unsigned  int layer)const
{
    checkLayer(layer); 
    double v(0) ;
    
    if( m_interpolate ){
        // interpolating
        // project using wcslib interface
        std::pair<double,double> p= pos.project(*m_wcs);
        double x( floor(p.first-0.5) ),
               y( floor(p.second-0.5) ),
               dx( p.first-x -1 ),
               dy( p.second-y-1 );
        unsigned int k = static_cast<unsigned int>(x + m_naxis1*(y+layer*m_naxis2));
        v = m_imageData[k];
        if( dx>0. && x < m_naxis1){ 
            v = v*(1.-dx) + m_imageData[k+1]*dx;
        }else if( x >1) {
            v = v*(1+dx) - m_imageData[k-1]*dx;
        }
        if( dy>0. && y < m_naxis2){ 
            v = v*(1.-dy) + m_imageData[k+m_naxis1]*dy;
        }else if( y > 1) {
            v = v*(1.+dy) - m_imageData[k-m_naxis1]*dy;
        }
    }else{
        unsigned int k = pixel_index(pos,layer);
        v = m_imageData[k];
    }
   
    return v;        
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
float &  SkyImage::operator[](const astro::SkyDir&  pixel)
{
    unsigned int k = pixel_index(pixel);
    return m_imageData[k];        

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
const float &  SkyImage::operator[](const astro::SkyDir&  pixel)const
{
    unsigned int k = pixel_index(pixel);
    return m_imageData[k];        

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SkyImage::operator()(const astro::SkyDir& s)const
{
    return pixelValue(s, m_layer);        
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SkyImage::getNeighbors(const astro::SkyDir& pos, std::vector<double>&neighbors)const
{
    int layer = 0; ///@todo: get neighbors on a different layer
    std::pair<double,double> p= pos.project(*m_wcs);
    if( p.first<0) p.first += m_naxis1;
    if(p.second<0) p.second += m_naxis2;
    unsigned int 
        i = static_cast<unsigned int>(p.first-0.5),
        j = static_cast<unsigned int>(p.second-0.5),
        k = i+m_naxis1*(j + layer*m_naxis2);
    if(i+1<(unsigned int)m_naxis1)neighbors.push_back(m_imageData[k+1]); 
    if(i>0) neighbors.push_back(m_imageData[k-1]);
    if(j+1<(unsigned int)m_naxis2)neighbors.push_back(m_imageData[k+m_naxis1]);
    if(j>0)neighbors.push_back(m_imageData[k-m_naxis1]);

}

// internal routine to convert a SkyDir to a pixel index
unsigned int SkyImage::pixel_index(const astro::SkyDir& pos, int layer) const
{
    // if not specified, use the data member
    if( layer<0 ) layer = m_layer;

    // project using wcslib interface, then adjust to be positive
    std::pair<double,double> p= pos.project(*m_wcs);
    if( p.first<0) p.first += m_naxis1;
    if(p.second<0) p.second += m_naxis2;
    unsigned int 
        i = static_cast<unsigned int>(p.first-0.5),
        j = static_cast<unsigned int>(p.second-0.5),
        k = i+m_naxis1*(j + layer*m_naxis2);
     if( k > m_pixelCount ) {
        throw std::range_error("SkyImage::pixel_index -- outside image hyper cube");
    }
    return k;
}


