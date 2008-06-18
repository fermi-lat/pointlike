/** @file likelihood_main.cxx


*/
#include "pointlike/SourceLikelihood.h"
#include "embed_python/Module.h"

#include "skymaps/PhotonMap.h"
#include "skymaps/DiffuseFunction.h"
#include "pointlike/TSmap.h"
#include "pointlike/Data.h"
#include "pointlike/DrawTS.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

// list of starting points, corresponding to data
class Points{public:
    std::string name; double ra, dec;
}points[]={
    {"one", 0, 0},
    {"two", 20, 0},
    {"three", 40, 0},
    {"",0,0}
};

using namespace astro;
using namespace pointlike;
using skymaps::PhotonMap;

/** @class AddPhoton

*/
class AddPhoton: public std::unary_function<astro::Photon, void> {
public:
    AddPhoton (skymaps::PhotonMap& map)
        : m_map(map)
    {}
    void operator()(const astro::Photon& gamma)
    {
        //int event_class = gamma.eventClass();

        m_map.addPhoton(gamma);
    }
    skymaps::PhotonMap& m_map;
};



int main(int argc , char** argv ){

  int rc(0);
  try{

    std::string python_path("../python");

    if( argc>1){
      python_path = argv[1];
    }
    
    embed_python::Module setup(python_path , "sourcefit_setup",  argc, argv);
    std::vector<double> ras, decs;
    std::vector<std::string> names;
    double fov = 10;
    double pixelSize = 0.1;

    setup.getList("name", names);
    setup.getList("ra", ras);
    setup.getList("dec", decs);
    setup.getValue("fov", fov, fov);
    setup.getValue("pixelsize", pixelSize, pixelSize);
    
    // use the  Data class to create the PhotonData object
    pointlike::Data healpixdata(setup);

    std::string path( ::getenv("EXTFILESSYS"));
    std::string diffusefile;
    setup.getValue("Diffuse.file", diffusefile);
    
//     skymaps::DiffuseFunction df1( path + "/galdiffuse/GP_gamma.fits");
    skymaps::DiffuseFunction df1( diffusefile);
    pointlike::DrawTS  draw(healpixdata, df1);

    double x = ras[0]; // For the moment, take the first ra and the first dec ...
    double y = decs[0];
    astro::SkyDir dirNew = astro::SkyDir(x, y, astro::SkyDir::EQUATORIAL);
//     pointlike::TSmap tsmap(healpixdata, df1);
//     tsmap.run(dirNew, 10, 8, true);

    std::cout << "Drawing region for : " << x 
	      <<" " << y << " " <<  pixelSize << " " << fov << std::endl;
    draw.region(dirNew, "testDraw.fits", pixelSize, fov, true);
    
  }catch(const std::exception& e){
    std::cerr << "Caught exception " << typeid(e).name() 
	      << " \"" << e.what() << "\"" << std::endl;
    rc=1;
  }
  return rc;
}

