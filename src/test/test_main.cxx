/** @file likelihood_main.cxx


*/
#include "pointlike/PointSourceLikelihood.h"
#include "embed_python/Module.h"

#include "skymaps/PhotonMap.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/Exposure.h"

#include "pointlike/SimpleTSmap.h"
#include "pointlike/ParamOptimization.h"
#include "PhotonList.h"

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

// test data file, generated using obsSim with 3 sources
std::string inputFile(  "../src/test/test_events.root" );

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



int main(int argc , char** argv )
{

    int rc(0);
    try{

#if 0 // test code for SimpleTSmap
        std::string path( ::getenv("EXTFILESSYS"));
        DiffuseFunction df1( path + "/galdiffuse/GP_gamma.fits");
        PhotonMap pmap("f:\\glast\\data\\SC2\\obssim\\allsky_noGRBs.fits","PHOTONMAP");
        SimpleTSmap tsmap(pmap, df1);
        tsmap.run();
    
#endif
        // use the python module to setup
        std::string python_path("../python");
        
        PointSourceLikelihood::setParameters(embed_python::Module(python_path , "test_pointlike_setup",  argc, argv));

        double  radius(10);

 
        skymaps::PhotonMap x;
        std::cout << "Loading data from file " << inputFile  <<std::endl;
        PhotonList photons(inputFile);

        std::for_each(photons.begin(), photons.end(), AddPhoton(x));
        std::cout << "photons found: "<< x.photonCount() 
            <<"  pixels created: " << x.pixelCount() <<std::endl;

        x.write("pointlike_test.fits");
        std::vector<astro::SkyDir> directions;



        for( int n=0; !points[n].name.empty(); ++n){
                        
            astro::SkyDir dir(points[n].ra,points[n].dec);
            PointSourceLikelihood like(x, points[n].name, dir);
            // test setting background expectation
#if 0
            std::vector<double> background(8,1E5);
            like.setBackgroundDensity(background);
#endif
            like.set_verbose(true);
            like.printSpectrum(); 

            if( like.TS()>10) {
               like.localize();
            }
            // check function value at peak
            double value = like(dir);
            directions.push_back(like.dir());

        };

        int minlevel(6), maxlevel(13);
        ParamOptimization so(x,directions,&std::cout,minlevel,maxlevel);
        so.compute();

    }catch(const std::exception& e){
        std::cerr << "Caught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        rc=1;
    }
     return rc;
}

/**
@page test Test Program output
Note that the signal fraction are not quite right: these were generated with DC1A, not DC2
@verbatim
Loading data from file ../src/test/test_events.root
photons found: 14396  pixels created: 3226

Spectrum of source one at ra, dec=0, 0
level events  sig fraction  TS
    6  4768  0.70 +/- 0.01  4641
    7  2923  0.64 +/- 0.02  2315
    8  1362  0.63 +/- 0.02  1051
    9   592  0.63 +/- 0.03   446
   10   227  0.72 +/- 0.05   246
   11   128  0.76 +/- 0.07   154
   12    48  0.80 +/- 0.11    61
   13    27  0.71 +/- 0.17    22
                            8934
      Searching for best position, skipping first 0 layers
Gradient   delta     ra        dec       error     Ts
   35797  0.0062    0.0000    0.0000    0.0032    8934.3773
    1861  0.0003    0.0062    0.0000    0.0030    8967.6785
    *** good fit ***

Spectrum of source two at ra, dec=20.000000, 0.000000
level events  sig fraction  TS
    6   503  0.73 +/- 0.03   524
    7   271  0.57 +/- 0.05   191
    8   144  0.73 +/- 0.07   149
    9    56  0.70 +/- 0.10    63
   10    30  0.65 +/- 0.15    24
   11    18  0.16 +/- 0.18     1
   12     4  0.97 +/- 0.70     6
   13     2  0.97 +/- 1.05     4
                             962
      Searching for best position, skipping first 0 layers
Gradient   delta     ra        dec       error     Ts
   17077  0.0334    20.0000   0.0000    0.0106    962.3884
    3204  0.0053    20.0190   -0.0254   0.0098    965.9061
     842  0.0013    20.0173   -0.0204   0.0095    966.1709
     204  0.0003    20.0174   -0.0217   0.0095    966.2021
    *** good fit ***

Spectrum of source three at ra, dec=40.000000, 0.000000
level events  sig fraction  TS
    6    52  0.83 +/- 0.09    74
    7    28  0.44 +/- 0.16    10
    8     9  1.00 +/- 0.44    23
    9     5  1.00 +/- 0.57    14
   10     1  0.00 +/- 13.73    -0
   11     3  0.62 +/- 0.37     3
   12     0
   13     0
                             125
      Searching for best position, skipping first 0 layers
Gradient   delta     ra        dec       error     Ts
    1851  0.0448    40.0000   0.0000    0.0372    124.7921
     651  0.0073    40.0379   0.0239    0.0253    130.3548
                                                 130.3201 --back off
     351  0.0040    40.0348   0.0219    0.0255    130.4250
                                                 130.3594 --back off
                                                 130.4144 --back off
    *** good fit ***
@endverbatim
*/
