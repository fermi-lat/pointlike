/** @file likelihood_main.cxx


*/
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
#include "embed_python/Module.h"

#include "skymaps/BinnedPhotonData.h"
#include "skymaps/PhotonBinner.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/Exposure.h"
#include "skymaps/IsotropicPowerLaw.h"

#include "facilities/commonUtilities.h"

#include "pointlike/SimpleTSmap.h"
#include "pointlike/ParamOptimization.h"
#include "pointlike/ConfidenceLevel.h"
#include "PhotonList.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <typeinfo>
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

// test data file, generated using obsSim with 3 sources
std::string inputFile(facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("pointlike"), 
                     "test_events.root" ));

std::string fitsFile(facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("pointlike"),
					 "test_pass8_ft1.fits"));
/** @class AddPhoton

*/
template<class M>
class AddPhoton: public std::unary_function<skymaps::Photon, void> {
public:
    AddPhoton (M& map)
        : m_map(map)
    {}
    void operator()(const skymaps::Photon& gamma)
    {
        //int event_class = gamma.eventClass();

        m_map.addPhoton(gamma);
    }

    M& m_map;
};



int main(int argc , char** argv )
{

    int rc(0);
    try{

        // use the python module to setup
        std::string python_path("../python");
        
        PointSourceLikelihood::setParameters(embed_python::Module(python_path , "test_pointlike_setup",  argc, argv));

        double  radius(10);


        std::cout << "Loading data from file " << inputFile  <<std::endl;
        PhotonList photons(inputFile);

        double bins_per_decade(5.0);
        skymaps::BinnedPhotonData& x= *new skymaps::BinnedPhotonData(bins_per_decade);
        std::for_each(photons.begin(), photons.end(), AddPhoton<skymaps::BinnedPhotonData>(x));
		//for(PhotonListIterator it=photons.begin();it!=photons.end();++it){
		//	skymaps::Band b = x.binner()(*it);
		//	std::cout << b.event_class() << " " << b.emin() << " " << b.emax() << " " << (b.event_class() + 10*static_cast<int>(b.emin()+.5)) << " " << (int)(b) << std::endl;
		//}

        x.info();

        std::cout << "Loading data from file " << fitsFile  <<std::endl;
		std::vector<std::string> ft1_list;
		ft1_list.push_back(fitsFile);
		Data test_data(ft1_list);
		test_data.info();


#if 0 // disable this since that file does not have GTI
        { 
            std::vector<std::string>ft1_list;
            ft1_list.push_back(facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("pointlike"),"r0236198321_ft1.fit"));
            Data test(ft1_list); 
            test.info();
        }
#endif

#if 0 // test code for SimpleTSmap
        std::cout << "\ntesting SimpleTSmap..." << std::endl;
        std::string path( ::getenv("EXTFILESSYS"));
        skymaps::DiffuseFunction df1( path + "/galdiffuse/GP_gamma.fits");
        PointSourceLikelihood::set_diffuse(&df1);
        PointSourceLikelihood::set_tolerance(0.1);
        PointSourceLikelihood::set_energy_range(500,500); // only the bin with these energies
//        skymaps::BinnedPhotonData pmap("f:\\glast\\data\\SC2\\obssim\\allsky_noGRBs.fits","PHOTONMAP");
        SimpleTSmap tsmap(x);
        tsmap.run(SkyDir(), 180.);
        tsmap.save("tsmap.txt");
        std::cout << "done" << std::endl;
#endif   


        std::vector<astro::SkyDir> directions;
#if 0 //?
        // test fitting with a diffuse component
        skymaps::SkySpectrum* diffuse(new skymaps::IsotropicPowerLaw());
        PointSourceLikelihood::set_diffuse(diffuse);
#endif
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
            std::vector<double> energies=like.energyList();

            if( like.TS()>10) {
               like.localize();
            }
            directions.push_back(like.dir());  // for PSF fits
            
            // check the ts evaluations
            double ts(like.TS()), ts1(like.TSmap(dir)), 
                ts2(like.TSmap(like.dir())); 
            std::cout << "full TS at maximum, starting point: " << ts << ", " 
                << ts1 << ", " << ts2 << std::endl;

            // check confidence level(s)
            double cl = ConfidenceLevel(*like[0], like.dir())();
            std::cout << "Confidence level: " << cl << std::endl;

        };

#if 1 // now refit the second one with the first as background
    PointSourceLikelihood like1(x, points[0].name, directions[0]);
    like1.maximize();
    PointSourceLikelihood like2(x, points[1].name, directions[1]);
    double ts_before(like2.maximize());
    like2.addBackgroundPointSource(&like1);
    double ts_after(like2.maximize());
    std::cout << "TS before, after background:"<<
        ts_before << ", " << ts_after << std::endl;
#endif




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

runing d:\Users\burnett\pointlike-dev\pointlike\v5r4\VC8debug\test_pointlike.exe
Loading data from file ../src/test/test_events.root
photons found: 14396  pixels created: 3226

Spectrum of source one at ra, dec=0, 0
level events   signal_fract  TS
    6  5306  0.74 +/- 0.01  8248
    7  3294  0.70 +/- 0.01  4340
    8  1601  0.83 +/- 0.01  3217
    9   712  0.80 +/- 0.02  1317
   10   248  0.85 +/- 0.03   591
   11   136  0.97 +/- 0.04   388
   12    52  1.00 +/- 0.11   162
   13    31  0.89 +/- 0.10    62
        TS sum  (alpha>0)  18326

Spectrum of source one at ra, dec=0.000000, 0.000000
level events   signal_fract  TS
    6  5306  0.74 +/- 0.01  8248
    7  3294  0.70 +/- 0.01  4340
    8  1601  0.83 +/- 0.01  3217
    9   712  0.80 +/- 0.02  1317
   10   248  0.85 +/- 0.03   591
   11   136  0.97 +/- 0.04   388
   12    52  1.00 +/- 0.11   162
   13    31  0.89 +/- 0.10    62
        TS sum  (alpha>0)  18326
      Searching for best position, start at level 7
Gradient   delta     ra        dec       error     Ts
   11852  0.0020    0.0000    0.0000    0.0031    10078.7
     398  0.0001    0.0020    -0.0001   0.0031    10083.7
    *** good fit ***

Spectrum of source one at ra, dec=0.002003, -0.000120
level events   signal_fract  TS
    6  5306  0.74 +/- 0.01  8248
    7  3292  0.70 +/- 0.01  4345
    8  1601  0.83 +/- 0.01  3217
    9   712  0.80 +/- 0.02  1317
   10   248  0.85 +/- 0.03   591
   11   136  0.97 +/- 0.04   389
   12    52  1.00 +/- 0.11   162
   13    31  0.89 +/- 0.11    63
        TS sum  (alpha>0)  18331

Spectrum of source two at ra, dec=20.000000, 0.000000
level events   signal_fract  TS
    6   564  0.74 +/- 0.03   900
    7   304  0.66 +/- 0.04   370
    8   162  0.92 +/- 0.04   411
    9    66  0.86 +/- 0.07   148
   10    34  0.87 +/- 0.10    73
   11    18  0.98 +/- 0.35    33
   12     4  1.00 +/- 0.60    16
   13     2  1.00 +/- 0.81     9
        TS sum  (alpha>0)  1957

Spectrum of source two at ra, dec=20.000000, 0.000000
level events   signal_fract  TS
    6   564  0.74 +/- 0.03   900
    7   304  0.66 +/- 0.04   370
    8   162  0.92 +/- 0.04   411
    9    66  0.86 +/- 0.07   148
   10    34  0.87 +/- 0.10    73
   11    18  0.99 +/- 0.35    33
   12     4  1.00 +/- 0.60    16
   13     2  1.00 +/- 0.81     9
        TS sum  (alpha>0)  1957
      Searching for best position, start at level 7
Gradient   delta     ra        dec       error     Ts
   18973  0.0319    20.0000   0.0000    0.0098    1057.6
    3478  0.0054    20.0183   -0.0261   0.0095    1061.4
     526  0.0008    20.0150   -0.0218   0.0093    1061.7
    *** good fit ***

Spectrum of source two at ra, dec=20.015227, -0.022581
level events   signal_fract  TS
    6   564  0.74 +/- 0.03   901
    7   305  0.65 +/- 0.04   366
    8   162  0.92 +/- 0.04   411
    9    66  0.85 +/- 0.07   145
   10    35  0.82 +/- 0.10    70
   11    18  1.00 +/- 0.34    39
   12     4  1.00 +/- 0.55    21
   13     2  1.00 +/- 0.78    10
        TS sum  (alpha>0)  1963

Spectrum of source three at ra, dec=40.000000, 0.000000
level events   signal_fract  TS
    6    65  0.67 +/- 0.08    99
    7    33  0.52 +/- 0.13    24
    8    14  0.66 +/- 0.16    27
    9     7  0.74 +/- 0.19    17
   10     1  0.98 +/- 1.34     3
   11     3  1.00 +/- 0.48     8
   12     0
   13     0
        TS sum  (alpha>0)  178

Spectrum of source three at ra, dec=40.000000, 0.000000
level events   signal_fract  TS
    6    65  0.67 +/- 0.08    99
    7    33  0.52 +/- 0.13    24
    8    14  0.66 +/- 0.16    27
    9     7  0.74 +/- 0.21    17
   10     1  0.99 +/- 1.37     3
   11     3  1.00 +/- 0.48     8
   12     0
   13     0
        TS sum  (alpha>0)  178
      Searching for best position, start at level 7
Gradient   delta     ra        dec       error     Ts
    1602  0.0385    40.0000   0.0000    0.0371    79.0
     388  0.0069    40.0248   0.0294    0.0319    80.3
      68  0.0011    40.0315   0.0279    0.0310    80.3
    *** good fit ***

Spectrum of source three at ra, dec=40.031499, 0.029067
level events   signal_fract  TS
    6    65  0.67 +/- 0.08   100
    7    33  0.52 +/- 0.13    24
    8    14  0.66 +/- 0.16    26
    9     7  0.75 +/- 0.21    17
   10     1  1.00 +/- 1.19     4
   11     3  0.82 +/- 0.36    10
   12     0
   13     0
        TS sum  (alpha>0)  180
Computing optimum sigma values
Level     sigma          error          k (~1)    alpha         photons
***************************************************************************
6         0              0              1         1              5935
7         0              0              1         1              3630
8         1              0              1         1              1777
9         1              0              1         1              785
10        1              0              1         1              284
11        1              0              1         1              157
12        1              0              1         1              56
13        1              0              1         1              33
Press any key to continue
@endverbatim
*/
