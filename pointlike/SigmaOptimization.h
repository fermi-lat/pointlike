/** @file SigmaOptimization.h 
    @brief declaration of the SigmaOptimization class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/CalData.h,v 1.3 2007/07/14 03:50:54 burnett Exp $
*/

#include "pointlike/PointSourceLikelihood.h"
#include "map_tools/PhotonMap.h"
#include <iostream>
#include <iomanip>

namespace pointlike {
/**
    @class SigmaOptimization
    @brief calculates the optimum resolution value, sigma, for each healpix level

    */
class SigmaOptimization {
public:
    /** @brief ctor.
        @param data is the PhotonMap data
        @param directions is an array of reference directions to calculate sigma values
        @param radius is the radius cone in degrees
    */
    SigmaOptimization(const map_tools::PhotonMap &data, const std::vector<astro::SkyDir>& directions, std::ostream* out=&std::cout, int minlevel=6, int maxlevel=13, double radius=7.0);
    
private:
    // Numerical Recipes algorithm for finding the minimum
    double goldensearch(const std::vector<astro::SkyDir>& directions, int num2look, int level, double radius);

    const map_tools::PhotonMap m_data;
    int m_minlevel;
    int m_maxlevel;
};

}