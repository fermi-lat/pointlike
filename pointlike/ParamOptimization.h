/** @file ParamOptimization.h 
    @brief declaration of the ParamOptimization class for optimizing point spread parameters

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/ParamOptimization.h,v 1.1 2007/11/01 21:49:59 mar0 Exp $
*/
#ifndef POINTLIKE_PARAMOPTIMIZATION_H
#define POINTLIKE_PARAMOPTIMIZATION_H
#include "pointlike/PointSourceLikelihood.h"
#include "map_tools/PhotonMap.h"
#include <iostream>
#include <iomanip>

namespace pointlike {
/**
    @class ParamOptimization
    @brief calculates the optimum resolution value, sigma, or the tail paramater, gamma, for each healpix level

    */
class ParamOptimization {

public:
    //parameters
    typedef enum {SIGMA, 
        GAMMA} Param;

    /** @brief ctor.
        @param data is the PhotonMap data
        @param directions is an array of reference directions to calculate sigma values
        @param radius is the radius cone in degrees
    */
    ParamOptimization(const map_tools::PhotonMap &data, const std::vector<astro::SkyDir>& directions, std::ostream* out=&std::cout, int minlevel=6, int maxlevel=13);
    
    //returns signal fraction for each energy bin
    std::vector<double> get_alphas() {return m_alphas;}

    /** @brief computes and stores parameter value which maximizes overall likelihood
        @param p parameter to optimize; SIGMA or GAMMA from the point spread function
    */
    void compute(Param p=SIGMA);

private:
    // Numerical Recipes algorithm for finding the minimum
    double goldensearch(int num2look, int level, bool sigma);
    std::vector<double> m_alphas;                 //signal fractions for each energy bin
    std::vector<pointlike::PointSourceLikelihood*> m_likes; 
    std::ostream * m_out;                         //where to send output
    const map_tools::PhotonMap m_data;            //points to skymap
    int m_minlevel;                               //minimum healpix level
    int m_maxlevel;                               //maximum healpix level
};

}
#endif
