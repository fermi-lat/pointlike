/** @file ParamOptimization.h 
    @brief declaration of the ParamOptimization class for optimizing point spread parameters

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/ParamOptimization.h,v 1.9 2008/06/28 06:17:25 mar0 Exp $
*/
#ifndef POINTLIKE_PARAMOPTIMIZATION_H
#define POINTLIKE_PARAMOPTIMIZATION_H
#include "pointlike/PointSourceLikelihood.h"
#ifdef OLD
#include "skymaps/PhotonMap.h"
#else
#include "skymaps/BinnedPhotonData.h"
#endif
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
 //   ParamOptimization(const skymaps::PhotonMap &data, const std::vector<astro::SkyDir>& directions, std::ostream* out=&std::cout, int minlevel=6, int maxlevel=13);
  
    ParamOptimization(const skymaps::BinnedPhotonData &data, const std::vector<astro::SkyDir>& directions, std::ostream* out=&std::cout);

    /** @brief computes and stores parameter value which maximizes overall likelihood
        @param p parameter to optimize; SIGMA or GAMMA from the point spread function
    */
    std::vector<double> compute();

    /** @brief fits the values of sigma to the parameterization described in skymaps::IParams
    */
    std::vector<double> fit_sigma();

private:
#ifdef OLD
    /// @brief Numerical Recipes algorithm for finding the minimum
    double goldensearch(int num2look, int band, bool sigma);
#endif

    /// @brief return the error estimates on the parameters
    std::vector<double> curvature(std::vector<double>& val);

    /// @brief returns the log likelihood of [sigma,gamma]
    static double loglike(std::vector<double> &params);

    //minimum bracketing
    static void brak(double& ax, double& bx, double &cx, double& fa, double& fb, double& fc, double func(const double));
    //brents's method
    static double brent(const double ax, const double bx, const double cx, double f(const double), const double tol, double& xmin);
    //minimization in 1D
    static double f1dim(const double x);
    //minimization on a line
    static void linmin(std::vector<double> &p, std::vector<double> &xi, double &fret, double func(std::vector<double>&));
    //powell's minimization method
    static void powell(std::vector<double> &p,std::vector<std::vector<double> > &xi, const double ftol,int &iter,double &fret,double func(std::vector<double>&));
    //returns chisq of sigma parameterization fit
    static double chisq(std::vector<double> &params);

    std::ostream * m_out;                         //where to send output
#ifdef OLD
    const skymaps::PhotonMap m_data;            //points to skymap
#else
    const skymaps::BinnedPhotonData m_data;
#endif

};

}
#endif
