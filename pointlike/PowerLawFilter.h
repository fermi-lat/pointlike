/** @file PowerLawFilter.h
@brief declare class PowerLawFilter

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/PowerLawFilter.h,v 1.1 2007/07/14 03:50:54 burnett Exp $
*/

#ifndef tools_PowerLawFilter_h
#define tools_PowerLawFilter_h

#include <vector>
#include <cmath>


namespace pointlike {

    /**
    @class DataPoint
    @brief Used to input a data point to PowerLawFilter class

    */
    class DataPoint {
    public:

        /** @brief ctor store a data point and the log of input values
        @param x_in x coordinate
        @param y_in y coordinate
        @param sigma_in measurement error
        */
        DataPoint(double x_in, double y_in, double sigma_in):
          m_x(x_in),
              m_y(y_in),
              m_sigma(sigma_in),
              m_ln_x(log(x_in)),
              m_ln_y(log(y_in)),
              m_ln_sigma(log(sigma_in))
          {}

          double x()const {return m_x;} 
          double y()const {return m_y;}
          double sigma()const {return m_sigma;}
          double ln_x()const {return m_ln_x;} 
          double ln_y()const {return m_ln_y;}
          double ln_sigma()const {return m_ln_sigma;} 

    private:
        double m_x;  
        double m_y; 
        double m_sigma; 
        double m_ln_x;  
        double m_ln_y; 
        double m_ln_sigma; 

    };


    /**
    @class PowerLawFilter
    @brief Fit a set of points to a power law

    */
    class PowerLawFilter {
    public:

        /** @brief ctor accepts input data and evaluates fit parameters.
        @param values A vector of pairs of (x, y) values
        */
        PowerLawFilter(const std::vector<std::pair<double, double> > & values);

        double slope()const {return m_slope;} 
        double constant()const {return m_constant;}
        double metric()const {return m_metric;}
        double chi_sq()const {return m_chi_sq;}

    private:
        double m_slope;  // m, from y = mx + b
        double m_constant; // b, from y = mx + b
        double m_metric; // How good is the fit?
        double m_chi_sq; // Chi squared?

    };



} // namespace pointlike

#endif
